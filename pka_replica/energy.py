"""Fast GFN2-xTB/ALPB microstate energetics (protocol v4).

Speed design (after profiling on methotrexate-class molecules):

  * ONE parent 3D geometry is generated (ETKDG + MMFF, best of 4) and every
    protonation state is derived from it by adding/removing H atoms on the
    same coordinates — no per-state embedding, no per-state MMFF.
  * All xtb jobs run through a ProcessPool (OMP_THREADS=1 per worker).
  * Rowan-style screening gate: single-flip states are computed first and
    calibrated; sites whose calibrated pKa falls outside [min-5, max+5]
    (stjames `reasonableness_buffer = 5`) are dropped before any multi-flip
    state is computed; multi-flip states are capped (16 cheapest) and use
    single-conformer solvated optimizations.
  * Neutral parents still get gas-phase optimization + solvated single point
    (avoids the zwitterion collapse); charged states are optimized directly
    in ALPB water (screens spurious intramolecular salt bridges).

Accuracy design:
  * site terms: functional-group-specific linear calibration (as Rowan);
  * state term: exact GFN2 many-body residual clipped to ±4 pKa units
    (larger values are geometry noise / ALPB artifacts).
"""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from .enumerator import (ALCOHOL_FLOOR, NON_GENUINE_ACID_N, Microstate,
                         acid_sites, base_sites)

LN10 = math.log(10.0)
RT_KCAL = 0.0019872041 * 298.15
HARTREE_KCAL = 627.509474
PKA_UNIT_KCAL = LN10 * RT_KCAL  # 1.3644 kcal/mol per pKa unit

# Empirical linear calibration pKa = a * raw + b (multistage protocol fit:
# AIMNet2 preopt + GFN2 opt + GFN2/ALPB single point; see calibrate_fit.py).
# raw = GFN2/ALPB deprotonation energy in pKa units, always in the
# deprotonation direction (base sites use the acidity of the conjugate acid).
CAL_TABLE: dict[str, tuple[float, float]] = {
    "O-acid": (0.553985, -46.321),
    "O-acid-fluoro": (0.161221, -11.303),  # refit 2026-10-24, 4-ref F ladder
    "N-acid": (0.489772, -37.585),
    "S-acid": (0.758134, -64.215),
    "N-base-ali": (0.398278, -24.542),   # refit 2026-10-22, 14 refs, slope positive
    "N-base-aro": (0.360147, -24.931),
    "N-base-aniline": (1.000000, -75.455),
    "N-base": (0.756512, -57.413),       # pooled fallback (legacy protocol)
}

# Literature default micro-pKas per SMARTS template (fragment fallback route).
FRAGMENT_PKA: dict[str, float] = {
    "ammonium_NH": 9.25,
    "carboxylic_acid": 4.3,
    "sulfonic_acid": -1.0,
    "phosphonic_acid": 2.1,
    "phosphoric_acid": 2.1,
    "boronic_acid": 8.8,
    "sulfonamide_NH": 10.0,
    "tetrazole_NH": 4.9,
    "imide_NH": 9.9,
    "hydroxamic_acid": 8.7,
    "phenol_enol_oxime": 10.0,
    "aromatic_NH": 14.0,
    "amide_NH": 15.0,
    "anilide_like_NH": 16.0,
    "thiol_sh": 10.3,
    "thiourea_NH": 21.0,
    "hydrazine_NH": 15.0,
    "alcohol_OH": 16.0,
    "amidine_guanidine_N": 12.5,
    "aromatic_N_acceptor": 5.2,
    "aniline_N": 4.9,
    "aliphatic_amine_N": 9.9,
    "nitrile_N": -9.0,
    "amide_N_protonation": -0.5,
}


def _flip_class(atom: Chem.Atom, delta: int, mol: Chem.Mol | None = None) -> str:
    sym = atom.GetSymbol().upper()
    if delta < 0:
        if sym == "N" and atom.GetFormalCharge() > 0:
            # ammonium N-H (zwitterion inputs): deprotonation of N+H IS the
            # conjugate-acid acidity of an amine — the N-base-ali domain,
            # NOT N-acid (fitted on sulfonamide/imide N-H, raws ~30-60)
            return "N-base-ali"
        if sym == "O" and mol is not None:
            # fluorinated carboxylic acids are FAR stronger than the generic
            # O-acid calibration domain (TFA 0.2 vs acetic 4.8) — route them
            # to their own class when any alpha carbon of the acid group
            # carries F (covers CH2F .. CF3 and perfluoro chains)
            for c in atom.GetNeighbors():
                if c.GetSymbol() != "C":
                    continue
                oxo = any(n.GetSymbol() == "O" and n.GetIdx() != atom.GetIdx()
                          for n in c.GetNeighbors())
                if not oxo:
                    continue
                for alpha in c.GetNeighbors():
                    if alpha.GetSymbol() == "C" and any(
                            n.GetSymbol() == "F" for n in alpha.GetNeighbors()):
                        return "O-acid-fluoro"
        return f"{sym}-acid"
    if delta > 0:
        if sym == "O":
            # carboxylate/sulfonate O- protonation = conjugate-acid acidity:
            # the O-acid calibration domain (same raws as acid dissociation)
            return "O-acid"
        if sym == "N":
            if atom.GetIsAromatic():
                return "N-base-aro"
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 6 and nb.GetIsAromatic():
                    return "N-base-aniline"
            return "N-base-ali"
        return f"{sym}-base"
    return "N-base"  # delta == 0: shouldn't occur; pooled fallback is safest


# ---------------------------------------------------------------------------
# xtb worker (parallel)
# ---------------------------------------------------------------------------

_XTB_EXE = shutil.which("xtb")

# AIMNet2 geometry engine (isayevlab/aimnetcentral, wB97M-D3-trained NNP).
# One calculator instance per worker thread; model load is cheap (~0.1 s).
import threading  # noqa: E402

_AIMNET_LOCAL = threading.local()
_AIMNET_LOCK = threading.Lock()  # warp's CPU tile arena is not thread-safe
AIMNET_AVAILABLE: bool | None = None


def _aimnet_calc():
    global AIMNET_AVAILABLE
    calc = getattr(_AIMNET_LOCAL, "calc", None)
    if calc is None:
        from aimnet.calculators import AIMNet2ASE

        calc = AIMNet2ASE("aimnet2", charge=0)
        _AIMNET_LOCAL.calc = calc
        AIMNET_AVAILABLE = True
    return calc


def aimnet_optimize(symbols: list[str], coords: list[list[float]], charge: int,
                    fmax: float = 0.1, steps: int = 400) -> list[list[float]] | None:
    """Gas-phase geometry optimization with the AIMNet2 neural potential.
    Returns optimized coordinates or None on failure.

    fmax is intentionally loose (0.1 eV/A): this stage only hands GFN2 a
    physically sensible starting geometry — GFN2 then converges to its own
    minimum, and an A/B on MTX confirmed site pKas shift <=0.05 vs fmax=0.05
    while the preopt gets ~20% cheaper.

    Charge is passed via atoms.info['charge'] (the documented highest-
    precedence channel) so a shared per-thread calculator keeps no stale
    state between calls.
    """
    try:
        from ase import Atoms
        from ase.optimize import BFGS

        Z = {"H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Si": 14, "P": 15,
             "S": 16, "Cl": 17, "As": 33, "Se": 34, "Br": 35, "I": 53}
        numbers = [Z[s.capitalize()] for s in symbols]
        atoms = Atoms(numbers=numbers, positions=coords)
        atoms.info["charge"] = float(charge)
        atoms.calc = _aimnet_calc()
        # warp's CPU runtime asserts under concurrent tile ops, so AIMNet2
        # optimization is serialized across worker threads (the xtb single
        # points around it still run in parallel)
        with _AIMNET_LOCK:
            BFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)
        return [list(map(float, p)) for p in atoms.get_positions()]
    except Exception:
        return None


def _xtb_job(job: dict, want_props: bool = False):
    """Run one energy job. `job` = {symbols, coords, charge, timeout}.

    Multistage protocol (Rowan-style; deterministic final geometry):
    1. AIMNet2 gas-phase preoptimization from the seed (fast NNP gets close);
    2. short GFN2-xTB gas-phase optimization — converges in a few steps to
       the *GFN2* minimum (sampling a different engine's minimum injects
       2-3 pKa of noise into the single-point differences);
    3. GFN2-xTB + ALPB-water single point (the solvation correction;
       optimizing in implicit solvent instead both shifts the calibration
       and collapses neutral amino acids into zwitterions).
    Returns (solvated single-point energy in Hartree or None, props dict).
    """
    exe = job.get("exe") or _XTB_EXE
    if not exe:
        return None, {}
    symbols, coords, charge = job["symbols"], job["coords"], job["charge"]
    timeout = job.get("timeout", 120)
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_STACKSIZE="1G")

    def write_xyz(path: str, syms, xyz) -> None:
        lines = [str(len(syms)), ""]
        for s, (x, y, z) in zip(syms, xyz):
            lines.append(f"{s} {x:.5f} {y:.5f} {z:.5f}")
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    with tempfile.TemporaryDirectory(prefix="pka_xtb_") as td:
        src = ""
        try:
            inp = os.path.join(td, "in.xyz")
            write_xyz(inp, symbols, coords)
            # 1) AIMNet2 preoptimization (best-effort; fall back to raw seed)
            pre = aimnet_optimize(symbols, coords, charge)
            if pre is not None:
                write_xyz(inp, symbols, pre)
            # 2) GFN2 gas-phase optimization (short: starts pre-converged)
            subprocess.run(
                [exe, inp, "--gfn", "2", "--opt", "--chrg", str(charge), "--norestart"],
                capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=timeout, cwd=td, env=env,
            )
            opt = os.path.join(td, "xtbopt.xyz")
            if not os.path.exists(opt):
                return None, {}
            # 3) solvated single point
            r = subprocess.run(
                [exe, opt, "--gfn", "2", "--alpb", "water", "--sp", "--chrg", str(charge)],
                capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=timeout, cwd=td, env=env,
            )
            src = r.stdout
        except Exception:
            return None, {}
        hits = re.findall(r"total energy\s+(-?\d+\.\d+)\s+Eh", src)
        energy = float(hits[-1]) if hits else None
        props: dict = {"total_energy_Eh": energy}
        m = re.search(r"HOMO-LUMO gap\s+(-?\d+\.\d+)\s+eV", src)
        if m:
            props["homo_lumo_gap_eV"] = float(m.group(1))
        m = re.search(r"-> Gsolv\s+(-?\d+\.\d+)\s+Eh", src)
        if m:
            props["Gsolv_kcal_mol"] = round(float(m.group(1)) * HARTREE_KCAL, 2)
        m = re.search(
            r"full:\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", src)
        if m:
            props["dipole_D"] = float(m.group(4))
        return energy, props


def _parent_geometry_from_mol(mol: Chem.Mol, n_embed: int = 4) -> tuple[list[str], list[list[float]]] | None:
    """ETKDG + MMFF on the *given* mol object (its atom ordering is what the
    enumeration flips refer to); return symbols and best-conformer coords (H included)."""
    m = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=n_embed, params=params))
    if not cids:
        return None
    try:
        scores = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=500)
        energies = {cid: e for cid, (_, e) in zip(cids, scores)}
        best = min(cids, key=lambda c: energies.get(c, 0.0))
    except Exception:
        best = cids[0]
    conf = m.GetConformer(best)
    symbols = [a.GetSymbol() for a in m.GetAtoms()]
    return symbols, [list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())]


def _parent_geometry(smiles: str, n_embed: int = 4):
    """SMILES convenience wrapper (atom order = parsed SMILES order)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return _parent_geometry_from_mol(mol, n_embed)


def _apply_flips_to_geometry(
    symbols: list[str], coords: list[list[float]],
    parent_mol_h: Chem.Mol, flips: dict[int, int],
) -> tuple[list[str], list[list[float]]] | None:
    """Derive a state's geometry from the parent by adding/removing H atoms.

    Deletes never reorder surviving atoms (index stability is what keeps
    element/coordinate pairing intact); additions append at the end.
    flips keys are heavy-atom indices of the parent molecule.
    """
    n_atoms = parent_mol_h.GetNumAtoms()
    # H atoms attached to each heavy atom (indices in the full parent list)
    h_map: dict[int, list[int]] = {}
    for atom in parent_mol_h.GetAtoms():
        if atom.GetSymbol() != "H":
            continue
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() != "H":
                h_map.setdefault(nb.GetIdx(), []).append(atom.GetIdx())
                break

    to_remove: set[int] = set()
    to_add: list[int] = []

    cx = sum(p[0] for p in coords) / len(coords)
    cy = sum(p[1] for p in coords) / len(coords)
    cz = sum(p[2] for p in coords) / len(coords)

    for idx, delta in flips.items():
        if idx >= n_atoms or parent_mol_h.GetAtomWithIdx(idx).GetSymbol() == "H":
            return None  # flips must address heavy atoms of the parent
        if delta < 0:
            hs = [h for h in h_map.get(idx, []) if h not in to_remove]
            if len(hs) < -delta:
                return None  # not enough H -> state impossible
            for h in hs[: -delta]:
                to_remove.add(h)
        elif delta > 0:
            to_add.extend([idx] * delta)

    syms = [s for i, s in enumerate(symbols) if i not in to_remove]
    xyz = [list(p) for i, p in enumerate(coords) if i not in to_remove]

    for idx in to_add:
        px, py, pz = coords[idx]
        r = math.sqrt(max((px - cx) ** 2 + (py - cy) ** 2 + (pz - cz) ** 2, 0.25))
        xyz.append([px + 1.09 * (px - cx) / r,
                    py + 1.09 * (py - cy) / r,
                    pz + 1.09 * (pz - cz) / r])
        syms.append("H")
    return syms, xyz


def calibrate(cls: str, raw_pka_units: float, lo: float = -2.0, hi: float = 16.0) -> float:
    """Functional-group-specific linear calibration, clamped to [lo, hi]
    (replica of Rowan's reasonableness screen)."""
    a, b = CAL_TABLE.get(cls, CAL_TABLE["N-base"])
    return max(lo, min(hi, a * raw_pka_units + b))


# ---------------------------------------------------------------------------
# quantum energy model
# ---------------------------------------------------------------------------

class QuantumEnergyModel:
    """Screened, parallel GFN2/ALPB microstate free energies (protocol v4)."""

    name = "xtb-gfn2-alpb"

    def __init__(self, n_confs: int = 2, max_workers: int | None = None,
                 timeout_s: int = 180, max_multi_states: int = 16):
        if not _XTB_EXE:
            raise RuntimeError("xtb executable not found on PATH")
        self.n_confs = n_confs
        self.max_workers = max_workers or max(1, (os.cpu_count() or 4) - 1)
        self.timeout_s = timeout_s
        self.max_multi_states = max_multi_states
        self.last_site_pka: dict[tuple[int, int], float] = {}

    def compute(
        self,
        mol: Chem.Mol,
        states: list[Microstate],
        pka_window: tuple[float, float] = (2.0, 12.0),
        buffer: float = 5.0,
        progress=None,
        with_multi: bool = True,
    ) -> list[float]:
        """Return dimensionless g_i for every state (input state ≈ 0 reference).

        progress: optional callback(done_jobs, total_jobs, phase).
        """
        parent = next((s for s in states if not s.flips), states[0] if states else None)
        if parent is None:
            return []
        # geometry MUST be generated from the same mol object the flip indices
        # refer to (canonical-SMILES reparsing silently permutes atom order)
        geo = _parent_geometry_from_mol(mol, n_embed=max(self.n_confs, 2))
        if geo is None:
            raise RuntimeError("3D 构象生成失败(ETKDG 嵌入无解),无法计算 pKa")
        symbols, coords = geo
        molH = Chem.AddHs(Chem.Mol(mol))
        z0 = Chem.GetFormalCharge(mol)

        def mk(syms, xyz, charge):
            return dict(symbols=syms, coords=xyz, charge=charge,
                        timeout=self.timeout_s, exe=_XTB_EXE)

        jobs: list[tuple[str, dict]] = []
        # 1) parent: single job — _parent_geometry_from_mol already embedded
        # n_confs conformers and kept the MMFF-best, so additional identical
        # parent jobs would deterministically return the same energy
        jobs.append(("parent|0", mk(symbols, coords, z0)))
        # 2) single-flip states
        single_states = [s for s in states if len(s.flips) == 1]
        for s in single_states:
            g = _apply_flips_to_geometry(symbols, coords, molH, s.flips)
            if g is not None:
                jobs.append((f"single|{s.smiles}", mk(g[0], g[1], s.charge)))

        total1 = len(jobs)
        done = 0
        results: dict[str, float | None] = {}
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futs = {pool.submit(_xtb_job, j): tag for tag, j in jobs}
            for fut, tag in futs.items():
                try:
                    results[tag] = fut.result()[0]
                except Exception:
                    # one broken worker (e.g. Windows temp-dir cleanup race)
                    # must not discard every other completed job
                    results[tag] = None
                done += 1
                if progress:
                    progress(done, total1, "single-flip screening")

        parent_es = [results.get(f"parent|{k}") for k in range(1)]
        parent_es = [e for e in parent_es if e is not None]
        if not parent_es:
            # silent all-zero energies would produce confidently wrong pKas
            # downstream; fail loudly instead
            raise RuntimeError(
                "xtb 母体能量计算失败(几何生成或全部任务超时/崩溃),无法给出 pKa")
        g_parent = parent_es[0]

        site_pka: dict[tuple[int, int], float] = {}
        site_raw: dict[tuple[int, int], float] = {}
        for s in single_states:
            r = results.get(f"single|{s.smiles}")
            if r is None:
                continue
            (idx, delta), = s.flips.items()
            if abs(delta) != 1:
                continue  # multi-unit flip at one atom is a state, not a site
            raw = (r - g_parent) * HARTREE_KCAL / PKA_UNIT_KCAL
            if delta == +1:
                raw = -raw
            cls = _flip_class(mol.GetAtomWithIdx(idx), delta, mol)
            site_raw[(idx, delta)] = raw
            site_pka[(idx, delta)] = calibrate(cls, raw)
        self.last_site_pka = dict(site_pka)
        # template name per site, from the enumeration tag ("acid:amide_NH@21")
        tpl_of: dict[tuple[int, int], str] = {}
        for s in single_states:
            if len(s.flips) == 1 and ":" in s.tag:
                (i0, d0), = s.flips.items()
                tpl_of[(i0, d0)] = s.tag.split(":", 1)[1].rsplit("@", 1)[0]

        def trusted(k: tuple[int, int], v: float) -> bool:
            # weakly-acidic N-H (amide/aniline/indole class) below pKa 10 is a
            # calibration artifact: letting it through used to put impossible
            # deprotonated states at the BOTTOM of the ensemble (vincristine
            # bug) and report -2 as "most acidic"
            if k[1] < 0 and tpl_of.get(k) in NON_GENUINE_ACID_N and v < 10.0:
                return False
            # aliphatic alcohols are >= ~13 in water; lower values (sterically
            # odd OH geometries) are artifacts (vinorelbine 2.2)
            if k[1] < 0 and tpl_of.get(k) == "alcohol_OH" and v < ALCOHOL_FLOOR:
                return False
            # a base with pKaH < 0 is not a base in water; its clamped value
            # understates the protonation cost and makes multi-protonated
            # combos artificially stable (direct-SMILES MTX +2-level bug)
            if k[1] > 0 and v < 0.0:
                return False
            return True

        # screening gate: drop sites outside window (+ buffer), keep the flip
        # direction that the window range actually covers
        lo, hi = pka_window[0] - buffer, pka_window[1] + buffer
        kept = {k: v for k, v in site_pka.items()
                if lo <= v <= hi and trusted(k, v)}
        # single source of truth for the UI: site_legend/n_active_sites must
        # show exactly the sites that participate in the ensemble
        self.last_kept_sites = dict(kept)

        # 3) multi-flip states over kept sites, capped by additive energy
        multi_states = [s for s in states if len(s.flips) > 1]
        ranked: list[tuple[float, Microstate]] = []
        for s in multi_states:
            add = 0.0
            ok = True
            for idx, delta in s.flips.items():
                if delta == 0:
                    continue
                sign = 1.0 if delta < 0 else -1.0
                unit = (idx, -1 if delta < 0 else +1)
                p = kept.get(unit)
                if p is None:
                    ok = False
                    break
                add += sign * abs(delta) * p
            if ok:
                ranked.append((add, s))
        ranked.sort(key=lambda t: t[0])
        chosen = [] if not with_multi else [s for _, s in ranked[: self.max_multi_states]]

        jobs2: list[tuple[str, dict]] = []
        for s in chosen:
            g = _apply_flips_to_geometry(symbols, coords, molH, s.flips)
            if g is not None:
                jobs2.append((f"multi|{s.smiles}", mk(g[0], g[1], s.charge)))
        total2 = total1 + len(jobs2)
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futs = {pool.submit(_xtb_job, j): tag for tag, j in jobs2}
            for fut, tag in futs.items():
                try:
                    results[tag] = fut.result()[0]
                except Exception:
                    results[tag] = None
                done += 1
                if progress:
                    progress(done, total2, "multi-state coupling")

        # assemble energies
        additive: dict[str, float | None] = {}
        raw_additive: dict[str, float] = {}
        for s in states:
            if not s.flips:
                additive[s.smiles] = 0.0
                raw_additive[s.smiles] = 0.0
                continue
            a, raw_sum, ok = 0.0, 0.0, True
            for idx, delta in s.flips.items():
                if delta == 0:
                    continue
                sign = 1.0 if delta < 0 else -1.0
                unit = (idx, -1 if delta < 0 else +1)
                p = kept.get(unit)  # no site_pka fallback: untrusted sites
                # must keep their states OUT of the ensemble entirely
                r = site_raw.get(unit)
                if p is None or r is None:
                    ok = False
                    break
                a += sign * abs(delta) * p
                raw_sum += sign * abs(delta) * r
            additive[s.smiles] = a if ok else None
            raw_additive[s.smiles] = raw_sum

        raw_multi: dict[str, float] = {}
        for s in chosen:
            r = results.get(f"multi|{s.smiles}")
            if r is not None:
                raw_multi[s.smiles] = (r - g_parent) * HARTREE_KCAL / PKA_UNIT_KCAL

        out: list[float] = []
        for s in states:
            if not s.flips:
                out.append(0.0)
                continue
            a = additive.get(s.smiles)
            if a is None:
                out.append(LN10 * 40.0)  # uncomputed multi state -> high energy
                continue
            if s.smiles in raw_multi:
                residual = max(-4.0, min(4.0, raw_multi[s.smiles] - raw_additive[s.smiles]))
                out.append(LN10 * (a + residual))
            else:
                out.append(LN10 * a)
        # re-reference to the lowest state
        if out:
            ref = min(out)
            out = [g - ref for g in out]
        return out


# ---------------------------------------------------------------------------
# fragment route (fast fallback)
# ---------------------------------------------------------------------------

@dataclass
class Site:
    atom_index: int
    kind: str
    template: str
    pka: float
    source: str


class FragmentBackend:
    """Site-pKa oracle from literature defaults per SMARTS template."""

    name = "fragment"

    def site_pkas(self, mol: Chem.Mol) -> tuple[dict[int, float], dict[int, float]]:
        acids = {idx: FRAGMENT_PKA[name] for idx, name in acid_sites(mol)}
        bases = {idx: FRAGMENT_PKA[name] for idx, name in base_sites(mol)}
        return acids, bases


class PkasolverBackend:
    """Optional pkasolver GNN oracle (not installed on this machine)."""

    name = "pkasolver"

    def __init__(self) -> None:
        from pkasolver.prediction import predict_pka  # type: ignore

        self._fn = predict_pka

    def site_pkas(self, mol: Chem.Mol):
        acids, bases = {}, {}
        try:
            rows = self._fn(Chem.Mol(mol))
        except Exception:
            return acids, bases
        for item in rows or []:
            try:
                idx, pka = int(item[0]), float(item[1])
            except Exception:
                continue
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetFormalCharge() > 0 or atom.GetTotalNumHs() == 0:
                bases[idx] = pka
            else:
                acids[idx] = pka
        return acids, bases


class EnergyModel:
    """Site-decomposition route: site pKas + through-space Coulomb coupling."""

    def __init__(self, backend: str = "auto", eps_eff: float = 21.0, coupling: bool = True):
        self.eps_eff = eps_eff
        self.coupling = coupling
        if backend == "fragment":
            self.oracle = FragmentBackend()
        elif backend == "pkasolver":
            self.oracle = PkasolverBackend()
        else:
            try:
                self.oracle = PkasolverBackend()
            except Exception:
                self.oracle = FragmentBackend()
        self.sites: list[Site] = []
        self._conf = None
        self.acid_pka: dict[int, float] = {}
        self.base_pka: dict[int, float] = {}
        self.acid_default: dict[int, float] = {}
        self.base_default: dict[int, float] = {}

    def prepare(self, mol: Chem.Mol) -> None:
        self.acid_pka, self.base_pka = self.oracle.site_pkas(mol)
        self.sites = []
        for idx, name in acid_sites(mol):
            p = self.acid_pka.get(idx, FRAGMENT_PKA.get(name, 14.0))
            self.sites.append(Site(idx, "acid", name, p,
                                   "pkasolver" if idx in self.acid_pka else "fragment"))
        for idx, name in base_sites(mol):
            p = self.base_pka.get(idx, FRAGMENT_PKA.get(name, 5.0))
            self.sites.append(Site(idx, "base", name, p,
                                   "pkasolver" if idx in self.base_pka else "fragment"))
        self._conf = _single_conformer(mol) if self.coupling else None
        self.acid_default = {s.atom_index: s.pka for s in self.sites if s.kind == "acid"}
        self.base_default = {s.atom_index: s.pka for s in self.sites if s.kind == "base"}

    def energy(self, ms: Microstate) -> float:
        g = 0.0
        ionized: list[tuple[int, int]] = []
        for idx, delta in ms.flips.items():
            if delta == -1:
                g += LN10 * self.acid_pka.get(idx, self.acid_default.get(idx, 14.0))
                ionized.append((idx, -1))
            elif delta == +1:
                g -= LN10 * self.base_pka.get(idx, self.base_default.get(idx, 5.0))
                ionized.append((idx, +1))
        if self._conf is not None and len(ionized) > 1:
            for i in range(len(ionized)):
                for j in range(i + 1, len(ionized)):
                    ia, qa = ionized[i]
                    ib, qb = ionized[j]
                    r = max(math.dist(self._conf[ia], self._conf[ib]), 1.0)
                    g += (qa * qb * 332.06371) / (self.eps_eff * r) / RT_KCAL
        return g

    def logp(self, smiles: str) -> float:
        return _logp(smiles)


def _single_conformer(mol: Chem.Mol) -> list[list[float]] | None:
    m = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    try:
        if AllChem.EmbedMolecule(m, params) != 0:
            return None
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:
        return None
    conf = m.GetConformer()
    return [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]


def _logp(smiles: str) -> float:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return 0.0
    try:
        return rdMolDescriptors.CalcCrippenDescriptors(m)[0]
    except Exception:
        return 0.0
