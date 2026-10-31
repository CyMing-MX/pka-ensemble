"""Workflow orchestration — mirrors Rowan's stjames MacropKaWorkflow / pKaWorkflow.

run_macropka(smiles)  -> the "Macroscopic pKa, logD" result object
run_micropka(smiles)  -> the site-level (microscopic) pKa result object

backend:
  "xtb"    — GFN2-xTB/ALPB physics route (needs xtb on PATH); parallel,
             geometry-reusing, Rowan-style screened (default)
  "sites"  — site-pKa decomposition (pkasolver if installed, else fragments)
  "auto"   — "xtb" when xtb is available, else "sites"

charge_window:
  "auto"   — derived from the molecule: the input charge widened by the
             number of active ionizable sites, clamped to ±3 (the replica
             equivalent of Rowan's min/max charge selector, which defaults
             to -2/+2 and allows -7..+7)
"""

from __future__ import annotations

import shutil

from rdkit import Chem

from . import thermo
from .energy import EnergyModel, QuantumEnergyModel, _logp
from .enumerator import (ALCOHOL_FLOOR, NON_GENUINE_ACID_N, Microstate,
                         _edit_charge, acid_sites, base_sites, canonical,
                         enumerate_microstates)


def _ensure_connected(mol: Chem.Mol, warnings: list[str]) -> Chem.Mol:
    """Salts / mixtures: keep the largest fragment by heavy-atom count and
    record a warning. A disconnected ion pair has no meaningful whole-pKa
    (Rowan outright rejects it; we compute the main component instead)."""
    if mol is None:
        return mol
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) <= 1:
        return mol
    best = max(frags, key=lambda m: (m.GetNumHeavyAtoms(), Chem.GetFormalCharge(m)))
    dropped = [Chem.MolToSmiles(m) for m in frags if m is not best]
    warnings.append(
        f"输入含 {len(frags)} 个独立片段(盐/混合物),已取最大片段 "
        f"{Chem.MolToSmiles(best)} 计算;被忽略: {', '.join(dropped)}。"
        "离子对的整体 pKa 无定义。")
    return best


def find_stable_tautomer(smiles: str, energy_fn=None) -> str:
    """QupKake-style tautomer search: RDKit enumerates tautomers, the most
    stable one is selected by GFN2-xTB + ALPB-water total energy (optional;
    falls back to RDKit's canonical ranking when xtb is unavailable)."""
    from rdkit.Chem.MolStandardize import rdMolStandardize

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    tautomers = list(rdMolStandardize.TautomerEnumerator().Enumerate(mol))
    if len(tautomers) <= 1:
        return canonical(tautomers[0] if len(tautomers) else mol)
    if energy_fn is None:
        return canonical(tautomers[0])
    scored = []
    for t in tautomers:
        try:
            e = energy_fn(canonical(t))
        except Exception:
            continue
        if e is None:
            continue
        scored.append((e, canonical(t)))
    if not scored:
        # every energy job failed: keep the input rather than silently
        # swapping to an arbitrary enumeration order (and claiming it
        # was "lower energy")
        return canonical(mol)
    return min(scored)[1]


def _mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"cannot parse SMILES: {smiles!r}")
    return mol


def _resolve_backend(backend: str) -> str:
    if backend != "auto":
        return backend
    return "xtb" if shutil.which("xtb") else "sites"


def suggest_charge_window(mol: Chem.Mol) -> tuple[int, int]:
    """Auto charge window: input charge widened by the count of ionizable
    sites, clamped to ±3."""
    z0 = Chem.GetFormalCharge(mol)
    n = len(acid_sites(mol)) + len(base_sites(mol))
    w = max(2, min(3, n))
    return z0 - w, z0 + w


def mol_properties(smiles: str) -> dict:
    """Curated RDKit physicochemical descriptors for the input state."""
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    mol = _mol(smiles)
    return {
        "MW": round(Descriptors.MolWt(mol), 2),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "logP_crippen": round(Crippen.MolLogP(mol), 2),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 1),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
    }


def run_macropka(
    smiles: str,
    min_pH: float = 0.0,
    max_pH: float = 14.0,
    min_charge: int | str = "auto",
    max_charge: int | str = "auto",
    backend: str = "auto",
    eps_eff: float = 21.0,
    coupling: bool = True,
    max_states: int = 200,
    n_confs: int = 2,
    pka_window: tuple[float, float] = (2.0, 12.0),
    progress=None,
    tautomer_search: bool = False,
) -> dict:
    """Replica of Rowan's macropka workflow result schema (stjames MacropKaWorkflow)."""
    warnings_out: list[str] = []
    mol = _ensure_connected(_mol(smiles), warnings_out)
    tautomer_switched = False
    if tautomer_search and shutil.which("xtb"):
        # QuantumEnergyModel intentionally NOT re-imported here: a conditional
        # import would shadow the module-level name and break the xtb route
        # below when tautomer_search is False (UnboundLocalError).
        from .energy import _XTB_EXE, _xtb_job, _parent_geometry_from_mol

        def xtb_energy(s):
            geo = _parent_geometry_from_mol(Chem.MolFromSmiles(s), n_embed=1)
            if geo is None:
                return None
            # return the bare energy: find_stable_tautomer min() compares
            # these directly (a (energy, props) tuple or None breaks it)
            return _xtb_job({"symbols": geo[0], "coords": geo[1],
                             "charge": Chem.GetFormalCharge(Chem.MolFromSmiles(s)),
                             "timeout": 120})[0]

        best = find_stable_tautomer(canonical(mol), energy_fn=xtb_energy)
        if best != canonical(mol):
            warnings_out.append(f"互变异构体搜索: 已切换到更低能互变体 {best}")
            mol = Chem.MolFromSmiles(best)
            tautomer_switched = True
    route = _resolve_backend(backend)

    if min_charge == "auto" or max_charge == "auto":
        amin, amax = suggest_charge_window(mol)
        if min_charge == "auto":
            min_charge = amin
        if max_charge == "auto":
            max_charge = amax
    min_charge, max_charge = int(min_charge), int(max_charge)

    states: list[Microstate] = enumerate_microstates(
        mol, min_charge=min_charge, max_charge=max_charge, max_states=max_states
    )
    if not states:
        raise ValueError(
            f"电荷窗 [{min_charge:+d}, {max_charge:+d}] 内没有任何微态"
            "(窗口未覆盖输入电荷 z0="
            f"{Chem.GetFormalCharge(mol):+d} 且无可达的质子化状态),请放宽 min_charge/max_charge")
    if not any(not s.flips for s in states):
        # window excluded the input state: macro transitions and logD lose
        # their anchor, pI breaks — refuse rather than return a truncated
        # result silently
        raise ValueError(
            f"电荷窗 [{min_charge:+d}, {max_charge:+d}] 不包含输入电荷 "
            f"z0={Chem.GetFormalCharge(mol):+d}:宏观 pKa/logD/pI 将失去参考态,"
            "请把窗口设为覆盖输入电荷")

    if route == "xtb":
        qm = QuantumEnergyModel(n_confs=n_confs)
        energies = qm.compute(mol, states, pka_window=pka_window, progress=progress)
        energy_backend = qm.name
        # only the sites that passed the screening gate (trusted + window) —
        # the raw table contains calibration artifacts the ensemble rejects,
        # which would make site_legend/n_active_sites contradict the pKa list
        kept_sites = getattr(qm, "last_kept_sites", qm.last_site_pka)
        # full-precision site pKas; the micropka stage reuses these instead of
        # recomputing every single-flip job (same protocol, half the cost).
        # After a tautomer switch the indices refer to the tautomer, not the
        # caller's SMILES — drop them so micropka recomputes consistently.
        site_pka_full = None if tautomer_switched else {
            f"{idx},{delta}": round(p, 4) for (idx, delta), p in kept_sites.items()}
        site_legend = [
            {"atom_index": idx, "kind": "acid" if d == -1 else "base",
             "pKa": round(p, 2)}
            for (idx, d), p in kept_sites.items()
        ]
    else:
        site_pka_full = None
        model = EnergyModel(backend="auto" if backend == "auto" else "fragment",
                            eps_eff=eps_eff, coupling=coupling)
        model.prepare(mol)
        energies = [model.energy(ms) for ms in states]
        energy_backend = model.oracle.name
        site_legend = [
            {"atom_index": s.atom_index, "kind": s.kind, "template": s.template,
             "pKa": round(s.pka, 2), "source": s.source}
            for s in model.sites
        ]

    charges = [ms.charge for ms in states]
    logps = [_logp(ms.smiles) for ms in states]
    for ms, g in zip(states, energies):
        ms.energy = g

    return {
        "object_data": {
            "initial_smiles": canonical(mol),
            # Order-preserving serialization for consumers that re-parse and
            # highlight atoms by index. MolToSmiles(canonical=False) is NOT
            # order-preserving once aromatic rings are rewritten (PubChem
            # Kekulé inputs silently permuted site highlights); a molblock is.
            "input_molblock": Chem.MolToMolBlock(mol, kekulize=False),
            "input_smiles_atomorder": Chem.MolToSmiles(mol, canonical=False),
            "input_query": canonical(mol),
            "min_pH": min_pH,
            "max_pH": max_pH,
            "min_charge": min_charge,
            "max_charge": max_charge,
            "charge_window_auto": True,
            "n_active_sites": len(site_legend),
            "method": "replica (open implementation of Rowan macropka)",
            "energy_backend": energy_backend,
        },
        "mol_properties": mol_properties(canonical(mol)),
        "microstates": [ms.to_dict() for ms in states],
        "site_legend": site_legend,
        "site_pka_full": site_pka_full,
        # keep transitions inside the physical pH axis (0-14): values outside
        # come either from placeholder charge levels (|pKa| ~ 40) or from
        # out-of-axis micro-pKa combos whose coupling accuracy is limited
        "pKa_values": [p for p in thermo.macro_pka_values(energies, charges)
                       if 0.0 <= p["pKa"] <= 14.0],
        "isoelectric_point": thermo.isoelectric_point(energies, charges, min_pH, max_pH),
        # logP of the INPUT state (identical to mol_properties.logP_crippen).
        # The old max-over-microstates heuristic broke for bases: Crippen
        # scores protonated cations higher than the neutral form, so the
        # "neutral logP" reference line sat above the whole logD curve.
        "logP_neutral": round(_logp(canonical(mol)), 3),
        "microstate_weights_by_pH": thermo.population_profile(
            energies, charges, min_pH, max_pH, 0.1
        ),
        "logD_by_pH": thermo.logd_profile(energies, charges, logps, min_pH, max_pH, 0.1),
        "warnings": warnings_out,
        "average_charge_by_pH": [
            [round(min_pH + 0.1 * i, 1),
             round(thermo.average_charge(energies, charges, round(min_pH + 0.1 * i, 1)), 4)]
            for i in range(int(round((max_pH - min_pH) / 0.1)))
        ],
    }


def run_micropka(
    smiles: str,
    pka_range: tuple[float, float] = (2.0, 12.0),
    deprotonate_elements: tuple[int, ...] = (7, 8, 16),
    protonate_elements: tuple[int, ...] = (7,),
    backend: str = "auto",
    progress=None,
    max_workers: int | None = None,
    site_pka: dict[tuple[int, int], float] | None = None,
) -> dict:
    """Replica of Rowan's (microscopic) pKa workflow: per-site pKa values.

    xtb route: single-flip states only, computed in parallel with the same
    geometry-reuse protocol; sites are screened against pka_range ± 5
    (stjames `reasonableness_buffer`).

    `site_pka` ({(atom_index, delta): pKa}) lets the caller reuse the single-
    flip energies a macropka run already computed for the same molecule —
    the server does this so the web UI pays for the quantum jobs once.
    """
    mol = _mol(smiles)
    route = _resolve_backend(backend)
    buffer = 5.0
    warnings_out: list[str] = []
    mol = _ensure_connected(mol, warnings_out)

    if site_pka is not None and site_pka:
        # wire format uses "idx,delta" string keys (JSON-safe); normalize to
        # the (atom_index, delta) tuples the lookups below expect
        norm: dict[tuple[int, int], float] = {}
        for k, v in site_pka.items():
            if isinstance(k, str):
                a, d = k.rsplit(",", 1)
                norm[(int(a), int(d))] = v
            else:
                norm[(int(k[0]), int(k[1]))] = v
        site_pka = norm

    conjugate_bases: list[dict] = []
    conjugate_acids: list[dict] = []

    # row-level display filters, shared by both backends:
    #   1. pH axis is 0-14 — values outside it are not physically meaningful
    #      to show (reference rows like pKa -2 get dropped);
    #   2. weakly-acidic N-H (amide/aniline/indole class) below pKa 10 is a
    #      calibration artifact, chemically those sites are >= 10.
    def _row_ok(name: str, pka: float) -> bool:
        if not (0.0 <= pka <= 14.0):
            return False
        if name in NON_GENUINE_ACID_N and pka < 10.0:
            return False
        if name == "alcohol_OH" and pka < ALCOHOL_FLOOR:
            return False
        return True

    if route == "xtb":
        from .energy import QuantumEnergyModel

        if not site_pka:
            # window must at least contain the parent and its single-flip
            # states (z0±1); [-2,+2] silently dropped sites for |z0|>=2 inputs
            z0 = Chem.GetFormalCharge(mol)
            states = enumerate_microstates(
                mol, min_charge=min(-2, z0 - 1), max_charge=max(2, z0 + 1),
                max_states=64)
            qm = QuantumEnergyModel(n_confs=2, max_workers=max_workers)
            # site-level pKa only needs single-flip states; skip the (expensive)
            # multi-flip coupling jobs entirely
            qm.compute(mol, states, pka_window=tuple(pka_range), progress=progress,
                       with_multi=False)
            site_pka = getattr(qm, "last_site_pka", {})

        for idx, name in acid_sites(mol):
            if mol.GetAtomWithIdx(idx).GetAtomicNum() not in deprotonate_elements:
                continue
            pka = site_pka.get((idx, -1))
            if pka is None or not (pka_range[0] - buffer <= pka <= pka_range[1] + buffer) or not _row_ok(name, pka):
                continue
            prod = _edit_charge(mol, idx, -1)
            conjugate_bases.append({
                "atom_index": idx, "template": name,
                "pKa": round(pka, 3),
                "conjugate": canonical(prod) if prod is not None else None,
            })
        for idx, name in base_sites(mol):
            if mol.GetAtomWithIdx(idx).GetAtomicNum() not in protonate_elements:
                continue
            pka = site_pka.get((idx, +1))
            if pka is None or not (pka_range[0] - buffer <= pka <= pka_range[1] + buffer) or not _row_ok(name, pka):
                continue
            prod = _edit_charge(mol, idx, +1)
            conjugate_acids.append({
                "atom_index": idx, "template": name,
                "pKa": round(pka, 3),
                "conjugate": canonical(prod) if prod is not None else None,
            })
        method = "replica-xtb (GFN2-xTB/ALPB water, Wagen-2024-style thermodynamic cycle)"
    else:
        model = EnergyModel(backend="auto" if backend == "auto" else "fragment")
        model.prepare(mol)

        for idx, name in acid_sites(mol):
            if mol.GetAtomWithIdx(idx).GetAtomicNum() not in deprotonate_elements:
                continue
            pka = model.acid_pka.get(idx, model.acid_default.get(idx))
            if pka is None or not (pka_range[0] - buffer <= pka <= pka_range[1] + buffer) or not _row_ok(name, pka):
                continue
            prod = _edit_charge(mol, idx, -1)
            conjugate_bases.append({
                "atom_index": idx, "template": name,
                "pKa": round(pka, 3),
                "conjugate": canonical(prod) if prod is not None else None,
            })
        for idx, name in base_sites(mol):
            if mol.GetAtomWithIdx(idx).GetAtomicNum() not in protonate_elements:
                continue
            pka = model.base_pka.get(idx, model.base_default.get(idx))
            if pka is None or not (pka_range[0] - buffer <= pka <= pka_range[1] + buffer) or not _row_ok(name, pka):
                continue
            prod = _edit_charge(mol, idx, +1)
            conjugate_acids.append({
                "atom_index": idx, "template": name,
                "pKa": round(pka, 3),
                "conjugate": canonical(prod) if prod is not None else None,
            })
        method = f"replica-{model.oracle.name}"

    acids = [b["pKa"] for b in conjugate_bases]
    bases = [a["pKa"] for a in conjugate_acids]
    return {
        "object_data": {
            "initial_smiles": canonical(mol),
            "input_molblock": Chem.MolToMolBlock(mol, kekulize=False),
            "input_smiles_atomorder": Chem.MolToSmiles(mol, canonical=False),
            "microscopic_pka_method": method,
            "warnings": warnings_out,
            "solvent": "water",
            "pka_range": list(pka_range),
            "deprotonate_elements": list(deprotonate_elements),
            "protonate_elements": list(protonate_elements),
        },
        "conjugate_bases": conjugate_bases,
        "conjugate_acids": conjugate_acids,
        "strongest_acid": round(min(acids), 3) if acids else None,
        "strongest_base": round(max(bases), 3) if bases else None,
    }
