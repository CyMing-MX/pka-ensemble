"""Microstate enumeration via SMARTS protonation rules.

Replicates the enumeration stage of Rowan's macropka workflow (which itself
follows Uni-pKa's enumerator design):

  * two template pools: acidic sites (can lose H+) and basic sites (can gain H+)
  * alternating deprotonation/protonation rounds (A-pool -> B-pool -> A-pool ...)
  * a "chemically reasonable" filter: RDKit sanitization rejects impossible
    valences / charge states, exactly like Uni-pKa's structure filter
  * states constrained to the user charge window [min_charge, max_charge]
  * deduplication by canonical (isomeric) SMILES

Every microstate carries `flips`: {atom_index: charge_delta} relative to the
input molecule, accumulated along the enumeration path. Atom indices are
stable across rounds because edits never reorder atoms. Rowan additionally
uses an AIMNet2-guided beam search and explicit tautomer enumeration;
proton-transfer tautomers appear here for free whenever an anion can be
re-protonated at a different site.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------------------
# SMARTS templates.
# Acidic templates match substructures containing an H-bearing atom that can
# lose H+; basic templates' first atom is the site that gains H+.
# Order matters: first match wins per atom (most specific first).
# ---------------------------------------------------------------------------

ACID_SMARTS: list[tuple[str, str]] = [
    # ammonium N-H (zwitterion inputs like [NH3+]...CC(=O)[O-] deprotonate
    # to the neutral amine — same ensemble as the neutral input form)
    ("ammonium_NH", "[NX4;+1;H2,H3]"),
    ("carboxylic_acid", "[CX3](=O)[OX2H1]"),
    # NOTE: these two used the $([SX4](=O)(=O)[OX2H1]) single-atom recursive
    # form, whose match contains ONLY the sulfur — the H-bearing-atom lookup
    # then found no H and silently dropped every sulfonic/sulfonamide site
    # (PFAS run caught it). Match through the O/N instead.
    ("sulfonic_acid", "[$([OX2H1][SX4](=O)(=O)),$([OX2H1][SX3]=O)]"),
    ("phosphonic_acid", "[PX4](=O)([OX2H1])[OX2H1]"),
    ("phosphoric_acid", "[PX4](=O)([OX2H1])[OX2H0]"),
    ("boronic_acid", "[BX3]([OX2H1])[OX2H1]"),
    ("sulfonamide_NH",
     "[$([NX3H1][SX4](=O)(=O)),$([NX3H2][SX4](=O)(=O)),$([NX4H2+][SX4](=O)(=O))]"),
    ("tetrazole_NH", "[nX3H1]1[nX2][nX2][nX2][nX2]1"),
    ("imide_NH", "[NX3H1]([CX3]=O)[CX3]=O"),
    ("hydroxamic_acid", "[CX3](=[OX1])[NX3][OX2H1]"),
    ("phenol_enol_oxime",
     "[$([OX2H1][c]),$([OX2H1][CX3]=[CX3]),$([OX2H1][NX2])]"),
    ("aromatic_NH", "[nX3H1]"),
    ("amide_NH", "[NX3H1][CX3]=[OX1]"),
    ("anilide_like_NH", "[NX3H1]"),
    ("thiol_sh", "[SX2H1]"),
    ("thiourea_NH", "[NX3H1][CX3]=[SX1]"),
    ("hydrazine_NH", "[NX3H1][NX3]"),
    ("alcohol_OH", "[OX2H1][#6,#1]"),
]

BASE_SMARTS: list[tuple[str, str]] = [
    ("amidine_guanidine_N", "[$([NX3;H0,H1][CX3]=[NX2]),$([NX3;H0,H1][CX3]([NX3])=[NX2])]"),
    ("aromatic_N_acceptor", "[nX2;r5,r6;!$([nH])]"),
    ("aniline_N", "[NX3;H0,H1,H2;$(Nc)]"),
    ("aliphatic_amine_N", "[NX3;H0,H1,H2;!$(N[C,S,P]=O);!$(N[a]);!$(N=*);!$([N+])]"),
    ("nitrile_N", "[NX1]#[CX2]"),  # off by default
    # protonating an anionic O (carboxylate/sulfonate/phosphate) gives the
    # conjugate acid — this is what makes zwitterion / salt inputs (e.g.
    # "[NH3+]CC(=O)[O-]", "CC(=O)[O-]") yield the same ensemble as their
    # neutral forms; neutral molecules have no -1 O so it stays inert there.
    # (NOTE: deliberately NO amide_N_protonation template — amide pKaH ~ -1
    # is never populated in pH 0-14, and routing it through N-base-ali once
    # let protonated-amide microstates dominate base ensembles, the MTX bug.)
    ("anion_O_protonation", "[O;-1;H0]"),
]

# Deprotonation must start from an H-bearing atom. N with +1 charge is
# allowed: a zwitterion input like [NH3+]CC(=O)[O-] legitimately
# deprotonates to the neutral amine (valence sanity is enforced downstream
# by _edit_charge's sanitization). Higher positive charges stay excluded.

# N-H acid templates that are NOT genuinely acidic in water: amide / aniline /
# indole-like N-H deprotonates at pKa >= ~10 by chemistry, so a calibrated
# value below 10 is a calibration artifact, not a real site. Genuinely acidic
# N-H templates (tetrazole ~4.9, sulfonamide ~10, imide ~9.9) are excluded.
NON_GENUINE_ACID_N = {"aromatic_NH", "amide_NH", "anilide_like_NH",
                      "hydrazine_NH", "thiourea_NH"}

# Aliphatic alcohol O-H is pKa >= ~13 in water (methanol 15.5); this
# calibration maps some sterically odd alcohols far lower (vinorelbine 2.2),
# which used to become the "most acidic site". Phenols/enols (7-10) are
# genuinely acidic and stay untouched.
ALCOHOL_FLOOR = 12.0


@dataclass
class Microstate:
    """One enumerated protonation microstate."""

    smiles: str                                 # canonical SMILES
    charge: int                                 # net formal charge
    flips: dict[int, int] = field(default_factory=dict)  # atom_idx -> charge delta vs input
    energy: float | None = None                 # dimensionless beta*G (set later)
    tag: str = "input"

    def to_dict(self) -> dict:
        return {
            "smiles": self.smiles,
            "charge": self.charge,
            "energy": None if self.energy is None else round(self.energy, 3),
            "tag": self.tag,
        }


def _matches(mol: Chem.Mol, smarts: str) -> list[tuple[int, ...]]:
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        return []
    return mol.GetSubstructMatches(patt, uniquify=False)


def _h_bearing_atom(mol: Chem.Mol, match: tuple[int, ...]) -> int | None:
    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetTotalNumHs() > 0 and atom.GetFormalCharge() <= 1:
            return idx
    return None


def _edit_charge(mol: Chem.Mol, idx: int, delta: int) -> Chem.Mol | None:
    """Add/remove one H and shift the formal charge by delta (+1/-1).

    Aromatic acceptors (fused azines, imidazoles...) sometimes fail full
    sanitization after explicit-H bookkeeping; a second attempt lets RDKit
    infer H counts implicitly before giving up.
    """
    def attempt(set_explicit_h: bool) -> Chem.Mol | None:
        rw = Chem.RWMol(mol)
        atom = rw.GetAtomWithIdx(idx)
        nh = atom.GetTotalNumHs() + delta
        if nh < 0:
            return None
        q = atom.GetFormalCharge() + delta
        if q > 2 or q < -2:
            return None
        atom.SetFormalCharge(q)
        if set_explicit_h:
            atom.SetNumExplicitHs(nh)
            atom.SetNoImplicit(True)
        out = rw.GetMol()
        try:
            Chem.SanitizeMol(out)
        except Exception:
            return None
        if not set_explicit_h and out.GetAtomWithIdx(idx).GetTotalNumHs() != nh:
            # implicit-H re-inference disagreed with the intended H count:
            # the state would no longer differ from its parent by exactly
            # delta H, corrupting every downstream flip-based energy
            return None
        return out

    out = attempt(set_explicit_h=True)
    if out is not None:
        return out
    return attempt(set_explicit_h=False)


def acid_sites(mol: Chem.Mol) -> list[tuple[int, str]]:
    """All deprotonatable atoms: (atom_index, template_name)."""
    out: list[tuple[int, str]] = []
    taken: set[int] = set()
    for name, smarts in ACID_SMARTS:
        for match in _matches(mol, smarts):
            idx = _h_bearing_atom(mol, match)
            if idx is not None and idx not in taken:
                taken.add(idx)
                out.append((idx, name))
    return out


def base_sites(mol: Chem.Mol, include_nitrile: bool = False) -> list[tuple[int, str]]:
    """All protonatable atoms: (atom_index, template_name)."""
    out: list[tuple[int, str]] = []
    taken: set[int] = set()
    for name, smarts in BASE_SMARTS:
        for match in _matches(mol, smarts):
            idx = match[0]
            if idx not in taken:
                taken.add(idx)
                out.append((idx, name))
    if not include_nitrile:
        out = [(i, n) for i, n in out if n != "nitrile_N"]
    return out


def net_charge(mol: Chem.Mol) -> int:
    return Chem.GetFormalCharge(mol)


def canonical(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol)


def enumerate_microstates(
    smiles_or_mol: str | Chem.Mol,
    min_charge: int = -2,
    max_charge: int = 2,
    max_states: int = 128,
    max_rounds: int = 12,
) -> list[Microstate]:
    """Enumerate the protonation ensemble within the charge window.

    Alternating rounds of deprotonation/protonation from the input state,
    deduplicated by canonical SMILES, mirroring Uni-pKa's A/B micropool logic.
    """
    mol = Chem.MolFromSmiles(smiles_or_mol) if isinstance(smiles_or_mol, str) else smiles_or_mol
    if mol is None:
        raise ValueError("invalid molecule")

    start = Microstate(smiles=canonical(mol), charge=net_charge(mol))
    by_smiles: dict[str, Microstate] = {start.smiles: start}
    frontier: list[tuple[Chem.Mol, dict[int, int]]] = [(mol, {})]

    for _round in range(max_rounds):
        if not frontier or len(by_smiles) >= max_states:
            break
        produced: list[tuple[Chem.Mol, dict[int, int]]] = []

        # -- deprotonation round (A2B) --------------------------------------
        for m, flips in frontier:
            if len(by_smiles) >= max_states:
                break
            z = net_charge(m)
            for idx, name in acid_sites(m):
                if z - 1 < min_charge:
                    break
                prod = _edit_charge(m, idx, -1)
                if prod is None:
                    continue
                s = canonical(prod)
                if s in by_smiles:
                    continue
                nf = dict(flips)
                nf[idx] = nf.get(idx, 0) - 1
                ms = Microstate(smiles=s, charge=net_charge(prod), flips=nf,
                                tag=f"acid:{name}@{idx}")
                by_smiles[s] = ms
                produced.append((prod, nf))

        # -- protonation round (B2A; also on fresh anions -> zwitterions) ---
        for m, flips in frontier + list(produced):
            if len(by_smiles) >= max_states:
                break
            z = net_charge(m)
            for idx, name in base_sites(m):
                if z + 1 > max_charge:
                    break
                prod = _edit_charge(m, idx, +1)
                if prod is None:
                    continue
                s = canonical(prod)
                if s in by_smiles:
                    continue
                nf = dict(flips)
                nf[idx] = nf.get(idx, 0) + 1
                ms = Microstate(smiles=s, charge=net_charge(prod), flips=nf,
                                tag=f"base:{name}@{idx}")
                by_smiles[s] = ms
                produced.append((prod, nf))

        frontier = produced

    states = [ms for ms in by_smiles.values() if min_charge <= ms.charge <= max_charge]
    states.sort(key=lambda ms: (ms.charge, ms.smiles))
    return states
