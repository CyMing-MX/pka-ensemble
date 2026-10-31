"""Refit the GFN2/ALPB linear calibration table for the AIMNet2 geometry
protocol (AIMNet2 gas-phase opt + GFN2/ALPB-water single point)."""

import collections
import sys

sys.path.insert(0, r"D:\zcode_rea\pka_replica")

from rdkit import Chem  # noqa: E402

from pka_replica.energy import (_aimnet_calc, _apply_flips_to_geometry,  # noqa: E402
                                _parent_geometry_from_mol, _xtb_job)
from pka_replica.enumerator import _edit_charge, canonical  # noqa: E402

_aimnet_calc()  # warm the model before forking threads

REFS = [
    # name, smiles, exp pKa, class, element, direction
    ("acetic acid", "CC(=O)O", 4.76, "O-acid", 8, "acid"),
    ("phenol", "c1ccccc1O", 9.99, "O-acid", 8, "acid"),
    ("methanol", "CO", 15.54, "O-acid", 8, "acid"),
    ("tetrazole", "c1nnn[nH]1", 4.90, "N-acid", 7, "acid"),
    ("pyrrole", "c1cc[nH]c1", 16.5, "N-acid", 7, "acid"),
    ("2-pyridone", "O=c1cccc[nH]1", 11.7, "N-acid", 7, "acid"),
    ("ethanethiol", "CCS", 10.60, "S-acid", 16, "acid"),
    ("thiophenol", "c1ccccc1S", 6.62, "S-acid", 16, "acid"),
    ("pyridine", "c1ccncc1", 5.23, "N-base-aro", 7, "base"),
    ("imidazole", "c1ncc[nH]1", 6.95, "N-base-aro", 7, "base"),
    ("methylamine", "CN", 10.64, "N-base-ali", 7, "base"),
    ("piperidine", "C1CCCCN1", 11.20, "N-base-ali", 7, "base"),
    ("aniline", "c1ccccc1N", 4.63, "N-base-aniline", 7, "base"),
]

from concurrent.futures import ThreadPoolExecutor  # noqa: E402

out = []


def job(item):
    name, smi, exp, cls, elem, kind = item
    mol = Chem.MolFromSmiles(smi)
    geo = _parent_geometry_from_mol(mol, n_embed=2)
    if geo is None:
        return None
    symbols, coords = geo
    molH = Chem.AddHs(Chem.Mol(mol))

    parent = _xtb_job({"symbols": symbols, "coords": coords, "charge": 0,
                       "timeout": 180})
    # site: first matching heavy atom of the requested element
    from pka_replica.enumerator import acid_sites, base_sites

    sites = (acid_sites if kind == "acid" else base_sites)(mol)
    sites = [(i, n) for i, n in sites if mol.GetAtomWithIdx(i).GetAtomicNum() == elem]
    if not sites or parent[0] is None:
        return None
    idx = sites[0][0]
    prod = _edit_charge(Chem.Mol(mol), idx, -1 if kind == "acid" else +1)
    g = _apply_flips_to_geometry(symbols, coords, molH, {idx: -1 if kind == "acid" else +1})
    conj = _xtb_job({"symbols": g[0], "coords": g[1],
                     "charge": Chem.GetFormalCharge(prod), "timeout": 180})
    if conj[0] is None:
        return None
    # raw in the unified deprotonation direction (positive = easier to
    # deprotonate): acid flips take (conj - parent); base flips take
    # (parent - conj), i.e. the acidity of the conjugate acid.
    raw = (conj[0] - parent[0]) * 627.509474 / 1.3644
    if kind == "base":
        raw = -raw
    return (name, cls, exp, raw)


with ThreadPoolExecutor(max_workers=6) as pool:
    for r in pool.map(job, REFS):
        if r is None:
            print("SKIP/FAIL")
            continue
        name, cls, exp, raw = r
        out.append((name, cls, exp, raw))
        print(f"{name:14s} {cls:16s} exp={exp:6.2f} raw={raw:9.3f}", flush=True)

by = collections.defaultdict(list)
for _, cls, exp, raw in out:
    by[cls].append((raw, exp))
print()
print("CAL_TABLE = {")
for cls, pts in sorted(by.items()):
    n = len(pts)
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] ** 2 for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    if n >= 2:
        a = (n * sxy - sx * sy) / (n * sxx - sx * sx)
        b = (sy - a * sx) / n
    else:
        a, b = 1.0, sy - sx
    resid = [abs(e - (a * r + b)) for r, e in pts]
    print(f'    "{cls}": ({a:.6f}, {b:.3f}),  # n={n} maxresid={max(resid):.2f}')
print("}")
