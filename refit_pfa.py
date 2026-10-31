"""Refit a dedicated calibration class for (partially) fluorinated carboxylic
acids. The generic O-acid fit (normal acids) extrapolates through the -2
clamp on perfluoro acids (PFAS run: all PFCAs clamped at -2 vs chem ~0.3).

Ladder references anchor the alpha-F effect at increasing fluorination.
Protocol matches production: AIMNet2 preopt -> GFN2 gas opt -> ALPB-water sp.
raw = (G_conjugate_base - G_parent) / 1.3644 in pKa units (deprotonation).
"""

import sys
import warnings

warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\zcode_rea\pka_replica")

import torch  # noqa: E402

torch.set_num_threads(2)

from rdkit import Chem  # noqa: E402

from pka_replica.energy import (_apply_flips_to_geometry,  # noqa: E402
                                _parent_geometry_from_mol, _xtb_job)
from pka_replica.enumerator import _edit_charge, acid_sites

# name, smiles, experimental pKa (25 C, water)
REFS = [
    ("fluoroacetic acid", "O=C(O)CF", 2.59),
    ("difluoroacetic acid", "O=C(O)C(F)F", 1.33),
    ("trifluoroacetic acid", "O=C(O)C(F)(F)F", 0.23),
    ("pentafluoropropionic acid", "O=C(O)C(F)(F)C(F)(F)F", 0.59),
]


def main():
    from concurrent.futures import ThreadPoolExecutor

    results = []

    def job(item):
        name, smi, exp = item
        mol = Chem.MolFromSmiles(smi)
        geo = _parent_geometry_from_mol(mol, n_embed=2)
        if geo is None:
            return None
        symbols, coords = geo
        molH = Chem.AddHs(Chem.Mol(mol))
        sites = [(i, n) for i, n in acid_sites(mol) if n == "carboxylic_acid"]
        if not sites:
            return None
        idx = sites[0][0]
        g_parent = _xtb_job({"symbols": symbols, "coords": coords,
                             "charge": 0, "timeout": 180})
        prod = _edit_charge(Chem.Mol(mol), idx, -1)
        g = _apply_flips_to_geometry(symbols, coords, molH, {idx: -1})
        g_conj = _xtb_job({"symbols": g[0], "coords": g[1],
                           "charge": Chem.GetFormalCharge(prod),
                           "timeout": 180})
        if g_parent[0] is None or g_conj[0] is None:
            return None
        raw = (g_conj[0] - g_parent[0]) * 627.509474 / 1.3644
        return (name, exp, raw)

    with ThreadPoolExecutor(max_workers=2) as pool:
        for r in pool.map(job, REFS):
            if r is None:
                print("SKIP/FAIL", flush=True)
                continue
            results.append(r)
            print(f"{r[0]:28s} exp={r[1]:6.2f} raw={r[2]:8.2f}", flush=True)

    n = len(results)
    if n < 3:
        print("too few points, keeping generic O-acid fit")
        return
    sx = sum(r[2] for r in results)
    sy = sum(r[1] for r in results)
    sxx = sum(r[2] ** 2 for r in results)
    sxy = sum(r[2] * r[1] for r in results)
    a = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    b = (sy - a * sx) / n
    resid = [abs(e - (a * raw + b)) for _, e, raw in results]
    print(f"\nn={n}  slope={a:.4f}  intercept={b:.3f}")
    print(f"residuals: mean={sum(resid)/n:.2f}  max={max(resid):.2f}")
    print(f'new CAL_TABLE["O-acid-fluoro"] = ({a:.6f}, {b:.3f})')
    # sanity: predictions for the references themselves
    for name, e, raw in results:
        print(f"  {name:28s} exp={e:5.2f} fit={a*raw+b:6.2f}")


if __name__ == "__main__":
    main()
