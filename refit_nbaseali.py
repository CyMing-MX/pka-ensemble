"""Refit the N-base-ali (aliphatic amine) calibration with a diverse
reference set. The old 2-point fit (methylamine/piperidine) produced a
physically wrong negative slope that extrapolated badly (+3.3 bias).

Protocol matches production: AIMNet2 preopt -> GFN2 gas opt -> ALPB-water sp.
raw = (G_conjugate_acid - G_base) / 1.3644 in pKa units.
"""

import sys
import warnings

warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\zcode_rea\pka_replica")

import torch  # noqa: E402

torch.set_num_threads(2)

from rdkit import Chem  # noqa: E402

from pka_replica.energy import _apply_flips_to_geometry, _parent_geometry_from_mol, _xtb_job  # noqa: E402
from pka_replica.enumerator import _edit_charge, base_sites  # noqa: E402

# name, smiles, experimental pKa of the conjugate acid (25 C, water)
REFS = [
    ("ammonia", "N", 9.25),
    ("methylamine", "CN", 10.64),
    ("ethylamine", "CCN", 10.63),
    ("propylamine", "CCCN", 10.53),
    ("dimethylamine", "CNC", 10.73),
    ("trimethylamine", "CN(C)C", 9.80),
    ("diethylamine", "CCNCC", 10.98),
    ("triethylamine", "CCN(CC)CC", 10.75),
    ("morpholine", "O1CCNCC1", 8.36),
    ("piperidine", "C1CCCCN1", 11.22),
    ("pyrrolidine", "C1CCCN1", 11.27),
    ("azetidine", "C1CNC1", 11.29),
    ("cyclohexylamine", "NC1CCCCC1", 10.64),
    ("benzylamine", "NCc1ccccc1", 9.34),
    ("ethanolamine", "NCCO", 9.50),
    ("trifluoroethylamine", "NCC(F)(F)F", 5.80),
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
        sites = [(i, n) for i, n in base_sites(mol)
                 if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]
        if not sites:
            return None
        idx = sites[0][0]
        g_parent = _xtb_job({"symbols": symbols, "coords": coords,
                             "charge": 0, "timeout": 180})
        prod = _edit_charge(Chem.Mol(mol), idx, +1)
        g = _apply_flips_to_geometry(symbols, coords, molH, {idx: +1})
        g_conj = _xtb_job({"symbols": g[0], "coords": g[1],
                           "charge": Chem.GetFormalCharge(prod),
                           "timeout": 180})
        if g_parent[0] is None or g_conj[0] is None:
            return None
        raw = (g_parent[0] - g_conj[0]) * 627.509474 / 1.3644
        return (name, exp, raw)

    with ThreadPoolExecutor(max_workers=2) as pool:
        for r in pool.map(job, REFS):
            if r is None:
                print("SKIP/FAIL", flush=True)
                continue
            results.append(r)
            print(f"{r[0]:22s} exp={r[1]:6.2f} raw={r[2]:8.2f}", flush=True)

    n = len(results)
    if n < 3:
        print("too few points, keeping old fit")
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
    print(f'new CAL_TABLE["N-base-ali"] = ({a:.6f}, {b:.3f})')


if __name__ == "__main__":
    main()
