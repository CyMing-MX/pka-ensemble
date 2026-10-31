"""Validation suite v4 (parallel + screened xtb protocol)."""

import sys
import time

sys.path.insert(0, r"D:\zcode_rea\pka_replica")

from pka_replica.workflow import run_macropka  # noqa: E402

TESTS = [
    ("acetic acid", "CC(=O)O", [4.76]),
    ("phenol", "c1ccccc1O", [9.99]),
    ("pyridine", "c1ccncc1", [5.23]),
    ("imidazole", "c1ncc[nH]1", [6.95]),
    ("aniline", "c1ccccc1N", [4.63]),
    ("glycine", "NCC(=O)O", [2.34, 9.60]),
    ("piperazine", "C1CNCCN1", [9.82, 5.35]),
    ("citric acid", "OC(=O)CC(O)(CC(=O)O)C(=O)O", [3.13, 4.76, 6.40]),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", [4.43]),
    ("4-aminopyridine", "Nc1ccncc1", [9.17, 3.6]),
]

backend = sys.argv[1] if len(sys.argv) > 1 else "xtb"
errs = []
t_all = time.time()
for name, smi, exp in TESTS:
    t0 = time.time()
    try:
        r = run_macropka(smi, backend=backend)
        dt = time.time() - t0
        pkas = sorted(
            [p["pKa"] for p in r["pKa_values"] if -1.5 < p["pKa"] < 15.5],
            reverse=True,
        )
        pred = [round(p, 2) for p in pkas]
        for e, p in zip(exp, pred):
            errs.append(abs(e - p))
        print(f"{name:15s} exp={exp}  got={pred}  pI={r['isoelectric_point']} "
              f"({dt:.1f}s, win={r['object_data']['min_charge']}..{r['object_data']['max_charge']})")
    except Exception as ex:
        import traceback; traceback.print_exc()
        print(f"{name:15s} ERROR {ex}")

if errs:
    print(f"\nMAE: {sum(errs)/len(errs):.2f} pKa over {len(errs)} pairs; total {time.time()-t_all:.0f}s")

# --- speed test: methotrexate ---
print("\n=== methotrexate (46 heavy atoms, ~9 ionizable sites) ===")
MTX = "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1"
t0 = time.time()
r = run_macropka(MTX, backend="xtb")
dt = time.time() - t0
pkas = sorted([p["pKa"] for p in r["pKa_values"]], reverse=True)
print(f"time: {dt:.1f}s   states: {len(r['microstates'])}   "
      f"window: {r['object_data']['min_charge']}..{r['object_data']['max_charge']}")
print("site pKas:", sorted([s['pKa'] for s in r['site_legend']], reverse=True))
print("macro pKas:", [round(p, 2) for p in pkas if -2 < p < 16])
print("pI:", r["isoelectric_point"])
