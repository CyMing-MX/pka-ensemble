"""Benchmark pka_replica against the QupKake/Uni-pKa benchmark SDFs
(Novartis 280 + literature 122 molecules, each with one experimental
site-specific pKa and its type basic/acidic).

Two fair comparison metrics for site-level datasets:
  best-match  — min |pred_site - exp| over our sites in the same direction
                (site *assignment* is not charged to the predictor)
  strongest   — the strongest-acid/base convention (Rowan display)
"""

import sys
import time
import warnings

warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\zcode_rea\pka_replica")

from rdkit import Chem  # noqa: E402

from pka_replica.workflow import run_micropka  # noqa: E402


def load_sdf(path):
    out = []
    for m in Chem.SDMolSupplier(path):
        if m is None:
            continue
        p = m.GetPropsAsDict()
        out.append({
            "smiles": p.get("smiles") or Chem.MolToSmiles(m),
            "exp": float(p["pka"]),
            "type": p.get("pka_type", "?"),
        })
    return out


def main():
    backend = sys.argv[1] if len(sys.argv) > 1 else "sites"
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    rows = []
    for path, label in [("novartis.sdf", "novartis"), ("literature.sdf", "literature")]:
        for r in load_sdf(path):
            r["set"] = label
            rows.append(r)
    if limit:
        rows = rows[:limit]
    print(f"benchmark rows: {len(rows)} (backend={backend})", flush=True)

    errs_best = []
    errs_strong = []
    per_set_best = {}
    per_set_strong = {}
    t0 = time.time()
    for n, r in enumerate(rows):
        try:
            micro = run_micropka(r["smiles"], backend=backend)
        except Exception:
            continue
        if r["type"] == "basic":
            site_preds = [a["pKa"] for a in micro["conjugate_acids"]]
            strongest = micro["strongest_base"]
        else:
            site_preds = [b["pKa"] for b in micro["conjugate_bases"]]
            strongest = micro["strongest_acid"]
        site_preds = [p for p in site_preds if -2 < p < 16]
        if not site_preds or strongest is None:
            continue

        err_best = min(abs(p - r["exp"]) for p in site_preds)
        err_strong = abs(strongest - r["exp"])
        errs_best.append(err_best)
        errs_strong.append(err_strong)
        per_set_best.setdefault(r["set"], []).append(err_best)
        per_set_strong.setdefault(r["set"], []).append(err_strong)
        if max(err_best, err_strong) > 3 and n % 10 == 0:
            print(f"    OUTLIER n={n} pred_best="
                  f"{min(site_preds, key=lambda p: abs(p - r['exp'])):.2f} "
                  f"pred_strong={strongest:.2f} exp={r['exp']:.2f} "
                  f"type={r['type']} {r['smiles'][:40]}", flush=True)
        if n % 25 == 0:
            print(f"  {n+1}/{len(rows)} done", flush=True)

    dt = time.time() - t0
    print(f"\n=== {backend} benchmark: {len(errs_best)} evaluated in {dt:.0f}s ===")
    for label in sorted(per_set_best):
        sb = per_set_best[label]
        ss = per_set_strong[label]
        print(f"  {label:10s} n={len(sb):3d}  bestMAE={sum(sb)/len(sb):.2f}  "
              f"strongMAE={sum(ss)/len(ss):.2f}")
    print(f"  {'TOTAL':10s} n={len(errs_best):3d}  bestMAE={sum(errs_best)/len(errs_best):.2f}  "
          f"strongMAE={sum(errs_strong)/len(errs_strong):.2f}")


if __name__ == "__main__":
    main()
