"""Benchmark pka_replica against the pKahub monoprotic dataset
(Sipos-Szabó et al., JCIM 2026 — pkahub.ttk.hu).

Monoprotic protocol (paper): molecules with a single experimental pKa and
an assigned charge-state transition 0→-1 (acid) or 1→0 (base).
Matching: closest predicted site in the matching direction
("best-match"); the strongest-acid/base convention is also reported.
"""

import csv
import glob
import sys
import time
import warnings
from collections import Counter, defaultdict

warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\zcode_rea\pka_replica")

from rdkit import Chem  # noqa: E402

from pka_replica.workflow import run_micropka  # noqa: E402


def load_pkahub(paths):
    entries = []
    seen_transitions = {}
    for path in paths:
        for row in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
            try:
                pka = float(row["pka_value"])
            except (KeyError, ValueError):
                continue
            if not (0.0 <= pka <= 14.0):
                continue
            pre, post = row["charge_state_pre"], row["charge_state_post"]
            molid = row["molid"]
            smiles = row["smiles"]
            if not smiles or Chem.MolFromSmiles(smiles) is None:
                continue
            key = (molid, pre, post)
            # average duplicates of the same transition (paper protocol)
            if key in seen_transitions:
                seen_transitions[key]["sum"] += pka
                seen_transitions[key]["n"] += 1
                continue
            ent = {"molid": molid, "pre": int(pre), "post": int(post),
                   "smiles": smiles, "source": row["data_source"],
                   "sum": pka, "n": 1}
            seen_transitions[key] = ent
    for ent in seen_transitions.values():
        ent["pka"] = ent["sum"] / ent["n"]
        del ent["sum"]
        entries.append(ent)
    return entries


def group_monoprotic(entries):
    """Group averaged transitions per molecule; keep molecules with exactly
    one transition of type 0→-1 (acid) or 1→0 (base) and no others."""
    per_mol = defaultdict(list)
    for e in entries:
        per_mol[e["molid"]].append(e)
    out = []
    for molid, es in per_mol.items():
        if len(es) != 1:
            continue
        e = es[0]
        if (e["pre"], e["post"]) == (0, -1):
            e["dir"] = "acid"
        elif (e["pre"], e["post"]) == (1, 0):
            e["dir"] = "base"
        else:
            continue
        e["molid"] = molid
        out.append(e)
    return out


def main():
    backend = sys.argv[1] if len(sys.argv) > 1 else "sites"
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    MAXW = int(sys.argv[3]) if len(sys.argv) > 3 else None

    entries = load_pkahub(glob.glob(r"D:\zcode_rea\pka_replica\pkahub\*.tsv"))
    mono = group_monoprotic(entries)
    # dedupe identical molecules (same canonical SMILES across sources)
    seen = {}
    for e in mono:
        m = Chem.MolFromSmiles(e["smiles"])
        can = Chem.MolToSmiles(m)
        if can not in seen:
            seen[can] = e
    mono = list(seen.values())
    if limit:
        mono = mono[:limit]
    print(f"monoprotic molecules: {len(mono)} "
          f"({Counter(e['dir'] for e in mono)})", flush=True)

    errs_best = []
    errs_strong = []
    t0 = time.time()
    skipped = 0
    for n, e in enumerate(mono):
        try:
            micro = run_micropka(e["smiles"], backend=backend,
                                 max_workers=MAXW)
        except Exception:
            skipped += 1
            continue
        if e["dir"] == "acid":
            preds = [b["pKa"] for b in micro["conjugate_bases"]]
            strongest = micro["strongest_acid"]
        else:
            preds = [a["pKa"] for a in micro["conjugate_acids"]]
            strongest = micro["strongest_base"]
        preds = [p for p in preds if -2 < p < 16]
        if not preds or strongest is None:
            skipped += 1
            continue
        errs_best.append(min(abs(p - e["pka"]) for p in preds))
        errs_strong.append(abs(strongest - e["pka"]))
        print(f"    m{n}: exp={e['pka']:.2f} "
              f"pred_best={min(preds, key=lambda p: abs(p - e['pka'])):.2f}",
              flush=True)
        if n % 500 == 0:
            print(f"  {n}/{len(mono)} bestMAE={sum(errs_best)/len(errs_best):.2f}",
                  flush=True)

    dt = time.time() - t0
    print(f"\n=== pKahub monoprotic ({backend}): {len(errs_best)} evaluated, "
          f"{skipped} skipped, {dt:.0f}s ===")
    print(f"  best-match MAE : {sum(errs_best)/len(errs_best):.2f}")
    print(f"  strongest  MAE : {sum(errs_strong)/len(errs_strong):.2f}")
    # error distribution like paper Figure 6
    import collections
    bins = collections.Counter()
    for e_ in errs_best:
        bins[int(min(e_, 4))] += 1
    print("  error bins (best):", dict(sorted(bins.items())))


if __name__ == "__main__":
    main()
