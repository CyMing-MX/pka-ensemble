"""Run our pKa tool on every drug in ACID-Analytes.xlsx and compare with
Chemicalize's pKa-A / pKa-B (most acidic / most basic apparent pKa).

Checkpointed: every molecule is appended to acid_analytes_result.tsv as soon
as it finishes, so an interrupted run can be inspected/restarted safely.
"""
import sys
import time
import faulthandler

sys.stdout.reconfigure(encoding="utf-8")
faulthandler.enable()

from openpyxl import load_workbook
from pka_replica import pubchem
from pka_replica.workflow import run_macropka, run_micropka
from pka_replica.server import parse_smiles_batch  # noqa: F401 (not used, kept for parity)

TSV = r"D:\zcode_rea\pka_replica\acid_analytes_result.tsv"

wb = load_workbook(r"D:\zcode_rea\ACID-Analytes.xlsx", read_only=True)
ws = wb["drug"]
rows = list(ws.iter_rows(values_only=True))
hdr = rows[0]
drugs = []
def to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def cell_float(v):
    if v is None:
        return None
    s = str(v).strip()
    if not s or s.lower() in ("n/a", "null", "none"):
        return None
    return to_float(s)


for r in rows[1:]:
    name = str(r[0] or "").strip()
    if not name:
        continue
    cas_field = str(r[1] or "")
    drugs.append({
        "name": name,
        "cas": [c.strip() for c in cas_field.replace("，", ",").split(",")
                if c.strip() and c.strip().lower() != "n/a"],
        "chem_A": cell_float(r[29]),
        "chem_B": cell_float(r[30]),
        "chem_pI": cell_float(r[31]),
    })
print(f"loaded {len(drugs)} drugs from xlsx", flush=True)


def to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def resolve_smiles(d):
    # name first, then each CAS
    for cand in [d["name"]] + d["cas"]:
        try:
            info = pubchem.resolve(cand)
            return info["smiles"], cand
        except Exception:
            continue
    return None, None


out = open(TSV, "a", encoding="utf-8")
if out.tell() == 0:
    out.write("name\tresolved_via\tsmiles\tchemA\tchemB\tourA\tourB\terrA\terrB\t"
              "chem_pI\tour_pI\terr_pI\tnote\n")
    out.flush()
else:
    out.seek(0)

# resume: names already successfully computed in a previous (killed) run
done_names = set()
retry_names = set()
for line in open(TSV, encoding="utf-8").read().splitlines()[1:]:
    f = line.split("\t")
    if len(f) < 13 or not f[0]:
        continue
    (our_a, our_b, note) = f[5], f[6], f[12]
    failed = ("计算失败" in note) or ("解析失败" in note)
    if failed:
        retry_names.add(f[0])
        done_names.discard(f[0])  # a stale failure must not shadow a success
    elif our_a or our_b or "无酸" in note or "无碱" in note:
        done_names.add(f[0])
        retry_names.discard(f[0])
print(f"resume: {len(done_names)} done, {len(retry_names)} to retry", flush=True)

done_A, done_B, done_pI = [], [], []
t_start = time.time()
for i, d in enumerate(drugs):
    if d["name"] in done_names and d["name"] not in retry_names:
        continue
    note = []
    smiles, via = resolve_smiles(d)
    if smiles is None:
        out.write(f"{d['name']}\t\t\t{d['chem_A']}\t{d['chem_B']}\t\t\t\t\t"
                  f"{d['chem_pI']}\t\t\tPubChem 解析失败\n")
        out.flush()
        print(f"[{i+1}/{len(drugs)}] {d['name']}: RESOLVE FAILED", flush=True)
        continue
    try:
        macro = run_macropka(smiles, backend="xtb")
        micro = run_micropka(smiles, backend="xtb",
                             site_pka=macro.get("site_pka_full"))
        our_a = micro.get("strongest_acid")
        our_b = micro.get("strongest_base")
        our_pi = macro.get("isoelectric_point")
        if macro.get("warnings"):
            note.append(";".join(str(w)[:40] for w in macro["warnings"]))
        if our_a is None:
            note.append("无酸性位点")
        if our_b is None:
            note.append("无碱性位点")
    except Exception as ex:
        out.write(f"{d['name']}\t{via}\t{smiles}\t{d['chem_A']}\t{d['chem_B']}\t\t\t\t\t"
                  f"{d['chem_pI']}\t\t\t计算失败: {type(ex).__name__}: {str(ex)[:60]}\n")
        out.flush()
        print(f"[{i+1}/{len(drugs)}] {d['name']}: CALC FAILED {type(ex).__name__}", flush=True)
        continue
    err_a = round(our_a - d["chem_A"], 2) if (our_a is not None and d["chem_A"] is not None) else ""
    err_b = round(our_b - d["chem_B"], 2) if (our_b is not None and d["chem_B"] is not None) else ""
    err_pi = round(our_pi - d["chem_pI"], 2) if (our_pi is not None and d["chem_pI"] is not None) else ""
    if err_a != "":
        done_A.append(abs(err_a))
    if err_b != "":
        done_B.append(abs(err_b))
    if err_pi != "":
        done_pI.append(abs(err_pi))
    out.write(f"{d['name']}\t{via}\t{smiles}\t{d['chem_A']}\t{d['chem_B']}\t"
              f"{our_a}\t{our_b}\t{err_a}\t{err_b}\t{d['chem_pI']}\t{our_pi}\t{err_pi}\t"
              f"{';'.join(note)}\n")
    out.flush()
    el = time.time() - t_start
    print(f"[{i+1}/{len(drugs)}] {d['name']}: A {our_a} (chem {d['chem_A']}, err {err_a}) | "
          f"B {our_b} (chem {d['chem_B']}, err {err_b}) | {el:.0f}s elapsed", flush=True)

import statistics as st
print("\n=== SUMMARY (abs error vs Chemicalize, all rows in TSV) ===", flush=True)
done_A, done_B, done_pI = [], [], []
for line in open(TSV, encoding="utf-8").read().splitlines()[1:]:
    f = line.split("\t")
    if len(f) < 13:
        continue
    for val, arr in ((f[7], done_A), (f[8], done_B), (f[11], done_pI)):
        try:
            arr.append(abs(float(val)))
        except ValueError:
            pass
for label, arr in (("pKa-A (most acidic)", done_A), ("pKa-B (most basic)", done_B),
                   ("pI", done_pI)):
    if arr:
        print(f"{label}: n={len(arr)} MAE={st.mean(arr):.2f} median={st.median(arr):.2f} "
              f"<=0.5: {sum(x<=0.5 for x in arr)}  <=1: {sum(x<=1 for x in arr)}", flush=True)
    else:
        print(f"{label}: n=0", flush=True)
out.close()
print("TSV written:", TSV, flush=True)
