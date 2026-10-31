"""PFAS-Analytes comparison: our pKa tool (xtb + fragment routes) vs
Chemicalize pKa-A for 11 perfluoro acids.

The strongest acid is read from site_pka_full (pre-display-filter) because
the 0-14 display screen intentionally hides superacid values like -3.
"""
import sys
import time
import faulthandler

sys.stdout.reconfigure(encoding="utf-8")
faulthandler.enable()

from openpyxl import load_workbook
from pka_replica import pubchem
from pka_replica.workflow import run_macropka, run_micropka

TSV = r"D:\zcode_rea\pka_replica\pfas_result.tsv"


def cell_float(v):
    if v is None:
        return None
    s = str(v).strip()
    if not s or s.lower() in ("n/a", "null", "none"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def acid_from_site_pka(site_pka):
    """min over acid entries ((idx,-1)) of the raw calibrated site table."""
    if not site_pka:
        return None
    vals = []
    for k, v in site_pka.items():
        d = int(str(k).rsplit(",", 1)[1]) if isinstance(k, str) else k[1]
        if d == -1:
            vals.append(float(v))
    return min(vals) if vals else None


wb = load_workbook(r"D:\zcode_rea\PFAS-Analytes.xlsx", read_only=True)
rows = list(wb["drug"].iter_rows(values_only=True))
subs = []
for r in rows[1:]:
    name = str(r[0] or "").strip()
    if not name:
        continue
    subs.append({"name": name, "cas": str(r[1] or "").strip(),
                 "chemA": cell_float(r[25])})
print(f"loaded {len(subs)} PFAS substances", flush=True)

done = set()
try:
    for line in open(TSV, encoding="utf-8").read().splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 6 and f[0] and "失败" not in f[5]:
            done.add(f[0])
except FileNotFoundError:
    pass

out = open(TSV, "a", encoding="utf-8")
if out.tell() == 0:
    out.write("name\tsmiles\tchemA\tourA_xtb\terr_xtb\tourA_frag\terr_frag\tnote\n")
    out.flush()

import statistics as st
errs_xtb, errs_frag = [], []
t0 = time.time()
for i, d in enumerate(subs):
    if d["name"] in done:
        continue
    note = ""
    try:
        info = pubchem.resolve(d["name"])
    except Exception:
        try:
            info = pubchem.resolve(d["cas"])
        except Exception as e:
            out.write(f"{d['name']}\t\t{d['chemA']}\t\t\t\t\t解析失败: {str(e)[:50]}\n")
            out.flush()
            print(f"[{i+1}/{len(subs)}] {d['name']}: RESOLVE FAILED", flush=True)
            continue
    smiles = info["smiles"]
    try:
        macro = run_macropka(smiles, backend="xtb")
        micro = run_micropka(smiles, backend="xtb",
                             site_pka=macro.get("site_pka_full"))
        our_xtb = micro.get("strongest_acid")
        if our_xtb is None:
            our_xtb = acid_from_site_pka(macro.get("site_pka_full"))
            if our_xtb is not None:
                note = "取自未筛位点表(超酸<0)"
        micro_f = run_micropka(smiles, backend="sites")
        our_frag = micro_f.get("strongest_acid")
    except Exception as ex:
        out.write(f"{d['name']}\t{smiles}\t{d['chemA']}\t\t\t\t\t计算失败: {type(ex).__name__}: {str(ex)[:50]}\n")
        out.flush()
        print(f"[{i+1}/{len(subs)}] {d['name']}: CALC FAILED {type(ex).__name__}", flush=True)
        continue
    ex_ = round(our_xtb - d["chemA"], 2) if (our_xtb is not None and d["chemA"] is not None) else ""
    ef_ = round(our_frag - d["chemA"], 2) if (our_frag is not None and d["chemA"] is not None) else ""
    if ex_ != "":
        errs_xtb.append(abs(ex_))
    if ef_ != "":
        errs_frag.append(abs(ef_))
    out.write(f"{d['name']}\t{smiles}\t{d['chemA']}\t{our_xtb}\t{ex_}\t{our_frag}\t{ef_}\t{note}\n")
    out.flush()
    print(f"[{i+1}/{len(subs)}] {d['name']}: xtb {our_xtb} (err {ex_}) | "
          f"frag {our_frag} (err {ef_}) | {time.time()-t0:.0f}s", flush=True)

# summary over the whole TSV
errs_xtb, errs_frag = [], []
for line in open(TSV, encoding="utf-8").read().splitlines()[1:]:
    f = line.split("\t")
    if len(f) < 8:
        continue
    for val, arr in ((f[4], errs_xtb), (f[6], errs_frag)):
        try:
            arr.append(abs(float(val)))
        except ValueError:
            pass
print("\n=== SUMMARY |our-chem| ===", flush=True)
for lab, arr in (("xtb route", errs_xtb), ("fragment route", errs_frag)):
    if arr:
        print(f"{lab}: n={len(arr)} MAE={st.mean(arr):.2f} median={st.median(arr):.2f} "
              f"<=0.5:{sum(x<=0.5 for x in arr)} <=1:{sum(x<=1 for x in arr)}", flush=True)
