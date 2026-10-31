"""Local web server for the pKa replica — a small Rowan-style macropka UI.

Run:  python -m pka_replica.server  (then open http://127.0.0.1:8765)
"""

from __future__ import annotations

import json
import math
import os
import re
import sys
import threading
import time
import urllib.parse
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from rdkit import Chem

from . import pubchem
from .depict import depict_mol, depict_smiles
from .enumerator import _edit_charge
from .workflow import run_macropka, run_micropka

def _resource_dir() -> str:
    """Data-file root: source tree normally, PyInstaller _MEIPASS when frozen."""
    base = getattr(sys, "_MEIPASS", None)
    if base and os.path.isdir(os.path.join(base, "web")):
        return base
    return os.path.dirname(__file__)


WEB_DIR = os.path.join(_resource_dir(), "web")
HL_COLOR = (0.99, 0.61, 0.16)  # warm orange highlight

# ---------------------------------------------------------------------------
# batch job store (txt/csv uploads -> sequential per-molecule computation)
# ---------------------------------------------------------------------------

_JOBS: dict = {}
_JOBS_LOCK = threading.Lock()
BATCH_MAX = 50

# session history of single-molecule predictions: re-opening or re-submitting
# the same molecule hits this cache instead of recomputing (server-lifetime)
_HISTORY: dict[str, dict] = {}      # id -> entry (ordered by insertion)
_HISTORY_ORDER: list[str] = []      # most recent first
_HIST_LOCK = threading.Lock()
_HIST_MAX = 50


def _history_put(entry: dict) -> None:
    """Store/refresh a history entry, most-recent-first, capped."""
    with _HIST_LOCK:
        hid = entry["id"]
        if hid in _HISTORY:
            _HISTORY_ORDER.remove(hid)
        _HISTORY[hid] = entry
        _HISTORY_ORDER.insert(0, hid)
        while len(_HISTORY_ORDER) > _HIST_MAX:
            old = _HISTORY_ORDER.pop()
            _HISTORY.pop(old, None)


def _history_summary(macro: dict) -> dict:
    pks = [p["pKa"] for p in macro.get("pKa_values", []) if -2 < p["pKa"] < 16]
    return {"pKa": [round(p, 2) for p in pks[:8]],
            "pI": macro.get("isoelectric_point")}


def _xlsx_response(rows: list[dict]) -> bytes:
    """Chemicalize-style long sheets: one row per molecule, one column per
    0.1-pH step over 0-14, 4-decimal cells. rows: [{name, logd, charge}]
    where logd/charge are [[pH, value], ...] series from the macro result."""
    import io

    from openpyxl import Workbook

    phs = [round(0.1 * i, 1) for i in range(141)]  # 0.0 .. 14.0

    def grid(series):
        vals = {}
        for p, v in (series or []):
            try:
                vals[round(float(p), 1)] = float(v)
            except (TypeError, ValueError):
                continue
        out = []
        for ph in phs:
            v = vals.get(ph)
            if v is None:  # profiles stop at 13.9: carry the tail into 14.0
                v = vals.get(round(ph - 0.1, 1))
            out.append(None if v is None else round(v, 4))
        return out

    wb = Workbook()
    ws_l = wb.active
    ws_l.title = "LogD"
    ws_c = wb.create_sheet("isoelectric")
    hdr = ["name", "pH"] + [f"{p:.1f}" for p in phs]
    ws_l.append(hdr)
    ws_c.append(hdr)
    for r in rows:
        ws_l.append([r["name"], "logD"] + grid(r.get("logd")))
        ws_c.append([r["name"], "Charge"] + grid(r.get("charge")))
    for ws in (ws_l, ws_c):
        for row in ws.iter_rows(min_row=2, min_col=3):
            for c in row:
                c.number_format = "0.0000"
        ws.column_dimensions["A"].width = 28
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def parse_smiles_batch(text: str) -> list[dict]:
    """Parse a txt/csv batch: entries separated by newlines; fields by comma,
    semicolon or tab. If a line has several fields, the first one that parses
    as a valid SMILES is taken and any other non-empty field becomes the name.
    """
    entries = []
    for raw in re.split(r"[\r\n]+", text or ""):
        line = raw.strip()
        if not line:
            continue
        fields = ([f.strip().strip('"') for f in re.split(r"[,;\t]", line)]
                  if re.search(r"[,;\t]", line) else [line])
        picked = None
        for f in fields:
            if f and Chem.MolFromSmiles(f) is not None:
                picked = f
                break
        name = None
        if picked is not None:
            name = next((f for f in fields if f and f != picked), None)
        entries.append({"smiles": picked, "name": name,
                        "valid": picked is not None, "raw": line,
                        "fields": [f for f in fields if f] if len(fields) > 1 else None})
    return entries


def _run_batch(job_id: str, backend: str) -> None:
    job = _JOBS[job_id]
    for item in job["items"]:
        if not item["valid"]:
            continue
        with _JOBS_LOCK:
            item["status"] = "running"
        try:
            macro = run_macropka(item["smiles"], backend=backend)
            macro["object_data"]["input_query"] = (item.get("resolved_from")
                                                   or item["smiles"])
            if item.get("resolved_info"):
                macro["resolved_info"] = item["resolved_info"]
            micro = run_micropka(item["smiles"], backend=backend,
                                 site_pka=macro.get("site_pka_full"))
            with _JOBS_LOCK:
                item["status"] = "done"
                item["summary"] = {
                    "pKa_values": macro["pKa_values"],
                    "isoelectric_point": macro["isoelectric_point"],
                    "energy_backend": macro["object_data"]["energy_backend"],
                    "n_states": len(macro["microstates"]),
                }
                item["macro"] = macro
                item["micro"] = micro
        except Exception as ex:
            with _JOBS_LOCK:
                item["status"] = "failed"
                item["error"] = f"{type(ex).__name__}: {ex}"
    with _JOBS_LOCK:
        job["status"] = "done"


def _job_payload(job: dict, include_full: bool = False) -> dict:
    with _JOBS_LOCK:
        items = []
        for i, it in enumerate(job["items"]):
            entry = {
                "index": i,
                "name": it.get("name"),
                "smiles": it.get("smiles"),
                "raw": it.get("raw"),
                "valid": it["valid"],
                "status": it["status"],
                "svg": it.get("svg"),
                "summary": it.get("summary"),
                "error": it.get("error"),
            }
            if include_full:
                entry["macro"] = it.get("macro")
                entry["micro"] = it.get("micro")
            items.append(entry)
        return {"status": job["status"], "items": items}


def _input_order_mol(object_data: dict) -> Chem.Mol | None:
    """Rebuild the mol in the atom order the site atom_index values refer to.

    MolToSmiles(canonical=False) is NOT order-preserving once aromatic rings
    get rewritten (PubChem Kekulé inputs silently permuted the highlights);
    a molblock round-trip preserves index order exactly. Falls back to the
    legacy SMILES fields for results produced by older code.
    """
    mb = object_data.get("input_molblock")
    if mb:
        try:
            m = Chem.MolFromMolBlock(mb)
            if m is not None:
                return m
        except Exception:
            pass
    s = object_data.get("input_smiles_atomorder") or object_data.get("initial_smiles")
    return Chem.MolFromSmiles(s) if s else None


def site_pair_images(initial_mol: Chem.Mol, atom_index: int, delta: int,
                     fallback_smiles: str | None = None) -> dict:
    """Rowan-style site row: parent structure with the site atom highlighted,
    and the conjugate structure with the same atom highlighted.

    `initial_mol` MUST be the mol in the order the site atom_index values
    refer to (see _input_order_mol). Fallbacks ensure the conjugate is always
    shown: unhighlighted depiction of the conjugate SMILES.
    """
    if initial_mol is None or atom_index >= initial_mol.GetNumAtoms():
        return {"parent_svg": None, "conjugate_svg": None}
    prod = _edit_charge(Chem.Mol(initial_mol), atom_index, delta)
    conj_svg = depict_mol(prod, width=210, height=145, highlight=[atom_index])
    if conj_svg is None and fallback_smiles:
        conj_svg = depict_smiles(fallback_smiles, width=210, height=145)
    return {
        "parent_svg": depict_mol(initial_mol, width=210, height=145,
                                 highlight=[atom_index]),
        "conjugate_svg": conj_svg,
    }


def enrich_micro_images(macro: dict, micro: dict) -> None:
    """Attach parent -> conjugate structure images to every site row of a
    micropka result (shared by /api/predict and the batch item endpoint).
    """
    ref = _input_order_mol(macro["object_data"])
    for row, delta in ([(r, -1) for r in micro["conjugate_bases"]] +
                       [(r, +1) for r in micro["conjugate_acids"]]):
        row["images"] = site_pair_images(
            ref, row["atom_index"], delta,
            fallback_smiles=row.get("conjugate"))


def charge_level_cards(macro: dict, max_charges: int = 8) -> list[dict]:
    """Chemicalize-style macrospecies: ONE card + ONE curve per net charge level.

    Users think in charge levels (-1..+3), not in the dozens of microstates
    (tautomers / protonation-site variants) inside each level — the previous
    per-microstate legend was unreadable. Returns, per charge: aggregate peak,
    the pH window where this level is the dominant one, and the dominant
    microstate (depicted).
    """
    states = macro["microstates"]
    prof = macro["microstate_weights_by_pH"]  # [[pH, [w...]], ...]
    if not states or not prof:
        return []
    by_charge: dict[int, list[int]] = {}
    for i, s in enumerate(states):
        by_charge.setdefault(s["charge"], []).append(i)
    charges = sorted(by_charge, reverse=True)[:max_charges]

    # aggregate weight per charge level, per pH row
    agg = {z: [sum(ws[i] for i in by_charge[z]) for _, ws in prof] for z in charges}

    out = []
    for z in charges:
        a = agg[z]
        peak = max(a)
        # pH window(s) where this level is the largest among all levels
        runs: list[list[float]] = []
        cur = None
        for k, (pH, _) in enumerate(prof):
            best = max(charges, key=lambda zz: agg[zz][k])
            if best == z:
                if cur and abs(pH - cur[1]) < 0.15:
                    cur[1] = pH
                else:
                    if cur:
                        runs.append(cur)
                    cur = [pH, pH]
            else:
                # another level dominates here — close any open run
                if cur:
                    runs.append(cur)
                cur = None
        if cur:
            runs.append(cur)
        lo = min(r[0] for r in runs) if runs else None
        hi = max(r[1] for r in runs) if runs else None
        # dominant microstate within this charge level (at the peak pH)
        k_peak = a.index(peak)
        dom_idx = max(by_charge[z], key=lambda i: prof[k_peak][1][i])
        st = states[dom_idx]
        out.append({
            "charge": z,
            "peak": round(peak, 4),
            "ph_lo": None if lo is None else round(max(0.0, lo - 0.05), 1),
            "ph_hi": None if hi is None else min(14.0, math.floor((hi + 0.05) * 10) / 10),
            "index": dom_idx,
            "smiles": st["smiles"],
            "svg": depict_smiles(st["smiles"], width=200, height=140),
        })
    return out


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):  # silence default logging
        pass

    def _send(self, code: int, body: bytes, ctype: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        try:
            self._do_get()
        except Exception as ex:
            # an unhandled exception would otherwise just drop the connection;
            # always answer with JSON so the UI can show something
            try:
                self._send(500, json.dumps(
                    {"ok": False, "error": f"{type(ex).__name__}: {ex}"},
                ).encode("utf-8"), "application/json")
            except Exception:
                pass

    def _do_get(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        qs = urllib.parse.parse_qs(parsed.query)
        if path in ("/", "/index.html"):
            with open(os.path.join(WEB_DIR, "index.html"), "rb") as fh:
                self._send(200, fh.read(), "text/html; charset=utf-8")
            return
        if path == "/api/batch/status":
            job = _JOBS.get(qs.get("id", [""])[0])
            if job is None:
                self._send(404, b'{"ok": false, "error": "unknown job"}',
                           "application/json")
                return
            payload = _job_payload(job)
            payload["ok"] = True
            self._send(200, json.dumps(payload).encode("utf-8"), "application/json")
            return
        if path == "/api/batch/item":
            job = _JOBS.get(qs.get("id", [""])[0])
            try:
                idx = int(qs.get("index", ["-1"])[0])
            except ValueError:
                idx = -1
            if job is None or not (0 <= idx < len(job["items"])):
                self._send(404, b'{"ok": false, "error": "unknown item"}',
                           "application/json")
                return
            it = job["items"][idx]
            if it["status"] != "done":
                self._send(409, b'{"ok": false, "error": "not ready"}',
                           "application/json")
                return
            # heavy RDKit rendering must stay OUTSIDE _JOBS_LOCK: rendering
            # while holding it stalls every /api/batch/status poll
            need_images = False
            with _JOBS_LOCK:
                # site-pair images are generated lazily, once, for molecules
                # actually viewed (batch computation itself skips them)
                if it["status"] == "done" and not it.get("_images_done"):
                    need_images = True
            if need_images:
                enrich_micro_images(it["macro"], it["micro"])
                with _JOBS_LOCK:
                    it["_images_done"] = True
            levels = charge_level_cards(it["macro"])
            payload = {"ok": True, "macro": it["macro"], "micro": it["micro"],
                       "charge_levels": levels,
                       "resolved_info": it.get("resolved_info")}
            self._send(200, json.dumps(payload).encode("utf-8"), "application/json")
            return
        if path == "/api/history":
            with _HIST_LOCK:
                items = [{k: e[k] for k in ("id", "query", "smiles", "ts", "svg",
                                            "summary")} for hid in _HISTORY_ORDER
                         for e in (_HISTORY[hid],)]
            self._send(200, json.dumps({"ok": True, "items": items}).encode("utf-8"),
                       "application/json")
            return
        if path == "/api/history/item":
            hid = qs.get("id", [""])[0]
            with _HIST_LOCK:
                entry = _HISTORY.get(hid)
            if entry is None:
                self._send(404, b'{"ok": false, "error": "unknown history id"}',
                           "application/json")
                return
            # same lazy-image pattern as /api/batch/item
            if not entry.get("_images_done"):
                enrich_micro_images(entry["macro"], entry["micro"])
                with _HIST_LOCK:
                    entry["_images_done"] = True
            payload = {"ok": True, "macro": entry["macro"], "micro": entry["micro"],
                       "charge_levels": entry["levels"],
                       "resolved_info": entry.get("resolved"),
                       "from_cache": True}
            self._send(200, json.dumps(payload).encode("utf-8"), "application/json")
            return
        if path == "/api/resolve":
            query = qs.get("q", [""])[0].strip()
            if not query:
                self._send(400, b'{"ok": false, "error": "empty query"}',
                           "application/json")
                return
            try:
                info = pubchem.resolve(query)
                self._send(200, json.dumps({"ok": True, **info}).encode("utf-8"),
                           "application/json")
            except pubchem.PubChemError as e:
                self._send(200, json.dumps({"ok": False,
                                            "error": str(e)}).encode("utf-8"),
                           "application/json")
            return
        self._send(404, b"not found", "text/plain")

    def do_POST(self) -> None:
        if self.path not in ("/api/predict", "/api/batch/start", "/api/resolve",
                             "/api/export"):
            self._send(404, b"not found", "text/plain")
            return
        try:
            n = int(self.headers.get("Content-Length", 0))
            req = json.loads(self.rfile.read(n).decode("utf-8"))
        except Exception:
            self._send(400, b"bad json", "text/plain")
            return
        if not isinstance(req, dict):
            self._send(400, json.dumps(
                {"ok": False, "error": "请求体必须是 JSON 对象"}).encode("utf-8"),
                "application/json")
            return
        try:
            self._do_post(req)
        except Exception as ex:
            try:
                self._send(500, json.dumps(
                    {"ok": False, "error": f"{type(ex).__name__}: {ex}"},
                ).encode("utf-8"), "application/json")
            except Exception:
                pass

    def _do_post(self, req: dict) -> None:
        backend = str(req.get("backend", "auto"))

        if self.path == "/api/export":
            # IUPAC-first naming for export rows: explicit name loses to
            # PubChem IUPAC; missing identities get a best-effort lookup
            # within a time budget so a big batch cannot hang the export
            deadline = time.time() + 20.0

            def row_name(explicit, macro, fallback):
                ri = (macro or {}).get("resolved_info") or {}
                if ri.get("iupac_name"):
                    return str(ri["iupac_name"])
                if time.time() < deadline:
                    smi = ((macro or {}).get("object_data") or {}).get(
                        "initial_smiles") or fallback
                    if smi:
                        try:
                            iu = pubchem.lookup_smiles(smi).get("iupac_name")
                            if iu:
                                return str(iu)
                        except Exception:
                            pass
                return str(explicit or fallback or "mol")

            rows = []
            if req.get("history"):
                with _HIST_LOCK:
                    entries = [_HISTORY[h] for h in _HISTORY_ORDER]
                for e in entries:
                    rows.append({
                        "name": row_name(e.get("query"), e.get("macro"),
                                         e.get("query")),
                        "logd": e["macro"].get("logD_by_pH"),
                        "charge": e["macro"].get("average_charge_by_pH"),
                    })
            if req.get("job_id"):
                with _JOBS_LOCK:
                    job = _JOBS.get(req["job_id"])
                    items = [dict(it) for it in job["items"]] if job else []
                for it in items:
                    if it.get("status") != "done" or "macro" not in it:
                        continue
                    rows.append({
                        "name": row_name(it.get("name"), it["macro"],
                                         it.get("smiles")),
                        "logd": it["macro"].get("logD_by_pH"),
                        "charge": it["macro"].get("average_charge_by_pH"),
                    })
            for rw in (req.get("rows") or []):
                rows.append({"name": str(rw.get("name") or "mol"),
                             "logd": rw.get("logd"),
                             "charge": rw.get("charge")})
            if not rows:
                # 409 so downloadExport's r.ok check fails and the UI shows
                # the message instead of saving a corrupt .xlsx
                self._send(409, json.dumps(
                    {"ok": False, "error": "没有可导出的数据"}
                ).encode("utf-8"), "application/json")
                return
            data = _xlsx_response(rows)
            fname = f"pka_export_{time.strftime('%Y%m%d_%H%M%S')}.xlsx"
            self.send_response(200)
            self.send_header("Content-Type",
                             "application/vnd.openxmlformats-officedocument"
                             ".spreadsheetml.sheet")
            self.send_header("Content-Disposition",
                             f'attachment; filename="{fname}"')
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return

        if self.path == "/api/batch/start":
            entries = parse_smiles_batch(str(req.get("text", "")))
            # try PubChem for entries that are not valid SMILES (CAS / names)
            if req.get("resolve", True):
                attempts = 0
                # bound the worst case: sequential network lookups with
                # retries inside the request thread would otherwise hang
                # the client for minutes during a PubChem outage
                deadline = time.time() + 25.0
                for e in entries:
                    if e["valid"] or attempts >= 30 or time.time() > deadline:
                        continue
                    # multi-field lines (e.g. "50-78-2,aspirin"): try each
                    # field separately — the raw line never resolves
                    for cand in (e.get("fields") or [e["raw"]]):
                        if attempts >= 30 or time.time() > deadline:
                            break
                        attempts += 1
                        try:
                            info = pubchem.resolve(cand)
                            e["smiles"] = info["smiles"]
                            e["name"] = e["raw"]
                            e["valid"] = True
                            e["resolved_from"] = cand
                            e["resolved_info"] = info
                            break
                        except Exception:
                            continue
            valid = [e for e in entries if e["valid"]]
            if not valid:
                self._send(200, json.dumps({
                    "ok": False,
                    "error": "没有解析到有效的 SMILES(支持换行/逗号/分号/制表符分隔)",
                }).encode("utf-8"), "application/json")
                return
            if len(valid) > BATCH_MAX:
                self._send(200, json.dumps({
                    "ok": False,
                    "error": f"批量上限 {BATCH_MAX} 个分子(当前 {len(valid)} 个)",
                }).encode("utf-8"), "application/json")
                return
            job_id = uuid.uuid4().hex[:12]
            items = []
            for e in valid:
                items.append({
                    "smiles": e["smiles"], "name": e["name"], "raw": e["raw"],
                    "valid": True, "status": "queued",
                    "resolved_from": e.get("resolved_from"),
                    "resolved_info": e.get("resolved_info"),
                    "svg": depict_smiles(e["smiles"], width=190, height=120),
                })
            with _JOBS_LOCK:
                _JOBS[job_id] = {"status": "running", "items": items,
                                 "started": time.time()}
            threading.Thread(target=_run_batch, args=(job_id, backend),
                             daemon=True).start()
            payload = {"ok": True, "job_id": job_id,
                       "invalid": [e["raw"] for e in entries if not e["valid"]],
                       "total": len(items)}
            self._send(200, json.dumps(payload).encode("utf-8"), "application/json")
            return

        smiles = str(req.get("smiles", "")).strip()
        try:
            def progress(done, total, phase):
                print(f"  [{done}/{total}] {phase}", flush=True)

            resolved = None
            orig_query = smiles  # what the user actually typed (name/CAS/SMILES)
            if Chem.MolFromSmiles(smiles) is None:
                # not a SMILES -> treat as CAS number or compound name
                try:
                    info = pubchem.resolve(smiles)
                    resolved = info
                    smiles = info["smiles"]
                except pubchem.PubChemError as e:
                    self._send(200, json.dumps({
                        "ok": False,
                        "error": f"输入既不是有效 SMILES,PubChem 也未能解析({e})",
                    }).encode("utf-8"), "application/json")
                    return

            # session history: same canonical molecule + SAME params ->
            # serve the cached result instead of recomputing
            key = None
            _m0 = Chem.MolFromSmiles(smiles)
            if _m0 is not None:
                key = Chem.MolToSmiles(_m0)
            params = (backend,
                      json.dumps(req.get("min_charge", "auto"), sort_keys=True),
                      json.dumps(req.get("max_charge", "auto"), sort_keys=True),
                      float(req.get("min_pH", 0.0)),
                      float(req.get("max_pH", 14.0)),
                      int(req.get("n_confs", 2)))
            if key is not None:
                with _HIST_LOCK:
                    hit = next((e for e in (_HISTORY[h] for h in _HISTORY_ORDER)
                                if e["key"] == key and e.get("params") == params),
                               None)
                    if hit is not None:
                        # re-promote under the same lock hold (Lock is NOT
                        # reentrant — a nested acquire here deadlocked the
                        # handler thread permanently)
                        _HISTORY_ORDER.remove(hit["id"])
                        _HISTORY_ORDER.insert(0, hit["id"])
                        images_done = hit.get("_images_done")
                if hit is not None and not images_done:
                    enrich_micro_images(hit["macro"], hit["micro"])
                    with _HIST_LOCK:
                        hit["_images_done"] = True
                if hit is not None:
                    payload = {"ok": True, "macro": hit["macro"],
                               "micro": hit["micro"], "charge_levels": hit["levels"],
                               "resolved": hit.get("resolved"),
                               "history_id": hit["id"], "from_cache": True,
                               "query": smiles}
                    self._send(200, json.dumps(payload).encode("utf-8"),
                               "application/json")
                    return

            macro = run_macropka(
                smiles,
                min_charge=req.get("min_charge", "auto"),
                max_charge=req.get("max_charge", "auto"),
                min_pH=float(req.get("min_pH", 0.0)),
                max_pH=float(req.get("max_pH", 14.0)),
                backend=backend,
                n_confs=int(req.get("n_confs", 2)),
                progress=progress,
            )
            macro["object_data"]["input_query"] = smiles
            if resolved:
                macro["resolved_info"] = resolved
            else:
                # direct SMILES input: best-effort PubChem identity so the
                # physchem card and export names can show the IUPAC name
                # (silent on failure; payload["resolved"] stays None so the
                # UI does not show a "resolved via PubChem" status)
                try:
                    macro["resolved_info"] = pubchem.lookup_smiles(smiles)
                except Exception:
                    pass
            micro = run_micropka(smiles, backend=backend,
                                 site_pka=macro.get("site_pka_full"))

            # Rowan-style site rows: parent -> conjugate structure images.
            enrich_micro_images(macro, micro)

            levels = charge_level_cards(macro)
            hid = uuid.uuid4().hex[:12]
            _history_put({
                "id": hid, "key": key, "backend": backend, "params": params,
                "query": orig_query,
                "display": (resolved or {}).get("iupac_name") or smiles,
                "smiles": macro["object_data"]["initial_smiles"],
                "ts": int(time.time()),
                "svg": depict_smiles(macro["object_data"]["initial_smiles"],
                                     width=190, height=120),
                "summary": _history_summary(macro),
                "macro": macro, "micro": micro, "levels": levels,
                "resolved": resolved, "_images_done": True,
            })
            payload = {
                "macro": macro,
                "micro": micro,
                "charge_levels": levels,
                "resolved": resolved,
                "history_id": hid,
                "ok": True,
            }
        except Exception as ex:
            import traceback

            traceback.print_exc()
            payload = {"ok": False, "error": f"{type(ex).__name__}: {ex}"}
        self._send(200, json.dumps(payload).encode("utf-8"), "application/json")


def main(port: int = 8765) -> None:
    # politeness: if the port is taken, walk upward instead of dying
    import socket

    for _ in range(10):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            if s.connect_ex(("127.0.0.1", port)) != 0:
                break  # port free
        port += 1
    server = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    url = f"http://127.0.0.1:{port}"
    print(f"pka_replica web UI on {url}  (Ctrl+C to stop)", flush=True)
    threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
