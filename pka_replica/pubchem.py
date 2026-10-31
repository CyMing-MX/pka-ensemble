"""PubChem PUG REST lookups: CAS number or compound name -> SMILES.

Results are cached in-process; network failures raise PubChemError with a
user-facing message.
"""

from __future__ import annotations

import re
import time

import requests

CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
_CACHE: dict[str, dict] = {}
_MISS = object()  # sentinel: transient lookup failures are retried next call


class PubChemError(Exception):
    pass


def looks_like_cas(query: str) -> bool:
    return bool(CAS_RE.match(query.strip()))


def is_smiles(query: str) -> bool:
    from rdkit import Chem

    return Chem.MolFromSmiles(query.strip()) is not None


def _cas_for_cid(cid) -> str | None:
    """CAS registry number for a CID: xrefs/RN first, CAS-shaped synonyms as
    fallback. Best effort (transient failures -> None); result cached so a
    session does not re-ask for the same compound."""
    if not cid:
        return None
    ck = f"cas:{cid}"
    hit = _CACHE.get(ck, _MISS)
    if hit is not _MISS:
        return None if hit is None else hit
    rns = []
    try:
        r = requests.get(f"{_BASE}/cid/{cid}/xrefs/RN/JSON", timeout=20,
                         headers={"User-Agent": "pka-replica/0.6"})
        if r.status_code == 200:
            rns = r.json().get("IdentifierList", {}).get("RN") or []
    except Exception:
        rns = []
    std = [x for x in rns if CAS_RE.match(x)]
    cas = std[0] if std else None
    if cas is None:
        try:
            r = requests.get(f"{_BASE}/cid/{cid}/synonyms/JSON", timeout=20,
                             headers={"User-Agent": "pka-replica/0.6"})
            if r.status_code == 200:
                info = r.json().get("InformationList", {}).get("Information") or [{}]
                syns = info[0].get("Synonym") or []
                cas = next((s for s in syns if CAS_RE.match(s)), None)
        except Exception:
            pass
    _CACHE[ck] = cas if cas is not None else _MISS
    return cas


def lookup_smiles(smiles: str) -> dict:
    """Best-effort identity lookup by SMILES: CID + IUPAC name etc.

    Used to fill the physchem card and export names for direct-SMILES input
    (name/CAS queries already carry this from resolve()). Raises
    PubChemError when not found or unreachable; successes are cached.
    """
    smiles = smiles.strip()
    ck = "smiles:" + smiles
    if ck in _CACHE:
        return dict(_CACHE[ck])
    url = (f"{_BASE}/smiles/{requests.utils.quote(smiles, safe='')}/property/"
           "IUPACName,SMILES,MolecularFormula,MolecularWeight,ExactMass,TPSA,"
           "XLogP,Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,"
           "HeavyAtomCount,Complexity,Title/JSON")
    r = None
    last = None
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=20,
                             headers={"User-Agent": "pka-replica/0.6"})
            break
        except Exception as e:
            last = e
            time.sleep(1.0)
    else:
        raise PubChemError(f"PubChem 网络请求失败(已重试 2 次): {last}") from last
    if r.status_code in (404, 400):
        raise PubChemError(f"PubChem 未收录该 SMILES")
    try:
        r.raise_for_status()
        props = r.json()["PropertyTable"]["Properties"][0]
    except PubChemError:
        raise
    except Exception as e:
        raise PubChemError(f"PubChem 返回格式异常: {e}") from e
    out = {"query": smiles,
           "smiles": props.get("SMILES") or smiles,
           "cid": props.get("CID"),
           "iupac_name": props.get("IUPACName"),
           "formula": props.get("MolecularFormula"),
           "mw": props.get("MolecularWeight"),
           "cas": _cas_for_cid(props.get("CID")),
           # PubChem CACTVS computed descriptors (their values, not ours)
           "xlogp": props.get("XLogP"),
           "tpsa": props.get("TPSA"),
           "charge": props.get("Charge"),
           "hbd": props.get("HBondDonorCount"),
           "hba": props.get("HBondAcceptorCount"),
           "rotb": props.get("RotatableBondCount"),
           "heavy_atoms": props.get("HeavyAtomCount"),
           "exact_mass": props.get("ExactMass"),
           "complexity": props.get("Complexity"),
           "title": props.get("Title")}
    _CACHE[ck] = out
    return dict(out)


def resolve(query: str) -> dict:
    """Resolve a CAS number or compound name to a SMILES via PubChem.

    Returns {"query", "smiles", "cid", "mw"}; raises PubChemError when the
    compound is not found or the service is unreachable.
    """
    query = query.strip()
    if query in _CACHE:
        return dict(_CACHE[query])
    url = (f"{_BASE}/name/{requests.utils.quote(query)}/property/"
           "IUPACName,SMILES,ConnectivitySMILES,MolecularFormula,MolecularWeight,"
           "ExactMass,TPSA,XLogP,Charge,HBondDonorCount,HBondAcceptorCount,"
           "RotatableBondCount,HeavyAtomCount,Complexity,Title/JSON")
    r = None
    for attempt in range(3):  # VPN/proxy resets are transient; retry twice
        try:
            r = requests.get(url, timeout=20,
                             headers={"User-Agent": "pka-replica/0.6"})
            break
        except Exception as e:
            last = e
            time.sleep(1.5)
    else:
        raise PubChemError(f"PubChem 网络请求失败(已重试 2 次): {last}") from last
    if r.status_code == 404:
        raise PubChemError(f"PubChem 未收录: {query}")
    if r.status_code == 400:
        # PubChem returns 400 for unknown names too
        raise PubChemError(f"PubChem 未收录或无法解析: {query}")
    try:
        r.raise_for_status()
    except Exception as e:
        # 429/503 etc. must arrive as PubChemError, or /api/resolve and the
        # CLI crash with a raw requests traceback instead of a usable message
        raise PubChemError(f"PubChem 服务暂时不可用(HTTP {r.status_code}): {query}") from e
    try:
        props = r.json()["PropertyTable"]["Properties"][0]
    except Exception as e:
        raise PubChemError(f"PubChem 返回格式异常: {e}") from e
    # newer PubChem API uses SMILES (isomeric) / ConnectivitySMILES; older
    # records may still expose IsomericSMILES / CanonicalSMILES
    smiles = (props.get("SMILES") or props.get("IsomericSMILES")
              or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES"))
    if not smiles:
        raise PubChemError(f"PubChem 无该化合物的 SMILES: {query}")
    out = {"query": query, "smiles": smiles,
           "cid": props.get("CID"),
           "iupac_name": props.get("IUPACName"),
           "formula": props.get("MolecularFormula"),
           "mw": props.get("MolecularWeight"),
           "cas": _cas_for_cid(props.get("CID")),
           "xlogp": props.get("XLogP"),
           "tpsa": props.get("TPSA"),
           "charge": props.get("Charge"),
           "hbd": props.get("HBondDonorCount"),
           "hba": props.get("HBondAcceptorCount"),
           "rotb": props.get("RotatableBondCount"),
           "heavy_atoms": props.get("HeavyAtomCount"),
           "exact_mass": props.get("ExactMass"),
           "complexity": props.get("Complexity"),
           "title": props.get("Title")}
    # a CAS-shaped query is itself the CAS number (covers xref outages)
    if out["cas"] is None and CAS_RE.match(query):
        out["cas"] = query
    _CACHE[query] = out
    return dict(out)
