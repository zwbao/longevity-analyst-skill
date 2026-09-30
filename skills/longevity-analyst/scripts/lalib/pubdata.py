"""Live public-data clients. Every answer is cached in work/evidence/pub_<source>_<key>.json with the URL and the time
it was retrieved, and every record carries a `ref` that findings may cite. Nothing here interprets a result.

Sources: GWAS Catalog REST v2 (associations by trait), myvariant.info (dbSNP/ClinVar/gnomAD per variant),
Ensembl REST (gene coordinates), EpiGraphDB (published Mendelian-randomisation estimates).
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import ssl
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from .common import EXIT_EXTERNAL, LAError, load_json, now_iso, write_json

try:
    import certifi
    _CTX = ssl.create_default_context(cafile=certifi.where())
except Exception:  # noqa: BLE001
    _CTX = ssl.create_default_context()

UA = {"User-Agent": "longevity-analyst/0.6 (research use)", "Accept": "application/json"}


SHARED = Path.home() / ".cache" / "longevity-analyst" / "pubdata"
SHARED_TTL_DAYS = 30


def _shared_get(key: str) -> Optional[Dict[str, Any]]:
    """A public-database answer fetched earlier on this machine (any member), still within the TTL."""
    p = SHARED / (hashlib.sha1(key.encode()).hexdigest() + ".json")
    if not p.exists():
        return None
    try:
        d = load_json(p)
        from datetime import datetime, timezone
        age = datetime.now(timezone.utc) - datetime.fromisoformat(d["retrieved_at"])
        return d if age.days < SHARED_TTL_DAYS else None
    except Exception:  # noqa: BLE001
        return None


def _shared_put(key: str, doc: Dict[str, Any]) -> None:
    SHARED.mkdir(parents=True, exist_ok=True)
    write_json(SHARED / (hashlib.sha1(key.encode()).hexdigest() + ".json"), doc)


def _http(url: str, data: Optional[bytes] = None, headers: Optional[Dict[str, str]] = None, tries: int = 8) -> Any:
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={**UA, **(headers or {})})
            with urllib.request.urlopen(req, timeout=150, context=_CTX) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(min(3 * (i + 1), 20))
    raise LAError(f"could not retrieve {url.split('?')[0]}: {last}", EXIT_EXTERNAL)


def _cache(ws: Path, source: str, key: str, url: str, records: List[Dict[str, Any]], extra: Optional[Dict] = None,
           retrieved_at: Optional[str] = None) -> Path:
    k = hashlib.sha1(key.encode()).hexdigest()[:12]
    p = ws / "work" / "evidence" / f"pub_{source}_{k}.json"
    doc = {"source": source, "query": key, "url": url, "retrieved_at": retrieved_at or now_iso(), "records": records, **(extra or {})}
    write_json(p, doc)
    if source in ("gwas", "myvariant", "ensembl", "mr") and not retrieved_at:
        _shared_put(f"{source}|{key}", doc)
    return p


def known_refs(ws: Path) -> Dict[str, Dict[str, Any]]:
    """Every citable public record retrieved in this workspace: ref -> record."""
    out: Dict[str, Dict[str, Any]] = {}
    for p in (ws / "work" / "evidence").glob("pub_*.json"):
        try:
            for r in load_json(p).get("records", []):
                if r.get("ref"):
                    out[r["ref"]] = r
        except Exception:  # noqa: BLE001
            continue
    return out


# ------------------------------------------------------------------ GWAS Catalog
_NUM = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _beta(text: Any) -> Optional[float]:
    """'1.7 ln nmol/L increase' -> 1.7; '0.3 unit decrease' -> -0.3; None when no number or no direction."""
    if text is None:
        return None
    t = str(text)
    m = _NUM.search(t)
    if not m:
        return None
    v = float(m.group(0))
    if re.search(r"decreas|lower|reduc", t, re.I):
        return -abs(v)
    if re.search(r"increas|higher|raise", t, re.I):
        return abs(v)
    return None


def gwas_trait(ws: Path, efo: str, pmax: float = 5e-8, max_pages: int = 3) -> Dict[str, Any]:
    """Associations for one trait, strongest first, down to pmax (the strongest 300 are enough to clump 40 loci).
    Records keep the reported effect text. A copy fetched on this machine within the TTL is reused with its original
    retrieval time."""
    key = f"trait:{efo}:{pmax}"
    hit = _shared_get(f"gwas|{key}")
    if hit:
        _cache(ws, "gwas", key, hit["url"], hit["records"], retrieved_at=hit["retrieved_at"])
        return {"efo": efo, "records": hit["records"], "cached_from": hit["retrieved_at"]}
    recs: List[Dict[str, Any]] = []
    base = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
    url = f"{base}?efo_id={efo}&size=100&sort=p_value&direction=asc"
    for page in range(max_pages):
        d = _http(f"{url}&page={page}")
        items = (d.get("_embedded") or {}).get("associations") or []
        stop = False
        for a in items:
            p = a.get("p_value")
            if p is None or p > pmax:
                stop = True
                break
            for sa in a.get("snp_allele") or []:
                rs = sa.get("rs_id") or ""
                ea = (sa.get("effect_allele") or "").upper()
                if not re.fullmatch(r"rs\d+", rs) or ea not in ("A", "C", "G", "T"):
                    continue
                b = _beta(a.get("beta"))
                orv = a.get("or_value")
                try:
                    orv = float(orv) if orv not in (None, "", "NR") else None
                except ValueError:
                    orv = None
                direction = (1 if b > 0 else -1) if b else ((1 if orv > 1 else -1) if orv and orv != 1 else 0)
                if not direction:
                    continue
                recs.append({"ref": f"gwas:{a.get('accession_id')}:{rs}-{ea}", "rsid": rs, "effect_allele": ea,
                             "p_value": p, "beta_text": a.get("beta"), "or": orv, "direction": direction,
                             "accession": a.get("accession_id"), "pubmed_id": str(a.get("pubmed_id") or ""),
                             "mapped_genes": a.get("mapped_genes") or [], "efo": efo,
                             "reported_trait": (a.get("reported_trait") or [""])[0]})
        if stop or len(items) < 100:
            break
    _cache(ws, "gwas", key, url, recs)
    return {"efo": efo, "records": recs}


# ------------------------------------------------------------------ myvariant.info
def myvariant_rsids(ws: Path, rsids: Iterable[str], assembly: str = "hg38") -> Dict[str, List[Dict[str, Any]]]:
    """rsID -> list of alleles with position (assembly), ref, alt, gnomAD genome EAS AF, ClinVar significance."""
    ids = sorted(set(rsids))
    key = f"rsids:{assembly}:{hashlib.sha1(','.join(ids).encode()).hexdigest()}"
    hit = _shared_get(f"myvariant|{key}")
    if hit:
        _cache(ws, "myvariant", key, hit["url"], hit["records"], retrieved_at=hit["retrieved_at"])
        out_c: Dict[str, List[Dict[str, Any]]] = {}
        for r in hit["records"]:
            out_c.setdefault(r["rsid"], []).append(r)
        return out_c
    out: Dict[str, List[Dict[str, Any]]] = {}
    recs: List[Dict[str, Any]] = []
    fields = "dbsnp.rsid,dbsnp.ref,dbsnp.alt,dbsnp.chrom,gnomad_genome.af.af_eas,gnomad_genome.af.af,clinvar.rcv.clinical_significance,clinvar.rcv.accession,clinvar.gene.symbol,chrom,vcf,hg19,hg38"
    for i in range(0, len(ids), 500):
        chunk = ids[i:i + 500]
        body = urllib.parse.urlencode({"q": ",".join(chunk), "scopes": "dbsnp.rsid", "fields": fields,
                                       "assembly": assembly, "size": 1000}).encode()
        res = _http("https://myvariant.info/v1/query", data=body,
                    headers={"Content-Type": "application/x-www-form-urlencoded"})
        for h in res:
            if h.get("notfound"):
                continue
            rs = h.get("query")
            v = h.get("vcf") or {}
            pos = (h.get(assembly) or {}).get("start") or v.get("position")
            eas = ((h.get("gnomad_genome") or {}).get("af") or {}).get("af_eas")
            cv = h.get("clinvar") or {}
            rcv = cv.get("rcv") or []
            rcv = rcv if isinstance(rcv, list) else [rcv]
            sig = sorted({str(x.get("clinical_significance")) for x in rcv if x.get("clinical_significance")})
            rec = {"ref": f"myvariant:{h.get('_id')}", "rsid": rs, "chrom": str(h.get("chrom") or ""), "pos": int(pos) if pos else None,
                   "ref_allele": (v.get("ref") or "").upper(), "alt_allele": (v.get("alt") or "").upper(),
                   "af_eas": float(eas) if isinstance(eas, (int, float)) else None, "clinvar_significance": sig,
                   "clinvar_rcv": [x.get("accession") for x in rcv if x.get("accession")][:5], "assembly": assembly}
            out.setdefault(rs, []).append(rec)
            recs.append(rec)
    _cache(ws, "myvariant", key, "https://myvariant.info/v1/query", recs)
    return out


def myvariant_hgvs(ws: Path, hgvs_ids: Iterable[str], assembly: str = "hg38") -> Dict[str, Dict[str, Any]]:
    """Annotate the member's own variants (chrN:g.POSREF>ALT) with ClinVar and gnomAD EAS AF."""
    ids = sorted(set(hgvs_ids))
    out: Dict[str, Dict[str, Any]] = {}
    recs = []
    fields = "clinvar.rcv.clinical_significance,clinvar.rcv.accession,clinvar.rcv.conditions.name,clinvar.gene.symbol,gnomad_genome.af.af_eas,dbsnp.rsid"
    for i in range(0, len(ids), 500):
        chunk = ids[i:i + 500]
        body = urllib.parse.urlencode({"ids": ",".join(chunk), "fields": fields, "assembly": assembly}).encode()
        res = _http("https://myvariant.info/v1/variant", data=body, headers={"Content-Type": "application/x-www-form-urlencoded"})
        for h in res:
            if h.get("notfound"):
                continue
            cv = h.get("clinvar") or {}
            rcv = cv.get("rcv") or []
            rcv = rcv if isinstance(rcv, list) else [rcv]
            rec = {"ref": f"clinvar:{h.get('_id')}" if rcv else f"myvariant:{h.get('_id')}", "hgvs": h.get("query") or h.get("_id"),
                   "gene": (cv.get("gene") or {}).get("symbol"), "rsid": (h.get("dbsnp") or {}).get("rsid"),
                   "clinvar_significance": sorted({str(x.get("clinical_significance")) for x in rcv if x.get("clinical_significance")}),
                   "clinvar_rcv": [x.get("accession") for x in rcv if x.get("accession")][:5],
                   "conditions": sorted({str(((x.get("conditions") or {}) if isinstance(x.get("conditions"), dict) else {}).get("name", ""))
                                         for x in rcv} - {""})[:5],
                   "af_eas": ((h.get("gnomad_genome") or {}).get("af") or {}).get("af_eas")}
            out[rec["hgvs"]] = rec
            recs.append(rec)
    _cache(ws, "clinvar", f"hgvs:{assembly}:{hashlib.sha1(','.join(ids).encode()).hexdigest()}", "https://myvariant.info/v1/variant", recs)
    return out


# ------------------------------------------------------------------ Ensembl
def gene_region(ws: Path, symbol: str, assembly: str = "GRCh38") -> Optional[Dict[str, Any]]:
    host = "https://rest.ensembl.org" if assembly == "GRCh38" else "https://grch37.rest.ensembl.org"
    url = f"{host}/lookup/symbol/homo_sapiens/{urllib.parse.quote(symbol)}?content-type=application/json"
    try:
        d = _http(url)
    except LAError:
        return None
    rec = {"ref": f"ensembl:{d.get('id')}", "symbol": symbol, "chrom": str(d.get("seq_region_name")), "start": d.get("start"),
           "end": d.get("end"), "assembly": assembly}
    _cache(ws, "ensembl", f"gene:{assembly}:{symbol}", url, [rec])
    return rec


# ------------------------------------------------------------------ EpiGraphDB MR
def mr(ws: Path, exposure: str, outcome: str, pval: float = 1e-5) -> List[Dict[str, Any]]:
    """Published Mendelian-randomisation estimates (MR-EvE, EpiGraphDB) for exposure -> outcome trait names."""
    url = "https://api.epigraphdb.org/mr?" + urllib.parse.urlencode({"exposure_trait": exposure, "outcome_trait": outcome,
                                                                      "pval_threshold": pval})
    d = _http(url)
    recs = []
    for r in d.get("results") or []:
        e, o, m = r.get("exposure") or {}, r.get("outcome") or {}, r.get("mr") or {}
        if m.get("b") is None or m.get("se") is None:
            continue
        recs.append({"ref": f"mr:{e.get('id')}->{o.get('id')}:{m.get('method')}", "exposure_id": e.get("id"),
                     "exposure": e.get("trait"), "outcome_id": o.get("id"), "outcome": o.get("trait"), "b": m.get("b"),
                     "se": m.get("se"), "pval": m.get("pval"), "method": m.get("method"), "moescore": m.get("moescore"),
                     "scale_note": "b is the published MR estimate from the outcome GWAS (log odds per SD of exposure for binary outcomes)"})
    _cache(ws, "mr", f"mr:{exposure}->{outcome}:{pval}", url, recs)
    return recs


def se_ci(b: float, se: float) -> List[float]:
    return [b - 1.96 * se, b + 1.96 * se]


def log_or(x: Optional[float]) -> Optional[float]:
    return math.log(x) if x and x > 0 else None
