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
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from .common import EXIT_EXTERNAL, LAError, load_json, now_iso, sha256_file, write_json

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


WRITTEN: Dict[str, str] = {}    # cache files written by this process -> sha256 (merged into state.evidence_ledger on save)
_DOWN: Dict[str, str] = {}      # host -> first error in this run (circuit breaker: a down source is not retried per item)


def _http(url: str, data: Optional[bytes] = None, headers: Optional[Dict[str, str]] = None, tries: int = 5) -> Any:
    host = urllib.parse.urlparse(url).netloc
    if host in _DOWN:
        raise LAError(f"{host} was unreachable earlier in this run: {_DOWN[host]}", EXIT_EXTERNAL)
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={**UA, **(headers or {})})
            with urllib.request.urlopen(req, timeout=90, context=_CTX) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            last = e
            if isinstance(e, urllib.error.HTTPError) and e.code in (400, 404):
                break                               # a definite answer, not an outage
            time.sleep(min(3 * (i + 1), 20))
    if not (isinstance(last, urllib.error.HTTPError) and last.code in (400, 404)):
        _DOWN[host] = str(last)[:120]
    raise LAError(f"could not retrieve {url.split('?')[0]}: {last}", EXIT_EXTERNAL)


def _cache(ws: Path, source: str, key: str, url: str, records: List[Dict[str, Any]], extra: Optional[Dict] = None,
           retrieved_at: Optional[str] = None) -> Path:
    k = hashlib.sha1(key.encode()).hexdigest()[:12]
    p = ws / "work" / "evidence" / f"pub_{source}_{k}.json"
    doc = {"source": source, "query": key, "url": url, "retrieved_at": retrieved_at or now_iso(), "records": records, **(extra or {})}
    write_json(p, doc)
    WRITTEN[str(p)] = sha256_file(p, limit=None)
    if source in ("gwas", "myvariant", "ensembl", "mr", "clinvar", "spdi") and not retrieved_at:
        _shared_put(f"{source}|{key}", doc)
    return p


def ledger(ws: Path) -> Dict[str, str]:
    """Cache files the harness wrote: the ledger saved in state plus files written in this process."""
    sp = ws / "state.json"
    led: Dict[str, str] = {}
    if sp.exists():
        try:
            led.update(load_json(sp).get("evidence_ledger") or {})
        except Exception:  # noqa: BLE001
            pass
    led.update({k: v for k, v in WRITTEN.items() if k.startswith(str(ws))})
    return led


def known_refs(ws: Path) -> Dict[str, Dict[str, Any]]:
    """Every citable public record retrieved in this workspace: ref -> record. A cache file counts only when its bytes
    match what the harness wrote (a hand-written or edited file is ignored)."""
    out: Dict[str, Dict[str, Any]] = {}
    led = ledger(ws)
    for p in (ws / "work" / "evidence").glob("pub_*.json"):
        if led.get(str(p)) != sha256_file(p, limit=None):
            continue
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


# ------------------------------------------------------------------ ClinVar (NCBI) + canonical SPDI
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
VARIATION = "https://api.ncbi.nlm.nih.gov/variation/v0/"
GRCH38_ACC = {"1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12", "5": "NC_000005.10",
              "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12", "10": "NC_000010.11",
              "11": "NC_000011.10", "12": "NC_000012.12", "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
              "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
              "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10"}


_NCBI_LAST = [0.0]


def _ncbi(url: str, data: Optional[bytes] = None, need: str = "result") -> Dict[str, Any]:
    """NCBI E-utilities at ≤3 requests/s; an error payload (rate limit, backend) is retried, then raised as LAError."""
    last = None
    for i in range(4):
        wait = 0.4 - (time.time() - _NCBI_LAST[0])
        if wait > 0:
            time.sleep(wait)
        _NCBI_LAST[0] = time.time()
        body = _http(url, data=data, headers={"Content-Type": "application/x-www-form-urlencoded"} if data else None, tries=3)
        if isinstance(body, dict) and need in body:
            return body
        last = str(body)[:400] if not isinstance(body, dict) else str(body.get("error") or body)[:400]
        if "cannot be transformed" in last:              # deterministic: the batch is too large, retrying won't help
            break
        time.sleep(2 * (i + 1))
    raise LAError(f"NCBI E-utilities returned no {need}: {last}", EXIT_EXTERNAL)


def clinvar_gene(ws: Path, gene: str) -> List[Dict[str, Any]]:
    """ClinVar germline records in a gene whose classification mentions pathogenic (incl. conflicting ones), with their
    canonical SPDI (GRCh38) and GRCh38/GRCh37 locations. Structural variants without an SPDI are left out."""
    key = f"clinvar:gene:{gene}"
    hit = _shared_get(f"clinvar|{key}")
    if hit:
        _cache(ws, "clinvar", key, hit["url"], hit["records"], retrieved_at=hit["retrieved_at"])
        return hit["records"]
    term = f"{gene}[gene] AND (clinsig_pathogenic[prop] OR clinsig_likely_pathogenic[prop] OR clinsig_has_conflicts[prop])"
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode({"db": "clinvar", "term": term, "retmode": "json", "retmax": 10000})
    ids = _ncbi(url, need="esearchresult")["esearchresult"].get("idlist", [])
    recs: List[Dict[str, Any]] = []
    oversized: List[str] = []

    def summaries(batch: List[str]) -> List[Dict[str, Any]]:
        body = urllib.parse.urlencode({"db": "clinvar", "id": ",".join(batch), "retmode": "json"}).encode()
        try:
            res = _ncbi(EUTILS + "esummary.fcgi", data=body)["result"]
        except LAError as e:
            if "max size" not in str(e) and "cannot be transformed" not in str(e):
                raise
            if len(batch) == 1:                          # one huge record (a multi-gene structural variant)
                oversized.append(batch[0])
                return []
            h = len(batch) // 2
            return summaries(batch[:h]) + summaries(batch[h:])
        return [res[u] for u in res.get("uids", [])]

    for i in range(0, len(ids), 300):
        for x in summaries(ids[i:i + 300]):
            vs = x.get("variation_set") or []
            if len(vs) != 1 or not vs[0].get("canonical_spdi"):
                continue                                 # haplotypes / structural variants: not matchable to one VCF allele
            spdi = vs[0]["canonical_spdi"]
            try:
                acc, pos, dele, ins = spdi.split(":")
            except ValueError:
                continue
            gc = x.get("germline_classification") or {}
            locs = {l.get("assembly_name"): {"chr": str(l.get("chr")), "start": int(l["start"]), "stop": int(l["stop"])}
                    for l in vs[0].get("variation_loc", []) if l.get("start") and str(l.get("start")).isdigit()}
            recs.append({"ref": f"clinvar:{x.get('accession')}", "gene": gene, "title": x.get("title"), "spdi": spdi,
                         "deleted": dele, "inserted": ins, "loc": locs, "classification": gc.get("description") or "",
                         "review_status": gc.get("review_status") or "",
                         "conditions": sorted({t.get("trait_name") for t in (gc.get("trait_set") or []) if t.get("trait_name")})[:5]})
    _cache(ws, "clinvar", key, EUTILS + "esearch.fcgi (" + term + ")", recs, extra={"oversized_skipped": oversized})
    return recs


def canonical_spdi(ws: Path, assembly: str, chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Any]:
    """The member's VCF allele as NCBI canonical SPDI on GRCh38 (repeat-aware, so dup/ins/delins/unaligned calls of the
    same change compare equal). A reference-base mismatch is returned as a warning (wrong assembly or corrupt call)."""
    key = f"spdi:{assembly}:{chrom}:{pos}:{ref}:{alt}"
    hit = _shared_get(f"spdi|{key}")
    if hit:
        _cache(ws, "spdi", key, hit["url"], hit["records"], retrieved_at=hit["retrieved_at"])
        return hit["records"][0]
    gcf = "GCF_000001405.40" if assembly == "hg38" else "GCF_000001405.25"
    url = VARIATION + f"vcf/chr{chrom}/{pos}/{ref}/{alt}/contextuals?assembly={gcf}"
    ctx = (_http(url, tries=3).get("data") or {}).get("spdis") or []
    rec = {"query": key, "spdi": None, "warning": None}
    if ctx:
        c = ctx[0]
        cs = f"{c['seq_id']}:{c['position']}:{c['deleted_sequence']}:{c['inserted_sequence']}"
        if assembly != "hg38":
            eq = (_http(VARIATION + f"spdi/{urllib.parse.quote(cs)}/all_equivalent_contextual", tries=3).get("data") or {}).get("spdis") or []
            g38 = [e for e in eq if e.get("seq_id") in GRCH38_ACC.values()]
            cs = f"{g38[0]['seq_id']}:{g38[0]['position']}:{g38[0]['deleted_sequence']}:{g38[0]['inserted_sequence']}" if g38 else None
        if cs:
            d = _http(VARIATION + f"spdi/{urllib.parse.quote(cs)}/canonical_representative", tries=3).get("data") or {}
            if d.get("warnings"):
                rec["warning"] = "; ".join(w.get("message", "") for w in d["warnings"])[:200]
            if d.get("seq_id"):
                rec["spdi"] = f"{d['seq_id']}:{d['position']}:{d['deleted_sequence']}:{d['inserted_sequence']}"
    _cache(ws, "spdi", key, url, [dict(rec, ref=f"spdi:{key}")])
    return rec


# ------------------------------------------------------------------ Ensembl
def gene_region(ws: Path, symbol: str, assembly: str = "GRCh38") -> Optional[Dict[str, Any]]:
    host = "https://rest.ensembl.org" if assembly == "GRCh38" else "https://grch37.rest.ensembl.org"
    url = f"{host}/lookup/symbol/homo_sapiens/{urllib.parse.quote(symbol)}?content-type=application/json"
    hit = _shared_get(f"ensembl|gene:{assembly}:{symbol}")
    if hit:
        _cache(ws, "ensembl", f"gene:{assembly}:{symbol}", url, hit["records"], retrieved_at=hit["retrieved_at"])
        return hit["records"][0]
    try:
        d = _http(url, tries=3)
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
