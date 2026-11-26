"""Age-direction vote for GTEx-style top-expressed gene lists (Host tool).

Uses curated literature gene sets (SenMayo / inflammaging / mito / ECM), not a
clock fitted on LongevityBench or GTEx eval. Presence + rank in the top-N list
drive votes; low feature count or near-tie → low_confidence.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

_SETS_REL = Path("data/reference/expression_age_gene_sets.json")
_CACHE: Dict[str, dict] = {}


def _load_sets(root: Path) -> dict:
    key = str(root)
    if key in _CACHE:
        return _CACHE[key]
    path = root / _SETS_REL
    if not path.is_file():
        raise FileNotFoundError(str(path))
    data = json.loads(path.read_text(encoding="utf-8"))
    _CACHE[key] = data
    return data


def parse_tissue(text: str) -> Optional[str]:
    m = re.search(
        r"two\s+([A-Za-z][A-Za-z \-/]+?)\s+samples?\s+comes?\s+from\s+an?\s+older",
        text,
        re.I,
    )
    if m:
        return m.group(1).strip()
    m = re.search(r"\b(Testis|Liver|Heart|Lung|Muscle|Blood|Nerve|Stomach|Breast|"
                  r"Prostate|Kidney|Skin|Adrenal|Pituitary|Vagina|Salivary|"
                  r"Blood Vessel|Adipose)\b", text, re.I)
    return m.group(1).title() if m else None


def _top_genes_block(block: str) -> List[str]:
    m = re.search(r"Top[- ]?\d+\s+expressed\s+genes\s*:\s*(.+)", block, re.I | re.S)
    if not m:
        return []
    raw = m.group(1)
    # stop before next Sample / GSEA
    raw = re.split(r"\n\s*(?:Sample\s+[bB]|GSEA|Transcriptomic)", raw, maxsplit=1)[0]
    genes = []
    for tok in re.split(r"[;\n]+", raw):
        g = tok.strip().split()[0] if tok.strip() else ""
        g = g.strip(" ,.")
        if g and re.match(r"^[A-Za-z][A-Za-z0-9_.\-]*$", g):
            genes.append(g.upper())
    return genes


def parse_top_lists(text: str) -> Tuple[List[str], List[str]]:
    """Return (genes_A, genes_B) ordered by decreasing expression rank."""
    parts = re.split(r"(?im)^\s*Sample\s+[bB]\s*:", text)
    if len(parts) < 2:
        # alternate: Sample a: ... Sample b:
        parts = re.split(r"(?im)^\s*Sample\s+[aA]\s*:", text)
        if len(parts) >= 3:
            return _top_genes_block(parts[1]), _top_genes_block(parts[2])
        return [], []
    # first part may include Sample a
    a_block = parts[0]
    if not re.search(r"Top[- ]?\d+\s+expressed", a_block, re.I):
        # try split on Sample a
        a_parts = re.split(r"(?im)^\s*Sample\s+[aA]\s*:", a_block)
        a_block = a_parts[-1] if a_parts else a_block
    return _top_genes_block(a_block), _top_genes_block(parts[1])


def _rank_map(genes: List[str]) -> Dict[str, int]:
    """1-based rank; first occurrence wins."""
    out: Dict[str, int] = {}
    for i, g in enumerate(genes, start=1):
        out.setdefault(g, i)
    return out


def _weight(sets: dict, gene: str, tissue: Optional[str]) -> Tuple[float, str]:
    """Signed weight: positive = rises with age."""
    rises = sets.get("rises_with_age") or {}
    falls = sets.get("falls_with_age") or {}
    w = 0.0
    direction = ""
    if gene in rises:
        w = float(rises[gene])
        direction = "rises"
    elif gene in falls:
        w = -float(falls[gene])
        direction = "falls"
    else:
        return 0.0, ""
    boosts = (sets.get("tissue_boost") or {}).get(tissue or "", {})
    if gene in boosts:
        # boost amplifies magnitude, keeps sign
        mag = abs(w) + float(boosts[gene])
        w = mag if w > 0 else -mag
    return w, direction


# Pathway name keywords → age direction for contrastive ssGSEA in the prompt.
# Positive weight: higher NES looks older. Negative: higher NES looks younger.
_PATHWAY_RISES = (
    "inflam", "interleukin", "complement", "interferon", "cytokine", "chemokine",
    "senescence", "nfkb", "nf-kb", "toll like", "tlr ", "innate immune",
    "matrix metalloprotein", "collagen", "extracellular matrix", "fibrosis",
    "wound healing", "platelet", "coagulation", "fibrin", "antigen present",
    "pd 1", "pd-1", "immune", "leukocyte", "neutrophil", "macrophage",
)
_PATHWAY_FALLS = (
    "mitochondri", "oxidative phosphoryl", "electron transport", "respiratory chain",
    "tca cycle", "citric acid", "fatty acid beta", "developmental lineage",
    "spermatogen", "sperm", "fertiliz", "meiosis", "germ cell",
    "dna replication", "dna strand elongation", "cell cycle mitotic",
    "ribosom", "translation initiation",
)

_GSEA_PAIR = re.compile(
    r"Sample-A shows\s+([\w\s\-]+?)\s+activity in\s+(.+?)\s*\(NES:\s*([+\-]?\d+(?:\.\d+)?)",
    re.I,
)
_GSEA_B = re.compile(
    r"while Sample-B shows\s+([\w\s\-]+?)\s+activity\s*\(NES:\s*([+\-]?\d+(?:\.\d+)?)",
    re.I,
)


def _pathway_weight(name: str) -> float:
    low = name.lower()
    for kw in _PATHWAY_RISES:
        if kw in low:
            return 1.2
    for kw in _PATHWAY_FALLS:
        if kw in low:
            return -1.0
    return 0.0


def parse_gsea_pairs(text: str) -> List[Dict[str, Any]]:
    """Parse contrastive ssGSEA lines already printed in LongevityBench prompts."""
    block_m = re.search(r"GSEA results:\s*(.+)", text, re.I | re.S)
    if not block_m:
        return []
    block = block_m.group(1)
    # Split on "Sample-A shows" keeping delimiters
    chunks = re.split(r"(?=Sample-A shows)", block)
    out: List[Dict[str, Any]] = []
    for ch in chunks:
        if "Sample-A shows" not in ch:
            continue
        ma = _GSEA_PAIR.search(ch)
        mb = _GSEA_B.search(ch)
        if not ma or not mb:
            continue
        pathway = ma.group(2).strip()
        nes_a, nes_b = float(ma.group(3)), float(mb.group(2))
        w = _pathway_weight(pathway)
        if w == 0.0 or nes_a == nes_b:
            continue
        if w > 0:
            points_to = "A" if nes_a > nes_b else "B"
        else:
            points_to = "A" if nes_a < nes_b else "B"
        out.append(
            {
                "pathway": pathway,
                "nes_A": nes_a,
                "nes_B": nes_b,
                "direction": "rises" if w > 0 else "falls",
                "points_to": points_to,
                "weight": round(abs(w) * abs(nes_a - nes_b), 4),
            }
        )
    return out


def _gene_vote_core(
    root: Path, text: str
) -> Tuple[Dict[str, Any], float, float, List[dict], int, int]:
    """Returns (meta_or_error, older_a, older_b, used, n_a, n_b)."""
    sets = _load_sets(root)
    tissue = parse_tissue(text)
    genes_a, genes_b = parse_top_lists(text)
    meta = {"tissue": tissue, "sets_sources": sets.get("sources") or []}
    if len(genes_a) < 10 or len(genes_b) < 10:
        return (
            {
                "status": "input_problem",
                "message": "need two Top-N expressed gene lists (Sample A/B)",
                "n_A": len(genes_a),
                "n_B": len(genes_b),
                **meta,
            },
            0.0,
            0.0,
            [],
            len(genes_a),
            len(genes_b),
        )
    ra, rb = _rank_map(genes_a), _rank_map(genes_b)
    n_a, n_b = len(genes_a), len(genes_b)
    older_a = older_b = 0.0
    used: List[dict] = []
    for gene in set(ra) | set(rb):
        w, direction = _weight(sets, gene, tissue)
        if w == 0.0:
            continue
        ea = (n_a + 1 - ra[gene]) / n_a if gene in ra else 0.0
        eb = (n_b + 1 - rb[gene]) / n_b if gene in rb else 0.0
        if ea == eb:
            continue
        if w > 0:
            points_to = "A" if ea > eb else "B"
            pts = abs(w) * abs(ea - eb)
        else:
            points_to = "A" if ea < eb else "B"
            pts = abs(w) * abs(ea - eb)
        if points_to == "A":
            older_a += pts
        else:
            older_b += pts
        used.append(
            {
                "feature": gene,
                "direction": direction,
                "rank_A": ra.get(gene),
                "rank_B": rb.get(gene),
                "points_to": points_to,
                "weight": round(abs(w), 3),
            }
        )
    used.sort(key=lambda x: -x["weight"])
    return ({"status": "ok", **meta, "min_features": int(sets.get("min_features") or 5),
             "share_abstain_band": float(sets.get("share_abstain_band") or 0.08)},
            older_a, older_b, used, n_a, n_b)


def compare_top_lists(
    root: Path,
    text: str,
) -> Dict[str, Any]:
    """Fuse literature gene-set vote with contrastive ssGSEA already in the prompt."""
    meta, older_a, older_b, used, n_a, n_b = _gene_vote_core(root, text)
    head = {
        "method": "age-direction vote (expression top-list + prompt ssGSEA)",
        "label_zh": "本库的方法：文献基因集 top-list 投票 + 题面 contrastive ssGSEA（不是 GTEx 拟合时钟）",
        "reference": "; ".join(meta.get("sets_sources") or []),
        "independent_check": "gene sets not fitted on LongevityBench; ssGSEA NES taken from the prompt text",
        "provided_by": "longevity-analyst age-evidence",
        "tissue": meta.get("tissue"),
    }
    if meta.get("status") == "input_problem":
        return {**head, **meta}

    gene_total = older_a + older_b
    gene_share = older_a / gene_total if gene_total else 0.5

    gsea = parse_gsea_pairs(text)
    gsea_a = sum(x["weight"] for x in gsea if x["points_to"] == "A")
    gsea_b = sum(x["weight"] for x in gsea if x["points_to"] == "B")
    gsea_total = gsea_a + gsea_b
    gsea_share = gsea_a / gsea_total if gsea_total else None

    # Fuse: prefer ssGSEA when it has aging-related pathways; else gene vote alone.
    if gsea_share is not None and len(gsea) >= 1:
        w_gsea, w_gene = 0.55, 0.45
        if len(used) < meta["min_features"]:
            w_gsea, w_gene = 0.75, 0.25
        share_a = w_gsea * gsea_share + w_gene * gene_share
        method_note = "fused"
    else:
        share_a = gene_share
        method_note = "gene_only"

    band = float(meta["share_abstain_band"])
    conflict = (
        gsea_share is not None
        and abs(gene_share - 0.5) >= band
        and abs(gsea_share - 0.5) >= band
        and (gene_share > 0.5) != (gsea_share > 0.5)
    )
    low = (
        (len(used) < meta["min_features"] and not gsea)
        or abs(share_a - 0.5) < band
        or conflict
    )
    looks = "A" if share_a > 0.5 else "B" if share_a < 0.5 else "tie"
    return {
        **head,
        "status": "ok",
        "outputs": {
            "looks_older": looks,
            "weighted_share_for_A_older": round(share_a, 3),
            "features_used": len(used),
            "features_given": n_a + n_b,
            "low_confidence": bool(low),
            "n_A": n_a,
            "n_B": n_b,
            "fusion": method_note,
            "gene_share_A": round(gene_share, 3),
            "gsea_share_A": round(gsea_share, 3) if gsea_share is not None else None,
            "gsea_pathways_used": len(gsea),
            "conflict_gene_vs_gsea": bool(conflict),
        },
        "strongest_features": used[:8],
        "gsea_votes": sorted(gsea, key=lambda x: -x["weight"])[:8],
        "caveats_zh": [
            "基因集来自文献先验；ssGSEA NES 来自题面，未在本库重算富集。",
            "top-N 列表只反映相对高表达，细胞组成混杂时把握要下调。",
            *(["基因票与 ssGSEA 冲突或份额接近 0.5：建议 Host 弃权 uncertain。"] if low else []),
        ],
    }
