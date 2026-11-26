"""Reference models fitted by this library on public data, with held-out validation (not published methods).

nhanes-age: how routine labs and body measurements change with age in US adults (NHANES 2011-2016, held out
2017-2018). For each sex and analyte the file stores the mean and SD of the transformed value at every age 18-80
and the correlation of residuals between analytes. Given any subset of a person's values, the engine scores every
age by a Gaussian likelihood and returns the posterior over age: a point estimate, a 90% interval and the
probability of each age band. It answers "what age do these values look like", for comparing people or reading
an unlabelled profile; the gap to the real age is NOT a validated biological age.

Standard library only (Cholesky on at most ~45 analytes).
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from . import aliases, panel, store
from .library import Library, LKError

MODEL_FILE = Path("data") / "reference" / "nhanes_age_reference.json"

# Factors from units seen on reports to the model's unit (value_in_model_unit = value * factor).
CONVERT: Dict[str, Dict[str, float]] = {
    "albumin": {"g/l": 1, "g/dl": 10},
    "alkaline_phosphatase": {"u/l": 1, "iu/l": 1},
    "urea": {"mg/dl": 1, "mmol/l": 2.801},  # model unit: BUN mg/dL; urea mmol/L x 2.801
    "creatinine": {"umol/l": 1, "mg/dl": 88.4},
    "glucose": {"mmol/l": 1, "mg/dl": 0.0555},
    "total_bilirubin": {"umol/l": 1, "mg/dl": 17.1},
    "uric_acid": {"umol/l": 1, "mg/dl": 59.48},
    "total_cholesterol": {"mmol/l": 1, "mg/dl": 0.02586},
    "triglycerides": {"mmol/l": 1, "mg/dl": 0.01129},
    "hdl": {"mmol/l": 1, "mg/dl": 0.02586},
    "ast": {"u/l": 1, "iu/l": 1}, "alt": {"u/l": 1, "iu/l": 1}, "ggt": {"u/l": 1, "iu/l": 1},
    "ldh": {"u/l": 1, "iu/l": 1},
    "total_protein": {"g/l": 1, "g/dl": 10}, "globulin": {"g/l": 1, "g/dl": 10},
    "sodium": {"mmol/l": 1, "meq/l": 1}, "potassium": {"mmol/l": 1, "meq/l": 1},
    "chloride": {"mmol/l": 1, "meq/l": 1}, "bicarbonate": {"mmol/l": 1, "meq/l": 1},
    "calcium": {"mmol/l": 1, "mg/dl": 0.2495}, "phosphorus": {"mmol/l": 1, "mg/dl": 0.3229},
    "wbc": {"10^9/l": 1, "10^3/ul": 1, "1000cells/ul": 1, "k/ul": 1},
    "rbc": {"10^12/l": 1, "10^6/ul": 1, "millioncells/ul": 1},
    "platelets": {"10^9/l": 1, "10^3/ul": 1, "1000cells/ul": 1},
    "lymphocyte_pct": {"%": 1}, "monocyte_pct": {"%": 1}, "neutrophil_pct": {"%": 1},
    "eosinophil_pct": {"%": 1}, "basophil_pct": {"%": 1}, "hematocrit": {"%": 1}, "rdw": {"%": 1},
    "hemoglobin": {"g/dl": 1, "g/l": 0.1}, "mchc": {"g/dl": 1, "g/l": 0.1},
    "mcv": {"fl": 1}, "mch": {"pg": 1}, "hba1c": {"%": 1},
    "crp": {"mg/l": 1, "mg/dl": 10},
    "bmi": {"kg/m^2": 1, "kg/m2": 1, "kg/㎡": 1},
    "waist": {"cm": 1, "m": 100}, "height": {"cm": 1, "m": 100}, "weight": {"kg": 1},
}
BANDS = [(18, 29), (30, 39), (40, 49), (50, 59), (60, 69), (70, 80)]


def _unit_key(unit: str) -> str:
    u = aliases.unit_of(unit or "")
    return "".join(str(u).split()).casefold().replace("µ", "u").replace("μ", "u").replace("**", "^")


def _index() -> Dict[str, str]:
    """folded name -> analyte, from both the English names and the library's names (Chinese report names)."""
    idx: Dict[str, str] = {}
    for analyte, names in aliases.ANALYTES.items():
        for n in names["lib"]:
            idx.setdefault(aliases._fold(n), analyte)
    return idx


_LIB_NAMES = _index()


def analyte_for(name: str) -> Optional[str]:
    return aliases.analyte_of(name) or _LIB_NAMES.get(aliases._fold(name))


def _transform(kind: str, x: float) -> Optional[float]:
    if kind == "none":
        return x
    if kind == "log":
        return math.log(x) if x > 0 else None
    if kind == "log1p":
        return math.log1p(x) if x >= 0 else None
    if kind == "logcrp":
        return math.log(max(x, 0.1))
    raise ValueError(kind)


def _cholesky(a: List[List[float]]) -> List[List[float]]:
    n = len(a)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = a[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                L[i][j] = math.sqrt(max(s, 1e-9))
            else:
                L[i][j] = s / L[j][j]
    return L


def _quad(L: List[List[float]], z: List[float]) -> float:
    """z' R^-1 z with R = L L'."""
    y = []
    for i in range(len(z)):
        y.append((z[i] - sum(L[i][k] * y[k] for k in range(i))) / L[i][i])
    return sum(v * v for v in y)


def load_model(lib: Library) -> Dict[str, Any]:
    path = lib.root / MODEL_FILE
    if not path.is_file():
        raise LKError(f"reference model {MODEL_FILE} is missing from this library", "library_broken")
    return json.loads(path.read_text(encoding="utf-8"))


def parse_bands(text: Any) -> Optional[List[Tuple[str, int, int]]]:
    """"18-30, 30-45, 75+" -> [(label, low, high)]; a shared edge belongs to the higher band (30 is in 30-45)."""
    if text in (None, "", []):
        return None
    parts = text if isinstance(text, list) else str(text).replace("，", ",").replace("；", ",").split(",")
    out = []
    for part in parts:
        p = str(part).strip()
        m = re.match(r"^(\d{1,3})\s*(?:-|–|~|to)\s*(\d{1,3})$", p) or re.match(r"^(\d{1,3})\s*\+$", p) or re.match(r"^<\s*(\d{1,3})$", p)
        if not m:
            raise LKError(f"band {p!r}: use 30-45, 75+ or <30", "bad_input", 3)
        if p.startswith("<"):
            out.append((p, 0, int(m.group(1)) - 1))
        elif p.endswith("+"):
            out.append((p, int(m.group(1)), 200))
        else:
            out.append((p, int(m.group(1)), int(m.group(2))))
    return out


def band_of(age: int, bands: List[Tuple[str, int, int]]) -> Optional[str]:
    """The band an age falls in; an edge shared by two bands (45 in 30-45 and 45-60) goes to the higher one."""
    hit = [(lo, lab) for lab, lo, hi in bands if lo <= age <= hi]
    return max(hit)[1] if hit else None


def estimate_age(lib: Library, rows: List[Dict[str, Any]], sex: Optional[str], bands: Any = None,
                 with_posterior: bool = False) -> Dict[str, Any]:
    meta = load_model(lib)
    head = {"method": "nhanes-age-reference",
            "label_zh": "自建参考模型（美国 NHANES 2011–2016 公开数据训练，2017–2018 留出验证；不是已发表方法）",
            "source": meta["fitted_on"], "held_out": meta["held_out"], "validation": meta["validation"],
            **lib.stamp()}
    if sex not in ("male", "female"):
        return {**head, "status": "input_problem", "problems": [{"key": "sex", "message_zh": "需要性别（男或女）：男女的参考曲线不同。"}]}
    m = meta["model"][sex]
    feats: List[str] = m["features"]
    used: Dict[str, Dict[str, Any]] = {}
    problems: List[Dict[str, str]] = []
    ignored: List[str] = []
    for r in rows:
        analyte = analyte_for(r["item"])
        if analyte not in feats:
            ignored.append(r["item"])
            continue
        try:
            value = float(str(r["value"]).replace(",", "."))
        except ValueError:
            problems.append({"key": analyte, "message_zh": f"{r['item']} 的值「{r['value']}」不是数。"})
            continue
        unit = _unit_key(r.get("unit", ""))
        table = CONVERT.get(analyte, {})
        factor = table.get(unit) if unit else (1.0 if len(set(table.values())) == 1 else None)
        if factor is None:
            problems.append({"key": analyte, "message_zh": f"{r['item']} 的单位「{r.get('unit') or '未写'}」不认识或缺失；可用："
                                                         f"{'、'.join(table)}。"})
            continue
        x = _transform(meta["features"][analyte]["transform"], value * factor)
        if x is None:
            problems.append({"key": analyte, "message_zh": f"{r['item']} 的值 {value} 不能用。"})
            continue
        lo = min(m["mean"][analyte][i] - 8 * m["sd"][analyte][i] for i in range(len(m["mean"][analyte])))
        hi = max(m["mean"][analyte][i] + 8 * m["sd"][analyte][i] for i in range(len(m["mean"][analyte])))
        if not lo <= x <= hi:
            problems.append({"key": analyte, "message_zh": f"{r['item']} {value} {r.get('unit', '')} 远超人群范围，多半是单位写错；这一项没有用。"})
            continue
        if analyte in used and abs(used[analyte]["x"] - x) > 1e-9:
            problems.append({"key": analyte, "message_zh": f"{r['item']} 出现两次且数值不同，这一项没有用。"})
            used.pop(analyte)
            continue
        used[analyte] = {"x": x, "from": f"{r['item']} {r['value']} {r.get('unit', '')}".strip()}
    if len(used) < 3:
        return {**head, "status": "input_problem", "problems": problems + [{"key": "", "message_zh": f"能用的项目只有 {len(used)} 个，至少要 3 个。"}],
                "ignored": ignored}
    keys = list(used)
    idx = [feats.index(k) for k in keys]
    R = [[m["corr"][i][j] for j in idx] for i in idx]
    L = _cholesky(R)
    ages = list(range(meta["ages"][0], meta["ages"][1] + 1))
    logp = []
    for i in range(len(ages)):
        z, logsd = [], 0.0
        for k in keys:
            sd = m["sd"][k][i]
            z.append(max(-4.0, min(4.0, (used[k]["x"] - m["mean"][k][i]) / sd)))
            logsd += math.log(sd)
        logp.append(-0.5 * _quad(L, z) - logsd)
    top = max(logp)
    p = [math.exp(v - top) for v in logp]
    total = sum(p)
    p = [v / total for v in p]
    cdf, acc = [], 0.0
    for v in p:
        acc += v
        cdf.append(acc)
    pick = lambda q: ages[next(i for i, c in enumerate(cdf) if c >= q)]  # noqa: E731
    estimate = sum(a * v for a, v in zip(ages, p))
    custom = parse_bands(bands)
    if custom:
        band_p: Dict[str, float] = {lab: 0.0 for lab, _lo, _hi in custom}
        for age, v in zip(ages, p):
            lab = band_of(age, custom)
            if lab:
                band_p[lab] += v
        bands_out = {k: round(v, 3) for k, v in band_p.items()}
    else:
        bands_out = {f"{a}-{b}" if b < 80 else f"{a}+": round(sum(v for age, v in zip(ages, p) if a <= age <= b), 3) for a, b in BANDS}
    extra = {"_posterior": {"ages": ages, "p": p}} if with_posterior else {}
    return {**head, **extra, "status": "ok",
            "outputs": {"age_estimate": round(estimate, 1), "interval_90": [pick(0.05), pick(0.95)],
                        "most_likely_age": ages[p.index(max(p))], "band_probability": bands_out,
                        "most_likely_band": max(bands_out, key=bands_out.get),
                        "n_values_used": len(keys)},
            "used": {k: used[k]["from"] for k in keys}, "ignored": ignored, "problems": problems,
            "caveats_zh": ["美国成人参考人群，没有在中国人群验证；80 表示 80 岁及以上。",
                           "这是“这些化验值像多少岁的人”，不是生物年龄：和实际年龄的差没有验证过与健康结局的关系。",
                           f"留出数据上（{meta['held_out']}，与 LongevityBench 同类的 21 项）平均误差约 "
                           f"{meta['validation'].get('longevitybench_like', {}).get('mae')} 岁，90% 区间覆盖约 "
                           f"{meta['validation'].get('longevitybench_like', {}).get('coverage_90')}。"]}


def estimate_from_payload(lib: Library, payload: Dict[str, Any]) -> Dict[str, Any]:
    payload = {**payload, "panel_text": store.resolve(payload)}
    if payload.get("panel_text") is not None:
        people = panel.parse_panel(str(payload["panel_text"]))
        if len(people) > 1:
            return {"status": "per_person", "people": [
                {"label": p["label"], **estimate_age(lib, p["rows"], payload.get("sex") or p["sex"], payload.get("bands"))}
                for p in people]}
        person = people[0]
        return estimate_age(lib, person["rows"], payload.get("sex") or person["sex"], payload.get("bands"))
    rows = payload.get("measurements") or []
    if not isinstance(rows, list):
        raise LKError("measurements must be a list of {item, value, unit}", "bad_input", 3)
    return estimate_age(lib, rows, payload.get("sex"), payload.get("bands"))
