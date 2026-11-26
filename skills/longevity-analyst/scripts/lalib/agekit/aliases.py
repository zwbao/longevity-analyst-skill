"""English lab-report names and units for the measurements the library's methods take.

Skills declare their own aliases (mostly the names on Chinese reports). English reports, NHANES-style labels and
LOINC-style long names often differ ("Serum glucose", "C-reactive protein", "Lymphocyte percent",
"White blood cell count (1000 cells/uL)"). Before a row is matched against a method's inputs, the engine tries the
row's own name; when the method does not know it, it looks the name up here and uses the first equivalent name the
method does know. Units are normalised the same way. Structural data, not judgement: a name maps to one analyte.

Unit conventions follow NHANES as documented in biolearn's clinical registry (BSD-3-Clause,
https://github.com/bio-learn/biolearn/blob/master/biolearn/clinical/registry.py).
"""

from __future__ import annotations

import re
import unicodedata
from typing import Dict, List, Optional

# analyte -> (names seen on English reports, names the library's skills use for it)
ANALYTES: Dict[str, Dict[str, List[str]]] = {
    "albumin": {"en": ["albumin", "serum albumin", "alb"], "lib": ["albumin", "alb", "白蛋白", "血清白蛋白", "albumin_g_l", "albumin_gl"]},
    "alkaline_phosphatase": {"en": ["alkaline phosphatase", "alk phos", "alp", "alkp"], "lib": ["alp", "alkaline phosphatase", "akp", "碱性磷酸酶", "alp_u_l"]},
    "creatinine": {"en": ["creatinine", "serum creatinine", "creat", "scr", "cr"], "lib": ["creatinine", "scr", "cr", "crea", "肌酐", "血肌酐", "creatinine_umol_l", "creat_umol"]},
    "glucose": {"en": ["serum glucose", "glucose", "fasting glucose", "fasting plasma glucose", "glucose fasting", "plasma glucose", "fbg", "fpg"],
                "lib": ["glucose", "glu", "fpg", "葡萄糖", "空腹血糖", "血糖", "glucose_mmol_l", "glucose_mmol"]},
    "crp": {"en": ["c-reactive protein", "c reactive protein", "crp", "hs-crp", "hscrp", "high-sensitivity c-reactive protein", "high sensitivity crp"],
            "lib": ["crp", "hs-crp", "hscrp", "c-反应蛋白", "超敏c反应蛋白", "crp_mg_dl"]},
    "lymphocyte_pct": {"en": ["lymphocyte percent", "lymphocytes percent", "lymphocyte percentage", "lymphocytes", "lymphs", "lymphocyte %", "lymph %"],
                       "lib": ["lymph%", "lymphocyte %", "ly%", "lymphocyte percent", "淋巴细胞%", "lymph_pct"]},
    "mcv": {"en": ["mean cell volume", "mean corpuscular volume", "mcv"], "lib": ["mcv", "mean cell volume", "平均红细胞体积", "mcv_fl"]},
    "rdw": {"en": ["red cell distribution width", "red blood cell distribution width", "rdw", "rdw-cv"], "lib": ["rdw", "rdw-cv", "red cell distribution width", "红细胞分布宽度", "rdw_pct"]},
    "wbc": {"en": ["white blood cell count", "white blood cells", "white cell count", "leukocytes", "leukocyte count", "wbc"],
            "lib": ["wbc", "white blood cell count", "leukocytes", "白细胞", "白细胞计数", "wbc_10e3", "wbc_10e9_l"]},
    "rbc": {"en": ["red blood cell count", "red blood cells", "erythrocytes", "rbc"], "lib": ["rbc", "red blood cell count", "红细胞", "红细胞计数", "rbc_10e12_l"]},
    "hemoglobin": {"en": ["hemoglobin", "haemoglobin", "hgb", "hb"], "lib": ["hgb", "hb", "hemoglobin", "血红蛋白", "hgb_g_dl"]},
    "hematocrit": {"en": ["hematocrit", "haematocrit", "hct"], "lib": ["hct", "hematocrit", "红细胞比容", "hct_pct"]},
    "platelets": {"en": ["platelet count", "platelets", "plt"], "lib": ["plt", "platelet count", "血小板", "血小板计数", "plt_10e9_l"]},
    "mch": {"en": ["mean cell hemoglobin", "mean corpuscular hemoglobin", "mch"], "lib": ["mch", "mean cell hemoglobin", "mch_pg"]},
    "mchc": {"en": ["mean cell hemoglobin concentration", "mean corpuscular hemoglobin concentration", "mchc", "mean cell hgb conc"], "lib": ["mchc", "mchc_g_dl"]},
    "neutrophil_pct": {"en": ["segmented neutrophils percent", "neutrophils percent", "neutrophil percent", "neutrophils", "neutrophil %"],
                       "lib": ["neut%", "ne%", "neutrophil percent", "中性粒细胞%", "neut_pct"]},
    "monocyte_pct": {"en": ["monocyte percent", "monocytes percent", "monocytes", "monocyte %"], "lib": ["mono%", "mo%", "monocyte percent", "单核细胞%", "mono_pct"]},
    "eosinophil_pct": {"en": ["eosinophils percent", "eosinophil percent", "eosinophils", "eosinophil %"], "lib": ["eo%", "eos%", "eosinophil percent", "嗜酸性粒细胞%", "eos_pct"]},
    "basophil_pct": {"en": ["basophils percent", "basophil percent", "basophils", "basophil %"], "lib": ["baso%", "ba%", "basophil percent", "嗜碱性粒细胞%", "baso_pct"]},
    "urea": {"en": ["blood urea nitrogen", "urea nitrogen", "bun", "urea"], "lib": ["bun", "urea", "blood urea nitrogen", "尿素氮", "urea_mmol_l"]},
    "uric_acid": {"en": ["uric acid", "urate", "ua"], "lib": ["uric acid", "ua", "尿酸", "uric_acid_umol_l"]},
    "total_cholesterol": {"en": ["total cholesterol", "cholesterol total", "cholesterol", "tc", "chol"],
                          "lib": ["tc", "chol", "total cholesterol", "总胆固醇", "胆固醇", "tc_mmol_l", "tc_mg_dl"]},
    "hdl": {"en": ["hdl", "hdl cholesterol", "hdl-c", "hdl-cholesterol", "high-density lipoprotein cholesterol", "direct hdl-cholesterol"],
            "lib": ["hdl", "hdl-c", "high-density lipoprotein", "高密度脂蛋白", "高密度脂蛋白胆固醇", "hdl_mmol_l", "hdl_mg_dl"]},
    "ldl": {"en": ["ldl", "ldl cholesterol", "ldl-c", "ldl-cholesterol", "low-density lipoprotein cholesterol"],
            "lib": ["ldl", "ldl-c", "低密度脂蛋白", "低密度脂蛋白胆固醇", "ldl_mmol_l"]},
    "triglycerides": {"en": ["triglycerides", "triglyceride", "tg", "trig"], "lib": ["tg", "triglycerides", "甘油三酯", "三酰甘油", "tg_mmol_l"]},
    "hba1c": {"en": ["glycohemoglobin", "hemoglobin a1c", "haemoglobin a1c", "hba1c", "a1c", "glycated hemoglobin", "glycated haemoglobin"],
              "lib": ["hba1c", "a1c", "glycohemoglobin", "糖化血红蛋白", "糖化血红蛋白a1c", "hba1c_pct"]},
    "total_bilirubin": {"en": ["total bilirubin", "bilirubin total", "bilirubin", "tbil"], "lib": ["tbil", "total bilirubin", "bilirubin", "总胆红素", "tbil_umol_l"]},
    "alt": {"en": ["alanine aminotransferase", "alt", "sgpt", "gpt"], "lib": ["alt", "gpt", "alanine aminotransferase", "谷丙转氨酶", "alt_u_l"]},
    "ast": {"en": ["aspartate aminotransferase", "ast", "sgot", "got"], "lib": ["ast", "got", "aspartate aminotransferase", "谷草转氨酶", "ast_u_l"]},
    "ggt": {"en": ["gamma glutamyl transferase", "gamma-glutamyl transferase", "ggt", "gamma gt"], "lib": ["ggt", "gamma glutamyl transferase", "谷氨酰转肽酶", "ggt_u_l"]},
    "ldh": {"en": ["lactate dehydrogenase", "ldh", "ld"], "lib": ["ldh", "ld", "lactate dehydrogenase", "乳酸脱氢酶", "ldh_u_l"]},
    "total_protein": {"en": ["total protein", "protein total", "tp"], "lib": ["tp", "total protein", "总蛋白", "total_protein_g_l"]},
    "globulin": {"en": ["globulin", "glob"], "lib": ["globulin", "glb", "glob", "球蛋白", "globulin_g_l"]},
    "calcium": {"en": ["calcium", "total calcium", "ca"], "lib": ["calcium", "ca", "钙", "总钙", "calcium_mmol_l"]},
    "phosphorus": {"en": ["phosphorus", "phosphate", "inorganic phosphorus"], "lib": ["phosphorus", "phos", "磷", "无机磷", "phosphorus_mmol_l"]},
    "sodium": {"en": ["sodium", "na"], "lib": ["sodium", "na", "钠", "sodium_mmol_l"]},
    "potassium": {"en": ["potassium", "k"], "lib": ["potassium", "k", "钾", "potassium_mmol_l"]},
    "chloride": {"en": ["chloride", "cl"], "lib": ["chloride", "cl", "氯", "chloride_mmol_l"]},
    "bicarbonate": {"en": ["bicarbonate", "total co2", "co2", "hco3"], "lib": ["bicarbonate", "hco3", "co2", "tco2", "bicarbonate_mmol_l"]},
    "cystatin_c": {"en": ["cystatin c", "cystatin-c", "cysc"], "lib": ["cystatin c", "cys-c", "cysc", "胱抑素c", "cystatin_c_mg_l"]},
    "ferritin": {"en": ["ferritin", "serum ferritin"], "lib": ["ferritin", "fer", "血清铁蛋白", "ferritin_ug_l"]},
    "iron": {"en": ["iron", "serum iron"], "lib": ["iron", "fe", "铁", "iron_umol_l"]},
    "tibc": {"en": ["total iron binding capacity", "tibc"], "lib": ["tibc", "total iron binding capacity", "tibc_umol_l"]},
    "transferrin_saturation": {"en": ["transferrin saturation", "tsat", "iron saturation"], "lib": ["tsat", "transferrin saturation", "铁饱和度", "tsat_pct"]},
    "folate": {"en": ["folate", "serum folate", "folic acid"], "lib": ["folate", "fol", "叶酸", "folate_nmol_l"]},
    "vitamin_b12": {"en": ["vitamin b12", "b12", "cobalamin"], "lib": ["vitamin b12", "b12", "vitb12", "维生素 b12", "b12_pmol_l"]},
    "ntprobnp": {"en": ["nt-probnp", "ntprobnp", "n-terminal pro-bnp", "n-terminal pro b-type natriuretic peptide"], "lib": ["nt-probnp", "ntprobnp", "ntprobnp_pg_ml"]},
    "fibrinogen": {"en": ["fibrinogen"], "lib": ["fibrinogen", "fib", "fbg", "fibrinogen_g_l"]},
    "cotinine": {"en": ["cotinine", "serum cotinine"], "lib": ["cotinine", "serum cotinine", "可替宁", "cotinine_ng_ml"]},
    "uacr": {"en": ["urine albumin creatinine ratio", "urine albumin-to-creatinine ratio", "albumin/creatinine ratio", "uacr", "acr"],
             "lib": ["uacr", "acr", "urine albumin-to-creatinine ratio", "尿白蛋白肌酐比", "uacr_mg_g", "acr_mg_g"]},
    "sbp": {"en": ["systolic blood pressure", "systolic bp", "sbp", "systolic", "systolic blood pres", "systolic: blood pres"], "lib": ["sbp", "systolic blood pressure", "收缩压", "高压", "sbp_mmhg"]},
    "dbp": {"en": ["diastolic blood pressure", "diastolic bp", "dbp", "diastolic", "diastolic blood pres", "diastolic: blood pres"], "lib": ["dbp", "diastolic blood pressure", "舒张压", "低压", "dbp_mmhg"]},
    "pulse": {"en": ["pulse", "heart rate", "resting heart rate", "pulse rate", "60 sec pulse", "60 sec. pulse"], "lib": ["pulse", "heart rate", "心率", "pulse_per_min"]},
    "bmi": {"en": ["bmi", "body mass index"], "lib": ["bmi", "body mass index", "身体质量指数", "bmi_kg_m2"]},
    "waist": {"en": ["waist circumference", "waist"], "lib": ["waist", "waist circumference", "腰围", "waist_cm"]},
    "height": {"en": ["height", "standing height", "body height"], "lib": ["height", "身高", "height_cm"]},
    "weight": {"en": ["weight", "body weight"], "lib": ["weight", "体重", "weight_kg"]},
}

# Units as English reports and NHANES print them -> the forms skillkit already understands.
UNITS: Dict[str, str] = {
    "1000 cells/ul": "10^3/uL", "1000 cells/µl": "10^3/uL", "thousand/ul": "10^3/uL", "x10^3/ul": "10^3/uL",
    "10^3 cells/ul": "10^3/uL", "k/µl": "10^3/uL", "10e3/ul": "10^3/uL", "x10e3/ul": "10^3/uL",
    "million cells/ul": "10^12/L", "million/ul": "10^12/L", "x10^6/ul": "10^12/L", "10^6/ul": "10^12/L", "10e6/ul": "10^12/L",
    "x10^9/l": "10^9/L", "x10^12/l": "10^12/L",
}


def _fold(text: str) -> str:
    s = unicodedata.normalize("NFKC", str(text)).casefold()
    s = re.sub(r"\s*\([^)]*\)\s*$", "", s.rstrip(" ."))  # drop a trailing "(unit)" or "(abbreviation)"
    return re.sub(r"[\s_\-·•:：,，/\\]+", "", s)


_EN_INDEX: Dict[str, str] = {}
for _analyte, _names in ANALYTES.items():
    for _name in _names["en"]:
        _EN_INDEX.setdefault(_fold(_name), _analyte)


_QUALIFIERS = {"serum", "plasma", "level", "levels", "s", "p", "the", "of", "in", "result", "frozen", "refrigerated", "refrig",
               "si", "direct", "1st", "2nd", "3rd", "rdg"}


def _tokens(text: str) -> str:
    """Word-order-free form: "Glucose, serum" and "Serum glucose" give the same key."""
    s = unicodedata.normalize("NFKC", str(text)).casefold()
    s = re.sub(r"\s*\([^)]*\)\s*$", "", s.rstrip(" ."))
    words = [w for w in re.split(r"[^0-9a-z%#]+", s) if w and w not in _QUALIFIERS]
    return " ".join(sorted(words))


_TOKEN_INDEX: Dict[str, str] = {}
_AMBIGUOUS = set()
for _analyte, _names in ANALYTES.items():
    for _name in _names["en"]:
        _key = _tokens(_name)
        if _key in _TOKEN_INDEX and _TOKEN_INDEX[_key] != _analyte:
            _AMBIGUOUS.add(_key)
        _TOKEN_INDEX.setdefault(_key, _analyte)
for _key in _AMBIGUOUS:
    del _TOKEN_INDEX[_key]


def analyte_of(name: str) -> Optional[str]:
    return _EN_INDEX.get(_fold(name)) or (_TOKEN_INDEX.get(_tokens(name)) if _tokens(name) else None)


def unit_of(unit: str) -> str:
    text = re.sub(r"\*\*", "^", str(unit or ""))   # kg/m**2 (NHANES codebook) -> kg/m^2
    text = re.sub(r"(?i)^iu/", "U/", text.strip())   # IU/L for enzymes is U/L
    key = re.sub(r"\s+", " ", unicodedata.normalize("NFKC", text)).casefold()
    return UNITS.get(key, text)


def translate_rows(rows: List[dict], known, normalize_unit=None) -> List[dict]:
    """Rows whose name the method does not know get the first equivalent name it does know.

    known(name) -> bool tells whether the method's own alias index matches a name (skillkit rules). A name in
    parentheses ("Glycohemoglobin (%)") is tried with and without it."""
    out = []
    for row in rows:
        row = dict(row)
        row["unit"] = unit_of(row.get("unit", ""))
        if not known(row["item"]):
            analyte = analyte_of(row["item"]) or analyte_of(re.sub(r"[(\[（【][^)\]）】]*[)\]）】]", " ", row["item"]))
            if analyte:
                for candidate in ANALYTES[analyte]["lib"]:
                    if known(candidate):
                        row["source_item"] = row["item"]
                        row["item"] = candidate
                        break
        out.append(row)
    return out
