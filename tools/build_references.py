"""Build skills/longevity-analyst/data/ref_population.json from public data (run once; the output is committed).

Inputs (downloaded by hand, all public):
  NHANES 2005-2010 XPT files (DEMO BIOPRO CRP CBC BMX BPX TCHOL HDL GHB VID, cycles D E F) - https://wwwn.cdc.gov/nchs/nhanes/
  GMHI 2020: 4347_final_relative_abundances.txt + Final_metadata_4347.csv - https://github.com/jaeyunsung/GMHI_2020
Usage: python3 tools/build_references.py <nhanes_dir> <gmhi_abundance_txt> <gmhi_metadata_csv>
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "skills" / "longevity-analyst" / "scripts"))

PCTS = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
AGE_BANDS = [(20, 29), (30, 39), (40, 49), (50, 59), (60, 69), (70, 79), (80, 150)]

# analyte key, label, NHANES file stem, column(s) per cycle, factor to SI unit, SI unit, names used to match lab rows
LABS = [
    ("glucose_nonfasting", "血糖（非空腹，全体受检）", "BIOPRO", {"*": "LBXSGL"}, 1 / 18.0, "mmol/L", []),
    ("hba1c", "糖化血红蛋白", "GHB", {"*": "LBXGH"}, 1.0, "%", ["糖化血红蛋白", "HbA1c"]),
    ("alt", "谷丙转氨酶", "BIOPRO", {"*": "LBXSATSI"}, 1.0, "U/L", ["谷丙转氨酶", "ALT"]),
    ("ast", "谷草转氨酶", "BIOPRO", {"*": "LBXSASSI"}, 1.0, "U/L", ["谷草转氨酶", "AST"]),
    ("ggt", "γ-谷氨酰转移酶", "BIOPRO", {"*": "LBXSGTSI"}, 1.0, "U/L", ["γ-谷氨酰转移酶", "GGT"]),
    ("alp", "碱性磷酸酶", "BIOPRO", {"*": "LBXSAPSI"}, 1.0, "U/L", ["碱性磷酸酶", "ALP"]),
    ("albumin", "白蛋白", "BIOPRO", {"*": "LBXSAL"}, 10.0, "g/L", ["白蛋白", "ALB"]),
    ("creatinine", "肌酐", "BIOPRO", {"*": "LBXSCR"}, 88.4, "umol/L", ["肌酐", "Cr", "CREA"]),
    ("urate", "尿酸", "BIOPRO", {"*": "LBXSUA"}, 59.48, "umol/L", ["尿酸", "UA"]),
    ("tc", "总胆固醇", "TCHOL", {"*": "LBXTC"}, 1 / 38.67, "mmol/L", ["总胆固醇", "TC"]),
    ("hdl", "高密度脂蛋白胆固醇", "HDL", {"*": "LBDHDD"}, 1 / 38.67, "mmol/L", ["高密度脂蛋白胆固醇", "HDL-C"]),
    ("tg_nonfasting", "甘油三酯（非空腹，全体受检）", "BIOPRO", {"*": "LBXSTR"}, 1 / 88.57, "mmol/L", []),
    ("crp", "C反应蛋白", "CRP", {"*": "LBXCRP"}, 10.0, "mg/L", ["超敏C反应蛋白", "C反应蛋白", "hs-CRP", "CRP"]),
    ("vitd", "25-羟基维生素D", "VID", {"D": "LBDVIDMS", "*": "LBXVIDMS"}, 1 / 2.496, "ng/mL", ["25-羟基维生素D", "25-OH-D", "维生素D"]),
    ("wbc", "白细胞计数", "CBC", {"*": "LBXWBCSI"}, 1.0, "10^9/L", ["白细胞计数", "WBC"]),
    ("lym_pct", "淋巴细胞百分比", "CBC", {"*": "LBXLYPCT"}, 1.0, "%", ["淋巴细胞百分比", "LYM%"]),
    ("mcv", "平均红细胞体积", "CBC", {"*": "LBXMCVSI"}, 1.0, "fL", ["平均红细胞体积", "MCV"]),
    ("rdw", "红细胞分布宽度", "CBC", {"*": "LBXRDW"}, 1.0, "%", ["红细胞分布宽度", "RDW-CV", "RDW"]),
    ("plt", "血小板计数", "CBC", {"*": "LBXPLTSI"}, 1.0, "10^9/L", ["血小板计数", "PLT"]),
    ("hgb", "血红蛋白", "CBC", {"*": "LBXHGB"}, 10.0, "g/L", ["血红蛋白", "HGB"]),
    ("bmi", "BMI", "BMX", {"*": "BMXBMI"}, 1.0, "kg/m2", ["BMI"]),
    ("waist", "腰围", "BMX", {"*": "BMXWAIST"}, 1.0, "cm", ["腰围"]),
    ("sbp", "收缩压", "BPX", {"*": ["BPXSY1", "BPXSY2"]}, 1.0, "mmHg", ["收缩压"]),
    # fasting morning subsample, weighted with its own WTSAF2YR
    ("ldl", "低密度脂蛋白胆固醇（空腹）", "TRIGLY", {"*": "LBDLDLSI"}, 1.0, "mmol/L", ["低密度脂蛋白胆固醇", "LDL-C"]),
    ("tg", "甘油三酯（空腹）", "TRIGLY", {"*": "LBDTRSI"}, 1.0, "mmol/L", ["甘油三酯", "TG"]),
    ("glucose", "空腹血糖", "GLU", {"*": "LBDGLUSI"}, 1.0, "mmol/L", ["空腹血糖", "GLU", "FPG"]),
    ("insulin", "空腹胰岛素", "GLU", {"*": "LBXIN"}, 1.0, "uIU/mL", ["空腹胰岛素", "INS"]),
]


def wpct(x, w, qs):
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = np.cumsum(w) / w.sum()
    return [float(np.interp(q / 100, c, x)) for q in qs]


def nhanes(dirp: Path):
    demo = []
    for s in ("D", "E", "F"):
        d = pd.read_sas(dirp / f"DEMO_{s}.xpt", format="xport")[["SEQN", "RIAGENDR", "RIDAGEYR", "WTMEC2YR"]]
        d["cycle"] = s
        demo.append(d)
    demo = pd.concat(demo)
    out = {}
    for key, label, stem, cols, fac, unit, names in LABS:
        frames = []
        for s in ("D", "E", "F"):
            f = dirp / f"{stem}_{s}.xpt"
            if not f.exists():
                continue
            t = pd.read_sas(f, format="xport")
            col = cols.get(s, cols.get("*"))
            if isinstance(col, list):
                have = [c for c in col if c in t.columns]
                if not have:
                    continue
                t["v"] = t[have].mean(axis=1)
            else:
                if col not in t.columns:
                    continue
                t["v"] = t[col]
            t["w_own"] = t["WTSAF2YR"] if "WTSAF2YR" in t.columns else np.nan
            frames.append(t[["SEQN", "v", "w_own"]])
        if not frames:
            continue
        m = demo.merge(pd.concat(frames), on="SEQN")
        fasting = m["w_own"].notna().any()
        m["w"] = m["w_own"] if fasting else m["WTMEC2YR"]
        m = m.dropna(subset=["v", "w"])
        m = m[m.w > 0]
        m["v"] = m.v * fac
        strata = {}
        for sex, code in (("male", 1), ("female", 2)):
            for lo, hi in AGE_BANDS:
                g = m[(m.RIAGENDR == code) & m.RIDAGEYR.between(lo, hi)]
                if len(g) < (60 if fasting else 100):
                    continue
                strata[f"{sex}:{lo}-{hi}"] = {"n": int(len(g)), "pct": [round(v, 4) for v in wpct(g.v.values, g.w.values, PCTS)],
                                              "sd": round(float(np.sqrt(np.cov(g.v.values, aweights=g.w.values))), 4)}
        out[key] = {"label_zh": label, "unit": unit, "names": names, "strata": strata}
    return out


def gmhi_ref(abund: Path, meta: Path):
    from lalib.native import gmhi_profile
    a = pd.read_csv(abund, sep="\t", index_col=0)
    md = pd.read_csv(meta, header=None, index_col=0, low_memory=False)
    md = md[~md.index.duplicated()].T
    vals = []
    for col in a.columns:
        sp = {k: float(v) for k, v in a[col].items() if v and v > 0}
        vals.append(gmhi_profile(sp)[0] if sp else np.nan)
    md = md.iloc[:len(vals)].copy()
    md["gmhi"] = vals
    ph = md["Phenotype"].astype(str)
    region = md["Geographical Region or Population"].astype(str)
    healthy = md[ph.str.lower() == "healthy"]["gmhi"].dropna().values
    nonhealthy = md[ph.str.lower() != "healthy"]["gmhi"].dropna().values
    ea = md[(ph.str.lower() == "healthy") & region.str.contains("China|Japan|Korea|Chinese|Hong Kong|Taiwan", case=False)]["gmhi"].dropna().values
    ref = {"source": "Gupta VK et al. Nat Commun 2020;11:4635 (GMHI), 4347 public stool metagenomes; GMHI recomputed here with the skill's GMHI.R-equivalent code",
           "doi": "10.1038/s41467-020-18476-8", "pcts": PCTS,
           "healthy": {"n": int(len(healthy)), "pct": [round(float(np.percentile(healthy, q)), 4) for q in PCTS]},
           "nonhealthy": {"n": int(len(nonhealthy)), "pct": [round(float(np.percentile(nonhealthy, q)), 4) for q in PCTS]}}
    if len(ea) >= 100:
        ref["healthy_east_asia"] = {"n": int(len(ea)), "pct": [round(float(np.percentile(ea, q)), 4) for q in PCTS]}
    return ref


def main():
    nh, ab, meta = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    ref = {"schema": "la-refpop/1",
           "note": "Population reference distributions. Labs: NHANES 2005-2010 US civilian non-institutionalised population, MEC-weighted percentiles by sex and age band; glucose and triglycerides there are non-fasting for most participants. Gut: GMHI in 4347 public metagenomes. A member's position is read off these tables; which table is appropriate for a Chinese member is stated in the report, not assumed.",
           "pcts": PCTS,
           "labs_nhanes": {"source": "NHANES 2005-2006, 2007-2008, 2009-2010 (CDC/NCHS), WTMEC2YR-weighted", "url": "https://wwwn.cdc.gov/nchs/nhanes/",
                           "analytes": nhanes(nh)},
           "gmhi": gmhi_ref(ab, meta)}
    out = ROOT / "skills" / "longevity-analyst" / "data" / "ref_population.json"
    out.write_text(json.dumps(ref, ensure_ascii=False, indent=1))
    print(out, "labs:", len(ref["labs_nhanes"]["analytes"]), "gmhi healthy n:", ref["gmhi"]["healthy"]["n"],
          "EA:", ref["gmhi"].get("healthy_east_asia", {}).get("n"))


if __name__ == "__main__":
    main()
