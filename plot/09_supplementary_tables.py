#!/usr/bin/env python3
"""
09_supplementary_tables.py
──────────────────────────
Manuscript-ready supplementary tables built from the 8 scatter_output TSVs:

  Table 1  per-phenotype SAIGE-Tractor-CCT-unique discoveries vs All by All
           meta-analysis  →  table1_cct_unique.{tsv,xlsx}

  Stable 4 per-phenotype SAIGE-Tractor-unique discoveries through Heterogeneous
           (HET) and Homogeneous (HOM) combiners vs All by All meta-analysis
           (i.e. loci flagged by HOM or HET but not seen in CCT-vs-META or
           seen but classified differently)
           → stable4_het_hom_unique.{tsv,xlsx}

  Stable 5 per-phenotype SAIGE-Tractor-unique ancestry-specific discoveries vs
           All by All within-global-ancestry GWAS
           (loci unique in any per-ancestry scatter table — anc1..anc5 —
            i.e. NOT seen in the corresponding ancestry's All by All GWAS)
           → stable5_per_ancestry_unique.{tsv,xlsx}

Only essential columns are kept for readability. Numeric p-values are written as
scientific-notation strings (1.23e-08) for clean Excel display.

All output goes to `manuscript_tables/`.
"""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
from pathlib import Path

from common import (
    ROOT, ANCS, PHENO_LABELS,
    SCATTER_TABLES, load_scatter, safe_float,
    NAME_SAIGE, NAME_ABA,
)

TABLES_DIR = ROOT / "manuscript_tables"
TABLES_DIR.mkdir(exist_ok=True)

def pfmt(p):
    """Format p-values: scientific notation for small p, fixed decimals otherwise."""
    v = safe_float(p)
    if v is None: return ""
    if v == 0: return "<1e-320"
    if v < 0.001: return f"{v:.2e}"
    return f"{v:.3g}"

def bfmt(b, d=4):
    v = safe_float(b)
    return "" if v is None else f"{v:+.{d}f}"

def afmt(a, d=4):
    v = safe_float(a)
    return "" if v is None else f"{v:.{d}f}"

def gene_short(g):
    if not isinstance(g, str): return ""
    return g.split(";")[0]

def pheno_label(p):
    p = str(p)
    return PHENO_LABELS.get(p, PHENO_LABELS.get(f"pheno_{p}", p))

# ═══════════════════════════════════════════════════════════════════════════
# TABLE 1 — per-phenotype CCT-vs-META unique discoveries
# ═══════════════════════════════════════════════════════════════════════════
def build_table1():
    df = load_scatter(SCATTER_TABLES["cct"])
    sub = df[(df["tophit_source"]=="SAIGE") & (df["locus_status"]=="SAIGE_only")].copy()
    sub["SAIGE_p_cct_c"] = pd.to_numeric(sub["SAIGE_P_cct_admixed_c"], errors="coerce")
    sub = sub.sort_values(["phenotype","SAIGE_p_cct_c"])

    rows = []
    for _, r in sub.iterrows():
        ph = str(r["phenotype"])
        rows.append({
            "Phenotype ID":       ph,
            "Phenotype":          pheno_label(ph),
            "Chr":                str(r.get("Chr","")),
            "Position":           str(r.get("Pos","")),
            "Ref / Alt":          f"{r.get('Ref','')} / {r.get('Alt','')}",
            "Gene":               gene_short(r.get("Gene","")),
            "Function":           r.get("Func",""),
            f"{NAME_SAIGE} CCT p-value (conditioned)":   pfmt(r.get("SAIGE_P_cct_admixed_c")),
            f"{NAME_SAIGE} BETA (all-ancestry, conditioned)": bfmt(r.get("SAIGE_BETA_c_ancALL")),
            f"{NAME_SAIGE} SE (all-ancestry, conditioned)":   bfmt(r.get("SAIGE_SE_c_ancALL")),
            f"{NAME_SAIGE} AF (all-ancestry)":            afmt(r.get("SAIGE_AF_Allele2_ancALL")),
            f"{NAME_ABA} meta-analysis p-value":          pfmt(r.get("ABA_pvalue")),
            f"{NAME_ABA} meta-analysis BETA":             bfmt(r.get("ABA_META_BETA")),
            f"{NAME_ABA} meta-analysis SE":               bfmt(r.get("ABA_META_SE")),
            f"{NAME_ABA} META AF":                        afmt(r.get("ABA_AF_Allele2")),
        })
    out = pd.DataFrame(rows)
    out.to_csv(TABLES_DIR / "table1_cct_unique.tsv", sep="\t", index=False)
    try:
        out.to_excel(TABLES_DIR / "table1_cct_unique.xlsx", index=False)
    except Exception as e:
        print(f"[09] (skipping xlsx for table1 — {e})")
    print(f"[09] Table 1: {len(out)} {NAME_SAIGE} CCT-unique loci")

# ═══════════════════════════════════════════════════════════════════════════
# STABLE 4 — per-phenotype HET/HOM-unique vs CCT (and vs META)
# ═══════════════════════════════════════════════════════════════════════════
def build_stable4():
    df_cct = load_scatter(SCATTER_TABLES["cct"])
    df_hom = load_scatter(SCATTER_TABLES["hom"])
    df_het = load_scatter(SCATTER_TABLES["het"])

    # Set of CCT-discovered loci to subtract (= what CCT already finds)
    def loci_keys(d, source):
        s = d[(d["tophit_source"]==source)]
        return set(s["phenotype"].astype(str) + "::"
                   + s["Chr"].astype(str) + ":" + s["Pos"].astype(str))
    cct_saige = loci_keys(df_cct, "SAIGE")

    # HOM-only and HET-only loci NOT in CCT
    def unique_to(combiner_df, combiner_name):
        sub = combiner_df[(combiner_df["tophit_source"]=="SAIGE")
                          & (combiner_df["locus_status"]=="SAIGE_only")].copy()
        sub["key"] = (sub["phenotype"].astype(str) + "::"
                      + sub["Chr"].astype(str) + ":" + sub["Pos"].astype(str))
        sub = sub[~sub["key"].isin(cct_saige)].copy()
        sub["combiner"] = combiner_name
        return sub

    hom_unique = unique_to(df_hom, "HOM")
    het_unique = unique_to(df_het, "HET")
    big = pd.concat([hom_unique, het_unique], ignore_index=True)
    big["pick_p"] = pd.to_numeric(big["SAIGE_pvalue"], errors="coerce")
    big = big.sort_values(["phenotype","combiner","pick_p"])

    rows = []
    for _, r in big.iterrows():
        ph = str(r["phenotype"])
        rows.append({
            "Phenotype ID":     ph,
            "Phenotype":        pheno_label(ph),
            "Combiner":         r["combiner"],
            "Chr":              str(r.get("Chr","")),
            "Position":         str(r.get("Pos","")),
            "Ref / Alt":        f"{r.get('Ref','')} / {r.get('Alt','')}",
            "Gene":             gene_short(r.get("Gene","")),
            f"{NAME_SAIGE} HOM p (conditioned)":   pfmt(r.get("SAIGE_P_hom_admixed_c")),
            f"{NAME_SAIGE} HET p (conditioned)":   pfmt(r.get("SAIGE_P_het_admixed_c")),
            f"{NAME_SAIGE} CCT p (conditioned)":   pfmt(r.get("SAIGE_P_cct_admixed_c")),
            f"{NAME_ABA} meta-analysis p":         pfmt(r.get("ABA_pvalue")),
            f"{NAME_SAIGE} β AFR (cond.)":         bfmt(r.get("SAIGE_BETA_c_anc1")),
            f"{NAME_SAIGE} β EAS (cond.)":         bfmt(r.get("SAIGE_BETA_c_anc2")),
            f"{NAME_SAIGE} β EUR (cond.)":         bfmt(r.get("SAIGE_BETA_c_anc3")),
            f"{NAME_SAIGE} β NatAm (cond.)":       bfmt(r.get("SAIGE_BETA_c_anc4")),
            f"{NAME_SAIGE} β SAS (cond.)":         bfmt(r.get("SAIGE_BETA_c_anc5")),
        })
    out = pd.DataFrame(rows)
    out.to_csv(TABLES_DIR / "stable4_het_hom_unique.tsv", sep="\t", index=False)
    try:
        out.to_excel(TABLES_DIR / "stable4_het_hom_unique.xlsx", index=False)
    except Exception as e:
        print(f"[09] (skipping xlsx for stable4 — {e})")
    print(f"[09] Stable 4: {len(out)} HOM/HET-unique loci (not in CCT)")

# ═══════════════════════════════════════════════════════════════════════════
# STABLE 5 — per-phenotype per-ancestry unique discoveries
# ═══════════════════════════════════════════════════════════════════════════
def build_stable5():
    """For each per-ancestry scatter table, list every SAIGE-Tractor unique
    locus and tag it with the local-ancestry test that flagged it."""
    rows = []
    anc_tables = [
        ("AFR",  "anc1", "anc_afr"),
        ("EAS",  "anc2", "anc_eas"),
        ("EUR",  "anc3", "anc_eur"),
        ("NatAm","anc4", "anc_amr"),
        ("SAS",  "anc5", "anc_sas"),
    ]
    for anc_name, suf, key in anc_tables:
        df = load_scatter(SCATTER_TABLES[key])
        sub = df[(df["tophit_source"]=="SAIGE") & (df["locus_status"]=="SAIGE_only")].copy()
        sub["pick_p"] = pd.to_numeric(sub["SAIGE_pvalue"], errors="coerce")
        sub = sub.sort_values(["phenotype","pick_p"])
        # The matching All-by-All per-ancestry column prefix
        aba_pref = "ABA_" + ("AMR" if anc_name=="NatAm" else anc_name)
        for _, r in sub.iterrows():
            ph = str(r["phenotype"])
            rows.append({
                "Phenotype ID":  ph,
                "Phenotype":     pheno_label(ph),
                "Ancestry":      anc_name,
                "Chr":           str(r.get("Chr","")),
                "Position":      str(r.get("Pos","")),
                "Ref / Alt":     f"{r.get('Ref','')} / {r.get('Alt','')}",
                "Gene":          gene_short(r.get("Gene","")),
                f"{NAME_SAIGE} {anc_name} p (conditioned)":
                                 pfmt(r.get(f"SAIGE_p.value_c_{suf}")),
                f"{NAME_SAIGE} {anc_name} BETA (conditioned)":
                                 bfmt(r.get(f"SAIGE_BETA_c_{suf}")),
                f"{NAME_SAIGE} {anc_name} SE (conditioned)":
                                 bfmt(r.get(f"SAIGE_SE_c_{suf}")),
                f"{NAME_SAIGE} {anc_name} AF (local)":
                                 afmt(r.get(f"SAIGE_AF_Allele2_{suf}")),
                f"{NAME_SAIGE} {anc_name} N$_{{haplo}}$":
                                 int(safe_float(r.get(f"SAIGE_N_haplo_{suf}")) or 0) or "",
                f"{NAME_ABA} {anc_name} p":
                                 pfmt(r.get(f"{aba_pref}_Pvalue")),
                f"{NAME_ABA} {anc_name} BETA":
                                 bfmt(r.get(f"{aba_pref}_BETA")),
                f"{NAME_ABA} {anc_name} SE":
                                 bfmt(r.get(f"{aba_pref}_SE")),
                f"{NAME_ABA} {anc_name} AF (global)":
                                 afmt(r.get(f"{aba_pref}_AF_Allele2")),
            })
    out = pd.DataFrame(rows)
    out.to_csv(TABLES_DIR / "stable5_per_ancestry_unique.tsv", sep="\t", index=False)
    try:
        out.to_excel(TABLES_DIR / "stable5_per_ancestry_unique.xlsx", index=False)
    except Exception as e:
        print(f"[09] (skipping xlsx for stable5 — {e})")
    print(f"[09] Stable 5: {len(out)} per-ancestry SAIGE-Tractor unique loci")

if __name__ == "__main__":
    build_table1()
    build_stable4()
    build_stable5()
    print(f"[09] all tables written to {TABLES_DIR}")
