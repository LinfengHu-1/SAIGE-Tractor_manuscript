#!/usr/bin/env python3
"""
08_precision_improvement.py
───────────────────────────
At every ABA META genome-wide-significant (p < 5e-8) locus, compare the
standard error and chi-square statistic between SAIGE-Tractor and All-by-All
at the *same variant*.  We compute SE reduction and chi-sq increase as
percentages and report median/mean by phenotype × ancestry stratum.

Two comparisons are reported:

  (1) **ancALL-vs-META** — SAIGE-Tractor's ancALL-conditioned estimate vs the
      ABA fixed-effect meta. Compares full-cohort estimates between the two
      pipelines using the same total N.
  (2) **Per-ancestry** — SAIGE-Tractor's per-local-ancestry conditioned
      estimate (anc1..5) vs the ABA per-global-ancestry estimate
      (AFR/EAS/EUR/AMR/SAS). Compares like-for-like ancestry strata.

For each pair (variant, comparison), we compute:
   SE_reduction_pct = 100 × (SE_ABA - SE_SAIGE) / SE_ABA
   chisq_ABA        = (BETA_ABA / SE_ABA)^2
   chisq_SAIGE      = (BETA_SAIGE / SE_SAIGE)^2
   chisq_increase_pct = 100 × (chisq_SAIGE - chisq_ABA) / chisq_ABA

We require both SAIGE and ABA estimates to be non-null and BETA != 0.

Outputs
─────────
  replication_output/precision_improvement_loci.tsv     – per-locus per-comparison
  replication_output/precision_improvement_summary.tsv  – mean/median by stratum
  replication_output/precision_improvement_summary.json – machine-readable summary

A one-liner suitable for the manuscript appears at the end of stdout.
"""
from __future__ import annotations
import json, math
import numpy as np
import pandas as pd

from common import (
    REPLICATION_DIR, ANCS, PHENO_LABELS,
    SCATTER_TABLES, load_scatter, safe_float,
)

GW = 5e-8

# Load all three combiner-vs-META scatter tables (same rows underneath; we
# pool them to get the maximum coverage of ABA-significant variants since each
# combiner-vs-META table independently clumps ABA top hits).
df_cct = load_scatter(SCATTER_TABLES["cct"])
df_hom = load_scatter(SCATTER_TABLES["hom"])
df_het = load_scatter(SCATTER_TABLES["het"])

all_rows = pd.concat([df_cct, df_hom, df_het], ignore_index=True)

# Keep only rows that are ABA meta tophits OR shared loci (so we anchor on
# ABA-significant variants). Then filter to p < 5e-8 in ABA META.
anchor = all_rows[all_rows["tophit_source"]=="ABA"].copy()
anchor["ABA_pvalue_num"] = pd.to_numeric(anchor["ABA_pvalue"], errors="coerce")
anchor = anchor[anchor["ABA_pvalue_num"] < GW]

# Dedup by phenotype + chr + pos (a locus may appear in all 3 tables)
anchor["_key"] = (anchor["phenotype"].astype(str) + "::"
                  + anchor["Chr"].astype(str) + ":" + anchor["Pos"].astype(str))
anchor = anchor.drop_duplicates(subset=["_key"])
print(f"[08] {len(anchor)} ABA-significant loci (p < 5e-8) to evaluate")

# ─── Per-locus comparison rows ────────────────────────────────────────────────
def per_anc_meta_col(anc_name):
    return f"ABA_{anc_name}_BETA", f"ABA_{anc_name}_SE", f"ABA_{anc_name}_Pvalue"

ABA_ANC_NAMES = ["AFR","EAS","EUR","AMR","SAS"]
SAIGE_SUFFIX_FOR_ABA = {"AFR":"anc1","EAS":"anc2","EUR":"anc3","AMR":"anc4","SAS":"anc5"}

rows = []
for _, r in anchor.iterrows():
    ph = str(r["phenotype"])
    label = PHENO_LABELS.get(ph, ph)
    chrom = r["Chr"]; pos = r["Pos"]
    gene  = str(r.get("Gene","")).split(";")[0]

    # ── (1) ancALL-vs-META ──
    b_meta = safe_float(r.get("ABA_META_BETA"))
    s_meta = safe_float(r.get("ABA_META_SE"))
    b_sgt  = safe_float(r.get("SAIGE_BETA_c_ancALL"))
    s_sgt  = safe_float(r.get("SAIGE_SE_c_ancALL"))
    if (b_meta is not None and s_meta is not None and s_meta > 0
        and b_sgt is not None and s_sgt is not None and s_sgt > 0):
        chi_meta = (b_meta/s_meta)**2
        chi_sgt  = (b_sgt /s_sgt )**2
        rows.append(dict(
            phenotype=ph, label=label, Gene=gene, Chr=chrom, Pos=pos,
            comparison="ancALL-vs-META",
            ancestry="ALL",
            ABA_BETA=b_meta, ABA_SE=s_meta, ABA_chisq=chi_meta,
            SAIGE_BETA=b_sgt, SAIGE_SE=s_sgt, SAIGE_chisq=chi_sgt,
            SE_reduction_pct       = 100*(s_meta - s_sgt)/s_meta,
            chisq_increase_pct     = 100*(chi_sgt - chi_meta)/chi_meta,
            same_direction         = (b_meta*b_sgt) > 0,
            ABA_p = safe_float(r.get("ABA_pvalue")),
            SAIGE_p_cct = safe_float(r.get("SAIGE_P_cct_admixed_c")),
        ))

    # ── (2) Per-ancestry ──
    for anc_name in ABA_ANC_NAMES:
        bc, sc, pc = per_anc_meta_col(anc_name)
        b_aba = safe_float(r.get(bc)); s_aba = safe_float(r.get(sc))
        suf   = SAIGE_SUFFIX_FOR_ABA[anc_name]
        b_sg  = safe_float(r.get(f"SAIGE_BETA_c_{suf}"))
        s_sg  = safe_float(r.get(f"SAIGE_SE_c_{suf}"))
        if (b_aba is not None and s_aba is not None and s_aba > 0
            and b_sg is not None and s_sg is not None and s_sg > 0):
            chi_a = (b_aba/s_aba)**2
            chi_s = (b_sg /s_sg )**2
            rows.append(dict(
                phenotype=ph, label=label, Gene=gene, Chr=chrom, Pos=pos,
                comparison="per-ancestry",
                ancestry=anc_name,
                ABA_BETA=b_aba, ABA_SE=s_aba, ABA_chisq=chi_a,
                SAIGE_BETA=b_sg, SAIGE_SE=s_sg, SAIGE_chisq=chi_s,
                SE_reduction_pct       = 100*(s_aba - s_sg)/s_aba,
                chisq_increase_pct     = 100*(chi_s - chi_a)/chi_a,
                same_direction         = (b_aba*b_sg) > 0,
                ABA_p = safe_float(r.get(pc)),
                SAIGE_p_cct = safe_float(r.get("SAIGE_P_cct_admixed_c")),
            ))

per_locus = pd.DataFrame(rows)
out_loci = REPLICATION_DIR / "precision_improvement_loci.tsv"
per_locus.to_csv(out_loci, sep="\t", index=False)
print(f"[08] wrote {out_loci} ({len(per_locus)} (locus,ancestry) comparison rows)")

# ─── Summary by stratum ─────────────────────────────────────────────────────
def summarize(grp):
    return pd.Series(dict(
        n_loci                          = len(grp),
        median_SE_reduction_pct         = grp["SE_reduction_pct"].median(),
        mean_SE_reduction_pct           = grp["SE_reduction_pct"].mean(),
        pct_SE_reduced                  = 100*(grp["SE_reduction_pct"] > 0).mean(),
        median_chisq_increase_pct       = grp["chisq_increase_pct"].median(),
        mean_chisq_increase_pct         = grp["chisq_increase_pct"].mean(),
        pct_chisq_increased             = 100*(grp["chisq_increase_pct"] > 0).mean(),
        pct_same_direction              = 100*grp["same_direction"].mean(),
    ))

# Overall: ancALL-vs-META and per-ancestry separately
summary_rows = []
for comp, sub in per_locus.groupby("comparison"):
    summary_rows.append(dict(
        stratum=f"OVERALL ({comp})",
        comparison=comp, ancestry="all",
        **summarize(sub).to_dict()
    ))
    if comp == "per-ancestry":
        for anc, sub2 in sub.groupby("ancestry"):
            summary_rows.append(dict(
                stratum=f"per-ancestry: {anc}",
                comparison=comp, ancestry=anc,
                **summarize(sub2).to_dict()
            ))

# Per-phenotype rows (ancALL-vs-META view, most directly comparable to old-draft bullet)
ancall = per_locus[per_locus["comparison"]=="ancALL-vs-META"]
for ph, sub in ancall.groupby("phenotype"):
    summary_rows.append(dict(
        stratum=f"phenotype (ancALL-vs-META): {PHENO_LABELS.get(ph, ph)} [{ph}]",
        comparison="ancALL-vs-META", ancestry="all",
        **summarize(sub).to_dict()
    ))

summary = pd.DataFrame(summary_rows)
out_summary = REPLICATION_DIR / "precision_improvement_summary.tsv"
summary.to_csv(out_summary, sep="\t", index=False)
print(f"[08] wrote {out_summary} ({len(summary)} summary rows)")

# JSON summary
jsummary = {
    "thresholds": {"ABA_GW_p": GW},
    "n_loci_anchor": int(len(anchor)),
    "overall_ancALL_vs_META": summarize(
        per_locus[per_locus["comparison"]=="ancALL-vs-META"]).to_dict(),
    "overall_per_ancestry": summarize(
        per_locus[per_locus["comparison"]=="per-ancestry"]).to_dict(),
    "per_ancestry_breakdown": {
        anc: summarize(per_locus[(per_locus["comparison"]=="per-ancestry")
                                  & (per_locus["ancestry"]==anc)]).to_dict()
        for anc in ABA_ANC_NAMES
        if not per_locus[(per_locus["comparison"]=="per-ancestry")
                          & (per_locus["ancestry"]==anc)].empty
    },
}
out_json = REPLICATION_DIR / "precision_improvement_summary.json"
with open(out_json, "w") as f:
    json.dump(jsummary, f, indent=2, default=str)
print(f"[08] wrote {out_json}")

# ─── One-line headline for the manuscript ──────────────────────────────────
o = jsummary["overall_ancALL_vs_META"]
print()
print("="*78)
print("MANUSCRIPT ONE-LINER")
print("="*78)
print(f"At {o['n_loci']:,} ABA-genome-wide-significant loci, SAIGE-Tractor "
      f"(ancALL-conditioned) reduced the standard error by a median of "
      f"{o['median_SE_reduction_pct']:.1f}% (mean {o['mean_SE_reduction_pct']:.1f}%) "
      f"and increased the chi-square statistic by a median of "
      f"{o['median_chisq_increase_pct']:.0f}% (mean {o['mean_chisq_increase_pct']:.0f}%) "
      f"vs All-by-All META. SE was reduced for {o['pct_SE_reduced']:.0f}% "
      f"of loci; chi-square was increased for {o['pct_chisq_increased']:.0f}%; "
      f"effect direction was concordant for {o['pct_same_direction']:.0f}%.")
print()
o2 = jsummary["overall_per_ancestry"]
print(f"At per-ancestry comparisons (same locus × ancestry strata, "
      f"n={o2['n_loci']:,}), SAIGE-Tractor's local-ancestry-conditioned "
      f"estimate reduced SE by a median of {o2['median_SE_reduction_pct']:.1f}% "
      f"and increased chi-square by a median of {o2['median_chisq_increase_pct']:.0f}%. "
      f"Per-ancestry breakdown:")
for anc, s in jsummary["per_ancestry_breakdown"].items():
    print(f"   {anc:5s}  n={int(s['n_loci']):>4d}   "
          f"median SE↓ {s['median_SE_reduction_pct']:+6.1f}%   "
          f"median χ²↑ {s['median_chisq_increase_pct']:+6.0f}%   "
          f"χ²↑ in {s['pct_chisq_increased']:.0f}% of loci")
