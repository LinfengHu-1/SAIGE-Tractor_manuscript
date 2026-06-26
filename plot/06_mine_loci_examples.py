#!/usr/bin/env python3
"""
06_mine_loci_examples.py
────────────────────────
Mine the scatter_output/CCT scatter table for manuscript-quality examples
across each of the four "advantage mechanisms" of SAIGE-Tractor:

   (1) Genuine SAIGE-only discovery — ABA META at the same variant ≥ 5e-6
   (2) Cross-ancestry effect heterogeneity (HET combiner advantage)
       — opposing/divergent β across ancestries; HET p ≫ HOM p
   (3) Single-ancestry signal diluted in global meta
       — one ancestry p < 1e-4 with concordant direction; others null
   (4) Local-ancestry conditioning unmasking
       — marginal CCT p ≫ conditioned CCT p (we report Δ -log10)

Each candidate row is enriched with per-ancestry stats and a one-line
"why interesting" note.  Output written to:
    replication_output/loci_candidates.tsv
    replication_output/loci_candidates_top.md   (boss-review shortlist)
"""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
import pandas as pd

from common import (
    REPLICATION_DIR, ANCS, PHENO_LABELS,
    SCATTER_TABLES, load_scatter, safe_float,
)

# Now uses the UNION of SAIGE-only loci across all 8 scatter tables (built
# by 02_replication_analysis.py).  Falls back to CCT-only if union missing.
_unified = REPLICATION_DIR / "unified_saige_only_loci.tsv"
if _unified.exists():
    print(f"[06] using unified SAIGE-only set from {_unified.name}")
    df_union_saige = pd.read_csv(_unified, sep="\t", low_memory=False)
    df_union_saige["phenotype"] = df_union_saige["phenotype"].astype(str)
else:
    print("[06] no unified set found — falling back to CCT-only (run 02 first)")
    df_union_saige = pd.DataFrame()

df = load_scatter(SCATTER_TABLES["cct"])
df["phenotype"] = df["phenotype"].astype(str)

def label(ph): return PHENO_LABELS.get(ph, ph)

# Pre-compute per-row metrics
def neg(p):
    p = safe_float(p)
    return None if (p is None or p <= 0) else -math.log10(p)

def per_anc_stats(r):
    out = {}
    for suf, name, _, _ in ANCS:
        out[name] = dict(
            af  = safe_float(r.get(f"SAIGE_AF_Allele2_{suf}")),
            beta= safe_float(r.get(f"SAIGE_BETA_c_{suf}")),
            se  = safe_float(r.get(f"SAIGE_SE_c_{suf}")),
            p   = safe_float(r.get(f"SAIGE_p.value_c_{suf}")),
            nh  = safe_float(r.get(f"SAIGE_N_haplo_{suf}")),
        )
    return out

# Annotate each SAIGE_only row
rows = []
# Use UNION SAIGE-only loci (across CCT/HOM/HET + 5 per-ancestry tables) — was
# CCT-only previously, which silently dropped ~370 per-ancestry-specific hits
# (e.g. RCOR1 in NatAm).
if not df_union_saige.empty:
    saige_only = df_union_saige.copy()
    saige_only["tests_flagging"] = saige_only["tests_flagging"].astype(str)
else:
    saige_only = df[(df["tophit_source"]=="SAIGE") & (df["locus_status"]=="SAIGE_only")].copy()
    saige_only["tests_flagging"] = "CCT"
shared = df[(df["tophit_source"]=="SAIGE") & (df["locus_status"]=="shared")].copy()
shared["tests_flagging"] = "CCT"

for src_label, src in [("SAIGE_only", saige_only), ("shared", shared)]:
    for _, r in src.iterrows():
        ph = r["phenotype"]
        sp_cct = safe_float(r.get("SAIGE_P_cct_admixed_c"))
        sp_het = safe_float(r.get("SAIGE_P_het_admixed_c"))
        sp_hom = safe_float(r.get("SAIGE_P_hom_admixed_c"))
        sp_marg= safe_float(r.get("SAIGE_P_cct_admixed"))
        ap     = safe_float(r.get("ABA_pvalue"))
        ap_afr = safe_float(r.get("ABA_AFR_Pvalue"))
        ap_eur = safe_float(r.get("ABA_EUR_Pvalue"))
        ap_eas = safe_float(r.get("ABA_EAS_Pvalue"))
        ap_amr = safe_float(r.get("ABA_AMR_Pvalue"))
        ap_sas = safe_float(r.get("ABA_SAS_Pvalue"))
        pa = per_anc_stats(r)

        # Mechanism scores
        ng_log_aba = neg(ap)
        ng_log_cct = neg(sp_cct)
        adv_oom = (ng_log_cct - ng_log_aba) if (ng_log_cct and ng_log_aba) else None
        het_minus_hom = (neg(sp_het) - neg(sp_hom)) if (neg(sp_het) and neg(sp_hom)) else None
        la_boost = (neg(sp_cct) - neg(sp_marg)) if (neg(sp_cct) and neg(sp_marg)) else None

        # Direction: opposing / divergent
        betas = {n: pa[n]["beta"] for n in pa if pa[n]["beta"] is not None}
        n_pos = sum(b > 0 for b in betas.values()); n_neg = sum(b < 0 for b in betas.values())
        beta_range = (max(betas.values()) - min(betas.values())) if betas else None
        opposing = n_pos >= 1 and n_neg >= 1

        # Significant per-ancestry breakdown (which anc drives SAIGE)
        sig_anc = [n for n in pa if pa[n]["p"] is not None and pa[n]["p"] < 1e-4]
        sig_anc_strong = [n for n in pa if pa[n]["p"] is not None and pa[n]["p"] < 5e-8]

        # ABA per-ancestry availability
        aba_ran = sum(p is not None and not (isinstance(p, float) and math.isnan(p))
                      for p in [ap_afr, ap_eur, ap_eas, ap_amr, ap_sas])

        # Genuineness
        is_genuine = (ap is None) or (ap >= 5e-6)

        rows.append(dict(
            phenotype=ph, label=label(ph),
            tests_flagging=r.get("tests_flagging",""),
            Gene=str(r.get("Gene","")).split(";")[0],
            FullGene=r.get("Gene",""),
            Chr=r.get("Chr",""), Pos=int(r.get("Pos",0)),
            Ref=r.get("Ref",""), Alt=r.get("Alt",""),
            locus_status=src_label,
            SAIGE_p_cct_c=sp_cct, SAIGE_p_het_c=sp_het, SAIGE_p_hom_c=sp_hom,
            SAIGE_p_cct_marg=sp_marg,
            ABA_p=ap,
            ABA_AFR_p=ap_afr, ABA_EUR_p=ap_eur, ABA_EAS_p=ap_eas,
            ABA_AMR_p=ap_amr, ABA_SAS_p=ap_sas,
            ABA_anc_ran=aba_ran,
            sig_anc_str=",".join(sig_anc),
            sig_anc_gw_str=",".join(sig_anc_strong),
            beta_AFR  =pa["AFR"]["beta"],   p_AFR  =pa["AFR"]["p"],   AF_AFR  =pa["AFR"]["af"],   Nhap_AFR  =pa["AFR"]["nh"],
            beta_EAS  =pa["EAS"]["beta"],   p_EAS  =pa["EAS"]["p"],   AF_EAS  =pa["EAS"]["af"],   Nhap_EAS  =pa["EAS"]["nh"],
            beta_EUR  =pa["EUR"]["beta"],   p_EUR  =pa["EUR"]["p"],   AF_EUR  =pa["EUR"]["af"],   Nhap_EUR  =pa["EUR"]["nh"],
            beta_NatAm=pa["NatAm"]["beta"], p_NatAm=pa["NatAm"]["p"], AF_NatAm=pa["NatAm"]["af"], Nhap_NatAm=pa["NatAm"]["nh"],
            beta_SAS  =pa["SAS"]["beta"],   p_SAS  =pa["SAS"]["p"],   AF_SAS  =pa["SAS"]["af"],   Nhap_SAS  =pa["SAS"]["nh"],
            n_anc_pos=n_pos, n_anc_neg=n_neg, beta_range=beta_range, opposing=opposing,
            advantage_OOM=adv_oom,
            het_minus_hom_OOM=het_minus_hom,
            la_conditioning_boost_OOM=la_boost,
            n_sig_anc=len(sig_anc),
            is_genuine_miss=is_genuine,
        ))

cand = pd.DataFrame(rows)
out_full = REPLICATION_DIR / "loci_candidates.tsv"
cand.to_csv(out_full, sep="\t", index=False)
print(f"[06] wrote {out_full} ({len(cand)} rows)")

# ─── Shortlists ──────────────────────────────────────────────────────────────
def fmt_p(v): return "—" if v is None or (isinstance(v,float) and pd.isna(v)) else (f"{v:.2e}")
def fmt_b(v): return "—" if v is None or (isinstance(v,float) and pd.isna(v)) else (f"{v:+.3f}")
def fmt_af(v): return "—" if v is None or (isinstance(v,float) and pd.isna(v)) else (f"{v:.3f}")
def fmt_n(v): return "—" if v is None or (isinstance(v,float) and pd.isna(v)) else (f"{int(v):,}")

# Mechanism 1: genuine SAIGE-only, sorted by SAIGE strength + ABA weakness
m1 = cand[(cand["locus_status"]=="SAIGE_only") & cand["is_genuine_miss"]].copy()
m1["score"] = (-np.log10(m1["SAIGE_p_cct_c"].clip(lower=1e-320))) - (
    -np.log10(m1["ABA_p"].fillna(1).clip(lower=1e-320)))
m1 = m1.sort_values("score", ascending=False)

# Mechanism 2: cross-ancestry heterogeneity (HET-dominant)
m2 = cand[
    (cand["locus_status"]=="SAIGE_only")
    & cand["het_minus_hom_OOM"].notna()
    & (cand["het_minus_hom_OOM"] > 2)
    & cand["opposing"]
].copy().sort_values("het_minus_hom_OOM", ascending=False)

# Mechanism 3: single-ancestry signal diluted in META
def single_anc_score(r):
    if r["n_sig_anc"] != 1: return -1
    # Stronger ancestry-specific p, weaker META p → better
    ap = r["ABA_p"] if pd.notna(r["ABA_p"]) else 1
    return (-math.log10(max(min(r["SAIGE_p_cct_c"] or 1, 1), 1e-320))
            ) - (-math.log10(max(ap or 1, 1e-320)))
m3 = cand[(cand["locus_status"]=="SAIGE_only") & (cand["n_sig_anc"]==1)].copy()
m3["score"] = m3.apply(single_anc_score, axis=1)
m3 = m3.sort_values("score", ascending=False)

# Mechanism 4: LA-conditioning unmasking
m4 = cand[
    cand["la_conditioning_boost_OOM"].notna()
    & (cand["la_conditioning_boost_OOM"] > 5)
].copy().sort_values("la_conditioning_boost_OOM", ascending=False)

# ─── Boss-review shortlist (markdown) ────────────────────────────────────────
out_md = REPLICATION_DIR / "loci_candidates_top.md"
with open(out_md, "w") as f:
    def section(title, df_, ntop, cols_blurb):
        f.write(f"\n## {title}\n\n")
        f.write(f"**Top {ntop} candidates** (full ranked list in `loci_candidates.tsv`)\n\n")
        f.write("| Gene | Pheno | Chr:Pos | SAIGE CCT p | ABA META p | Flagged by | Key per-ancestry detail |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for _, r in df_.head(ntop).iterrows():
            per_anc = []
            for n in ["AFR","EAS","EUR","NatAm","SAS"]:
                b = r[f"beta_{n}"]; p = r[f"p_{n}"]
                if p is not None and not pd.isna(p) and p < 0.05:
                    per_anc.append(f"{n} β={fmt_b(b)} p={fmt_p(p)}")
            blurb = " · ".join(per_anc) if per_anc else "—"
            f.write(f"| **{r['Gene']}** | {r['label']} | {r['Chr']}:{r['Pos']} | {fmt_p(r['SAIGE_p_cct_c'])} | {fmt_p(r['ABA_p'])} | {r['tests_flagging']} | {blurb} |\n")

    f.write("# Manuscript loci shortlist (mined from current scatter_output)\n")
    f.write("\nGenerated by `scripts/06_mine_loci_examples.py` for boss review before final manuscript example selection.\n")
    f.write("\nAll loci are from the **CCT-vs-META** scatter table after local-ancestry conditioning. ")
    f.write("'Genuine miss' = ABA META p ≥ 5×10⁻⁶ or variant not tested in META. ")
    f.write("'Opposing' = at least one ancestry with β>0 and one with β<0 in SAIGE-Tractor.\n")

    section("Mechanism 1 — Genuine SAIGE-Tractor-only discoveries (most extreme advantage)",
            m1, 30, [])
    section("Mechanism 2 — Cross-ancestry effect heterogeneity (HET ≫ HOM, opposing directions)",
            m2, 30, [])
    section("Mechanism 3 — Single-ancestry signal diluted in global meta",
            m3, 30, [])
    section("Mechanism 4 — Local-ancestry conditioning unmasking (Δ -log10p > 5)",
            m4, 30, [])

    # Per-source breakdown: top 15 SAIGE-only hits per scatter table
    f.write("\n## Per-source breakdown — top 15 SAIGE-only loci flagged by each test\n\n")
    f.write("(The same locus may be flagged by multiple tests; the 'flagged by' column "
            "shows the full set. Loci flagged by a per-ancestry test but missed by CCT "
            "are the ones most often dropped by the old report.)\n")
    for src in ["CCT","HOM","HET","AFR","EAS","EUR","AMR","SAS"]:
        sub = cand[(cand["locus_status"]=="SAIGE_only")
                   & cand["tests_flagging"].fillna("").apply(
                        lambda s: src in [x.strip() for x in s.split(",")])
                   ].copy()
        if sub.empty: continue
        sub = sub.sort_values("SAIGE_p_cct_c").head(15)
        section(f"Flagged by **{src}** test (n={len(sub)})", sub, 15, [])

print(f"[06] wrote {out_md}")
print(f"[06] mechanism counts:  genuine_only={len(m1)}  het_opposing={len(m2)}  single_anc={len(m3)}  la_boost={len(m4)}")
