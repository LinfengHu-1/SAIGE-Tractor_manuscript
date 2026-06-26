#!/usr/bin/env python3
"""
02_replication_analysis.py
──────────────────────────
Per-locus replication lookup for every SAIGE-Tractor-only and ABA-only
discovery, where "SAIGE-only" and "ABA-only" are taken as the **union**
across all 8 scatter tables (CCT/HOM/HET-vs-META plus the five
per-ancestry tables anc1..anc5 vs AFR/EAS/EUR/AMR/SAS).

  ┌─────────────────────────┬─────────────────────────────────────────────────────┐
  │ Source                   │ Lookup performed                                    │
  ├─────────────────────────┼─────────────────────────────────────────────────────┤
  │ SAIGE-Tractor unique     │ ABA META + per-ancestry sumstats: exact CHR:POS    │
  │   top hit (union)        │ + best ±500 kb window p-value & best direction     │
  │                          │                                                    │
  │ ABA unique top hit       │ SAIGE-T sumstat: exact CHR:POS + best ±500 kb      │
  │   (union)                │ window over P_cct_admixed_c (and per-ancestry      │
  │                          │ p.value_c_anc{i}) + concordance                    │
  │                          │                                                    │
  │ Shared top hit           │ Records the matched p-values already present in    │
  │                          │ the scatter table for a replication-rate panel     │
  └─────────────────────────┴─────────────────────────────────────────────────────┘

Outputs (under replication_output/)
  replication_saige_only_in_aba.tsv   – one row per SAIGE-only locus
  replication_aba_only_in_saige.tsv   – one row per ABA-only locus
  replication_shared_summary.tsv      – per-phenotype: shared, replicated/total
  replication_summary.json            – aggregate numbers for the report
  unified_saige_only_loci.tsv         – master union table of SAIGE-only hits
  unified_aba_only_loci.tsv           – master union table of ABA-only hits

Each output row carries a `tests_flagging` column listing which scatter tables
flagged that locus as SAIGE_only/ABA_only (e.g. "CCT,HOM,EUR" or "AMR" alone).
"""
from __future__ import annotations
import json, math
import numpy as np
import pandas as pd
from common import (
    REPLICATION_DIR, ANCS, PHENO_LABELS,
    SCATTER_TABLES, load_scatter,
    saige_to_canonical, canonical_to_aba_file_id,
    load_aba_sumstat, load_saige_chr,
    safe_float,
)

WINDOW       = 500_000
P_NOMINAL    = 0.05
P_SUGGESTIVE = 5e-6
P_GW         = 5e-8

SAIGE_KEEP = (
    ["CHR","POS","Allele1","Allele2",
     "BETA_c_ancALL","SE_c_ancALL","p.value_c_ancALL",
     "P_cct_admixed_c","P_het_admixed_c","P_hom_admixed_c",
     "AF_Allele2_ancALL"]
    + [c for s,_,_,_ in ANCS for c in (f"BETA_c_{s}", f"SE_c_{s}",
                                        f"p.value_c_{s}", f"AF_Allele2_{s}")]
)
ABA_KEEP = ["CHR","POS","Allele1","Allele2","BETA","SE","Pvalue","AF_Allele2"]

_SAIGE_CACHE: dict[tuple[str,str], pd.DataFrame] = {}
_ABA_CACHE  : dict[tuple[str,str], pd.DataFrame] = {}

def saige_chr(saige_pheno, chrom):
    key = (saige_pheno, str(chrom))
    if key not in _SAIGE_CACHE:
        try:
            df = load_saige_chr(saige_pheno, chrom, usecols=SAIGE_KEEP)
        except ValueError:
            df = load_saige_chr(saige_pheno, chrom)
        _SAIGE_CACHE[key] = df
    return _SAIGE_CACHE[key]

def aba_load(canonical_pheno, anc):
    key = (canonical_pheno, anc)
    if key not in _ABA_CACHE:
        try:
            df = load_aba_sumstat(canonical_pheno, anc)
            if not df.empty:
                df = df[[c for c in ABA_KEEP if c in df.columns]]
        except Exception:
            df = pd.DataFrame()
        _ABA_CACHE[key] = df
    return _ABA_CACHE[key]

def best_in_window(df, chrom, pos, p_col="Pvalue", beta_col="BETA"):
    if df is None or df.empty: return (None, None, None)
    chrom = str(chrom).replace("chr","")
    sub = df[(df["CHR"].astype(str)==chrom) &
             (df["POS"].between(pos-WINDOW, pos+WINDOW))]
    if sub.empty: return (None, None, None)
    sub = sub.dropna(subset=[p_col])
    if sub.empty: return (None, None, None)
    i = sub[p_col].idxmin()
    return (float(sub.loc[i, p_col]),
            int(sub.loc[i, "POS"]),
            None if beta_col not in sub.columns else safe_float(sub.loc[i, beta_col]))

def exact_lookup(df, chrom, pos, a1=None, a2=None, p_col="Pvalue", beta_col="BETA"):
    if df is None or df.empty: return (None, None)
    chrom = str(chrom).replace("chr","")
    sub = df[(df["CHR"].astype(str)==chrom) & (df["POS"]==pos)]
    if sub.empty: return (None, None)
    if a1 and a2 and "Allele1" in sub.columns and "Allele2" in sub.columns:
        m = ((sub["Allele1"]==a1) & (sub["Allele2"]==a2)) | \
            ((sub["Allele1"]==a2) & (sub["Allele2"]==a1))
        if m.any(): sub = sub[m]
    sub = sub.dropna(subset=[p_col])
    if sub.empty: return (None, None)
    r = sub.iloc[0]
    return (float(r[p_col]),
            None if beta_col not in sub.columns else safe_float(r[beta_col]))

def concord_label(beta_ref, beta_test):
    a = safe_float(beta_ref); b = safe_float(beta_test)
    if a is None or b is None: return "NA"
    if a*b > 0: return "concordant"
    if a*b < 0: return "discordant"
    return "NA"

def replication_call(p):
    if p is None: return "no_signal"
    if p < P_GW:        return "genome_wide"
    if p < P_SUGGESTIVE:return "suggestive"
    if p < P_NOMINAL:   return "nominal"
    return "no_signal"

# ─── Build the unified SAIGE-only / ABA-only enumeration ─────────────────────
# Iterate every scatter table and union rows tagged SAIGE_only / ABA_only.
# Dedup by (phenotype, Chr, Pos); keep the row with the smallest SAIGE p-value
# (for SAIGE-only) or smallest ABA p-value (for ABA-only).  Carry along a
# `tests_flagging` column listing which scatter tables tagged the locus.
TABLE_NAMES = ["CCT","HOM","HET","AFR","EAS","EUR","AMR","SAS"]
TABLE_KEYS  = ["cct","hom","het","anc_afr","anc_eas","anc_eur","anc_amr","anc_sas"]

print("[02] loading all 8 scatter tables …")
all_tabs = {nm: load_scatter(SCATTER_TABLES[k]) for nm, k in zip(TABLE_NAMES, TABLE_KEYS)}
for nm in TABLE_NAMES:
    all_tabs[nm]["phenotype"] = all_tabs[nm]["phenotype"].astype(str)
    all_tabs[nm]["__source_table__"] = nm

def unify(status: str, p_col_for_dedup: str):
    """status ∈ {SAIGE_only, ABA_only}; p_col is the metric used to pick best row."""
    frames = []
    for nm, t in all_tabs.items():
        side = "SAIGE" if status == "SAIGE_only" else "ABA"
        sub = t[(t["locus_status"]==status) & (t["tophit_source"]==side)].copy()
        frames.append(sub)
    big = pd.concat(frames, ignore_index=True)
    # Key by (phenotype, Chr, Pos). Some rows may have NA Pos — drop.
    big = big.dropna(subset=["Chr","Pos"])
    big["Pos"] = pd.to_numeric(big["Pos"], errors="coerce").astype("Int64")
    big = big.dropna(subset=["Pos"])
    big["_key"] = big["phenotype"].astype(str) + "::" + big["Chr"].astype(str) + ":" + big["Pos"].astype(str)
    # Build "tests_flagging" per key
    tests = big.groupby("_key")["__source_table__"].apply(
        lambda s: ",".join(sorted(set(s)))).to_dict()
    # Pick best row per key by p_col_for_dedup
    big["_p_pick"] = pd.to_numeric(big[p_col_for_dedup], errors="coerce")
    big = big.sort_values("_p_pick")
    uniq = big.drop_duplicates(subset=["_key"], keep="first").copy()
    uniq["tests_flagging"] = uniq["_key"].map(tests)
    uniq = uniq.drop(columns=["_key","_p_pick","__source_table__"])
    return uniq

saige_only = unify("SAIGE_only", "SAIGE_pvalue")
aba_only   = unify("ABA_only",   "ABA_pvalue")
print(f"[02] unified SAIGE-only loci: {len(saige_only)}  (was 158 in CCT-only)")
print(f"[02] unified ABA-only   loci: {len(aba_only)}  (was 46 in CCT-only)")

# Dump master union tables for transparency
saige_only.to_csv(REPLICATION_DIR / "unified_saige_only_loci.tsv", sep="\t", index=False)
aba_only.to_csv  (REPLICATION_DIR / "unified_aba_only_loci.tsv",   sep="\t", index=False)
print(f"[02] wrote unified_saige_only_loci.tsv and unified_aba_only_loci.tsv")

# Helper: count what each row sees in tests_flagging
def sources_summary(df, label):
    cnt = {}
    for s in df["tests_flagging"]:
        for t in str(s).split(","):
            cnt[t] = cnt.get(t,0)+1
    print(f"[02] {label}: per-table membership = {sorted(cnt.items(), key=lambda x: -x[1])}")
sources_summary(saige_only, "SAIGE-only union")
sources_summary(aba_only,   "ABA-only union")

# ─── SAIGE-only loci: replicate in ABA ────────────────────────────────────────
print(f"[02] {len(saige_only)} SAIGE-only loci to check in ABA …")
rows_s = []
for _, r in saige_only.iterrows():
    canon = saige_to_canonical(str(r["phenotype"]))
    chrom = str(r["Chr"]).replace("chr","")
    pos   = int(r["Pos"])
    a1, a2 = str(r.get("Ref","")), str(r.get("Alt",""))
    saige_p = safe_float(r.get("SAIGE_P_cct_admixed_c") or r.get("SAIGE_pvalue"))
    saige_beta_all = safe_float(r.get("SAIGE_BETA_c_ancALL"))

    rec = dict(phenotype=canon,
               label=PHENO_LABELS.get(canon, canon),
               tests_flagging=r["tests_flagging"],
               Gene=str(r.get("Gene","")).split(";")[0],
               Chr=chrom, Pos=pos, Allele1=a1, Allele2=a2,
               SAIGE_p=saige_p,
               SAIGE_beta_ancALL=saige_beta_all,
               SAIGE_p_cct_c=safe_float(r.get("SAIGE_P_cct_admixed_c")),
               SAIGE_p_het_c=safe_float(r.get("SAIGE_P_het_admixed_c")),
               SAIGE_p_hom_c=safe_float(r.get("SAIGE_P_hom_admixed_c")))

    df_meta = aba_load(canon, "META")
    p_meta_exact, b_meta_exact = exact_lookup(df_meta, chrom, pos, a1, a2)
    p_meta_win,   pos_meta_win, b_meta_win   = best_in_window(df_meta, chrom, pos)
    rec.update(ABA_META_p_exact=p_meta_exact,
               ABA_META_beta_exact=b_meta_exact,
               ABA_META_p_window=p_meta_win,
               ABA_META_pos_window=pos_meta_win,
               ABA_META_beta_window=b_meta_win,
               ABA_META_concord_exact=concord_label(saige_beta_all, b_meta_exact),
               ABA_META_concord_window=concord_label(saige_beta_all, b_meta_win),
               ABA_META_replication_class=replication_call(p_meta_win))

    for _suf, name, _aba_col_prefix, anc_csv in ANCS:
        anc_up = name.replace("NatAm","AMR")
        df_anc = aba_load(canon, anc_up)
        p_e, b_e = exact_lookup(df_anc, chrom, pos, a1, a2)
        p_w, pos_w, b_w = best_in_window(df_anc, chrom, pos)
        rec[f"ABA_{anc_up}_p_exact"]  = p_e
        rec[f"ABA_{anc_up}_p_window"] = p_w
        rec[f"ABA_{anc_up}_beta_window"] = b_w
        rec[f"ABA_{anc_up}_replication_class"] = replication_call(p_w)
    rows_s.append(rec)

rep_s = pd.DataFrame(rows_s).sort_values(["phenotype","Chr","Pos"])
out_s = REPLICATION_DIR / "replication_saige_only_in_aba.tsv"
rep_s.to_csv(out_s, sep="\t", index=False)
print(f"[02] wrote {out_s}")

# ─── ABA-only loci: replicate in SAIGE-T ──────────────────────────────────────
print(f"[02] {len(aba_only)} ABA-only loci to check in SAIGE-Tractor …")
rows_a = []
all_saige_files = {p.name.replace("merged_saigetractor_","").replace(".txt.gz","")
                   for p in (REPLICATION_DIR.parent / "saige_sumstat").iterdir()
                   if p.name.startswith("merged_saigetractor_")}

def saige_filename(canon):
    if canon in all_saige_files: return canon
    if f"pheno_{canon}" in all_saige_files: return f"pheno_{canon}"
    return None

def _skipped_row(r, canon, reason):
    return dict(phenotype=canon,
                label=PHENO_LABELS.get(canon, canon),
                tests_flagging=r.get("tests_flagging",""),
                Gene=str(r.get("Gene","")).split(";")[0],
                Chr=str(r["Chr"]).replace("chr",""),
                Pos=int(r["Pos"]),
                Allele1=str(r.get("Ref","")), Allele2=str(r.get("Alt","")),
                ABA_p=safe_float(r.get("ABA_pvalue")),
                ABA_beta=safe_float(r.get("ABA_META_BETA") or r.get("ABA_BETA")),
                SAIGE_cct_p_exact=None, SAIGE_cct_p_window=None,
                SAIGE_cct_replication_class=reason)

for canon, grp in aba_only.groupby(aba_only["phenotype"].astype(str).map(saige_to_canonical)):
    sf = saige_filename(canon)
    if sf is None:
        for _, r in grp.iterrows():
            rows_a.append(_skipped_row(r, canon, "saige_sumstat_missing"))
        continue
    for _, r in grp.iterrows():
        chrom = str(r["Chr"]).replace("chr","")
        pos   = int(r["Pos"])
        a1, a2 = str(r.get("Ref","")), str(r.get("Alt",""))
        aba_p = safe_float(r.get("ABA_pvalue"))
        aba_beta = safe_float(r.get("ABA_META_BETA") or r.get("ABA_BETA"))

        df_sg = saige_chr(sf, chrom)
        p_exact, b_exact = exact_lookup(df_sg, chrom, pos, a1, a2,
                                        p_col="P_cct_admixed_c",
                                        beta_col="BETA_c_ancALL")
        p_win, pos_win, b_win = best_in_window(df_sg, chrom, pos,
                                               p_col="P_cct_admixed_c",
                                               beta_col="BETA_c_ancALL")
        per_anc = {}
        for s,n,_,_ in ANCS:
            pcol = f"p.value_c_{s}"
            bcol = f"BETA_c_{s}"
            if df_sg is not None and not df_sg.empty and pcol in df_sg.columns:
                pw, _, bw = best_in_window(df_sg, chrom, pos, p_col=pcol, beta_col=bcol)
                per_anc[(n,"p_window")] = pw
                per_anc[(n,"beta_window")] = bw
                per_anc[(n,"rep_class")] = replication_call(pw)
            else:
                per_anc[(n,"p_window")] = None
                per_anc[(n,"beta_window")] = None
                per_anc[(n,"rep_class")] = "no_data"

        rec = dict(phenotype=canon,
                   label=PHENO_LABELS.get(canon, canon),
                   tests_flagging=r.get("tests_flagging",""),
                   Gene=str(r.get("Gene","")).split(";")[0],
                   Chr=chrom, Pos=pos, Allele1=a1, Allele2=a2,
                   ABA_p=aba_p, ABA_beta=aba_beta,
                   SAIGE_cct_p_exact=p_exact,
                   SAIGE_cct_beta_exact=b_exact,
                   SAIGE_cct_p_window=p_win,
                   SAIGE_cct_pos_window=pos_win,
                   SAIGE_cct_beta_window=b_win,
                   SAIGE_cct_concord_exact=concord_label(aba_beta, b_exact),
                   SAIGE_cct_concord_window=concord_label(aba_beta, b_win),
                   SAIGE_cct_replication_class=replication_call(p_win))
        for (n, k), v in per_anc.items():
            rec[f"SAIGE_{n}_{k}"] = v
        rows_a.append(rec)

rep_a = pd.DataFrame(rows_a).sort_values(["phenotype","Chr","Pos"])
out_a = REPLICATION_DIR / "replication_aba_only_in_saige.tsv"
rep_a.to_csv(out_a, sep="\t", index=False)
print(f"[02] wrote {out_a}")

# ─── Shared loci summary (unchanged) ────────────────────────────────────────
print("[02] Building shared-loci replication summary (CCT only)…")
df_cct = all_tabs["CCT"]
shared = df_cct[df_cct["locus_status"]=="shared"].copy()
shared["_sp"] = pd.to_numeric(shared["SAIGE_pvalue"], errors="coerce")
shared["_ap"] = pd.to_numeric(shared["ABA_pvalue"],   errors="coerce")
shared_summary_rows = []
for canon_raw, grp in shared.groupby(shared["phenotype"].astype(str).map(saige_to_canonical)):
    s = grp[grp["tophit_source"]=="SAIGE"]
    a = grp[grp["tophit_source"]=="ABA"]
    n_shared = len(s)
    shared_summary_rows.append(dict(
        phenotype=canon_raw,
        label=PHENO_LABELS.get(canon_raw, canon_raw),
        n_shared=n_shared,
        median_log10p_saige = float(-np.log10(s["_sp"].median())) if not s["_sp"].dropna().empty else None,
        median_log10p_aba   = float(-np.log10(a["_ap"].median())) if not a["_ap"].dropna().empty else None,
    ))
shared_summary = pd.DataFrame(shared_summary_rows)
out_sh = REPLICATION_DIR / "replication_shared_summary.tsv"
shared_summary.to_csv(out_sh, sep="\t", index=False)
print(f"[02] wrote {out_sh}")

# ─── Aggregate JSON ──────────────────────────────────────────────────────────
def class_counts(df, col):
    if df.empty or col not in df.columns: return {}
    return df[col].fillna("no_signal").value_counts().to_dict()

agg = dict(
    n_saige_only_union = len(rep_s),
    n_aba_only_union   = len(rep_a),
    n_shared_cct       = int(shared["tophit_source"].eq("SAIGE").sum()),
    saige_only_in_ABA_META   = class_counts(rep_s, "ABA_META_replication_class"),
    aba_only_in_SAIGE_cct    = class_counts(rep_a, "SAIGE_cct_replication_class"),
    saige_only_concord_window= class_counts(rep_s, "ABA_META_concord_window"),
    aba_only_concord_window  = class_counts(rep_a, "SAIGE_cct_concord_window"),
    window_kb = WINDOW//1000,
)
# per-test-source membership for the union sets
agg["saige_only_per_test_count"] = {
    t: int((rep_s["tests_flagging"].str.contains(t, na=False)).sum())
    for t in TABLE_NAMES
}
agg["aba_only_per_test_count"] = {
    t: int((rep_a["tests_flagging"].str.contains(t, na=False)).sum())
    for t in TABLE_NAMES
}
out_json = REPLICATION_DIR / "replication_summary.json"
with open(out_json, "w") as f:
    json.dump(agg, f, indent=2, default=str)
print(f"[02] wrote {out_json}")
print(f"[02] done.  Summary:\n{json.dumps(agg, indent=2, default=str)}")
