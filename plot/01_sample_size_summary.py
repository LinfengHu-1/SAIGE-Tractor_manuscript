#!/usr/bin/env python3
"""
01_sample_size_summary.py
─────────────────────────
Per-phenotype × ancestry sample-size table directly from the summary statistics.

Sources
─────────
- SAIGE-Tractor: read N_haplo_anc{1..5} from every SAIGE-Tractor sumstat in
  saige_sumstat/, take the MAX across all variants, divide by 2 and round to
  the nearest integer (haplotypes → individual-equivalent count).
- All by All: per-ancestry sample sizes are NOT present in the All by All
  per-ancestry sumstats (which only carry CHR/POS/Allele/BETA/SE/Pvalue/AF),
  so they are read from the All by All phenotype CSV
  (All_by_All_phenotypes_v7_Analyzed_phenotypes.csv) where ancestries reported
  there were the ones the All by All pipeline ran for each phenotype.

Output (replication_output/sample_size_table.tsv) — 15 columns:
    phenotype  label  AFR_ABA  AFR_TRACTOR  EAS_ABA  EAS_TRACTOR
    EUR_ABA  EUR_TRACTOR  NatAm_ABA  NatAm_TRACTOR
    SAS_ABA  SAS_TRACTOR  MID_ABA  TOTAL_ABA  TOTAL_TRACTOR

A long-form table (one row per phenotype × ancestry) and a JSON snapshot are
also written for downstream use.
"""
from __future__ import annotations
import json
import math
import numpy as np
import pandas as pd
from common import (
    SAIGE_SS, ABA_SS, REPLICATION_DIR, PHENO_LABELS,
    list_saige_phenos, saige_to_canonical, load_aba_pheno_csv,
    aba_sample_size,
)

# ─── Per-ancestry mapping ────────────────────────────────────────────────────
# (display_name, saige_n_haplo_column,    aba_csv_ancestry_code)
ANCS = [
    ("AFR",   "N_haplo_anc1", "afr"),
    ("EAS",   "N_haplo_anc2", "eas"),
    ("EUR",   "N_haplo_anc3", "eur"),
    ("NatAm", "N_haplo_anc4", "amr"),
    ("SAS",   "N_haplo_anc5", "sas"),
]

def _div2(x):
    if x is None or pd.isna(x): return None
    return int(round(float(x) / 2))

# ─── SAIGE-Tractor: per-ancestry max(N_haplo) per phenotype ─────────────────
def saige_max_n(saige_pheno: str) -> dict[str, int | None]:
    """Stream the SAIGE-Tractor sumstat in chunks; return per-ancestry
    individual-equivalent N. Auto-detects column scheme:
      - Continuous traits:  use max(N_haplo_anc{i}) / 2.
      - Binary traits:      use (max(N_case_anc{i}) + max(N_ctrl_anc{i})) / 2.
    """
    path = SAIGE_SS / f"merged_saigetractor_{saige_pheno}.txt.gz"
    if not path.exists():
        return {name: None for name, _, _ in ANCS}
    # Peek the header to determine which N columns exist.
    header = pd.read_csv(path, sep="\t", compression="gzip", nrows=0).columns.tolist()
    has_nhap = any(c.startswith("N_haplo_anc") for c in header)
    has_case = any(c.startswith("N_case_anc")  for c in header)
    if has_nhap:
        cols = [f"N_haplo_anc{i}" for i in range(1, 6)]
    elif has_case:
        cols = [f"N_case_anc{i}" for i in range(1, 6)] + \
               [f"N_ctrl_anc{i}" for i in range(1, 6)]
    else:
        return {name: None for name, _, _ in ANCS}
    cols = [c for c in cols if c in header]
    maxes = {c: 0.0 for c in cols}
    rdr = pd.read_csv(path, sep="\t", compression="gzip",
                      usecols=cols, chunksize=200_000, low_memory=False)
    for ch in rdr:
        for c in cols:
            v = pd.to_numeric(ch[c], errors="coerce")
            if not v.empty:
                m = v.max()
                if pd.notna(m) and m > maxes[c]:
                    maxes[c] = float(m)
    out = {}
    for i, (name, _, _) in enumerate(ANCS, start=1):
        if has_nhap:
            tot_haps = maxes.get(f"N_haplo_anc{i}", 0)
        else:
            tot_haps = maxes.get(f"N_case_anc{i}", 0) + maxes.get(f"N_ctrl_anc{i}", 0)
        out[name] = _div2(tot_haps) if tot_haps > 0 else None
    return out

# ─── All by All sample sizes from the CSV ───────────────────────────────────
df_csv = load_aba_pheno_csv()

def aba_n_per_anc(canon_pheno: str) -> dict[str, int | None]:
    out = {}
    for name, _, anc_csv in ANCS:
        _, _, n_total = aba_sample_size(df_csv, canon_pheno, anc_csv)
        out[name] = None if n_total is None else int(n_total)
    _, _, n_mid = aba_sample_size(df_csv, canon_pheno, "mid")
    out["MID"] = None if n_mid is None else int(n_mid)
    return out

# ─── Driver ─────────────────────────────────────────────────────────────────
saige_phenos = list_saige_phenos()
print(f"[01] found {len(saige_phenos)} SAIGE-Tractor phenotypes in {SAIGE_SS}")

rows = []
long_rows = []
json_payload = {}

for saige_p in sorted(saige_phenos):
    canon = saige_to_canonical(saige_p)
    label = PHENO_LABELS.get(canon, PHENO_LABELS.get(saige_p, canon))
    print(f"[01]   reading {saige_p} …")
    s_n = saige_max_n(saige_p)
    a_n = aba_n_per_anc(canon)
    # Totals (sum of available per-ancestry, skipping None)
    total_aba     = sum(v for k,v in a_n.items() if v is not None)
    total_tractor = sum(v for k,v in s_n.items() if v is not None)
    row = {
        "phenotype":      canon,
        "label":          label,
        "AFR_ABA":        a_n["AFR"],
        "AFR_TRACTOR":    s_n["AFR"],
        "EAS_ABA":        a_n["EAS"],
        "EAS_TRACTOR":    s_n["EAS"],
        "EUR_ABA":        a_n["EUR"],
        "EUR_TRACTOR":    s_n["EUR"],
        "NatAm_ABA":      a_n["NatAm"],
        "NatAm_TRACTOR":  s_n["NatAm"],
        "SAS_ABA":        a_n["SAS"],
        "SAS_TRACTOR":    s_n["SAS"],
        "MID_ABA":        a_n["MID"],
        "TOTAL_ABA":      total_aba if total_aba > 0 else None,
        "TOTAL_TRACTOR":  total_tractor if total_tractor > 0 else None,
    }
    rows.append(row)
    json_payload[canon] = {"label": label, **{k: v for k,v in row.items() if k not in ("phenotype","label")}}
    for name, _, _ in ANCS:
        long_rows.append(dict(phenotype=canon, label=label, ancestry=name,
                              ABA_N=a_n[name], TRACTOR_N=s_n[name]))
    long_rows.append(dict(phenotype=canon, label=label, ancestry="MID",
                          ABA_N=a_n["MID"], TRACTOR_N=None))

wide_df = pd.DataFrame(rows)
long_df = pd.DataFrame(long_rows)

wide_path = REPLICATION_DIR / "sample_size_table.tsv"
long_path = REPLICATION_DIR / "sample_size_table_long.tsv"
json_path = REPLICATION_DIR / "sample_size_table.json"

wide_df.to_csv(wide_path, sep="\t", index=False)
long_df.to_csv(long_path, sep="\t", index=False)
with open(json_path, "w") as f:
    json.dump(json_payload, f, indent=2, default=str)

print(f"\n[01] wrote {wide_path}")
print(f"[01] wrote {long_path}")
print(f"[01] wrote {json_path}")
print(f"[01] {len(rows)} phenotypes processed.")
print(f"[01] NOTE: TRACTOR columns are max(N_haplo_anc{{i}}) ÷ 2 (rounded), giving")
print( "          individual-equivalent counts comparable to ABA sample sizes.")
print(f"          ABA per-ancestry counts come from the All by All phenotype CSV")
print( "          because the All by All per-ancestry sumstats don't carry N columns.")
