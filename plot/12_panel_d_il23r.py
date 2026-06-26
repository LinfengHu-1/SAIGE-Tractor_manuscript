#!/usr/bin/env python3
"""
12_panel_d_il23r.py
───────────────────
Standalone clean bar chart for IL23R / Crohn's disease.

Per-ancestry sample size analysed at chr1:67,240,275 — global-ancestry
meta-analysis (single EUR stratum only) vs SAIGE-Tractor (every ancestry,
EUR expanded by recovering haplotypes from admixed individuals).

No surrounding card / narrative / stats strip — just the chart itself, so
it can be dropped directly into a manually-composed Figure 1.

Output:  manuscript_figures/panel_d_il23r.{pdf,png}
"""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from common import FIGS_DIR, REPLICATION_DIR, safe_float, apply_manuscript_style

apply_manuscript_style()
plt.rcParams.update({"font.family": "DejaVu Sans"})

# Colour palette — Okabe-Ito
CB_AFR    = "#D55E00"
CB_EAS    = "#56B4E9"
CB_EUR    = "#0072B2"
CB_NATAM  = "#E69F00"
CB_SAS    = "#CC79A7"

# Load sample-size table
ss = pd.read_csv(REPLICATION_DIR / "sample_size_table.tsv", sep="\t")
row = ss[ss["phenotype"].astype(str) == "GI_522.11"].iloc[0]

ancs = [("AFR",   "AFR_ABA",   "AFR_TRACTOR",   CB_AFR),
        ("EAS",   "EAS_ABA",   "EAS_TRACTOR",   CB_EAS),
        ("EUR",   "EUR_ABA",   "EUR_TRACTOR",   CB_EUR),
        ("NatAm", "NatAm_ABA", "NatAm_TRACTOR", CB_NATAM),
        ("SAS",   "SAS_ABA",   "SAS_TRACTOR",   CB_SAS)]

aba_vals   = []
saige_vals = []
cols       = []
labels     = []
for name, ab_col, sg_col, col in ancs:
    ab = safe_float(row[ab_col])
    sg = safe_float(row[sg_col])
    ab = 0 if (ab is None or math.isnan(ab)) else ab
    sg = 0 if (sg is None or math.isnan(sg)) else sg
    aba_vals.append(ab); saige_vals.append(sg)
    cols.append(col); labels.append(name)

# ───────────────────────────────────────────────────────────────────────────
# Figure — horizontal bars, two per ancestry (faint = global, solid = SAIGE)
# ───────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5.5))

ypos = np.arange(len(ancs))[::-1]
bh   = 0.38
gap  = 0.04

ax.barh(ypos + bh/2 + gap, aba_vals,    height=bh,
         color=cols, alpha=0.28, edgecolor="white", linewidth=0.8,
         label="global-ancestry meta-analysis")
ax.barh(ypos - bh/2 - gap, saige_vals,  height=bh,
         color=cols, alpha=1.00, edgecolor="white", linewidth=0.8,
         label="SAIGE-Tractor")

max_val = max(saige_vals + aba_vals)
pad     = max_val * 0.015
for yy, ab, sg in zip(ypos, aba_vals, saige_vals):
    ax.text(ab + pad, yy + bh/2 + gap,
             f"{int(ab):,}" if ab > 0 else "not tested",
             ha="left", va="center", fontsize=11,
             color="#555")
    ax.text(sg + pad, yy - bh/2 - gap, f"{int(sg):,}",
             ha="left", va="center", fontsize=12, color="#111",
             weight="bold")

ax.set_yticks(ypos)
ax.set_yticklabels(labels, fontsize=14, weight="bold")
for tl, c in zip(ax.get_yticklabels(), cols):
    tl.set_color(c)

ax.set_xlim(0, max_val * 1.30)
ax.set_xlabel("Individuals analysed at this variant", fontsize=13)

# Custom legend with neutral grey swatches showing the two alpha levels
legend_handles = [
    Rectangle((0, 0), 1, 1, fc="#777", alpha=0.28,
              label="global-ancestry meta-analysis"),
    Rectangle((0, 0), 1, 1, fc="#777", alpha=1.00,
              label="SAIGE-Tractor"),
]
ax.legend(handles=legend_handles, loc="lower right",
           fontsize=11, frameon=False)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.grid(axis="x", color="#eee", linewidth=0.8)
ax.set_axisbelow(True)
ax.tick_params(axis="y", which="both", left=False)

out_pdf = FIGS_DIR / "panel_d_il23r.pdf"
out_png = FIGS_DIR / "panel_d_il23r.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=200)
plt.close(fig)
print(f"[panel D · IL23R] wrote {out_pdf}")
print(f"[panel D · IL23R] wrote {out_png}")
