#!/usr/bin/env python3
"""
12_panel_d_hpr.py
─────────────────
Standalone clean see-saw / forest-bar chart for HPR / TXNL4B and LDL
cholesterol at chr16:72,080,103 (rs217181).

Five per-ancestry β estimates with 95% CI whiskers, oriented so the
opposing-direction AFR-positive / EUR-negative pattern is immediately
obvious.  No surrounding card / narrative / stats strip — just the chart
itself, so it can be dropped directly into a manually-composed Figure 1.

Output:  manuscript_figures/panel_d_hpr.{pdf,png}
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import (FIGS_DIR, SCATTER_TABLES, load_scatter,
                     safe_float, apply_manuscript_style)

apply_manuscript_style()
plt.rcParams.update({"font.family": "DejaVu Sans"})

CB_AFR    = "#D55E00"
CB_EAS    = "#56B4E9"
CB_EUR    = "#0072B2"
CB_NATAM  = "#E69F00"
CB_SAS    = "#CC79A7"

cct = load_scatter(SCATTER_TABLES["cct"])
sub = cct[(cct["phenotype"].astype(str) == "3028288") &
          (cct["Gene"].astype(str).str.contains("HPR", na=False))]
row = sub.iloc[0]

ancs = [("AFR",   "anc1", CB_AFR),
        ("EAS",   "anc2", CB_EAS),
        ("EUR",   "anc3", CB_EUR),
        ("NatAm", "anc4", CB_NATAM),
        ("SAS",   "anc5", CB_SAS)]

betas, ses, cols, labels = [], [], [], []
for name, suf, col in ancs:
    b  = safe_float(row[f"SAIGE_BETA_c_{suf}"])
    se = safe_float(row[f"SAIGE_SE_c_{suf}"])
    betas.append(b if b is not None else 0.0)
    ses.append(se if se is not None else 0.0)
    cols.append(col); labels.append(name)

# ───────────────────────────────────────────────────────────────────────────
# Figure — vertical bars sticking up / down from the zero line, with 95% CI
# ───────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5.5))

xs = np.arange(len(ancs))
ax.axhline(0, color="#37474f", lw=1.5, zorder=4)
ax.bar(xs, betas, width=0.62, color=cols, edgecolor="white", linewidth=1.0,
        zorder=5)

# 95% CI whiskers
for x, b, se in zip(xs, betas, ses):
    lo = b - 1.96 * se; hi = b + 1.96 * se
    ax.plot([x, x], [lo, hi], color="#37474f", lw=1.6, zorder=6)
    ax.plot([x - 0.12, x + 0.12], [lo, lo], color="#37474f", lw=1.6,
             zorder=6)
    ax.plot([x - 0.12, x + 0.12], [hi, hi], color="#37474f", lw=1.6,
             zorder=6)

# Numeric β labels, placed outside the CI so they never collide with the
# whiskers
for x, b, se in zip(xs, betas, ses):
    lo = b - 1.96 * se; hi = b + 1.96 * se
    if b >= 0:
        ax.text(x, hi + 0.012, f"{b:+.3f}",
                 ha="center", va="bottom", fontsize=12,
                 color="#111", weight="bold")
    else:
        ax.text(x, lo - 0.012, f"{b:+.3f}",
                 ha="center", va="top", fontsize=12,
                 color="#111", weight="bold")

ax.set_xticks(xs)
ax.set_xticklabels(labels, fontsize=14, weight="bold")
for tl, c in zip(ax.get_xticklabels(), cols):
    tl.set_color(c)

# y-limits — symmetric, generous so the +0.108 label fits cleanly
top    = max([b + 1.96 * se for b, se in zip(betas, ses)])
bottom = min([b - 1.96 * se for b, se in zip(betas, ses)])
pad    = max(abs(top), abs(bottom)) * 0.30
ax.set_ylim(bottom - pad, top + pad)
ax.set_ylabel("β  per local-ancestry stratum  (95 % CI)", fontsize=13)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.grid(axis="y", color="#eee", linewidth=0.8)
ax.set_axisbelow(True)
ax.tick_params(axis="x", which="both", bottom=False)

out_pdf = FIGS_DIR / "panel_d_hpr.pdf"
out_png = FIGS_DIR / "panel_d_hpr.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=200)
plt.close(fig)
print(f"[panel D · HPR] wrote {out_pdf}")
print(f"[panel D · HPR] wrote {out_png}")
