#!/usr/bin/env python3
"""
10_beta_scatterplots.py
───────────────────────
Beta concordance scatterplots between All by All and SAIGE-Tractor at All by
All meta-analysis top hits.

  Fig 12  beta_concordance_overall   — 1 panel:   All by All BETA (meta) vs
                                        SAIGE-Tractor BETA (mega) at ancALL.
  Fig 13  beta_concordance_per_anc   — 5 panels:  per-ancestry BETA concordance
                                        (AFR/EAS/EUR/NatAm/SAS).

Anchor set = every (phenotype × variant) row in the All by All scatter rows
(`tophit_source == "ABA"`) across the CCT-vs-META scatter table — i.e. All by
All meta-analysis top hits. Both BETAs must be present.

Variant-level filter: only variants where **both** BETA estimates are
**strictly positive** are plotted. This removes the sign-driven y = x trend
(which is mostly an artifact of effect-direction concordance at GW-significant
loci) and lets the reader judge effect-magnitude concordance on its own merits.

Manuscript style: very large dots, very large axis/title fonts (target ≥ 22 pt
for axis labels), full method names spelled out, no jargon shortcuts.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as sstats

from common import (
    FIGS_DIR, SCATTER_TABLES, load_scatter, safe_float,
    apply_manuscript_style,
    NAME_SAIGE, NAME_ABA, NAME_ABA_META, NAME_SAIGE_MEGA,
    COLOR_SAIGE, COLOR_ABA, ANC_COLORS,
)
apply_manuscript_style()
# Override defaults further — Figs 12/13 are headline concordance figures and
# need to be readable at half-column-width on a printed page.
plt.rcParams.update({
    "font.size":         18,
    "axes.titlesize":    24,
    "axes.labelsize":    22,
    "xtick.labelsize":   18,
    "ytick.labelsize":   18,
    "legend.fontsize":   18,
    "figure.titlesize":  26,
})

def savefig(fig, name):
    fig.savefig(FIGS_DIR / f"{name}.pdf")
    fig.savefig(FIGS_DIR / f"{name}.png")
    plt.close(fig)
    print(f"[10] wrote {name}.pdf + .png")

# Load the CCT scatter table (which carries all the ancestry-resolved columns)
df = load_scatter(SCATTER_TABLES["cct"])
# Anchor on All by All top hits (shared + ABA_only).
aba_rows = df[df["tophit_source"] == "ABA"].copy()

def _scatter_with_diag(ax, x, y, color, alpha=0.75, point_size=120):
    """Render a beta-vs-beta scatter with diagonal y=x and Pearson r in a box.

    Variants are **pre-filtered to both x > 0 and y > 0** so the plot shows
    effect-size magnitude concordance rather than the (almost trivially
    perfect) sign-driven y = x trend.  Plot extent is the positive quadrant
    only.
    """
    finite = (~np.isnan(x)) & (~np.isnan(y))
    pos    = (x > 0) & (y > 0)
    keep   = finite & pos
    n_all  = int(finite.sum())
    x = x[keep]; y = y[keep]
    if len(x) == 0:
        ax.text(0.5, 0.5, "no positive-β\nvariants",
                transform=ax.transAxes, ha="center", va="center", fontsize=18)
        return None
    ax.scatter(x, y, s=point_size, color=color, alpha=alpha,
               edgecolor="white", linewidth=1.0)
    lim = max(np.concatenate([x, y])) * 1.10
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.plot([0, lim], [0, lim], ls="--", lw=2.0, color="#333", label="y = x")
    r, _ = sstats.pearsonr(x, y)
    # within ±20% magnitude band
    ratio = y / x
    within20 = int(((ratio >= 0.80) & (ratio <= 1.20)).sum())
    pct20 = 100 * within20 / len(x)
    ax.text(0.05, 0.96,
            f"N (positive-β) = {len(x):,}\n"
            f"N (all variants) = {n_all:,}\n"
            f"Pearson r = {r:.3f}\n"
            f"within ±20%: {pct20:.0f}%",
            transform=ax.transAxes, va="top", ha="left",
            fontsize=17, family="monospace",
            bbox=dict(facecolor="white", edgecolor="#cfd8dc", alpha=0.96, pad=8))
    return r

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 12 — Overall meta vs mega concordance
# ═══════════════════════════════════════════════════════════════════════════
def fig12_overall():
    x = pd.to_numeric(aba_rows["ABA_META_BETA"],      errors="coerce").to_numpy()
    y = pd.to_numeric(aba_rows["SAIGE_BETA_c_ancALL"], errors="coerce").to_numpy()
    fig, ax = plt.subplots(figsize=(11, 11))
    r = _scatter_with_diag(ax, x, y, color=COLOR_SAIGE, point_size=160)
    ax.set_xlabel(f"{NAME_ABA} β (meta-analysis)")
    ax.set_ylabel(f"{NAME_SAIGE} β (mega-analysis)")
    ax.set_title(f"Effect-size magnitude concordance\n"
                 f"{NAME_ABA} top hits — positive-β variants only",
                 loc="left", weight="bold", pad=18, fontsize=24)
    ax.legend(loc="lower right", framealpha=0.95, fontsize=18)
    fig.tight_layout()
    savefig(fig, "fig12_beta_concordance_overall")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 13 — Per-ancestry concordance (5 panels)
# ═══════════════════════════════════════════════════════════════════════════
def fig13_per_ancestry():
    names     = ["AFR","EAS","EUR","NatAm","SAS"]
    aba_names = ["AFR","EAS","EUR","AMR","SAS"]
    sufs      = ["anc1","anc2","anc3","anc4","anc5"]
    fig, axes = plt.subplots(2, 3, figsize=(24, 16))
    axes_flat = axes.flatten()
    for ax, name, aba_name, suf in zip(axes_flat[:5], names, aba_names, sufs):
        x = pd.to_numeric(aba_rows[f"ABA_{aba_name}_BETA"], errors="coerce").to_numpy()
        y = pd.to_numeric(aba_rows[f"SAIGE_BETA_c_{suf}"],  errors="coerce").to_numpy()
        _scatter_with_diag(ax, x, y, color=ANC_COLORS.get(name, "#444"),
                           point_size=180)
        ax.set_xlabel(f"{NAME_ABA} β")
        ax.set_ylabel(f"{NAME_SAIGE} β")
        ax.set_title(name, loc="center", weight="bold", fontsize=28)
    # Use the 6th cell as a key/legend caption
    ax6 = axes_flat[5]
    ax6.axis("off")
    ax6.text(0.02, 0.92,
             "Per-ancestry effect-size\nmagnitude concordance",
             fontsize=24, weight="bold", color="#1a237e",
             transform=ax6.transAxes, va="top")
    ax6.text(0.02, 0.65,
             f"x-axis: {NAME_ABA} β at {NAME_ABA} top hits\n"
             f"        (per global-ancestry meta-analysis\n"
             f"         estimate at the variant)\n\n"
             f"y-axis: {NAME_SAIGE} β at the same variant\n"
             f"        (per-local-ancestry strata,\n"
             f"         conditioned on local-ancestry dosage)\n\n"
             f"Filter: only variants where BOTH β > 0\n"
             f"        (removes the sign-driven y = x trend)\n\n"
             f"Dashed line = y = x (perfect concordance)\n"
             f"r = Pearson correlation on positive subset\n"
             f"within ±20% = fraction of variants with\n"
             f"  SAIGE β / ABA β between 0.8 and 1.2",
             fontsize=16, color="#333",
             transform=ax6.transAxes, va="top", family="monospace")
    fig.suptitle(f"Per-ancestry effect-size magnitude concordance — "
                 f"{NAME_SAIGE} vs {NAME_ABA}",
                 x=0.04, ha="left", y=0.995, weight="bold", fontsize=26)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    savefig(fig, "fig13_beta_concordance_per_ancestry")

if __name__ == "__main__":
    fig12_overall()
    fig13_per_ancestry()
    print(f"[10] done. Outputs in {FIGS_DIR}")
