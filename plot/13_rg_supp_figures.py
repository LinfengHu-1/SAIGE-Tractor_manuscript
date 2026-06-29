#!/usr/bin/env python3
"""
13_rg_supp_figures.py
─────────────────────
Build Supplementary Figures 12–14 for the rg / h² (step 3) section from Wei's
source data in `Manuscript/Wei_supplementary_tables.xlsx`. Wei provided the raw
numbers but not the rendered panels; this script produces production-ready,
high-resolution, colourblind-safe versions.

  SFig 12  — Precision, bias, and calibration of the cross-ancestry genetic-
             correlation estimator as a function of sample size.
             (sheet: SFig12_calibration_summary)
  SFig 13  — Computational cost: matrix-free STadmix vs explicit-GRM radmix,
             runtime and peak memory vs N.
             (sheet: SFig13_cost)
  SFig 14  — Genotype centering determines the estimand across five admixture
             cohorts, under centered vs uncentered simulating truth.
             (sheet: SFig14_centering_raw)

Style: Okabe-Ito colourblind-safe palette, large fonts (apply_manuscript_style),
vector-editable PDF + 400-dpi PNG.

Output: manuscript_figures/SFig12_*, SFig13_*, SFig14_* (.pdf + .png)
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from common import ROOT, FIGS_DIR, apply_manuscript_style
apply_manuscript_style()

XLSX = ROOT / "Manuscript" / "Wei_supplementary_tables.xlsx"

# ── Okabe-Ito colourblind-safe palette ───────────────────────────────────────
CB_BLUE   = "#0072B2"   # deep blue   — centered (STadmix, our default)
CB_ORANGE = "#E69F00"   # orange      — no-center (uncentered)
CB_GREEN  = "#009E73"   # bluish green — radmix (MLE)
CB_SKY    = "#56B4E9"   # sky blue
CB_PURPLE = "#CC79A7"   # reddish purple
CB_VERM   = "#D55E00"   # vermilion
GREY      = "#666666"

# true-rg colour ramp (for SFig 12)
RG_COLORS = {0.0: CB_SKY, 0.6: CB_BLUE, 1.0: CB_PURPLE}

PNG_DPI = 400


def _read(sheet: str) -> pd.DataFrame:
    """Read a sheet; row 0 is a title, row 1 is the header."""
    return pd.read_excel(XLSX, sheet_name=sheet, skiprows=1, engine="openpyxl")


def savefig(fig, name):
    fig.savefig(FIGS_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGS_DIR / f"{name}.png", bbox_inches="tight", dpi=PNG_DPI)
    plt.close(fig)
    print(f"[13] wrote {name}.pdf + .png")


# ═══════════════════════════════════════════════════════════════════════════
# SFig 12 — precision, bias, calibration vs sample size
# ═══════════════════════════════════════════════════════════════════════════
def sfig12():
    df = _read("SFig12_calibration_summary")
    df = df.dropna(subset=["N"]).copy()
    df["N"] = df["N"].astype(int)
    df["true_rg"] = df["true_rg"].astype(float)
    Ns = sorted(df["N"].unique())
    rgs = sorted(df["true_rg"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(21, 6.6))

    # ── Panel A: precision (empirical SD) vs N, log-log, with 1/sqrt(N) guide ─
    axA = axes[0]
    for rg in rgs:
        sub = df[df["true_rg"] == rg].sort_values("N")
        axA.plot(sub["N"], sub["empirical_SD"], "-o",
                 color=RG_COLORS.get(rg, GREY), lw=2.4, ms=10,
                 label=f"true $r_g$ = {rg:g}")
    # 1/sqrt(N) reference anchored to the largest-SD series' first point
    anchor = df[df["true_rg"] == rgs[0]].sort_values("N").iloc[0]
    refN = np.array(Ns, dtype=float)
    ref = anchor["empirical_SD"] * np.sqrt(anchor["N"] / refN)
    axA.plot(refN, ref, ls="--", lw=2.0, color=GREY, label=r"$\propto 1/\sqrt{N}$")
    axA.set_xscale("log"); axA.set_yscale("log")
    axA.set_xlabel("Sample size $N$")
    axA.set_ylabel("Empirical SD of $\\hat{r}_g$")
    axA.set_title("a  Precision improves with sample size",
                  loc="left", weight="bold")
    axA.set_xticks(Ns); axA.set_xticklabels([f"{n//1000}k" for n in Ns])
    axA.legend(frameon=False, loc="upper right")
    axA.grid(True, which="both", color="#eee", lw=0.6); axA.set_axisbelow(True)

    # ── Panel B: bias vs N, with 2.5–97.5% replicate spread whiskers ─────────
    axB = axes[1]
    x_idx = {n: i for i, n in enumerate(Ns)}
    width = 0.22
    for j, rg in enumerate(rgs):
        sub = df[df["true_rg"] == rg].sort_values("N")
        xs = np.array([x_idx[n] for n in sub["N"]]) + (j - (len(rgs)-1)/2) * width
        lo = (sub["bias"] - sub["bias_lo"]).clip(lower=0).values
        hi = (sub["bias_hi"] - sub["bias"]).clip(lower=0).values
        axB.errorbar(xs, sub["bias"], yerr=[lo, hi], fmt="o",
                     color=RG_COLORS.get(rg, GREY), ms=10, lw=2.0,
                     capsize=5, elinewidth=2.0, label=f"true $r_g$ = {rg:g}")
    axB.axhline(0, color="#222", lw=1.4, ls="-")
    axB.set_xticks(range(len(Ns))); axB.set_xticklabels([f"{n//1000}k" for n in Ns])
    axB.set_xlabel("Sample size $N$")
    axB.set_ylabel("Bias  $\\hat{r}_g - r_g$")
    axB.set_title("b  Approximately unbiased\n(points = mean; whiskers = 2.5–97.5% of replicates)",
                  loc="left", weight="bold", fontsize=15)
    axB.legend(frameon=False, loc="upper right", ncol=1)
    axB.grid(True, axis="y", color="#eee", lw=0.6); axB.set_axisbelow(True)

    # ── Panel C: calibration — reported SE vs empirical SD, identity line ────
    axC = axes[2]
    mx = max(df["empirical_SD"].max(), df["mean_reported_SE"].max()) * 1.08
    axC.plot([0, mx], [0, mx], ls="--", lw=2.0, color=GREY, label="identity")
    for rg in rgs:
        sub = df[df["true_rg"] == rg]
        sizes = 60 + (np.log10(sub["N"]) - np.log10(min(Ns))) * 240
        axC.scatter(sub["empirical_SD"], sub["mean_reported_SE"],
                    s=sizes, color=RG_COLORS.get(rg, GREY), alpha=0.85,
                    edgecolor="white", linewidth=1.2, label=f"true $r_g$ = {rg:g}")
    axC.set_xlim(0, mx); axC.set_ylim(0, mx); axC.set_aspect("equal")
    axC.set_xlabel("Empirical SD across replicates")
    axC.set_ylabel("Mean reported jackknife SE")
    cov_lo, cov_hi = df["coverage95"].min(), df["coverage95"].max()
    axC.set_title(f"c  Well-calibrated SE\n(95% CI coverage {cov_lo:.2f}–{cov_hi:.2f}; "
                  "marker size $\\propto N$)",
                  loc="left", weight="bold", fontsize=15)
    axC.legend(frameon=False, loc="lower right")
    axC.grid(True, color="#eee", lw=0.6); axC.set_axisbelow(True)

    fig.suptitle("Supplementary Figure 12 — Precision, bias, and calibration of the "
                 "cross-ancestry genetic-correlation estimator vs sample size",
                 x=0.012, ha="left", y=1.02, weight="bold", fontsize=18)
    fig.tight_layout()
    savefig(fig, "SFig12_rg_calibration_vs_N")


# ═══════════════════════════════════════════════════════════════════════════
# SFig 13 — computational cost: STadmix (matrix-free) vs radmix (explicit GRM)
# ═══════════════════════════════════════════════════════════════════════════
def sfig13():
    df = _read("SFig13_cost")
    df = df.dropna(subset=["N"]).copy()
    df["N"] = df["N"].astype(int)
    tools = {
        "STadmix":               dict(c=CB_BLUE,   m="o", lbl="STadmix (matrix-free RHE)"),
        "radmix (explicit GRM)": dict(c=CB_ORANGE, m="s", lbl="radmix (explicit GRM)"),
    }
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.4))

    # ── Panel A: wall-clock runtime ──
    axA = axes[0]
    for tool, st in tools.items():
        sub = df[df["tool"] == tool].sort_values("N")
        axA.plot(sub["N"], sub["wall_s"], "-"+st["m"], color=st["c"],
                 lw=2.6, ms=11, label=st["lbl"])
    axA.set_xscale("log")
    axA.set_xlabel("Sample size $N$")
    axA.set_ylabel("Wall-clock time per chromosome (s)")
    axA.set_title("a  Runtime", loc="left", weight="bold")
    Ns = sorted(df["N"].unique())
    axA.set_xticks(Ns); axA.set_xticklabels([f"{n//1000}k" for n in Ns])
    axA.legend(frameon=False, loc="upper left")
    axA.grid(True, which="both", color="#eee", lw=0.6); axA.set_axisbelow(True)

    # ── Panel B: peak memory, with quadratic extrapolation for radmix ──
    axB = axes[1]
    for tool, st in tools.items():
        sub = df[df["tool"] == tool].sort_values("N")
        axB.plot(sub["N"], sub["peak_GB"], "-"+st["m"], color=st["c"],
                 lw=2.6, ms=11, label=st["lbl"])
    # Extrapolate radmix explicit-GRM memory (∝ N²) beyond the measured range
    rad = df[df["tool"] == "radmix (explicit GRM)"].sort_values("N")
    coef = np.polyfit(rad["N"], rad["peak_GB"], 2)          # a N² + b N + c
    xmax_meas = rad["N"].max()
    xext = np.array([xmax_meas, 50000, 100000, 245000], dtype=float)
    yext = np.polyval(coef, xext)
    axB.plot(xext, yext, ls=":", lw=2.4, color=CB_ORANGE,
             label="radmix, quadratic extrapolation ($\\propto N^2$)")
    axB.scatter([245000], [np.polyval(coef, 245000)], marker="*", s=420,
                color=CB_VERM, edgecolor="white", linewidth=1.4, zorder=12)
    axB.annotate(f"All of Us scale\n($N\\approx$245k): ~{np.polyval(coef,245000):.0f} GB",
                 xy=(245000, np.polyval(coef, 245000)),
                 xytext=(0.42, 0.78), textcoords="axes fraction",
                 fontsize=13, color=CB_VERM, weight="bold",
                 ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", color=CB_VERM, lw=1.8))
    # typical compute-node RAM ceilings
    for ceil, lab in [(128, "128 GB node")]:
        axB.axhline(ceil, ls="--", lw=1.5, color=GREY)
        axB.text(1050, ceil*1.04, lab, fontsize=12, color=GREY, va="bottom")
    axB.set_xscale("log"); axB.set_yscale("log")
    axB.set_xlabel("Sample size $N$")
    axB.set_ylabel("Peak memory (GB)")
    axB.set_title("b  Peak memory", loc="left", weight="bold")
    xt = Ns + [50000, 245000]
    axB.set_xticks(xt); axB.set_xticklabels([f"{n//1000}k" for n in xt], rotation=0)
    axB.legend(frameon=False, loc="lower right", fontsize=13)
    axB.grid(True, which="both", color="#eee", lw=0.6); axB.set_axisbelow(True)

    fig.suptitle("Supplementary Figure 13 — Computational cost: matrix-free STadmix "
                 "vs explicit-GRM radmix",
                 x=0.012, ha="left", y=1.02, weight="bold", fontsize=18)
    fig.tight_layout()
    savefig(fig, "SFig13_rg_compute_cost")


# ═══════════════════════════════════════════════════════════════════════════
# SFig 14 — centering determines the estimand across admixture cohorts
# ═══════════════════════════════════════════════════════════════════════════
COHORT_LABEL = {
    "eur50": "Balanced 2-way\n(50% / 50%)",
    "eur80": "Unbalanced 2-way\n(80% / 20%)",
    "bal3":  "Balanced 3-way\n(1/3 each)",
    "imb3":  "Unbalanced 3-way\n(0.6 / 0.2 / 0.2)",
    "imb3b": "Unbalanced 3-way\n(0.6 / 0.3 / 0.1)",
}
COHORT_ORDER = ["eur50", "eur80", "bal3", "imb3", "imb3b"]
BASIS_LABEL = {"cen": "Simulated under centered truth",
               "unc": "Simulated under uncentered truth"}
BASIS_ORDER = ["cen", "unc"]
STYLE_CENTERED = "centered (RHE in STadmix)"
STYLE_NOCENTER = "no-center (RHE in STadmix)"


def sfig14():
    df = _read("SFig14_centering_raw")
    df = df.dropna(subset=["cohort"]).copy()
    df["true_rg"] = df["true_rg"].astype(float)
    df["rg"] = pd.to_numeric(df["rg"], errors="coerce")
    df = df.dropna(subset=["rg"])
    trs = sorted(df["true_rg"].unique())

    nrows, ncols = len(COHORT_ORDER), len(BASIS_ORDER)
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 22), squeeze=False)

    # Evenly-spaced categorical x positions for the true_rg levels so the two
    # grouped boxes never collide (true_rg = 0.9 and 1.0 are only 0.1 apart on
    # the natural scale).
    xcat = np.arange(len(trs), dtype=float)
    box_w = 0.34
    off = 0.20

    for i, cohort in enumerate(COHORT_ORDER):
        for j, basis in enumerate(BASIS_ORDER):
            ax = axes[i][j]
            sub = df[(df["cohort"] == cohort) & (df["basis"] == basis)]
            # identity reference: truth value at each categorical position
            ax.plot(xcat, trs, ls="--", lw=1.8, color=GREY, zorder=1)
            ax.scatter(xcat, trs, s=28, color=GREY, zorder=2)
            for style, scol, soff in [(STYLE_CENTERED, CB_BLUE, -off),
                                      (STYLE_NOCENTER, CB_ORANGE, +off)]:
                data = [sub[(sub["style"] == style) & (sub["true_rg"] == t)]["rg"].values
                        for t in trs]
                positions = xcat + soff
                bp = ax.boxplot(data, positions=positions, widths=box_w,
                                patch_artist=True, manage_ticks=False,
                                showfliers=False, zorder=3)
                for box in bp["boxes"]:
                    box.set(facecolor=scol, alpha=0.65, edgecolor=scol, linewidth=1.4)
                for med in bp["medians"]:
                    med.set(color="black", linewidth=2.0)
                for w in bp["whiskers"] + bp["caps"]:
                    w.set(color=scol, linewidth=1.6)
            ax.set_xlim(-0.6, len(trs) - 0.4)
            ax.set_ylim(-0.45, 1.5)
            ax.set_xticks(xcat); ax.set_xticklabels([f"{t:g}" for t in trs])
            ax.grid(True, axis="y", color="#eee", lw=0.6); ax.set_axisbelow(True)
            if i == 0:
                ax.set_title(BASIS_LABEL[basis], weight="bold", fontsize=16, pad=10)
            if i == nrows - 1:
                ax.set_xlabel("True $r_g$")
            if j == 0:
                ax.set_ylabel(f"{COHORT_LABEL[cohort]}\n\nestimated $\\hat{{r}}_g$",
                              fontsize=13)

    legend_handles = [
        Patch(facecolor=CB_BLUE, alpha=0.65, edgecolor=CB_BLUE,
              label="centered (STadmix, default)"),
        Patch(facecolor=CB_ORANGE, alpha=0.65, edgecolor=CB_ORANGE,
              label="no-center (uncentered)"),
        Line2D([], [], ls="--", color=GREY, lw=1.8, label="identity ($\\hat{r}_g = r_g$)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
               frameon=False, fontsize=15, bbox_to_anchor=(0.5, -0.012))
    fig.suptitle("Supplementary Figure 14 — Genotype centering determines the estimand "
                 "across admixture cohorts\n"
                 "Centered STadmix recovers the truth under either basis; the uncentered "
                 "coding recovers it only when the truth is itself uncentered, otherwise "
                 "saturating toward 1.\n"
                 "(3-way cohorts pool all three ancestry pairs)",
                 x=0.012, ha="left", y=1.005, weight="bold", fontsize=17)
    fig.tight_layout(rect=(0, 0.02, 1, 0.99))
    savefig(fig, "SFig14_rg_centering_across_cohorts")


if __name__ == "__main__":
    sfig12()
    sfig13()
    sfig14()
    print(f"[13] done — outputs in {FIGS_DIR}")
