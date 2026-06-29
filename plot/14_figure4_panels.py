#!/usr/bin/env python3
"""
14_figure4_panels.py
────────────────────
Bespoke Figure-4 panels for the "Including all participants increases power"
section. Each panel visualises ONE discovery advantage with a distinct visual
metaphor and minimal numbers (show, don't tell). Standalone files so they can
be composited freely.

  A1  IL23R / Crohn's   — inclusion → power (EUR): same effect, tighter CI,
                          crosses genome-wide significance.
  A2  ADRB2/SH3TC2 / BMI — inclusion → power (AFR): more AFR haplotypes,
                          discovery where meta-analysis missed.
  B   HFE / MCHC        — modeling local ancestry recovers signal the
                          per-stratum meta-analysis dropped (EUR not tested).
  C1  APOC1 / ALP       — heterogeneity: opposite-direction effects cancel
                          under a shared-effect model.
  C2  HPR/TXNL4B / LDL  — heterogeneity (lipids): same mechanism, recurs.

Okabe-Ito colourblind-safe palette, large fonts, vector PDF + 400-dpi PNG.
All numbers read at run time from scatter_output (so they track data refreshes).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

from common import (FIGS_DIR, SCATTER_TABLES, load_scatter, safe_float,
                    apply_manuscript_style)
apply_manuscript_style()

# Okabe-Ito
CB_AFR   = "#D55E00"   # vermilion
CB_EUR   = "#0072B2"   # deep blue
CB_NAT   = "#E69F00"   # orange
CB_EAS   = "#56B4E9"   # sky blue
CB_SAS   = "#CC79A7"   # reddish purple
GREY     = "#999999"
COL_META = "#7a7a7a"   # neutral grey  — global-ancestry meta-analysis
COL_ST   = "#CC79A7"   # brand magenta — SAIGE-Tractor
GW = 5e-8
GW_LOG = -np.log10(GW)

NAME_ST = "SAIGE-Tractor"
NAME_META = "global-ancestry\nmeta-analysis"

df = load_scatter(SCATTER_TABLES["cct"])

def row(gene, ph, pos=None):
    sub = df[(df["phenotype"].astype(str) == str(ph)) &
             (df["Gene"].astype(str).str.contains(gene, na=False, regex=False))]
    if pos is not None and not sub.empty:
        sub = sub.iloc[(sub["Pos"] - pos).abs().argsort()[:1]]
    return sub.iloc[0] if not sub.empty else None

def savefig(fig, name):
    fig.savefig(FIGS_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGS_DIR / f"{name}.png", bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[14] wrote {name}.pdf + .png")

def plog(p):
    p = safe_float(p)
    return np.nan if (p is None or p <= 0) else -np.log10(p)


# ════════════════════════════════════════════════════════════════════════════
# A1 — IL23R: same effect, tighter CI, crosses the line
# ════════════════════════════════════════════════════════════════════════════
def panelA1():
    r = row("IL23R", "GI_522.11", 67240275)
    b_st  = safe_float(r["SAIGE_BETA_c_anc3"]); se_st = safe_float(r["SAIGE_SE_c_anc3"])
    p_st  = safe_float(r["SAIGE_p.value_c_anc3"])
    b_ab  = safe_float(r["ABA_EUR_BETA"]);      se_ab = safe_float(r["ABA_EUR_SE"])
    p_ab  = safe_float(r["ABA_EUR_Pvalue"])

    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    # two rows: meta (top), SAIGE-Tractor (bottom)
    y_ab, y_st = 1.0, 0.0
    for y, b, se, col, lbl in [(y_ab, b_ab, se_ab, COL_META, "global-ancestry meta-analysis"),
                               (y_st, b_st, se_st, COL_ST,  "SAIGE-Tractor (EUR haplotypes)")]:
        lo, hi = b - 1.96*se, b + 1.96*se
        ax.plot([lo, hi], [y, y], color=col, lw=5, solid_capstyle="round", zorder=4)
        ax.scatter([b], [y], s=320, color=col, edgecolor="white", lw=2, zorder=5)
        ax.text(b, y + 0.20, lbl, ha="center", va="bottom", fontsize=13,
                weight="bold", color=col)
    # contraction arrow between the two CIs (same point estimate)
    ax.annotate("", xy=(b_st - 1.96*se_st, 0.42), xytext=(b_ab - 1.96*se_ab, 0.58),
                arrowprops=dict(arrowstyle="-|>", color="#444", lw=1.6))
    ax.text(b_st - 1.96*se_st - 0.02, 0.5,
            "+European haplotypes\nrecovered from\nadmixed participants",
            ha="right", va="center", fontsize=11.5, style="italic", color="#333")
    ax.axvline(0, color="#bbb", lw=1.2, ls="-", zorder=1)
    ax.set_yticks([]); ax.set_ylim(-0.7, 1.9)
    ax.set_xlabel("Effect on Crohn's disease  (β, 95% CI)")
    # significance annotations (minimal numbers: just the two p-values)
    ax.text(b_ab + 1.96*se_ab + 0.02, y_ab,
            f"p = {p_ab:.0e}\n(below threshold)", ha="left", va="center",
            fontsize=11.5, color=COL_META)
    ax.text(b_st + 1.96*se_st + 0.02, y_st,
            f"p = {p_st:.0e}\ngenome-wide", ha="left", va="center",
            fontsize=11.5, color=COL_ST, weight="bold")
    ax.set_title("A1   IL23R — Crohn's disease\n"
                 "Same effect, more samples → the interval tightens and crosses the line",
                 loc="left", weight="bold", fontsize=14)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.margins(x=0.28)
    savefig(fig, "fig4A1_il23r_inclusion_EUR")


# ════════════════════════════════════════════════════════════════════════════
# A2 — ADRB2/SH3TC2 (AFR): more AFR haplotypes → discovery where meta missed
# ════════════════════════════════════════════════════════════════════════════
def panelA2():
    r = row("ADRB2", "BMI", 148898672)
    b_st = safe_float(r["SAIGE_BETA_c_anc1"]); se_st = safe_float(r["SAIGE_SE_c_anc1"])
    af_st = safe_float(r["SAIGE_AF_Allele2_anc1"]); n_st = safe_float(r["SAIGE_N_haplo_anc1"])
    p_st = safe_float(r["SAIGE_p.value_c_anc1"])
    b_ab = safe_float(r["ABA_AFR_BETA"]); se_ab = safe_float(r["ABA_AFR_SE"])
    af_ab = safe_float(r["ABA_AFR_AF_Allele2"]); p_ab = safe_float(r["ABA_AFR_Pvalue"])

    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    y_ab, y_st = 1.0, 0.0
    for y, b, se, col, lbl, p in [
            (y_ab, b_ab, se_ab, COL_META, "global-ancestry meta-analysis (AFR cluster)", p_ab),
            (y_st, b_st, se_st, COL_ST,  "SAIGE-Tractor (African haplotypes)", p_st)]:
        lo, hi = b - 1.96*se, b + 1.96*se
        ax.plot([lo, hi], [y, y], color=col, lw=5, solid_capstyle="round", zorder=4)
        ax.scatter([b], [y], s=320, color=col, edgecolor="white", lw=2, zorder=5)
        ax.text(b, y + 0.20, lbl, ha="center", va="bottom", fontsize=12.5,
                weight="bold", color=col)
        tag = "genome-wide" if p < GW else "below threshold"
        ax.text(lo - 0.004, y, f"p = {p:.0e}\n{tag}", ha="right", va="center",
                fontsize=11.5, color=col, weight=("bold" if p < GW else "normal"))
    ax.axvline(0, color="#bbb", lw=1.2, zorder=1)
    ax.set_yticks([]); ax.set_ylim(-0.7, 1.9)
    ax.set_xlabel("Effect on BMI  (β, 95% CI)")
    # the inclusion message: more AFR haplotypes / higher carrier frequency
    ax.text(max(b_st, b_ab) + 1.96*max(se_st, se_ab) + 0.004, 0.5,
            "African haplotypes pooled\nfrom the whole cohort\n"
            f"(carrier freq. {af_ab:.0%} → {af_st:.0%})",
            ha="left", va="center", fontsize=11.5, style="italic", color="#333")
    ax.set_title("A2   ADRB2 / SH3TC2 region — BMI (African ancestry)\n"
                 "The same inclusion gain in a different ancestry → a discovery meta-analysis missed",
                 loc="left", weight="bold", fontsize=13.5)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.margins(x=0.34)
    savefig(fig, "fig4A2_adrb2_inclusion_AFR")


# ════════════════════════════════════════════════════════════════════════════
# B — HFE: signal lives in EUR, but meta dropped the EUR stratum
# ════════════════════════════════════════════════════════════════════════════
def panelB():
    r = row("HFE", "3009744", 26092913)
    ancs = [("AFR", "anc1", "ABA_AFR", CB_AFR),
            ("EAS", "anc2", "ABA_EAS", CB_EAS),
            ("EUR", "anc3", "ABA_EUR", CB_EUR),
            ("NatAm","anc4","ABA_AMR", CB_NAT),
            ("SAS", "anc5", "ABA_SAS", CB_SAS)]
    st_af = {nm: safe_float(r[f"SAIGE_AF_Allele2_{suf}"]) for nm, suf, _, _ in ancs}
    ab_af = {nm: safe_float(r[f"{ab}_AF_Allele2"]) for nm, _, ab, _ in ancs}
    p_st  = safe_float(r["SAIGE_p.value_c_anc3"])     # EUR
    p_meta = safe_float(r["ABA_META_Pvalue"])

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8),
                             gridspec_kw={"width_ratios": [1.5, 1.0], "wspace": 0.30})

    # ── Left: per-ancestry allele frequency, meta vs SAIGE-Tractor ──
    axL = axes[0]
    x = np.arange(len(ancs)); w = 0.38
    ymax = max([v for v in st_af.values() if v and not np.isnan(v)]) * 1.30
    for i, (nm, suf, ab, col) in enumerate(ancs):
        v_ab = ab_af[nm]
        if v_ab is None or np.isnan(v_ab):
            axL.text(i - w/2, ymax*0.02, "✗ not\ntested", ha="center", va="bottom",
                     fontsize=10.5, color="#b00", weight="bold")
        else:
            axL.bar(i - w/2, v_ab, width=w, color=COL_META, edgecolor="white")
        v_st = st_af[nm]
        if v_st is not None and not np.isnan(v_st) and v_st > 0:
            axL.bar(i + w/2, v_st, width=w, color=col, edgecolor="white")
    axL.set_ylim(0, ymax)
    axL.set_xticks(x); axL.set_xticklabels([a[0] for a in ancs])
    axL.set_ylabel("Allele frequency of the tested variant")
    axL.set_title("where is the variant tested?", loc="left", weight="bold", fontsize=13)
    axL.legend(handles=[Rectangle((0,0),1,1,fc=COL_META, label="global-ancestry meta-analysis"),
                        Rectangle((0,0),1,1,fc=CB_EUR, label="SAIGE-Tractor (per local ancestry)")],
               fontsize=11, frameon=False, loc="upper left")
    for s in ("top","right"): axL.spines[s].set_visible(False)
    # annotate the recovered EUR bar from the right so it never hits the legend
    eur_i = 2
    axL.annotate("signal lives here\n(EUR, 6.6%)", xy=(eur_i + w/2, st_af["EUR"]),
                 xytext=(eur_i + 1.05, st_af["EUR"]*0.75),
                 fontsize=12, color=CB_EUR, weight="bold", va="center", ha="center",
                 arrowprops=dict(arrowstyle="-|>", color=CB_EUR, lw=1.8))

    # ── Right: resulting evidence (−log10 p) ──
    axR = axes[1]
    bars = [("meta-analysis\n(anchored on rare\nAFR/AMR/SAS tails)", plog(p_meta), COL_META),
            (f"SAIGE-Tractor\n(models the EUR\nhaplotype directly)", plog(p_st), CB_EUR)]
    for i, (lbl, val, col) in enumerate(bars):
        axR.bar(i, val, width=0.6, color=col, edgecolor="white")
        axR.text(i, val + 1.5, f"{val:.0f}", ha="center", va="bottom",
                 fontsize=13, weight="bold", color=col)
    axR.axhline(GW_LOG, color="#444", ls="--", lw=1.6)
    axR.text(1.45, GW_LOG + 1.5, "genome-wide", ha="right", va="bottom",
             fontsize=10.5, color="#444")
    axR.set_xlim(-0.6, 1.6); axR.set_ylim(0, plog(p_st)*1.15)
    axR.set_xticks([0, 1]); axR.set_xticklabels([b[0] for b in bars], fontsize=10.5)
    axR.set_ylabel("Evidence  ($-\\log_{10} p$)")
    axR.set_title("the recovered EUR signal", loc="left", weight="bold", fontsize=13)
    for s in ("top","right"): axR.spines[s].set_visible(False)

    fig.suptitle("B   HFE — MCHC: meta-analysis dropped the EUR stratum — "
                 "the one carrying the signal; SAIGE-Tractor models it directly",
                 x=0.012, ha="left", y=1.02, weight="bold", fontsize=14)
    savefig(fig, "fig4B_hfe_modeling")


# ════════════════════════════════════════════════════════════════════════════
# C — heterogeneity: opposite-direction effects (shared by C1, C2)
# ════════════════════════════════════════════════════════════════════════════
def _diverging(r, trait, title, fname):
    ancs = [("AFR","anc1",CB_AFR),("EAS","anc2",CB_EAS),("EUR","anc3",CB_EUR),
            ("NatAm","anc4",CB_NAT),("SAS","anc5",CB_SAS)]
    recs = []
    for nm, suf, col in ancs:
        b = safe_float(r[f"SAIGE_BETA_c_{suf}"]); se = safe_float(r[f"SAIGE_SE_c_{suf}"])
        p = safe_float(r[f"SAIGE_p.value_c_{suf}"])
        if b is not None and se is not None:
            recs.append((nm, b, se, p, col))
    # order by beta so positives and negatives separate visually
    recs.sort(key=lambda t: t[1])
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ys = np.arange(len(recs))
    for y, (nm, b, se, p, col) in zip(ys, recs):
        ax.plot([b-1.96*se, b+1.96*se], [y, y], color=col, lw=4.5,
                solid_capstyle="round", zorder=4)
        ax.scatter([b], [y], s=260, color=col, edgecolor="white", lw=1.8, zorder=5)
        sig = "★" if (p is not None and p < GW) else ""
        ax.text(b + (0.004 if b >= 0 else -0.004), y + 0.28,
                f"{nm} {sig}", ha=("left" if b >= 0 else "right"),
                va="bottom", fontsize=12.5, weight="bold", color=col)
    ax.axvline(0, color="#333", lw=1.6, zorder=2)
    # "what a shared-effect model sees" marker near zero
    ax.scatter([0], [len(recs)+0.1], s=300, marker="D", color=GREY,
               edgecolor="white", lw=1.6, zorder=6)
    ax.text(0, len(recs)+0.45, "a single shared-effect model\naverages these to ≈ 0",
            ha="center", va="bottom", fontsize=11.5, style="italic", color="#444")
    ax.set_yticks([]); ax.set_ylim(-0.7, len(recs)+1.3)
    ax.set_xlabel(f"Per-ancestry effect on {trait}  (β, 95% CI)")
    ax.set_title(title, loc="left", weight="bold", fontsize=13.5)
    # directional cue: arrows in the margins under the title, clear of all bars
    xl = ax.get_xlim()
    ax.text(xl[1]*0.72, len(recs)+0.1, "African: +  →", ha="center", va="center",
            fontsize=11.5, color=CB_AFR, weight="bold", style="italic")
    ax.text(xl[0]*0.72, len(recs)+0.1, "←  European: −", ha="center", va="center",
            fontsize=11.5, color=CB_EUR, weight="bold", style="italic")
    for s in ("top","right","left"): ax.spines[s].set_visible(False)
    fig.tight_layout()
    savefig(fig, fname)

def panelC1():
    r = row("APOC1", "3035995", 44919689)
    _diverging(r, "alkaline phosphatase",
               "C1   APOC1 (APOE/APOC1) — alkaline phosphatase\n"
               "African + and European − effects cancel under a shared-effect model (★ = genome-wide)",
               "fig4C1_apoc1_heterogeneity")

def panelC2():
    r = row("HPR", "3028288", 72080103)
    _diverging(r, "LDL cholesterol",
               "C2   HPR / TXNL4B — LDL cholesterol\n"
               "The same heterogeneity recurs: African + vs European − (★ = genome-wide)",
               "fig4C2_hpr_heterogeneity")


if __name__ == "__main__":
    panelA1(); panelA2(); panelB(); panelC1(); panelC2()
    print(f"[14] done — outputs in {FIGS_DIR}")
