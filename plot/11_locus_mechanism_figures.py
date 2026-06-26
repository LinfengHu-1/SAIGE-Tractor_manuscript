#!/usr/bin/env python3
"""
11_locus_mechanism_figures.py
─────────────────────────────
Manuscript figures that visualise the **mechanisms** by which SAIGE-Tractor
outperforms All by All meta-analysis at specific case-study loci.

We use the canonical forest-plot layout common in Nature/Nature Genetics tier
manuscripts: per-ancestry β±95%CI from both methods stacked vertically, with
the combined-test summaries (HET, HOM, CCT) appended at the bottom and the
All by All meta-analysis estimate added for direct visual comparison. A side
table lists exact β, SE, p, and AF for each row.

Figures produced (one per mechanism):
  Fig 14a  LA-refinement                — VEGFA / triglycerides
  Fig 14b  Effective-N rescue           — IL23R / Crohn's
  Fig 14c  Single-ancestry signal       — RCOR1 / platelet count
  Fig 15   Cross-ancestry heterogeneity — HPR/TXNL4B / LDL  (opposing β)
  Fig 16   Global-MAF filtering         — HLA-DRA / multiple sclerosis
  Fig 17   LA-conditioning unmasking    — OR51B5 / sickle-cell anemia

Plus a 2x3 "mechanism overview" panel combining all six.
"""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from common import (
    FIGS_DIR, ANCS, PHENO_LABELS,
    SCATTER_TABLES, load_scatter, safe_float,
    apply_manuscript_style,
    NAME_SAIGE, NAME_ABA, NAME_ABA_META, NAME_SAIGE_MEGA,
    COLOR_SAIGE, COLOR_ABA, ANC_COLORS,
)
apply_manuscript_style()

df_cct = load_scatter(SCATTER_TABLES["cct"])
# Also load the per-ancestry tables for fallback lookups (e.g. RCOR1 is only
# flagged in the AMR per-ancestry scatter, not in the CCT-vs-META table).
_ALL_TABLES = {k: load_scatter(SCATTER_TABLES[k]) for k in
                ["cct","hom","het","anc_afr","anc_eas","anc_eur","anc_amr","anc_sas"]}

# Map ancestry name → (anc suffix in SAIGE, ABA per-ancestry column prefix)
ANC_KEYS = [
    ("AFR",   "anc1", "ABA_AFR"),
    ("EAS",   "anc2", "ABA_EAS"),
    ("EUR",   "anc3", "ABA_EUR"),
    ("NatAm", "anc4", "ABA_AMR"),
    ("SAS",   "anc5", "ABA_SAS"),
]

def savefig(fig, name):
    fig.savefig(FIGS_DIR / f"{name}.pdf")
    fig.savefig(FIGS_DIR / f"{name}.png")
    plt.close(fig)
    print(f"[11] wrote {name}.pdf + .png")

def find_locus(gene_pattern, phenotype, chrom=None, pos=None):
    """Return the scatter row for a (gene, phenotype) locus, searching CCT
    first then falling back to HOM, HET, and the 5 per-ancestry tables. If
    multiple rows match, return the one closest to (chrom, pos)."""
    for key in ["cct","hom","het","anc_afr","anc_eas","anc_eur","anc_amr","anc_sas"]:
        d = _ALL_TABLES[key]
        hits = d[d["Gene"].astype(str).str.contains(gene_pattern, na=False)
                 & (d["phenotype"].astype(str) == str(phenotype))]
        if hits.empty: continue
        if chrom is not None and pos is not None:
            hits = hits.copy()
            hits["_d"] = abs(pd.to_numeric(hits["Pos"], errors="coerce") - pos)
            hits = hits.sort_values("_d")
        return hits.iloc[0]
    return None

def collect_per_anc(r):
    """Pull per-ancestry β/SE/p/AF/N_haplo for both methods from a row."""
    out = []
    for name, suf, aba in ANC_KEYS:
        rec = dict(
            ancestry = name,
            saige_beta = safe_float(r.get(f"SAIGE_BETA_c_{suf}")),
            saige_se   = safe_float(r.get(f"SAIGE_SE_c_{suf}")),
            saige_p    = safe_float(r.get(f"SAIGE_p.value_c_{suf}")),
            saige_af   = safe_float(r.get(f"SAIGE_AF_Allele2_{suf}")),
            saige_nhap = safe_float(r.get(f"SAIGE_N_haplo_{suf}")),
            aba_beta   = safe_float(r.get(f"{aba}_BETA")),
            aba_se     = safe_float(r.get(f"{aba}_SE")),
            aba_p      = safe_float(r.get(f"{aba}_Pvalue")),
            aba_af     = safe_float(r.get(f"{aba}_AF_Allele2")),
        )
        out.append(rec)
    return out

def pfmt(p):
    v = safe_float(p)
    if v is None: return "—"
    if v == 0: return "<1e-320"
    if v < 0.001: return f"{v:.2e}"
    return f"{v:.3g}"

def bfmt(b):
    v = safe_float(b)
    return "—" if v is None else f"{v:+.3f}"

# ═══════════════════════════════════════════════════════════════════════════
# Core forest-plot helper
# ═══════════════════════════════════════════════════════════════════════════
def forest_plot(r, title, mechanism_label, mechanism_color="#1a237e",
                point_size_scale=300, show_aba_meta=True, x_pad=0.10,
                figsize=(13.5, 8.5), filename=None,
                summary_only_for=None):
    """Forest plot of per-ancestry β±95%CI for one variant.

    Parameters
    ----------
    r : pd.Series
        A row from the CCT-vs-META scatter table.
    title : str
        Top-line title (gene — trait — chr:pos).
    mechanism_label : str
        Short subtitle describing the discovery mechanism.
    summary_only_for : list[str] or None
        If given, only show SAIGE-Tractor combined summaries (CCT/HOM/HET) for
        these labels (default: all).
    """
    per_anc = collect_per_anc(r)
    # Build display rows
    display_rows = []
    for rec in per_anc:
        label = rec["ancestry"]
        # SAIGE side
        display_rows.append(dict(
            label=label, method=NAME_SAIGE,
            beta=rec["saige_beta"], se=rec["saige_se"], p=rec["saige_p"],
            af=rec["saige_af"], n=rec["saige_nhap"], n_lbl="N$_{hap}$",
            color=ANC_COLORS.get(label,"#444"),
        ))
        # ABA side
        display_rows.append(dict(
            label=label, method=NAME_ABA,
            beta=rec["aba_beta"], se=rec["aba_se"], p=rec["aba_p"],
            af=rec["aba_af"], n=None, n_lbl="",
            color=ANC_COLORS.get(label,"#444"),
        ))

    # Combined summaries (SAIGE-Tractor only)
    summaries = [
        ("HOM (homogeneous)", safe_float(r.get("SAIGE_P_hom_admixed_c")),
         safe_float(r.get("SAIGE_BETA_c_ancALL")), safe_float(r.get("SAIGE_SE_c_ancALL")), "#1a237e"),
        ("HET (heterogeneous)", safe_float(r.get("SAIGE_P_het_admixed_c")),
         None, None, "#e65100"),
        ("CCT (combined)", safe_float(r.get("SAIGE_P_cct_admixed_c")),
         safe_float(r.get("SAIGE_BETA_c_ancALL")), safe_float(r.get("SAIGE_SE_c_ancALL")), "#b71c1c"),
    ]
    if show_aba_meta:
        summaries.append((f"{NAME_ABA} META", safe_float(r.get("ABA_pvalue")),
                          safe_float(r.get("ABA_META_BETA")), safe_float(r.get("ABA_META_SE")),
                          COLOR_ABA))

    # ── Compute β±95%CI ranges for axis scaling ──
    betas, lo, hi = [], [], []
    for d in display_rows:
        b, se = d["beta"], d["se"]
        if b is not None and se is not None and se > 0:
            betas.append(b); lo.append(b-1.96*se); hi.append(b+1.96*se)
    for _, _, b, se, _ in summaries:
        if b is not None and se is not None and se > 0:
            betas.append(b); lo.append(b-1.96*se); hi.append(b+1.96*se)
    if not betas:
        print(f"[11]  forest_plot: nothing to plot for {title}")
        return
    x_min = min(lo) - x_pad * max(1e-3, max(hi)-min(lo))
    x_max = max(hi) + x_pad * max(1e-3, max(hi)-min(lo))

    # ── Layout: forest column + side table column ───────────────────────────
    fig = plt.figure(figsize=figsize)
    gs  = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.10)
    ax  = fig.add_subplot(gs[0, 0])
    axT = fig.add_subplot(gs[0, 1]); axT.axis("off")

    # y positions: per-ancestry rows (SAIGE then ABA pair) then a divider then summaries
    n_per_anc = len(display_rows)
    n_sum     = len(summaries)
    # Add gap between per-anc and summaries
    gap = 1.0
    n_total = n_per_anc + n_sum + 1
    y_positions = list(range(n_total, 0, -1))
    y_per_anc = y_positions[:n_per_anc]
    y_sum     = y_positions[n_per_anc+1 : n_per_anc+1+n_sum]

    # ── Plot per-ancestry rows ──
    for y, d in zip(y_per_anc, display_rows):
        b, se = d["beta"], d["se"]
        if b is None or se is None or se <= 0:
            # No estimate available — draw "NA" tick
            ax.text(x_min, y, "  not tested",
                    ha="left", va="center", color="#999",
                    fontsize=12, style="italic")
            continue
        ci_lo, ci_hi = b-1.96*se, b+1.96*se
        marker = "o" if d["method"] == NAME_SAIGE else "s"
        ax.plot([ci_lo, ci_hi], [y, y],
                color=d["color"], lw=2.4, solid_capstyle="round", alpha=0.9)
        ax.scatter([b], [y], s=point_size_scale, color=d["color"],
                   marker=marker, edgecolor="white", linewidth=1.5, zorder=10)

    # Divider line between per-ancestry and summary rows
    div_y = n_sum + 0.5 + 0.5
    ax.axhline(div_y, color="#bbb", ls="--", lw=1.0)

    # ── Plot SAIGE-Tractor combined summaries ──
    for y, (lbl, p, b, se, col) in zip(y_sum, summaries):
        if b is None or se is None or se <= 0:
            # HET doesn't have a single β/SE — show only the p-value bar
            ax.text(0, y, f"  p = {pfmt(p)}",
                    ha="left", va="center", fontsize=13,
                    color=col, weight="bold")
            continue
        ci_lo, ci_hi = b-1.96*se, b+1.96*se
        ax.plot([ci_lo, ci_hi], [y, y], color=col, lw=3.5,
                solid_capstyle="round")
        ax.scatter([b], [y], s=point_size_scale*1.2, color=col,
                   marker="D", edgecolor="white", linewidth=1.5, zorder=10)

    # ── Y-axis labels ──
    yticks = []; yticklabels = []
    for y, d in zip(y_per_anc, display_rows):
        method_short = "SAIGE-T" if d["method"] == NAME_SAIGE else "ABA"
        yticks.append(y); yticklabels.append(f"{d['label']:6s}   {method_short}")
    for y, (lbl, _, _, _, _) in zip(y_sum, summaries):
        yticks.append(y); yticklabels.append(lbl)
    ax.set_yticks(yticks); ax.set_yticklabels(yticklabels, fontsize=12)
    ax.axvline(0, color="#999", lw=1)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0.4, n_total + 0.6)
    ax.set_xlabel(f"Effect-size β (95% CI)")
    # Reserve extra top space (~25% of axes height) for title + banner without overlap
    fig.suptitle(title, x=0.04, ha="left", y=0.98, weight="bold", fontsize=17)
    # Subtitle banner with mechanism — placed BELOW suptitle, above the axes
    fig.text(0.04, 0.93, "  " + mechanism_label + "  ",
             ha="left", va="center",
             fontsize=12, color="white", weight="bold",
             bbox=dict(facecolor=mechanism_color, edgecolor="none", pad=6))

    # Custom legend below the axes (so it never overlaps any row)
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([],[], marker="o", color="w", markerfacecolor="#666",
               markersize=12, markeredgecolor="white", label=f"{NAME_SAIGE} (per-ancestry)"),
        Line2D([],[], marker="s", color="w", markerfacecolor="#666",
               markersize=12, markeredgecolor="white", label=f"{NAME_ABA} (per-ancestry)"),
        Line2D([],[], marker="D", color="w", markerfacecolor="#444",
               markersize=12, markeredgecolor="white", label="Combined / summary estimate"),
    ]
    ax.legend(handles=legend_handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.10), ncol=3, fontsize=11,
              frameon=False)
    fig.subplots_adjust(top=0.88, bottom=0.16)

    # ── Side table ──
    axT.set_xlim(0, 1); axT.set_ylim(0.4, n_total + 0.6)
    headers = ["", "β", "SE", "p", "AF"]
    col_x   = [0.00, 0.36, 0.52, 0.66, 0.88]
    for x, h in zip(col_x, headers):
        axT.text(x, n_total + 0.4, h, fontsize=12, weight="bold", color="#1a237e")
    for y, d in zip(y_per_anc, display_rows):
        axT.text(col_x[1], y, bfmt(d["beta"]),    fontsize=11, va="center")
        axT.text(col_x[2], y, bfmt(d["se"]),      fontsize=11, va="center")
        axT.text(col_x[3], y, pfmt(d["p"]),       fontsize=11, va="center", family="monospace")
        af_s = "—" if d["af"] is None else f"{d['af']:.3f}"
        axT.text(col_x[4], y, af_s,               fontsize=11, va="center")
    for y, (lbl, p, b, se, col) in zip(y_sum, summaries):
        axT.text(col_x[1], y, bfmt(b),  fontsize=11, va="center", weight="bold", color=col)
        axT.text(col_x[2], y, bfmt(se), fontsize=11, va="center", color=col)
        axT.text(col_x[3], y, pfmt(p),  fontsize=11, va="center", family="monospace", weight="bold", color=col)
    if filename:
        savefig(fig, filename)


# ═══════════════════════════════════════════════════════════════════════════
# Build the per-mechanism figures
# ═══════════════════════════════════════════════════════════════════════════
def fig14a_VEGFA():
    r = find_locus("VEGFA", "3022192", chrom=6, pos=43791136)
    if r is None: print("[11] VEGFA not found"); return
    forest_plot(
        r,
        title=f"VEGFA — triglycerides (chr6:43,791,136)",
        mechanism_label=f"MECHANISM 1 · Local-ancestry refinement: single-ancestry signal diluted in {NAME_ABA} global meta",
        mechanism_color="#2e7d32",
        filename="fig14a_mechanism_VEGFA_triglycerides",
    )

def fig14b_IL23R():
    r = find_locus("IL23R", "GI_522.11", chrom=1, pos=67240275)
    if r is None: print("[11] IL23R not found"); return
    forest_plot(
        r,
        title=f"IL23R — Crohn's disease (chr1:67,240,275)",
        mechanism_label=f"MECHANISM 2 · Effective-N rescue: {NAME_SAIGE} recovers EUR haplotypes from admixed individuals",
        mechanism_color="#0d47a1",
        filename="fig14b_mechanism_IL23R_crohns",
    )

def fig14c_RCOR1():
    r = find_locus("RCOR1", "3024929", chrom=14, pos=102663104)
    if r is None: print("[11] RCOR1 not found"); return
    forest_plot(
        r,
        title=f"RCOR1 — platelet count (chr14:102,663,104)",
        mechanism_label=f"MECHANISM 3 · Single-ancestry signal: NatAm-specific effect at high local AF, invisible to combined CCT",
        mechanism_color="#ef6c00",
        filename="fig14c_mechanism_RCOR1_platelet",
    )

def fig15_HPR():
    r = find_locus("HPR;TXNL4B", "3028288", chrom=16, pos=72080103)
    if r is None: print("[11] HPR/TXNL4B not found"); return
    forest_plot(
        r,
        title=f"HPR / TXNL4B — LDL cholesterol (chr16:72,080,103)",
        mechanism_label=f"MECHANISM 4 · Cross-ancestry heterogeneity: opposite-direction β (AFR + vs EUR −) cancels in HOM/META, captured by HET/CCT",
        mechanism_color="#b71c1c",
        figsize=(14.5, 9),
        filename="fig15_mechanism_HPR_LDL",
    )

def fig16_HLA_DRA():
    r = find_locus("HLA-DRA", "NS_326.1", chrom=6, pos=32445768)
    if r is None: print("[11] HLA-DRA not found"); return
    forest_plot(
        r,
        title=f"HLA-DRA — multiple sclerosis (chr6:32,445,768)",
        mechanism_label=f"MECHANISM 5 · Variant filtered out of {NAME_ABA} per-ancestry strata; recovered per-haplotype by {NAME_SAIGE}",
        mechanism_color="#4a148c",
        figsize=(14.5, 9),
        filename="fig16_mechanism_HLA_DRA_MS",
    )

def fig17_OR51B5():
    r = find_locus("OR51B5", "282.5", chrom=11, pos=5498240)
    if r is None: print("[11] OR51B5 not found"); return
    forest_plot(
        r,
        title=f"OR51B5 / HBB locus — sickle-cell anemia (chr11:5,498,240)",
        mechanism_label=f"MECHANISM 6 · Local-ancestry conditioning unmasks the AFR-specific Hb-locus signal at much greater effective N",
        mechanism_color="#37474f",
        figsize=(14.5, 9),
        filename="fig17_mechanism_OR51B5_SCA",
    )

# ═══════════════════════════════════════════════════════════════════════════
# Overview panel: 6 mechanisms side-by-side mini-forests
# ═══════════════════════════════════════════════════════════════════════════
def fig_mechanism_overview():
    """A 2x3 grid summarising all six mechanisms in compact forest format.
    Useful as the main-text 'mechanism summary' figure of the manuscript."""
    cases = [
        ("VEGFA", "3022192", 6, 43791136,
         "Single-ancestry signal diluted in meta",  "#2e7d32"),
        ("IL23R", "GI_522.11", 1, 67240275,
         "Effective-N rescue in admixed cohort",    "#0d47a1"),
        ("RCOR1", "3024929", 14, 102663104,
         "NatAm-specific signal invisible to CCT",  "#ef6c00"),
        ("HPR;TXNL4B", "3028288", 16, 72080103,
         "Opposing β across ancestries (HET advantage)", "#b71c1c"),
        ("HLA-DRA", "NS_326.1", 6, 32445768,
         "Variant filtered out of per-ancestry META", "#4a148c"),
        ("OR51B5", "282.5", 11, 5498240,
         "LA-conditioning unmasks AFR-specific Hb signal", "#37474f"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(22, 16))
    for ax, (gene, ph, chrom, pos, mech, col) in zip(axes.flatten(), cases):
        r = find_locus(gene, ph, chrom=chrom, pos=pos)
        if r is None:
            ax.axis("off"); continue
        per_anc = collect_per_anc(r)
        # Mini forest: per-ancestry SAIGE β±95%CI, with All by All META overlaid
        y_pos = np.arange(len(per_anc), 0, -1)
        all_betas = []
        for y, d in zip(y_pos, per_anc):
            if d["saige_beta"] is not None and d["saige_se"] and d["saige_se"] > 0:
                lo = d["saige_beta"] - 1.96*d["saige_se"]
                hi = d["saige_beta"] + 1.96*d["saige_se"]
                ax.plot([lo, hi], [y, y],
                        color=ANC_COLORS.get(d["ancestry"], "#444"), lw=2.4,
                        solid_capstyle="round")
                ax.scatter([d["saige_beta"]], [y], s=160,
                           color=ANC_COLORS.get(d["ancestry"], "#444"),
                           marker="o", edgecolor="white", linewidth=1.2, zorder=10)
                all_betas += [lo, hi]
            else:
                ax.text(0.02, y, "n.t.", color="#999", fontsize=10,
                        va="center", style="italic", transform=ax.get_yaxis_transform())
        # Overlay All by All META as a horizontal short marker on the diagonal
        b_meta = safe_float(r.get("ABA_META_BETA"))
        se_meta= safe_float(r.get("ABA_META_SE"))
        cct_p  = safe_float(r.get("SAIGE_P_cct_admixed_c"))
        meta_p = safe_float(r.get("ABA_pvalue"))
        if b_meta is not None and se_meta and se_meta > 0:
            ax.axvline(b_meta, ls=":", color=COLOR_ABA, lw=2.5,
                       label=f"{NAME_ABA_META} β = {bfmt(b_meta)}")
            all_betas += [b_meta - 1.96*se_meta, b_meta + 1.96*se_meta]
        ax.axvline(0, color="#bbb", lw=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([d["ancestry"] for d in per_anc], fontsize=12, weight="bold")
        ax.set_ylim(0.4, len(per_anc)+0.6)
        ax.set_xlabel("β (95% CI)", fontsize=12)
        if all_betas:
            x_pad = max(0.05, 0.10*(max(all_betas)-min(all_betas)))
            ax.set_xlim(min(all_betas)-x_pad, max(all_betas)+x_pad)
        title_short = gene.replace(";", "/")
        phenotype_lbl = PHENO_LABELS.get(str(ph), str(ph))
        ax.set_title(f"{title_short} — {phenotype_lbl}", loc="left",
                     weight="bold", fontsize=14)
        # Mechanism subtitle and p-values placed BELOW the x-axis with room
        ax.text(0.02, -0.36,
                f"{mech}",
                transform=ax.transAxes, fontsize=12,
                color=col, va="top", ha="left", weight="bold")
        ax.text(0.02, -0.50,
                f"{NAME_SAIGE} CCT p = {pfmt(cct_p)}      "
                f"{NAME_ABA} meta-analysis p = {pfmt(meta_p)}",
                transform=ax.transAxes, fontsize=11,
                color="#333", va="top", ha="left", family="monospace")
        ax.legend(loc="lower right", fontsize=10, framealpha=0.95)
    fig.suptitle(f"Discovery mechanisms of {NAME_SAIGE} — per-ancestry effect estimates at six case-study loci",
                 x=0.04, ha="left", y=0.99, weight="bold", fontsize=18)
    fig.subplots_adjust(hspace=0.95, wspace=0.30, top=0.93, bottom=0.08, left=0.06, right=0.97)
    savefig(fig, "fig_mechanism_overview_6loci")

if __name__ == "__main__":
    fig14a_VEGFA()
    fig14b_IL23R()
    fig14c_RCOR1()
    fig15_HPR()
    fig16_HLA_DRA()
    fig17_OR51B5()
    fig_mechanism_overview()
    print(f"[11] done — outputs in {FIGS_DIR}")
