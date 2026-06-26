#!/usr/bin/env python3
"""
05_extra_figures.py
───────────────────
Additional manuscript figures requested after the first pass:

  Fig 9   meta_anchored_scatter         ABA META tophits × 3 panels (HOM/CCT/HET p)
  Fig 10  bmi_haplotype_bars            BMI ancestry haplotype counts —
                                         local (SAIGE-T median N_haplo) vs
                                         global (ABA N × 2), MID included
  Fig 11  method_overlap                Locus overlap between META, CCT, HOM, HET
                                         (Venn3 for HOM/CCT/META + UpSet-like bar
                                          for all 4 sets) and near-miss breakdown
                                          for the unique-locus categories.
"""
from __future__ import annotations
import math, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle, Circle

from common import (
    FIGS_DIR, REPLICATION_DIR, SAIGE_SS, ANCS, PHENO_LABELS,
    SCATTER_TABLES, load_scatter, load_aba_pheno_csv, safe_float,
    load_saige_chr, open_saige_sumstat,
    apply_manuscript_style, NAME_SAIGE, NAME_ABA, NAME_ABA_META, NAME_SAIGE_MEGA,
    COLOR_SAIGE, COLOR_ABA, COLOR_SHARED, ANC_COLORS as _ANC_COLORS,
)

apply_manuscript_style()

def savefig(fig, name):
    fig.savefig(FIGS_DIR / f"{name}.pdf")
    fig.savefig(FIGS_DIR / f"{name}.png")
    plt.close(fig)
    print(f"[05] wrote {name}.pdf + .png")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 9 — META-anchored scatter (3 combiners)
# ═══════════════════════════════════════════════════════════════════════════
def _meta_anchored_panel(ax, df, y_col, y_label, color_meta, color_saige):
    """Plot every All by All META tophit (locus_status in {shared, ABA_only})
       against the SAIGE-Tractor combiner p-value at that variant."""
    sub = df[df["tophit_source"]=="ABA"].copy()
    sub["x"] = -np.log10(pd.to_numeric(sub["ABA_pvalue"], errors="coerce").clip(1e-320))
    sub["y"] = -np.log10(pd.to_numeric(sub[y_col],        errors="coerce").clip(1e-320))
    sub = sub.dropna(subset=["x","y"])
    if sub.empty:
        ax.text(0.5,0.5,"no data",transform=ax.transAxes,ha="center"); return
    maxv = max(sub["x"].max(), sub["y"].max(), 9.0) * 1.05
    diff = sub["y"] - sub["x"]
    colors = np.where(diff > 1, color_saige, np.where(diff < -1, color_meta, "#666"))
    ax.scatter(sub["x"], sub["y"], c=colors, s=42, alpha=0.75, edgecolor="white", linewidth=0.4)
    ax.plot([0,maxv],[0,maxv], ls="--", c="#888", lw=1.2)
    ax.axvline(7.3, c="#c44", lw=1.0, alpha=0.5)
    ax.axhline(7.3, c="#c44", lw=1.0, alpha=0.5)
    ax.set_xlim(0, maxv); ax.set_ylim(0, maxv)
    ax.set_xlabel(f"{NAME_ABA_META}  $-\\log_{{10}}(p)$")
    ax.set_ylabel(y_label)
    ax.set_aspect("equal")
    n_s = int((diff >  1).sum()); n_m = int((diff < -1).sum())
    ax.text(0.04, 0.97,
            f"{NAME_SAIGE} stronger: {n_s}\n{NAME_ABA} stronger:   {n_m}\nN top hits:           {len(sub)}",
            transform=ax.transAxes, va="top", fontsize=12, family="monospace",
            bbox=dict(facecolor="white", edgecolor="#cfd8dc", alpha=0.95, pad=6))

def fig9_meta_anchored_scatter():
    df_cct = load_scatter(SCATTER_TABLES["cct"])
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    _meta_anchored_panel(axes[0], df_cct, "SAIGE_P_hom_admixed_c",
                         f"{NAME_SAIGE_MEGA}   HOM  $-\\log_{{10}}(p)$",
                         COLOR_ABA, COLOR_SAIGE)
    axes[0].set_title(f"a   HOM vs All by All",
                       loc="left", weight="bold")
    _meta_anchored_panel(axes[1], df_cct, "SAIGE_P_het_admixed_c",
                         f"{NAME_SAIGE}   HET  $-\\log_{{10}}(p)$",
                         COLOR_ABA, COLOR_SAIGE)
    axes[1].set_title(f"b   HET vs All by All",
                       loc="left", weight="bold")
    _meta_anchored_panel(axes[2], df_cct, "SAIGE_P_cct_admixed_c",
                         f"{NAME_SAIGE}   CCT  $-\\log_{{10}}(p)$",
                         COLOR_ABA, COLOR_SAIGE)
    axes[2].set_title(f"c   CCT (combined) vs All by All",
                       loc="left", weight="bold")
    fig.suptitle(f"All by All meta-analysis top hits comparison with {NAME_SAIGE} (CCT, HOM, HET)",
                 x=0.06, ha="left", y=1.02, weight="bold")
    savefig(fig, "fig9_meta_anchored_scatter")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 10 — BMI haplotype bars
# ═══════════════════════════════════════════════════════════════════════════
def fig10_bmi_haplotype_bars():
    """BMI haplotype counts per ancestry: local (SAIGE-Tractor median N_haplo
    per tested variant) vs global (All by All BMI sample size × 2 = haplotype
    equivalent). Manuscript-style: large fonts, no label overlap, fully spelled-out
    method names."""
    path = SAIGE_SS / "merged_saigetractor_BMI.txt.gz"
    if not path.exists():
        print("[05] BMI SAIGE-Tractor sumstat missing — skipping fig10")
        return
    cols = [f"N_haplo_{s}" for s,_,_,_ in ANCS]
    rdr = pd.read_csv(path, sep="\t", compression="gzip",
                      usecols=cols, chunksize=500_000, low_memory=False)
    medians = {c: [] for c in cols}
    for ch in rdr:
        for c in cols:
            v = pd.to_numeric(ch[c], errors="coerce")
            v = v[v > 0]
            if not v.empty: medians[c].append(v.median())
    saige_med = {c: (float(np.median(v)) if v else None) for c, v in medians.items()}

    csv = load_aba_pheno_csv()
    bmi = csv[csv["phenoname"]=="BMI"].set_index("ancestry")

    ancestries = ["AFR","EAS","EUR","NatAm","SAS","MID"]
    aba_keys   = ["AFR","EAS","EUR","AMR","SAS","MID"]
    saige_keys = [f"N_haplo_anc{i+1}" for i in range(5)] + [None]

    local_hap, global_hap = [], []
    for anc, ak, sk in zip(ancestries, aba_keys, saige_keys):
        if ak in bmi.index:
            nc = float(bmi.loc[ak,"n_cases"])    if pd.notna(bmi.loc[ak,"n_cases"])    else 0
            nk = float(bmi.loc[ak,"n_controls"]) if pd.notna(bmi.loc[ak,"n_controls"]) else 0
            global_hap.append((nc + nk) * 2.0)
        else:
            global_hap.append(np.nan)
        local_hap.append(saige_med[sk] if (sk and saige_med.get(sk)) else np.nan)

    def fmt_n(n):
        if pd.isna(n): return "—"
        if n >= 1e6:   return f"{n/1e6:.2f}M"
        if n >= 1e3:   return f"{n/1e3:.1f}k"
        return f"{n:.0f}"

    # ── Use gridspec with an extra row below the axes for ratio labels ──
    fig = plt.figure(figsize=(16, 9))
    gs  = fig.add_gridspec(2, 1, height_ratios=[10, 1.2], hspace=0.08,
                            left=0.12, right=0.97, top=0.90, bottom=0.05)
    ax  = fig.add_subplot(gs[0, 0])
    ax_ratio = fig.add_subplot(gs[1, 0]); ax_ratio.axis("off")

    x = np.arange(len(ancestries))
    bar_w = 0.38
    ax.bar(x - bar_w/2, local_hap,  width=bar_w, color=COLOR_SAIGE,
           edgecolor="white", linewidth=1.4,
           label=f"{NAME_SAIGE} — median N$_{{haplo}}$ per tested variant (local ancestry)")
    ax.bar(x + bar_w/2, global_hap, width=bar_w, color=COLOR_ABA,
           edgecolor="white", linewidth=1.4,
           label=f"{NAME_ABA} — sample size × 2 (global ancestry, haplotype equivalent)")

    # Counts above each bar
    y_max = max([v for v in local_hap + global_hap if pd.notna(v)])
    for i, (loc, glo) in enumerate(zip(local_hap, global_hap)):
        if pd.notna(loc):
            ax.text(x[i]-bar_w/2, loc + y_max*0.015, fmt_n(loc),
                    ha="center", va="bottom", fontsize=14,
                    color="#7d1742", weight="bold")
        if pd.notna(glo):
            ax.text(x[i]+bar_w/2, glo + y_max*0.015, fmt_n(glo),
                    ha="center", va="bottom", fontsize=14,
                    color="#1a3e8a", weight="bold")

    ax.set_xticks(x); ax.set_xticklabels(ancestries, fontsize=15, weight="bold")
    ax.set_ylabel("Number of haplotypes (BMI)")
    ax.set_title(f"{NAME_SAIGE} local-ancestry haplotype count vs {NAME_ABA} global-ancestry haplotype count\n(All of Us, body mass index)",
                 loc="left", weight="bold", pad=16)
    ax.legend(loc="upper right", framealpha=0.96, edgecolor="#cfd8dc")
    ax.grid(axis="y", color="#e8e8e8", linewidth=0.9); ax.set_axisbelow(True)
    ax.set_ylim(0, y_max * 1.20)
    ax.set_xlim(-0.55, len(ancestries) - 0.45)

    # ── Ratio row (lives in its own axes below) ────────────────────────────
    ax_ratio.set_xlim(ax.get_xlim())
    ax_ratio.set_ylim(0, 1)
    # Row label is placed FAR to the left, outside the bar area, with right-alignment
    ax_ratio.text(-0.52, 0.55, "Local / global\nratio:",
                  ha="right", va="center", fontsize=13, weight="bold", color="#333")
    for i, (loc, glo) in enumerate(zip(local_hap, global_hap)):
        if pd.notna(loc) and pd.notna(glo) and glo > 0:
            r = loc / glo
            col = "#1b5e20" if r >= 1.10 else ("#0d47a1" if r >= 0.9
                   else ("#e65100" if r >= 0.5 else "#b71c1c"))
            ax_ratio.text(x[i], 0.55, f"{r:.2f}×",
                          ha="center", va="center", fontsize=15,
                          weight="bold", color=col)
        elif pd.isna(loc) and pd.notna(glo):
            ax_ratio.text(x[i], 0.55, "no local\n(global only)",
                          ha="center", va="center", fontsize=11,
                          color="#666", style="italic")

    savefig(fig, "fig10_bmi_haplotype_bars")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 11 — Method overlap + near-miss breakdown
# ═══════════════════════════════════════════════════════════════════════════
def _loci_in_method(df, method_is_saige: bool, near_miss_only=False):
    """Return set of (phenotype, locus_id) strings that the method discovered.
       Optionally returns the per-locus 'other' p-value for near-miss tabs."""
    # 'shared' rows count toward both methods; tophit_source rows show who finds it
    if method_is_saige:
        sub = df[(df["tophit_source"]=="SAIGE")]
    else:
        sub = df[(df["tophit_source"]=="ABA")]
    keys = (sub["phenotype"].astype(str) + "::" + sub["locus_id"].astype(str)).tolist()
    return set(keys)

def fig11_method_overlap():
    """CCT-vs-META locus overlap and near-miss breakdown.

    Restricted to one scatter table (CCT-vs-META) so all counts add up
    consistently across panels.  See `fig11b_combiner_supplement` for the
    HET/HOM comparison kept as a supplementary figure.
    """
    df_cct = load_scatter(SCATTER_TABLES["cct"])

    # ── Pull the four locus pools (all from one scatter table) ─────────────
    s_rows = df_cct[df_cct["tophit_source"] == "SAIGE"]
    a_rows = df_cct[df_cct["tophit_source"] == "ABA"]
    n_shared      = int((s_rows["locus_status"] == "shared").sum())
    n_saige_only  = int((s_rows["locus_status"] == "SAIGE_only").sum())
    n_meta_only   = int((a_rows["locus_status"] == "ABA_only").sum())
    n_saige_total = n_shared + n_saige_only
    n_meta_total  = n_shared + n_meta_only
    n_total_loci  = n_shared + n_saige_only + n_meta_only

    # ── Near-miss breakdown (for SAIGE-only and META-only loci) ───────────
    def breakdown(sub, other_p_col):
        other = pd.to_numeric(sub[other_p_col], errors="coerce")
        n  = len(sub)
        gw = int((other < 5e-8).sum())
        sg = int(((other >= 5e-8) & (other < 5e-6)).sum())
        nm = int(((other >= 5e-6) & (other < 0.01)).sum())
        ge = int((other.isna() | (other >= 0.01)).sum())
        return n, gw, sg, nm, ge

    n_so, so_gw, so_sg, so_nm, so_ge = breakdown(
        df_cct[(df_cct["tophit_source"]=="SAIGE") & (df_cct["locus_status"]=="SAIGE_only")],
        "ABA_pvalue")
    n_mo, mo_gw, mo_sg, mo_nm, mo_ge = breakdown(
        df_cct[(df_cct["tophit_source"]=="ABA") & (df_cct["locus_status"]=="ABA_only")],
        "SAIGE_pvalue")
    rows_c   = np.array([[so_gw, so_sg, so_nm, so_ge],
                          [mo_gw, mo_sg, mo_nm, mo_ge]], dtype=float)
    labels_c = [f"SAIGE-Tractor (CCT) unique\n(n = {n_so})",
                f"{NAME_ABA} meta-analysis unique\n(n = {n_mo})"]
    pct = rows_c / rows_c.sum(axis=1, keepdims=True).clip(min=1) * 100

    cats = [
        "Genome-wide significant in the other method (p < 5×10⁻⁸)",
        "Near-miss in the other method (5×10⁻⁸ ≤ p < 5×10⁻⁶)",
        "Weak signal in the other method (5×10⁻⁶ ≤ p < 0.01)",
        "Genuine miss (p ≥ 0.01 or variant not tested)",
    ]
    cat_colors = ["#2e7d32", "#fb8c00", "#fdd835", "#b71c1c"]

    # ── Layout ─────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(15, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.5], wspace=1.2)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # ── Panel a: TWO VERTICAL stacked bars — totals on top, no bleed ──────
    x = np.array([0, 1])
    bar_w = 0.55
    # Bar 0 = SAIGE-Tractor CCT (shared base + CCT-only top)
    ax_a.bar(x[0], n_shared,     width=bar_w, color="#7fb069",
             edgecolor="white", label="Shared (found by both methods)")
    ax_a.bar(x[0], n_saige_only, width=bar_w, bottom=n_shared,
             color="#d8527c", edgecolor="white",
             label="SAIGE-Tractor (CCT) only")
    # Bar 1 = All by All META  (shared base + META-only top)
    ax_a.bar(x[1], n_shared,    width=bar_w, color="#7fb069", edgecolor="white")
    ax_a.bar(x[1], n_meta_only, width=bar_w, bottom=n_shared,
             color="#5b8def", edgecolor="white",
             label=f"{NAME_ABA} meta-analysis only")

    # Inside-bar segment labels (centred)
    def _annot(xi, y_bottom, h, text):
        if h >= 25:
            ax_a.text(xi, y_bottom + h/2, text, ha="center", va="center",
                      color="white", fontsize=10.5, weight="bold")
    _annot(x[0], 0, n_shared, f"{n_shared}")
    _annot(x[0], n_shared, n_saige_only, f"{n_saige_only}")
    _annot(x[1], 0, n_shared, f"{n_shared}")
    _annot(x[1], n_shared, n_meta_only, f"{n_meta_only}")
    # Totals above each bar
    ax_a.text(x[0], n_saige_total + 12, f"{n_saige_total}\nloci",
              ha="center", va="bottom", fontsize=11,
              color="#7d1742", weight="bold")
    ax_a.text(x[1], n_meta_total + 12,  f"{n_meta_total}\nloci",
              ha="center", va="bottom", fontsize=11,
              color="#1a3e8a", weight="bold")

    ax_a.set_xticks(x)
    ax_a.set_xticklabels([f"{NAME_SAIGE}\n(CCT)", f"{NAME_ABA}\nmeta-analysis"], fontsize=11)
    ax_a.set_ylabel("Number of independent loci")
    ax_a.set_title("B. Locus discovery — SAIGE-Tractor CCT vs All-by-All META",
                   loc="left", weight="bold", fontsize=12, pad=14)
    overlap_pct_of_cct  = 100 * n_shared / n_saige_total if n_saige_total else 0
    overlap_pct_of_meta = 100 * n_shared / n_meta_total  if n_meta_total  else 0
    ax_a.text(0.02, 0.985,
              f"Total unique loci (union) = {n_total_loci}\n"
              f"Shared = {n_shared}  ({overlap_pct_of_cct:.0f}% of CCT, "
              f"{overlap_pct_of_meta:.0f}% of META)",
              transform=ax_a.transAxes, fontsize=9.5, color="#333",
              va="top", ha="left",
              bbox=dict(facecolor="#f5f7fa", edgecolor="#cfd8dc", pad=4))
    ax_a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16),
                ncol=1, fontsize=10, frameon=False,
                title="Locus category", title_fontsize=10)
    for spine in ("top","right"): ax_a.spines[spine].set_visible(False)
    ax_a.grid(axis="y", color="#eee", linewidth=0.8); ax_a.set_axisbelow(True)
    ax_a.set_ylim(0, max(n_saige_total, n_meta_total) * 1.18)
    ax_a.set_xlim(-0.7, 1.7)

    # ── Panel b: 2-row near-miss breakdown ────────────────────────────────
    y = np.arange(len(labels_c))
    bottoms = np.zeros(len(labels_c))
    for j, cat in enumerate(cats):
        ax_b.barh(y, pct[:, j], left=bottoms, color=cat_colors[j], label=cat,
                  height=0.55, edgecolor="white", linewidth=0.6)
        for yi, val in enumerate(pct[:, j]):
            if val >= 5:
                ax_b.text(bottoms[yi] + val/2, yi,
                          f"{int(rows_c[yi, j])}  ({val:.0f}%)",
                          ha="center", va="center", color="white",
                          fontsize=10, weight="bold")
        bottoms += pct[:, j]
    ax_b.set_yticks(y); ax_b.set_yticklabels(labels_c, fontsize=11)
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 100)
    ax_b.set_xlabel("Percent of unique loci")
    ax_b.set_title("C. Where each method's unique loci stand in the other method",
                   loc="left", weight="bold", fontsize=12, pad=12)
    ax_b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14),
                ncol=1, fontsize=10, frameon=False,
                title="Status in the other method (replication tier)",
                title_fontsize=10)
    for spine in ("top","right"): ax_b.spines[spine].set_visible(False)

    fig.suptitle(f"Locus discovery and replication tiers — {NAME_SAIGE} CCT vs {NAME_ABA_META}",
                 x=0.06, ha="left", y=1.02, weight="bold")
    savefig(fig, "fig11_method_overlap")

    # ─── JSON dump (clean CCT-vs-META only) ────────────────────────────────
    out = {
        "comparison": "CCT-vs-META (single scatter table; counts add up across panels)",
        "panel_a": {
            "shared": n_shared,
            "saige_tractor_cct_only": n_saige_only,
            "aba_meta_only": n_meta_only,
            "saige_tractor_cct_total": n_saige_total,
            "aba_meta_total": n_meta_total,
            "total_unique_loci": n_total_loci,
        },
        "panel_b": {
            "saige_only_n_loci": n_so,
            "saige_only_breakdown_in_meta": dict(
                genome_wide=int(so_gw), near_miss=int(so_sg),
                weak=int(so_nm), genuine=int(so_ge)),
            "meta_only_n_loci": n_mo,
            "meta_only_breakdown_in_cct": dict(
                genome_wide=int(mo_gw), near_miss=int(mo_sg),
                weak=int(mo_nm), genuine=int(mo_ge)),
        },
    }
    out_path = REPLICATION_DIR / "method_overlap_counts.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[05] wrote {out_path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 11b (SUPPLEMENT) — HET / HOM combiner panel (kept separate so the
# main-text Fig 11 has no inconsistent totals)
# ═══════════════════════════════════════════════════════════════════════════
def fig11_supp_combiners():
    """Supplementary: vertical stacked-bar comparison for each of CCT, HOM, HET
    against ABA META.  All within-panel counts add up. META totals differ
    slightly across panels because each scatter table was built by independently
    clumping META top hits against the SAIGE-side combiner seed (noted in the
    caption)."""
    df_cct = load_scatter(SCATTER_TABLES["cct"])
    df_hom = load_scatter(SCATTER_TABLES["hom"])
    df_het = load_scatter(SCATTER_TABLES["het"])

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.0), sharey=True)
    y_max_global = 0
    for ax, df, name in [
        (axes[0], df_cct, "CCT"),
        (axes[1], df_hom, "HOM"),
        (axes[2], df_het, "HET"),
    ]:
        s = df[df["tophit_source"]=="SAIGE"]
        a = df[df["tophit_source"]=="ABA"]
        sh = int((s["locus_status"]=="shared").sum())
        so = int((s["locus_status"]=="SAIGE_only").sum())
        mo = int((a["locus_status"]=="ABA_only").sum())
        st = sh + so; mt = sh + mo
        y_max_global = max(y_max_global, st, mt)
        x = np.array([0, 1]); bar_w = 0.55
        ax.bar(x[0], sh, width=bar_w, color="#7fb069", edgecolor="white",
               label="Shared")
        ax.bar(x[0], so, width=bar_w, bottom=sh, color="#d8527c",
               edgecolor="white", label=f"SAIGE-Tractor ({name}) only")
        ax.bar(x[1], sh, width=bar_w, color="#7fb069", edgecolor="white")
        ax.bar(x[1], mo, width=bar_w, bottom=sh, color="#5b8def",
               edgecolor="white", label=f"{NAME_ABA} meta-analysis only")
        for xi, segs in zip(x, [(sh, so), (sh, mo)]):
            b = 0
            for seg in segs:
                if seg >= 25:
                    ax.text(xi, b + seg/2, str(seg), ha="center", va="center",
                            color="white", fontsize=10, weight="bold")
                b += seg
        ax.text(x[0], st + 10, f"{st}", ha="center", va="bottom",
                fontsize=10, color="#7d1742", weight="bold")
        ax.text(x[1], mt + 10, f"{mt}", ha="center", va="bottom",
                fontsize=10, color="#1a3e8a", weight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels([f"SAIGE-T\n({name})", f"{NAME_ABA}\nmeta-analysis"], fontsize=10)
        ax.set_title(f"{name}-vs-META", loc="left", weight="bold", fontsize=12)
        ax.set_xlim(-0.7, 1.7)
        if ax is axes[0]:
            ax.set_ylabel("Number of independent loci")
            ax.legend(loc="upper left", fontsize=8.5, frameon=True, framealpha=0.95)
        for spine in ("top","right"): ax.spines[spine].set_visible(False)
        ax.grid(axis="y", color="#eee", linewidth=0.8); ax.set_axisbelow(True)
    for ax in axes:
        ax.set_ylim(0, y_max_global * 1.18)

    fig.suptitle(f"Locus overlap with {NAME_ABA_META} for all three "
                 f"{NAME_SAIGE} combiners.  Within each panel, counts add up; "
                 f"meta-analysis totals differ slightly between panels because each "
                 f"scatter table independently clumps {NAME_ABA} hits against the "
                 f"{NAME_SAIGE} seed.",
                 fontsize=12, x=0.06, ha="left", y=1.04, wrap=True)
    savefig(fig, "supp_fig_hom_het_overlap")

# ───────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    fig9_meta_anchored_scatter()
    fig10_bmi_haplotype_bars()
    fig11_method_overlap()
    fig11_supp_combiners()
    print("[05] done.")
