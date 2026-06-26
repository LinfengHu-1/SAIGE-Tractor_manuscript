#!/usr/bin/env python3
"""
12_panel_b_workflow.py
──────────────────────
Standalone Panel B — SAIGE-Tractor method workflow with the four method
strengths integrated as in-line visualisations at each pipeline stage.

  1. INPUT       — Inclusive cohort  (solid + striped/admixed persons)
  2. STORAGE     — Scalable & compressed  (raw stack → tiny stack + clock)
  3. ANALYSIS    — HOM + HET tests, real math  (5 score vectors flowing into
                    two complementary aggregations: sum-then-square vs
                    square-then-sum; then CCT combination)
  4. OUTPUT      — Method products  (per-ancestry β/SE/p table, combined
                    p_CCT, and the cross-ancestry r_g / h² matrix)

Design constraints from user feedback
- No background fills on stage cards — only thin borders + title labels —
  so each piece can be lifted/edited independently.
- Less prose, more "show": equations and visuals carry the load.
- Some persons in Stage 1 must look ADMIXED (multiple colours in one body),
  not just different solid colours.
- Stage 3 visualisation must reflect the actual mathematical difference
  between HOM (1-df) and HET (K-df) score tests.
- Stage 4 must demonstrate the workflow's downstream products (not just
  one strength).

Output:  manuscript_figures/panel_b_workflow.{pdf,png}
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (FancyBboxPatch, Rectangle, Circle,
                                  Ellipse, Polygon)
from common import FIGS_DIR, apply_manuscript_style

apply_manuscript_style()
plt.rcParams.update({"font.family": "DejaVu Sans"})

CB_AFR    = "#D55E00"
CB_EAS    = "#56B4E9"
CB_EUR    = "#0072B2"
CB_NATAM  = "#E69F00"
CB_SAS    = "#CC79A7"
CB_GRAY   = "#999999"

COL_HOM   = "#0072B2"
COL_HET   = "#E69F00"
COL_CCT   = "#D55E00"

STG_EC = ["#0d1b4d", "#7a3e00", "#7d1742", "#1f6f4a"]


def person_solid(ax, cx, cy, *, body_col=CB_GRAY, head_col="#cfd8dc",
                  scale=1.0):
    head_r = 0.22 * scale
    body_w = 0.50 * scale
    body_h = 0.60 * scale
    ax.add_patch(Polygon([(cx - body_w/2, cy),
                           (cx + body_w/2, cy),
                           (cx + body_w/2 * 0.72, cy + body_h),
                           (cx - body_w/2 * 0.72, cy + body_h)],
                          closed=True, facecolor=body_col,
                          edgecolor="#37474f", lw=0.5, zorder=4))
    ax.add_patch(Circle((cx, cy + body_h + head_r * 0.55), head_r,
                          facecolor=head_col, edgecolor="#37474f",
                          lw=0.5, zorder=4))


def person_admixed(ax, cx, cy, *, body_cols, head_col="#cfd8dc", scale=1.0):
    """Trapezoidal body split vertically into stripes (one per ancestry
    contributing to this person), to signify local-ancestry admixture."""
    head_r = 0.22 * scale
    body_w = 0.50 * scale
    body_h = 0.60 * scale
    # Build the trapezoid corners
    bx0, bx1 = cx - body_w/2,            cx + body_w/2
    tx0, tx1 = cx - body_w/2 * 0.72,     cx + body_w/2 * 0.72
    n = len(body_cols)
    # Interpolate body left and right edges as linear functions of y
    ys = np.linspace(0, 1, 32)
    left_at  = lambda t: bx0 + t * (tx0 - bx0)
    right_at = lambda t: bx1 + t * (tx1 - bx1)
    for i, col in enumerate(body_cols):
        frac0 = i / n
        frac1 = (i + 1) / n
        # Build a parallelogram-ish strip between vertical x slices
        # interpolated through the trapezoid
        verts = []
        for t in ys:
            x_lo = left_at(t)  * (1 - frac0) + right_at(t) * frac0
            verts.append((x_lo, cy + t * body_h))
        for t in ys[::-1]:
            x_hi = left_at(t)  * (1 - frac1) + right_at(t) * frac1
            verts.append((x_hi, cy + t * body_h))
        ax.add_patch(Polygon(verts, closed=True, facecolor=col,
                              edgecolor="none", lw=0, zorder=4))
    # Outline the body
    ax.add_patch(Polygon([(bx0, cy), (bx1, cy), (tx1, cy + body_h),
                           (tx0, cy + body_h)], closed=True,
                          facecolor="none", edgecolor="#37474f", lw=0.5,
                          zorder=5))
    ax.add_patch(Circle((cx, cy + body_h + head_r * 0.55), head_r,
                          facecolor=head_col, edgecolor="#37474f",
                          lw=0.5, zorder=4))


def stage_outline(ax, x, y, w, h, title, ec):
    """Thin-border stage card with title above — NO background fill."""
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                  boxstyle="round,pad=0.04,rounding_size=0.15",
                                  fc="none", ec=ec, lw=1.4, zorder=1))
    ax.text(x + w/2, y + h + 0.20, title,
            ha="center", va="bottom", fontsize=13.5, weight="bold",
            color=ec, zorder=2)


# ── STAGE 1 ─────────────────────────────────────────────────────────────────
def stage_input(ax, x0, y0, w, h):
    stage_outline(ax, x0, y0, w, h, "1. Inclusive cohort", STG_EC[0])

    # Crowd: ~50% solid (single-ancestry), ~30% 2-colour (admixed),
    # ~20% 3-colour (heavily admixed) — explicitly demonstrates the
    # admixed individuals SAIGE-Tractor recovers that conventional global-
    # ancestry pipelines exclude.
    layout = [
        (CB_AFR,),                       (CB_EUR,),         (CB_AFR, CB_EUR),    (CB_NATAM,),       (CB_EUR, CB_NATAM, CB_AFR),
        (CB_EUR,),                       (CB_EUR, CB_AFR),  (CB_AFR, CB_NATAM),  (CB_EUR, CB_SAS),  (CB_EAS,),
        (CB_NATAM, CB_EUR),              (CB_AFR,),         (CB_EUR,),           (CB_AFR, CB_EUR),  (CB_EUR, CB_NATAM),
        (CB_SAS,),                       (CB_EUR, CB_AFR),  (CB_NATAM, CB_EUR),  (CB_AFR,),         (CB_EAS, CB_EUR, CB_SAS),
    ]
    n_rows, n_cols = 4, 5
    body_y_bot = y0 + 1.00
    body_y_top = y0 + h - 1.20
    body_x_lf  = x0 + 0.40
    body_x_rt  = x0 + w * 0.66
    dxs = np.linspace(body_x_lf, body_x_rt, n_cols)
    dys = np.linspace(body_y_bot, body_y_top, n_rows)
    for r, py in enumerate(dys):
        for c, px in enumerate(dxs):
            jitter = (r % 2) * (dxs[1] - dxs[0]) * 0.18
            cols = layout[r * n_cols + c]
            if len(cols) == 1:
                person_solid(ax, px + jitter, py,
                              body_col=cols[0], scale=0.85)
            else:
                person_admixed(ax, px + jitter, py,
                                body_cols=list(cols), scale=0.85)
    # Funnel arrow
    ax.annotate("",
                xy=(x0 + w - 0.20, y0 + h/2 - 0.30),
                xytext=(x0 + w * 0.80, y0 + h/2 - 0.30),
                arrowprops=dict(arrowstyle="-|>", color=STG_EC[0],
                                 lw=2.6, mutation_scale=24),
                zorder=5)
    # Small in-stage caption + key swatches
    ax.text(x0 + 0.40, y0 + 0.55,
            "■ single ancestry      ▌▌▌ admixed (recovered)",
            ha="left", va="bottom", fontsize=9.5, color="#444")


# ── STAGE 2 ─────────────────────────────────────────────────────────────────
def stage_storage(ax, x0, y0, w, h):
    stage_outline(ax, x0, y0, w, h,
                   "2. Scalable & compressed storage", STG_EC[1])
    # Raw stack
    raw_x = x0 + 0.55
    raw_y = y0 + 2.30
    raw_w = 0.90
    raw_h = 1.85
    n_disks = 6
    disk_h = raw_h / n_disks
    for k in range(n_disks):
        ax.add_patch(Ellipse((raw_x + raw_w/2, raw_y + k * disk_h),
                              raw_w, 0.24, facecolor=CB_GRAY,
                              edgecolor="white", lw=1.0, zorder=4 + k))
        if k < n_disks - 1:
            ax.add_patch(Rectangle((raw_x, raw_y + k * disk_h),
                                    raw_w, disk_h,
                                    facecolor=CB_GRAY, edgecolor="white",
                                    lw=0.8, zorder=3 + k))
    ax.text(raw_x + raw_w/2, raw_y + raw_h + 0.35, "raw\ngenotypes",
            ha="center", va="bottom", fontsize=10, color="#222",
            linespacing=1.0)
    # Arrow + compress label
    arr_y = raw_y + raw_h/2
    arr_x_left  = raw_x + raw_w + 0.20
    arr_x_right = arr_x_left + 1.05
    ax.annotate("",
                xy=(arr_x_right, arr_y),
                xytext=(arr_x_left, arr_y),
                arrowprops=dict(arrowstyle="-|>", color=STG_EC[1],
                                 lw=2.4, mutation_scale=20),
                zorder=8)
    ax.text((arr_x_left + arr_x_right)/2, arr_y + 0.28, "compress",
            ha="center", va="bottom", fontsize=10.5, color=STG_EC[1],
            style="italic")
    # Tiny stack
    tiny_x = arr_x_right + 0.15
    tiny_y = raw_y + 0.90
    tiny_w = 0.45
    tiny_h = 0.45
    n_tiny = 3
    disk_h_t = tiny_h / n_tiny
    for k in range(n_tiny):
        ax.add_patch(Ellipse((tiny_x + tiny_w/2, tiny_y + k * disk_h_t),
                              tiny_w, 0.12, facecolor=CB_EAS,
                              edgecolor="white", lw=0.8, zorder=5 + k))
        if k < n_tiny - 1:
            ax.add_patch(Rectangle((tiny_x, tiny_y + k * disk_h_t),
                                    tiny_w, disk_h_t,
                                    facecolor=CB_EAS, edgecolor="white",
                                    lw=0.7, zorder=4 + k))
    ax.text(tiny_x + tiny_w/2, tiny_y + tiny_h + 0.35,
            r"DS$_k$, ANC$_k$",
            ha="center", va="bottom", fontsize=10.5, color="#222",
            weight="bold")
    # Runtime tag
    tag_w, tag_h = 3.30, 0.60
    tag_x = x0 + w/2 - tag_w/2
    tag_y = y0 + 1.30
    ax.add_patch(FancyBboxPatch((tag_x, tag_y), tag_w, tag_h,
                                  boxstyle="round,pad=0.02,rounding_size=0.10",
                                  fc="white", ec=STG_EC[1], lw=1.4,
                                  zorder=20))
    sw_cx = tag_x + 0.32
    sw_cy = tag_y + tag_h/2
    ax.add_patch(Circle((sw_cx, sw_cy), 0.20,
                          facecolor="white", edgecolor=STG_EC[1],
                          lw=1.6, zorder=21))
    ax.plot([sw_cx, sw_cx + 0.14],
             [sw_cy, sw_cy + 0.10],
             color=STG_EC[1], lw=1.8, zorder=22)
    ax.add_patch(Circle((sw_cx, sw_cy), 0.03,
                          facecolor=STG_EC[1], zorder=23))
    ax.text(sw_cx + 0.36, sw_cy,
            "reruns in hours, not days",
            ha="left", va="center", fontsize=11.5, color=STG_EC[1],
            weight="bold", zorder=24)


# ── STAGE 3 — actual statistical difference between HOM and HET ─────────────
def stage_testing(ax, x0, y0, w, h):
    stage_outline(ax, x0, y0, w, h,
                   "3. Two complementary score tests", STG_EC[2])
    # Row of 5 per-ancestry score vectors at the top of the stage
    ancs = [("AFR",   CB_AFR),
            ("EAS",   CB_EAS),
            ("EUR",   CB_EUR),
            ("NatAm", CB_NATAM),
            ("SAS",   CB_SAS)]
    n = len(ancs)
    chip_y = y0 + h - 1.50
    chip_w = 0.62; chip_h = 0.48
    chip_gap = 0.10
    total_w = n * chip_w + (n - 1) * chip_gap
    chips_x0 = x0 + w/2 - total_w/2
    chip_centres_x = []
    for i, (name, col) in enumerate(ancs):
        cx = chips_x0 + i * (chip_w + chip_gap)
        chip_centres_x.append(cx + chip_w/2)
        ax.add_patch(FancyBboxPatch((cx, chip_y), chip_w, chip_h,
                                      boxstyle="round,pad=0.01,rounding_size=0.05",
                                      fc=col, ec="white", lw=1.0,
                                      zorder=5))
        ax.text(cx + chip_w/2, chip_y + chip_h/2,
                f"U$_{{{name}}}$",
                ha="center", va="center", fontsize=10.5,
                color="white", weight="bold", zorder=6)
    # HOM vs HET test boxes, side by side
    box_y      = y0 + 1.60
    box_h      = 1.65
    box_w      = (w - 0.80) / 2 - 0.15
    box_x_hom  = x0 + 0.40
    box_x_het  = x0 + w - 0.40 - box_w
    # HOM box
    ax.add_patch(FancyBboxPatch((box_x_hom, box_y), box_w, box_h,
                                  boxstyle="round,pad=0.02,rounding_size=0.08",
                                  fc="none", ec=COL_HOM, lw=1.6,
                                  zorder=5))
    ax.text(box_x_hom + box_w/2, box_y + box_h - 0.18,
            "HOM",
            ha="center", va="top", fontsize=13, weight="bold",
            color=COL_HOM, zorder=6)
    ax.text(box_x_hom + box_w/2, box_y + box_h/2 + 0.05,
            r"$T_{\mathrm{hom}}\,=\,\frac{(\,\sum_k U_k\,)^{2}}{\sum_k V_k}$",
            ha="center", va="center", fontsize=15, color="#111", zorder=6)
    ax.text(box_x_hom + box_w/2, box_y + 0.22,
            r"$\sim\,\chi^{2}_{1}$  ·  shared effect",
            ha="center", va="center", fontsize=10.5, color=COL_HOM,
            style="italic", zorder=6)
    # HET box
    ax.add_patch(FancyBboxPatch((box_x_het, box_y), box_w, box_h,
                                  boxstyle="round,pad=0.02,rounding_size=0.08",
                                  fc="none", ec=COL_HET, lw=1.6,
                                  zorder=5))
    ax.text(box_x_het + box_w/2, box_y + box_h - 0.18,
            "HET",
            ha="center", va="top", fontsize=13, weight="bold",
            color=COL_HET, zorder=6)
    ax.text(box_x_het + box_w/2, box_y + box_h/2 + 0.05,
            r"$T_{\mathrm{het}}\,=\,\sum_k \frac{U_k^{\,2}}{V_k}$",
            ha="center", va="center", fontsize=15, color="#111", zorder=6)
    ax.text(box_x_het + box_w/2, box_y + 0.22,
            r"$\sim\,\chi^{2}_{K}$  ·  effects may differ",
            ha="center", va="center", fontsize=10.5, color=COL_HET,
            style="italic", zorder=6)
    # Arrows from each ancestry chip into both boxes (curved)
    for cx in chip_centres_x:
        # to HOM (left box)
        ax.annotate("",
                      xy=(box_x_hom + box_w * 0.85, box_y + box_h - 0.03),
                      xytext=(cx, chip_y),
                      arrowprops=dict(arrowstyle="-", color=COL_HOM,
                                       lw=0.6, alpha=0.55,
                                       connectionstyle="arc3,rad=-0.10"),
                      zorder=3)
        # to HET (right box)
        ax.annotate("",
                      xy=(box_x_het + box_w * 0.15, box_y + box_h - 0.03),
                      xytext=(cx, chip_y),
                      arrowprops=dict(arrowstyle="-", color=COL_HET,
                                       lw=0.6, alpha=0.55,
                                       connectionstyle="arc3,rad=0.10"),
                      zorder=3)
    # CCT combination strip at the very bottom
    cct_w = 4.00
    cct_h = 0.65
    cct_x = x0 + w/2 - cct_w/2
    cct_y = y0 + 0.50
    ax.add_patch(FancyBboxPatch((cct_x, cct_y), cct_w, cct_h,
                                  boxstyle="round,pad=0.02,rounding_size=0.10",
                                  fc="none", ec=COL_CCT, lw=1.8, zorder=10))
    ax.text(cct_x + cct_w/2, cct_y + cct_h/2,
            r"$p_{\mathrm{CCT}}\;=\;\mathrm{CCT}(\,p_{\mathrm{hom}},\,p_{\mathrm{het}}\,)$",
            ha="center", va="center", fontsize=13, color=COL_CCT,
            weight="bold", zorder=11)
    # Two arrows from each test box into CCT
    ax.annotate("", xy=(cct_x + cct_w * 0.30, cct_y + cct_h + 0.04),
                  xytext=(box_x_hom + box_w/2, box_y),
                  arrowprops=dict(arrowstyle="-|>", color=COL_HOM,
                                   lw=1.6, mutation_scale=14,
                                   connectionstyle="arc3,rad=-0.18"),
                  zorder=6)
    ax.annotate("", xy=(cct_x + cct_w * 0.70, cct_y + cct_h + 0.04),
                  xytext=(box_x_het + box_w/2, box_y),
                  arrowprops=dict(arrowstyle="-|>", color=COL_HET,
                                   lw=1.6, mutation_scale=14,
                                   connectionstyle="arc3,rad=0.18"),
                  zorder=6)


# ── STAGE 4 — three concrete outputs ────────────────────────────────────────
def stage_output(ax, x0, y0, w, h):
    stage_outline(ax, x0, y0, w, h,
                   "4. Method outputs", STG_EC[3])

    # ─── Output A: per-ancestry β/SE/p table (top half) ────────────────────
    tbl_x  = x0 + 0.30
    tbl_y  = y0 + h - 3.30
    tbl_w  = w - 0.60
    tbl_h  = 2.30
    col_w  = tbl_w / 4
    row_h  = tbl_h / 6   # 1 header + 5 ancestry rows
    # Header
    headers = ["ancestry", r"$\widehat{\beta}$", "SE", "p"]
    for j, hd in enumerate(headers):
        hx = tbl_x + j * col_w
        ax.add_patch(Rectangle((hx, tbl_y + 5 * row_h), col_w, row_h,
                                facecolor="#37474f", edgecolor="white",
                                lw=0.6, zorder=4))
        ax.text(hx + col_w/2, tbl_y + 5.5 * row_h, hd,
                ha="center", va="center", fontsize=10.5,
                color="white", weight="bold", zorder=5)
    # 5 ancestry rows with placeholder symbolic content
    rows = [
        ("AFR",   CB_AFR,   r"$\widehat{\beta}_{\mathrm{AFR}}$",   "SE", "p"),
        ("EAS",   CB_EAS,   r"$\widehat{\beta}_{\mathrm{EAS}}$",   "SE", "p"),
        ("EUR",   CB_EUR,   r"$\widehat{\beta}_{\mathrm{EUR}}$",   "SE", "p"),
        ("NatAm", CB_NATAM, r"$\widehat{\beta}_{\mathrm{NatAm}}$", "SE", "p"),
        ("SAS",   CB_SAS,   r"$\widehat{\beta}_{\mathrm{SAS}}$",   "SE", "p"),
    ]
    for i, (name, col, b_sym, se_sym, p_sym) in enumerate(rows):
        ry = tbl_y + (4 - i) * row_h
        # zebra striping
        if i % 2 == 0:
            ax.add_patch(Rectangle((tbl_x, ry), tbl_w, row_h,
                                    facecolor="#f5f7fa", edgecolor="none",
                                    zorder=2))
        # Ancestry column with colored swatch + name
        ax.add_patch(Rectangle((tbl_x + 0.08, ry + row_h * 0.20),
                                0.18, row_h * 0.60,
                                facecolor=col, edgecolor="white", lw=0.6,
                                zorder=3))
        ax.text(tbl_x + 0.34, ry + row_h/2, name,
                ha="left", va="center", fontsize=10.5, weight="bold",
                color=col, zorder=4)
        # β, SE, p columns — symbolic
        for j, sym in enumerate([b_sym, se_sym, p_sym], start=1):
            hx = tbl_x + j * col_w
            ax.text(hx + col_w/2, ry + row_h/2, sym,
                    ha="center", va="center", fontsize=10.5,
                    color="#222", zorder=4)
        # Row divider
        ax.plot([tbl_x, tbl_x + tbl_w], [ry, ry],
                 color="#cfd8dc", lw=0.4, zorder=3)
    # Table outer border
    ax.add_patch(Rectangle((tbl_x, tbl_y), tbl_w, tbl_h,
                            facecolor="none", edgecolor="#37474f",
                            lw=0.8, zorder=5))
    ax.text(tbl_x, tbl_y + tbl_h + 0.12,
            "per-ancestry estimates at every variant",
            ha="left", va="bottom", fontsize=10, color=STG_EC[3],
            style="italic")

    # ─── Output B: combined p_CCT callout (middle) ─────────────────────────
    cct_y = tbl_y - 0.85
    cct_w = tbl_w
    cct_h = 0.55
    ax.add_patch(FancyBboxPatch((tbl_x, cct_y), cct_w, cct_h,
                                  boxstyle="round,pad=0.02,rounding_size=0.10",
                                  fc="none", ec=COL_CCT, lw=1.4, zorder=4))
    ax.text(tbl_x + cct_w/2, cct_y + cct_h/2,
            r"+ combined  $p_{\mathrm{hom}}\;\;p_{\mathrm{het}}\;\;p_{\mathrm{CCT}}$",
            ha="center", va="center", fontsize=11.5, color=COL_CCT,
            weight="bold", zorder=5)

    # ─── Output C: r_g / h² matrix (bottom) ────────────────────────────────
    m_names = ["AFR", "EUR", "NatAm"]
    m_cols  = [CB_AFR, CB_EUR, CB_NATAM]
    cell    = 0.50
    m_w     = 3 * cell
    m_x     = x0 + w/2 - m_w/2 + 0.15
    m_y     = y0 + 0.55
    for i in range(3):
        for j in range(3):
            cx = m_x + j * cell
            cy = m_y + (2 - i) * cell
            if i == j:
                fc = m_cols[i]; txt = r"$h^{2}$"; tc = "white"
            else:
                fc = "#cfd8dc"; txt = r"$r_g$"; tc = "#222"
            ax.add_patch(Rectangle((cx, cy), cell, cell,
                                    facecolor=fc, edgecolor="white", lw=1.0,
                                    zorder=4))
            ax.text(cx + cell/2, cy + cell/2, txt,
                    ha="center", va="center", fontsize=10.5,
                    color=tc, weight="bold", zorder=5)
    for k, name in enumerate(m_names):
        ax.text(m_x + k * cell + cell/2, m_y + 3 * cell + 0.08,
                name, ha="center", va="bottom", fontsize=8.5,
                color=m_cols[k], weight="bold")
        ax.text(m_x - 0.05, m_y + (2 - k) * cell + cell/2,
                name, ha="right", va="center", fontsize=8.5,
                color=m_cols[k], weight="bold")
    # Single tight label below the matrix
    ax.text(m_x + m_w/2, m_y - 0.15,
            r"within-individual  $r_g$ / $h^{2}$",
            ha="center", va="top", fontsize=10, color=STG_EC[3],
            weight="bold")


# ───────────────────────────────────────────────────────────────────────────
# Build figure
# ───────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(22, 8.5))
ax.set_xlim(0, 22); ax.set_ylim(0, 8.5)
ax.set_aspect("equal")
ax.set_xticks([]); ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

n_stages  = 4
left_pad  = 0.40
right_pad = 0.40
gap       = 0.55
stage_h   = 6.40
stage_y   = 0.70
stage_w   = (22 - left_pad - right_pad - (n_stages - 1) * gap) / n_stages

xs = [left_pad + i * (stage_w + gap) for i in range(n_stages)]

stage_input  (ax, xs[0], stage_y, stage_w, stage_h)
stage_storage(ax, xs[1], stage_y, stage_w, stage_h)
stage_testing(ax, xs[2], stage_y, stage_w, stage_h)
stage_output (ax, xs[3], stage_y, stage_w, stage_h)

# Inter-stage arrows
arrow_y = stage_y + stage_h/2 - 0.10
for i in range(n_stages - 1):
    x_left  = xs[i] + stage_w
    x_right = xs[i+1]
    ax.annotate("",
                  xy=(x_right - 0.03, arrow_y),
                  xytext=(x_left + 0.03, arrow_y),
                  arrowprops=dict(arrowstyle="-|>", color="#37474f",
                                   lw=3.0, mutation_scale=28),
                  zorder=10)

out_pdf = FIGS_DIR / "panel_b_workflow.pdf"
out_png = FIGS_DIR / "panel_b_workflow.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=180)
plt.close(fig)
print(f"[panel B] wrote {out_pdf}")
print(f"[panel B] wrote {out_png}")
