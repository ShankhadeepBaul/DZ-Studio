"""
dz_figs.py  --  the science and the plotting
=============================================

This module draws the two "figure" tabs in DZ Studio: the KDE + pie panels and
the MDS map. 

What lives here, top to bottom:

  1. Configuration    -- two small "settings" objects (KDEConfig, MDSConfig)
                         that hold every knob the figures expose. The GUI fills
                         these in from its widgets; nothing here talks to the
                         GUI directly.
  2. Data loading     -- reading an Excel sheet, listing its samples, and
                         pulling out the age column for each chosen sample.
  3. KDE figure       -- the kernel-density plot with the histogram-scaled
                         y-axis, provenance-coloured fill, peak labels, and the
                         provenance pie chart.
  4. MDS figure       -- turns sample-to-sample dissimilarity into a 2-D map
                         (multidimensional scaling), following Vermeesch (2013).

Design note for readers: this file does no Qt and no file dialogs. It takes
plain inputs (age arrays, a config object) and returns a matplotlib Figure.
That separation is deliberate -- if a NUMBER looks wrong, the bug is in here;
if the WINDOW looks wrong, the bug is in dz_app.py. It also means you can
import these functions into a notebook and use them without the GUI at all.

Run `python dz_figs.py` to render a self-test figure from synthetic data.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.patheffects as pe
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from matplotlib.collections import PolyCollection
from scipy.stats import gaussian_kde, ks_2samp
from scipy.signal import find_peaks
from sklearn.manifold import MDS

import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")


# ======================================================================
#  Configuration
#
#  Two dataclasses that hold every setting the figures expose. Think of
#  them as the old notebook "SETTINGS block" turned into objects: the GUI
#  builds one of these from its widgets and hands it to the figure
#  functions. Every field has a sensible default, so you can also just
#  write `KDEConfig()` in a script and tweak one or two values.
# ======================================================================

DEFAULT_AGE_BINS = [
    (0, 32), (32, 37), (40, 50), (60, 250), (250, 700),
    (900, 1300), (1300, 1600), (1600, 1800), (1800, 2500), (2500, 4000),
]
DEFAULT_BIN_LABELS = [
    "< 32 Ma",
    "F2: Paleogene (32–37 Ma)",
    "F1: Paleogene (40–50 Ma)",
    "E: Cordilleran arc (60–250 Ma)",
    "D: Appalachian–Peri-Gondwana (250–700 Ma)",
    "C3: Grenville (0.9–1.3 Ga)",
    "C2: Midcontinent gr.–rhy. (1.3–1.6 Ga)",
    "C1: Yavapai–Mazatzal (1.6–1.8 Ga)",
    "B: Trans-Hudson–Penokean (1.8–2.5 Ga)",
    "A: Archean (> 2.5 Ga)",
]
DEFAULT_BIN_COLORS = [
    "#D0D0D0", "#F8B6CC", "#E8527A", "#5BA3D9", "#5FB876",
    "#9D6FB8", "#F4C536", "#E08542", "#A87A56", "#5D3A2A",
]


@dataclass
class KDEConfig:
    """All settings for the KDE + pie figure.

    Every field has a default, so `KDEConfig()` gives the standard figure and
    you override only what you want. The GUI builds one of these from its
    widgets. To add a new option: add a field here, use it inside
    make_kde_figure()/_plot_panel(), and add a widget in dz_app.py that sets
    it. Ages and bandwidths are in Ma; panel widths and font sizes are in the
    usual matplotlib units.
    """
    # 3. age ranges
    young_range: tuple = (0, 100)
    old_range: tuple = (100, 4000)
    pie_range: tuple = (0, 4000)
    # 4. bandwidths
    young_bandwidth: float = 3
    old_bandwidth: float = 40
    # 5. histogram bin widths
    young_bin_width: float = 3
    old_bin_width: float = 40
    single_bandwidth: float = 25       # used only when show_split is False
    single_bin_width: float = 25
    # 6. provenance bins
    age_bins: list = field(default_factory=lambda: list(DEFAULT_AGE_BINS))
    bin_labels: list = field(default_factory=lambda: list(DEFAULT_BIN_LABELS))
    bin_colors: list = field(default_factory=lambda: list(DEFAULT_BIN_COLORS))
    # 7. sizes & layout
    fig_width: float = 12
    young_panel_width: float = 2
    old_panel_width: float = 4.5
    pie_panel_width: float = 3
    row_height: float = 2
    vertical_spacing: float = 0.55
    horizontal_spacing: float = 0.10
    # 8. font sizes
    sample_name_size: float = 7
    n_label_size: float = 8
    axis_tick_size: float = 7.5
    axis_label_size: float = 9
    legend_size: float = 7.5
    peak_label_size: float = 6.5
    pie_inside_label_size: float = 8.5
    pie_outside_label_size: float = 7.5
    legend_ncol: int = 5
    # 9. peak labelling
    label_peaks: bool = True
    peak_threshold: float = 0.30
    # 10. pie labels
    pie_inside_threshold: float = 90.0
    pie_min_percent: float = 5.0
    pie_radius: float = 0.9     # keep < ~1.1 so outside labels have room
    pie_label_dist: float = 0.18   # gap beyond the rim for the % labels
    show_pie: bool = True
    show_split: bool = True        # two panels (young/old) vs one full-range
    show_legend: bool = True


@dataclass
class MDSConfig:
    """All settings for the MDS figure.

    Same idea as KDEConfig: defaults give the standard map, the GUI fills it
    from widgets, add a field + use it in make_mds_figure() + add a widget to
    extend it. `metric` picks the dissimilarity ("ks" is the Vermeesch
    default). The normalize/axis_limit/equal_aspect fields only reframe the
    finished map so separate figures share axes -- they do not change the
    statistics.
    """
    age_threshold: float = 70
    metric: str = "ks"                 # 'ks' | 'cross_correlation' | 'likeness'
    kde_bandwidth: float = 40
    mds_type: str = "non-metric"       # 'non-metric' | 'metric'
    seed: int = 42
    n_init: int = 10
    normalize_scale: bool = True       # scale each map so farthest point = radius 1
    axis_limit: float = 1.15           # fixed square axis range when normalising
    equal_aspect: bool = True          # 1 unit on x == 1 unit on y
    fig_width: float = 8.5
    fig_height: float = 8.5
    point_size: float = 140
    point_color: str = "#3B6FB6"
    label_font_size: float = 10
    show_labels: bool = True
    show_neighbor_lines: bool = True
    show_title: bool = True
    show_stress: bool = True
    title: str = ""                    # blank = automatic
    xlabel: str = ""                   # blank = "MDS dimension 1"
    ylabel: str = ""                   # blank = "MDS dimension 2"
    per_sample_colors: bool = False    # colour points individually instead
    colors: dict = field(default_factory=dict)


def apply_rc(cfg: KDEConfig):
    """The notebook's rcParams block."""
    mpl.rcParams.update({
        "font.family":     "serif",
        "axes.labelsize":  cfg.axis_label_size,
        "axes.titlesize":  cfg.axis_label_size,
        "axes.linewidth":  0.7,
        "xtick.labelsize": cfg.axis_tick_size,
        "ytick.labelsize": cfg.axis_tick_size,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "pdf.fonttype":    42,
        "ps.fonttype":     42,
        "svg.fonttype":    "none",
        "savefig.dpi":     300,
        "savefig.bbox":    "tight",
    })


# ======================================================================
#  Data loading
#
#  Small helpers for getting ages out of an Excel workbook: list the
#  sheets, list the sample names in a sheet, and pull the age column for a
#  chosen set of samples (optionally dropping grains younger than a
#  threshold). Everything downstream works on plain NumPy age arrays.
# ======================================================================

def read_sheet(path, sheet_name):
    if str(path).lower().endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_excel(path, sheet_name=sheet_name)


def sheet_names(path):
    if str(path).lower().endswith(".csv"):
        return []
    return pd.ExcelFile(path).sheet_names


def all_samples(df, id_col):
    """Sample names in first-appearance order."""
    seen, out = set(), []
    for s in df[id_col].astype(str):
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def load_ages(df, sample_list, id_col, age_col, threshold=None):
    """
    Notebook's load_ages, but taking an already-read DataFrame.
    threshold: if given, drop grains younger than it (the MDS behaviour).
    """
    ages, labels = [], []
    for s in sample_list:
        a = df.loc[df[id_col].astype(str) == s, age_col]
        a = pd.to_numeric(a, errors="coerce").dropna().to_numpy(float)
        if threshold is not None:
            a = a[a > threshold]
            if a.size < 3:
                continue
        if a.size == 0:
            continue
        ages.append(a)
        labels.append(s)
    return ages, labels


# ======================================================================
#  KDE figure
#
#  Builds the kernel-density-estimate panels: one row per sample, with a
#  young-age panel, an old-age panel (or a single combined panel), and an
#  optional provenance pie. The KDE curve is scaled to sit on a histogram
#  count axis, the area under it is coloured by provenance bin, and the
#  main peaks are labelled. This is the plotting logic adapted from
#  Zircon_KDE_main.ipynb.
# ======================================================================

def _kde(data, x_grid, bandwidth):
    """Evaluate a Gaussian KDE of `data` on `x_grid` at a fixed bandwidth in Ma.

    scipy's gaussian_kde expects `bw_method` as a multiple of the data's
    standard deviation, not an absolute width, so we divide the requested
    bandwidth (in Ma) by the sample's std to convert. Returns None if there
    are too few grains to estimate a density.
    """
    if len(data) < 3:
        return None
    s = np.std(data, ddof=1)
    if s == 0:
        return None
    kde = gaussian_kde(data, bw_method=bandwidth / s)
    return kde(x_grid)


def _plot_panel(ax, a, age_range, bandwidth, bin_width, cfg: KDEConfig):
    """Draw one KDE panel (one age window of one sample) onto `ax`.

    The layout choices worth knowing if you want to change this:

    * The KDE curve is rescaled so its peak matches the tallest histogram bar
      (`kde * h_max / kde.max()`). That's what lets a smooth density and a
      discrete count share one y-axis -- the y-axis is really the histogram
      count, and the KDE is drawn to fit it. Remove that scaling if you want a
      true probability-density y-axis instead.
    * The area under the curve is filled in segments, one per provenance
      age-bin (`cfg.age_bins` / `cfg.bin_colors`), which is what colours the
      peaks by source. Gaps between bins are allowed and simply stay unfilled.
    * The histogram outline is drawn as a single PolyCollection rather than
      matplotlib's ax.bar(). It looks identical (open black bars) but is one
      artist instead of one rectangle per bar -- a big speed-up on figures with
      many bins. If you edit the bar look, edit the vertices here.
    """
    a = np.asarray(a, dtype=float)
    a = a[(a >= age_range[0]) & (a <= age_range[1])]

    edges = np.arange(age_range[0], age_range[1] + bin_width, bin_width)
    counts, _ = np.histogram(a, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])

    x_grid = np.linspace(age_range[0], age_range[1], 600)
    kde = _kde(a, x_grid, bandwidth)
    h_max = counts.max() if counts.size and counts.max() > 0 else 1

    # rescale the density so its peak equals the tallest bar (shared y-axis)
    if kde is not None and kde.max() > 0:
        kde_scaled = kde * (h_max / kde.max())
    else:
        kde_scaled = None

    # coloured KDE fill segmented by age bin
    if kde_scaled is not None:
        for (lo, hi), col in zip(cfg.age_bins, cfg.bin_colors):
            lo_p, hi_p = max(lo, age_range[0]), min(hi, age_range[1])
            if hi_p <= lo_p:
                continue
            mask = (x_grid >= lo_p) & (x_grid <= hi_p)
            if mask.sum() < 2:
                continue
            ax.fill_between(x_grid[mask], 0, kde_scaled[mask],
                            color=col, alpha=0.95, linewidth=0)
        ax.plot(x_grid, kde_scaled, color="black", lw=0.7)

    hw = bin_width * 0.95 / 2.0
    verts = [[(c - hw, 0), (c + hw, 0), (c + hw, h), (c - hw, h)]
             for c, h in zip(centers, counts)]
    ax.add_collection(PolyCollection(
        verts, facecolors="none", edgecolors="black",
        linewidths=0.5, zorder=3))

    # peak labels in tidy stacked rows
    max_label_y = h_max
    if cfg.label_peaks and kde is not None and kde_scaled is not None:
        min_gap_ma = bandwidth
        label_pad_ma = bandwidth * 6
        dx = (age_range[1] - age_range[0]) / 599
        min_dist = max(1, int(round(min_gap_ma / dx)))
        pks, _ = find_peaks(kde,
                            height=cfg.peak_threshold * float(kde.max()),
                            distance=min_dist)
        row_right_x = []
        row0_y = h_max * 1.10
        row_step = h_max * 0.16
        for p in pks:
            x_p = x_grid[p]
            row = None
            for ri, right in enumerate(row_right_x):
                if x_p - right >= label_pad_ma:
                    row = ri
                    row_right_x[ri] = x_p
                    break
            if row is None:
                row_right_x.append(x_p)
                row = len(row_right_x) - 1
            label_y = row0_y + row * row_step
            max_label_y = max(max_label_y, label_y + row_step * 0.4)
            ax.plot([x_p, x_p], [kde_scaled[p], label_y - h_max * 0.015],
                    color="0.55", lw=0.4, ls=":", zorder=3.5)
            ax.text(x_p, label_y, f"{x_p:.0f}",
                    fontsize=cfg.peak_label_size, ha="center", va="bottom",
                    color="#222", zorder=4)

    ax.set_xlim(*age_range)
    ax.set_ylim(0, max(h_max * 1.15, max_label_y * 1.05) if h_max > 0 else 1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="both", length=2)


def _plot_pie(ax, a, cfg: KDEConfig):
    a = np.asarray(a, dtype=float)
    a = a[(a > cfg.pie_range[0]) & (a <= cfg.pie_range[1])]
    counts = np.array([((a >= lo) & (a < hi)).sum()
                       for lo, hi in cfg.age_bins], dtype=float)
    total = counts.sum()
    if total == 0:
        ax.axis("off")
        return
    pcts = 100.0 * counts / total

    wedges, _ = ax.pie(
        counts, colors=cfg.bin_colors,
        startangle=90, counterclock=False,
        radius=cfg.pie_radius,
        wedgeprops={"linewidth": 0.4, "edgecolor": "white"},
    )

    R = cfg.pie_radius
    LIM = 1.5
    gap = cfg.pie_label_dist            # how far beyond the rim the label sits
    # label ring: just outside the rim, but never past the viewport.  Clamp the
    # RADIUS (not x and y separately) so every label stays on its slice's
    # radial line and therefore stays angularly aligned at any pie size.
    lr = min(R + gap, LIM - 0.12)
    for w, p in zip(wedges, pcts):
        if p < cfg.pie_min_percent:
            continue
        ang = 0.5 * (w.theta1 + w.theta2)
        cx, cy = np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))
        ha = "left" if cx >= 0 else "right"
        va = "bottom" if cy > 0.35 else "top" if cy < -0.35 else "center"
        ax.text(lr * cx, lr * cy, f"{p:.0f}%",
                ha=ha, va=va, fontsize=cfg.pie_outside_label_size,
                color="#1A1A1A")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Fixed viewport, independent of the radius.  The pie now genuinely fills
    # more of the cell as pie_radius grows (radius 1.0 ~ 70% of the cell,
    # radius 1.4 fills it).  Previously the limits scaled with R too, which
    # cancelled the radius out and made anything past 1 look identical.
    LIM = 1.5
    ax.set_xlim(-LIM, LIM)
    ax.set_ylim(-LIM, LIM)


def make_kde_figure(ages, labels, cfg: KDEConfig, fig: Figure = None) -> Figure:
    """Draw the full KDE + pie figure and return it as a matplotlib Figure.

    `ages` is a list of age arrays (one per sample) and `labels` their names;
    `cfg` is a KDEConfig holding every style choice. Pass an existing `fig` to
    draw into it (the GUI does this to reuse its canvas), or leave it None to
    get a fresh figure. Nothing is saved to disk here -- the caller decides
    what to do with the returned figure.
    """
    apply_rc(cfg)
    n = len(ages)
    if n == 0:
        raise ValueError("no samples selected")

    if fig is None:
        fig = Figure()
    fig.clear()
    fig.set_size_inches(cfg.fig_width, cfg.row_height * n + 1.4)

    # KDE columns: two (young/old) when split, otherwise one full-range panel
    if cfg.show_split:
        kde_cols = [
            dict(rng=cfg.young_range, bw=cfg.young_bandwidth,
                 bin=cfg.young_bin_width, w=cfg.young_panel_width),
            dict(rng=cfg.old_range, bw=cfg.old_bandwidth,
                 bin=cfg.old_bin_width, w=cfg.old_panel_width),
        ]
    else:
        full = (cfg.young_range[0], cfg.old_range[1])
        kde_cols = [
            dict(rng=full, bw=cfg.single_bandwidth, bin=cfg.single_bin_width,
                 w=cfg.young_panel_width + cfg.old_panel_width),
        ]

    nrows = n + (1 if cfg.show_legend else 0)
    heights = ([0.7] if cfg.show_legend else []) + [1.0] * n
    widths = [c["w"] for c in kde_cols]
    if cfg.show_pie:
        widths.append(cfg.pie_panel_width)
    ncols = len(widths)
    gs = fig.add_gridspec(
        nrows=nrows, ncols=ncols,
        width_ratios=widths,
        height_ratios=heights,
        hspace=cfg.vertical_spacing, wspace=cfg.horizontal_spacing,
    )

    off = 0
    if cfg.show_legend:
        ax_leg = fig.add_subplot(gs[0, :])
        ax_leg.axis("off")
        handles = [Patch(facecolor=c, edgecolor="none", label=l)
                   for c, l in zip(cfg.bin_colors, cfg.bin_labels)]
        ax_leg.legend(handles=handles, loc="center", ncol=cfg.legend_ncol,
                      frameon=False, fontsize=cfg.legend_size,
                      handlelength=1.2, handleheight=1.0,
                      columnspacing=1.5, borderpad=0.2)
        off = 1

    for i, (a, lab) in enumerate(zip(ages, labels)):
        a = np.asarray(a, dtype=float)
        kde_axes = []
        for ci, col in enumerate(kde_cols):
            ax = fig.add_subplot(gs[i + off, ci])
            _plot_panel(ax, a, col["rng"], col["bw"], col["bin"], cfg)
            lo, hi = col["rng"]
            n_in = int(((a >= lo) & (a <= hi)).sum())
            ax.text(0.99, 1.05, f"n = {n_in}", transform=ax.transAxes,
                    ha="right", va="bottom",
                    fontsize=cfg.n_label_size, fontweight="bold")
            if i == n - 1:
                ax.set_xlabel("Age (Ma)")
            else:
                ax.set_xticklabels([])
            kde_axes.append(ax)

        # sample name at the top-left of the first KDE panel
        kde_axes[0].text(0.03, 1.05, lab, transform=kde_axes[0].transAxes,
                         ha="left", va="bottom",
                         fontsize=cfg.sample_name_size, fontweight="bold")

        if cfg.show_pie:
            ax_p = fig.add_subplot(gs[i + off, len(kde_cols)])
            _plot_pie(ax_p, a, cfg)

    return fig


# ======================================================================
#  MDS figure   (multidimensional scaling, after Vermeesch 2013)
#
#  Turns a set of samples into a 2-D "map" where similar samples plot
#  close together. The steps are: (1) measure how different every pair of
#  samples is -- pairwise_dissim() -- using the K-S statistic (or KDE
#  cross-correlation / likeness); (2) run_mds() places the samples in 2-D
#  so those dissimilarities are reproduced as distances, seeded from the
#  classical (cmdscale) solution the way R's MASS::isoMDS does; (3)
#  make_mds_figure() draws the points, nearest-neighbour lines, and labels.
#
#  Only the DISTANCES between points carry meaning -- the map can be freely
#  rotated or mirrored and still be the same result, so never read the axes
#  themselves as physical quantities.
# ======================================================================

METRIC_LABEL = {
    "ks":                "K–S",
    "cross_correlation": "Cross-correlation",
    "likeness":          "Likeness",
}


def pairwise_dissim(ages, metric, bw):
    """Build the sample-by-sample dissimilarity matrix that MDS then maps.

    For every pair of samples it measures how different their age
    distributions are, returning an n x n matrix of values from 0 (identical)
    to 1 (completely different). `metric` chooses how the difference is
    measured: "ks" is the Kolmogorov-Smirnov statistic (the Vermeesch default,
    needs no bandwidth); "cross_correlation" and "likeness" compare the two
    KDEs and use the `bw` bandwidth. This matrix -- not the plot -- is the real
    result; the map is just a way to look at it.
    """
    n = len(ages)
    D = np.zeros((n, n))

    if metric == "ks":
        for i in range(n):
            for j in range(i + 1, n):
                D[i, j] = D[j, i] = ks_2samp(ages[i], ages[j]).statistic
        return D

    x = np.linspace(min(a.min() for a in ages),
                    max(a.max() for a in ages), 2000)
    kdes = [gaussian_kde(a, bw_method=bw / np.std(a, ddof=1))(x) for a in ages]

    if metric == "cross_correlation":
        for i in range(n):
            for j in range(i + 1, n):
                D[i, j] = D[j, i] = 1.0 - np.corrcoef(kdes[i], kdes[j])[0, 1]
    elif metric == "likeness":
        kn = [k / k.sum() for k in kdes]
        for i in range(n):
            for j in range(i + 1, n):
                D[i, j] = D[j, i] = 0.5 * np.sum(np.abs(kn[i] - kn[j]))
    else:
        raise ValueError(f"unknown metric: {metric}")
    return D


def classical_mds(D):
    """Classical (metric) MDS embedding from the dissimilarity matrix, via
    eigen-decomposition of the double-centred matrix. This is R's cmdscale(),
    which MASS::isoMDS() uses as its default starting configuration (y =
    cmdscale(d, k)). Deterministic: no random seed."""
    D = np.asarray(D, float)
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    L = np.sqrt(np.clip(vals[:2], 0, None))
    return vecs[:, :2] * L


def run_mds(D, mds_type, seed, n_init=10):
    is_metric = (mds_type == "metric")
    kwargs = dict(n_components=2, dissimilarity="precomputed",
                  metric=is_metric, max_iter=500)
    if not is_metric:
        kwargs["normalized_stress"] = "auto"
        # Follow Vermeesch's provenance package (MASS::isoMDS), which seeds the
        # non-metric fit with the classical (cmdscale) configuration rather than
        # random restarts. Deterministic; the dissimilarities are unchanged, so
        # the structure matches the R output up to rotation/reflection.
        init = classical_mds(D)
        kwargs.update(n_init=1, random_state=0)
        m = MDS(**kwargs)
        coords = m.fit_transform(D, init=init)
    else:
        kwargs.update(n_init=n_init, random_state=seed)
        m = MDS(**kwargs)
        coords = m.fit_transform(D)
    return coords, m.stress_


def make_mds_figure(ages, labels, cfg: MDSConfig, fig: Figure = None):
    """Draw the MDS map and return (figure, dissimilarity matrix, labels, stress).

    Given the per-sample age arrays and an MDSConfig, this computes the
    dissimilarity matrix, runs MDS to place the samples in 2-D, and plots them
    with optional nearest-neighbour lines and labels. It returns the figure
    plus the dissimilarity matrix and the stress value, so the GUI can offer
    the matrix as a CSV export and show the stress in the caption. "Stress" is
    a goodness-of-fit number: near 0 means the 2-D map reproduces the
    dissimilarities well.
    """
    mpl.rcParams.update({
        "font.family": "serif", "font.size": 10,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.dpi": 300, "savefig.bbox": "tight",
    })

    D = pairwise_dissim(ages, cfg.metric, cfg.kde_bandwidth)
    coords, stress = run_mds(D, cfg.mds_type, cfg.seed, cfg.n_init)

    # normalise to a common scale so every figure shares the same axes.
    # MDS scale is arbitrary, so this only rescales/recentres.
    if cfg.normalize_scale:
        coords = coords - coords.mean(axis=0)
        radius = np.sqrt((coords ** 2).sum(axis=1)).max()
        if radius > 0:
            coords = coords / radius

    if fig is None:
        fig = Figure()
    fig.clear()
    fig.set_size_inches(cfg.fig_width, cfg.fig_height)
    ax = fig.add_subplot(111)

    xs, ys = coords[:, 0], coords[:, 1]
    if cfg.normalize_scale:
        ax.set_xlim(-cfg.axis_limit, cfg.axis_limit)
        ax.set_ylim(-cfg.axis_limit, cfg.axis_limit)
    else:
        xpad = 0.25 * (xs.max() - xs.min() + 1e-9)
        ypad = 0.25 * (ys.max() - ys.min() + 1e-9)
        ax.set_xlim(xs.min() - xpad, xs.max() + xpad)
        ax.set_ylim(ys.min() - ypad, ys.max() + ypad)
    if cfg.equal_aspect:
        ax.set_aspect("equal", "box")

    if cfg.show_neighbor_lines:
        seen = set()
        for kind, color, lw, alpha in [("nn1", "0.55", 1.0, 0.65),
                                       ("nn2", "0.80", 0.7, 0.50)]:
            for i in range(len(labels)):
                d = D[i].copy()
                d[i] = np.inf
                j = int(np.argmin(d))
                if kind == "nn2":
                    d[j] = np.inf
                    j = int(np.argmin(d))
                key = (min(i, j), max(i, j))
                if key in seen:
                    continue
                seen.add(key)
                ax.plot([coords[i, 0], coords[j, 0]],
                        [coords[i, 1], coords[j, 1]],
                        color=color, lw=lw, alpha=alpha, zorder=1)

    if cfg.per_sample_colors:
        cols = [cfg.colors.get(l, cfg.point_color) for l in labels]
    else:
        cols = cfg.point_color
    ax.scatter(xs, ys, s=cfg.point_size, c=cols,
               edgecolor="black", linewidth=0.9, zorder=3)

    if cfg.show_labels:
        for (x, y), lab in zip(coords, labels):
            ax.annotate(lab, (x, y), xytext=(8, 8),
                        textcoords="offset points",
                        fontsize=cfg.label_font_size, fontweight="bold",
                        zorder=5)

    ax.set_xlabel(cfg.xlabel or "MDS dimension 1")
    ax.set_ylabel(cfg.ylabel or "MDS dimension 2")
    if cfg.show_title:
        auto = (f"Detrital-zircon MDS  ·  {METRIC_LABEL.get(cfg.metric, cfg.metric)}"
                f"  ·  grains > {cfg.age_threshold:g} Ma")
        ax.set_title(cfg.title or auto)
    ax.grid(True, ls=":", color="0.8", alpha=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    if cfg.show_stress:
        fig.text(0.5, 0.015,
                 f"{cfg.mds_type} MDS   ·   stress = {stress:.3f}"
                 f"   ·   n = {len(labels)}",
                 ha="center", fontsize=8, color="0.4")
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    return fig, D, labels, stress


# ======================================================================
if __name__ == "__main__":
    rng = np.random.default_rng(3)
    rows = []
    for name, mix in {
        "WBS Fluvial": [(1050, 60, .18), (1750, 70, .25), (2720, 60, .57)],
        "WBS Eolian":  [(44, 2, .06), (1080, 45, .40), (1720, 55, .34), (2700, 70, .20)],
        "TSA 1":       [(35, 2, .04), (1090, 40, .45), (1700, 50, .31), (2650, 90, .20)],
        "TSA 2":       [(1085, 45, .48), (1715, 55, .32), (2710, 75, .20)],
        "TSA 3":       [(120, 15, .05), (1070, 55, .42), (1730, 60, .34), (2690, 85, .19)],
        "TSA 4":       [(1095, 50, .44), (1705, 65, .35), (2670, 80, .21)],
    }.items():
        for mu, sd, w in mix:
            rows.append(pd.DataFrame({"Sample_ID": name,
                                      "BestAge": rng.normal(mu, sd, int(120 * w))}))
    df = pd.concat(rows, ignore_index=True)
    df = df[df.BestAge > 1]
    with pd.ExcelWriter("example_data.xlsx") as w:
        df.to_excel(w, sheet_name="ZrUPb", index=False)

    S = all_samples(df, "Sample_ID")
    ages, labels = load_ages(df, S, "Sample_ID", "BestAge")
    f = make_kde_figure(ages, labels, KDEConfig())
    f.savefig("_selftest_kde.pdf")
    print("KDE figure ok  ", [f"{l}:{len(a)}" for l, a in zip(labels, ages)])

    ages2, labels2 = load_ages(df, S, "Sample_ID", "BestAge", threshold=70)
    f2, D, lab, stress = make_mds_figure(ages2, labels2, MDSConfig())
    f2.savefig("_selftest_mds.pdf")
    print(f"MDS figure ok   stress={stress:.4f}")
