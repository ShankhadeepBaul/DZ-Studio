"""
dz_mda.py  --  maximum depositional age
=======================================

Computes the maximum depositional age (MDA) of a sediment from its youngest
detrital-zircon grains -- the youngest age population puts an upper bound on
when the rock was deposited. This module holds the maths; dz_app.py's MDA tab
is the interface to it.

It follows the standard, published methods (matching detritalPy / Sharman &
Malkowski 2020), verified to reproduce their numbers exactly over thousands of
random test cases. Three MDAs are reported per sample:

    YSG     youngest single grain (by age + 1 sigma)
    YGC1s   youngest cluster of >= 2 grains overlapping at 1 sigma
    YGC2s   youngest cluster of >= 3 grains overlapping at 2 sigma

For the clusters, grains are sorted by (age + error) and the cluster is the
youngest run of grains whose error bars overlap the anchor grain; the MDA is
their inverse-variance weighted mean, reported with an MSWD (a scatter check:
~1 means the grains agree within their errors).

Convention: everything is computed internally at 1 sigma. The loader converts
a 2-sigma error column to 1 sigma on the way in, and results are reported back
at whichever sigma the caller asks for. Grains without a valid age AND error
are dropped, because an MDA is undefined without a real uncertainty per grain.

Run `python dz_mda.py` for a self-test with hand-checkable numbers.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


# ----------------------------------------------------------------------
#  Data preparation
#  Turning a spreadsheet column into clean 1-sigma errors and pulling the
#  ages + errors for one sample.
# ----------------------------------------------------------------------

def to_one_sigma(err, err_sigma):
    """Convert an error column to 1-sigma. err_sigma is 1 or 2."""
    e = np.asarray(err, float)
    return e / 2.0 if int(err_sigma) == 2 else e


def sample_ages_errors(df, id_col, age_col, err_col, sample, err_sigma=2):
    """
    Ages and 1-sigma errors for one sample, cleaned.

    Grains with a missing/zero/negative age or error are dropped, because MDA
    (unlike a KDE) is meaningless without a real uncertainty on each grain.
    Returns (ages, sig1) sorted youngest-first is NOT applied here; callers
    sort as needed.
    """
    import pandas as pd
    m = df[id_col].astype(str) == str(sample)
    a = pd.to_numeric(df.loc[m, age_col], errors="coerce").to_numpy(float)
    if err_col is None or err_col not in df.columns:
        e = np.full(a.shape, np.nan)
    else:
        e = pd.to_numeric(df.loc[m, err_col], errors="coerce").to_numpy(float)
    e = to_one_sigma(np.abs(e), err_sigma)
    ok = np.isfinite(a) & np.isfinite(e) & (a > 0) & (e > 0)
    return a[ok], e[ok]


# ----------------------------------------------------------------------
#  Building blocks
#  The pieces the methods are made of: the inverse-variance weighted mean
#  (with its MSWD), the youngest-cluster search, and a probability-density
#  curve. Read these to see exactly how a cluster is chosen and averaged.
# ----------------------------------------------------------------------

def weighted_mean(ages, sig1):
    """
    Inverse-variance weighted mean.

    Returns dict: mean, se1 (1-sigma of the mean), mswd, n.
    MSWD = (1/(n-1)) * sum( ((age-mean)/sigma)^2 )  over the grains used.
    """
    a = np.asarray(ages, float)
    s = np.asarray(sig1, float)
    n = a.size
    if n == 0:
        return {}
    w = 1.0 / s ** 2
    m = float(np.sum(w * a) / np.sum(w))
    se1 = float(np.sqrt(1.0 / np.sum(w)))
    mswd = float(np.sum(((a - m) / s) ** 2) / (n - 1)) if n > 1 else np.nan
    return {"mean": m, "se1": se1, "mswd": mswd, "n": n}


def youngest_cluster_indices(ages, sig1, n_min=3, k=2.0, contiguous=True):
    """
    Youngest grain cluster, following detritalPy / Sharman & Malkowski (2020).

    The grains are sorted by (age + k*sigma).  The cluster is anchored on the
    grain with the smallest (age + k*sigma); a later grain j joins the cluster
    if its lower bound falls below the anchor's upper bound:

        age_j - k*sig_j  <  age_anchor + k*sig_anchor

    With contiguous=True (the detritalPy default) the cluster is the run of
    consecutive grains from the anchor up to the first one that fails that
    test, and it is accepted only if that run is at least n_min long; otherwise
    the next anchor is tried.  With contiguous=False, any >= n_min grains that
    overlap the anchor qualify.

    Returns indices into the ORIGINAL arrays.  [] if no cluster is found.

    (YC1s is n_min=2, k=1; YC2s is n_min=3, k=2.)
    """
    a = np.asarray(ages, float)
    s = np.asarray(sig1, float) * k          # k-sigma half-widths
    n = a.size
    if n < n_min:
        return []
    order = np.argsort(a + s)                 # sort by age + k*sigma
    ao, so = a[order], s[order]
    tops = ao + so
    bottoms = ao - so

    # Walk candidate anchors from youngest (smallest age+k*sigma) upward. For
    # each anchor i, `overlaps` marks which later grains reach back below the
    # anchor's top. Contiguous mode stops at the first non-overlapping grain;
    # the run counts only if it is at least n_min long.
    for i in range(n):
        overlaps = bottoms[i:] < tops[i]      # each grain vs the anchor top
        if contiguous:
            false_idx = np.flatnonzero(~overlaps)
            if false_idx.size == 0:           # everything from i overlaps
                if n - i >= n_min:
                    return list(order[i:])
                continue
            elif false_idx[0] >= n_min:
                return list(order[i:i + false_idx[0]])
        else:
            if int(overlaps.sum()) >= n_min:
                return list(order[i:][overlaps])
    return []


def pdp(ages, sig1, x):
    """Probability density plot (each grain a gaussian of its own 1-sigma)."""
    a = np.asarray(ages, float)
    s = np.asarray(sig1, float)
    if a.size == 0:
        return np.zeros_like(x)
    s = np.where(s <= 0, np.nan, s)
    d = (x[:, None] - a[None, :]) / s[None, :]
    y = np.nansum(np.exp(-0.5 * d ** 2) / (s[None, :] * np.sqrt(2 * np.pi)),
                  axis=1) / a.size
    return y


def youngest_pdp_peak(ages, sig1):
    a = np.asarray(ages, float)
    if a.size == 0:
        return np.nan
    lo = max(0.0, a.min() - 50)
    hi = a.min() + 400
    x = np.arange(lo, hi, 0.25)
    y = pdp(a, sig1, x)
    if y.size < 3:
        return float(a.min())
    peaks = np.flatnonzero((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])) + 1
    peaks = peaks[y[peaks] >= 0.05 * y.max()]
    return float(x[peaks[0]]) if peaks.size else float(a.min())


# ----------------------------------------------------------------------
#  The methods
#  compute() runs every MDA method on one sample and returns their ages,
#  uncertainties, MSWD and cluster membership. table_row() formats those
#  into one tidy row for the results table.
# ----------------------------------------------------------------------

@dataclass
class MDAParams:
    n_min: int = 3            # grains needed in the cluster (YC2s convention)
    k_sigma: float = 2.0      # overlap test width
    contiguous: bool = True   # detritalPy default: cluster grains adjacent
    exclude_youngest: bool = False
    report_sigma: int = 2     # 1 or 2, for the numbers handed back


def compute(ages, sig1, params: MDAParams) -> dict:
    """
    Run every method on one sample.  Returns a dict of results; each value is
    itself a dict with at least 'age' and 'err' (err at params.report_sigma),
    plus method-specific extras (mswd, n_cluster, cluster indices...).

    If exclude_youngest is on, the single youngest grain is removed before ANY
    method runs, so YSG, the cluster search and YPP all ignore it.
    """
    a = np.asarray(ages, float)
    s = np.asarray(sig1, float)
    rs = params.report_sigma

    excluded = None
    if params.exclude_youngest and a.size > 0:
        j = int(np.argmin(a))
        excluded = (float(a[j]), float(s[j]))
        keep = np.ones(a.size, bool)
        keep[j] = False
        a, s = a[keep], s[keep]

    out = {"n_used": int(a.size), "excluded_grain": excluded,
           "report_sigma": rs}
    if a.size == 0:
        return out

    # --- YSG ---------------------------------------------------------
    # detritalPy selects the youngest grain by (age + 1-sigma), not youngest
    # age, then reports that grain's age and error.
    j = int(np.argmin(a + s))
    out["YSG"] = {"age": float(a[j]), "err": float(s[j] * rs), "n": 1}

    # --- clusters ----------------------------------------------------
    def cluster_result(n_min, k, label):
        idx = youngest_cluster_indices(a, s, n_min=n_min, k=k,
                                       contiguous=params.contiguous)
        if not idx:
            return {"age": np.nan, "err": np.nan, "mswd": np.nan,
                    "n": 0, "idx": []}
        wm = weighted_mean(a[idx], s[idx])
        return {"age": wm["mean"], "err": wm["se1"] * rs,
                "mswd": wm["mswd"], "n": wm["n"], "idx": idx}

    # YC1s: >=2 grains at 1 sigma
    out["YC1s"] = cluster_result(2, 1.0, "YC1s")
    # YC2s: user's main definition (default n_min/k come from params)
    out["YC2s"] = cluster_result(params.n_min, params.k_sigma, "YC2s")

    # --- YPP ---------------------------------------------------------
    out["YPP"] = {"age": youngest_pdp_peak(a, s), "err": np.nan, "n": int(a.size)}

    return out


METHODS = ["YSG", "YC1s", "YC2s"]
METHOD_LABEL = {
    "YSG": "YSG — youngest single grain",
    "YC1s": "YGC1σ — youngest ≥2 grains overlapping at 1σ",
    "YC2s": "YGC2σ — youngest ≥3 grains overlapping at 2σ",
}
METHOD_SHORT = {"YSG": "YSG", "YC1s": "YGC1σ", "YC2s": "YGC2σ"}


def _fmt(age, err):
    if age is None or not np.isfinite(age):
        return "—"
    if err is None or not np.isfinite(err):
        return f"{age:.1f}"
    return f"{age:.1f} ± {err:.1f}"


def table_row(sample, ages, sig1, params: MDAParams) -> dict:
    """
    One row for the summary table: the three reported MDAs as 'age ± error'
    (at params.report_sigma).  MSWD and cluster size live on the figure, not
    here, so the table stays readable.
    """
    r = compute(ages, sig1, params)
    rs = params.report_sigma
    row = {"Sample": sample, "n": r.get("n_used", 0)}
    labels = [("YSG", "YSG"), ("YC1s", "YGC1σ"), ("YC2s", "YGC2σ")]
    for key, col in labels:
        m = r.get(key, {})
        row[f"{col} (±{rs}σ)"] = _fmt(m.get("age"), m.get("err"))
    return row


# ----------------------------------------------------------------------
#  Self-test
#  Runs the methods on a tiny hand-picked dataset so you can confirm the
#  numbers by eye. Executed when you run `python dz_mda.py` directly.
# ----------------------------------------------------------------------

if __name__ == "__main__":
    # three grains tightly clustered at ~40, plus older grains and one young
    # outlier at 34 that should NOT drag a >=3 cluster.
    ages = np.array([34.0, 39.5, 40.0, 40.5, 55.0, 60.0, 200.0, 1000.0])
    err2 = np.array([1.0,  1.2,  1.0,  1.1,  2.0,  2.0,  5.0,   20.0])  # 2-sigma
    sig1 = to_one_sigma(err2, 2)

    print("== include youngest ==")
    r = compute(ages, sig1, MDAParams(n_min=3, k_sigma=2.0,
                                      exclude_youngest=False, report_sigma=2))
    for m in METHODS:
        d = r[m]
        print(f"  {m:5s} age={d['age']:7.2f}  ±{d.get('err', float('nan')):.2f}"
              f"  n={d.get('n','-')}  mswd={d.get('mswd', float('nan')):.2f}"
              if np.isfinite(d.get('mswd', np.nan)) else
              f"  {m:5s} age={d['age']:7.2f}  n={d.get('n','-')}")
    print("  YC2s cluster grains:", [round(ages[i], 1) for i in r['YC2s']['idx']])

    print("== exclude youngest (drops the 34 outlier) ==")
    r2 = compute(ages, sig1, MDAParams(exclude_youngest=True, report_sigma=2))
    print("  excluded grain:", r2["excluded_grain"])
    print("  YSG now:", round(r2["YSG"]["age"], 2))
    print("  YC2s cluster grains:", [round(ages[i], 1) for i in
          youngest_cluster_indices(*(lambda a,s:(a,s))(
              np.delete(ages, np.argmin(ages)),
              np.delete(sig1, np.argmin(ages))), n_min=3, k=2.0)])

    # overlap sanity: 39.5,40,40.5 at 2sigma -> bars ~[38.3,40.7],[38,42],[38.9,42.1]
    # common overlap [38.9,40.7] non-empty -> cluster of 3. good.
    print("== weighted mean of the 3 clustered grains ==")
    idx = r["YC2s"]["idx"]
    print("  ", weighted_mean(ages[idx], sig1[idx]))


# ----------------------------------------------------------------------
#  Rank plot
#  The figure people put in papers: the youngest grains ranked left to
#  right with their 1-sigma and 2-sigma error bars, the chosen cluster
#  highlighted, and a band at the MDA. Drawing only; the maths is above.
# ----------------------------------------------------------------------

def make_rank_plot(fig, sample, ages, sig1, params: MDAParams,
                   method="YC2s", n_show=25, age_max=None,
                   point_color="#4C72B0", cluster_color="#C44E52", fs=9,
                   title="", xlabel="", ylabel="", show_legend=True,
                   label_in="Grains in cluster",
                   label_out="Grains not in cluster"):
    """
    Ranked-age (weighted-mean) plot for a single sample.

    Grains sorted youngest->oldest on x; age +/- (report_sigma) bars on y.
    The grains in the chosen cluster are highlighted; a horizontal band marks
    the MDA +/- its uncertainty.  This is the figure people put in the paper.
    """
    import matplotlib as mpl
    mpl.rcParams.update({"font.family": "serif", "pdf.fonttype": 42,
                         "ps.fonttype": 42, "svg.fonttype": "none"})
    fig.clear()
    ax = fig.add_subplot(111)

    a = np.asarray(ages, float)
    s = np.asarray(sig1, float)
    if a.size == 0:
        ax.text(0.5, 0.5, f"{sample}: no usable grains\n"
                "(need age and error > 0)", ha="center", va="center",
                color="#888", transform=ax.transAxes)
        ax.axis("off")
        return

    res = compute(a, s, params)
    md = res.get(method, {})

    # work on the youngest grains for a readable plot; optionally drop grains
    # older than age_max (they are irrelevant to the MDA and only squash the
    # y-axis)
    order = np.argsort(a)
    a_s, s_s = a[order], s[order]
    if age_max is not None:
        keep = a_s <= age_max
        if keep.sum() >= 2:
            a_s, s_s = a_s[keep], s_s[keep]
    show = min(n_show, a_s.size)
    xs = np.arange(1, show + 1)
    ays, sys = a_s[:show], s_s[:show]
    rs = params.report_sigma

    # which of the shown grains are in the cluster?
    in_cluster = np.zeros(show, bool)
    if method in ("YC1s", "YC2s") and md.get("idx"):
        # md['idx'] indexes the (possibly youngest-excluded) array used inside
        # compute(); rebuild that array's ages to match back by value+order
        cl_ages = set(np.round(
            (a if not params.exclude_youngest
             else np.delete(a, np.argmin(a)))[md["idx"]], 6))
        in_cluster = np.array([round(v, 6) in cl_ages for v in ays])

    base = ~in_cluster
    from matplotlib.lines import Line2D
    from matplotlib.colors import to_rgb

    def _lighten(c, f=0.55):
        r, g, b = to_rgb(c)
        return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)

    def draw_grains(mask, color, z):
        if not mask.any():
            return
        xg, yg, sg = xs[mask], ays[mask], sys[mask]
        soft = _lighten(color, 0.55)
        # 2σ: soft, thin whisker with small caps (the outer extent)
        _, caps, bars = ax.errorbar(
            xg, yg, yerr=2 * sg, fmt="none", ecolor=soft,
            elinewidth=1.2, capsize=2.5, capthick=1.2, zorder=z)
        # 1σ: bold, rounded segment sitting on top (the inner, darker bar)
        _, _, bars1 = ax.errorbar(
            xg, yg, yerr=1 * sg, fmt="none", ecolor=color,
            elinewidth=3.2, capsize=0, zorder=z + 0.1)
        for b in bars1:
            b.set_capstyle("round")
        # marker with a white edge so it reads on top of the bar
        ax.plot(xg, yg, "o", ms=5.2, color=color, mec="white", mew=0.9,
                zorder=z + 0.2)

    draw_grains(base, point_color, 3)
    draw_grains(in_cluster, cluster_color, 4)

    # MDA line + uncertainty band
    if np.isfinite(md.get("age", np.nan)):
        mda = md["age"]
        err = md.get("err", np.nan)
        ax.axhline(mda, color=cluster_color, lw=1.1, ls=(0, (5, 2)), zorder=2)
        if np.isfinite(err):
            ax.axhspan(mda - err, mda + err, color=cluster_color,
                       alpha=0.12, zorder=1)

    ax.set_xlabel(xlabel or "Grain rank (youngest → oldest)", fontsize=fs)
    ax.set_ylabel(ylabel or "Age (Ma)", fontsize=fs)
    ax.tick_params(labelsize=fs - 1)
    ax.margins(x=0.04)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    short = {"YSG": "YSG", "YC1s": "YGC1σ", "YC2s": "YGC2σ"}.get(method, method)
    auto = f"{sample}   ·   {short}"
    if np.isfinite(md.get("age", np.nan)):
        auto += f" = {md['age']:.1f}"
        if np.isfinite(md.get("err", np.nan)):
            auto += f" ± {md['err']:.1f} ({rs}σ)"
        if np.isfinite(md.get("mswd", np.nan)):
            auto += f",  MSWD = {md['mswd']:.2f},  n = {md.get('n','-')}"
    ax.set_title(title or auto, fontsize=fs)

    if show_legend:
        grey = "#555555"
        handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=cluster_color,
                   markeredgecolor="white", markersize=7, label=label_in),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=point_color,
                   markeredgecolor="white", markersize=7, label=label_out),
            Line2D([0], [0], color=grey, lw=3.2, solid_capstyle="round",
                   label="1σ"),
            Line2D([0], [0], color=_lighten(grey, 0.55), lw=1.2, label="2σ"),
        ]
        ax.legend(handles=handles, fontsize=fs - 2, frameon=False,
                  loc="lower right", ncol=2, handlelength=1.4,
                  columnspacing=1.2)
    fig.tight_layout()
