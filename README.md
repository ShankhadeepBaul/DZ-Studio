# DZ Studio

A small desktop app for making detrital-zircon figures — kernel-density (KDE)
plots with provenance pies, multidimensional-scaling (MDS) maps, and
maximum-depositional-age (MDA) plots — without writing any code.

It began as a set of Jupyter notebooks and was wrapped in a window so that
people who don't program can pick samples, adjust settings with sliders and
dropdowns, and export publication-ready figures. The plotting and the
statistics live in plain, documented Python modules, so anyone who *does*
program can read, reuse, or modify them.

---

## What it does

- **KDE + pie** — one row per sample: a young-age panel, an old-age panel (or a
  single combined panel), the KDE curve scaled onto a histogram count axis,
  the area coloured by provenance age-bin, the main peaks labelled, and an
  optional provenance pie chart.
- **MDS** — turns sample-to-sample dissimilarity into a 2-D map where similar
  samples plot together, following Vermeesch (2013). K–S, cross-correlation and
  likeness metrics; metric and non-metric scaling; the dissimilarity matrix
  exports to CSV.
- **MDA** — youngest single grain (YSG) and youngest grain clusters (YGC1σ,
  YGC2σ) per sample, plus a ranked-grain plot, following the detritalPy /
  Sharman & Malkowski (2020) definitions.

---

## Install and run (from source)

Requires Python 3.10 or newer.

```
pip install -r requirements.txt
python dz_app.py
```

That opens the app. Open your Excel file, pick the sample and age columns,
tick some samples, and the figure appears.

## Run without Python (built app)

Non-programmers can use a pre-built app that bundles everything:

- Download the zip for your OS from the **Releases** page, unzip it, and run
  `DZ Studio` inside. No Python needed.
- Or build it yourself once: run `build_windows.bat` (Windows) or
  `build_macos.sh` (macOS) from a machine that has Python. The result is a
  self-contained folder/app you can share.

Unsigned-app note: Windows may say "unknown publisher" (More info → Run anyway)
and macOS may block the first launch (right-click → Open → Open). Normal for
un-notarised software.

---

## The three tabs in detail

**KDE + pie.** Young/old age ranges, the two KDE bandwidths and bin widths (or
a single combined panel with its own bandwidth), peak-label threshold, pie
range, pie size and label placement, column widths, row height, spacing and
every font size are controls on the right-hand panel.

**MDS.** Age threshold, dissimilarity metric, KDE bandwidth (for the
correlation metrics), metric vs non-metric scaling, seed and restarts, point
size and colour, labels and neighbour lines. A "Consistent axes" group scales
every figure to the same square axes so separate panels line up — this only
reframes the plot, it does not change the statistics. Stress is reported in the
caption (near 0 = the map reproduces the dissimilarities well).

**MDA.** Choose the error column and whether it is 1σ or 2σ, and whether to
report at 1σ or 2σ. The table gives YSG, YGC1σ and YGC2σ as `age ± error` for
every ticked sample and exports to CSV. The rank plot shows the youngest grains
with nested 1σ/2σ error bars, the chosen cluster highlighted, and a band at the
MDA. "Exclude the single youngest grain" drops an obvious outlier before every
method runs.

---

## Methods and credits

The statistics follow published methods; this app is an interface to them, not
a new method.

- **MDS** follows **Vermeesch, P. (2013)**, *Multi-sample comparison of detrital
  age distributions*, Chemical Geology 341, 140–146. Non-metric MDS is seeded
  from the classical (`cmdscale`) configuration the way R's `MASS::isoMDS` does,
  so the result is deterministic and reproduces the structure of the R
  `provenance` package (up to the usual free rotation/reflection of any MDS
  map).
- **MDA** follows the youngest-grain-cluster definitions of **Sharman, G.R. &
  Malkowski, M.A. (2020)** and the **detritalPy** implementation
  (**Sharman et al., 2018**). The YSG / YGC1σ / YGC2σ results were checked
  against detritalPy over thousands of random samples and match exactly.
- The **K–S** dissimilarity is the Kolmogorov–Smirnov statistic (`scipy.stats`).

If you use figures from this tool in a publication, please cite the underlying
methods above.

---

## Project files (`.dzp`)

**Save project** writes a small `.dzp` file (plain JSON) holding the data-file
path, the chosen columns, the sample order, the provenance-bin scheme and every
setting on every tab. **Load project** restores all of it. This is how you
regenerate the exact same figure months later — during revision, say — instead
of trying to remember which bandwidth you used. The provenance-bin scheme can
also be saved and loaded on its own, so one colour scheme can serve a whole lab.

---

## Preview vs export

The **preview** is a smooth image view, like a PDF reader: the figure is
rendered once and zoom/scroll never re-render it. Mouse wheel scrolls,
Ctrl+wheel zooms, the View dropdown fits width/page or a fixed percentage, and
F11 gives the figure the whole window. The preview only changes how big the
image looks — never the figure itself.

**Export** (Ctrl+E) rebuilds the figure from the currently ticked samples, then
lets you pick format (PDF / SVG / PNG), printed size (with journal-column
presets) and, for PNG, resolution. Text is kept live and editable in
Illustrator (`pdf.fonttype = 42`, `svg.fonttype = none`).

---

## For developers — how the code is laid out

The project is deliberately split so that the maths and the interface never mix:

| File | What it is | Edit it to… |
|---|---|---|
| `dz_figs.py` | KDE and MDS: config objects, data loading, and the plotting/stats | change how the KDE or MDS looks or is computed |
| `dz_mda.py`  | MDA maths and the rank plot | change the MDA methods or their figure |
| `dz_app.py`  | the PySide6 window; no maths of its own | change the interface, controls, or export |

The rule of thumb: **if a number is wrong, it's in `dz_figs.py` / `dz_mda.py`;
if a window or control misbehaves, it's in `dz_app.py`.**

Each module starts with a docstring that walks through its sections, and the
functions carry docstrings explaining inputs, outputs and the tricky bits
(for example, exactly how an MDA cluster is chosen, or why MDS axes carry no
meaning). You can import the figure functions into a notebook and use them
without the GUI at all — that is the quickest way to experiment with a change.

Every setting is a field on a config dataclass (`KDEConfig`, `MDSConfig` in
`dz_figs.py`; `MDAParams` in `dz_mda.py`) with a sensible default, so adding a
new option is usually: add a field to the config, read it in the plotting
function, and add one widget in `dz_app.py` that sets it.

Contributions and forks are welcome under the MIT licence.

---

## Requirements

See `requirements.txt`. In short: numpy, pandas, scipy, scikit-learn,
matplotlib, openpyxl, and PySide6. `example_data.xlsx` is a synthetic dataset
you can open to try the app immediately.
