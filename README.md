# DZ Studio

A small desktop app for making detrital-zircon figures such as kernel density (KDE)
plots with provenance pies, multidimensional scaling (MDS) maps, and
maximum depositional age (MDA) plots, without writing any code.

It began as a set of Jupyter notebooks and was wrapped in a window so that
people who don't program can pick samples, adjust settings with sliders and
dropdowns, and export publication-ready figures. The plotting and the
statistics live in plain, documented Python modules, so anyone who *does*
program can read, reuse, or modify them.

---

## What it does

- **KDE + pie** - one row per sample: a young-age panel, an old-age panel (or a
  single combined panel), the KDE curve scaled onto a histogram count axis,
  the area coloured by provenance age-bin, the main peaks labelled, and an
  optional provenance pie chart.
- **MDS** - turns sample-to-sample dissimilarity into a 2-D map where similar
  samples plot together, following Vermeesch (2013). K–S, cross-correlation, and
  likeness metrics; metric and non-metric scaling; the dissimilarity matrix
  exports to CSV.
- **MDA** - youngest single grain (YSG) and youngest grain clusters (YGC1σ,
  YGC2σ) per sample, plus a ranked-grain plot, following the detritalPy/
  Sharman & Malkowski (2020) definitions.

---
## Run without Python (pre-built GUI) (Preferred)

Non-programmers can use a pre-built GUI that bundles everything

**Go to the releases page and download the respective zip according to your Operating system**

- Download the zip for your OS from the **Releases** page, unzip it, and run
  `DZ Studio` inside. No Python needed.
- Or build it yourself once: run `build_windows.bat` (Windows) or
  `build_macos.sh` (macOS) from a machine that has Python. The result is a
  self-contained folder/app you can share.

Unsigned-app note: Windows may say "unknown publisher" (More info → Run anyway)
and macOS may block the first launch (right-click → Open → Open). 

---
## Install and run (from source) (If you have Python installed)

Download the repository and unzip everything into one folder. ## Requires Python 3.10 or newer ##

```
pip install -r requirements.txt
python dz_app.py
```

That opens the app. Open your Excel file, pick the sample, age, and error columns,
tick some samples, and the figure appears.

---

## The three tabs in detail

**KDE + pie** -- Young/old age ranges, the two KDE bandwidths and bin widths (or
a single combined panel with its own bandwidth), peak-label threshold, pie
range, pie size and label placement, column widths, row height, spacing and
every font size are controls on the right-hand panel.
You can choose the number of provenance groups, rename and color them accordingly from the 
bottom left panel. This will appear in the legend; you can turn it off as well.

**MDS** -- Age threshold, dissimilarity metric, KDE bandwidth (for the
correlation metrics), metric vs non-metric scaling, seed and restarts, point
size and colour, labels and neighbour lines. A "Consistent axes" group scales
every figure to the same square axes if you plan to plot multiple figures and 
maintain consistent scaling for all of them. This only reframes the plot; it 
does not change the statistics. Stress is reported in the caption.

**MDA** -- Choose the error column and whether it is 1σ or 2σ, and whether to
report at 1σ or 2σ. The table at the bottom right gives YSG, YGC1σ and YGC2σ as 
`age ± error` for every ticked sample and exports to CSV. Choose the desired method 
from the drop-down menu in the middle panel on the right-hand side. The rank plot shows 
the youngest grains with nested 1σ/2σ error bars; the grains falling into the chosen
cluster are highlighted. The MDA appears as a band. 
"Exclude the single youngest grain" drops an obvious outlier before every method runs.

---

## Methods and credits

The statistics follow published methods; this app is an interface to them, not
a new method.

- **MDS** follows **Vermeesch, P. (2013)**, *Multi-sample comparison of detrital
  age distributions*, Chemical Geology 341, 140–146. 
- **MDA** follows the youngest-grain-cluster definitions of **Sharman, G.R. &
  Malkowski, M.A. (2020)** and the **detritalPy** implementation
  (**Sharman et al., 2018**).
- The **K–S** dissimilarity is the Kolmogorov–Smirnov statistic.

If you use figures from this tool in a publication, please also cite the underlying
methods above.

---

## Project files (`.dzp`)

**Save project** writes a small `.dzp` file (plain JSON) holding the data-file
path, the chosen columns, the sample order, the provenance-bin scheme and every
setting on every tab. **Load project** restores all of it. 

---

## Preview vs export

The **preview** is a smooth image view, like a PDF reader: the figure is
rendered once and zoom/scroll never re-render it. Mouse wheel scrolls,
Ctrl+wheel zooms, the View dropdown fits width/page or a fixed percentage, and
F11 gives the figure the whole window. You can untick and stop the live preview if you want 
and use F5 or click redraw to regenerate the plot after you make any changes. 

**Export** (Ctrl+E) rebuilds the figure from the currently ticked samples, then
lets you pick format (PDF / SVG / PNG), size and resolution. Text is kept live and editable in
Illustrator.

---

## For developers - how the code is laid out

The project is deliberately split so that the maths and the interface never mix:

| File | What it is | Edit it to… |
|---|---|---|
| `dz_figs.py` | KDE and MDS: config objects, data loading, and the plotting/stats | change how the KDE or MDS looks or is computed |
| `dz_mda.py`  | MDA maths and the rank plot | change the MDA methods or their figure |
| `dz_app.py`  | the PySide6 window; no maths of its own | change the interface, controls, or export |

The rule of thumb: **if a number is wrong, it's in `dz_figs.py` / `dz_mda.py`;
if a window or control misbehaves, it's in `dz_app.py`.**

Each module starts with a docstring that walks through its sections, and the
functions carry docstrings with explanations. You can import the figure functions into a notebook and use them
without the GUI at all.

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
