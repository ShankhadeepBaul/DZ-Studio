"""
dz_app.py  --  the DZ Studio desktop app
=========================================

This is the graphical front end. It contains no geochronology and no plotting maths of its own: it reads an Excel file, lets you pick samples and set options with widgets, then calls the real work in dz_figs.py (KDE + pie,
MDS) and dz_mda.py (maximum depositional age) and shows the result.

That split is the thing to understand before reading this file: if a NUMBER is wrong, look in dz_figs.py / dz_mda.py; if the WINDOW or a control misbehaves,
look here.

Layout of the file, top to bottom:

  - helpers        small widget factories (labelled spin boxes, text fields,
                   collapsible groups) reused everywhere.
  - ColorButton    the little colour swatch used for bins and MDS points.
  - SampleList     the tick-list of samples on the left, with reordering.
  - BinEditor      the provenance-bin table (age ranges, labels, colours).
  - Tab            a base class: each tab is a figure preview on the left and a
                   scrollable column of settings on the right, with shared
                   redraw/zoom/export behaviour.
  - KDETab / MDSTab / MDATab   the three tabs, each exposing its figure's
                   options and calling the matching function in dz_figs/dz_mda.
  - MainWindow     wires it all together: the menu, the data dock, saving and
                   loading .dzp project files, and image export.

Run:  python dz_app.py
"""

from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.figure import Figure

from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import Qt

import dz_figs as F
import dz_mda as MDA

APP = "DZ Studio"
VERSION = "0.2"


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------

def lineedit(placeholder=""):
    e = QtWidgets.QLineEdit()
    e.setPlaceholderText(placeholder)
    e.setClearButtonEnabled(True)
    return e


def spin(lo, hi, val, step=1.0, dec=0, suffix=""):
    w = QtWidgets.QSpinBox() if dec == 0 else QtWidgets.QDoubleSpinBox()
    w.setRange(lo, hi)
    w.setSingleStep(step)
    if dec:
        w.setDecimals(dec)
    w.setValue(val)
    if suffix:
        w.setSuffix(suffix)
    return w


def form(rows):
    f = QtWidgets.QFormLayout()
    f.setLabelAlignment(Qt.AlignRight)
    f.setContentsMargins(6, 6, 6, 6)
    f.setSpacing(4)
    for lab, w in rows:
        f.addRow(w) if lab is None else f.addRow(lab, w)
    box = QtWidgets.QWidget()
    box.setLayout(f)
    return box


def group(title, w, collapsed=False):
    g = QtWidgets.QGroupBox(title)
    g.setCheckable(True)
    g.setChecked(not collapsed)
    lay = QtWidgets.QVBoxLayout(g)
    lay.setContentsMargins(4, 4, 4, 4)
    lay.addWidget(w)
    g.toggled.connect(w.setVisible)
    w.setVisible(not collapsed)
    return g


class ColorButton(QtWidgets.QPushButton):
    changed = QtCore.Signal(str)
    _n = 0

    def __init__(self, color):
        super().__init__()
        self.setFixedSize(30, 18)
        self.setCursor(Qt.PointingHandCursor)
        ColorButton._n += 1
        self.setObjectName(f"swatch{ColorButton._n}")
        self._c = color
        self._paint()
        self.clicked.connect(self._pick)

    def _paint(self):
        self.setStyleSheet(
            f"QPushButton#{self.objectName()} {{"
            f" background-color: {self._c};"
            f" border: 1px solid rgba(0,0,0,0.35);"
            f" border-radius: 3px; }}"
            f"QPushButton#{self.objectName()}:hover {{"
            f" border: 1px solid rgba(0,0,0,0.75); }}")
        self.setToolTip(f"{self._c}   (click to change)")

    def color(self):
        return self._c

    def _pick(self):
        c = QtWidgets.QColorDialog.getColor(
            QtGui.QColor(self._c), self.window(), "Select colour")
        if c.isValid():
            self._c = c.name()
            self._paint()
            self.changed.emit(self._c)


# ----------------------------------------------------------------------
#                   sample list: tick + drag to reorder  
# ----------------------------------------------------------------------

class SampleList(QtWidgets.QWidget):
    """The tick-list of samples on the left of the window. Lets the user
    choose which samples appear in the figures and drag them into the order
    they should be drawn. Emits `changed` whenever the selection or order
    changes so the current tab can redraw."""
    changed = QtCore.Signal()

    def __init__(self):
        super().__init__()
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4)
        v.addWidget(QtWidgets.QLabel(
            "Tick to include. Drag to set top-to-bottom order."))
        self.list = QtWidgets.QListWidget()
        self.list.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
        self.list.itemChanged.connect(lambda *_: self.changed.emit())
        self.list.model().rowsMoved.connect(lambda *_: self.changed.emit())
        v.addWidget(self.list, 1)
        row = QtWidgets.QHBoxLayout()
        for t, fn in (("All", lambda: self._all(True)),
                      ("None", lambda: self._all(False))):
            b = QtWidgets.QPushButton(t)
            b.clicked.connect(fn)
            row.addWidget(b)
        v.addLayout(row)

    def populate(self, names, counts):
        self.list.blockSignals(True)
        self.list.clear()
        for s in names:
            it = QtWidgets.QListWidgetItem(f"{s}   (n={counts[s]})")
            it.setData(Qt.UserRole, s)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked)
            self.list.addItem(it)
        self.list.blockSignals(False)
        self.changed.emit()

    def _all(self, on):
        self.list.blockSignals(True)
        for i in range(self.list.count()):
            self.list.item(i).setCheckState(Qt.Checked if on else Qt.Unchecked)
        self.list.blockSignals(False)
        self.changed.emit()

    def selected(self):
        return [self.list.item(i).data(Qt.UserRole)
                for i in range(self.list.count())
                if self.list.item(i).checkState() == Qt.Checked]

    def set_state(self, order, checked):
        lw = self.list
        lw.blockSignals(True)
        cur = {lw.item(i).data(Qt.UserRole): lw.item(i).text()
               for i in range(lw.count())}
        lw.clear()
        for s in order:
            if s not in cur:
                continue
            it = QtWidgets.QListWidgetItem(cur[s])
            it.setData(Qt.UserRole, s)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if s in checked else Qt.Unchecked)
            lw.addItem(it)
        lw.blockSignals(False)


# ----------------------------------------------------------------------
# provenance bin editor  (this replaces AGE_BINS / BIN_LABELS / BIN_COLORS)
# ----------------------------------------------------------------------

class BinEditor(QtWidgets.QWidget):
    """The provenance-bin table: each row is an age range with a label and a
    colour, used to colour the KDE fill and slice the pie. Add or remove rows
    to use as many or as few bins as you like. Emits `changed` on any edit."""
    changed = QtCore.Signal()

    def __init__(self):
        super().__init__()
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4)
        self.t = QtWidgets.QTableWidget(0, 4)
        self.t.setHorizontalHeaderLabels(["From", "To", "Label", "Colour"])
        self.t.horizontalHeader().setSectionResizeMode(
            2, QtWidgets.QHeaderView.Stretch)
        self.t.setColumnWidth(0, 52)
        self.t.setColumnWidth(1, 52)
        self.t.setColumnWidth(3, 48)
        self.t.itemChanged.connect(lambda *_: self.changed.emit())
        v.addWidget(self.t, 1)
        row = QtWidgets.QHBoxLayout()
        for t, fn, tip in (
                ("＋ Add bin", self.add, "Add a provenance bin (row)"),
                ("− Remove bin", self.rm, "Remove the selected bin"),
                ("Reset", self.reset, "Restore the default bin scheme"),
                ("Load", self.load, "Load a saved bin scheme"),
                ("Save", self.save, "Save this bin scheme to a file")):
            b = QtWidgets.QPushButton(t)
            b.setToolTip(tip)
            b.clicked.connect(fn)
            row.addWidget(b)
        v.addLayout(row)
        self.reset()

    def _row(self, lo, hi, lab, col):
        r = self.t.rowCount()
        self.t.blockSignals(True)
        self.t.insertRow(r)
        for c, val in ((0, lo), (1, hi), (2, lab)):
            self.t.setItem(r, c, QtWidgets.QTableWidgetItem(str(val)))
        cb = ColorButton(col)
        cb.changed.connect(lambda *_: self.changed.emit())
        self.t.setCellWidget(r, 3, cb)
        self.t.blockSignals(False)

    def reset(self):
        self.t.blockSignals(True)
        self.t.setRowCount(0)
        self.t.blockSignals(False)
        for (lo, hi), lab, col in zip(F.DEFAULT_AGE_BINS,
                                      F.DEFAULT_BIN_LABELS,
                                      F.DEFAULT_BIN_COLORS):
            self._row(lo, hi, lab, col)
        self.changed.emit()

    def add(self):
        self._row(0, 100, "new bin", "#999999")
        self.changed.emit()

    def rm(self):
        r = self.t.currentRow()
        if r >= 0:
            self.t.removeRow(r)
            self.changed.emit()

    def values(self):
        bins, labels, colors = [], [], []
        for r in range(self.t.rowCount()):
            try:
                lo = float(self.t.item(r, 0).text())
                hi = float(self.t.item(r, 1).text())
                lab = self.t.item(r, 2).text()
                col = self.t.cellWidget(r, 3).color()
            except Exception:
                continue
            bins.append((lo, hi))
            labels.append(lab)
            colors.append(col)
        return bins, labels, colors

    def save(self):
        p, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save bin scheme", "age_bins.json", "JSON (*.json)")
        if p:
            b, l, c = self.values()
            Path(p).write_text(json.dumps(
                {"bins": b, "labels": l, "colors": c}, indent=2))

    def load(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Load bin scheme", "", "JSON (*.json)")
        if not p:
            return
        d = json.loads(Path(p).read_text())
        self.set_values(d["bins"], d["labels"], d["colors"])

    def set_values(self, bins, labels, colors):
        self.t.blockSignals(True)
        self.t.setRowCount(0)
        self.t.blockSignals(False)
        for (lo, hi), lab, col in zip(bins, labels, colors):
            self._row(lo, hi, lab, col)
        self.changed.emit()


# ----------------------------------------------------------------------
# base tab
# ----------------------------------------------------------------------

def render_figure_to_pixmap(fig, dpi):
    """Render a matplotlib Figure straight to a QPixmap at the given DPI."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fig.set_dpi(dpi)
    c = FigureCanvasAgg(fig)
    c.draw()
    w, h = c.get_width_height()
    buf = c.buffer_rgba()
    img = QtGui.QImage(bytes(buf), w, h,
                       QtGui.QImage.Format_RGBA8888).copy()
    return QtGui.QPixmap.fromImage(img)


class FigureView(QtWidgets.QGraphicsView):
    """
    Smooth preview.  The figure is rendered once to a pixmap; zoom and scroll
    are pure Qt view transforms and never re-render, so there is nothing to
    glitch.  Mouse wheel scrolls; Ctrl+wheel zooms toward the cursor.
    """
    zoomChanged = QtCore.Signal(float)   # emits effective dpi

    def __init__(self, base_dpi=96, render_dpi=150):
        super().__init__()
        self.base_dpi = base_dpi
        self.render_dpi = render_dpi
        self._scene = QtWidgets.QGraphicsScene(self)
        self.setScene(self._scene)
        self.item = None
        self.fit_mode = "width"          # "width" | "page" | None(=fixed zoom)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.AnchorViewCenter)
        self.setRenderHint(QtGui.QPainter.SmoothPixmapTransform, True)
        self.setRenderHint(QtGui.QPainter.Antialiasing, True)
        self.setBackgroundBrush(QtGui.QColor("#e9e9ee"))
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setAlignment(Qt.AlignCenter)

    def set_pixmap(self, pm):
        self._scene.clear()
        self.item = self._scene.addPixmap(pm)
        self.item.setTransformationMode(Qt.SmoothTransformation)
        self._scene.setSceneRect(QtCore.QRectF(pm.rect()))
        if self.fit_mode:
            self.apply_fit()
        else:
            self._emit()

    def _scale_now(self):
        return self.transform().m11()

    def _set_scale(self, s):
        self.setTransform(QtGui.QTransform().scale(s, s))

    def apply_fit(self):
        if self.item is None:
            return
        pm = self.item.pixmap()
        pw, ph = pm.width(), pm.height()
        if pw == 0 or ph == 0:
            return
        vp = self.viewport().size()
        aw, ah = vp.width() - 4, vp.height() - 4
        if self.fit_mode == "width":
            s = aw / pw
        elif self.fit_mode == "page":
            s = min(aw / pw, ah / ph)
        else:
            return
        self._set_scale(max(s, 0.01))
        self._emit()

    def set_zoom_percent(self, pct):
        self.fit_mode = None
        self._set_scale(self.base_dpi / self.render_dpi * pct / 100.0)
        self._emit()

    def set_fit(self, mode):
        self.fit_mode = mode
        self.apply_fit()

    def wheelEvent(self, e):
        if e.modifiers() & Qt.ControlModifier:
            self.fit_mode = None
            f = 1.15 if e.angleDelta().y() > 0 else 1 / 1.15
            self.scale(f, f)
            self._emit()
            e.accept()
        else:
            super().wheelEvent(e)        # normal vertical scroll

    def resizeEvent(self, e):
        super().resizeEvent(e)
        if self.fit_mode:
            self.apply_fit()

    def effective_dpi(self):
        return self._scale_now() * self.render_dpi

    def _emit(self):
        self.zoomChanged.emit(self.effective_dpi())


class Tab(QtWidgets.QWidget):
    """Base class for the three tabs. Handles everything they share.

    Layout: a smooth image preview on the left, a scrollable column of settings
    on the right (built with self.add(...) in each subclass, closed with
    self.done()).

    How drawing works, which matters if you add or change a tab:

      * Each subclass implements draw(), which builds self.fig by calling the
        relevant dz_figs/dz_mda function. It should NOT save or display -- just
        populate self.fig.
      * redraw() calls draw() then renders self.fig once to an image the
        preview shows. Zooming and scrolling transform that image; they never
        re-run draw(), so the preview is fast and never flickers.
      * The figure's true print size (in inches) is fixed and is what gets
        exported. The preview only changes how big that image looks on screen,
        never the figure itself.

    So: to add a control, add its widget in the subclass __init__, read it in
    that subclass's config()/params, and the shared redraw()/export() here take
    care of the rest.
    """
    ZOOMS = ["Fit width", "Fit page", "50 %", "75 %", "100 %",
             "150 %", "200 %", "300 %"]
    RENDER_DPI = 150

    def __init__(self, app):
        super().__init__()
        self.app = app
        self.fig = Figure(figsize=(9, 7), facecolor="white")

        self.view = FigureView(render_dpi=self.RENDER_DPI)
        self.view.zoomChanged.connect(self._show_zoom)

        self.zoom = QtWidgets.QComboBox()
        self.zoom.addItems(self.ZOOMS)
        self.zoom.setCurrentText("Fit width")
        self.zoom.currentIndexChanged.connect(self._zoom_changed)
        self.zoom_lbl = QtWidgets.QLabel("")
        self.zoom_lbl.setStyleSheet("color:#666;")
        hint = QtWidgets.QLabel("wheel = scroll · Ctrl+wheel = zoom")
        hint.setStyleSheet("color:#999;")

        bar = QtWidgets.QHBoxLayout()
        bar.setContentsMargins(4, 2, 6, 2)
        bar.addWidget(QtWidgets.QLabel("View:"))
        bar.addWidget(self.zoom)
        bar.addWidget(self.zoom_lbl)
        bar.addStretch(1)
        bar.addWidget(hint)

        self.panel = QtWidgets.QWidget()
        self.pl = QtWidgets.QVBoxLayout(self.panel)
        self.pl.setContentsMargins(2, 2, 2, 2)
        self.pl.setSpacing(4)
        self.settings_scroll = QtWidgets.QScrollArea()
        self.settings_scroll.setWidgetResizable(True)
        self.settings_scroll.setWidget(self.panel)
        self.settings_scroll.setMinimumWidth(300)
        self.settings_scroll.setMaximumWidth(360)

        left = QtWidgets.QWidget()
        lv = QtWidgets.QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.addLayout(bar)
        lv.addWidget(self.view, 1)

        sp = QtWidgets.QSplitter()
        sp.addWidget(left)
        sp.addWidget(self.settings_scroll)
        sp.setStretchFactor(0, 1)
        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(sp)

    def add(self, w):
        self.pl.addWidget(w)

    def done(self):
        self.pl.addStretch(1)

    def _zoom_changed(self):
        t = self.zoom.currentText()
        if t == "Fit width":
            self.view.set_fit("width")
        elif t == "Fit page":
            self.view.set_fit("page")
        else:
            self.view.set_zoom_percent(float(t.rstrip(" %")))

    def _show_zoom(self, eff_dpi):
        w_in, h_in = self.fig.get_size_inches()
        self.zoom_lbl.setText(
            f"{w_in:.1f} × {h_in:.1f} in   ·   {eff_dpi / 96 * 100:.0f}%")

    def render_current(self):
        """Re-render the pixmap from the current figure and show it."""
        try:
            pm = render_figure_to_pixmap(self.fig, self.RENDER_DPI)
            self.view.set_pixmap(pm)
        except Exception:
            traceback.print_exc()

    def redraw(self):
        try:
            self.draw()
        except Exception:
            self.fig.clear()
            ax = self.fig.add_subplot(111)
            ax.text(0.5, 0.5, traceback.format_exc(limit=2), ha="center",
                    va="center", fontsize=7, color="crimson", wrap=True)
            ax.axis("off")
        self.render_current()

    # back-compat name used by a couple of callers
    def rescale(self):
        self.view.apply_fit()


# ----------------------------------------------------------------------
# TAB 1 -- KDE + pie
# ----------------------------------------------------------------------

class KDETab(Tab):
    """The "KDE + pie" tab. Exposes the bandwidths, age ranges, bin colours,
    pie and label options, then calls dz_figs.make_kde_figure()."""

    def __init__(self, app):
        super().__init__(app)
        R = app.request_redraw
        c = F.KDEConfig()

        self.y0 = spin(0, 4000, c.young_range[0], 5, 0)
        self.y1 = spin(0, 4000, c.young_range[1], 5, 0)
        self.o0 = spin(0, 4000, c.old_range[0], 10, 0)
        self.o1 = spin(0, 4600, c.old_range[1], 100, 0)
        self.p0 = spin(0, 4000, c.pie_range[0], 10, 0)
        self.p1 = spin(0, 4600, c.pie_range[1], 100, 0)

        self.ybw = spin(0.5, 200, c.young_bandwidth, 1, 1, " Ma")
        self.obw = spin(0.5, 400, c.old_bandwidth, 5, 1, " Ma")
        self.ybin = spin(0.5, 200, c.young_bin_width, 1, 1, " Ma")
        self.obin = spin(0.5, 400, c.old_bin_width, 5, 1, " Ma")
        self.sbw = spin(0.5, 400, c.single_bandwidth, 1, 1, " Ma")
        self.sbw.setToolTip("Bandwidth for the single full-range panel "
                            "(used only when Split is off).")
        self.sbin = spin(0.5, 400, c.single_bin_width, 1, 1, " Ma")

        self.lockbin = QtWidgets.QCheckBox("Bin width follows bandwidth")
        self.lockbin.setChecked(True)
        self.ybw.valueChanged.connect(self._sync)
        self.obw.valueChanged.connect(self._sync)
        self.sbw.valueChanged.connect(self._sync)

        self.figw = spin(4, 30, c.fig_width, 0.5, 1, " in")
        self.wy = spin(0.5, 20, c.young_panel_width, 0.5, 1)
        self.wo = spin(0.5, 20, c.old_panel_width, 0.5, 1)
        self.wp = spin(0.5, 20, c.pie_panel_width, 0.5, 1)
        self.rh = spin(0.5, 8, c.row_height, 0.1, 2, " in")
        self.vs = spin(0, 3, c.vertical_spacing, 0.05, 2)
        self.hs = spin(0, 3, c.horizontal_spacing, 0.05, 2)

        self.lp = QtWidgets.QCheckBox("Label peak ages")
        self.lp.setChecked(c.label_peaks)
        self.pt = spin(0, 1, c.peak_threshold, 0.05, 2)
        self.pt.setToolTip("0 = label every peak.  1 = label only the tallest.")

        self.pin = spin(0, 100, c.pie_inside_threshold, 5, 1, " %")
        self.pmin = spin(0, 100, c.pie_min_percent, 1, 1, " %")
        self.split = QtWidgets.QCheckBox("Split into young + old panels")
        self.split.setChecked(c.show_split)
        self.split.setToolTip("Off = one panel spanning the whole age range "
                              "(young min to old max), using the old-panel "
                              "bandwidth and bin width.")

        self.showpie = QtWidgets.QCheckBox("Show pie column")
        self.showpie.setChecked(c.show_pie)
        self.prad = spin(0.3, 1.25, c.pie_radius, 0.05, 2)
        self.prad.setToolTip("Pie size within its cell. Labels sit just "
                             "outside the rim, so this caps around 1.25. "
                             "For bigger pies raise 'Pie col width' / "
                             "'Row height'.")
        self.pdist = spin(0.0, 0.8, c.pie_label_dist, 0.02, 2)
        self.pdist.setToolTip("Gap between the pie rim and the % labels, "
                             "which sit outside the pie with leader lines. "
                             "Larger = labels further out.")

        self.leg = QtWidgets.QCheckBox("Legend across the top")
        self.leg.setChecked(True)
        self.lncol = spin(1, 10, c.legend_ncol, 1, 0)

        self.fs_name = spin(3, 24, c.sample_name_size, 0.5, 1)
        self.fs_n = spin(3, 24, c.n_label_size, 0.5, 1)
        self.fs_tick = spin(3, 24, c.axis_tick_size, 0.5, 1)
        self.fs_axis = spin(3, 24, c.axis_label_size, 0.5, 1)
        self.fs_leg = spin(3, 24, c.legend_size, 0.5, 1)
        self.fs_peak = spin(3, 24, c.peak_label_size, 0.5, 1)
        self.fs_pin = spin(3, 24, c.pie_inside_label_size, 0.5, 1)
        self.fs_pout = spin(3, 24, c.pie_outside_label_size, 0.5, 1)

        self.add(group("Age ranges (Ma)", form([
            (None, self.split),
            ("Young panel from", self.y0), ("Young panel to", self.y1),
            ("Old panel from", self.o0), ("Old panel to", self.o1),
            ("Pie from", self.p0), ("Pie to", self.p1),
        ])))
        self.add(group("Bandwidth & bins", form([
            ("Young bandwidth", self.ybw),
            ("Old bandwidth", self.obw),
            (None, self.lockbin),
            ("Young bin width", self.ybin),
            ("Old bin width", self.obin),
            ("Single-panel bandwidth", self.sbw),
            ("Single-panel bin width", self.sbin),
        ])))
        self.add(group("Peak labels", form([
            (None, self.lp), ("Threshold", self.pt),
        ])))
        self.add(group("Pie", form([
            (None, self.showpie),
            ("Radius", self.prad),
            ("Label distance", self.pdist),
            ("Don't label below", self.pmin),
        ])))
        self.add(group("Layout", form([
            ("Figure width", self.figw),
            ("Young col width", self.wy),
            ("Old col width", self.wo),
            ("Pie col width", self.wp),
            ("Row height", self.rh),
            ("Row spacing", self.vs),
            ("Column spacing", self.hs),
            (None, self.leg), ("Legend columns", self.lncol),
        ]), collapsed=True))
        self.add(group("Font sizes", form([
            ("Sample name", self.fs_name), ("n =", self.fs_n),
            ("Ticks", self.fs_tick), ("Axis label", self.fs_axis),
            ("Legend", self.fs_leg), ("Peak label", self.fs_peak),
            ("Pie inside", self.fs_pin), ("Pie outside", self.fs_pout),
        ]), collapsed=True))
        self.done()

        for w in self.__dict__.values():
            if isinstance(w, (QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox)):
                w.valueChanged.connect(R)
            elif isinstance(w, QtWidgets.QCheckBox):
                w.stateChanged.connect(R)

    def _sync(self):
        if self.lockbin.isChecked():
            self.ybin.setValue(self.ybw.value())
            self.obin.setValue(self.obw.value())
            self.sbin.setValue(self.sbw.value())

    def config(self) -> F.KDEConfig:
        bins, labels, colors = self.app.bins.values()
        return F.KDEConfig(
            young_range=(self.y0.value(), self.y1.value()),
            old_range=(self.o0.value(), self.o1.value()),
            pie_range=(self.p0.value(), self.p1.value()),
            young_bandwidth=self.ybw.value(), old_bandwidth=self.obw.value(),
            young_bin_width=self.ybin.value(), old_bin_width=self.obin.value(),
            single_bandwidth=self.sbw.value(),
            single_bin_width=self.sbin.value(),
            age_bins=bins, bin_labels=labels, bin_colors=colors,
            fig_width=self.figw.value(),
            young_panel_width=self.wy.value(),
            old_panel_width=self.wo.value(),
            pie_panel_width=self.wp.value(),
            row_height=self.rh.value(),
            vertical_spacing=self.vs.value(),
            horizontal_spacing=self.hs.value(),
            sample_name_size=self.fs_name.value(),
            n_label_size=self.fs_n.value(),
            axis_tick_size=self.fs_tick.value(),
            axis_label_size=self.fs_axis.value(),
            legend_size=self.fs_leg.value(),
            peak_label_size=self.fs_peak.value(),
            pie_inside_label_size=self.fs_pin.value(),
            pie_outside_label_size=self.fs_pout.value(),
            legend_ncol=int(self.lncol.value()),
            label_peaks=self.lp.isChecked(),
            peak_threshold=self.pt.value(),
            pie_inside_threshold=self.pin.value(),
            pie_min_percent=self.pmin.value(),
            pie_radius=self.prad.value(),
            pie_label_dist=self.pdist.value(),
            show_pie=self.showpie.isChecked(),
            show_split=self.split.isChecked(),
            show_legend=self.leg.isChecked(),
        )

    def draw(self):
        df, samples = self.app.df, self.app.samples()
        if df is None or not samples:
            self.fig.clear()
            ax = self.fig.add_subplot(111)
            ax.text(0.5, 0.5, "Open a file and tick some samples",
                    ha="center", va="center", color="#888")
            ax.axis("off")
            return
        ages, labels = F.load_ages(df, samples, self.app.id_col(),
                                   self.app.age_col())
        F.make_kde_figure(ages, labels, self.config(), fig=self.fig)


# ----------------------------------------------------------------------
# TAB 2 -- MDS
# ----------------------------------------------------------------------

class MDSTab(Tab):
    """The MDS tab. Exposes the dissimilarity metric, scaling options and
    consistent-axis controls, then calls dz_figs.make_mds_figure(). Also
    offers the underlying dissimilarity matrix as a CSV export."""

    def __init__(self, app):
        super().__init__(app)
        R = app.request_redraw
        c = F.MDSConfig()

        self.thr = spin(0, 4000, c.age_threshold, 5, 0, " Ma")
        self.thr.setToolTip("Grains younger than this are dropped before "
                            "the dissimilarity is computed.")
        self.metric = QtWidgets.QComboBox()
        self._metric_keys = {"K–S": "ks",
                             "Cross-correlation": "cross_correlation",
                             "Likeness": "likeness"}
        self.metric.addItems(list(self._metric_keys))
        self.metric.setCurrentText("Cross-correlation")
        self.bw = spin(1, 400, c.kde_bandwidth, 5, 1, " Ma")
        self.bw.setToolTip("Only used by cross-correlation and likeness.")
        self.mtype = QtWidgets.QComboBox()
        self.mtype.addItems(["Non-metric", "Metric"])
        self.seed = spin(0, 9999, c.seed, 1, 0)
        self.ninit = spin(1, 50, c.n_init, 1, 0)

        self.fw = spin(3, 20, c.fig_width, 0.5, 1, " in")
        self.fh = spin(3, 20, c.fig_height, 0.5, 1, " in")
        self.ps = spin(10, 600, c.point_size, 10, 0)
        self.lfs = spin(4, 24, c.label_font_size, 0.5, 1)
        self.showlab = QtWidgets.QCheckBox("Sample labels")
        self.showlab.setChecked(True)
        self.showlab.setToolTip("Untick for a clean figure to label in "
                                "Illustrator.")
        self.showlines = QtWidgets.QCheckBox("Nearest-neighbour lines")
        self.showlines.setChecked(True)
        self.showtitle = QtWidgets.QCheckBox("Show title")
        self.showtitle.setChecked(True)
        self.showstress = QtWidgets.QCheckBox("Show stress caption")
        self.showstress.setChecked(True)
        self.pcol = ColorButton(c.point_color)
        self.pcol.changed.connect(R)

        # consistent axes across figures (Vermeesch statistics unchanged)
        self.norm = QtWidgets.QCheckBox("Same axes for every figure")
        self.norm.setChecked(c.normalize_scale)
        self.norm.setToolTip("Scale each map so its farthest point sits at "
                             "radius 1 and use a fixed square axis range, so "
                             "separate figures share identical axes. This only "
                             "rescales the plot; it does not change the MDS.")
        self.axlim = spin(0.5, 3.0, c.axis_limit, 0.05, 2)
        self.aspect = QtWidgets.QCheckBox("Equal aspect (square units)")
        self.aspect.setChecked(c.equal_aspect)

        # editable text (blank = automatic default)
        self.e_title = lineedit("automatic")
        self.e_xlab = lineedit("MDS dimension 1")
        self.e_ylab = lineedit("MDS dimension 2")

        self.matbtn = QtWidgets.QPushButton("Export dissimilarity matrix (CSV)")
        self.matbtn.clicked.connect(self.export_matrix)

        self.add(group("Dissimilarity", form([
            ("Drop grains younger than", self.thr),
            ("Metric", self.metric),
            ("KDE bandwidth", self.bw),
        ])))
        self.add(group("Scaling", form([
            ("MDS type", self.mtype),
            ("Seed", self.seed),
            ("Restarts", self.ninit),
        ])))
        self.add(group("Style", form([
            ("Fig width", self.fw), ("Fig height", self.fh),
            ("Point size", self.ps), ("Point colour", self.pcol),
            ("Label font", self.lfs),
            (None, self.showlab), (None, self.showlines),
            (None, self.showtitle), (None, self.showstress),
        ])))
        self.add(group("Consistent axes", form([
            (None, self.norm),
            ("Axis limit", self.axlim),
            (None, self.aspect),
        ])))
        self.add(group("Text & labels", form([
            ("Title", self.e_title),
            ("X-axis", self.e_xlab),
            ("Y-axis", self.e_ylab),
        ]), collapsed=True))
        self.add(group("Export", form([(None, self.matbtn)])))
        self.done()

        for w in (self.thr, self.bw, self.seed, self.ninit, self.fw, self.fh,
                  self.ps, self.lfs, self.axlim):
            w.valueChanged.connect(R)
        for w in (self.metric, self.mtype):
            w.currentIndexChanged.connect(R)
        for w in (self.showlab, self.showlines, self.showtitle, self.showstress,
                  self.norm, self.aspect):
            w.stateChanged.connect(R)
        for w in (self.e_title, self.e_xlab, self.e_ylab):
            w.textChanged.connect(R)

        self._D = None
        self._labels = []

    def config(self) -> F.MDSConfig:
        return F.MDSConfig(
            age_threshold=self.thr.value(),
            metric=self._metric_keys[self.metric.currentText()],
            kde_bandwidth=self.bw.value(),
            mds_type=self.mtype.currentText().lower(),
            seed=int(self.seed.value()),
            n_init=int(self.ninit.value()),
            fig_width=self.fw.value(), fig_height=self.fh.value(),
            point_size=self.ps.value(), point_color=self.pcol.color(),
            label_font_size=self.lfs.value(),
            show_labels=self.showlab.isChecked(),
            show_neighbor_lines=self.showlines.isChecked(),
            show_title=self.showtitle.isChecked(),
            show_stress=self.showstress.isChecked(),
            normalize_scale=self.norm.isChecked(),
            axis_limit=self.axlim.value(),
            equal_aspect=self.aspect.isChecked(),
            title=self.e_title.text(),
            xlabel=self.e_xlab.text(),
            ylabel=self.e_ylab.text(),
        )

    def draw(self):
        df, samples = self.app.df, self.app.samples()
        self.fig.clear()
        if df is None or len(samples) < 3:
            ax = self.fig.add_subplot(111)
            ax.text(0.5, 0.5, "MDS needs at least three samples",
                    ha="center", va="center", color="#888")
            ax.axis("off")
            return
        cfg = self.config()
        ages, labels = F.load_ages(df, samples, self.app.id_col(),
                                   self.app.age_col(),
                                   threshold=cfg.age_threshold)
        if len(labels) < 3:
            ax = self.fig.add_subplot(111)
            ax.text(0.5, 0.5,
                    f"Only {len(labels)} samples survive the >{cfg.age_threshold:g} Ma\n"
                    "filter with 3+ grains. Lower the threshold.",
                    ha="center", va="center", color="crimson")
            ax.axis("off")
            return
        _, self._D, self._labels, stress = F.make_mds_figure(
            ages, labels, cfg, fig=self.fig)
        self.app.status(f"MDS stress = {stress:.3f}  "
                        f"({'excellent' if stress < .10 else 'acceptable' if stress < .20 else 'UNRELIABLE — do not over-read this map'})")

    def export_matrix(self):
        if self._D is None:
            return
        p, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save matrix", "dissimilarity.csv", "CSV (*.csv)")
        if p:
            pd.DataFrame(self._D, index=self._labels,
                         columns=self._labels).to_csv(p)
            self.app.status(f"Saved {p}")


# ----------------------------------------------------------------------
#  TAB 3  --  MDA
# ----------------------------------------------------------------------

class MDATab(Tab):
    """The MDA tab. Lets the user pick the error column and sigma convention,
    computes YSG / YGC1s / YGC2s for every ticked sample into a table, and
    draws the ranked-grain plot for one chosen sample via dz_mda."""

    def __init__(self, app):
        super().__init__(app)
        R = app.request_redraw

        # error column + sigma (MDA is the only thing that needs the error)
        self.errcol = QtWidgets.QComboBox()
        self.errcol.currentIndexChanged.connect(R)
        self.errsig = QtWidgets.QComboBox()
        self.errsig.addItems(["2σ", "1σ"])
        self.errsig.currentIndexChanged.connect(R)
        self.repsig = QtWidgets.QComboBox()
        self.repsig.addItems(["2σ", "1σ"])
        self.repsig.currentIndexChanged.connect(R)
        self._cols_seen = None

        # exclude youngest applies to every method
        self.excl = QtWidgets.QCheckBox("Exclude the single youngest grain")
        self.excl.setToolTip("Use when the youngest grain looks like an "
                             "outlier. It is removed before every method runs.")

        # which sample to draw, which method, plot extent.  Methods are the
        # standard set with fixed definitions (YGC1σ = >=2 grains at 1σ,
        # YGC2σ = >=3 grains at 2σ, both contiguous, per detritalPy).
        self.plotsample = QtWidgets.QComboBox()
        self.plotsample.currentIndexChanged.connect(R)
        self._method_keys = {"YSG": "YSG", "YGC1σ": "YC1s", "YGC2σ": "YC2s"}
        self.method = QtWidgets.QComboBox()
        self.method.addItems(list(self._method_keys))
        self.method.setCurrentText("YGC2σ")
        self.method.currentIndexChanged.connect(R)
        self.nshow = spin(5, 200, 25, 5, 0, " grains")
        self.agemax = spin(0, 4000, 300, 50, 0, " Ma")
        self.agemax_on = QtWidgets.QCheckBox("Only show grains younger than")
        self.agemax_on.setChecked(True)
        self.fw = spin(3, 20, 6.5, 0.5, 1, " in")
        self.fh = spin(2, 16, 4.5, 0.5, 1, " in")

        # editable text + legend
        self.showleg = QtWidgets.QCheckBox("Show legend")
        self.showleg.setChecked(True)
        self.e_title = lineedit("automatic")
        self.e_xlab = lineedit("Grain rank (youngest → oldest)")
        self.e_ylab = lineedit("Age (Ma)")
        self.e_in = lineedit("Grains in cluster")
        self.e_out = lineedit("Grains not in cluster")

        self.add(group("Uncertainty", form([
            ("Error column", self.errcol),
            ("Error is", self.errsig),
            ("Report at", self.repsig),
            (None, self.excl),
        ])))
        self.add(group("Rank plot", form([
            ("Sample to plot", self.plotsample),
            ("Method", self.method),
            ("Show youngest", self.nshow),
            (None, self.agemax_on),
            ("Age cut-off", self.agemax),
            ("Fig width", self.fw), ("Fig height", self.fh),
        ])))
        self.add(group("Text & labels", form([
            ("Title", self.e_title),
            ("X-axis", self.e_xlab),
            ("Y-axis", self.e_ylab),
            (None, self.showleg),
            ("Legend: in", self.e_in),
            ("Legend: out", self.e_out),
        ]), collapsed=True))

        # results table for ALL selected samples
        self.table = QtWidgets.QTableWidget()
        self.table.setMinimumHeight(160)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.add(QtWidgets.QLabel("Results (all ticked samples):"))
        self.add(self.table)
        self.csvbtn = QtWidgets.QPushButton("Export results table (CSV)")
        self.csvbtn.clicked.connect(self.export_csv)
        self.add(self.csvbtn)
        self.done()

        for w in (self.nshow, self.agemax, self.fw, self.fh):
            w.valueChanged.connect(R)
        for w in (self.excl, self.agemax_on, self.showleg):
            w.stateChanged.connect(R)
        for w in (self.e_title, self.e_xlab, self.e_ylab, self.e_in, self.e_out):
            w.textChanged.connect(R)

        self._df_cache = None

    def _sync_columns(self):
        """Populate the error-column dropdown from the current sheet."""
        df = self.app.df
        if df is None:
            return
        cols = [str(c) for c in df.columns]
        if cols == self._cols_seen:
            return
        self._cols_seen = cols
        self.errcol.blockSignals(True)
        self.errcol.clear()
        self.errcol.addItems(cols)
        # guess a column with "err" in the name
        guess = next((c for c in cols if "err" in c.lower()), None)
        if guess:
            self.errcol.setCurrentText(guess)
        self.errcol.blockSignals(False)

    def _method_key(self):
        return self._method_keys[self.method.currentText()]

    def _params(self):
        # standard, fixed definitions -- no user-editable thresholds
        return MDA.MDAParams(
            n_min=3, k_sigma=2.0, contiguous=True,
            exclude_youngest=self.excl.isChecked(),
            report_sigma=1 if self.repsig.currentText() == "1σ" else 2,
        )

    def _err_sigma(self):
        return 1 if self.errsig.currentText() == "1σ" else 2

    def _refresh_sample_list(self, samples):
        cur = self.plotsample.currentText()
        self.plotsample.blockSignals(True)
        self.plotsample.clear()
        self.plotsample.addItems(samples)
        if cur in samples:
            self.plotsample.setCurrentText(cur)
        self.plotsample.blockSignals(False)

    def _fill_table(self, rows):
        if not rows:
            self.table.clear()
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            return
        cols = list(rows[0].keys())
        self.table.setColumnCount(len(cols))
        self.table.setHorizontalHeaderLabels(cols)
        self.table.setRowCount(len(rows))
        for i, r in enumerate(rows):
            for j, c in enumerate(cols):
                v = r[c]
                txt = "" if v is None or (isinstance(v, float) and np.isnan(v)) \
                    else str(v)
                self.table.setItem(i, j, QtWidgets.QTableWidgetItem(txt))
        self.table.resizeColumnsToContents()
        self._rows = rows

    def draw(self):
        df = self.app.df
        samples = self.app.samples()
        self.fig.clear()
        if df is None or not samples:
            ax = self.fig.add_subplot(111)
            ax.text(0.5, 0.5, "Open a file and tick some samples",
                    ha="center", va="center", color="#888")
            ax.axis("off")
            return

        self._sync_columns()
        self._refresh_sample_list(samples)
        params = self._params()
        errcol = self.errcol.currentText() or None
        esig = self._err_sigma()

        # build the results table for all ticked samples
        rows = []
        for s in samples:
            a, e = MDA.sample_ages_errors(df, self.app.id_col(),
                                          self.app.age_col(), errcol, s, esig)
            if a.size == 0:
                rows.append({"Sample": s, "n": 0})
                continue
            rows.append(MDA.table_row(s, a, e, params))
        self._fill_table(rows)

        # rank plot for the chosen sample
        ps = self.plotsample.currentText() or samples[0]
        a, e = MDA.sample_ages_errors(df, self.app.id_col(),
                                      self.app.age_col(), errcol, ps, esig)
        self.fig.set_size_inches(self.fw.value(), self.fh.value())
        amax = self.agemax.value() if self.agemax_on.isChecked() else None
        MDA.make_rank_plot(self.fig, ps, a, e, params,
                           method=self._method_key(),
                           n_show=int(self.nshow.value()), age_max=amax,
                           title=self.e_title.text(),
                           xlabel=self.e_xlab.text(),
                           ylabel=self.e_ylab.text(),
                           show_legend=self.showleg.isChecked(),
                           label_in=self.e_in.text() or "Grains in cluster",
                           label_out=self.e_out.text() or "Grains not in cluster")

        if errcol is None:
            self.app.status("No error column selected — MDA needs one")
        else:
            self.app.status(f"MDA computed for {len(samples)} samples "
                            f"(error column '{errcol}', {esig}σ)")

    def export_csv(self):
        rows = getattr(self, "_rows", None)
        if not rows:
            return
        p, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save MDA table", "mda_results.csv", "CSV (*.csv)")
        if p:
            pd.DataFrame(rows).to_csv(p, index=False)
            self.app.status(f"Saved {Path(p).name}")


class ExportDialog(QtWidgets.QDialog):
    """
    One dialog for all exporting.  Pick the format, the printed size and (for
    PNG) the resolution.  Size presets are the usual AGU/GSA column widths.
    The preview zoom has nothing to do with any of this.
    """
    PRESETS = {
        "Current figure size": None,
        "Single column (3.5 in)": 3.5,
        "Two column (7.5 in)": 7.5,
        "Full page (11 in)": 11.0,
    }

    def __init__(self, parent, w_in, h_in):
        super().__init__(parent)
        self.setWindowTitle("Export figure")
        self.setMinimumWidth(340)
        self._w0, self._h0 = w_in, h_in
        self._ratio = h_in / w_in if w_in else 1.0

        self.fmt = QtWidgets.QComboBox()
        self.fmt.addItems(["PDF  (vector, editable in Illustrator)",
                           "SVG  (vector, editable in Illustrator)",
                           "PNG  (raster image)"])
        self.fmt.currentIndexChanged.connect(self._fmt_changed)

        self.preset = QtWidgets.QComboBox()
        self.preset.addItems(list(self.PRESETS))
        self.preset.currentIndexChanged.connect(self._preset_changed)

        self.w = spin(1, 60, w_in, 0.25, 2, " in")
        self.h = spin(1, 80, h_in, 0.25, 2, " in")
        self.lock = QtWidgets.QCheckBox("Keep aspect ratio")
        self.lock.setChecked(True)
        self.w.valueChanged.connect(self._w_changed)

        self.dpi = spin(72, 1200, 300, 50, 0, " dpi")
        self.tight = QtWidgets.QCheckBox("Trim white margins")
        self.tight.setChecked(True)

        self.note = QtWidgets.QLabel()
        self.note.setWordWrap(True)
        self.note.setStyleSheet("color:#666;")

        btns = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Save | QtWidgets.QDialogButtonBox.Cancel)
        btns.button(QtWidgets.QDialogButtonBox.Save).setText("Export…")
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)

        v = QtWidgets.QVBoxLayout(self)
        v.addWidget(form([
            ("Format", self.fmt),
            ("Size preset", self.preset),
            ("Width", self.w),
            ("Height", self.h),
            (None, self.lock),
            ("Resolution", self.dpi),
            (None, self.tight),
        ]))
        v.addWidget(self.note)
        v.addWidget(btns)
        self._fmt_changed()

    def fmt_ext(self):
        return ["pdf", "svg", "png"][self.fmt.currentIndex()]

    def _fmt_changed(self):
        is_png = self.fmt_ext() == "png"
        self.dpi.setEnabled(is_png)
        self.note.setText(
            "PNG is a raster image. 300 dpi for review, 600 for print."
            if is_png else
            "Vector format: text stays live and editable in Illustrator; "
            "resolution does not apply.")

    def _preset_changed(self):
        target = self.PRESETS[self.preset.currentText()]
        if target is None:
            self.w.setValue(self._w0)
            self.h.setValue(self._h0)
        else:
            self.w.setValue(target)

    def _w_changed(self, val):
        if self.lock.isChecked():
            self.h.blockSignals(True)
            self.h.setValue(val * self._ratio)
            self.h.blockSignals(False)

    def values(self):
        return (self.fmt_ext(), self.w.value(), self.h.value(),
                int(self.dpi.value()), self.tight.isChecked())


# ----------------------------------------------------------------------
# main window
# ----------------------------------------------------------------------

class MainWindow(QtWidgets.QMainWindow):
    """The application window: the menu bar, the dataset dock (file, sheet,
    sample and age columns), the three tabs, image export, and saving/loading
    .dzp project files."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{APP} {VERSION}")
        self.resize(1500, 950)
        self.df = None
        self.path = None

        self.file_lbl = QtWidgets.QLabel("no file loaded")
        self.file_lbl.setWordWrap(True)
        self.sheet = QtWidgets.QComboBox()
        self.c_id = QtWidgets.QComboBox()
        self.c_age = QtWidgets.QComboBox()
        ob = QtWidgets.QPushButton("Open Excel file…")
        ob.clicked.connect(self.open_file)
        ab = QtWidgets.QPushButton("Apply columns")
        ab.clicked.connect(self.rebuild)
        self.sheet.currentIndexChanged.connect(self._sheet_changed)

        self.slist = SampleList()
        self.slist.changed.connect(self.request_redraw)
        self.bins = BinEditor()
        self.bins.changed.connect(self.request_redraw)

        left = QtWidgets.QWidget()
        lv = QtWidgets.QVBoxLayout(left)
        lv.setContentsMargins(4, 4, 4, 4)
        lv.addWidget(group("Data", form([
            (None, ob), (None, self.file_lbl),
            ("Sheet", self.sheet),
            ("Sample column", self.c_id),
            ("Age column", self.c_age),
            (None, ab),
        ])))
        lv.addWidget(group("Samples", self.slist), 2)
        lv.addWidget(group("Provenance bins", self.bins), 3)

        self.dock = QtWidgets.QDockWidget("Dataset")
        self.dock.setWidget(left)
        self.dock.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable |
                              QtWidgets.QDockWidget.DockWidgetClosable)
        self.dock.setMinimumWidth(360)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock)

        self.tabs = QtWidgets.QTabWidget()
        # macOS centres the tab bar by default; Windows left-aligns it. Pin it
        # to the left with a stylesheet so both platforms look the same.
        self.tabs.setDocumentMode(True)
        self.tabs.tabBar().setExpanding(False)
        self.tabs.setStyleSheet(
            "QTabWidget::tab-bar { alignment: left; }")
        self.kde = KDETab(self)
        self.mds = MDSTab(self)
        self.mda = MDATab(self)
        self.tabs.addTab(self.kde, "KDE + pie")
        self.tabs.addTab(self.mds, "MDS")
        self.tabs.addTab(self.mda, "MDA")
        self.tabs.currentChanged.connect(self.request_redraw)
        self.setCentralWidget(self.tabs)

        tb = self.addToolBar("main")
        tb.setMovable(False)
        exp = QtGui.QAction("Export…", self)
        exp.setShortcut("Ctrl+E")
        exp.triggered.connect(self.export)
        tb.addAction(exp)
        tb.addSeparator()
        for t, fn in (("Save project", self.save_project),
                      ("Load project", self.load_project)):
            a = QtGui.QAction(t, self)
            a.triggered.connect(fn)
            tb.addAction(a)
        tb.addSeparator()
        hide = QtGui.QAction("Hide panels (F11)", self)
        hide.setShortcut("F11")
        hide.setCheckable(True)
        hide.setToolTip("Collapse the side panels to give the figure the "
                        "whole window. Press F11 again to bring them back.")
        hide.toggled.connect(self._hide_panels)
        tb.addAction(hide)
        tb.addSeparator()
        self.live = QtWidgets.QCheckBox("Live preview")
        self.live.setChecked(True)
        tb.addWidget(self.live)
        a = QtGui.QAction("Redraw (F5)", self)
        a.setShortcut("F5")
        a.triggered.connect(self._redraw_now)
        tb.addAction(a)

        self.statusBar().showMessage("Open an Excel file to begin")
        self._t = QtCore.QTimer(self)
        self._t.setSingleShot(True)
        self._t.setInterval(300)
        self._t.timeout.connect(self._redraw_now)

    def _hide_panels(self, on):
        self.dock.setVisible(not on)
        for t in (self.kde, self.mds, self.mda):
            t.settings_scroll.setVisible(not on)
            if t.view.fit_mode:
                QtCore.QTimer.singleShot(0, t.view.apply_fit)

    # -- plumbing
    def status(self, m):
        self.statusBar().showMessage(m, 8000)

    def id_col(self):
        return self.c_id.currentText()

    def age_col(self):
        return self.c_age.currentText()

    def samples(self):
        return self.slist.selected() if self.df is not None else []

    def request_redraw(self, *_):
        if not self.live.isChecked():
            self.statusBar().showMessage("Live preview off — press F5 to update")
            return
        self._t.start()

    def _redraw_now(self):
        w = self.tabs.currentWidget()
        t0 = time.perf_counter()
        w.redraw()
        self.statusBar().showMessage(
            f"redrawn in {(time.perf_counter()-t0)*1000:.0f} ms", 3000)

    # -- data
    def open_file(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open DZ data", "", "Data (*.xlsx *.xls *.csv)")
        if not p:
            return
        self.path = p
        self.file_lbl.setText(Path(p).name)
        self.sheet.blockSignals(True)
        self.sheet.clear()
        names = F.sheet_names(p)
        if names:
            self.sheet.addItems(names)
            pick = next((n for n in names if "zrupb" in n.lower()
                         .replace(" ", "").replace("_", "")), names[0])
            self.sheet.setCurrentText(pick)
        self.sheet.blockSignals(False)
        self._sheet_changed()

    def _sheet_changed(self, *_):
        if not self.path:
            return
        try:
            self.df = F.read_sheet(self.path, self.sheet.currentText() or 0)
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Read error", str(e))
            return
        cols = [str(c) for c in self.df.columns]
        for box, want in ((self.c_id, "Sample_ID"), (self.c_age, "BestAge")):
            box.blockSignals(True)
            box.clear()
            box.addItems(cols)
            if want in cols:
                box.setCurrentText(want)
            box.blockSignals(False)
        self.rebuild()

    def rebuild(self):
        if self.df is None:
            return
        try:
            names = F.all_samples(self.df, self.id_col())
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Column error", str(e))
            return
        counts = {s: int((self.df[self.id_col()].astype(str) == s).sum())
                  for s in names}
        self.slist.populate(names, counts)
        self.status(f"{len(names)} samples, {len(self.df)} analyses")
        self.request_redraw()

    # -- export
    def export(self):
        tab = self.tabs.currentWidget()

        # ALWAYS rebuild the figure from the current ticked samples first, so
        # the export matches what is selected even if Live preview is off or a
        # redraw is still pending.  This is the fix for "it exported all
        # samples": tab.fig could otherwise be stale.
        tab.redraw()
        fig = tab.fig
        w_in, h_in = fig.get_size_inches()

        dlg = ExportDialog(self, w_in, h_in)
        if dlg.exec() != QtWidgets.QDialog.Accepted:
            return
        fmt, out_w, out_h, dpi, keep = dlg.values()

        p, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save figure as", f"figure.{fmt}",
            f"{fmt.upper()} (*.{fmt})")
        if not p:
            return
        if not p.lower().endswith("." + fmt):
            p += "." + fmt

        old_dpi = fig.get_dpi()
        try:
            fig.set_size_inches(out_w, out_h)
            kw = dict(facecolor="white")

            if keep:
                # Compute the tight bounding box ourselves using an Agg
                # renderer (the same path PNG uses and which works in the
                # frozen app), then hand savefig an explicit Bbox.  Relying on
                # savefig's own bbox_inches="tight" can yield a BLANK page in
                # a PyInstaller build, because the crop is computed with the
                # vector backend's renderer, which can return an empty box.
                from matplotlib.backends.backend_agg import FigureCanvasAgg
                agg = FigureCanvasAgg(fig)
                agg.draw()
                bbox = fig.get_tightbbox(agg.get_renderer()).padded(0.1)
                kw["bbox_inches"] = bbox

            if fmt == "png":
                kw["dpi"] = dpi

            try:
                fig.savefig(p, **kw)
            except Exception:
                # last-ditch fallback: save without the tight crop rather than
                # produce nothing
                kw.pop("bbox_inches", None)
                fig.savefig(p, **kw)
        finally:
            fig.set_size_inches(w_in, h_in)
            fig.set_dpi(old_dpi)
            tab.render_current()

        n = len(self.samples())
        note = "text stays editable in Illustrator" if fmt in ("pdf", "svg") \
            else f"{dpi} dpi"
        self.status(f"Saved {Path(p).name}  —  {n} samples, "
                    f"{out_w:.1f} × {out_h:.1f} in, {note}")

    # -- project
    def _state(self, obj):
        out = {}
        for k, w in vars(obj).items():
            if isinstance(w, (QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox)):
                out[k] = w.value()
            elif isinstance(w, QtWidgets.QCheckBox):
                out[k] = w.isChecked()
            elif isinstance(w, QtWidgets.QComboBox):
                out[k] = w.currentText()
            elif isinstance(w, QtWidgets.QLineEdit):
                out[k] = w.text()
            elif isinstance(w, ColorButton):
                out[k] = w.color()
        return out

    def _restore(self, obj, st):
        for k, v in st.items():
            w = getattr(obj, k, None)
            if isinstance(w, (QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox)):
                w.setValue(v)
            elif isinstance(w, QtWidgets.QCheckBox):
                w.setChecked(bool(v))
            elif isinstance(w, QtWidgets.QComboBox):
                w.setCurrentText(str(v))
            elif isinstance(w, QtWidgets.QLineEdit):
                w.setText(str(v))

    def save_project(self):
        p, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save project", "figure.dzp", "DZ project (*.dzp)")
        if not p:
            return
        if not p.lower().endswith(".dzp"):
            p += ".dzp"
        b, l, c = self.bins.values()
        state = {
            "version": VERSION,
            "path": self.path, "sheet": self.sheet.currentText(),
            "id_col": self.id_col(), "age_col": self.age_col(),
            "err_col": self.mda.errcol.currentText(),
            "order": [self.slist.list.item(i).data(Qt.UserRole)
                      for i in range(self.slist.list.count())],
            "checked": self.samples(),
            "bins": b, "bin_labels": l, "bin_colors": c,
            "kde": self._state(self.kde), "mds": self._state(self.mds),
            "mda": self._state(self.mda),
        }
        Path(p).write_text(json.dumps(state, indent=2))
        self.status("Project saved")

    def load_project(self):
        p, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Load project", "", "DZ project (*.dzp)")
        if not p:
            return
        st = json.loads(Path(p).read_text())

        if st.get("path") and Path(st["path"]).exists():
            self.path = st["path"]
            self.file_lbl.setText(Path(self.path).name)
            self.sheet.blockSignals(True)
            self.sheet.clear()
            names = F.sheet_names(self.path)
            if names:
                self.sheet.addItems(names)
                self.sheet.setCurrentText(st.get("sheet", names[0]))
            self.sheet.blockSignals(False)
            self.df = F.read_sheet(self.path, st.get("sheet") or 0)
            cols = [str(c) for c in self.df.columns]
            for box, key in ((self.c_id, "id_col"), (self.c_age, "age_col")):
                box.blockSignals(True)
                box.clear()
                box.addItems(cols)
                box.setCurrentText(st.get(key, cols[0]))
                box.blockSignals(False)
            self.rebuild()
            self.slist.set_state(st.get("order", []), set(st.get("checked", [])))
        else:
            QtWidgets.QMessageBox.warning(
                self, "Data not found",
                "This project links to a data file that isn't on this "
                "computer. Open the spreadsheet first, then load the project.")
        if st.get("bins"):
            self.bins.set_values([tuple(b) for b in st["bins"]],
                                 st["bin_labels"], st["bin_colors"])
        self._restore(self.kde, st.get("kde", {}))
        self._restore(self.mds, st.get("mds", {}))
        # tolerate old projects that stored raw metric/type keys
        mstate = st.get("mds", {})
        legacy = {"ks": "K–S", "cross_correlation": "Cross-correlation",
                  "likeness": "Likeness", "non-metric": "Non-metric",
                  "metric": "Metric"}
        for attr in ("metric", "mtype"):
            val = mstate.get(attr)
            if val in legacy:
                getattr(self.mds, attr).setCurrentText(legacy[val])
        self._restore(self.mda, st.get("mda", {}))
        if st.get("err_col"):
            self.mda._sync_columns()
            self.mda.errcol.setCurrentText(st["err_col"])
        if "pcol" in st.get("mds", {}):
            self.mds.pcol._c = st["mds"]["pcol"]
            self.mds.pcol._paint()
        self.status("Project loaded")
        self.request_redraw()


def main():
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName(APP)
    MainWindow().show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
