from __future__ import annotations
import datetime
import json
import math
import os
import re
import sys
import textwrap
import time
from pathlib import Path
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from PySide6.QtCore import QSettings, QSize, Qt, QThread, QTimer, QUrl, Signal
from PySide6.QtGui import QAction, QCursor, QDesktopServices, QGuiApplication, QFont, QFontDatabase, QIcon, QKeySequence, QPixmap, QShortcut
from PySide6.QtWidgets import (QApplication, QButtonGroup, QCheckBox, QComboBox, QDialog, QDoubleSpinBox,
                               QFileDialog, QFrame, QGridLayout, QHBoxLayout, QHeaderView, QLabel, QLineEdit,
                               QMainWindow, QMenu, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QRadioButton,
                               QScrollArea, QSizePolicy, QSlider, QSpinBox, QSplitter, QStackedWidget, QTableView,
                               QTabWidget, QTextEdit, QVBoxLayout, QWidget)
import core_calcolo as engine
import gene_plots as gplot
import pca_coloring as pcol
import report_pdf
from help_content import HELP_CONTENT
from qt_models import DataFrameModel
from qt_plots import HAS_MPLCURSORS, PlotPanel, confidence_ellipse
from qt_theme import (PALETTE, Card, Chip, ColorSwatch, DropZone, NavButton, RefreshCombo, Toast, apply_theme,
                      fade_in, field, open_file_dialog, repolish, MONO)

def _fmt_sci(value) -> str:
    if value is None or pd.isna(value):
        return ""
    try:
        return f"{float(value):.2e}"
    except (TypeError, ValueError):
        return str(value)

def _round_or_blank(v, n=3):
    try:
        return "" if pd.isna(v) else round(float(v), n)
    except (TypeError, ValueError):
        return ""

class Worker(QThread):
    log = Signal(str)
    done = Signal(object)
    failed = Signal(str, object)

    def __init__(self, fn, title: str, parent=None):
        super().__init__(parent)
        self._fn, self._title = fn, title

    def run(self):
        try:
            self.done.emit(self._fn(self.log.emit))
        except Exception as exc:  
            self.failed.emit(self._title, exc)

class ResultPage(QWidget):

    def __init__(self, content: QWidget, empty_text: str):
        super().__init__()
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.stack = QStackedWidget()
        empty = QWidget()
        el = QVBoxLayout(empty)
        el.addStretch()
        big = QLabel("◔")
        big.setAlignment(Qt.AlignCenter)
        big.setStyleSheet(f"font-size: 46pt; color: {PALETTE['border_strong']};")
        msg = QLabel(empty_text)
        msg.setAlignment(Qt.AlignCenter)
        msg.setObjectName("muted")
        msg.setWordWrap(True)
        el.addWidget(big)
        el.addWidget(msg)
        el.addStretch()
        self.stack.addWidget(empty)
        self.stack.addWidget(content)
        lay.addWidget(self.stack)

    def set_ready(self, ready: bool):
        self.stack.setCurrentIndex(1 if ready else 0)

def _scroll(widget: QWidget, width: int | None = None) -> QScrollArea:
    sa = QScrollArea()
    sa.setWidgetResizable(True)
    sa.setFrameShape(QFrame.NoFrame)
    sa.setWidget(widget)
    sa.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
    if width:
        sa.setFixedWidth(width)
    return sa

class DEAWindow(QMainWindow):
    LABEL_MODE_DISPLAY = ["No labels", "All significant genes", "Top N by padj", "Top N by |log2FC|",
                          "Specific genes (list)"]
    LABEL_MODE_MAP = {"No labels": "none", "All significant genes": "all_sig", "Top N by padj": "top_padj",
                      "Top N by |log2FC|": "top_lfc", "Specific genes (list)": "custom"}
    ORA_DIRECTION_DISPLAY = ["All DEGs (both directions)",
                             "Only genes higher in the 1st group (log2FC > 0)",
                             "Only genes higher in the 2nd/reference group (log2FC < 0)"]
    ORA_DIRECTION_MAP = {ORA_DIRECTION_DISPLAY[0]: "all", ORA_DIRECTION_DISPLAY[1]: "up",
                         ORA_DIRECTION_DISPLAY[2]: "down"}
    SETTINGS_FORMAT = "de-gui-settings"      
    SETTINGS_VERSION = 1

    PAGES = [
        ("setup", "1", "Setup & Parameters", "WORKFLOW"),
        ("pca", "2", "PCA", "RESULTS"),
        ("volcano", "3", "Volcano plot", None),
        ("table", "4", "Gene table", None),
        ("gsea", "5", "GSEA / Pathways", None),
        ("ora", "6", "ORA / Over-rep.", None),
        ("genes", "7", "Gene counts", None),
        ("log", "≡", "Log", "SYSTEM"),]

    def __init__(self):
        super().__init__()
        self.setWindowTitle("DEA Explorer")
        self.resize(1360, 900)
        self.setMinimumSize(1100, 720)

        self.results: engine.RunResults | None = None
        self.run_params: dict | None = None
        self._pending_params: dict | None = None
        self._pending_group_restore = None
        self.available_groups: list[str] = []
        self.active_tag = ""
        self.active_key = "setup"
        self._meta_cache = (None, None, None)
        self._contrast_combos: list[QComboBox] = []
        self._direction_labels: dict[str, QLabel] = {}
        self._result_pages: list[ResultPage] = []
        self._workers: list[Worker] = []
        self._pca_colorbar = self._gene_colorbar = self._ora_colorbar = None
        self._gene_plot_frame = None
        self._run_start: float | None = None

        self._build_ui()
        self._build_shortcuts()
        geo = QSettings("DEAExplorer", "Main").value("geometry")
        if geo is not None:
            self.restoreGeometry(geo)
        self.toast = Toast(self.centralWidget())
        self._validate_setup()
        self._update_contrast_preview()
        self.statusBar().showMessage("Ready.")
        QTimer.singleShot(400, self._check_environment)

    def _build_ui(self):
        root = QWidget()
        root.setObjectName("root")
        self.setCentralWidget(root)
        outer = QVBoxLayout(root)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        outer.addWidget(self._build_header())

        body = QHBoxLayout()
        body.setSpacing(0)
        body.addWidget(self._build_sidebar())
        self.pages = QStackedWidget()
        self.pages.setObjectName("page")
        body.addWidget(self.pages, 1)
        outer.addLayout(body, 1)

        self._page_widgets: dict[str, QWidget] = {}
        builders = {"setup": self._build_setup_page, "pca": self._build_pca_page, "volcano": self._build_volcano_page,
                    "table": self._build_table_page, "gsea": self._build_gsea_page, "ora": self._build_ora_page,
                    "genes": self._build_genes_page, "log": self._build_log_page}
        for key, *_ in self.PAGES:
            w = builders[key]()
            self._page_widgets[key] = w
            self.pages.addWidget(w)

        self.run_elapsed = QLabel("")
        self.statusBar().addPermanentWidget(self.run_elapsed)
        self.statusBar().setSizeGripEnabled(False)

        self.metadata_dz.pathChanged.connect(lambda *_: (self._validate_setup(), self._refresh_pca_color_columns(), self._refresh_gene_plot_columns()))
        self.counts_dz.pathChanged.connect(lambda *_: self._validate_setup())
        self._switch_page("setup", animate=False)

    def _build_header(self) -> QFrame:
        h = QFrame()
        h.setObjectName("header")
        h.setFixedHeight(66)
        lay = QHBoxLayout(h)
        lay.setContentsMargins(22, 0, 22, 0)
        lay.setSpacing(12)
        script_dir = Path(__file__).resolve().parent
        logo = script_dir / "GUI_logo.png"
        title = QLabel("DEA Explorer")
        title.setObjectName("appTitle")
        badge = QLabel("BETA")
        badge.setObjectName("badge")
        sub = QLabel("Differential Expression & Gene Set Enrichment Analysis")
        sub.setObjectName("appSubtitle")
        lay.addWidget(title)
        lay.addWidget(badge, alignment=Qt.AlignVCenter)
        lay.addSpacing(18)
        lay.addWidget(sub)
        lay.addStretch()
        help_btn = QPushButton("?  Help")
        help_btn.setObjectName("headerBtn")
        help_btn.setCursor(Qt.PointingHandCursor)
        help_btn.setToolTip("Explanation of the current page (F1).")
        help_btn.clicked.connect(lambda: self._show_help(self.active_key))
        lay.addWidget(help_btn)
        if logo.exists():
            pm = QPixmap(str(logo))
            if not pm.isNull():
                lbl = QLabel()
                lbl.setPixmap(pm.scaledToHeight(44, Qt.SmoothTransformation))
                lay.addSpacing(10)
                lay.addWidget(lbl)
        return h

    def _build_sidebar(self) -> QFrame:
        sb = QFrame()
        sb.setObjectName("sidebar")
        sb.setFixedWidth(232)
        lay = QVBoxLayout(sb)
        lay.setContentsMargins(0, 6, 0, 8)
        lay.setSpacing(0)
        self._nav_group = QButtonGroup(self)
        self._nav_group.setExclusive(True)
        self._nav_buttons: dict[str, NavButton] = {}
        for key, badge, title, section in self.PAGES:
            if section:
                s = QLabel(section)
                s.setObjectName("sideSection")
                lay.addWidget(s)
            b = NavButton(badge, title)
            b.clicked.connect(lambda _=False, k=key: self._switch_page(k))
            self._nav_group.addButton(b)
            self._nav_buttons[key] = b
            lay.addWidget(b)
        lay.addStretch()
        foot = QLabel("Ctrl+R run  ·  F1 help\nCtrl+1…7 switch page\nCtrl+H heatmap  ·  Ctrl+E export")
        foot.setObjectName("sideFoot")
        lay.addWidget(foot)
        return sb

    def _build_shortcuts(self):
        QShortcut(QKeySequence("F1"), self, activated=lambda: self._show_help(self.active_key))
        QShortcut(QKeySequence("Ctrl+R"), self, activated=lambda: self.run_button.isEnabled() and self._run_analysis())
        for i, (key, *_rest) in enumerate(self.PAGES[:7], start=1):
            QShortcut(QKeySequence(f"Ctrl+{i}"), self, activated=lambda k=key: self._switch_page(k))
        QShortcut(QKeySequence("Ctrl+8"), self, activated=lambda: self._switch_page("log"))
        QShortcut(QKeySequence("Ctrl+H"), self, activated=self._open_heatmap)
        QShortcut(QKeySequence("Ctrl+E"), self, activated=self._export_current_figure)

    def _export_current_figure(self):
        panels = {"pca": "pca_panel", "volcano": "volcano_panel", "gsea": "gsea_panel", "ora": "ora_panel",
                  "genes": "gene_panel"}
        panel = getattr(self, panels.get(self.active_key, ""), None)
        if panel is None or self.results is None:
            self.notify("Open a chart page (PCA, Volcano, GSEA, ORA, Gene counts) after running an analysis.")
            return
        panel._open_export_dialog()

    def _gene_to_counts(self, gene: str | None):
        if not gene:
            return
        self.gene_entry.setText(str(gene))
        self._switch_page("genes")
        self._draw_gene_plot()

    def _switch_page(self, key: str, animate: bool = True):
        self.active_key = key
        self._nav_buttons[key].setChecked(True)
        w = self._page_widgets[key]
        self.pages.setCurrentWidget(w)
        if animate:
            fade_in(w)

    @staticmethod
    def _spin(lo, hi, val, step=1, decimals=None, width=84):
        if decimals is None:
            s = QSpinBox()
            s.setRange(int(lo), int(hi))
            s.setSingleStep(int(step))
            s.setValue(int(val))
        else:
            s = QDoubleSpinBox()
            s.setDecimals(decimals)
            s.setRange(lo, hi)
            s.setSingleStep(step)
            s.setValue(val)
        s.setMinimumWidth(width)
        return s

    @staticmethod
    def _combo(items, current=None, editable=False, width=None):
        c = QComboBox()
        c.setEditable(editable)
        c.addItems(list(items))
        if current is not None:
            c.setCurrentText(current)
        if width:
            c.setMinimumWidth(width)
        return c

    def _debounce(self, fn, ms=220):
        t = QTimer(self)
        t.setSingleShot(True)
        t.setInterval(ms)
        t.timeout.connect(fn)
        return lambda *_: t.start()

    def _page_title(self, title: str, subtitle: str) -> QWidget:
        w = QWidget()
        l = QVBoxLayout(w)
        l.setContentsMargins(0, 0, 0, 0)
        l.setSpacing(2)
        t = QLabel(title)
        t.setObjectName("h1")
        s = QLabel(subtitle)
        s.setObjectName("muted")
        l.addWidget(t)
        l.addWidget(s)
        return w

    def _contrast_bar(self, kind: str) -> QWidget:
        w = QWidget()
        l = QVBoxLayout(w)
        l.setContentsMargins(0, 0, 0, 0)
        l.setSpacing(6)
        row = QHBoxLayout()
        lab = QLabel("CONTRAST")
        lab.setObjectName("fieldLabel")
        combo = QComboBox()
        combo.setMinimumWidth(300)
        combo.currentTextChanged.connect(self._on_contrast_selected)
        self._contrast_combos.append(combo)
        row.addWidget(lab)
        row.addWidget(combo)
        row.addStretch()
        l.addLayout(row)
        d = QLabel("")
        d.setObjectName("direction")
        d.setWordWrap(True)
        d.setVisible(False)
        self._direction_labels[kind] = d
        l.addWidget(d)
        return w

    def _result_layout(self, kind: str, title: str, subtitle: str, main: QWidget, inspector: QWidget | None,
                       empty_text: str = "Run an analysis from the Setup page to see this view.",
                       with_contrast: bool = True) -> ResultPage:
        content = QWidget()
        v = QVBoxLayout(content)
        v.setContentsMargins(24, 18, 24, 18)
        v.setSpacing(12)
        head = QHBoxLayout()
        head.addWidget(self._page_title(title, subtitle))
        head.addStretch()
        v.addLayout(head)
        if with_contrast:
            v.addWidget(self._contrast_bar(kind))
        row = QHBoxLayout()
        row.setSpacing(16)
        row.addWidget(main, 1)
        if inspector is not None:
            row.addWidget(_scroll(inspector, 322))
        v.addLayout(row, 1)
        page = ResultPage(content, empty_text)
        self._result_pages.append(page)
        return page

    @staticmethod
    def _inspector(*cards: QWidget) -> QWidget:
        w = QWidget()
        l = QVBoxLayout(w)
        l.setContentsMargins(0, 0, 6, 0)
        l.setSpacing(12)
        for c in cards:
            l.addWidget(c)
        l.addStretch()
        return w

    @staticmethod
    def _table_view(model: DataFrameModel) -> QTableView:
        t = QTableView()
        t.setModel(model)
        t.setAlternatingRowColors(True)
        t.setSelectionBehavior(QTableView.SelectRows)
        t.setEditTriggers(QTableView.NoEditTriggers)
        t.setSortingEnabled(True)
        t.verticalHeader().setVisible(False)
        t.verticalHeader().setDefaultSectionSize(30)
        t.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        t.horizontalHeader().setSortIndicatorShown(True)
        t.horizontalHeader().setSortIndicator(-1, Qt.AscendingOrder)
        t.setShowGrid(False)
        return t

    def _build_setup_page(self) -> QWidget:
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(0, 0, 0, 0)
        inner = QWidget()
        col = QVBoxLayout(inner)
        col.setContentsMargins(28, 22, 28, 24)
        col.setSpacing(14)
        col.addWidget(self._page_title("Setup & parameters", "Load the data, choose the two groups to compare and set the thresholds. Then run the analysis."))

        st = Card()
        sl = QHBoxLayout()
        sl.setSpacing(10)
        self.chip_counts = Chip("○  Counts file not selected")
        self.chip_meta = Chip("○  Metadata file not selected")
        self.chip_groups = Chip("○  Groups not yet extracted")
        for c in (self.chip_counts, self.chip_meta, self.chip_groups):
            sl.addWidget(c)
        sl.addStretch()
        st.body().addLayout(sl)
        col.addWidget(st)

        files = Card("Input files", "Drag & drop a file onto a box, or click it to browse.")
        dz = QHBoxLayout()
        dz.setSpacing(14)
        self.counts_dz = DropZone("Counts file", "csv · tsv · txt · xlsx · xls · json", "+", self._pick_counts)
        self.metadata_dz = DropZone("Metadata file", "csv · tsv · txt · xlsx · xls · json", "+", self._pick_metadata)
        dz.addWidget(self.counts_dz)
        dz.addWidget(self.metadata_dz)
        files.body().addLayout(dz)

        cols = QGridLayout()
        cols.setHorizontalSpacing(14)
        self.gene_col = QLineEdit("1")
        self.sample_col = QLineEdit("1")
        self.condition_col = QLineEdit("3")
        for e in (self.gene_col, self.sample_col, self.condition_col):
            e.setMaximumWidth(150)
        cols.addWidget(field("Gene column", self.gene_col,
                             "Column number (1 = first column) or exact name of the column with the gene "
                             "identifiers in the counts file."), 0, 0)
        cols.addWidget(field("Metadata: sample column", self.sample_col,
                             "Column of the metadata file with the sample IDs (must match the counts columns)."), 0, 1)
        cols.addWidget(field("Metadata: condition column", self.condition_col,
                             "Column of the metadata file with the experimental group of each sample."), 0, 2)
        reset = QPushButton("Reset fields")
        reset.setToolTip("Clears all the fields in this section so you can start over.")
        reset.clicked.connect(self._reset_setup_fields)
        cols.addWidget(reset, 0, 3, alignment=Qt.AlignBottom)
        cols.setColumnStretch(4, 1)
        files.body().addLayout(cols)
        hint = QLabel("Columns accept a number (1 = first column) or the exact column name.")
        hint.setObjectName("cardHint")
        files.body().addWidget(hint)
        col.addWidget(files)

        meth = Card("Statistical method", "Standard RNA-seq count data require DESeq2.")
        mr = QHBoxLayout()
        mr.setSpacing(22)
        self.rb_rna = QRadioButton("RNA-seq → DESeq2 (Wald test)")
        self.rb_rna.setChecked(True)
        self.rb_rna.setToolTip("Raw RNA-seq counts (integers, not normalized). Uses DESeq2.")
        self.rb_micro = QRadioButton("Microarray → limma (eBayes)")
        self.rb_micro.setToolTip("Microarray intensities (typically log-transformed). Uses limma.")
        self.species = self._combo(engine.SPECIES_CHOICES.keys(), width=230)
        self.species.setToolTip(
            "Species of the data: decides how Ensembl/Entrez IDs become gene symbols and which MSigDB sets are "
            "used for GSEA/ORA. For mouse, rat or zebrafish the annotation package must be installed in R "
            "(org.Mm.eg.db, org.Rn.eg.db, org.Dr.eg.db).")
        mr.addWidget(self.rb_rna)
        mr.addWidget(self.rb_micro)
        mr.addSpacing(10)
        mr.addWidget(field("Organism", self.species))
        mr.addStretch()
        meth.body().addLayout(mr)
        col.addWidget(meth)

        grp = Card("Groups to compare")
        top = QHBoxLayout()
        self.extract_btn = QPushButton("Extract groups from metadata")
        self.extract_btn.setToolTip("Reads the metadata file and lists every condition with its sample count.")
        self.extract_btn.clicked.connect(self._extract_groups)
        top.addWidget(self.extract_btn)
        top.addStretch()
        grp.body().addLayout(top)

        g = QGridLayout()
        g.setHorizontalSpacing(14)
        self.test_combo = self._combo([], width=250)
        self.test_combo.setToolTip("log2FC > 0 (and positive NES in GSEA) means HIGHER in this group.")
        self.ref_combo = self._combo([], width=250)
        self.ref_combo.setToolTip("The baseline: the group the test group is compared against (denominator).")
        self.test_combo.setPlaceholderText("Select a group…")
        self.ref_combo.setPlaceholderText("Select a group…")
        self.swap_btn = QPushButton("Swap")
        g.addWidget(field("Test group (numerator)", self.test_combo), 0, 0)
        g.addWidget(field("Reference group (denominator)", self.ref_combo), 0, 1)
        g.addWidget(self.swap_btn, 0, 2, alignment=Qt.AlignBottom)
        g.setColumnStretch(3, 1)
        grp.body().addLayout(g)
        self.contrast_preview = QLabel("")
        self.contrast_preview.setObjectName("direction")
        self.contrast_preview.setWordWrap(True)
        grp.body().addWidget(self.contrast_preview)
        self.groups_info = QLabel("")
        self.groups_info.setObjectName("cardHint")
        self.groups_info.setWordWrap(True)
        grp.body().addWidget(self.groups_info)

        self.pairwise = QCheckBox("Compare ALL pairs of groups (pairwise)")
        self.pairwise.setToolTip("Analyzes every possible pair of groups instead of one two-group comparison.")
        grp.body().addWidget(self.pairwise)
        self.gene_list = QLineEdit()
        self.gene_list.setPlaceholderText("e.g. TP53, BRCA1, MYC   (empty = all genes)")
        self.gene_list.setToolTip(
            "Optional: gene symbols separated by commas. If filled, ONLY these genes are used for the analysis "
            "(PCA included); ORA/GSEA are skipped. Leave empty to analyze all genes.")
        grp.body().addWidget(field("Genes to analyze (optional)", self.gene_list))
        col.addWidget(grp)

        self.swap_btn.clicked.connect(self._swap_groups)
        self.test_combo.currentTextChanged.connect(lambda *_: self._update_contrast_preview())
        self.ref_combo.currentTextChanged.connect(lambda *_: self._update_contrast_preview())
        self.pairwise.toggled.connect(lambda *_: self._toggle_pairwise())
        self._toggle_pairwise()

        par = Card("Analysis parameters")
        pg = QGridLayout()
        pg.setHorizontalSpacing(18)
        pg.setVerticalSpacing(12)
        self.padj_cutoff = self._spin(0.001, 1.0, 0.05, 0.005, 3)
        self.padj_cutoff.setToolTip("Significance threshold on the multiple-testing corrected p-value. Typical: 0.05.")
        self.lfc_cutoff = self._spin(0.0, 10.0, 1.0, 0.25, 2)
        self.lfc_cutoff.setToolTip("Minimum expression change (log2) to call a gene up/down. Typical: 1 (= 2-fold).")
        self.gsea_category = self._combo(["H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"], "H", editable=True)
        self.gsea_category.setToolTip("MSigDB collection for GSEA/ORA. 'H' (Hallmark) is a good generic start.")
        self.gsea_subcategory = QLineEdit()
        self.gsea_subcategory.setPlaceholderText("e.g. CP:KEGG")
        self.gsea_subcategory.setToolTip("Optional: restricts the category to a subcategory (e.g. CP:KEGG in C2).")
        pg.addWidget(field("q-value (padj) cutoff", self.padj_cutoff), 0, 0)
        pg.addWidget(field("log2FoldChange cutoff", self.lfc_cutoff), 0, 1)
        pg.addWidget(field("MSigDB category", self.gsea_category), 0, 2)
        pg.addWidget(field("Subcategory (optional)", self.gsea_subcategory), 0, 3)

        self.min_counts = QSlider(Qt.Horizontal)
        self.min_counts.setRange(0, 500)
        self.min_counts.setValue(10)
        self.min_counts_lbl = QLabel("10")
        self.min_counts_lbl.setMinimumWidth(34)
        self.min_counts_lbl.setStyleSheet(f"font-weight: 700; color: {PALETTE['ink']};")
        self.min_counts.valueChanged.connect(lambda v: self.min_counts_lbl.setText(str(v)))
        mcw = QWidget()
        mcl = QHBoxLayout(mcw)
        mcl.setContentsMargins(0, 0, 0, 0)
        mcl.addWidget(self.min_counts, 1)
        mcl.addWidget(self.min_counts_lbl)
        mcf = field("Minimum counts filter (RNA-seq)", mcw,
                    "Removes genes whose total count across the compared samples is below this value. "
                    "Default 10. Not applied to microarray data.")
        pg.addWidget(mcf, 1, 0, 1, 2)

        self.run_ora = QCheckBox("Run ORA (over-representation)")
        self.run_ora.setChecked(True)
        self.run_ora.setToolTip("Also runs ORA on the significant genes (same cutoffs and MSigDB category).")
        self.ora_dir = self._combo(self.ORA_DIRECTION_DISPLAY)
        self.ora_dir.setToolTip("Which significant genes are tested: all DEGs, only up or only down.")
        pg.addWidget(self.run_ora, 2, 0, 1, 1, alignment=Qt.AlignBottom)
        pg.addWidget(field("ORA gene list", self.ora_dir), 2, 1, 1, 3)
        par.body().addLayout(pg)
        col.addWidget(par)

        act = Card()
        ar = QHBoxLayout()
        ar.setSpacing(10)
        self.run_button = QPushButton("Run analysis")
        self.run_button.setObjectName("primary")
        self.run_button.setToolTip("Runs DESeq2/limma + ORA + GSEA with the parameters above (Ctrl+R).")
        self.run_button.clicked.connect(self._run_analysis)
        self.progress = QProgressBar()
        self.progress.setTextVisible(False)
        self.progress.setFixedWidth(220)
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        save_b = QPushButton("Save settings")
        save_b.setToolTip("Saves files, columns, groups, genes, cutoffs and GSEA/ORA options to a JSON file.")
        save_b.clicked.connect(self._save_settings)
        load_b = QPushButton("Load settings")
        load_b.setToolTip("Loads a settings JSON and fills every field to repeat the same analysis.")
        load_b.clicked.connect(self._load_settings)
        self.report_button = QPushButton("PDF report")
        self.report_button.setToolTip("Detailed PDF of the finished analysis: parameters, PCA, volcano, top DEGs, "
                                      "GSEA, ORA and log. Enabled once an analysis has finished.")
        self.report_button.setEnabled(False)
        self.report_button.clicked.connect(self._export_pdf_report)
        ar.addWidget(self.run_button)
        ar.addWidget(self.progress)
        ar.addStretch()
        ar.addWidget(save_b)
        ar.addWidget(load_b)
        ar.addWidget(self.report_button)
        act.body().addLayout(ar)
        col.addWidget(act)
        col.addStretch()

        redraw = self._debounce(self._on_cutoffs_changed)
        self.padj_cutoff.valueChanged.connect(redraw)
        self.lfc_cutoff.valueChanged.connect(redraw)

        v.addWidget(_scroll(inner))
        return page

    def _pick_counts(self):
        p = open_file_dialog(self, "Select the counts file")
        if p:
            self.counts_dz.set_path(p)

    def _pick_metadata(self):
        p = open_file_dialog(self, "Select the metadata file")
        if p:
            self.metadata_dz.set_path(p)

    @property
    def counts_path(self) -> str:
        return self.counts_dz.path()

    @property
    def metadata_path(self) -> str:
        return self.metadata_dz.path()

    def _method(self) -> str:
        return "RNAseq" if self.rb_rna.isChecked() else "Microarray"

    def _species_value(self) -> str:
        return engine.SPECIES_CHOICES.get(self.species.currentText(), engine.DEFAULT_SPECIES)

    @staticmethod
    def _parse_col(value: str):
        value = value.strip()
        return int(value) if value.isdigit() else value

    def _gene_list(self) -> list[str]:
        return [g.strip() for g in re.split(r"[,;\n]+", self.gene_list.text()) if g.strip()]

    def _reset_setup_fields(self):
        self.counts_dz.set_path("")
        self.metadata_dz.set_path("")
        self.gene_col.setText("1")
        self.sample_col.setText("1")
        self.condition_col.setText("2")
        self.gene_list.clear()
        self.test_combo.clear()
        self.ref_combo.clear()
        self.available_groups = []
        self.groups_info.clear()
        self._validate_setup()
        self.statusBar().showMessage("Fields reset.")

    def _validate_setup(self):
        c, m = self.counts_path.strip(), self.metadata_path.strip()
        if c and os.path.isfile(c):
            self.chip_counts.set_state("ok", "●  Counts file selected")
        elif c:
            self.chip_counts.set_state("err", "●  Counts file path not found")
        else:
            self.chip_counts.set_state("idle", "○  Counts file not selected")
        if m and os.path.isfile(m):
            self.chip_meta.set_state("ok", "●  Metadata file selected")
        elif m:
            self.chip_meta.set_state("err", "●  Metadata file path not found")
        else:
            self.chip_meta.set_state("idle", "○  Metadata file not selected")
        if not self.run_button.property("running"):
            self.run_button.setEnabled(bool(c and m))

    def _toggle_pairwise(self):
        on = self.pairwise.isChecked()
        for w in (self.test_combo, self.ref_combo, self.swap_btn):
            w.setEnabled(not on)
        self._update_contrast_preview()

    def _update_contrast_preview(self):
        if self.pairwise.isChecked():
            txt = ('Pairwise mode: every pair is analysed as "A vs B" (A first alphabetically). '
                   "log2FC > 0 = higher in A.")
        else:
            t, r = self.test_combo.currentText(), self.ref_combo.currentText()
            if not t or not r:
                txt = "Choose the test group and the reference group (after extracting the groups)."
            elif t == r:
                txt = "●  The test group and the reference group must be different."
            else:
                txt = (f'Contrast: {t} vs {r}   →   log2FC > 0 = higher in "{t}"   |   '
                       f'log2FC < 0 = higher in "{r}" (reference)')
        self.contrast_preview.setText(txt)

    def _swap_groups(self):
        t, r = self.test_combo.currentText(), self.ref_combo.currentText()
        self.test_combo.setCurrentText(r)
        self.ref_combo.setCurrentText(t)

    def _extract_groups(self):
        if not self.metadata_path:
            QMessageBox.critical(self, "Error", "Please select the metadata file first.")
            return
        self.statusBar().showMessage("Extracting groups…")
        self._busy(True)
        meta, sc, cc = self.metadata_path, self._parse_col(self.sample_col.text()), self._parse_col(
            self.condition_col.text())
        w = Worker(lambda cb: engine.list_groups(meta, sample_col=sc, condition_col=cc, log_callback=cb),
                   "Group extraction failed", self)
        w.log.connect(self._log)
        w.done.connect(self._on_groups_ready)
        w.failed.connect(self._on_error)
        w.finished.connect(lambda: (self._busy(False), self.statusBar().showMessage("Ready.")))
        self._start(w)

    def _on_groups_ready(self, data: dict):
        self.available_groups = data.get("groups", [])
        counts = data.get("counts", {})
        for c in (self.test_combo, self.ref_combo):
            c.blockSignals(True)
            c.clear()
            c.addItems(self.available_groups)
            c.setCurrentIndex(-1)
            c.blockSignals(False)
        self.groups_info.setText("Samples per group:  " + "   ·   ".join(
            f"{g} (n={counts.get(g, '?')})" for g in self.available_groups))
        pending, self._pending_group_restore = self._pending_group_restore, None
        if pending:
            t, r = pending
            missing = [g for g in (t, r) if g and g not in self.available_groups]
            if missing:
                QMessageBox.warning(self, "Settings", "These groups from the settings file are not in the "
                                                      "metadata: " + ", ".join(missing))
            else:
                self.test_combo.setCurrentText(t or "")
                self.ref_combo.setCurrentText(r or "")
        self._update_contrast_preview()
        self._log(f">>> Groups found: {', '.join(self.available_groups)}")
        self.statusBar().showMessage(f"{len(self.available_groups)} groups found.")
        if self.available_groups:
            self.chip_groups.set_state("ok", f"●  {len(self.available_groups)} groups extracted")
        else:
            self.chip_groups.set_state("err", "●  No groups found in metadata")

    def _collect_settings(self) -> dict:
        return {
            "counts_path": self.counts_path, "metadata_path": self.metadata_path, "method": self._method(),
            "gene_col": self.gene_col.text(), "sample_col": self.sample_col.text(),
            "condition_col": self.condition_col.text(), "pairwise_all": self.pairwise.isChecked(),
            "test_group": self.test_combo.currentText(), "ref_group": self.ref_combo.currentText(),
            "gene_list": self._gene_list(), "species": self._species_value(),
            "padj_cutoff": self.padj_cutoff.value(), "lfc_cutoff": self.lfc_cutoff.value(),
            "min_counts": self.min_counts.value(), "gsea_category": self.gsea_category.currentText(),
            "gsea_subcategory": self.gsea_subcategory.text(), "run_ora": self.run_ora.isChecked(),
            "ora_direction": self.ORA_DIRECTION_MAP.get(self.ora_dir.currentText(), "all"),
            "top_n_pathways": self.gsea_topn.value(), "top_n_ora": self.ora_topn.value(),}

    def _save_settings(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save analysis settings", "analysis_settings.json",
                                              "JSON settings (*.json)")
        if not path:
            return
        payload = {"format": self.SETTINGS_FORMAT, "version": self.SETTINGS_VERSION,
                   "saved": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                   "settings": self._collect_settings()}
        try:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2, ensure_ascii=False)
        except OSError as e:
            QMessageBox.critical(self, "Error", f"Unable to save the settings: {e}")
            return
        self._log(f">>> Settings saved to {path}")
        self.notify("Settings saved")

    def _load_settings(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load analysis settings", "",
                                              "JSON settings (*.json);;All files (*)")
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as fh:
                payload = json.load(fh)
            if not isinstance(payload, dict) or payload.get("format") != self.SETTINGS_FORMAT:
                raise ValueError("this is not a settings file created by this application.")
            if int(payload.get("version", 0)) > self.SETTINGS_VERSION:
                raise ValueError("the file was created by a newer version of the application.")
            s = payload["settings"]
            if not isinstance(s, dict):
                raise ValueError("invalid 'settings' section.")
        except (OSError, ValueError, KeyError, TypeError) as e:
            QMessageBox.critical(self, "Error", f"Unable to load the settings: {e}")
            return

        warnings: list[str] = []

        def put(setter, key, cast=None):
            if key not in s or s[key] is None:
                return
            try:
                setter(cast(s[key]) if cast else s[key])
            except (ValueError, TypeError):
                warnings.append(f"invalid value for '{key}' (ignored)")

        put(self.counts_dz.set_path, "counts_path", str)
        put(self.metadata_dz.set_path, "metadata_path", str)
        (self.rb_micro if s.get("method") == "Microarray" else self.rb_rna).setChecked(True)
        inv_sp = {v: k for k, v in engine.SPECIES_CHOICES.items()}
        if s.get("species") in inv_sp:
            self.species.setCurrentText(inv_sp[s["species"]])
        put(self.gene_col.setText, "gene_col", str)
        put(self.sample_col.setText, "sample_col", str)
        put(self.condition_col.setText, "condition_col", str)
        put(self.pairwise.setChecked, "pairwise_all", bool)
        put(self.padj_cutoff.setValue, "padj_cutoff", float)
        put(self.lfc_cutoff.setValue, "lfc_cutoff", float)
        put(self.min_counts.setValue, "min_counts", int)
        put(self.gsea_category.setCurrentText, "gsea_category", str)
        put(self.gsea_subcategory.setText, "gsea_subcategory", str)
        put(self.run_ora.setChecked, "run_ora", bool)
        put(self.gsea_topn.setValue, "top_n_pathways", int)
        put(self.ora_topn.setValue, "top_n_ora", int)
        self.gene_list.setText(", ".join(str(g) for g in (s.get("gene_list") or [])))
        inv = {v: k for k, v in self.ORA_DIRECTION_MAP.items()}
        self.ora_dir.setCurrentText(inv.get(s.get("ora_direction", "all"), self.ORA_DIRECTION_DISPLAY[0]))

        for c in (self.test_combo, self.ref_combo):
            c.clear()
        self.available_groups = []
        self._pending_group_restore = None
        for label, key in (("counts file", "counts_path"), ("metadata file", "metadata_path")):
            p = s.get(key)
            if p and not os.path.isfile(p):
                warnings.append(f"the {label} was not found: {p}")
        meta = self.metadata_path
        if meta and os.path.isfile(meta):
            if s.get("test_group") or s.get("ref_group"):
                self._pending_group_restore = (s.get("test_group", ""), s.get("ref_group", ""))
            self._extract_groups()
        else:
            warnings.append("groups were not restored (metadata file unavailable): select the file and press "
                            "'Extract groups from metadata'.")
        self._validate_setup()
        self._update_contrast_preview()
        self._log(f">>> Settings loaded from {path}")
        self.notify("Settings loaded")
        if warnings:
            QMessageBox.warning(self, "Settings loaded with warnings", "\n".join("- " + w for w in warnings))

    def _busy(self, on: bool):
        if on:
            self.progress.setRange(0, 0)
        else:
            self.progress.setRange(0, 1)
            self.progress.setValue(0)

    def _start(self, worker: Worker):
        self._workers.append(worker)
        worker.finished.connect(lambda w=worker: self._workers.remove(w) if w in self._workers else None)
        worker.start()

    def _run_analysis(self):
        if not self.counts_path or not self.metadata_path:
            QMessageBox.critical(self, "Error", "Please select both the counts file and the metadata file.")
            return
        contrast = None
        if not self.pairwise.isChecked():
            t, r = self.test_combo.currentText(), self.ref_combo.currentText()
            if not t or not r or t == r:
                QMessageBox.critical(self, "Error", "Choose a test group and a DIFFERENT reference group, or "
                                                    "enable the pairwise comparison of all pairs.")
                return
            contrast = [t, r]  
        gene_list = self._gene_list()
        if gene_list and len({g.upper() for g in gene_list}) < 2:
            QMessageBox.critical(self, "Error",
                                 "Enter at least 2 genes, or leave the field empty to analyze all genes.")
            return

        kw = dict(
            counts_path=self.counts_path, metadata_path=self.metadata_path, method=self._method(),
            gene_col=self._parse_col(self.gene_col.text()), sample_col=self._parse_col(self.sample_col.text()),
            condition_col=self._parse_col(self.condition_col.text()), contrast=contrast,
            pairwise_all=self.pairwise.isChecked(), padj_cutoff=self.padj_cutoff.value(),
            lfc_cutoff=self.lfc_cutoff.value(), min_counts=self.min_counts.value(),
            gsea_category=self.gsea_category.currentText(),
            gsea_subcategory=(self.gsea_subcategory.text() or None), run_ora=self.run_ora.isChecked(),
            ora_direction=self.ORA_DIRECTION_MAP.get(self.ora_dir.currentText(), "all"),
            gene_list=gene_list or None, species=self._species_value())

        self._pending_params = {
            "counts_file": self.counts_path, "metadata_file": self.metadata_path, "method": self._method(),
            "gene_col": self.gene_col.text(), "sample_col": self.sample_col.text(),
            "condition_col": self.condition_col.text(), "contrast": list(contrast) if contrast else None,
            "pairwise_all": self.pairwise.isChecked(), "species": self._species_value(),
            "gene_list": list(gene_list), "padj_cutoff": self.padj_cutoff.value(),
            "lfc_cutoff": self.lfc_cutoff.value(), "min_counts": self.min_counts.value(),
            "gsea_category": self.gsea_category.currentText(), "gsea_subcategory": self.gsea_subcategory.text(),
            "run_ora": self.run_ora.isChecked(), "ora_direction": self.ora_dir.currentText(),
            "top_n_pathways": self.gsea_topn.value(), "top_n_ora": self.ora_topn.value(),
            "started": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),}

        self.run_button.setProperty("running", True)
        self.run_button.setEnabled(False)
        self.run_button.setText("Running…")
        self._busy(True)
        self._run_start = time.time()
        self._tick()
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._tick)
        self._timer.start(1000)

        w = Worker(lambda cb: engine.run_dea(log_callback=cb, **kw), "Analysis failed", self)
        w.log.connect(self._log)
        w.done.connect(self._on_results_ready)
        w.failed.connect(self._on_error)
        w.finished.connect(self._stop_run_ui)
        self._start(w)

    def _tick(self):
        if self._run_start is None:
            return
        mm, ss = divmod(int(time.time() - self._run_start), 60)
        self.run_elapsed.setText(f"Elapsed {mm:02d}:{ss:02d}")
        self.statusBar().showMessage("Analysis running (DESeq2/limma + ORA + GSEA)…")

    def _stop_run_ui(self):
        self._busy(False)
        t = getattr(self, "_timer", None)
        if t is not None:
            t.stop()
        self._run_start = None
        self.run_elapsed.setText("")
        self.run_button.setProperty("running", False)
        self.run_button.setText("Run analysis")
        self._validate_setup()
        self.statusBar().showMessage("Ready.")

    def _on_error(self, title: str, err):
        self._log(f"[ERROR] {title}: {err}")
        QMessageBox.critical(self, title, str(err))

    def _on_results_ready(self, results):
        self.results = results
        self.run_params = dict(self._pending_params or {})
        if self._run_start is not None:
            mm, ss = divmod(int(time.time() - self._run_start), 60)
            self.run_params["duration"] = f"{mm:02d}:{ss:02d} (mm:ss)"
        self.report_button.setEnabled(True)
        tags = list(results.per_contrast.keys())
        for c in self._contrast_combos:
            c.blockSignals(True)
            c.clear()
            c.addItems(tags)
            c.blockSignals(False)
        self.active_tag = tags[0] if tags else ""
        self._sync_contrast_combos()
        for p in self._result_pages:
            p.set_ready(True)
        self._log(">>> Analysis completed successfully.")
        self._refresh_pca_color_columns()
        self._refresh_gene_plot_columns()
        self._refresh_all_views()
        self._switch_page("pca")
        self.notify("Analysis completed")

    def _sync_contrast_combos(self):
        for c in self._contrast_combos:
            c.blockSignals(True)
            c.setCurrentText(self.active_tag)
            c.blockSignals(False)

    def _on_contrast_selected(self, text: str):
        if not text or text == self.active_tag:
            return
        self.active_tag = text
        self._sync_contrast_combos()
        self._refresh_all_views()

    def _active_result(self):
        return None if self.results is None else self.results.per_contrast.get(self.active_tag)

    @staticmethod
    def _dir_names(res):
        hi = res.group_high if res is not None and res.group_high else "group 1 (log2FC>0)"
        lo = res.group_low if res is not None and res.group_low else "group 2 (log2FC<0)"
        return hi, lo

    def _contrast_title(self, tag: str) -> str:
        res = self.results.per_contrast.get(tag) if self.results else None
        return res.label if res is not None else str(tag).replace("_", " ")

    def _regulation_options(self, res):
        hi, lo = self._dir_names(res)
        return ["All", f"Higher in {hi} only", f"Higher in {lo} only", "Not significant only"]

    @staticmethod
    def _direction_text(kind: str, res) -> str:
        if res is None:
            return ""
        if not (res.group_high and res.group_low):
            return "The direction of the comparison cannot be determined automatically for this contrast."
        hi, lo = res.group_high, res.group_low
        if kind == "gsea":
            return (f'Contrast {hi} vs {lo}:  NES > 0 = enriched among genes higher in "{hi}"   |   '
                    f'NES < 0 = among genes higher in "{lo}" (reference)')
        if kind == "ora":
            return (f'Contrast {hi} vs {lo}:  UP = genes higher in "{hi}"   |   '
                    f'DOWN = genes higher in "{lo}" (reference)')
        return (f'Contrast {hi} vs {lo}:  log2FC > 0 = higher in "{hi}"   |   '
                f'log2FC < 0 = higher in "{lo}" (reference)')

    def _update_direction_labels(self):
        res = self._active_result()
        for kind, lab in self._direction_labels.items():
            txt = self._direction_text(kind, res)
            lab.setText(txt)
            lab.setVisible(bool(txt))
        if res is not None:
            opts = self._regulation_options(res)
            idx = max(self.table_reg.currentIndex(), 0)
            self.table_reg.blockSignals(True)
            self.table_reg.clear()
            self.table_reg.addItems(opts)
            self.table_reg.setCurrentIndex(min(idx, len(opts) - 1))
            self.table_reg.blockSignals(False)

    def _refresh_all_views(self):
        if self.results is None:
            return
        self._update_direction_labels()
        self._draw_pca()
        self._draw_volcano()
        self._populate_gene_table()
        self._draw_gsea()
        self._draw_ora()
        self._draw_gene_plot()

    def _on_cutoffs_changed(self):
        if self.results is not None:
            self._draw_volcano()
            self._populate_gene_table()
            self._draw_gsea()
            self._draw_ora()

    def _build_pca_page(self) -> QWidget:
        self.pca_labels = QCheckBox("Show sample labels")
        self.pca_labels.setToolTip("Writes the ID of each sample next to its point.")
        self.pca_ellipses = QCheckBox("Confidence ellipses per group")
        self.pca_ellipses.setToolTip("Ellipse around each group (needs ≥ 3 samples per group).")
        self.pca_colorblind = QCheckBox("Colour-blind-friendly palette")
        self.pca_colorblind.setToolTip("Okabe-Ito palette, readable with colour vision deficiency.")
        self.pca_title = QLineEdit()
        self.pca_title.setPlaceholderText("PCA of Count Data")
        disp = Card("Display")
        for w in (self.pca_labels, self.pca_ellipses, self.pca_colorblind):
            disp.body().addWidget(w)
        disp.body().addWidget(field("Custom title", self.pca_title))

        self.pca_color_combo = RefreshCombo(self._refresh_pca_color_columns)
        self.pca_color_combo.addItem(pcol.DEFAULT_COLUMN)
        self.pca_color_combo.setToolTip(
            "Colors the points by any metadata column. Numeric columns get a continuous scale + colour bar; text "
            "columns one color per group. Missing values are grey.")
        self.pca_kind = self._combo(pcol.KIND_CHOICES)
        self.pca_kind.setToolTip("'Auto' decides from the values. Force 'Categorical' for numeric codes that are "
                                 "really groups, or 'Numeric' for a continuous scale.")
        self.pca_cmap = self._combo(pcol.CMAPS.keys(), "Red → Blue")
        colr = Card("Colour points by")
        colr.body().addWidget(field("Metadata column", self.pca_color_combo))
        colr.body().addWidget(field("Type", self.pca_kind))
        colr.body().addWidget(field("Scale (numeric)", self.pca_cmap))

        expl = ("PCA (Principal Component Analysis): each point is a sample. Samples with similar colour/condition "
                "that sit close together have similar expression profiles. PC1 and PC2 report the percentage of "
                "variance explained. Use the toolbar's magnifier to zoom and the annotation tools to add arrows, "
                "boxes or text before exporting the figure for publication.")
        self.pca_panel = PlotPanel(expl, redraw_callback=self._draw_pca)

        for w in (self.pca_labels, self.pca_ellipses, self.pca_colorblind):
            w.toggled.connect(lambda *_: self._draw_pca())
        self.pca_title.editingFinished.connect(self._draw_pca)
        self.pca_color_combo.activated.connect(lambda *_: self._draw_pca())
        self.pca_kind.activated.connect(lambda *_: self._draw_pca())
        self.pca_cmap.activated.connect(lambda *_: self._draw_pca())
        return self._result_layout("pca", "PCA", "Global structure of the samples — one point per sample.",
                                   self.pca_panel, self._inspector(disp, colr), with_contrast=False)

    def _pca_meta_source(self):
        rp = self.run_params or {}
        path = rp.get("metadata_file") or self.metadata_path
        sample_col = self._parse_col(str(rp.get("sample_col", self.sample_col.text())))
        return path, sample_col

    def _read_meta_cached(self, path: str) -> pd.DataFrame:
        mtime = os.path.getmtime(path)
        if self._meta_cache[0] != path or self._meta_cache[1] != mtime:
            self._meta_cache = (path, mtime, pcol.read_metadata(path))
        return self._meta_cache[2]

    def _refresh_pca_color_columns(self):
        combo = getattr(self, "pca_color_combo", None)
        if combo is None:
            return
        cols = []
        path = self.metadata_path.strip()
        if path and os.path.isfile(path):
            try:
                cols = list(self._read_meta_cached(path).columns)
            except Exception:  # noqa: BLE001
                cols = []
        items = [pcol.DEFAULT_COLUMN] + cols
        current = combo.currentText()
        combo.blockSignals(True)
        combo.clear()
        combo.addItems(items)
        combo.setCurrentText(current if current in items else pcol.DEFAULT_COLUMN)
        combo.blockSignals(False)

    def pca_color_settings(self) -> dict:
        return {"column": self.pca_color_combo.currentText(), "kind": self.pca_kind.currentText(),
                "colorblind": self.pca_colorblind.isChecked(), "cmap": self.pca_cmap.currentText()}

    def _draw_pca(self):
        if self.results is None:
            return
        panel = self.pca_panel
        if self._pca_colorbar is not None:
            try:
                self._pca_colorbar.remove()
            except Exception:  # noqa: BLE001
                pass
            self._pca_colorbar = None
        panel.clear()
        df = self.results.pca_data
        if "condition" not in df.columns:
            panel.redraw()
            return
        meta_path, sample_col = self._pca_meta_source()
        try:
            spec = pcol.spec_from_settings(df, meta_path, sample_col, self.pca_color_settings())
        except Exception as e:  # noqa: BLE001
            QMessageBox.warning(self, "PCA colors", f"Cannot color by '{self.pca_color_combo.currentText()}':\n{e}"
                                                    "\n\nFalling back to the condition colors.")
            self.pca_color_combo.setCurrentText(pcol.DEFAULT_COLUMN)
            spec = pcol.condition_spec(df, self.pca_colorblind.isChecked())

        show_ellipses = self.pca_ellipses.isChecked()

        def ellipse_fn(x, y, ax, color):
            if show_ellipses:
                confidence_ellipse(x, y, ax, n_std=1.5, facecolor=color, alpha=0.15, edgecolor=color, linewidth=1.2)

        artists, extra_handles = pcol.draw_pca_points(panel.ax, df, spec, size=80, alpha=0.85, ellipse_fn=ellipse_fn)
        if self.pca_labels.isChecked() and "SampleID" in df.columns:
            for _, r in df.iterrows():
                panel.ax.annotate(str(r["SampleID"]), xy=(r["PC1"], r["PC2"]), xytext=(4, 4),
                                  textcoords="offset points", fontsize=7)
        panel.ax.set_xlabel(self.results.pca_x_label)
        panel.ax.set_ylabel(self.results.pca_y_label)
        panel.ax.set_title(self.pca_title.text().strip() or "PCA of Count Data")

        if spec["mode"] == "numeric":
            if panel.legend_enabled():
                self._pca_colorbar = pcol.add_colorbar(panel.figure, panel.ax, spec)
            if extra_handles and panel.legend_enabled():
                from matplotlib.lines import Line2D
                handles = [Line2D([], [], marker="o", ls="", markerfacecolor=c, markeredgecolor="black", label=lbl)
                           for lbl, c in extra_handles]
                panel.place_legend(handles=handles)
        else:
            panel.place_legend(title=spec["column"])
        panel.apply_common_style()

        if HAS_MPLCURSORS and artists:
            all_sc = [a for a, _ in artists]
            all_sub = [s for _, s in artists]
            col_name = spec["column"]

            def fmt(idx_tuple):
                ai, pi = idx_tuple
                row = all_sub[ai].iloc[pi]
                val = row.get("_val")
                val = "n.a." if pd.isna(val) else val
                return f"{row.get('SampleID', '?')}\n{col_name}: {val}"
            panel.enable_hover(all_sc, lambda idx: fmt(idx) if isinstance(idx, tuple) else str(idx))
        panel.redraw()

    def _build_volcano_page(self) -> QWidget:
        self.volcano_mode = self._combo(self.LABEL_MODE_DISPLAY)
        self.volcano_mode.setToolTip("Which genes get their name written on the chart.")
        self.volcano_topn = self._spin(1, 200, 15)
        self.volcano_font = self._spin(6, 20, 8)
        self.volcano_custom = QLineEdit()
        self.volcano_custom.setPlaceholderText("TP53, MYC, EGFR")
        self.volcano_leader = QCheckBox("Leader lines")
        self.volcano_leader.setChecked(True)
        self.volcano_leader.setToolTip("Thin line connecting each point to its label.")
        lab = Card("Gene labels")
        lab.body().addWidget(field("Mode", self.volcano_mode))
        r = QHBoxLayout()
        r.addWidget(field("Top N", self.volcano_topn))
        r.addWidget(field("Font size", self.volcano_font))
        lab.body().addLayout(r)
        lab.body().addWidget(field("Specific genes", self.volcano_custom, "Comma-separated, e.g. TP53, MYC"))
        lab.body().addWidget(self.volcano_leader)

        self.sw_up = ColorSwatch(PALETTE["up"])
        self.sw_down = ColorSwatch(PALETTE["down"])
        self.sw_ns = ColorSwatch(PALETTE["ns"])
        sw = QHBoxLayout()
        for name, s in (("Up", self.sw_up), ("Down", self.sw_down), ("NS", self.sw_ns)):
            sw.addWidget(QLabel(name))
            sw.addWidget(s)
            sw.addSpacing(8)
        sw.addStretch()
        self.volcano_size = self._spin(4, 80, 18, 1, 0)
        self.volcano_alpha = self._spin(0.1, 1.0, 0.7, 0.05, 2)
        self.volcano_title = QLineEdit()
        app = Card("Appearance")
        app.body().addLayout(sw)
        pr = QHBoxLayout()
        pr.addWidget(field("Point size", self.volcano_size))
        pr.addWidget(field("Alpha", self.volcano_alpha))
        app.body().addLayout(pr)
        app.body().addWidget(field("Custom title", self.volcano_title))
        note = QLabel("Cutoffs set in Setup update this chart live — no need to rerun R.")
        note.setObjectName("cardHint")
        note.setWordWrap(True)
        app.body().addWidget(note)

        expl = ("Volcano plot: each point is a gene. X axis = log2 Fold Change; Y axis = -log10(padj). "
                "Colours follow the cutoffs set in Setup or the custom selections on the right.")
        self.volcano_panel = PlotPanel(expl, redraw_callback=self._draw_volcano)

        self._volcano_df = None
        self.volcano_panel.canvas.mpl_connect("button_press_event", self._on_volcano_click)
        trig = self._debounce(self._draw_volcano, 250)
        for w in (self.volcano_topn, self.volcano_font, self.volcano_size, self.volcano_alpha):
            w.valueChanged.connect(trig)
        self.volcano_mode.activated.connect(lambda *_: self._draw_volcano())
        self.volcano_leader.toggled.connect(lambda *_: self._draw_volcano())
        for e in (self.volcano_custom, self.volcano_title):
            e.editingFinished.connect(self._draw_volcano)
        for s in (self.sw_up, self.sw_down, self.sw_ns):
            s.colorChanged.connect(lambda *_: self._draw_volcano())
        return self._result_layout("volcano", "Volcano plot", "Effect size against statistical significance.",
                                   self.volcano_panel, self._inspector(lab, app))

    def _on_volcano_click(self, ev):
        p, df = self.volcano_panel, self._volcano_df
        if (df is None or df.empty or not ev.dblclick or ev.inaxes is not p.ax or p.annotation_mode != "none"
                or getattr(p.toolbar, "mode", "")):
            return
        pts = p.ax.transData.transform(df[["log2FoldChange", "neglog10_padj"]].to_numpy(dtype=float))
        dist = np.hypot(pts[:, 0] - ev.x, pts[:, 1] - ev.y)
        i = int(np.argmin(dist))
        if dist[i] <= 14:
            self._gene_to_counts(str(df.iloc[i]["GeneID"]))

    def _compute_label_subset(self, df: pd.DataFrame, mode: str) -> pd.DataFrame:
        if mode == "none" or df.empty:
            return df.iloc[0:0]
        sig = df[df["status"] != "Not Significant"].copy()
        if mode == "all_sig":
            return sig
        n = self.volcano_topn.value()
        if mode == "top_padj":
            return sig.sort_values("padj").head(n)
        if mode == "top_lfc":
            sig["_abs_lfc"] = sig["log2FoldChange"].abs()
            return sig.sort_values("_abs_lfc", ascending=False).head(n)
        if mode == "custom":
            genes = [g.strip() for g in self.volcano_custom.text().split(",") if g.strip()]
            return df[df["GeneID"].astype(str).isin(genes)] if genes else df.iloc[0:0]
        return df.iloc[0:0]

    def _draw_volcano(self):
        if self.results is None:
            return
        panel = self.volcano_panel
        panel.clear()
        tag = self.active_tag
        if not tag or tag not in self.results.per_contrast:
            panel.redraw()
            return
        res = self.results.per_contrast[tag]
        padj_c, lfc_c = self.padj_cutoff.value(), self.lfc_cutoff.value()
        df = engine.categorize_volcano(res.res_tbl, padj_c, lfc_c)
        df = df.dropna(subset=["padj"]).copy()
        df["neglog10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
        self._volcano_df = df
        color_map = {"Upregulated": self.sw_up.color(), "Downregulated": self.sw_down.color(),
                     "Not Significant": self.sw_ns.color()}
        size, alpha = self.volcano_size.value(), self.volcano_alpha.value()
        scatters = {}
        for status, color in color_map.items():
            sub = df[df["status"] == status]
            if sub.empty:
                continue
            panel.ax.scatter(sub["log2FoldChange"], sub["neglog10_padj"], c=color, s=size, alpha=alpha,
                             label=status, edgecolor="none")
            scatters[status] = sub
        panel.ax.axvline(lfc_c, linestyle="--", color="#8a8a8a", linewidth=1)
        panel.ax.axvline(-lfc_c, linestyle="--", color="#8a8a8a", linewidth=1)
        panel.ax.axhline(-math.log10(padj_c), linestyle="--", color="#8a8a8a", linewidth=1)
        panel.ax.set_xlabel("log2 Fold Change")
        panel.ax.set_ylabel("-log10(padj)")
        panel.ax.set_title(self.volcano_title.text().strip() or f"Volcano Plot — {self._contrast_title(tag)}")
        panel.place_legend()

        mode = self.LABEL_MODE_MAP.get(self.volcano_mode.currentText(), "none")
        label_df = self._compute_label_subset(df, mode)
        fontsize, leader = self.volcano_font.value(), self.volcano_leader.isChecked()
        for _, row in label_df.iterrows():
            x, y, gene = row["log2FoldChange"], row["neglog10_padj"], str(row.get("GeneID", ""))
            if leader:
                panel.ax.annotate(gene, xy=(x, y), xytext=(8, 8), textcoords="offset points", fontsize=fontsize,
                                  ha="left", va="bottom", arrowprops=dict(arrowstyle="-", lw=0.4, color="#666666"))
            else:
                panel.ax.text(x, y, gene, fontsize=fontsize, ha="left", va="bottom")
        panel.apply_common_style()

        if HAS_MPLCURSORS and scatters:
            try:
                collections = panel.ax.collections[:len(scatters)]
                frames = list(scatters.values())

                def fmt(idx_tuple, frames=frames):
                    ai, pi = idx_tuple
                    return str(frames[ai].iloc[pi]["GeneID"])
                panel.enable_hover(collections, lambda idx: fmt(idx) if isinstance(idx, tuple) else str(idx))
            except Exception:  # noqa: BLE001
                pass
        panel.redraw()

    def _build_table_page(self) -> QWidget:
        content = QWidget()
        v = QVBoxLayout(content)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(10)
        bar = Card()
        r = QHBoxLayout()
        r.setSpacing(10)
        self.table_filter = QLineEdit()
        self.table_filter.setPlaceholderText("Search gene…")
        self.table_filter.setClearButtonEnabled(True)
        self.table_filter.setMinimumWidth(220)
        self.table_reg = QComboBox()
        self.table_reg.addItems(["All", "Upregulated only", "Downregulated only", "Not Significant only"])
        self.table_reg.setMinimumWidth(240)
        exp = QPushButton("Export CSV…")
        exp.clicked.connect(self._export_gene_table)
        hm = QPushButton("Heatmap…")
        hm.setObjectName("primary")
        hm.setStyleSheet("padding: 7px 18px; font-size: 9.5pt;")
        hm.setToolTip("Interactive clustered heatmap of the differentially expressed genes in a separate window.")
        hm.clicked.connect(lambda *_: self._open_heatmap())
        self.hm_sel = QPushButton("Heatmap of selection…")
        self.hm_sel.setEnabled(False)
        self.hm_sel.setToolTip("Select rows in the table (Ctrl/Shift+click) and open them as a custom-gene heatmap.")
        self.hm_sel.clicked.connect(lambda *_: self._open_heatmap(self._selected_genes()))
        self.table_count = QLabel("")
        self.table_count.setObjectName("muted")
        r.addWidget(self.table_filter)
        r.addWidget(self.table_reg)
        r.addWidget(exp)
        r.addWidget(hm)
        r.addWidget(self.hm_sel)
        r.addStretch()
        r.addWidget(self.table_count)
        bar.body().addLayout(r)
        v.addWidget(bar)
        self.gene_model = DataFrameModel(pd.DataFrame(columns=["GeneID", "Regulation", "Significant",
                                                               "log2FoldChange", "pvalue", "padj", "stat"]))
        self.gene_view = self._table_view(self.gene_model)
        self.gene_view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.gene_view.setSelectionMode(QTableView.ExtendedSelection)
        self.gene_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.gene_view.customContextMenuRequested.connect(self._gene_table_menu)
        self.gene_view.doubleClicked.connect(lambda idx: self._gene_to_counts(self.gene_model.gene_at(idx.row())))
        self.gene_view.selectionModel().selectionChanged.connect(
            lambda *_: self.hm_sel.setEnabled(bool(self.gene_view.selectionModel().selectedRows())))
        QShortcut(QKeySequence(QKeySequence.Copy), self.gene_view, activated=self._copy_selected_rows).setContext(
            Qt.WidgetShortcut)
        v.addWidget(self.gene_view, 1)
        self.table_filter.textChanged.connect(self._debounce(self._populate_gene_table, 180))
        self.table_reg.activated.connect(lambda *_: self._populate_gene_table())
        return self._result_layout("table", "Gene table", "Every tested gene with its statistics for the active contrast.",
                                   content, None)

    def _selected_rows(self) -> list[int]:
        return sorted({i.row() for i in self.gene_view.selectionModel().selectedRows()})

    def _selected_genes(self) -> list[str]:
        genes = [self.gene_model.gene_at(r) for r in self._selected_rows()]
        return [g for g in genes if g]

    def _copy_selected_rows(self):
        rows = self._selected_rows()
        if rows:
            QGuiApplication.clipboard().setText(self.gene_model.rows_as_text(rows))
            self.notify(f"{len(rows)} rows copied")

    def _add_to_volcano_labels(self, genes: list[str]):
        if not genes:
            return
        cur = [g.strip() for g in self.volcano_custom.text().split(",") if g.strip()]
        self.volcano_custom.setText(", ".join(dict.fromkeys(cur + genes)))
        self.volcano_mode.setCurrentText("Specific genes (list)")
        self._draw_volcano()
        self.notify(f"{len(genes)} gene(s) labelled on the volcano plot")

    def _gene_table_menu(self, pos):
        genes = self._selected_genes()
        idx = self.gene_view.indexAt(pos)
        if not genes and idx.isValid():
            genes = [g for g in [self.gene_model.gene_at(idx.row())] if g]
        m = QMenu(self)
        if genes:
            first = genes[0]
            m.addAction(f"Open {first} in Gene counts", lambda: self._gene_to_counts(first))
            m.addAction(f"Label on volcano plot ({len(genes)})", lambda: self._add_to_volcano_labels(genes))
            m.addAction(f"Heatmap of selection ({len(genes)})", lambda: self._open_heatmap(genes))
            m.addSeparator()
            m.addAction("Copy gene names", lambda: QGuiApplication.clipboard().setText("\n".join(genes)))
            m.addAction("Copy rows (Ctrl+C)", self._copy_selected_rows)
        else:
            m.addAction("Select one or more rows first").setEnabled(False)
        m.exec(QCursor.pos())

    def _populate_gene_table(self):
        tag = self.active_tag
        if not tag or self.results is None or tag not in self.results.per_contrast:
            self.table_count.setText("")
            return
        res = self.results.per_contrast[tag]
        ranked = res.ranked_tbl
        if ranked is None or ranked.empty:
            self.gene_model.set_data(pd.DataFrame(columns=self.gene_model._raw.columns))
            self.table_count.setText("No genes to show.")
            return
        ranked = engine.annotate_regulation(ranked, res.group_high, res.group_low,
                                            self.padj_cutoff.value(), self.lfc_cutoff.value())
        hi, lo = self._dir_names(res)
        filt = self.table_filter.text().strip().lower()
        if filt:
            ranked = ranked[ranked["GeneID"].astype(str).str.lower().str.contains(filt, regex=False)]
        sel = self.table_reg.currentIndex()
        if sel == 1:
            ranked = ranked[ranked["Regulation"] == f"Higher in {hi}"]
        elif sel == 2:
            ranked = ranked[ranked["Regulation"] == f"Higher in {lo}"]
        elif sel == 3:
            ranked = ranked[ranked["Regulation"] == "Not significant"]

        cols = ["GeneID", "Regulation", "Significant", "log2FoldChange", "pvalue", "padj", "stat"]
        raw = ranked.reindex(columns=cols).copy()
        shown = pd.DataFrame({
            "GeneID": raw["GeneID"].astype(str), "Regulation": raw["Regulation"].astype(str),
            "Significant": raw["Significant"].astype(str),
            "log2FoldChange": raw["log2FoldChange"].map(_round_or_blank),
            "pvalue": raw["pvalue"].map(_fmt_sci), "padj": raw["padj"].map(_fmt_sci),
            "stat": raw["stat"].map(_round_or_blank)})
        colors = [PALETTE["up"] if r == f"Higher in {hi}" else PALETTE["down"] if r == f"Higher in {lo}"
                  else PALETTE["muted"] for r in raw["Regulation"].astype(str)]
        self.gene_model.set_data(raw, shown, colors)
        self.table_count.setText(f"{len(raw):,} genes shown")

    def _export_gene_table(self):
        tag = self.active_tag
        if not tag or self.results is None or tag not in self.results.per_contrast:
            QMessageBox.information(self, "Info", "No results to export.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export gene table", "gene_table.csv", "CSV (*.csv)")
        if path:
            res = self.results.per_contrast[tag]
            out = engine.annotate_regulation(res.ranked_tbl, res.group_high, res.group_low,
                                             self.padj_cutoff.value(), self.lfc_cutoff.value())
            out.to_csv(path, index=False)
            self.notify(f"Table exported to {os.path.basename(path)}")

    def _open_heatmap(self, preset_genes: list[str] | None = None):
        """Heatmap clusterizzata del contrasto attivo, in una finestra a parte."""
        res = self._active_result()
        if self.results is None or res is None:
            QMessageBox.information(self, "Heatmap", "Run an analysis first.")
            return
        nm = self.results.norm_matrix
        if nm is None:
            QMessageBox.information(self, "Heatmap", "The normalized expression matrix is not available for this "
                                                     "run (see the Log page for the reason).")
            return
        if res.ranked_tbl is None or res.ranked_tbl.empty:
            QMessageBox.information(self, "Heatmap", "No genes to show for this contrast.")
            return
        try:
            import heatmap_window_qt as hmap  # finestra heatmap PySide6
        except ImportError as e:
            QMessageBox.warning(self, "Heatmap", f"Cannot load heatmap_window_qt.py:\n{e}")
            return
        pca = self.results.pca_data
        cond = pd.Series(pca["condition"].astype(str).values, index=pca["SampleID"].astype(str).values)
        meta, sample_col = None, 1
        try:
            meta_path, sample_col = self._pca_meta_source()
            if meta_path and os.path.isfile(meta_path):
                meta = self._read_meta_cached(meta_path)
        except Exception:  # noqa: BLE001
            meta = None

        open_gene = self._gene_to_counts

        ctx = {"norm_mat": nm, "norm_label": self.results.norm_label, "ranked": res.ranked_tbl, "cond": cond,
               "group_high": res.group_high, "group_low": res.group_low, "padj": self.padj_cutoff.value(),
               "lfc": self.lfc_cutoff.value(), "title": res.label, "meta": meta, "sample_col": sample_col,
               "on_gene": open_gene, "preset_genes": preset_genes or None}
        self._heatmap = hmap.HeatmapWindow(self, ctx)
        self._heatmap.show()
        self._heatmap.raise_()
        self._heatmap.activateWindow()

    @staticmethod
    def _plot_table_stack(panel: PlotPanel, view: QTableView, table_height: int = 280) -> QScrollArea:
        inner = QWidget()
        lay = QVBoxLayout(inner)
        lay.setContentsMargins(0, 0, 8, 0)
        lay.setSpacing(14)
        lay.addWidget(panel)
        view.setMinimumHeight(table_height)
        lay.addWidget(view)
        sa = QScrollArea()
        sa.setWidgetResizable(True)
        sa.setFrameShape(QFrame.NoFrame)
        sa.setWidget(inner)
        sa.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        return sa

    @staticmethod
    def _fit_plot_height(panel: PlotPanel, n_rows: int, per_row: int = 30, base: int = 230, min_h: int = 480):
        panel.canvas.setMinimumHeight(max(min_h, base + per_row * int(n_rows)))

    def _build_gsea_page(self) -> QWidget:
        self.gsea_topn = self._spin(3, 100, 15)
        self.gsea_stars = QCheckBox("Significance asterisks")
        self.gsea_stars.setChecked(True)
        self.gsea_stars.setToolTip("* padj<0.05   ** padj<0.01   *** padj<0.001 next to each bar.")
        self.gsea_title = QLineEdit()
        exp = QPushButton("Export pathway table…")
        exp.clicked.connect(self._export_gsea_table)
        opts = Card("Display")
        opts.body().addWidget(field("Top N pathways", self.gsea_topn))
        opts.body().addWidget(self.gsea_stars)
        opts.body().addWidget(field("Custom title", self.gsea_title))
        opts.body().addWidget(exp)

        expl = ("GSEA (Gene Set Enrichment Analysis): shows the most modulated biological pathways. NES positive = "
                'enriched among genes higher in the FIRST group of the contrast ("A" in "A vs B"), negative = among '
                "genes higher in the reference group (B). Only pathways with padj below the chosen threshold (or "
                "the top N by |NES| if none is significant) are shown.")
        self.gsea_panel = PlotPanel(expl, figsize=(8.0, 4.5), redraw_callback=self._draw_gsea, min_height=480)
        self.gsea_model = DataFrameModel(pd.DataFrame(columns=["pathway", "NES", "padj", "size"]), align_left=(0,))
        self.gsea_view = self._table_view(self.gsea_model)
        self.gsea_view.horizontalHeader().setSectionResizeMode(0, QHeaderView.Interactive)
        self.gsea_view.horizontalHeader().resizeSection(0, 360)
        split = self._plot_table_stack(self.gsea_panel, self.gsea_view)

        trig = self._debounce(self._draw_gsea, 200)
        self.gsea_topn.valueChanged.connect(trig)
        self.gsea_stars.toggled.connect(lambda *_: self._draw_gsea())
        self.gsea_title.editingFinished.connect(self._draw_gsea)
        return self._result_layout("gsea", "GSEA — pathways", "Gene Set Enrichment Analysis on the ranked gene list.",
                                   split, self._inspector(opts))

    def _draw_gsea(self):
        if self.results is None:
            return
        panel = self.gsea_panel
        panel.clear()
        self._fit_plot_height(panel, 0)
        self.gsea_model.set_data(pd.DataFrame(columns=["pathway", "NES", "padj", "size"]))
        tag = self.active_tag
        if not tag or tag not in self.results.per_contrast:
            panel.redraw()
            return
        gsea_tbl = self.results.per_contrast[tag].gsea_tbl
        if gsea_tbl is None or gsea_tbl.empty:
            panel.ax.text(0.5, 0.5, "GSEA not available for this contrast.", ha="center", va="center")
            panel.redraw()
            return
        plot_data, msg = engine.filter_top_pathways(gsea_tbl, top_n=self.gsea_topn.value(),
                                                    padj_cutoff=self.padj_cutoff.value())
        self._fit_plot_height(panel, len(plot_data))
        colors = plot_data["Direction"].map({"Up-regulated": PALETTE["up"], "Down-regulated": PALETTE["down"]})
        bars = panel.ax.barh(plot_data["pathway_clean"], plot_data["NES"], color=colors)
        panel.ax.set_xlabel("Normalized Enrichment Score (NES)")
        panel.ax.set_title(self.gsea_title.text().strip() or f"Top Modulated Pathways — {self._contrast_title(tag)}")
        panel.ax.axvline(0, color="black", linewidth=0.8)
        panel.ax.tick_params(axis="y", labelsize=8)

        if self.gsea_stars.isChecked():
            for bar, padj in zip(bars, plot_data["padj"]):
                stars = ""
                if pd.notna(padj):
                    stars = "***" if padj < 0.001 else "**" if padj < 0.01 else "*" if padj < 0.05 else ""
                if stars:
                    width = bar.get_width()
                    offset = 0.05 * (abs(width) if width != 0 else 1)
                    panel.ax.text(width + offset if width >= 0 else width - offset,
                                  bar.get_y() + bar.get_height() / 2, stars, va="center",
                                  ha="left" if width >= 0 else "right", fontsize=9, fontweight="bold")
        xmin, xmax = panel.ax.get_xlim()
        panel.ax.set_xlim(left=xmin * 1.15 if xmin < 0 else xmin, right=xmax * 1.15 if xmax > 0 else xmax)

        hi, lo = self._dir_names(self.results.per_contrast[tag])
        panel.place_legend(handles=[mpatches.Patch(color=PALETTE["up"], label=f"NES > 0: higher in {hi}"),
                                    mpatches.Patch(color=PALETTE["down"], label=f"NES < 0: higher in {lo}")])
        panel.apply_common_style()
        panel.redraw()
        self.statusBar().showMessage(msg)

        d = plot_data.sort_values("padj")
        raw = pd.DataFrame({"pathway": d.get("pathway_clean", ""), "NES": d.get("NES"), "padj": d.get("padj"),
                            "size": d.get("size", "")})
        shown = pd.DataFrame({"pathway": raw["pathway"].astype(str), "NES": raw["NES"].map(_round_or_blank),
                              "padj": raw["padj"].map(_fmt_sci), "size": raw["size"].astype(str)})
        self.gsea_model.set_data(raw, shown)

    def _export_gsea_table(self):
        tag = self.active_tag
        if not tag or self.results is None or tag not in self.results.per_contrast:
            QMessageBox.information(self, "Info", "No GSEA results to export.")
            return
        tbl = self.results.per_contrast[tag].gsea_tbl
        if tbl is None:
            QMessageBox.information(self, "Info", "GSEA is not available for this contrast.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export GSEA table", "gsea_table.csv", "CSV (*.csv)")
        if path:
            tbl.to_csv(path, index=False)
            self.notify(f"GSEA table exported to {os.path.basename(path)}")

    def _build_ora_page(self) -> QWidget:
        self.ora_topn = self._spin(3, 100, 15)
        self.ora_type = self._combo(["Bar", "Dot"])
        self.ora_type.setToolTip("Bar: -log10 padj as bar length with k/K, fold enrichment and padj. "
                                 "Dot: Gene Ratio on X, size = genes, colour = padj.")
        self.ora_title = QLineEdit()
        exp = QPushButton("Export ORA table…")
        exp.clicked.connect(self._export_ora_table)
        opts = Card("Display")
        opts.body().addWidget(field("Top N pathways", self.ora_topn))
        opts.body().addWidget(field("Chart type", self.ora_type))
        opts.body().addWidget(field("Custom title", self.ora_title))
        opts.body().addWidget(exp)

        expl = ("ORA (Over-Representation Analysis): tests whether the significant genes (padj and log2FC cutoffs "
                "from Setup) are over-represented in each pathway, against the background of all genes tested. "
                "k/K = genes in the list / genes of the pathway in the universe; FE = fold enrichment. Bar colour = "
                "direction of the supporting genes (UP / DOWN if at least 80% agree, otherwise MIXED).")
        self.ora_panel = PlotPanel(expl, figsize=(8.0, 4.5), redraw_callback=self._draw_ora, min_height=480)
        self._ora_cols = ["pathway", "genes", "k/K", "FE", "padj", "direction", "N_UP", "N_DOWN"]
        self.ora_model = DataFrameModel(pd.DataFrame(columns=self._ora_cols), align_left=(0,))
        self.ora_view = self._table_view(self.ora_model)
        self.ora_view.horizontalHeader().setSectionResizeMode(0, QHeaderView.Interactive)
        self.ora_view.horizontalHeader().resizeSection(0, 300)
        split = self._plot_table_stack(self.ora_panel, self.ora_view)

        self.ora_topn.valueChanged.connect(self._debounce(self._draw_ora, 200))
        self.ora_type.activated.connect(lambda *_: self._draw_ora())
        self.ora_title.editingFinished.connect(self._draw_ora)
        return self._result_layout("ora", "ORA — over-representation", "Are the significant genes enriched in a pathway?",
                                   split, self._inspector(opts))

    def _draw_ora(self):
        if self.results is None:
            return
        panel = self.ora_panel
        if self._ora_colorbar is not None:
            try:
                self._ora_colorbar.remove()
            except Exception: 
                pass
            self._ora_colorbar = None
        panel.clear()
        self._fit_plot_height(panel, 0)
        self.ora_model.set_data(pd.DataFrame(columns=self._ora_cols))
        tag = self.active_tag
        if not tag or tag not in self.results.per_contrast:
            panel.redraw()
            return
        res = self.results.per_contrast[tag]
        ora_tbl = res.ora_tbl
        if ora_tbl is None or ora_tbl.empty:
            if res.ora_meta is None:
                text = "ORA was not run for this analysis."
            else:
                text = "ORA not available for this contrast."
                note = (res.ora_meta.get("note") or "").strip()
                if note:
                    text += "\n\n" + textwrap.fill(note, 70)
            panel.ax.text(0.5, 0.5, text, ha="center", va="center")
            panel.redraw()
            return

        cutoff = self.padj_cutoff.value()
        plot_data, msg = engine.filter_top_ora(ora_tbl, top_n=self.ora_topn.value(), padj_cutoff=cutoff)
        self._fit_plot_height(panel, len(plot_data))
        n_in = (res.ora_meta or {}).get("n_input_genes")
        title = self.ora_title.text().strip() or f"Top Enriched Pathways (ORA) — {self._contrast_title(tag)}"
        dir_colors = {"UP": PALETTE["up"], "DOWN": PALETTE["down"], "MIXED": PALETTE["ns"]}

        if self.ora_type.currentText() == "Dot":
            from matplotlib.colors import LinearSegmentedColormap
            d = plot_data.sort_values("GeneRatioNum", kind="stable")
            cmap = LinearSegmentedColormap.from_list("ora_padj", [PALETTE["up"], PALETTE["down"]])
            sc = panel.ax.scatter(d["GeneRatioNum"], d["pathway_clean"], s=20 + 8 * d["Count"],
                                  c=d["p.adjust"], cmap=cmap)
            self._ora_colorbar = panel.figure.colorbar(sc, ax=panel.ax, pad=0.02)
            self._ora_colorbar.set_label("padj")
            panel.ax.set_xlabel("Gene Ratio")
        else:
            colors = plot_data["Direction"].map(dir_colors).fillna(PALETTE["ns"])
            bars = panel.ax.barh(plot_data["pathway_clean"], plot_data["neglog10"], color=colors, height=0.7)
            xmax = max(float(plot_data["neglog10"].max()), -np.log10(cutoff))
            if not np.isfinite(xmax) or xmax <= 0:
                xmax = 1.0
            for bar, (_, r) in zip(bars, plot_data.iterrows()):
                label = (f"{int(r['Count'])}/{int(r['SetSize'])} genes ({int(r['N_UP'])} UP, {int(r['N_DOWN'])} DOWN)"
                         f" | FE {r['FoldEnrichment']:.1f} | padj {r['p.adjust']:.1e}")
                panel.ax.text(bar.get_width() + 0.01 * xmax, bar.get_y() + bar.get_height() / 2, label,
                              va="center", ha="left", fontsize=7)
            panel.ax.axvline(-np.log10(cutoff), linestyle="--", color=PALETTE["err"], linewidth=1)
            panel.ax.set_xlim(0, xmax * 1.9)
            panel.ax.set_xlabel(r"$-\log_{10}$(padj)")
            hi, lo = self._dir_names(res)
            panel.place_legend(handles=[
                mpatches.Patch(color=PALETTE["up"], label=f"UP: higher in {hi} (>= 80% of genes)"),
                mpatches.Patch(color=PALETTE["down"], label=f"DOWN: higher in {lo} (>= 80% of genes)"),
                mpatches.Patch(color=PALETTE["ns"], label="MIXED")])
        panel.ax.tick_params(axis="y", labelsize=8)
        panel.ax.set_title(title)
        panel.apply_common_style()
        panel.redraw()
        self.statusBar().showMessage(msg + (f" Genes in input: {n_in}." if n_in else ""))

        d = plot_data.sort_values("p.adjust", kind="stable")
        raw = pd.DataFrame({
            "pathway": d.get("pathway_clean", ""), "genes": d["Count"].astype(int),
            "k/K": d["Count"].astype(float) / d["SetSize"].astype(float),
            "FE": d["FoldEnrichment"].astype(float), "padj": d["p.adjust"].astype(float),
            "direction": d.get("Direction", ""), "N_UP": d["N_UP"].astype(int), "N_DOWN": d["N_DOWN"].astype(int)})
        shown = pd.DataFrame({
            "pathway": raw["pathway"].astype(str), "genes": raw["genes"].astype(str),
            "k/K": d["Count"].astype(int).astype(str) + "/" + d["SetSize"].astype(int).astype(str),
            "FE": raw["FE"].round(2).astype(str), "padj": raw["padj"].map(_fmt_sci),
            "direction": raw["direction"].astype(str), "N_UP": raw["N_UP"].astype(str),
            "N_DOWN": raw["N_DOWN"].astype(str)})
        self.ora_model.set_data(raw, shown)

    def _export_ora_table(self):
        tag = self.active_tag
        if not tag or self.results is None or tag not in self.results.per_contrast:
            QMessageBox.information(self, "Info", "No ORA results to export.")
            return
        tbl = self.results.per_contrast[tag].ora_tbl
        if tbl is None or tbl.empty:
            QMessageBox.information(self, "Info", "ORA is not available for this contrast.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export ORA table", "ora_table.csv", "CSV (*.csv)")
        if path:
            tbl.to_csv(path, index=False)
            self.notify(f"ORA table exported to {os.path.basename(path)}")

    def _build_genes_page(self) -> QWidget:
        self.gene_entry = QLineEdit()
        self.gene_entry.setPlaceholderText("e.g. TP53")
        self.gene_entry.setToolTip("Gene symbol as in the analysis (case-insensitive). Genes removed by the "
                                   "low-counts filter are not available.")
        self.gene_var = RefreshCombo(self._refresh_gene_plot_columns)
        self.gene_var.setToolTip("Metadata column to compare with expression. Numeric → dot plot with Spearman; "
                                 "text → box plot with Mann-Whitney / Kruskal-Wallis.")
        self.gene_kind = self._combo(pcol.KIND_CHOICES)
        plot_btn = QPushButton("Plot")
        plot_btn.setObjectName("primary")
        plot_btn.clicked.connect(self._draw_gene_plot)
        c1 = Card("Gene & variable")
        c1.body().addWidget(field("Gene", self.gene_entry))
        c1.body().addWidget(field("Metadata (Y axis)", self.gene_var))
        c1.body().addWidget(field("Type", self.gene_kind))
        c1.body().addWidget(plot_btn)

        self.gene_color = RefreshCombo(self._refresh_gene_plot_columns)
        self.gene_color.addItem(gplot.NO_COLOR)
        self.gene_color.setToolTip("Colors the points (one point = one sample) by any metadata column.")
        self.gene_color_kind = self._combo(pcol.KIND_CHOICES)
        self.gene_cmap = self._combo(pcol.CMAPS.keys(), "Red → Blue")
        self.gene_colorblind = QCheckBox("Colour-blind-friendly palette")
        c2 = Card("Colour points by")
        c2.body().addWidget(field("Metadata column", self.gene_color))
        c2.body().addWidget(field("Type", self.gene_color_kind))
        c2.body().addWidget(field("Scale (numeric)", self.gene_cmap))
        c2.body().addWidget(self.gene_colorblind)

        self.gene_labels = QCheckBox("Show sample labels")
        self.gene_trend = QCheckBox("Linear trend line (numeric metadata)")
        self.gene_title = QLineEdit()
        exp = QPushButton("Export plot data…")
        exp.clicked.connect(self._export_gene_plot_data)
        c3 = Card("Display")
        c3.body().addWidget(self.gene_labels)
        c3.body().addWidget(self.gene_trend)
        c3.body().addWidget(field("Custom title", self.gene_title))
        c3.body().addWidget(exp)

        main = QWidget()
        mv = QVBoxLayout(main)
        mv.setContentsMargins(0, 0, 0, 0)
        mv.setSpacing(8)
        self.gene_msg = QLabel("Run an analysis first, then type a gene and press 'Plot'.")
        self.gene_msg.setObjectName("muted")
        self.gene_msg.setWordWrap(True)
        expl = ("Gene counts: each point is a SAMPLE; X axis = expression of the chosen gene (log2 of the "
                "size-factor-normalized counts). With a numeric metadata column the chart is a dot plot with "
                "Spearman's rho; with a text column it is a box plot per group with Mann-Whitney (2 groups) or "
                "Kruskal-Wallis (>2). The normalized values are exploratory only: the differential expression "
                "statistics come from the raw counts model.")
        self.gene_panel = PlotPanel(expl, redraw_callback=self._draw_gene_plot)
        mv.addWidget(self.gene_msg)
        mv.addWidget(self.gene_panel, 1)

        self.gene_entry.returnPressed.connect(self._draw_gene_plot)
        for w in (self.gene_var, self.gene_kind, self.gene_color, self.gene_color_kind, self.gene_cmap):
            w.activated.connect(lambda *_: self._draw_gene_plot())
        for w in (self.gene_colorblind, self.gene_labels, self.gene_trend):
            w.toggled.connect(lambda *_: self._draw_gene_plot())
        self.gene_title.editingFinished.connect(self._draw_gene_plot)
        return self._result_layout("genes", "Gene counts", "Expression of a single gene across all samples.",
                                   main, self._inspector(c1, c2, c3), with_contrast=False)

    def _refresh_gene_plot_columns(self):
        if not hasattr(self, "gene_var") or not hasattr(self, "gene_color"):
            return
        cols, cols_y = [], []
        path = self.metadata_path.strip()
        if path and os.path.isfile(path):
            try:
                meta = self._read_meta_cached(path)
                cols = list(meta.columns)
                sname = pcol.resolve_column(meta, self._parse_col(self.sample_col.text()))
                cols_y = [c for c in cols if c != sname]
            except Exception:  # noqa: BLE001
                cols, cols_y = [], []
        cur_y, cur_c = self.gene_var.currentText(), self.gene_color.currentText()
        self.gene_var.blockSignals(True)
        self.gene_color.blockSignals(True)
        self.gene_var.clear()
        self.gene_var.addItems(cols_y)
        self.gene_var.setCurrentText(cur_y if cur_y in cols_y else (cols_y[0] if cols_y else ""))
        items_c = [gplot.NO_COLOR] + cols
        self.gene_color.clear()
        self.gene_color.addItems(items_c)
        self.gene_color.setCurrentText(cur_c if cur_c in items_c else gplot.NO_COLOR)
        self.gene_var.blockSignals(False)
        self.gene_color.blockSignals(False)

    def _draw_gene_plot(self):
        panel = getattr(self, "gene_panel", None)
        if panel is None:
            return
        if self._gene_colorbar is not None:
            try:
                self._gene_colorbar.remove()
            except Exception:  
                pass
            self._gene_colorbar = None
        panel.clear()
        self._gene_plot_frame = None
        if self.results is None:
            panel.redraw()
            return
        nm = self.results.norm_matrix
        if nm is None:
            self.gene_msg.setText("The normalized expression matrix is not available for this run "
                                  "(see the Log page for the reason).")
            panel.redraw()
            return
        gene = self.gene_entry.text().strip()
        if not gene:
            self.gene_msg.setText("Type a gene symbol and press 'Plot'.")
            panel.redraw()
            return
        var = self.gene_var.currentText()
        if not var:
            self.gene_msg.setText("Select a metadata column for the Y axis.")
            panel.redraw()
            return

        meta_path, sample_col = self._pca_meta_source()
        try:
            meta = self._read_meta_cached(meta_path)
            frame = gplot.build_frame(nm, meta, sample_col, gene, var, self.gene_kind.currentText())
        except Exception as e:  # noqa: BLE001
            self.gene_msg.setText(f"Warning: {e}")
            panel.redraw()
            return
        try:
            spec = gplot.color_spec_for(frame, meta, sample_col, self.gene_color.currentText(),
                                        self.gene_color_kind.currentText(), self.gene_colorblind.isChecked(),
                                        self.gene_cmap.currentText())
        except Exception as e:  # noqa: BLE001
            QMessageBox.warning(self, "Gene counts colors",
                                f"Cannot color by '{self.gene_color.currentText()}':\n{e}\n\n"
                                "Falling back to a single color.")
            self.gene_color.setCurrentText(gplot.NO_COLOR)
            spec = gplot.color_spec_for(frame, meta, sample_col, gplot.NO_COLOR)

        artists, stats, extra = gplot.draw_gene_plot(
            panel.ax, frame, spec, self.results.norm_label, size=60, alpha=0.85,
            title=self.gene_title.text().strip() or None,
            show_sample_labels=self.gene_labels.isChecked(), show_trend=self.gene_trend.isChecked())
        if panel.legend_enabled():
            self._gene_colorbar = gplot.add_legend_or_colorbar(panel.figure, panel.ax, spec, extra,
                                                               legend_fn=panel.place_legend)
        panel.apply_common_style()
        self._gene_plot_frame = (frame, stats)

        note = f"{frame['gene']} vs {var}: {len(frame['df'])} samples plotted"
        if not spec.get("uniform"):
            note += f" | colored by '{spec['column']}'"
            if spec.get("n_missing"):
                note += f" ({spec['n_missing']} without a value, grey)"
        if spec.get("n_unmatched"):
            note += f" | Warning: {spec['n_unmatched']} sample IDs not found in the metadata"
        self.gene_msg.setText(note)

        if HAS_MPLCURSORS and artists:
            all_sc = [a for a, _ in artists]
            all_sub = [sdf for _, sdf in artists]
            col_name = spec["column"]

            def fmt(idx_tuple):
                ai, pi = idx_tuple
                row = all_sub[ai].iloc[pi]
                txt = f"{row.get('SampleID', '?')}\nexpr: {row['expr']:.2f} | {var}: {row['val']}"
                if not spec.get("uniform"):
                    cv = row.get("_val")
                    txt += f"\n{col_name}: {'n.a.' if pd.isna(cv) else cv}"
                return txt
            panel.enable_hover(all_sc, lambda idx: fmt(idx) if isinstance(idx, tuple) else str(idx))
        panel.redraw()

    def _export_gene_plot_data(self):
        if not self._gene_plot_frame:
            QMessageBox.information(self, "Info", "Draw a gene plot first.")
            return
        frame, _stats = self._gene_plot_frame
        path, _ = QFileDialog.getSaveFileName(self, "Save plot data", f"{frame['gene']}_vs_{frame['variable']}.csv",
                                              "CSV (*.csv)")
        if not path:
            return
        out = frame["df"].rename(columns={"expr": f"expression ({self.results.norm_label})",
                                          "val": frame["variable"]})
        try:
            out.to_csv(path, index=False)
        except OSError as e:
            QMessageBox.critical(self, "Error", f"Unable to save the file: {e}")
            return
        self.notify(f"Plot data saved to {os.path.basename(path)}")

    def _build_log_page(self) -> QWidget:
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(24, 18, 24, 18)
        v.setSpacing(12)
        v.addWidget(self._page_title("Log", "Everything the R pipeline printed: messages, warnings and errors."))
        bar = QHBoxLayout()
        clear = QPushButton("Clear log")
        clear.setObjectName("danger")
        clear.clicked.connect(lambda: self.log_text.clear())
        save = QPushButton("Save log to file…")
        save.clicked.connect(self._save_log)
        bar.addWidget(clear)
        bar.addWidget(save)
        bar.addStretch()
        v.addLayout(bar)
        self.log_text = QPlainTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setLineWrapMode(QPlainTextEdit.WidgetWidth)
        self.log_text.setStyleSheet(
            f"QPlainTextEdit {{ background: #0e1f31; color: #d5e4f3; border: none; border-radius: 12px; "
            f"padding: 14px; font-family: {MONO}; font-size: 9.5pt; selection-background-color: #1f6fb2; }}")
        v.addWidget(self.log_text, 1)
        return page

    def _log(self, msg: str):
        self.log_text.appendPlainText(msg)

    def _save_log(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save log", "log.txt", "Text (*.txt)")
        if path:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(self.log_text.toPlainText())
            self.notify(f"Log saved to {os.path.basename(path)}")

    def _export_pdf_report(self):
        if self.results is None or not self.results.per_contrast:
            QMessageBox.information(self, "Info", "Run an analysis first: there are no results to report.")
            return
        default = f"DEA_report_{datetime.datetime.now():%Y%m%d_%H%M}.pdf"
        path, _ = QFileDialog.getSaveFileName(self, "Save PDF report", default, "PDF (*.pdf)")
        if not path:
            return
        results = self.results
        params = dict(self.run_params or {})
        params["pca_color"] = self.pca_color_settings()  # stessa colorazione della PCA a schermo
        if self._gene_plot_frame:
            params["gene_plot"] = {
                "gene": self._gene_plot_frame[0]["gene"], "variable": self.gene_var.currentText(),
                "kind": self.gene_kind.currentText(), "color_by": self.gene_color.currentText(),
                "color_kind": self.gene_color_kind.currentText(), "cmap": self.gene_cmap.currentText(),
                "colorblind": self.gene_colorblind.isChecked(), "labels": self.gene_labels.isChecked(),
                "trend": self.gene_trend.isChecked()}
        log_text = self.log_text.toPlainText()
        self.report_button.setEnabled(False)
        self.statusBar().showMessage("Generating PDF report… this may take a few seconds.")
        self._busy(True)

        def task(_cb):
            report_pdf.build_report(path, results, params, log_text)
            return path

        w = Worker(task, "PDF report failed", self)
        w.done.connect(self._on_report_done)
        w.failed.connect(lambda title, err: (self.statusBar().showMessage("PDF report failed."),
                                             QMessageBox.critical(self, "Error",
                                                                  f"Could not create the PDF report:\n{err}")))
        w.finished.connect(lambda: (self.report_button.setEnabled(True), self._busy(False)))
        self._start(w)

    def _on_report_done(self, path: str):
        self.statusBar().showMessage(f"PDF report saved to: {path}")
        box = QMessageBox(self)
        box.setWindowTitle("Report saved")
        box.setText(f"Report saved to:\n{path}\n\nOpen it now?")
        box.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        if box.exec() == QMessageBox.Yes:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(path).resolve())))

    def notify(self, text: str):
        self.toast.show_message(text)
        self.statusBar().showMessage(text)

    def _check_environment(self):
        problems = engine.check_environment()
        if problems:
            QMessageBox.warning(self, "Environment check",
                                "Some issues were detected:\n\n" + "\n".join(f"- {p}" for p in problems))

    def _show_help(self, tab_key: str):
        content = HELP_CONTENT.get(tab_key)
        if not content:
            return
        dlg = QDialog(self)
        dlg.setWindowTitle(content["title"])
        dlg.resize(700, 600)
        v = QVBoxLayout(dlg)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)
        head = QFrame()
        head.setObjectName("header")
        hl = QHBoxLayout(head)
        hl.setContentsMargins(22, 14, 22, 14)
        t = QLabel(content["title"])
        t.setObjectName("appTitle")
        t.setStyleSheet("font-size: 14pt;")
        hl.addWidget(t)
        v.addWidget(head)
        tabs = QTabWidget()
        for label, key in (("What it means", "meaning"), ("How to use it", "usage")):
            te = QTextEdit()
            te.setReadOnly(True)
            te.setFrameShape(QFrame.NoFrame)
            te.setPlainText(content[key])
            te.setStyleSheet("QTextEdit { padding: 14px; font-size: 10.5pt; background: white; border: none; }")
            tabs.addTab(te, label)
        body = QWidget()
        bl = QVBoxLayout(body)
        bl.setContentsMargins(16, 14, 16, 14)
        bl.addWidget(tabs)
        close = QPushButton("Close")
        close.setObjectName("primary")
        close.clicked.connect(dlg.accept)
        row = QHBoxLayout()
        row.addStretch()
        row.addWidget(close)
        bl.addLayout(row)
        v.addWidget(body, 1)
        dlg.exec()

    def closeEvent(self, e):
        QSettings("DEAExplorer", "Main").setValue("geometry", self.saveGeometry())
        hm = getattr(self, "_heatmap", None)
        if hm is not None:
            try:
                hm.close()
            except RuntimeError:  
                pass
        super().closeEvent(e)

    def resizeEvent(self, e):
        super().resizeEvent(e)
        if hasattr(self, "toast") and self.toast.isVisible():
            self.toast.move(self.centralWidget().width() - self.toast.width() - 24,
                            self.centralWidget().height() - self.toast.height() - 52)

def launch() -> int:
    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("DEA Explorer")
    apply_theme(app)
    script_dir = Path(__file__).resolve().parent
    for name in ("logo.png", "GUI_logo.png"):
        p = script_dir / name
        if p.exists():
            app.setWindowIcon(QIcon(str(p)))
            break
    win = DEAWindow()
    win.show()
    return app.exec()

if __name__ == "__main__":
    sys.exit(launch())
