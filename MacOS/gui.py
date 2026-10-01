from __future__ import annotations
import math
import os
import sys
import threading
import time
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext, colorchooser, simpledialog
import ctypes
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import transforms
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

try:
    import mplcursors
    HAS_MPLCURSORS = True
except ImportError:
    HAS_MPLCURSORS = False

try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
    HAS_DND = True
except ImportError:
    HAS_DND = False

import core_calcolo as engine


# =============================================================================
# Visual theme (palette, fonts, ttk style)
# =============================================================================
PALETTE = {
    "bg": "#f4f8fc",            # general background, very light blue-gray
    "bg_card": "#ffffff",       # panel/box background (crisp white cards)
    "accent": "#123c63",        # main academic blue (deeper, more premium)
    "accent_light": "#2f74ad",
    "accent_pale": "#d7e6f2",
    "accent_pale_2": "#eaf1f8",
    "accent_soft": "#5b9bd5",   # secondary accent used for hover/highlights
    "text": "#1b2b3a",
    "muted": "#8a97a3",
    "border": "#d7e3ee",
    "border_strong": "#b7cee2",
    "up": "#c0392b",
    "down": "#2471a3",
    "ns": "#999b9d",
    "ok": "#1e8449",
    "warn": "#b9770e",
    "err": "#c0392b",
    "drop_active": "#eaf6ee",
}

# Color-blind-friendly palette (Okabe & Ito, 2008) — useful for publication
# figures that are also readable by people with color vision deficiencies.
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
             "#0072B2", "#D55E00", "#CC79A7", "#000000"]

FONT_BASE = ("Trebuchet MS", 10)
FONT_BOLD = ("Trebuchet MS", 10, "bold")
FONT_SMALL = ("Trebuchet MS", 9)
FONT_HEADER = ("Trebuchet MS", 16, "bold")
FONT_SUBHEADER = ("Trebuchet MS", 10)
IS_MAC = sys.platform == "darwin"
FONT_MONO = ("Menlo", 10) if IS_MAC else ("Consolas", 9)
FONT_EMOJI = "Apple Color Emoji" if IS_MAC else "Segoe UI Emoji"

class ColorSwatch(tk.Label):
    """Quadratino colorato cliccabile (tk.Button su macOS ignora il colore di sfondo)."""

    def __init__(self, master, command=None, **kw):
        kw.pop("bd", None)
        kw.setdefault("cursor", "hand2")
        kw.setdefault("relief", "solid")
        super().__init__(master, borderwidth=1, **kw)
        self._command = command
        self.bind("<Button-1>", lambda _e: self._command() if self._command else None)


LEGEND_POSITIONS = ["best", "upper right", "upper left", "lower left", "lower right",
                     "center left", "center right", "upper center", "lower center", "outside"]

# Figure size presets for scientific journals (width, height) in inches.
EXPORT_PRESETS = {
    "Custom": None,
    "Single column (89 mm)": (89 / 25.4, 89 / 25.4 * 0.75),
    "1.5 column (140 mm)": (140 / 25.4, 140 / 25.4 * 0.62),
    "Double column (183 mm)": (183 / 25.4, 183 / 25.4 * 0.5),
    "Square 10x10 cm": (10 / 2.54, 10 / 2.54),
    "Slide 16:9 (25.4x14.3 cm)": (10.0, 5.63),
}

HELP_CONTENT = {
    "setup": {
        "title": "Help — Setup & Parameters",
        "meaning": (
            "This tab prepares the Differential Expression Analysis (DEA): the "
            "statistical comparison of expression levels for each gene between two "
            "conditions (e.g. 'Treated' vs 'Control', or 'Responder' vs "
            "'Non-responder').\n\n"
            "What you need:\n"
            "- The COUNTS file: a gene (rows) x sample (columns) table. For RNA-seq "
            "these must be raw counts (integers, e.g. the number of reads aligned to "
            "each gene), NOT already normalized values (TPM/FPKM), because DESeq2 "
            "applies its own statistical normalization.\n"
            "- The METADATA file: a table linking each sample (same ID column as the "
            "counts) to its experimental condition (e.g. treatment group).\n\n"
            "Statistical method:\n"
            "- DESeq2 (Wald test): models counts with a negative binomial "
            "distribution, suited to the typical noise of sequencing data (variance "
            "growing faster than the mean). It is the correct choice for raw RNA-seq "
            "counts.\n"
            "- limma (eBayes): designed for continuous, already log-transformed "
            "intensity data such as microarrays; it uses Bayesian moderation of "
            "variance to stabilize the estimate when there are few samples per "
            "group.\n\n"
            "Parameters:\n"
            "- padj cutoff: threshold on the multiple-testing CORRECTED p-value "
            "(typically via the Benjamini-Hochberg/FDR method). It is needed because "
            "when testing thousands of genes in parallel, some small p-values would "
            "arise purely by chance: padj keeps the expected proportion of false "
            "positives under control. Example: with padj<0.05 you accept that, on "
            "average, 5% of genes declared significant might not really be so.\n"
            "- log2FoldChange cutoff: threshold on the magnitude of change. log2FC=1 "
            "corresponds to a doubling of expression; log2FC=-1 to a halving. A gene "
            "can be statistically significant but with a biologically negligible "
            "change: this cutoff is used to discard those cases.\n"
            "- MSigDB category (for the subsequent GSEA): 'H' (Hallmark) groups a few "
            "dozen well-characterized and coherent gene sets (e.g. the KRAS signaling "
            "pathway, EMT, inflammatory response); other categories (C2 = curated "
            "pathways such as KEGG/Reactome, C5 = Gene Ontology, etc.) are more "
            "granular or specific.\n"
            "- ORA (Over-Representation Analysis): runs alongside the GSEA, with the "
            "same MSigDB category and the same padj / log2FC cutoffs. You can turn it "
            "off and choose which genes it uses: all DEGs (up + down), only "
            "up-regulated or only down-regulated ones."
        ),
        "usage": (
            "1. Press 'Browse...' to select the counts file and the metadata file "
            "(accepts csv/tsv/txt/xlsx/xls/json).\n"
            "2. Check/correct the column indices or names (genes, samples, "
            "condition): you can enter either a number (1 = first column) or the "
            "exact column name.\n"
            "3. Choose the correct statistical method based on the data type "
            "(RNA-seq -> DESeq2, microarray -> limma).\n"
            "4. Press 'Extract groups from metadata': the list of conditions found "
            "will appear, with the corresponding number of samples (n=...).\n"
            "5. Check EXACTLY 2 groups to compare, or enable the 'Compare ALL pairs "
            "of groups (pairwise)' checkbox if you want all possible contrasts at "
            "once (in that case the individual group checkboxes are disabled).\n"
            "6. Set the padj cutoff, log2FC cutoff, and MSigDB category for the "
            "GSEA and ORA (tick/untick 'Run ORA' and pick which genes it should use).\n"
            "7. The 'Preparation status' panel at the top and the red borders on "
            "required fields show you what is still missing; the 'Run analysis' "
            "button is enabled only once both counts and metadata have been "
            "selected.\n"
            "8. Press 'Run analysis': the progress bar and timer show that the "
            "computation (DESeq2/limma + GSEA in R) is running; when finished you "
            "are automatically switched to the PCA tab.\n"
            "9. 'Reset fields' clears the whole file/columns/groups section so you "
            "can start over."
        ),
    },
    "pca": {
        "title": "Help — PCA",
        "meaning": (
            "PCA (Principal Component Analysis) is a statistical technique that "
            "reduces the thousands of genes measured for each sample to a few "
            "'principal components' (PC1, PC2, ...), linear combinations of genes "
            "chosen to capture as much of the data's variability as possible in the "
            "fewest dimensions.\n\n"
            "In practice: each point in the chart is a sample (not a gene). If two "
            "samples have a similar expression profile, their points will be close; "
            "if they are very different, they will be far apart. PC1 and PC2 report "
            "in parentheses the percentage of total variance they manage to explain: "
            "the higher it is, the more faithfully that 2D chart represents the true "
            "differences between samples (rather than losing information by "
            "compressing too many dimensions into two).\n\n"
            "Biological purpose: to check whether samples group according to the "
            "expected experimental condition (e.g. all 'responders' close together, "
            "separated from 'non-responders'). If instead the groups do not "
            "separate, or a sample is isolated far from all the others, this can "
            "indicate: a batch effect (samples processed on different days/batches), "
            "a low-quality sample, or simply that condition has no strong "
            "transcriptomic effect."
        ),
        "usage": (
            "- 'Show sample labels': writes the ID of each sample next to its "
            "point.\n"
            "- 'Confidence ellipses per group': draws an ellipse around the samples "
            "of each condition (requires at least 3 samples per group) to better "
            "visualize separation.\n"
            "- 'Color-blind-friendly palette': switches the group colors to the "
            "Okabe-Ito palette, readable even by people with color vision "
            "deficiencies.\n"
            "- 'Custom title' + 'Apply': replaces the chart's automatic title.\n"
            "- Zoom/pan: use the matplotlib toolbar above the chart (the "
            "magnifying glass and the hand) to zoom into an area or move around.\n\n"
            "How to add manual annotations (arrows, lines, boxes, text):\n"
            "1. In the 'Manual annotations' bar, choose the tool: Arrow, Line, Box "
            "or Text (or 'None' to disable).\n"
            "2. For Arrow/Line/Box: click on the chart at the starting point and "
            "drag to the end point, then release the mouse button.\n"
            "3. For Text: simply click where you want the text, then type the "
            "content in the window that appears.\n"
            "4. You can change color (colored button next to 'Color') and line "
            "width before drawing the next annotation.\n"
            "5. 'Undo last' removes only the last annotation added; 'Clear all' "
            "removes all of them.\n"
            "Note: if you change an option that redraws the chart from scratch "
            "(e.g. enabling/disabling ellipses), manual annotations are cleared: add "
            "them last, once the chart is already set up.\n\n"
            "To export the figure in high resolution (for a paper/poster): press "
            "'Export image', choose a size preset (single/1.5/"
            "double column), the file format (PDF/SVG/EPS are vector formats, ideal "
            "for crisp text and lines; PNG/TIFF are raster), the DPI (300 is the "
            "minimum standard required by most journals) and press 'Save...'."
        ),
    },
    "volcano": {
        "title": "Help — Volcano Plot",
        "meaning": (
            "The volcano plot combines, for each gene, TWO pieces of information at "
            "the same time: how much its expression changed (X axis) and how "
            "statistically reliable that change is (Y axis).\n\n"
            "- X axis = log2FoldChange: expression change between the two "
            "conditions on a log2 scale. Example: log2FC = 2 means the gene is "
            "expressed 4 times more (2^2) in the reference condition; log2FC = -1 "
            "means expression is halved.\n"
            "- Y axis = -log10(padj): the HIGHER a point is, the lower (hence more "
            "significant) its padj is. -log10 is used only to make very small "
            "p-values readable on a chart (e.g. padj = 0.00001 becomes 5 on the Y "
            "axis).\n\n"
            "The dashed horizontal and vertical lines mark the cutoffs chosen in the "
            "Setup tab: a gene is colored as 'Upregulated' or 'Downregulated' only if "
            "it passes BOTH thresholds (significance AND effect size); otherwise it "
            "stays gray ('Not Significant'). The typical 'butterfly/volcano' shape of "
            "the chart arises because genes with more extreme changes (tails on the "
            "right/left) also tend to be the most significant (points higher up), "
            "hence the name."
        ),
        "usage": (
            "- 'Contrast': if you ran multiple comparisons (e.g. with the pairwise "
            "option), choose here which one to display.\n"
            "- Gene labels: choose the mode ('None', 'All significant genes', 'Top N "
            "by padj', 'Top N by |log2FC|', 'Specific genes (list)'). With 'Specific "
            "genes' type the names separated by commas in the field below (e.g. "
            "TP53, MYC, EGFR). Adjust 'Top N' and 'Font size' as needed; 'Leader "
            "lines to label' helps readability when there are many labels close "
            "together.\n"
            "- Appearance: customize the Up/Down/NS color (click the colored "
            "square), point size and transparency, and the title.\n"
            "- After changing a setting, press 'Apply changes' to regenerate the "
            "chart (there is no need to rerun the whole analysis in R).\n\n"
            "Manual annotations (arrows/lines/boxes/text): same procedure described "
            "for PCA — choose the tool in the 'Manual annotations' bar, drag (or "
            "click for text) directly on the chart. Remember: pressing 'Apply "
            "changes' redraws the chart and clears manual annotations, so add them "
            "last.\n\n"
            "For high-resolution export, use 'Export image' as "
            "in the PCA tab."
        ),
    },
    "table": {
        "title": "Help — Gene Table",
        "meaning": (
            "This tab shows, in tabular form, the same statistical results as the "
            "volcano plot: a complete (or filtered) list of genes with the values "
            "produced by DESeq2/limma.\n\n"
            "Columns:\n"
            "- GeneID: gene identifier (symbol, e.g. KRAS, or Ensembl ID if the "
            "symbol is not available).\n"
            "- Regulation: direction of change (Up-regulated / Down-regulated) "
            "based on the sign of log2FoldChange, regardless of whether it is "
            "significant or not.\n"
            "- Significant: indicates whether the gene passes BOTH the padj cutoff "
            "AND the log2FC cutoff set in the Setup tab (i.e. whether it is colored "
            "in the volcano plot).\n"
            "- log2FoldChange: the magnitude of the change on a log2 scale (see the "
            "Volcano Plot help for details).\n"
            "- pvalue: the 'raw' statistical significance, BEFORE multiple-testing "
            "correction.\n"
            "- padj: the corrected p-value (e.g. Benjamini-Hochberg/FDR method): "
            "this is the value to use to decide whether a gene is significant, "
            "never the raw pvalue, because testing thousands of genes together "
            "would underestimate the false positive rate.\n"
            "- stat: the test statistic (e.g. the Wald statistic for DESeq2, "
            "moderated t for limma) from which the pvalue is derived; the sign "
            "indicates direction, the absolute value the strength of the signal."
        ),
        "usage": (
            "- 'Filter (search GeneID)': type even part of a gene name to filter "
            "the table in real time (e.g. typing 'KRAS' also shows any related "
            "genes containing that string).\n"
            "- Regulation dropdown: shows only Upregulated genes, only "
            "Downregulated, or only Not Significant.\n"
            "- Click a column header (e.g. 'padj') to sort the table by that "
            "column; click again to reverse ascending/descending order.\n"
            "- 'Export CSV...': saves the currently filtered/sorted table to a CSV "
            "file, useful for further analysis in Excel/R/Python or as "
            "supplementary material for a paper.\n"
            "- The number in the bottom right ('N genes shown') updates based on "
            "the applied filters."
        ),
    },
    "gsea": {
        "title": "Help — GSEA / Pathway",
        "meaning": (
            "GSEA (Gene Set Enrichment Analysis) answers a different question than "
            "the single-gene one: instead of asking 'did this gene change "
            "significantly?', it asks 'does an entire group of genes belonging to "
            "the same biological pathway move coherently in one direction, even if "
            "each individual gene on its own is not strong enough to be "
            "significant?'.\n\n"
            "How it works, briefly: all genes are ranked based on a statistic "
            "derived from the comparison between the two conditions; for each gene "
            "set (e.g. 'HALLMARK_KRAS_SIGNALING', a group of genes typically "
            "activated/modulated when the KRAS signaling pathway is hyperactive) it "
            "is checked whether its genes cluster unusually toward the top or "
            "bottom of this ranking, compared to what would be expected by pure "
            "chance.\n\n"
            "- NES (Normalized Enrichment Score): the enrichment score, normalized "
            "for gene set size (so pathways with more or fewer genes are "
            "comparable). Positive NES = the gene set is shifted toward the most "
            "UP-regulated genes; negative NES = toward the most DOWN-regulated "
            "genes.\n"
            "- padj: as for single genes, this is the enrichment p-value corrected "
            "for the number of gene sets tested together.\n"
            "- Asterisks (*, **, ***) indicate padj<0.05, <0.01, <0.001 "
            "respectively, following the most common graphical convention in "
            "publications.\n\n"
            "Important point: a significant NES does NOT imply that the individual "
            "genes of that pathway are significant DEGs one by one (nor that they "
            "are mutated); it means that the overall pathway pattern, summed over "
            "many genes that agree even if weak individually, is statistically "
            "robust."
        ),
        "usage": (
            "- 'Top N pathways': how many pathways to show in the bar chart (in "
            "order of significance/|NES|).\n"
            "- 'Significance asterisks': show/hide the * ** *** symbols next to "
            "each bar.\n"
            "- 'Custom title' + 'Apply': replaces the automatic title.\n"
            "- 'Export pathway table CSV...': saves the full list of pathways "
            "shown (name, NES, padj, gene set size) to a CSV file.\n"
            "- The table at the bottom shows the same data as the chart; click "
            "the headers to sort (e.g. by ascending padj).\n"
            "- In the chart, red bars are 'Up-regulated' pathways (positive NES), "
            "blue bars are 'Down-regulated' (negative NES).\n\n"
            "Manual annotations and high-resolution export: identical to those "
            "described for the PCA tab — use the 'Manual annotations' bar to "
            "highlight a pathway of particular interest with an arrow or box "
            "before exporting the figure with 'Export image'."
        ),
    },
    "ora": {
        "title": "Help — ORA / Over-representation",
        "meaning": (
            "ORA (Over-Representation Analysis) asks a simpler question than GSEA: "
            "'among the genes I declared differentially expressed, are there more "
            "genes of a given pathway than I would expect by chance?'.\n\n"
            "How it works, briefly: the significant genes (padj and log2FC below/above "
            "the cutoffs set in the Setup tab) are taken as a list; for each gene set "
            "a hypergeometric test compares how many of them fall in that pathway "
            "with how many would be expected given the size of the pathway. The "
            "background (universe) is the set of genes actually tested in the DEA, "
            "not the whole genome.\n\n"
            "- k/K: genes of the list that belong to the pathway (k) over the genes "
            "of the pathway present in the universe (K).\n"
            "- FE (fold enrichment): (k / list size) divided by (K / universe size). "
            "FE = 3 means the pathway is 3 times more represented than expected.\n"
            "- padj: hypergeometric p-value corrected for the number of pathways "
            "tested (Benjamini-Hochberg).\n"
            "- Direction: ORA itself has no direction (the list is unordered). It is "
            "derived afterwards from the genes supporting each pathway: UP / DOWN if "
            "at least 80% of them have log2FC > 0 / < 0, MIXED otherwise. N_UP and "
            "N_DOWN give the exact counts.\n\n"
            "ORA and GSEA are complementary, not alternatives: ORA depends on the "
            "cutoffs used to build the gene list and ignores the magnitude of the "
            "changes, while GSEA uses the whole ranking of genes."
        ),
        "usage": (
            "- 'Top N pathways': how many pathways to show (significant ones first; "
            "if none passes the padj cutoff, the top N by p-value are shown as "
            "exploratory).\n"
            "- 'Chart type': 'Bar' shows significance as bar length (-log10 padj), "
            "with k/K, FE and padj written next to each bar; 'Dot' shows the Gene "
            "Ratio on the X axis, dot size = number of genes, color = padj.\n"
            "- 'Custom title' + 'Apply': replaces the automatic title.\n"
            "- 'Export ORA table CSV...': saves the complete ORA table "
            "(all pathways tested, including the gene lists).\n"
            "- Which genes enter the ORA (all DEGs / only up / only down) is chosen "
            "in the Setup tab and requires re-running the analysis; the padj cutoff "
            "of this chart follows the Setup tab and updates without rerunning R.\n"
            "- If the tab says ORA is not available, the reason is shown (typically "
            "fewer than 5 genes passing the cutoffs, or ORA disabled in Setup).\n\n"
            "Manual annotations and high-resolution export: identical to those "
            "described for the PCA tab."
        ),
    },
    "log": {
        "title": "Help — Log",
        "meaning": (
            "This tab shows the text log of everything the R pipeline "
            "(DESeq2/limma + GSEA) printed during execution: informational "
            "messages, warnings, and any errors. It is the same kind of output you "
            "would see running the R script from the command line, useful for "
            "understanding exactly what was done (e.g. how many genes were "
            "filtered out for low counts, which experimental design was used, how "
            "many gene sets were tested in the GSEA) or for diagnosing an error."
        ),
        "usage": (
            "- The log fills in automatically during 'Group extraction' and "
            "'Run analysis'.\n"
            "- 'Clear log': empties the text area (does not delete anything on "
            "disk).\n"
            "- 'Save log to file...': exports the current content to a .txt file, "
            "useful for attaching to a bug report or for your own internal project "
            "documentation."
        ),
    },
}

if sys.platform == "win32":
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass

def apply_theme(root: tk.Tk) -> ttk.Style:
    """Configure a consistent ttk theme (academic palette, clean borders,
    colored tabs, readable tables, validation styles)."""
    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass

    root.configure(bg=PALETTE["bg"])

    style.configure(".", background=PALETTE["bg"], foreground=PALETTE["text"], font=FONT_BASE)
    style.configure("TFrame", background=PALETTE["bg"])
    style.configure("TLabel", background=PALETTE["bg"], foreground=PALETTE["text"], font=FONT_BASE)
    style.configure("Header.TLabel", background=PALETTE["accent"], foreground="white", font=FONT_HEADER)
    style.configure("Subheader.TLabel", background=PALETTE["accent"], foreground=PALETTE["accent_pale"],
                     font=FONT_SUBHEADER)
    style.configure("Muted.TLabel", background=PALETTE["bg"], foreground=PALETTE["muted"], font=FONT_SMALL)
    style.configure("CardMuted.TLabel", background=PALETTE["bg_card"], foreground=PALETTE["muted"],
                     font=FONT_SMALL)
    style.configure("Ok.TLabel", background=PALETTE["bg_card"], foreground=PALETTE["ok"], font=FONT_BOLD)
    style.configure("Warn.TLabel", background=PALETTE["bg_card"], foreground=PALETTE["warn"], font=FONT_BOLD)
    style.configure("Err.TLabel", background=PALETTE["bg_card"], foreground=PALETTE["err"], font=FONT_BOLD)

    style.configure("TLabelframe", background=PALETTE["bg_card"], bordercolor=PALETTE["border_strong"],
                     relief="solid", borderwidth=1)
    style.configure("TLabelframe.Label", background=PALETTE["bg_card"], foreground=PALETTE["accent"],
                     font=("Trebuchet MS", 11, "bold"))

    style.configure("TButton", background=PALETTE["accent_pale_2"], foreground=PALETTE["accent"],
                     padding=(12, 7), font=FONT_BASE, borderwidth=0, relief="flat")
    style.map("TButton",
              background=[("active", PALETTE["accent_pale"]), ("disabled", "#eceff2")],
              foreground=[("disabled", "#9aa4ad")])

    style.configure("Accent.TButton", background=PALETTE["accent"], foreground="white",
                     padding=(16, 9), font=FONT_BOLD, borderwidth=0, relief="flat")
    style.map("Accent.TButton",
              background=[("active", PALETTE["accent_light"]), ("disabled", "#a9bccb")])

    style.configure("Danger.TButton", background="#f8ebe9", foreground=PALETTE["err"],
                     padding=6, font=FONT_BASE, borderwidth=0)
    style.map("Danger.TButton", background=[("active", "#f1d7d3")])

    style.configure("Help.TButton", background=PALETTE["accent_pale"], foreground=PALETTE["accent"],
                     padding=6, font=FONT_BOLD, borderwidth=0)
    style.map("Help.TButton", background=[("active", PALETTE["accent_light"])],
              foreground=[("active", "white")])

    style.configure("Toolbutton", background=PALETTE["bg_card"], foreground=PALETTE["text"],
                     padding=5, font=FONT_SMALL, borderwidth=1, relief="solid",
                     bordercolor=PALETTE["border"])
    style.map("Toolbutton",
              background=[("selected", PALETTE["accent"]), ("active", PALETTE["accent_pale"])],
              foreground=[("selected", "white")])

    style.configure("TNotebook", background=PALETTE["bg"], borderwidth=0, tabmargins=(8, 8, 8, 0))
    style.configure("TNotebook.Tab", background=PALETTE["accent_pale_2"], foreground=PALETTE["accent"],
                     padding=(16, 9), font=FONT_BASE, borderwidth=0)
    style.map("TNotebook.Tab",
              background=[("selected", PALETTE["accent"]), ("active", PALETTE["accent_pale"])],
              foreground=[("selected", "white")])

    style.configure("Treeview", background="white", fieldbackground="white", foreground=PALETTE["text"],
                     rowheight=25, font=FONT_BASE, borderwidth=0)
    style.configure("Treeview.Heading", background=PALETTE["accent"], foreground="white", font=FONT_BOLD,
                     relief="flat", padding=(6, 6))
    style.map("Treeview.Heading", background=[("active", PALETTE["accent_light"])])
    style.map("Treeview", background=[("selected", PALETTE["accent_pale"])],
              foreground=[("selected", PALETTE["text"])])

    style.configure("Horizontal.TProgressbar", background=PALETTE["accent"],
                     troughcolor=PALETTE["accent_pale_2"], borderwidth=0)

    style.configure("TCheckbutton", background=PALETTE["bg_card"], font=FONT_BASE)
    style.configure("TRadiobutton", background=PALETTE["bg_card"], font=("Trebuchet MS", 11, "bold"), foreground="#008026")
    style.configure("TEntry", padding=4, fieldbackground="white")
    style.configure("Invalid.TEntry", padding=4, fieldbackground="#fdecea")
    style.configure("TSpinbox", padding=4)
    style.configure("TCombobox", padding=4)

    style.configure("TPanedwindow", background=PALETTE["bg"])
    style.configure("TScrollbar", background=PALETTE["accent_pale_2"], troughcolor=PALETTE["bg"],
                     bordercolor=PALETTE["border"])

    return style


# =============================================================================
# Tooltip: small contextual help that pops up on hover
# =============================================================================
class Tooltip:
    """Shows a small text box when the mouse hovers over a widget for a
    moment. Used throughout the GUI to explain what a field does, without
    cluttering the interface with fixed text."""

    def __init__(self, widget, text: str, delay: int = 500, wraplength: int = 280):
        self.widget = widget
        self.text = text
        self.delay = delay
        self.wraplength = wraplength
        self._after_id = None
        self._tip = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _event=None):
        self._unschedule()
        self._after_id = self.widget.after(self.delay, self._show)

    def _unschedule(self):
        if self._after_id:
            self.widget.after_cancel(self._after_id)
            self._after_id = None

    def _show(self):
        if self._tip or not self.text:
            return
        x = self.widget.winfo_rootx() + 12
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 8
        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        try:
            self._tip.wm_attributes("-topmost", True)
        except tk.TclError:
            pass
        self._tip.wm_geometry(f"+{x}+{y}")
        tk.Label(self._tip, text=self.text, justify="left", background="#fffee0",
                 relief="solid", borderwidth=1, wraplength=self.wraplength,
                 font=FONT_SMALL, padx=6, pady=4).pack()

    def _hide(self, _event=None):
        self._unschedule()
        if self._tip:
            self._tip.destroy()
            self._tip = None


def add_tip(widget, text: str):
    """Shortcut: attaches a Tooltip to a widget and returns it (useful for
    chaining calls while building the UI)."""
    Tooltip(widget, text)
    return widget


# =============================================================================
# DropZone: a card-style control that lets the user either drag a file
# in from the OS file explorer, or click anywhere on the card to open the
# usual "Browse..." dialog. Falls back cleanly to click-only selection if
# tkinterdnd2 is not installed (HAS_DND == False).
# =============================================================================
class DropZone(tk.Frame):
    """A single-file picker presented as a rounded-looking card. Reflects
    three states visually: empty, file selected & found, path set but
    file missing on disk."""

    def __init__(self, parent, title: str, subtitle: str, path_var: tk.StringVar,
                 pick_command, icon: str = "📄", accepted_ext=None):
        super().__init__(parent, bg=PALETTE["bg_card"], highlightthickness=2,
                          highlightbackground=PALETTE["border_strong"], bd=0, cursor="hand2")
        self.path_var = path_var
        self.pick_command = pick_command
        self.accepted_ext = accepted_ext or (".csv", ".tsv", ".txt", ".xlsx", ".xls", ".json")
        self._bg_idle = PALETTE["bg_card"]
        self._bg_hover = PALETTE["accent_pale_2"]
        self._bg_drop = PALETTE["drop_active"]

        self.icon_lbl = tk.Label(self, text=icon, font=(FONT_EMOJI, 28),
                                  bg=self._bg_idle, fg=PALETTE["accent"])
        self.icon_lbl.pack(pady=(16, 4))

        self.title_lbl = tk.Label(self, text=title, font=FONT_BOLD, bg=self._bg_idle,
                                   fg=PALETTE["accent"])
        self.title_lbl.pack()

        self.subtitle_lbl = tk.Label(self, text=subtitle, font=FONT_SMALL, bg=self._bg_idle,
                                      fg=PALETTE["muted"], wraplength=230, justify="center")
        self.subtitle_lbl.pack(pady=(2, 10))

        self.file_lbl = tk.Label(self, text="No file selected — drag it here or click",
                                  font=FONT_SMALL, bg=self._bg_idle, fg=PALETTE["muted"],
                                  wraplength=230, justify="center")
        self.file_lbl.pack(pady=(0, 8), padx=10)

        self.browse_btn = ttk.Button(self, text="Sfoglia...", command=self._trigger_pick)
        self.browse_btn.pack(pady=(0, 16))

        self._bg_widgets = [self, self.icon_lbl, self.title_lbl, self.subtitle_lbl, self.file_lbl]
        for w in (self, self.icon_lbl, self.title_lbl, self.subtitle_lbl, self.file_lbl):
            w.bind("<Button-1>", lambda _e: self._trigger_pick())
            w.bind("<Enter>", self._on_enter)
            w.bind("<Leave>", self._on_leave)

        self._register_dnd()
        self.path_var.trace_add("write", lambda *a: self._refresh())
        self._refresh()

    # ------------------------------------------------------------------
    def _register_dnd(self):
        if not HAS_DND:
            add_tip(self, "Drag & drop is not available (install 'tkinterdnd2' to "
                          "enable it). You can still click the card to select the file.")
            return
        try:
            self.drop_target_register(DND_FILES)
            self.dnd_bind("<<Drop>>", self._on_drop)
            self.dnd_bind("<<DragEnter>>", lambda e: self._set_bg(self._bg_drop))
            self.dnd_bind("<<DragLeave>>", lambda e: self._on_leave())
        except Exception:
            pass

    def _trigger_pick(self):
        self.pick_command()

    def _on_enter(self, _event=None):
        self._set_bg(self._bg_hover)

    def _on_leave(self, _event=None):
        self._set_bg(self._bg_idle)

    def _set_bg(self, color):
        self.configure(bg=color)
        for w in self._bg_widgets[1:]:
            w.configure(bg=color)

    def _on_drop(self, event):
        try:
            paths = self.tk.splitlist(event.data)
        except Exception:
            paths = [event.data]
        if paths:
            path = paths[0].strip("{}")
            self.path_var.set(path)
        self._on_leave()

    def _refresh(self):
        path = self.path_var.get().strip()
        if not path:
            self.file_lbl.configure(text="No file selected — drag it here or click",
                                     fg=PALETTE["muted"])
            self.icon_lbl.configure(fg=PALETTE["accent"])
            self.configure(highlightbackground=PALETTE["border_strong"])
        elif os.path.isfile(path):
            self.file_lbl.configure(text=f"✓  {os.path.basename(path)}", fg=PALETTE["ok"])
            self.icon_lbl.configure(fg=PALETTE["ok"])
            self.configure(highlightbackground=PALETTE["ok"])
        else:
            self.file_lbl.configure(text=f"✗  File not found: {os.path.basename(path)}",
                                     fg=PALETTE["err"])
            self.icon_lbl.configure(fg=PALETTE["err"])
            self.configure(highlightbackground=PALETTE["err"])


# =============================================================================
# ScrollableFrame: container with a vertical scrollbar always on the RIGHT.
# Used to wrap the content of EVERY tab in the program, so that no control
# is ever cut off by the window edge.
# =============================================================================
class ScrollableFrame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.canvas = tk.Canvas(self, borderwidth=0, highlightthickness=0, bg=PALETTE["bg"])
        self.vbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)

        self.canvas.configure(yscrollcommand=self.vbar.set)
        # The scrollbar always stays on the right, the canvas takes the rest.
        self.vbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self._window_id = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")

        self.inner.bind("<Configure>", self._on_inner_configure)
        self.canvas.bind("<Configure>", self._on_canvas_configure)
        # Mouse wheel scrolling only works while the mouse is over THIS tab:
        # it binds/unbinds bind_all on every Enter/Leave, so the other
        # (hidden) tabs don't intercept the event.
        self.canvas.bind("<Enter>", self._bind_wheel)
        self.canvas.bind("<Leave>", self._unbind_wheel)

    def _on_inner_configure(self, _event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def _on_canvas_configure(self, event):
        self.canvas.itemconfig(self._window_id, width=event.width)

    def _bind_wheel(self, _event):
        self.canvas.bind_all("<MouseWheel>", self._on_wheel)
        self.canvas.bind_all("<Button-4>", self._on_wheel)
        self.canvas.bind_all("<Button-5>", self._on_wheel)

    def _unbind_wheel(self, _event):
        self.canvas.unbind_all("<MouseWheel>")
        self.canvas.unbind_all("<Button-4>")
        self.canvas.unbind_all("<Button-5>")

    def _on_wheel(self, event):
        if event.num == 5 or event.delta < 0:
            self.canvas.yview_scroll(1, "units")
        elif event.num == 4 or event.delta > 0:
            self.canvas.yview_scroll(-1, "units")


def confidence_ellipse(x, y, ax, n_std=1.5, **kwargs):
    """Draws a confidence ellipse (approximated at n_std standard
    deviations) around the centroid of the points (x, y). Standard
    technique for highlighting sample clustering in a PCA chart, widely
    used in publication figures. Requires at least 3 points."""
    if len(x) < 3:
        return None
    cov = np.cov(x, y)
    if cov[0, 0] <= 0 or cov[1, 1] <= 0:
        return None
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    pearson = max(min(pearson, 0.999), -0.999)
    radius_x = np.sqrt(1 + pearson)
    radius_y = np.sqrt(1 - pearson)
    ellipse = mpatches.Ellipse((0, 0), width=radius_x * 2, height=radius_y * 2, **kwargs)

    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_x, mean_y = np.mean(x), np.mean(y)

    transf = (transforms.Affine2D()
              .rotate_deg(45)
              .scale(scale_x, scale_y)
              .translate(mean_x, mean_y))
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


# =============================================================================
# Reusable widget: interactive "publication-ready" chart
#   - zoom/pan via the standard matplotlib toolbar
#   - hover on points (if mplcursors is available)
#   - common style controls: grid, legend + position, minimalist axes,
#     serif/sans font, panel label (e.g. "A")
#   - manual annotations: arrow / line / box / text, drawn with the mouse
#   - high-resolution export with journal format presets
# =============================================================================
class InteractivePlotPanel(ttk.Frame):
    """A panel with an embedded matplotlib chart, a style bar, an
    annotation bar to fine-tune the figure by hand, and a text box
    explaining how to read the chart."""

    def __init__(self, parent, explanation: str, figsize=(6.4, 5.0), redraw_callback=None):
        super().__init__(parent)
        self.redraw_callback = redraw_callback or (lambda: None)

        self.figure = Figure(figsize=figsize, dpi=100)
        self.figure.patch.set_facecolor("white")
        self.ax = self.figure.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar_frame = ttk.Frame(self)
        toolbar_frame.pack(side=tk.TOP, fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame, pack_toolbar=False)
        self.toolbar.update()
        bg_color = PALETTE["accent_pale_2"]  # Same light color as the box at the bottom
        button_bg = "#cae2f9"                # Button background
        toolbar_frame.config(style="TFrame") # Make sure the standard ttk style is used
        self.toolbar.config(background=bg_color)
        if hasattr(self.toolbar, '_message_label'):
            self.toolbar._message_label.config(
                background=bg_color, 
                foreground=PALETTE["text"], 
                font=FONT_SMALL
            )

        for child in self.toolbar.winfo_children():
            if isinstance(child, tk.Button):
                child.config(
                    background=button_bg,
                    activebackground="#e4e6eb",  # Slightly darker color on mouse hover
                    relief="flat",               # Removes the raised 3D border
                    borderwidth=1,
                    highlightthickness=0,        # Removes the extra focus border
                    padx=6,
                    pady=4
                )
                child.bind("<Enter>", lambda e, btn=child: btn.config(background="#e4e6eb"))
                child.bind("<Leave>", lambda e, btn=child: btn.config(background=button_bg))
                
            elif isinstance(child, (tk.Frame, tk.Label)):
                child.config(background=bg_color)

        self.toolbar.pack(side=tk.TOP, fill=tk.X)

        # --- "publication-ready" style controls state ---
        self.grid_on = tk.BooleanVar(value=True)
        self.minimal_spines = tk.BooleanVar(value=True)
        self.legend_on = tk.BooleanVar(value=True)
        self.legend_pos = tk.StringVar(value="best")
        self.font_style = tk.StringVar(value="sans")
        self.panel_label = tk.StringVar(value="")
        self.panel_label_size = tk.IntVar(value=14)

        self._build_style_toolbar()

        # --- manual annotations state ---
        self.annotation_mode = tk.StringVar(value="none")   # none|arrow|line|rect|text
        self.annotation_color = "#2c3e50"
        self.annotation_lw = tk.DoubleVar(value=1.5)
        self.custom_artists: list = []
        self._drag_start = None

        self._build_annotation_toolbar()

        box = tk.Text(self, height=4, wrap="word", bg=PALETTE["accent_pale_2"], relief=tk.FLAT,
                      font=FONT_SMALL, fg=PALETTE["text"], padx=8, pady=6)
        box.insert("1.0", explanation)
        box.configure(state="disabled")
        box.pack(side=tk.BOTTOM, fill=tk.X, padx=4, pady=(2, 6))
        self.explanation_box = box

        self.canvas.mpl_connect("button_press_event", self._on_press)
        self.canvas.mpl_connect("button_release_event", self._on_release)

    # ------------------------------------------------------------------
    def _build_style_toolbar(self):
        bar = ttk.Frame(self)
        bar.pack(side=tk.TOP, fill=tk.X, pady=(4, 0))

        ttk.Label(bar, text="Figure style (for publication):", font=FONT_SMALL).pack(
            side=tk.LEFT, padx=(8, 6))

        cb_grid = ttk.Checkbutton(bar, text="Grid", variable=self.grid_on,
                                   command=self.redraw_callback)
        cb_grid.pack(side=tk.LEFT, padx=3)
        add_tip(cb_grid, "Shows/hides the grid lines behind the chart.")

        cb_spines = ttk.Checkbutton(bar, text="Minimalist axes", variable=self.minimal_spines,
                                     command=self.redraw_callback)
        cb_spines.pack(side=tk.LEFT, padx=3)
        add_tip(cb_spines, "Hides the top and right borders of the chart, as required by "
                            "the typical style of scientific journals.")

        cb_legend = ttk.Checkbutton(bar, text="Legend", variable=self.legend_on,
                                     command=self.redraw_callback)
        cb_legend.pack(side=tk.LEFT, padx=(10, 3))
        add_tip(cb_legend, "Shows/hides the chart legend.")

        pos_combo = ttk.Combobox(bar, textvariable=self.legend_pos, state="readonly", width=12,
                                  values=LEGEND_POSITIONS)
        pos_combo.pack(side=tk.LEFT, padx=3)
        pos_combo.bind("<<ComboboxSelected>>", lambda e: self.redraw_callback())
        add_tip(pos_combo, "Legend position. 'outside' places it outside the chart area, "
                            "useful when it overlaps data points.")

        ttk.Label(bar, text="Font:", font=FONT_SMALL).pack(side=tk.LEFT, padx=(10, 2))
        font_combo = ttk.Combobox(bar, textvariable=self.font_style, state="readonly", width=8,
                                   values=["sans", "serif"])
        font_combo.pack(side=tk.LEFT, padx=2)
        font_combo.bind("<<ComboboxSelected>>", lambda e: self.redraw_callback())
        add_tip(font_combo, "Many journals require a serif font (e.g. Times) for figures; "
                             "'sans' is the most common style for presentations and posters.")

        ttk.Label(bar, text="Panel label:", font=FONT_SMALL).pack(side=tk.LEFT, padx=(10, 2))
        label_entry = ttk.Entry(bar, textvariable=self.panel_label, width=4)
        label_entry.pack(side=tk.LEFT, padx=2)
        label_entry.bind("<KeyRelease>", lambda e: self.redraw_callback())
        add_tip(label_entry, "E.g. 'A', 'B', 'C': adds a bold letter in the top-left corner, "
                              "useful for multi-panel figures in a paper.")

        ttk.Spinbox(bar, from_=8, to=28, textvariable=self.panel_label_size, width=4,
                    command=self.redraw_callback).pack(side=tk.LEFT, padx=2)

    def apply_common_style(self):
        """Applies grid/minimalist axes/font/panel label. Must be called
        inside every drawing function, after plotting the data but before
        panel.redraw()."""
        ax = self.ax
        ax.grid(self.grid_on.get(), alpha=0.25)
        show_side_spines = not self.minimal_spines.get()
        for side in ("top", "right"):
            ax.spines[side].set_visible(show_side_spines)

        family = "serif" if self.font_style.get() == "serif" else "sans-serif"
        text_objs = [ax.title, ax.xaxis.label, ax.yaxis.label]
        text_objs += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        legend = ax.get_legend()
        if legend:
            text_objs += list(legend.get_texts())
            if legend.get_title():
                text_objs.append(legend.get_title())
        for obj in text_objs:
            try:
                obj.set_family(family)
            except Exception:
                pass

        label = self.panel_label.get().strip()
        if label:
            ax.text(-0.12, 1.06, label, transform=ax.transAxes,
                     fontsize=self.panel_label_size.get(), fontweight="bold",
                     va="bottom", ha="left")

    def place_legend(self, **kwargs):
        """Shows/hides/positions the legend according to the style
        controls. Call this instead of ax.legend(...) directly in drawing
        functions, to respect the user's choice."""
        if not self.legend_on.get():
            leg = self.ax.get_legend()
            if leg:
                leg.remove()
            return
        pos = self.legend_pos.get()
        if pos == "outside":
            self.ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0,
                            frameon=False, **kwargs)
        else:
            self.ax.legend(loc=pos, frameon=False, **kwargs)

    # ------------------------------------------------------------------
    def _build_annotation_toolbar(self):
        bar = ttk.Frame(self)
        bar.pack(side=tk.TOP, fill=tk.X, pady=(2, 0))

        ttk.Label(bar, text="Manual annotations:", font=FONT_SMALL).pack(side=tk.LEFT, padx=(8, 4))
        modes = [("None", "none"), ("Arrow", "arrow"), ("Line", "line"),
                 ("Box", "rect"), ("Text", "text")]
        for label, val in modes:
            rb = ttk.Radiobutton(bar, text=label, value=val, variable=self.annotation_mode,
                                  style="Toolbutton")
            rb.pack(side=tk.LEFT, padx=1)
        add_tip(bar, "Choose a tool and then drag (or click, for text) directly on the "
                     "chart to add the chosen annotation.")

        ttk.Label(bar, text="Color:", font=FONT_SMALL).pack(side=tk.LEFT, padx=(10, 2))
        self.color_swatch = ColorSwatch(bar, width=2, bg=self.annotation_color, relief="solid",
                                       bd=1, command=self._pick_color, cursor="hand2")
        self.color_swatch.pack(side=tk.LEFT, padx=2)

        ttk.Label(bar, text="Width:", font=FONT_SMALL).pack(side=tk.LEFT, padx=(10, 2))
        ttk.Spinbox(bar, from_=0.5, to=6.0, increment=0.5, textvariable=self.annotation_lw,
                    width=4).pack(side=tk.LEFT)

        btn_undo = ttk.Button(bar, text="Undo last", command=self._undo_annotation)
        btn_undo.pack(side=tk.LEFT, padx=(10, 2))
        add_tip(btn_undo, "Removes the last manual annotation added.")
        btn_clear = ttk.Button(bar, text="Clear all", command=self._clear_annotations)
        btn_clear.pack(side=tk.LEFT, padx=2)
        add_tip(btn_clear, "Removes all manual annotations from the current chart.")

        btn_export = ttk.Button(bar, text="Export image",
                                 style="Accent.TButton", command=self._open_export_dialog)
        btn_export.pack(side=tk.RIGHT, padx=8)
        add_tip(btn_export, "Saves this chart in high resolution, with journal format "
                             "presets, ready to use directly in a paper or poster.")

    # ------------------------------------------------------------------
    def clear(self):
        self.ax.clear()
        self.custom_artists = []
        self.canvas.draw_idle()

    def redraw(self):
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def enable_hover(self, artist, formatter):
        """Adds a hover tooltip on points, if mplcursors is available.
        formatter(index) -> string to display."""
        if not HAS_MPLCURSORS:
            return
        cursor = mplcursors.cursor(artist, hover=True)

        @cursor.connect("add")
        def _(sel):
            sel.annotation.set_text(formatter(sel.index))
            sel.annotation.get_bbox_patch().set(fc="lightyellow", alpha=0.95)

    # ------------------------------------------------------------------
    # Manual annotations (arrows, lines, boxes, text) drawn with the mouse
    # ------------------------------------------------------------------
    def _pick_color(self):
        color = colorchooser.askcolor(color=self.annotation_color,
                                       title="Annotation color")[1]
        if color:
            self.annotation_color = color
            self.color_swatch.configure(bg=color)

    def _on_press(self, event):
        if event.inaxes != self.ax or self.annotation_mode.get() == "none":
            return
        if getattr(self.toolbar, "mode", ""):
            # the zoom/pan toolbar is active: don't draw annotations
            return
        if self.annotation_mode.get() == "text":
            self._add_text_annotation(event.xdata, event.ydata)
            return
        self._drag_start = (event.xdata, event.ydata)

    def _on_release(self, event):
        if self._drag_start is None:
            return
        x0, y0 = self._drag_start
        self._drag_start = None
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return
        x1, y1 = event.xdata, event.ydata
        if x0 == x1 and y0 == y1:
            return

        mode = self.annotation_mode.get()
        color = self.annotation_color
        lw = self.annotation_lw.get()
        artist = None

        if mode == "arrow":
            artist = self.ax.annotate(
                "", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=14),
            )
        elif mode == "line":
            (artist,) = self.ax.plot([x0, x1], [y0, y1], color=color, lw=lw,
                                      solid_capstyle="round")
        elif mode == "rect":
            x_min, y_min = min(x0, x1), min(y0, y1)
            w, h = abs(x1 - x0), abs(y1 - y0)
            artist = mpatches.Rectangle((x_min, y_min), w, h, fill=False,
                                         edgecolor=color, linewidth=lw)
            self.ax.add_patch(artist)

        if artist is not None:
            self.custom_artists.append(artist)
            self.redraw()

    def _add_text_annotation(self, x, y):
        text = simpledialog.askstring("Add text", "Annotation text:", parent=self)
        if text:
            artist = self.ax.text(
                x, y, text, fontsize=10, color=self.annotation_color, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=self.annotation_color, alpha=0.9),
            )
            self.custom_artists.append(artist)
            self.redraw()

    def _undo_annotation(self):
        if not self.custom_artists:
            return
        artist = self.custom_artists.pop()
        try:
            artist.remove()
        except Exception:
            pass
        self.redraw()

    def _clear_annotations(self):
        for artist in self.custom_artists:
            try:
                artist.remove()
            except Exception:
                pass
        self.custom_artists = []
        self.redraw()

    # ------------------------------------------------------------------
    # High-resolution export, with journal format presets
    # ------------------------------------------------------------------
    def _open_export_dialog(self):
        dlg = tk.Toplevel(self)
        dlg.title("Export image")
        dlg.configure(bg=PALETTE["bg_card"])
        dlg.resizable(False, False)
        dlg.transient(self.winfo_toplevel())
        dlg.grab_set()

        frm = ttk.Frame(dlg, padding=16)
        frm.pack(fill=tk.BOTH, expand=True)

        cur_w, cur_h = self.figure.get_size_inches()
        aspect = float(cur_w) / float(cur_h) if cur_h else 1.0

        ttk.Label(frm, text="Size preset (scientific journals):").grid(
            row=0, column=0, sticky="w", pady=5)
        preset_var = tk.StringVar(value="Custom")
        preset_combo = ttk.Combobox(frm, textvariable=preset_var, state="readonly", width=26,
                                     values=list(EXPORT_PRESETS.keys()))
        preset_combo.grid(row=0, column=1, sticky="w", pady=5)
        add_tip(preset_combo, "Sets standard width/height for single, 1.5 or double column, "
                               "typical of scientific journals.")

        ttk.Label(frm, text="File format:").grid(row=1, column=0, sticky="w", pady=5)
        fmt_var = tk.StringVar(value="PNG")
        ttk.Combobox(frm, textvariable=fmt_var, state="readonly", width=10,
                     values=["PNG", "PDF", "SVG", "TIFF", "EPS"]).grid(row=1, column=1, sticky="w", pady=5)
        ttk.Label(frm, text="PDF/SVG/EPS are vector formats: perfect quality at\nany size, "
                            "ideal for figures with text/lines.", style="Muted.TLabel",
                  justify="left").grid(row=1, column=2, rowspan=1, sticky="w", padx=(10, 0))

        ttk.Label(frm, text="Resolution (DPI):").grid(row=2, column=0, sticky="w", pady=5)
        dpi_var = tk.IntVar(value=300)
        ttk.Combobox(frm, textvariable=dpi_var, width=10,
                     values=[150, 300, 600, 1200]).grid(row=2, column=1, sticky="w", pady=5)

        ttk.Label(frm, text="Width (inches):").grid(row=3, column=0, sticky="w", pady=5)
        width_var = tk.DoubleVar(value=round(float(cur_w), 2))
        ttk.Spinbox(frm, from_=2.0, to=20.0, increment=0.5, textvariable=width_var,
                    width=10).grid(row=3, column=1, sticky="w", pady=5)

        ttk.Label(frm, text="Height (inches):").grid(row=4, column=0, sticky="w", pady=5)
        height_var = tk.DoubleVar(value=round(float(cur_h), 2))
        height_spin = ttk.Spinbox(frm, from_=2.0, to=20.0, increment=0.5, textvariable=height_var, width=10)
        height_spin.grid(row=4, column=1, sticky="w", pady=5)

        lock_var = tk.BooleanVar(value=False)
        lock_cb = ttk.Checkbutton(frm, text="Lock aspect ratio (width/height)", variable=lock_var)
        lock_cb.grid(row=5, column=0, columnspan=2, sticky="w", pady=(2, 8))
        add_tip(lock_cb, "If enabled, changing the width automatically adjusts the height "
                          "to keep the current aspect ratio.")

        def on_preset_change(_event=None):
            preset = EXPORT_PRESETS.get(preset_var.get())
            if preset:
                width_var.set(round(preset[0], 2))
                height_var.set(round(preset[1], 2))
        preset_combo.bind("<<ComboboxSelected>>", on_preset_change)

        def on_width_change(*_args):
            if lock_var.get():
                try:
                    height_var.set(round(width_var.get() / aspect, 2))
                except (ZeroDivisionError, tk.TclError):
                    pass
        width_var.trace_add("write", on_width_change)

        transparent_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(frm, text="Transparent background", variable=transparent_var).grid(
            row=6, column=0, columnspan=2, sticky="w", pady=(0, 10))

        ttk.Label(frm, text="Tip: most journals require at least\n"
                            "300 DPI (often 600 for color TIFF/EPS).",
                  style="Muted.TLabel", justify="left").grid(row=7, column=0, columnspan=2, sticky="w")

        def do_export():
            ext_map = {"PNG": ".png", "PDF": ".pdf", "SVG": ".svg", "TIFF": ".tiff", "EPS": ".eps"}
            ext = ext_map[fmt_var.get()]
            path = filedialog.asksaveasfilename(defaultextension=ext,
                                                 filetypes=[(fmt_var.get(), f"*{ext}")])
            if not path:
                return
            orig_size = self.figure.get_size_inches()
            try:
                self.figure.set_size_inches(float(width_var.get()), float(height_var.get()))
                self.figure.savefig(path, dpi=int(dpi_var.get()),
                                     transparent=transparent_var.get(),
                                     bbox_inches="tight", pad_inches=0.15)
            except Exception as e:
                messagebox.showerror("Export error", str(e), parent=dlg)
                return
            finally:
                self.figure.set_size_inches(*orig_size)
                self.canvas.draw_idle()
            dlg.destroy()
            messagebox.showinfo("Export complete", f"Image saved to:\n{path}")

        btns = ttk.Frame(frm)
        btns.grid(row=8, column=0, columnspan=3, pady=(14, 0), sticky="e")
        ttk.Button(btns, text="Cancel", command=dlg.destroy).pack(side=tk.RIGHT, padx=4)
        ttk.Button(btns, text="Save...", style="Accent.TButton", command=do_export).pack(side=tk.RIGHT, padx=4)


# =============================================================================
# Main application
# =============================================================================
class DEAApp(ttk.Frame):

    LABEL_MODE_DISPLAY = [
        "No labels",
        "All significant genes",
        "Top N by padj",
        "Top N by |log2FC|",
        "Specific genes (list)",
    ]
    LABEL_MODE_MAP = {
        "No labels": "none",
        "All significant genes": "all_sig",
        "Top N by padj": "top_padj",
        "Top N by |log2FC|": "top_lfc",
        "Specific genes (list)": "custom",
    }

    REGULATION_FILTER_OPTIONS = ["All", "Upregulated only", "Downregulated only", "Not Significant only"]

    ORA_DIRECTION_DISPLAY = ["All DEGs (up + down)", "Up-regulated only", "Down-regulated only"]
    ORA_DIRECTION_MAP = {
        "All DEGs (up + down)": "all",
        "Up-regulated only": "up",
        "Down-regulated only": "down",
    }

    def __init__(self, root):
        super().__init__(root, style="TFrame")
        self.root = root
        self.pack(fill=tk.BOTH, expand=True)

        self.counts_path = tk.StringVar()
        self.metadata_path = tk.StringVar()
        self.gene_col = tk.StringVar(value="1")
        self.sample_col = tk.StringVar(value="1")
        self.condition_col = tk.StringVar(value="3")
        self.method = tk.StringVar(value="RNAseq")
        self.pairwise_all = tk.BooleanVar(value=False)
        self.padj_cutoff = tk.DoubleVar(value=0.05)
        self.lfc_cutoff = tk.DoubleVar(value=1.0)
        self.gsea_category = tk.StringVar(value="H")
        self.gsea_subcategory = tk.StringVar(value="")
        self.min_counts_var = tk.IntVar(value=0)
        self.top_n_pathways = tk.IntVar(value=15)

        # --- "publication-ready" chart controls ---
        self.pca_show_labels = tk.BooleanVar(value=False)
        self.pca_show_ellipses = tk.BooleanVar(value=False)
        self.pca_colorblind = tk.BooleanVar(value=False)
        self.pca_title_var = tk.StringVar(value="")

        self.volcano_label_mode_display = tk.StringVar(value=self.LABEL_MODE_DISPLAY[0])
        self.volcano_label_topn = tk.IntVar(value=15)
        self.volcano_label_fontsize = tk.IntVar(value=8)
        self.volcano_label_custom = tk.StringVar(value="")
        self.volcano_leader_lines = tk.BooleanVar(value=True)
        self.volcano_color_up = tk.StringVar(value=PALETTE["up"])
        self.volcano_color_down = tk.StringVar(value=PALETTE["down"])
        self.volcano_color_ns = tk.StringVar(value=PALETTE["ns"])
        self.volcano_point_size = tk.DoubleVar(value=18)
        self.volcano_point_alpha = tk.DoubleVar(value=0.7)
        self.volcano_title_var = tk.StringVar(value="")

        self.gsea_title_var = tk.StringVar(value="")
        self.gsea_show_asterisks = tk.BooleanVar(value=True)

        self.run_ora = tk.BooleanVar(value=True)
        self.ora_direction_display = tk.StringVar(value=self.ORA_DIRECTION_DISPLAY[0])
        self.top_n_ora = tk.IntVar(value=15)
        self.ora_chart_type = tk.StringVar(value="Bar")
        self.ora_title_var = tk.StringVar(value="")
        self._ora_colorbar = None

        self.table_regulation_filter = tk.StringVar(value=self.REGULATION_FILTER_OPTIONS[0])

        self.available_groups: list[str] = []
        self.group_vars: dict[str, tk.BooleanVar] = {}
        self.results: engine.RunResults | None = None
        self.active_tag: tk.StringVar = tk.StringVar(value="")

        self._run_start_time: float | None = None
        self._timer_job = None

        self._build_layout()
        self._check_environment()
        self._validate_setup()

    # ------------------------------------------------------------------
    def _check_environment(self):
        problems = engine.check_environment()
        if problems:
            messagebox.showwarning(
                "Environment check",
                "Some issues were detected:\n\n" + "\n".join(f"- {p}" for p in problems)
            )

    # ------------------------------------------------------------------
    def _show_help(self, tab_key: str):
        """Opens a window showing the explanation and usage instructions
        for the given tab (from HELP_CONTENT)."""
        content = HELP_CONTENT.get(tab_key)
        if not content:
            return

        dlg = tk.Toplevel(self)
        dlg.title(content["title"])
        dlg.geometry("640x560")
        dlg.configure(bg=PALETTE["bg_card"])
        dlg.transient(self.winfo_toplevel())

        header = tk.Frame(dlg, bg=PALETTE["accent"])
        header.pack(side=tk.TOP, fill=tk.X)
        tk.Label(header, text=content["title"], font=FONT_HEADER, fg="white",
                 bg=PALETTE["accent"]).pack(anchor="w", padx=16, pady=12)

        notebook = ttk.Notebook(dlg)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        def make_text_tab(text: str) -> ttk.Frame:
            tab = ttk.Frame(notebook)
            txt = scrolledtext.ScrolledText(tab, wrap="word", font=FONT_BASE,
                                             bg="white", fg=PALETTE["text"], padx=10, pady=10)
            txt.insert("1.0", text)
            txt.configure(state="disabled")
            txt.pack(fill=tk.BOTH, expand=True)
            return tab

        notebook.add(make_text_tab(content["meaning"]), text="What it means")
        notebook.add(make_text_tab(content["usage"]), text="How to use it")

        btns = ttk.Frame(dlg)
        btns.pack(fill=tk.X, padx=10, pady=(0, 10))
        ttk.Button(btns, text="Close", style="Accent.TButton", command=dlg.destroy).pack(side=tk.RIGHT)

    def _build_help_button(self, parent, tab_key: str) -> ttk.Button:
        """Creates a '? Help' button that opens the help window for the
        given tab, ready to be packed/gridded wherever needed in the
        caller's layout."""
        btn = ttk.Button(parent, text="? Help", style="Help.TButton",
                         command=lambda: self._show_help(tab_key))
        add_tip(btn, "Opens the explanation of what this tab shows and how to use its controls.")
        return btn

    # ------------------------------------------------------------------
    def _build_layout(self):
        self._build_header()

        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=(8, 0))

        self.tab_setup = ttk.Frame(self.notebook)
        self.tab_pca = ttk.Frame(self.notebook)
        self.tab_volcano = ttk.Frame(self.notebook)
        self.tab_table = ttk.Frame(self.notebook)
        self.tab_gsea = ttk.Frame(self.notebook)
        self.tab_ora = ttk.Frame(self.notebook)
        self.tab_log = ttk.Frame(self.notebook)

        self.notebook.add(self.tab_setup, text="1 · Setup & Parameters")
        self.notebook.add(self.tab_pca, text="2 · PCA")
        self.notebook.add(self.tab_volcano, text="3 · Volcano Plot")
        self.notebook.add(self.tab_table, text="4 · Gene Table")
        self.notebook.add(self.tab_gsea, text="5 · GSEA / Pathway")
        self.notebook.add(self.tab_ora, text="6 · ORA / Over-representation")
        self.notebook.add(self.tab_log, text="Log")

        # Every tab is wrapped in a ScrollableFrame: the scrollbar is
        # always on the right, in ALL tabs.
        self._build_setup_tab(self._make_scrollable(self.tab_setup))
        self._build_pca_tab(self._make_scrollable(self.tab_pca))
        self._build_volcano_tab(self._make_scrollable(self.tab_volcano))
        self._build_table_tab(self._make_scrollable(self.tab_table))
        self._build_gsea_tab(self._make_scrollable(self.tab_gsea))
        self._build_ora_tab(self._make_scrollable(self.tab_ora))
        self._build_log_tab(self._make_scrollable(self.tab_log))

        # Status bar
        status_frame = tk.Frame(self, bg=PALETTE["accent_pale_2"])
        status_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.status_var = tk.StringVar(value="Ready.")
        ttk.Label(status_frame, textvariable=self.status_var, anchor="w",
                  background=PALETTE["accent_pale_2"], foreground=PALETTE["accent"],
                  font=FONT_SMALL, padding=(10, 5)).pack(fill=tk.X)

    def _make_scrollable(self, tab: ttk.Frame) -> ttk.Frame:
        """Wraps a tab's content in a ScrollableFrame (vertical scrollbar
        on the right) and returns the inner frame in which to build the
        controls."""
        sf = ScrollableFrame(tab)
        sf.pack(fill=tk.BOTH, expand=True)
        return sf.inner

    def _build_header(self):
        header = tk.Frame(self, bg=PALETTE["accent"])
        header.pack(side=tk.TOP, fill=tk.X)
        try:
            logo_icona = tk.PhotoImage(file=engine.resource_path("logo.png"))
            logo_piccolo = logo_icona.subsample(3, 3)
        except tk.TclError:
            logo_piccolo = None

        inner = tk.Frame(header, bg=PALETTE["accent"])
        inner.pack(fill=tk.X, padx=24, pady=(16, 14))

        if logo_piccolo is not None:
            logo_header = tk.Label(inner, image=logo_piccolo, bg=PALETTE["accent"])
            logo_header.image = logo_piccolo  # Keeps Python from garbage-collecting the image
            logo_header.pack(side=tk.RIGHT, padx=15, pady=5)

        title_row = tk.Frame(inner, bg=PALETTE["accent"])
        title_row.pack(fill="x", anchor="w")
        tk.Label(title_row, text="DEA Explorer", font=("Trebuchet MS", 21, "bold"),
                 fg="white", bg=PALETTE["accent"]).pack(side=tk.LEFT)
        badge = tk.Label(title_row, text=" BETA ", font=("Trebuchet MS", 9, "bold"),
                          fg=PALETTE["accent"], bg=PALETTE["accent_soft"], padx=2)
        badge.pack(side=tk.LEFT, padx=(10, 0), pady=(4, 0))

        tk.Label(inner, text="Differential Expression & Gene Set Enrichment Analysis",
                 font=("Monotype Corsiva", 20), fg=PALETTE["accent_pale"],
                 bg=PALETTE["accent"]).pack(fill="x", anchor="w", pady=(3, 0))

        accent_line = tk.Frame(self, bg=PALETTE["accent_soft"], height=3)
        accent_line.pack(side=tk.TOP, fill=tk.X)

    # ------------------------------------------------------------------
    def _build_setup_tab(self, f):
        pad = {"padx": 8, "pady": 6}

        # --- Header row with the Help button ---
        top_bar = ttk.Frame(f)
        top_bar.pack(fill=tk.X, padx=10, pady=(10, 0))
        self._build_help_button(top_bar, "setup").pack(anchor="center", pady=(10, 5))

        # --- Preparation status: quick checklist of what's missing ---
        status_frame = ttk.LabelFrame(f, text="Preparation status")
        status_frame.pack(fill=tk.X, padx=10, pady=(10, 8))
        self.chk_counts_lbl = ttk.Label(status_frame, text="○ Counts file not selected",
                                         style="CardMuted.TLabel")
        self.chk_counts_lbl.grid(row=0, column=0, sticky="w", padx=10, pady=3)
        self.chk_metadata_lbl = ttk.Label(status_frame, text="○ Metadata file not selected",
                                           style="CardMuted.TLabel")
        self.chk_metadata_lbl.grid(row=1, column=0, sticky="w", padx=10, pady=3)
        self.chk_groups_lbl = ttk.Label(status_frame, text="○ Groups not yet extracted",
                                         style="CardMuted.TLabel")
        self.chk_groups_lbl.grid(row=2, column=0, sticky="w", padx=10, pady=3)

        # --- File input: drag-and-drop cards + column mapping ---
        file_frame = ttk.LabelFrame(f, text="Input files")
        file_frame.pack(fill=tk.X, padx=10, pady=8)

        dropzones_row = ttk.Frame(file_frame)
        dropzones_row.pack(fill=tk.X, padx=10, pady=(12, 4))
        dropzones_row.columnconfigure(0, weight=1)
        dropzones_row.columnconfigure(1, weight=1)

        self.counts_dropzone = DropZone(
            dropzones_row, title="Counts file", icon="🧬",
            subtitle="csv · tsv · txt · xlsx · xls · json",
            path_var=self.counts_path, pick_command=self._pick_counts)
        self.counts_dropzone.grid(row=0, column=0, sticky="nsew", padx=(0, 8))

        self.metadata_dropzone = DropZone(
            dropzones_row, title="Metadata file", icon="🗂️",
            subtitle="csv · tsv · txt · xlsx · xls · json",
            path_var=self.metadata_path, pick_command=self._pick_metadata)
        self.metadata_dropzone.grid(row=0, column=1, sticky="nsew", padx=(8, 0))

        if not HAS_DND:
            ttk.Label(file_frame,
                      text="Tip: install 'tkinterdnd2' (pip install tkinterdnd2) to "
                           "be able to drag files directly into the boxes above.",
                      style="Muted.TLabel", wraplength=800).pack(fill=tk.X, padx=12, pady=(2, 8))

        # Hidden entries kept for backward-compatible validation styling
        # (the visible file selection now happens via the DropZone cards).
        self.counts_entry = None
        self.metadata_entry = None

        cols_frame = ttk.Frame(file_frame)
        cols_frame.pack(fill=tk.X, padx=10, pady=(6, 10))

        gc_entry = ttk.Entry(cols_frame, textvariable=self.gene_col, width=10)
        ttk.Label(cols_frame, text="Gene column (index or name):").grid(row=0, column=0, sticky="w", **pad)
        gc_entry.grid(row=0, column=1, sticky="w", **pad)
        add_tip(gc_entry, "Column number (1 = first column) or exact name of the column "
                          "that contains the gene identifiers in the counts file.")

        sc_entry = ttk.Entry(cols_frame, textvariable=self.sample_col, width=10)
        ttk.Label(cols_frame, text="Metadata sample column (index or name):").grid(
            row=1, column=0, sticky="w", **pad)
        sc_entry.grid(row=1, column=1, sticky="w", **pad)
        add_tip(sc_entry, "Column of the metadata file containing the sample IDs (must "
                          "match the columns of the counts file).")

        cc_entry = ttk.Entry(cols_frame, textvariable=self.condition_col, width=10)
        ttk.Label(cols_frame, text="Metadata condition column (index or name):").grid(
            row=2, column=0, sticky="w", **pad)
        cc_entry.grid(row=2, column=1, sticky="w", **pad)
        add_tip(cc_entry, "Column of the metadata file containing the experimental "
                          "group/condition of each sample (e.g. 'Drug' / 'Control').")

        reset_btn = ttk.Button(cols_frame, text="↺ Reset fields", command=self._reset_setup_fields)
        reset_btn.grid(row=0, column=2, rowspan=3, padx=(30, 8), pady=6, sticky="ns")
        add_tip(reset_btn, "Clears all the fields in this section so you can start over.")

        # --- Statistical method ---
        method_frame = ttk.LabelFrame(f, text="Statistical method (chosen by the user)")
        method_frame.pack(fill=tk.X, padx=10, pady=8)
        rb1 = ttk.Radiobutton(method_frame, text="RNA-seq -> DESeq2 (Wald test)",
                               variable=self.method, value="RNAseq")
        rb1.pack(side=tk.LEFT, padx=12, pady=8)
        add_tip(rb1, "Choose this option if your data are raw RNA-seq counts (integers, "
                     "not normalized). The analysis will use DESeq2.")
        rb2 = ttk.Radiobutton(method_frame, text="Microarray -> limma (eBayes)",
                               variable=self.method, value="Microarray")
        rb2.pack(side=tk.LEFT, padx=12, pady=8)
        add_tip(rb2, "Choose this option if your data are microarray intensity values "
                     "(typically already log-transformed). The analysis will use limma.")

        # --- Groups ---
        group_frame = ttk.LabelFrame(f, text="Groups to compare")
        group_frame.pack(fill=tk.X, padx=10, pady=8)
        extract_btn = ttk.Button(group_frame, text="Extract groups from metadata",
                                  command=self._extract_groups)
        extract_btn.grid(row=0, column=0, **pad)
        add_tip(extract_btn, "Reads the metadata file and shows below all the conditions "
                              "found, with the number of samples per condition.")
        self.groups_container = ttk.Frame(group_frame)
        self.groups_container.grid(row=1, column=0, columnspan=3, sticky="w", **pad)
        pairwise_cb = ttk.Checkbutton(group_frame, text="Compare ALL pairs of groups (pairwise)",
                                      variable=self.pairwise_all, command=self._toggle_pairwise)
        pairwise_cb.grid(row=2, column=0, sticky="w", **pad)
        add_tip(pairwise_cb, "If enabled, automatically analyzes every possible pair of "
                              "groups instead of one specific two-group comparison.")

        # --- Analysis parameters ---
        param_frame = ttk.LabelFrame(f, text="Analysis parameters (customizable)")
        param_frame.pack(fill=tk.X, padx=10, pady=8)
# --- MINIMUM COUNTS SLIDER ---

        row_idx = 3  

        # Descriptive label (column 0)
        counts_lbl = ttk.Label(param_frame, text="Minimum counts filter:")
        counts_lbl.grid(row=row_idx, column=0, sticky="w", padx=5, pady=6)
        add_tip(counts_lbl, "Removes genes whose total read/count sum across all samples does not "
                            "reach this value, improving the statistical power of the test.")

        # Inner sub-frame (column 1)
        slider_subframe = ttk.Frame(param_frame)
        slider_subframe.grid(row=row_idx, column=1, sticky="ew", padx=5, pady=6)

        # Dynamic numeric label (right-aligned in the sub-frame)
        val_lbl = ttk.Label(slider_subframe, text=f"{self.min_counts_var.get()}", font=FONT_BOLD, width=4, anchor="e")
        val_lbl.pack(side="right", padx=(5, 0))

        # Runs AUTOMATICALLY every time self.min_counts_var changes
        def _update_label_on_write(*args):
            try:
                val_lbl.configure(text=str(self.min_counts_var.get()))
            except Exception:
                pass

        # Hook the variable to the auto-update function
        self.min_counts_var.trace_add("write", _update_label_on_write)

        # The slider itself
        # Note: no 'command' parameter, to avoid conflicts with the variable trace
        counts_slider = ttk.Scale(
            slider_subframe, 
            from_=0, 
            to=500, 
            variable=self.min_counts_var, 
            orient="horizontal"
        )
        counts_slider.pack(side="left", fill="x", expand=True)

        ttk.Label(param_frame, text="q-value (padj) cutoff:").grid(row=0, column=0, sticky="w", **pad)
        padj_spin = ttk.Spinbox(param_frame, from_=0.001, to=1.0, increment=0.005,
                                 textvariable=self.padj_cutoff, width=8)
        padj_spin.grid(row=0, column=1, sticky="w", **pad)
        add_tip(padj_spin, "Statistical significance threshold (multiple-testing corrected "
                            "p-value). Typical value: 0.05.")

        ttk.Label(param_frame, text="log2FoldChange cutoff:").grid(row=0, column=2, sticky="w", **pad)
        lfc_spin = ttk.Spinbox(param_frame, from_=0.0, to=10.0, increment=0.25,
                                textvariable=self.lfc_cutoff, width=8)
        lfc_spin.grid(row=0, column=3, sticky="w", **pad)
        add_tip(lfc_spin, "Minimum expression change threshold (in log2) for a gene to be "
                           "considered up/down-regulated. Typical value: 1 (= 2-fold).")

        ttk.Label(param_frame, text="MSigDB category (GSEA):").grid(row=1, column=0, sticky="w", **pad)
        cat_combo = ttk.Combobox(param_frame, textvariable=self.gsea_category, width=8,
                                  values=["H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"])
        cat_combo.grid(row=1, column=1, sticky="w", **pad)
        add_tip(cat_combo, "MSigDB gene set collection to use for the GSEA. 'H' (Hallmark) "
                            "is a good generic starting point.")

        ttk.Label(param_frame, text="Subcategory (optional, e.g. CP:KEGG):").grid(
            row=1, column=2, sticky="w", **pad)
        subcat_entry = ttk.Entry(param_frame, textvariable=self.gsea_subcategory, width=14)
        subcat_entry.grid(row=1, column=3, sticky="w", **pad)
        add_tip(subcat_entry, "Optional: restricts the chosen category to a specific "
                               "subcategory (e.g. CP:KEGG inside C2).")

        ora_cb = ttk.Checkbutton(param_frame, text="Run ORA (over-representation)",
                                  variable=self.run_ora)
        ora_cb.grid(row=2, column=0, columnspan=2, sticky="w", **pad)
        add_tip(ora_cb, "Also runs the Over-Representation Analysis on the significant genes "
                         "(same padj / log2FC cutoffs and same MSigDB category as the GSEA).")

        ttk.Label(param_frame, text="ORA gene list:").grid(row=2, column=2, sticky="w", **pad)
        ora_dir_combo = ttk.Combobox(param_frame, textvariable=self.ora_direction_display,
                                      values=self.ORA_DIRECTION_DISPLAY, state="readonly", width=22)
        ora_dir_combo.grid(row=2, column=3, sticky="w", **pad)
        add_tip(ora_dir_combo, "Which significant genes are tested: all DEGs (up + down together), "
                                "only the up-regulated or only the down-regulated ones.")

        # --- Execution ---
        run_frame = ttk.Frame(f)
        run_frame.pack(fill=tk.X, padx=10, pady=10)
        self.run_button = ttk.Button(run_frame, text="Run analysis", style="Accent.TButton",
                                      command=self._run_analysis, state="disabled")
        self.run_button.pack(side=tk.LEFT, padx=6, pady=6)
        add_tip(self.run_button, "Runs DESeq2/limma and GSEA with the parameters set above. "
                                  "Enabled once the files have been selected.")
        self.progress = ttk.Progressbar(run_frame, mode="indeterminate", length=250)
        self.progress.pack(side=tk.LEFT, padx=14)

        info = ("Tip: select the files, press 'Extract groups from metadata' to see the "
                "available conditions, check EXACTLY 2 groups (or enable pairwise "
                "comparison), set the parameters and press 'Run analysis'.")
        ttk.Label(f, text=info, wraplength=800, style="Muted.TLabel").pack(fill=tk.X, padx=14, pady=(2, 12))

        # Live validation of required fields
        self.counts_path.trace_add("write", lambda *a: self._validate_setup())
        self.metadata_path.trace_add("write", lambda *a: self._validate_setup())

    def _reset_setup_fields(self):
        self.counts_path.set("")
        self.metadata_path.set("")
        self.gene_col.set("1")
        self.sample_col.set("1")
        self.condition_col.set("2")
        for child in self.groups_container.winfo_children():
            child.destroy()
        self.group_vars = {}
        self.available_groups = []
        self._validate_setup()
        self.status_var.set("Fields reset.")

    def _validate_setup(self):
        counts_ok = bool(self.counts_path.get().strip())
        metadata_ok = bool(self.metadata_path.get().strip())

        if hasattr(self, "chk_counts_lbl"):
            if counts_ok and os.path.isfile(self.counts_path.get()):
                self.chk_counts_lbl.configure(text="✓ Counts file selected", style="Ok.TLabel")
            elif counts_ok:
                self.chk_counts_lbl.configure(text="✗ Counts file path not found", style="Err.TLabel")
            else:
                self.chk_counts_lbl.configure(text="○ Counts file not selected", style="CardMuted.TLabel")

        if hasattr(self, "chk_metadata_lbl"):
            if metadata_ok and os.path.isfile(self.metadata_path.get()):
                self.chk_metadata_lbl.configure(text="✓ Metadata file selected", style="Ok.TLabel")
            elif metadata_ok:
                self.chk_metadata_lbl.configure(text="✗ Metadata file path not found", style="Err.TLabel")
            else:
                self.chk_metadata_lbl.configure(text="○ Metadata file not selected", style="CardMuted.TLabel")

        if hasattr(self, "run_button"):
            self.run_button.configure(state="normal" if (counts_ok and metadata_ok) else "disabled")

    def _toggle_pairwise(self):
        state = "disabled" if self.pairwise_all.get() else "normal"
        for child in self.groups_container.winfo_children():
            child.configure(state=state)

    # ------------------------------------------------------------------
    def _build_pca_tab(self, f):
        top_bar = ttk.Frame(f)
        top_bar.pack(fill=tk.X, padx=8, pady=(8, 0))
        self._build_help_button(top_bar, "pca").pack(anchor="center", pady=(10, 5))

        top = ttk.Frame(f)
        top.pack(fill=tk.X, padx=8, pady=(4, 4))
        cb_labels = ttk.Checkbutton(top, text="Show sample labels", variable=self.pca_show_labels)
        cb_labels.pack(side=tk.LEFT)
        add_tip(cb_labels, "Writes the ID of each sample next to its point in the chart.")

        cb_ellipse = ttk.Checkbutton(top, text="Confidence ellipses per group",
                                      variable=self.pca_show_ellipses)
        cb_ellipse.pack(side=tk.LEFT, padx=(14, 0))
        add_tip(cb_ellipse, "Draws an ellipse around the samples of each group (requires "
                             "at least 3 samples per group): makes the separation between "
                             "conditions visually clearer, widely used in PCA figures in papers.")

        cb_cb = ttk.Checkbutton(top, text="Color-blind-friendly palette", variable=self.pca_colorblind)
        cb_cb.pack(side=tk.LEFT, padx=(14, 0))
        add_tip(cb_cb, "Uses the Okabe-Ito palette, readable even by people with a color "
                       "vision deficiency (recommended for published figures).")

        ttk.Label(top, text="  Custom title:").pack(side=tk.LEFT, padx=(16, 4))
        title_entry = ttk.Entry(top, textvariable=self.pca_title_var, width=32)
        title_entry.pack(side=tk.LEFT)
        add_tip(title_entry, "Overrides the chart's default title.")
        ttk.Button(top, text="Apply", style="Accent.TButton",
                   command=self._draw_pca).pack(side=tk.LEFT, padx=12)

        expl = ("PCA (Principal Component Analysis): each point is a sample. Samples "
                "with a similar color/condition that are close together indicate "
                "similar expression profiles. PC1 and PC2 show the percentage of "
                "variance explained: the higher it is, the more that component "
                "summarizes the dataset's variability. Use the toolbar's magnifying "
                "glass to zoom, and the 'Manual annotations' bar to add arrows/boxes/"
                "text before exporting the figure for publication.")
        self.pca_panel = InteractivePlotPanel(f, expl, redraw_callback=self._draw_pca)
        self.pca_panel.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))

    def _build_volcano_tab(self, f):
        top_bar = ttk.Frame(f)
        top_bar.pack(fill=tk.X, padx=8, pady=(8, 0))
        self._build_help_button(top_bar, "volcano").pack(anchor="center", pady=(10, 5))

        top = ttk.Frame(f)
        top.pack(fill=tk.X, padx=8, pady=(4, 4))
        ttk.Label(top, text="Contrast:", font=FONT_BOLD).pack(side=tk.LEFT)
        self.volcano_contrast_combo = ttk.Combobox(top, textvariable=self.active_tag,
                                                     state="readonly", width=28)
        self.volcano_contrast_combo.pack(side=tk.LEFT, padx=6)
        self.volcano_contrast_combo.bind("<<ComboboxSelected>>", lambda e: self._refresh_all_views())
        ttk.Label(top, text="(after changing the cutoffs in the Setup tab, press 'Apply "
                             "changes' below to update without rerunning R)",
                  style="Muted.TLabel").pack(side=tk.LEFT, padx=10)

        direction_bar = ttk.Frame(f)
        direction_bar.pack(fill=tk.X, padx=8, pady=(0, 4))
        self.direction_label_var = tk.StringVar(value="")
        ttk.Label(direction_bar, textvariable=self.direction_label_var,
                  font=FONT_BOLD, foreground=PALETTE["accent"]).pack(side=tk.LEFT)

        controls = ttk.Frame(f)
        controls.pack(fill=tk.X, padx=8, pady=(0, 4))
        controls.columnconfigure(0, weight=1)
        controls.columnconfigure(1, weight=1)

        # --- Gene labels ---
        label_frame = ttk.LabelFrame(controls, text="Gene labels (for publication)")
        label_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 6), pady=4)

        ttk.Label(label_frame, text="Mode:").grid(row=0, column=0, sticky="w", padx=6, pady=4)
        mode_combo = ttk.Combobox(label_frame, textvariable=self.volcano_label_mode_display, state="readonly",
                                   width=26, values=self.LABEL_MODE_DISPLAY)
        mode_combo.grid(row=0, column=1, columnspan=3, sticky="w", padx=6, pady=4)
        add_tip(mode_combo, "Choose which genes to label with their name directly on the chart.")

        ttk.Label(label_frame, text="Top N:").grid(row=1, column=0, sticky="w", padx=6, pady=4)
        ttk.Spinbox(label_frame, from_=1, to=200, textvariable=self.volcano_label_topn,
                    width=6).grid(row=1, column=1, sticky="w", padx=6, pady=4)

        ttk.Label(label_frame, text="Font size:").grid(row=1, column=2, sticky="w", padx=6, pady=4)
        ttk.Spinbox(label_frame, from_=6, to=20, textvariable=self.volcano_label_fontsize,
                    width=6).grid(row=1, column=3, sticky="w", padx=6, pady=4)

        ttk.Label(label_frame, text="Specific genes (comma-separated):").grid(
            row=2, column=0, columnspan=2, sticky="w", padx=6, pady=4)
        custom_entry = ttk.Entry(label_frame, textvariable=self.volcano_label_custom, width=34)
        custom_entry.grid(row=3, column=0, columnspan=4, sticky="we", padx=6, pady=(0, 4))
        add_tip(custom_entry, "Use this field with the 'Specific genes (list)' mode to label "
                               "only the genes of interest to you, e.g.: TP53, MYC, EGFR")

        leader_cb = ttk.Checkbutton(label_frame, text="Leader lines to label",
                                     variable=self.volcano_leader_lines)
        leader_cb.grid(row=4, column=0, columnspan=4, sticky="w", padx=6, pady=4)
        add_tip(leader_cb, "Draws a thin line connecting the point to the label, useful "
                            "when there are many labels close together.")

        # --- Appearance ---
        style_frame = ttk.LabelFrame(controls, text="Chart appearance")
        style_frame.grid(row=0, column=1, sticky="nsew", pady=4)

        ttk.Label(style_frame, text="Up color:").grid(row=0, column=0, sticky="w", padx=6, pady=4)
        self.volcano_up_swatch = ColorSwatch(
            style_frame, width=2, bg=self.volcano_color_up.get(), relief="solid", bd=1, cursor="hand2",
            command=lambda: self._pick_series_color(self.volcano_color_up, self.volcano_up_swatch))
        self.volcano_up_swatch.grid(row=0, column=1, sticky="w", padx=4)

        ttk.Label(style_frame, text="Down color:").grid(row=0, column=2, sticky="w", padx=6, pady=4)
        self.volcano_down_swatch = ColorSwatch(
            style_frame, width=2, bg=self.volcano_color_down.get(), relief="solid", bd=1, cursor="hand2",
            command=lambda: self._pick_series_color(self.volcano_color_down, self.volcano_down_swatch))
        self.volcano_down_swatch.grid(row=0, column=3, sticky="w", padx=4)

        ttk.Label(style_frame, text="NS color:").grid(row=0, column=4, sticky="w", padx=6, pady=4)
        self.volcano_ns_swatch = ColorSwatch(
            style_frame, width=2, bg=self.volcano_color_ns.get(), relief="solid", bd=1, cursor="hand2",
            command=lambda: self._pick_series_color(self.volcano_color_ns, self.volcano_ns_swatch))
        self.volcano_ns_swatch.grid(row=0, column=5, sticky="w", padx=4)

        ttk.Label(style_frame, text="Point size:").grid(row=1, column=0, sticky="w", padx=6, pady=4)
        ttk.Spinbox(style_frame, from_=4, to=80, textvariable=self.volcano_point_size,
                    width=6).grid(row=1, column=1, sticky="w", padx=4)

        ttk.Label(style_frame, text="Transparency:").grid(row=1, column=2, sticky="w", padx=6, pady=4)
        ttk.Spinbox(style_frame, from_=0.1, to=1.0, increment=0.05, textvariable=self.volcano_point_alpha,
                    width=6).grid(row=1, column=3, sticky="w", padx=4)

        ttk.Label(style_frame, text="Custom title:").grid(row=2, column=0, sticky="w", padx=6, pady=4)
        ttk.Entry(style_frame, textvariable=self.volcano_title_var, width=26).grid(
            row=2, column=1, columnspan=5, sticky="we", padx=6, pady=4)

        ttk.Button(controls, text="Apply changes", style="Accent.TButton",
                   command=self._draw_volcano).grid(row=1, column=0, columnspan=2, pady=(6, 2))

        expl = ("Volcano plot: each point is a gene. X axis = log2 Fold Change (how much "
                "and in which direction expression changes); Y axis = -log10(padj) "
                "(statistical significance, higher = more significant). Colors follow "
                "the cutoffs set in the Setup tab, but you can customize them above. "
                "Use the 'Gene labels' controls to show only the names you want in the "
                "figure, and the 'Manual annotations' bar to add arrows/boxes/lines: "
                "add them last, since they are cleared when the chart is regenerated.")
        self.volcano_panel = InteractivePlotPanel(f, expl, redraw_callback=self._draw_volcano)
        self.volcano_panel.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))

    def _pick_series_color(self, var: tk.StringVar, button: tk.Button):
        color = colorchooser.askcolor(color=var.get(), title="Choose color")[1]
        if color:
            var.set(color)
            button.configure(bg=color)

    def _build_table_tab(self, f):
        top_bar = ttk.Frame(f)
        top_bar.pack(fill=tk.X, padx=8, pady=(6, 0))
        self._build_help_button(top_bar, "table").pack(anchor="center", pady=(10, 5))

        top = ttk.Frame(f)
        top.pack(fill=tk.X, padx=8, pady=6)
        ttk.Label(top, text="Filter (search GeneID):").pack(side=tk.LEFT)
        self.table_filter_var = tk.StringVar()
        entry = ttk.Entry(top, textvariable=self.table_filter_var, width=25)
        entry.pack(side=tk.LEFT, padx=6)
        entry.bind("<KeyRelease>", lambda e: self._populate_gene_table())
        add_tip(entry, "Type to filter the table in real time by gene name.")

        reg_combo = ttk.Combobox(top, textvariable=self.table_regulation_filter, state="readonly",
                                  width=22, values=self.REGULATION_FILTER_OPTIONS)
        reg_combo.pack(side=tk.LEFT, padx=10)
        reg_combo.bind("<<ComboboxSelected>>", lambda e: self._populate_gene_table())
        add_tip(reg_combo, "Shows only genes with a specific regulation direction.")

        ttk.Button(top, text="Export CSV...", command=self._export_gene_table).pack(side=tk.LEFT, padx=10)

        self.table_count_lbl = ttk.Label(top, text="", style="Muted.TLabel")
        self.table_count_lbl.pack(side=tk.RIGHT, padx=10)

        table_container = ttk.Frame(f)
        table_container.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))

        columns = ("GeneID", "Regulation", "Significant", "log2FoldChange", "pvalue", "padj", "stat")
        self.gene_tree = ttk.Treeview(table_container, columns=columns, show="headings", height=20)
        for c in columns:
            self.gene_tree.heading(c, text=c, command=lambda cc=c: self._sort_tree(self.gene_tree, cc, False))
            self.gene_tree.column(c, width=110, anchor="center")
        self.gene_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.gene_tree.tag_configure("up", foreground=PALETTE["up"])
        self.gene_tree.tag_configure("down", foreground=PALETTE["down"])
        self.gene_tree.tag_configure("ns", foreground=PALETTE["muted"])

        vsb = ttk.Scrollbar(table_container, orient="vertical", command=self.gene_tree.yview)
        self.gene_tree.configure(yscrollcommand=vsb.set)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

    def _build_gsea_tab(self, f):
        top_bar = ttk.Frame(f)
        top_bar.pack(fill=tk.X, padx=8, pady=(6, 0))
        self._build_help_button(top_bar, "gsea").pack(anchor="center", pady=(10, 5))

        top = ttk.Frame(f)
        top.pack(fill=tk.X, padx=8, pady=6)
        ttk.Label(top, text="Top N pathways:").pack(side=tk.LEFT)
        ttk.Spinbox(top, from_=3, to=100, increment=1, textvariable=self.top_n_pathways,
                    width=6, command=self._draw_gsea).pack(side=tk.LEFT, padx=6)

        asterisk_cb = ttk.Checkbutton(top, text="Significance asterisks",
                                       variable=self.gsea_show_asterisks, command=self._draw_gsea)
        asterisk_cb.pack(side=tk.LEFT, padx=(14, 0))
        add_tip(asterisk_cb, "Adds * (padj<0.05), ** (padj<0.01), *** (padj<0.001) next to "
                              "each bar, as in the figures of many publications.")

        ttk.Label(top, text="  Custom title:").pack(side=tk.LEFT, padx=(10, 4))
        ttk.Entry(top, textvariable=self.gsea_title_var, width=26).pack(side=tk.LEFT)
        ttk.Button(top, text="Apply", style="Accent.TButton", command=self._draw_gsea).pack(
            side=tk.LEFT, padx=10)
        ttk.Button(top, text="Export pathway table CSV...",
                   command=self._export_gsea_table).pack(side=tk.LEFT, padx=10)

        paned = ttk.Panedwindow(f, orient=tk.VERTICAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))

        expl = ("GSEA (Gene Set Enrichment Analysis): shows the most modulated "
                "biological pathways. NES (Normalized Enrichment Score) positive = "
                "pathway activated in the reference group, negative = repressed. Only "
                "pathways with padj below the chosen threshold (or the top N by |NES| "
                "if none is significant) are shown. Use the 'Manual annotations' bar "
                "to highlight pathways of interest before exporting the figure.")
        self.gsea_panel = InteractivePlotPanel(paned, expl, figsize=(6.4, 4.0), redraw_callback=self._draw_gsea)
        paned.add(self.gsea_panel, weight=2)

        table_frame = ttk.Frame(paned)
        columns_g = ("pathway", "NES", "padj", "size")
        self.gsea_tree = ttk.Treeview(table_frame, columns=columns_g, show="headings", height=8)
        for c in columns_g:
            self.gsea_tree.heading(c, text=c, command=lambda cc=c: self._sort_tree(self.gsea_tree, cc, False))
            self.gsea_tree.column(c, width=120, anchor="center")
        self.gsea_tree.pack(fill=tk.BOTH, expand=True)
        paned.add(table_frame, weight=1)

    def _build_ora_tab(self, f):
        top_bar = ttk.Frame(f)
        top_bar.pack(fill=tk.X, padx=8, pady=(6, 0))
        self._build_help_button(top_bar, "ora").pack(anchor="center", pady=(10, 5))

        top = ttk.Frame(f)
        top.pack(fill=tk.X, padx=8, pady=6)
        ttk.Label(top, text="Top N pathways:").pack(side=tk.LEFT)
        ttk.Spinbox(top, from_=3, to=100, increment=1, textvariable=self.top_n_ora,
                    width=6, command=self._draw_ora).pack(side=tk.LEFT, padx=6)

        ttk.Label(top, text="  Chart type:").pack(side=tk.LEFT, padx=(10, 4))
        type_combo = ttk.Combobox(top, textvariable=self.ora_chart_type, values=["Bar", "Dot"],
                                   state="readonly", width=6)
        type_combo.pack(side=tk.LEFT)
        type_combo.bind("<<ComboboxSelected>>", lambda e: self._draw_ora())
        add_tip(type_combo, "Bar: significance as bar length (-log10 padj), with k/K, fold enrichment "
                             "and padj next to each bar. Dot: Gene Ratio on X, size = genes, color = padj.")

        ttk.Label(top, text="  Custom title:").pack(side=tk.LEFT, padx=(10, 4))
        ttk.Entry(top, textvariable=self.ora_title_var, width=26).pack(side=tk.LEFT)
        ttk.Button(top, text="Apply", style="Accent.TButton", command=self._draw_ora).pack(
            side=tk.LEFT, padx=10)
        ttk.Button(top, text="Export ORA table CSV...",
                   command=self._export_ora_table).pack(side=tk.LEFT, padx=10)

        paned = ttk.Panedwindow(f, orient=tk.VERTICAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))

        expl = ("ORA (Over-Representation Analysis): tests whether the significant genes "
                "(padj and log2FC cutoffs from the Setup tab) are over-represented in each "
                "pathway, against the background of all genes tested in the DEA. k/K = genes "
                "in the list / genes of the pathway in the universe; FE = fold enrichment. "
                "Bar colour shows the direction of the genes supporting the pathway (UP / "
                "DOWN if at least 80% of them agree, otherwise MIXED). Only pathways with "
                "padj below the chosen threshold (or the top N by p-value if none is "
                "significant) are shown.")
        self.ora_panel = InteractivePlotPanel(paned, expl, figsize=(6.4, 4.0), redraw_callback=self._draw_ora)
        paned.add(self.ora_panel, weight=2)

        table_frame = ttk.Frame(paned)
        columns_o = ("pathway", "count", "k_K", "FE", "padj", "direction", "N_UP", "N_DOWN")
        headings_o = {"pathway": "pathway", "count": "genes", "k_K": "k/K", "FE": "FE",
                      "padj": "padj", "direction": "direction", "N_UP": "N_UP", "N_DOWN": "N_DOWN"}
        self.ora_tree = ttk.Treeview(table_frame, columns=columns_o, show="headings", height=8)
        for c in columns_o:
            self.ora_tree.heading(c, text=headings_o[c],
                                   command=lambda cc=c: self._sort_tree(self.ora_tree, cc, False))
            self.ora_tree.column(c, width=120 if c != "pathway" else 260, anchor="center")
        self.ora_tree.pack(fill=tk.BOTH, expand=True)
        paned.add(table_frame, weight=1)

    def _build_log_tab(self, f):
        top_bar = ttk.Frame(f)
        top_bar.pack(fill=tk.X, padx=8, pady=(8, 0))
        self._build_help_button(top_bar, "log").pack(anchor="center", pady=(10, 5))

        top = ttk.Frame(f)
        top.pack(fill=tk.X, padx=8, pady=(4, 4))
        ttk.Button(top, text="Clear log", style="Danger.TButton",
                   command=self._clear_log).pack(side=tk.LEFT)
        ttk.Button(top, text="Save log to file...", command=self._save_log).pack(side=tk.LEFT, padx=8)

        self.log_text = scrolledtext.ScrolledText(f, wrap="word", state="disabled",
                                                    font=FONT_MONO, bg="white", fg=PALETTE["text"],
                                                    height=28)
        self.log_text.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

    def _clear_log(self):
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.configure(state="disabled")

    def _save_log(self):
        path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text", "*.txt")])
        if path:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(self.log_text.get("1.0", "end"))
            self.status_var.set(f"Log saved to: {path}")

    # ------------------------------------------------------------------
    def _log(self, msg: str):
        self.log_text.configure(state="normal")
        self.log_text.insert("end", msg + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def _pick_counts(self):
        path = filedialog.askopenfilename(title="Select the counts file",
                                           filetypes=[("Data", "*.csv *.tsv *.txt *.xlsx *.xls *.json"),
                                                      ("All files", "*.*")])
        if path:
            self.counts_path.set(path)

    def _pick_metadata(self):
        path = filedialog.askopenfilename(title="Select the metadata file",
                                           filetypes=[("Data", "*.csv *.tsv *.txt *.xlsx *.xls *.json"),
                                                      ("All files", "*.*")])
        if path:
            self.metadata_path.set(path)

    # ------------------------------------------------------------------
    def _extract_groups(self):
        if not self.metadata_path.get():
            messagebox.showerror("Error", "Please select the metadata file first.")
            return
        self.status_var.set("Extracting groups...")
        self.progress.start(10)

        def task():
            try:
                data = engine.list_groups(
                    self.metadata_path.get(),
                    sample_col=self._parse_col(self.sample_col.get()),
                    condition_col=self._parse_col(self.condition_col.get()),
                    log_callback=self._threadsafe_log,
                )
                self.root.after(0, lambda: self._on_groups_ready(data))
            except engine.RPipelineError as e:
                self.root.after(0, lambda: self._on_error("Group extraction failed", e))
            finally:
                self.root.after(0, lambda: (self.progress.stop(), self.status_var.set("Ready.")))

        threading.Thread(target=task, daemon=True).start()

    def _threadsafe_log(self, line: str):
        self.root.after(0, lambda: self._log(line))

    def _on_groups_ready(self, data: dict):
        self.available_groups = data.get("groups", [])
        for child in self.groups_container.winfo_children():
            child.destroy()
        self.group_vars = {}
        counts = data.get("counts", {})
        for i, g in enumerate(self.available_groups):
            var = tk.BooleanVar(value=False)
            n = counts.get(g, "?")
            cb = ttk.Checkbutton(self.groups_container, text=f"{g} (n={n})", variable=var)
            cb.grid(row=i // 4, column=i % 4, sticky="w", padx=8, pady=2)
            self.group_vars[g] = var
        self._log(f">>> Groups found: {', '.join(self.available_groups)}")
        self.status_var.set(f"{len(self.available_groups)} groups found.")
        if hasattr(self, "chk_groups_lbl"):
            if self.available_groups:
                self.chk_groups_lbl.configure(
                    text=f"✓ {len(self.available_groups)} groups extracted: {', '.join(self.available_groups)}",
                    style="Ok.TLabel")
            else:
                self.chk_groups_lbl.configure(text="✗ No groups found in metadata", style="Err.TLabel")

    def _on_error(self, title: str, err: Exception):
        self._log(f"[ERROR] {title}: {err}")
        messagebox.showerror(title, str(err))

    @staticmethod
    def _parse_col(value: str):
        value = value.strip()
        if value.isdigit():
            return int(value)
        return value

    # ------------------------------------------------------------------
    def _run_analysis(self):
        if not self.counts_path.get() or not self.metadata_path.get():
            messagebox.showerror("Error", "Please select both the counts file and the metadata file.")
            return

        selected = [g for g, v in self.group_vars.items() if v.get()]
        contrast = None
        if not self.pairwise_all.get():
            if len(selected) != 2:
                messagebox.showerror(
                    "Error",
                    "Select EXACTLY 2 groups to compare, or enable the pairwise "
                    "comparison of all pairs."
                )
                return
            contrast = selected

        # Lettura delle variabili Tk nel thread principale (il task gira in un thread a parte)
        run_ora = self.run_ora.get()
        ora_direction = self.ORA_DIRECTION_MAP.get(self.ora_direction_display.get(), "all")

        self.run_button.configure(state="disabled")
        self.progress.start(10)
        self._run_start_time = time.time()
        self._tick_timer()

        def task():
            try:
                results = engine.run_dea(
                    counts_path=self.counts_path.get(),
                    metadata_path=self.metadata_path.get(),
                    method=self.method.get(),
                    gene_col=self._parse_col(self.gene_col.get()),
                    sample_col=self._parse_col(self.sample_col.get()),
                    condition_col=self._parse_col(self.condition_col.get()),
                    contrast=contrast,
                    pairwise_all=self.pairwise_all.get(),
                    padj_cutoff=self.padj_cutoff.get(),
                    lfc_cutoff=self.lfc_cutoff.get(),
                    min_counts=self.min_counts_var.get(),
                    gsea_category=self.gsea_category.get(),
                    gsea_subcategory=(self.gsea_subcategory.get() or None),
                    run_ora=run_ora,
                    ora_direction=ora_direction,
                    log_callback=self._threadsafe_log,
                )
                self.root.after(0, lambda: self._on_results_ready(results))
            except engine.RPipelineError as e:
                self.root.after(0, lambda: self._on_error("Analysis failed", e))
            finally:
                self.root.after(0, self._stop_run_ui)

        threading.Thread(target=task, daemon=True).start()

    def _tick_timer(self):
        if self._run_start_time is None:
            return
        elapsed = int(time.time() - self._run_start_time)
        mm, ss = divmod(elapsed, 60)
        self.status_var.set(f"Analysis running (DESeq2/limma + ORA + GSEA)... elapsed time: {mm:02d}:{ss:02d}")
        self._timer_job = self.root.after(1000, self._tick_timer)

    def _stop_run_ui(self):
        self.progress.stop()
        self._validate_setup()
        self._run_start_time = None
        if self._timer_job:
            self.root.after_cancel(self._timer_job)
            self._timer_job = None
        self.status_var.set("Ready.")

    def _on_results_ready(self, results: engine.RunResults):
        self.results = results
        tags = list(results.per_contrast.keys())
        self.volcano_contrast_combo["values"] = tags
        if tags:
            self.active_tag.set(tags[0])
        self._log(">>> Analysis completed successfully.")
        self._refresh_all_views()
        self.notebook.select(self.tab_pca)

    # ------------------------------------------------------------------
    def _refresh_all_views(self):
        if self.results is None:
            return
        self._draw_pca()
        self._draw_volcano()
        self._populate_gene_table()
        self._draw_gsea()
        self._draw_ora()

    def _draw_pca(self):
        if self.results is None:
            return
        panel = self.pca_panel
        panel.clear()
        df = self.results.pca_data
        if "condition" not in df.columns:
            panel.redraw()
            return
        conditions = sorted(df["condition"].unique())
        colors = OKABE_ITO if self.pca_colorblind.get() else plt.cm.tab10.colors
        artists = []
        for i, cond in enumerate(conditions):
            sub = df[df["condition"] == cond]
            color = colors[i % len(colors)]
            sc = panel.ax.scatter(sub["PC1"], sub["PC2"], label=str(cond),
                                   color=color, s=80, alpha=0.85,
                                   edgecolor="black", linewidth=0.5)
            artists.append((sc, sub))
            if self.pca_show_ellipses.get():
                confidence_ellipse(sub["PC1"].values, sub["PC2"].values, panel.ax,
                                    n_std=1.5, facecolor=color, alpha=0.15,
                                    edgecolor=color, linewidth=1.2)
            if self.pca_show_labels.get() and "SampleID" in sub.columns:
                for _, r in sub.iterrows():
                    panel.ax.annotate(str(r["SampleID"]), xy=(r["PC1"], r["PC2"]),
                                       xytext=(4, 4), textcoords="offset points", fontsize=7)
        panel.ax.set_xlabel(self.results.pca_x_label)
        panel.ax.set_ylabel(self.results.pca_y_label)
        title = self.pca_title_var.get().strip() or "PCA of Count Data"
        panel.ax.set_title(title)
        panel.place_legend(title="Condition")
        panel.apply_common_style()

        if HAS_MPLCURSORS:
            all_sc = [a for a, _ in artists]
            all_sub = [s for _, s in artists]

            def fmt(idx_tuple):
                artist_idx, point_idx = idx_tuple
                sample_id = all_sub[artist_idx].iloc[point_idx].get("SampleID", "?")
                return str(sample_id)
            panel.enable_hover(all_sc, lambda idx: fmt(idx) if isinstance(idx, tuple) else str(idx))
        panel.redraw()

    def _compute_label_subset(self, df: pd.DataFrame, mode: str) -> pd.DataFrame:
        """Selects the subset of genes to label in the volcano plot,
        according to the mode chosen by the user in the 'Volcano Plot'
        tab."""
        if mode == "none" or df.empty:
            return df.iloc[0:0]
        if mode == "all_sig":
            return df[df["status"] != "Not Significant"]
        if mode == "top_padj":
            sig = df[df["status"] != "Not Significant"].copy()
            sig = sig.sort_values("padj")
            return sig.head(self.volcano_label_topn.get())
        if mode == "top_lfc":
            sig = df[df["status"] != "Not Significant"].copy()
            sig["_abs_lfc"] = sig["log2FoldChange"].abs()
            sig = sig.sort_values("_abs_lfc", ascending=False)
            return sig.head(self.volcano_label_topn.get())
        if mode == "custom":
            genes = [g.strip() for g in self.volcano_label_custom.get().split(",") if g.strip()]
            if not genes:
                return df.iloc[0:0]
            return df[df["GeneID"].astype(str).isin(genes)]
        return df.iloc[0:0]

    def _draw_volcano(self):
        if self.results is None:
            return
        panel = self.volcano_panel
        panel.clear()
        tag = self.active_tag.get()
        if not tag or tag not in self.results.per_contrast:
            panel.redraw()
            return
        res = self.results.per_contrast[tag]
        if hasattr(self, "direction_label_var"):
            if res.group_high and res.group_low:
                self.direction_label_var.set(
                    f"log2FC positivo = più espresso in \"{res.group_high}\"   |   "
                    f"log2FC negativo = più espresso in \"{res.group_low}\""
                )
            else:
                self.direction_label_var.set(
                    "Direzione del log2FC non determinabile automaticamente per questo contrasto."
                )
        df = engine.categorize_volcano(res.res_tbl, self.padj_cutoff.get(), self.lfc_cutoff.get())
        df = df.dropna(subset=["padj"]).copy()
        df["neglog10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))

        color_map = {
            "Upregulated": self.volcano_color_up.get(),
            "Downregulated": self.volcano_color_down.get(),
            "Not Significant": self.volcano_color_ns.get(),
        }
        size = self.volcano_point_size.get()
        alpha = self.volcano_point_alpha.get()
        scatters = {}
        for status, color in color_map.items():
            sub = df[df["status"] == status]
            if sub.empty:
                continue
            sc = panel.ax.scatter(sub["log2FoldChange"], sub["neglog10_padj"],
                                   c=color, s=size, alpha=alpha, label=status, edgecolor="none")
            scatters[status] = sub

        panel.ax.axvline(self.lfc_cutoff.get(), linestyle="--", color="#8a8a8a", linewidth=1)
        panel.ax.axvline(-self.lfc_cutoff.get(), linestyle="--", color="#8a8a8a", linewidth=1)
        panel.ax.axhline(-math.log10(self.padj_cutoff.get()), linestyle="--", color="#8a8a8a", linewidth=1)
        panel.ax.set_xlabel("log2 Fold Change")
        panel.ax.set_ylabel("-log10(padj)")
        title = self.volcano_title_var.get().strip() or f"Volcano Plot — {tag.replace('_', ' ')}"
        panel.ax.set_title(title)
        panel.place_legend()

        # --- gene labels, according to the chosen mode ---
        mode = self.LABEL_MODE_MAP.get(self.volcano_label_mode_display.get(), "none")
        label_df = self._compute_label_subset(df, mode)
        fontsize = self.volcano_label_fontsize.get()
        use_leader = self.volcano_leader_lines.get()
        for _, row in label_df.iterrows():
            x, y = row["log2FoldChange"], row["neglog10_padj"]
            gene = str(row.get("GeneID", ""))
            if use_leader:
                panel.ax.annotate(gene, xy=(x, y), xytext=(8, 8), textcoords="offset points",
                                   fontsize=fontsize, ha="left", va="bottom",
                                   arrowprops=dict(arrowstyle="-", lw=0.4, color="#666666"))
            else:
                panel.ax.text(x, y, gene, fontsize=fontsize, ha="left", va="bottom")

        panel.apply_common_style()

        if HAS_MPLCURSORS and scatters:
            try:
                collections = panel.ax.collections[:len(scatters)]
                sub_frames = list(scatters.values())

                def fmt(idx_tuple, frames=sub_frames):
                    artist_idx, point_idx = idx_tuple
                    return str(frames[artist_idx].iloc[point_idx]["GeneID"])
                panel.enable_hover(collections, lambda idx: fmt(idx) if isinstance(idx, tuple) else str(idx))
            except Exception:
                pass
        panel.redraw()

    def _populate_gene_table(self):
        for row in self.gene_tree.get_children():
            self.gene_tree.delete(row)
        tag = self.active_tag.get()
        if not tag or self.results is None or tag not in self.results.per_contrast:
            if hasattr(self, "table_count_lbl"):
                self.table_count_lbl.configure(text="")
            return
        ranked = self.results.per_contrast[tag].ranked_tbl
        if ranked is None or ranked.empty:
            if hasattr(self, "table_count_lbl"):
                self.table_count_lbl.configure(text="No genes to show.")
            return
        filt = self.table_filter_var.get().strip().lower()
        if filt:
            ranked = ranked[ranked["GeneID"].astype(str).str.lower().str.contains(filt)]

        reg_filter = self.table_regulation_filter.get()
        if reg_filter == "Upregulated only" and "Regolazione" in ranked.columns:
            ranked = ranked[ranked["Regolazione"].astype(str).str.contains("Up", case=False, na=False)]
        elif reg_filter == "Downregulated only" and "Regolazione" in ranked.columns:
            ranked = ranked[ranked["Regolazione"].astype(str).str.contains("Down", case=False, na=False)]
        elif reg_filter == "Not Significant only" and "Significativo" in ranked.columns:
            ranked = ranked[~ranked["Significativo"].astype(str).str.contains("Sì|Si|Yes|True", case=False, na=False)]

        for _, r in ranked.iterrows():
            reg_text = str(r.get("Regolazione", "")).lower()
            sig_text = str(r.get("Significativo", ""))
            if sig_text.strip().upper() == "SI":   # R writes SI/NO/NA: show it in English
                sig_text = "YES"
            row_tag = "up" if "up" in reg_text else ("down" if "down" in reg_text else "ns")
            self.gene_tree.insert("", "end", values=(
                r.get("GeneID", ""), r.get("Regolazione", ""), sig_text,
                round(r.get("log2FoldChange", float("nan")), 3) if pd.notna(r.get("log2FoldChange")) else "",
                _fmt_sci(r.get("pvalue")), _fmt_sci(r.get("padj")),
                round(r.get("stat", float("nan")), 3) if pd.notna(r.get("stat")) else "",
            ), tags=(row_tag,))

        if hasattr(self, "table_count_lbl"):
            self.table_count_lbl.configure(text=f"{len(ranked)} genes shown")

    def _draw_gsea(self):
        if self.results is None:
            return
        panel = self.gsea_panel
        panel.clear()
        for row in self.gsea_tree.get_children():
            self.gsea_tree.delete(row)
        tag = self.active_tag.get()
        if not tag or tag not in self.results.per_contrast:
            panel.redraw()
            return
        gsea_tbl = self.results.per_contrast[tag].gsea_tbl
        if gsea_tbl is None or gsea_tbl.empty:
            panel.ax.text(0.5, 0.5, "GSEA not available for this contrast.",
                           ha="center", va="center")
            panel.redraw()
            return

        plot_data, msg = engine.filter_top_pathways(
            gsea_tbl, top_n=self.top_n_pathways.get(), padj_cutoff=self.padj_cutoff.get())

        colors = plot_data["Direction"].map({"Up-regulated": PALETTE["up"], "Down-regulated": PALETTE["down"]})
        bars = panel.ax.barh(plot_data["pathway_clean"], plot_data["NES"], color=colors)
        panel.ax.set_xlabel("Normalized Enrichment Score (NES)")
        title = self.gsea_title_var.get().strip() or f"Top Modulated Pathways — {tag.replace('_', ' ')}"
        panel.ax.set_title(title)
        panel.ax.axvline(0, color="black", linewidth=0.8)
        panel.ax.tick_params(axis="y", labelsize=8)

        if self.gsea_show_asterisks.get():
            for bar, padj in zip(bars, plot_data["padj"]):
                stars = ""
                if pd.notna(padj):
                    if padj < 0.001:
                        stars = "***"
                    elif padj < 0.01:
                        stars = "**"
                    elif padj < 0.05:
                        stars = "*"
                if stars:
                    width = bar.get_width()
                    offset = 0.05 * (abs(width) if width != 0 else 1)
                    x_pos = width + offset if width >= 0 else width - offset
                    ha = "left" if width >= 0 else "right"
                    panel.ax.text(x_pos, bar.get_y() + bar.get_height() / 2, stars,
                                   va="center", ha=ha, fontsize=9, fontweight="bold")

        xmin, xmax = panel.ax.get_xlim()
        new_xmin = xmin*1.15 if xmin < 0 else xmin
        new_xmax = xmax*1.15 if xmax > 0 else xmax
        panel.ax.set_xlim(left=new_xmin, right=new_xmax)
        
        up_patch = mpatches.Patch(color=PALETTE["up"], label="Up-regulated")
        down_patch = mpatches.Patch(color=PALETTE["down"], label="Down-regulated")
        panel.place_legend(handles=[up_patch, down_patch])
        panel.apply_common_style()
        panel.redraw()
        self.status_var.set(msg)

        for _, r in plot_data.sort_values("padj").iterrows():
            self.gsea_tree.insert("", "end", values=(
                r.get("pathway_clean", ""), round(r.get("NES", float("nan")), 3),
                _fmt_sci(r.get("padj")), r.get("size", "")))

    def _draw_ora(self):
        if self.results is None:
            return
        panel = self.ora_panel
        # la colorbar del dot plot vive in un asse a parte: va rimossa prima di ax.clear()
        if self._ora_colorbar is not None:
            try:
                self._ora_colorbar.remove()
            except Exception:
                pass
            self._ora_colorbar = None
        panel.clear()
        for row in self.ora_tree.get_children():
            self.ora_tree.delete(row)
        tag = self.active_tag.get()
        if not tag or tag not in self.results.per_contrast:
            panel.redraw()
            return
        res = self.results.per_contrast[tag]
        ora_tbl = res.ora_tbl
        if ora_tbl is None or ora_tbl.empty:
            import textwrap
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

        cutoff = self.padj_cutoff.get()
        plot_data, msg = engine.filter_top_ora(
            ora_tbl, top_n=self.top_n_ora.get(), padj_cutoff=cutoff)
        n_in = (res.ora_meta or {}).get("n_input_genes")
        title = self.ora_title_var.get().strip() or f"Top Enriched Pathways (ORA) — {tag.replace('_', ' ')}"
        dir_colors = {"UP": PALETTE["up"], "DOWN": PALETTE["down"], "MISTO": PALETTE["ns"]}

        if self.ora_chart_type.get() == "Dot":
            from matplotlib.colors import LinearSegmentedColormap
            d = plot_data.sort_values("GeneRatioNum", kind="stable")
            cmap = LinearSegmentedColormap.from_list("ora_padj", [PALETTE["up"], PALETTE["down"]])
            sc = panel.ax.scatter(d["GeneRatioNum"], d["pathway_clean"],
                                   s=20 + 8 * d["Count"], c=d["p.adjust"], cmap=cmap)
            self._ora_colorbar = panel.figure.colorbar(sc, ax=panel.ax, pad=0.02)
            self._ora_colorbar.set_label("padj")
            panel.ax.set_xlabel("Gene Ratio")
            panel.ax.tick_params(axis="y", labelsize=8)
            panel.ax.set_title(title)
        else:
            colors = plot_data["Direzione"].map(dir_colors).fillna(PALETTE["ns"])
            bars = panel.ax.barh(plot_data["pathway_clean"], plot_data["neglog10"],
                                  color=colors, height=0.7)
            xmax = max(float(plot_data["neglog10"].max()), -np.log10(cutoff))
            if not np.isfinite(xmax) or xmax <= 0:
                xmax = 1.0
            for bar, (_, r) in zip(bars, plot_data.iterrows()):
                label = (f"{int(r['Count'])}/{int(r['SetSize'])} genes "
                         f"({int(r['N_UP'])} UP, {int(r['N_DOWN'])} DOWN) | "
                         f"FE {r['FoldEnrichment']:.1f} | padj {r['p.adjust']:.1e}")
                panel.ax.text(bar.get_width() + 0.01 * xmax, bar.get_y() + bar.get_height() / 2,
                               label, va="center", ha="left", fontsize=7)
            panel.ax.axvline(-np.log10(cutoff), linestyle="--", color=PALETTE["err"], linewidth=1)
            panel.ax.set_xlim(0, xmax * 1.9)
            panel.ax.set_xlabel(r"$-\log_{10}$(padj)")
            panel.ax.tick_params(axis="y", labelsize=8)
            panel.ax.set_title(title)
            handles = [
                mpatches.Patch(color=PALETTE["up"], label="UP (>= 80% of genes)"),
                mpatches.Patch(color=PALETTE["down"], label="DOWN (>= 80% of genes)"),
                mpatches.Patch(color=PALETTE["ns"], label="MIXED"),
            ]
            panel.place_legend(handles=handles)

        panel.apply_common_style()
        panel.redraw()
        self.status_var.set(msg + (f" Genes in input: {n_in}." if n_in else ""))

        for _, r in plot_data.sort_values("p.adjust", kind="stable").iterrows():
            self.ora_tree.insert("", "end", values=(
                r.get("pathway_clean", ""), int(r["Count"]),
                f"{int(r['Count'])}/{int(r['SetSize'])}",
                round(float(r["FoldEnrichment"]), 2), _fmt_sci(r.get("p.adjust")),
                "MIXED" if r.get("Direzione", "") == "MISTO" else r.get("Direzione", ""),
                int(r["N_UP"]), int(r["N_DOWN"])))

    # ------------------------------------------------------------------
    def _sort_tree(self, tree: ttk.Treeview, col: str, reverse: bool):
        data = [(tree.set(k, col), k) for k in tree.get_children("")]

        def key(v):
            try:
                return float(v[0])
            except (ValueError, TypeError):
                return v[0]
        data.sort(key=key, reverse=reverse)
        for index, (_, k) in enumerate(data):
            tree.move(k, "", index)
        tree.heading(col, command=lambda: self._sort_tree(tree, col, not reverse))

    def _export_gene_table(self):
        tag = self.active_tag.get()
        if not tag or self.results is None or tag not in self.results.per_contrast:
            messagebox.showinfo("Info", "No results to export.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                             filetypes=[("CSV", "*.csv")])
        if path:
            self.results.per_contrast[tag].ranked_tbl.to_csv(path, index=False)
            self.status_var.set(f"Table exported to: {path}")

    def _export_gsea_table(self):
        tag = self.active_tag.get()
        if not tag or self.results is None or tag not in self.results.per_contrast:
            messagebox.showinfo("Info", "No GSEA results to export.")
            return
        gsea_tbl = self.results.per_contrast[tag].gsea_tbl
        if gsea_tbl is None:
            messagebox.showinfo("Info", "GSEA is not available for this contrast.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if path:
            gsea_tbl.to_csv(path, index=False)
            self.status_var.set(f"GSEA table exported to: {path}")

    def _export_ora_table(self):
        tag = self.active_tag.get()
        if not tag or self.results is None or tag not in self.results.per_contrast:
            messagebox.showinfo("Info", "No ORA results to export.")
            return
        ora_tbl = self.results.per_contrast[tag].ora_tbl
        if ora_tbl is None or ora_tbl.empty:
            messagebox.showinfo("Info", "ORA is not available for this contrast.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if path:
            ora_tbl.to_csv(path, index=False)
            self.status_var.set(f"ORA table exported to: {path}")


# --- small local helpers ---
def _fmt_sci(value) -> str:
    if value is None or pd.isna(value):
        return ""
    try:
        return f"{float(value):.2e}"
    except (TypeError, ValueError):
        return str(value)


def launch():
    root = TkinterDnD.Tk() if HAS_DND else tk.Tk()
    root.title("DEA Explorer")
    try:
        logo_icona = tk.PhotoImage(file=engine.resource_path("logo.png"))
        root.iconphoto(False, logo_icona)
    except tk.TclError:
        pass
    root.geometry("1200x860")
    root.minsize(1000, 700)
    apply_theme(root)
    app = DEAApp(root)
    root.mainloop()


if __name__ == "__main__":
    launch()
