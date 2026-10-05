from __future__ import annotations

import datetime
import os
import textwrap

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle, Patch

import core_calcolo as engine

A4 = (8.27, 11.69)
C = {
    "accent": "#123c63", "accent_light": "#2f74ad", "pale": "#d7e6f2",
    "pale2": "#eaf1f8", "text": "#1b2b3a", "muted": "#6c7a89",
    "up": "#c0392b", "down": "#2471a3", "ns": "#999b9d", "grid": "#d7e3ee",
}
TAB10 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
         "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
LOG_LINES_PER_PAGE = 84
LOG_MAX_PAGES = 30


def _fmt_sci(v) -> str:
    if v is None or pd.isna(v):
        return ""
    try:
        return f"{float(v):.2e}"
    except (TypeError, ValueError):
        return str(v)


def _fmt_num(v, nd=3) -> str:
    if v is None or pd.isna(v):
        return ""
    try:
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return str(v)


def _trunc(s, n) -> str:
    s = str(s)
    return s if len(s) <= n else s[: n - 1] + "…"


def _dir_names(res):
    hi = res.group_high if getattr(res, "group_high", None) else "group 1 (log2FC>0)"
    lo = res.group_low if getattr(res, "group_low", None) else "group 2 (log2FC<0)"
    return hi, lo


def _contrast_label(res) -> str:
    if res.group_high and res.group_low:
        return f"{res.group_high} vs {res.group_low}"
    return str(res.tag).replace("_", " ")


class _Writer:
    def __init__(self, pdf: PdfPages, generated: str):
        self.pdf = pdf
        self.generated = generated
        self.n = 0

    def new(self, title: str, subtitle: str | None = None) -> Figure:
        fig = Figure(figsize=A4)
        fig.add_artist(Rectangle((0, 0.94), 1, 0.06, transform=fig.transFigure,
                                 facecolor=C["accent"], edgecolor="none"))
        fig.add_artist(Rectangle((0, 0.936), 1, 0.004, transform=fig.transFigure,
                                 facecolor=C["accent_light"], edgecolor="none"))
        fig.text(0.07, 0.97, _trunc(title, 70), color="white", fontsize=15,
                 fontweight="bold", va="center")
        if subtitle:
            fig.text(0.07, 0.915, subtitle, fontsize=8.5, color=C["muted"], va="center")
        return fig

    def save(self, fig: Figure):
        self.n += 1
        fig.text(0.07, 0.02, f"DEA Explorer  ·  generated {self.generated}",
                 fontsize=7, color=C["muted"])
        fig.text(0.93, 0.02, f"Page {self.n}", fontsize=7, color=C["muted"], ha="right")
        self.pdf.savefig(fig)


def _section(fig: Figure, text: str, y: float, x: float = 0.07) -> float:
    fig.text(x, y, text, fontsize=11, fontweight="bold", color=C["accent"], va="top")
    fig.add_artist(Rectangle((x, y - 0.022), 0.86, 0.0012, transform=fig.transFigure,
                             facecolor=C["accent_light"], edgecolor="none"))
    return y - 0.035


def _paragraph(fig: Figure, text: str, y: float, width: int = 100, fs: float = 8.5,
               color: str | None = None, x: float = 0.07, step: float = 0.0165,
               bold: bool = False) -> float:
    for block in str(text).split("\n"):
        lines = textwrap.wrap(block, width=width) or [""]
        for ln in lines:
            fig.text(x, y, ln, fontsize=fs, color=color or C["text"], va="top",
                     fontweight="bold" if bold else "normal")
            y -= step
    return y - 0.006


def _kv(fig: Figure, items: list[tuple[str, str]], y: float, x: float = 0.07,
        label_w: float = 0.24, width: int = 70, step: float = 0.0185) -> float:
    for k, v in items:
        fig.text(x, y, k, fontsize=8.5, color=C["muted"], va="top")
        lines = textwrap.wrap(str(v), width=width) or [""]
        for i, ln in enumerate(lines):
            fig.text(x + label_w, y - i * step, ln, fontsize=8.5, color=C["text"],
                     va="top", fontweight="bold")
        y -= step * len(lines) + 0.003
    return y - 0.004


def _table(fig: Figure, rows: list[list[str]], headers: list[str], widths: list[float],
           top: float, row_h: float = 0.0165, x: float = 0.07, w: float = 0.86,
           fs: float = 7) -> float:
    if not rows:
        rows = [["—"] + [""] * (len(headers) - 1)]
    n = len(rows) + 1
    h = n * row_h
    ax = fig.add_axes([x, top - h, w, h])
    ax.axis("off")
    tbl = ax.table(cellText=rows, colLabels=headers, colWidths=widths,
                   loc="center", cellLoc="left", bbox=[0, 0, 1, 1])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fs)
    for (r, _c), cell in tbl.get_celld().items():
        cell.set_edgecolor(C["grid"])
        cell.set_linewidth(0.5)
        if r == 0:
            cell.set_facecolor(C["accent"])
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 0:
            cell.set_facecolor(C["pale2"])
    return top - h - 0.015


def _style_axes(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=7)


def _note_page_text(fig: Figure, text: str, y: float = 0.5):
    fig.text(0.5, y, text, ha="center", va="center", fontsize=10, color=C["muted"])


def _contrast_stats(res, padj_cut: float, lfc_cut: float) -> dict:
    df = engine.categorize_volcano(res.res_tbl, padj_cut, lfc_cut)
    n_up = int((df["status"] == "Upregulated").sum())
    n_down = int((df["status"] == "Downregulated").sum())
    gsea_sig = ora_sig = None
    if res.gsea_tbl is not None and not res.gsea_tbl.empty and "padj" in res.gsea_tbl:
        gsea_sig = int((res.gsea_tbl["padj"] < padj_cut).sum())
    if res.ora_tbl is not None and not res.ora_tbl.empty and "p.adjust" in res.ora_tbl:
        ora_sig = int((res.ora_tbl["p.adjust"] < padj_cut).sum())
    return {
        "df": df,
        "n_rows": len(df),
        "n_tested": int(df["pvalue"].notna().sum()) if "pvalue" in df else len(df),
        "n_padj": int(df["padj"].notna().sum()),
        "n_sig": int((df["padj"] < padj_cut).sum()),
        "n_up": n_up, "n_down": n_down,
        "gsea_sig": gsea_sig, "ora_sig": ora_sig,
    }


def _page_cover(w: _Writer, results, params: dict, stats: dict, generated: str):
    fig = w.new("Differential Expression Analysis Report",
                f"Generated on {generated}")
    y = 0.885
    y = _section(fig, "Input data", y)
    pca = results.pca_data
    n_samples = len(pca) if pca is not None else 0
    conds = (sorted(pca["condition"].astype(str).unique())
             if pca is not None and "condition" in pca.columns else [])
    cond_txt = ", ".join(
        f"{c} (n={int((pca['condition'].astype(str) == c).sum())})" for c in conds) or "n/a"
    y = _kv(fig, [
        ("Counts file", os.path.basename(params.get("counts_file", "")) or "n/a"),
        ("Metadata file", os.path.basename(params.get("metadata_file", "")) or "n/a"),
        ("Samples", f"{n_samples}"),
        ("Conditions", cond_txt),
        ("Gene / sample / cond. columns",
         f"{params.get('gene_col', '?')} / {params.get('sample_col', '?')} / "
         f"{params.get('condition_col', '?')}"),
    ], y)

    y = _section(fig, "Analysis parameters", y - 0.005)
    method = params.get("method", "")
    method_txt = {"RNAseq": "DESeq2 (Wald test) — RNA-seq raw counts",
                  "Microarray": "limma (eBayes) — microarray"}.get(method, method or "n/a")
    contrast = params.get("contrast")
    if params.get("pairwise_all"):
        comp = "All pairwise comparisons between groups"
    elif contrast:
        comp = " vs ".join(contrast)
    else:
        comp = "n/a"
    cat = params.get("gsea_category", "") or "n/a"
    if params.get("gsea_subcategory"):
        cat += f" / {params['gsea_subcategory']}"
    ora_txt = "Not run"
    if params.get("run_ora"):
        ora_txt = f"Run on: {params.get('ora_direction', 'all DEGs')}"
    items = [
        ("Statistical method", method_txt),
        ("Comparison", comp),
        ("padj cutoff", f"{params.get('padj_cutoff', '?')}"),
        ("|log2FC| cutoff", f"{params.get('lfc_cutoff', '?')}"),
        ("Minimum counts filter", f"{params.get('min_counts', 0)}"),
        ("MSigDB category (GSEA/ORA)", cat),
        ("ORA", ora_txt),
    ]
    if params.get("started"):
        items.append(("Analysis started", params["started"]))
    if params.get("duration"):
        items.append(("Run time", params["duration"]))
    y = _kv(fig, items, y)

    y = _section(fig, "Results at a glance", y - 0.005)
    rows = []
    for tag, res in list(results.per_contrast.items())[:40]:
        s = stats[tag]
        rows.append([
            _trunc(_contrast_label(res), 34), str(s["n_tested"]), str(s["n_sig"]),
            str(s["n_up"]), str(s["n_down"]),
            "n/a" if s["gsea_sig"] is None else str(s["gsea_sig"]),
            "n/a" if s["ora_sig"] is None else str(s["ora_sig"]),
        ])
    y = _table(fig, rows,
               ["Contrast (A vs B)", "Genes tested", "padj < cutoff", "DEG higher in A", "DEG higher in B",
                "GSEA sig.", "ORA sig."],
               [0.34, 0.12, 0.12, 0.10, 0.10, 0.11, 0.11], y, fs=7)
    _paragraph(fig, "DEG = padj below the cutoff AND |log2FC| at or above the cutoff. "
                    "'GSEA sig.' / 'ORA sig.' = pathways with padj below the cutoff.",
               y, width=120, fs=7, color=C["muted"])
    w.save(fig)


def _page_guide(w: _Writer):
    fig = w.new("How to read this report")
    y = 0.90
    sections = [
        ("PCA",
         "Each point is a sample. Samples with similar expression profiles lie close "
         "together. The percentage on each axis is the share of total variance explained by "
         "that component. Check that samples cluster by condition; an isolated sample or "
         "clustering by something other than the condition may indicate a batch effect or "
         "a low-quality sample."),
        ("Volcano plot",
         "Each point is a gene. The x axis is the log2 fold change, the y axis is "
         "-log10(padj): the higher a point, the stronger the statistical evidence. A "
         "positive log2FC means higher expression in the first group of the contrast "
         "(\"A\" in \"A vs B\"), a negative one higher expression in the second group "
         "(B, the reference). Dashed lines mark the "
         "cutoffs; genes beyond both are the differentially expressed genes (DEGs)."),
        ("Top DEG tables",
         "The most significant genes higher in each group of the contrast (sorted by padj). padj is the "
         "p-value corrected for multiple testing (Benjamini-Hochberg / FDR); 'stat' is the "
         "test statistic of the differential expression test."),
        ("GSEA",
         "Gene Set Enrichment Analysis uses the whole ranked gene list. The Normalized "
         "Enrichment Score (NES) is positive when a pathway is enriched among genes "
         "higher in the first group of the contrast (\"A\" in \"A vs B\") and negative when "
         "enriched among genes higher in the second group (the reference). Asterisks mark padj < 0.05 (*), < 0.01 (**), < 0.001 (***)."),
        ("ORA",
         "Over-Representation Analysis tests whether the selected DEGs are over-represented "
         "in a pathway (hypergeometric test). The background is the set of genes actually "
         "tested. 'Count/SetSize' is the number of input genes in the pathway over the "
         "pathway size, 'FE' is the fold enrichment. Direction is derived afterwards: UP (higher in "
         "group A) / DOWN (higher in group B) if at least 80% of the pathway's input genes "
         "go in that direction, MIXED otherwise."),
        ("Reproducibility",
         "All parameters used for the run are listed on the first page, and the full "
         "processing log is attached at the end of the report."),
    ]
    for title, text in sections:
        y = _section(fig, title, y)
        y = _paragraph(fig, text, y, width=108, fs=8.5)
        y -= 0.012
    w.save(fig)


def _page_pca(w: _Writer, results):
    df = results.pca_data
    fig = w.new("Principal Component Analysis", "Exploratory view of sample similarity")
    if df is None or df.empty or "condition" not in df.columns:
        _note_page_text(fig, "PCA data not available.")
        w.save(fig)
        return
    ax = fig.add_axes([0.11, 0.42, 0.80, 0.46])
    conds = sorted(df["condition"].astype(str).unique())
    for i, cond in enumerate(conds):
        sub = df[df["condition"].astype(str) == cond]
        ax.scatter(sub["PC1"], sub["PC2"], s=55, alpha=0.85, edgecolor="black",
                   linewidth=0.5, color=TAB10[i % len(TAB10)], label=str(cond))
        if len(df) <= 40 and "SampleID" in sub.columns:
            for _, r in sub.iterrows():
                ax.annotate(str(r["SampleID"]), (r["PC1"], r["PC2"]), xytext=(4, 4),
                            textcoords="offset points", fontsize=6)
    ax.set_xlabel(results.pca_x_label, fontsize=8)
    ax.set_ylabel(results.pca_y_label, fontsize=8)
    ax.legend(title="Condition", fontsize=7, title_fontsize=8, frameon=False)
    ax.grid(alpha=0.25)
    _style_axes(ax)

    rows = [[_trunc(c, 40), str(int((df["condition"].astype(str) == c).sum()))] for c in conds]
    y = _section(fig, "Samples per condition", 0.37)
    _table(fig, rows, ["Condition", "Samples"], [0.7, 0.3], y, w=0.5)
    w.save(fig)


def _page_volcano(w: _Writer, res, s: dict, params: dict):
    padj_cut = float(params.get("padj_cutoff", 0.05))
    lfc_cut = float(params.get("lfc_cutoff", 1.0))
    label = _contrast_label(res)
    if res.group_high and res.group_low:
        sub = (f"Positive log2FC = higher expression in \"{res.group_high}\"  |  "
               f"negative log2FC = higher expression in \"{res.group_low}\"")
    else:
        sub = "log2FC direction cannot be determined automatically for this contrast."
    fig = w.new(f"Contrast: {label}", sub)

    y = _kv(fig, [
        ("Genes in results table", str(s["n_rows"])),
        ("Genes tested (non-NA p-value)", str(s["n_tested"])),
        ("padj < %.3g (any fold change)" % padj_cut, str(s["n_sig"])),
        (f'DEGs — higher in "{_dir_names(res)[0]}"', str(s["n_up"])),
        (f'DEGs — higher in "{_dir_names(res)[1]}"', str(s["n_down"])),
    ], 0.895, label_w=0.30)

    df = s["df"].dropna(subset=["padj"]).copy()
    ax = fig.add_axes([0.11, 0.09, 0.82, 0.62])
    if df.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No genes with a valid padj to plot.", ha="center", va="center",
                color=C["muted"])
        w.save(fig)
        return
    df["nl"] = -np.log10(df["padj"].clip(lower=1e-300))
    for status, color in (("Not Significant", C["ns"]), ("Downregulated", C["down"]),
                          ("Upregulated", C["up"])):
        sub_df = df[df["status"] == status]
        if sub_df.empty:
            continue
        ax.scatter(sub_df["log2FoldChange"], sub_df["nl"], s=9, c=color, alpha=0.7,
                   edgecolor="none", label=f"{status} ({len(sub_df)})")
    ax.axvline(lfc_cut, ls="--", lw=0.9, color="#8a8a8a")
    ax.axvline(-lfc_cut, ls="--", lw=0.9, color="#8a8a8a")
    ax.axhline(-np.log10(padj_cut), ls="--", lw=0.9, color="#8a8a8a")
    top = df[df["status"] != "Not Significant"].sort_values(
        ["nl", "log2FoldChange"], ascending=[False, False], key=None).head(10)
    for i, (_, r) in enumerate(top.iterrows()):
        ax.annotate(str(r.get("GeneID", "")), (r["log2FoldChange"], r["nl"]),
                    xytext=(6, 6 + 7 * (i % 2)), textcoords="offset points", fontsize=6.5,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#666666"))
    ax.set_xlabel("log2 Fold Change", fontsize=8)
    ax.set_ylabel("-log10(padj)", fontsize=8)
    ax.set_title(f"Volcano plot — {label}", fontsize=10)
    ax.legend(fontsize=7, frameon=False)
    _style_axes(ax)
    w.save(fig)


def _page_top_degs(w: _Writer, res, s: dict, params: dict, top_n: int = 20):
    label = _contrast_label(res)
    fig = w.new(f"Top DEGs — {label}",
                f"Up to {top_n} most significant genes per direction (sorted by padj)")
    df = s["df"]
    cols = ["GeneID", "log2FoldChange", "pvalue", "padj"] + (["stat"] if "stat" in df else [])
    headers = ["Gene", "log2FC", "p-value", "padj"] + (["stat"] if "stat" in df else [])
    widths = [0.30, 0.17, 0.18, 0.18] + ([0.17] if "stat" in df else [])
    if "stat" not in df:
        widths = [0.34, 0.22, 0.22, 0.22]

    def rows_for(status):
        sub = df[df["status"] == status].sort_values(
            ["padj", "log2FoldChange"], ascending=[True, status == "Downregulated"]).head(top_n)
        out = []
        for _, r in sub.iterrows():
            row = [_trunc(r.get("GeneID", ""), 24), _fmt_num(r.get("log2FoldChange")),
                   _fmt_sci(r.get("pvalue")), _fmt_sci(r.get("padj"))]
            if "stat" in df:
                row.append(_fmt_num(r.get("stat")))
            out.append(row)
        return out

    _hi, _lo = _dir_names(res)
    y = _section(fig, f"Higher in {_trunc(_hi, 40)} (log2FC > 0)", 0.895)
    y = _table(fig, rows_for("Upregulated"), headers, widths, y)
    y = _section(fig, f"Higher in {_trunc(_lo, 40)} (log2FC < 0)", y - 0.01)
    _table(fig, rows_for("Downregulated"), headers, widths, y)
    w.save(fig)


def _page_gsea(w: _Writer, res, params: dict):
    label = _contrast_label(res)
    tbl = res.gsea_tbl
    if tbl is None or tbl.empty:
        fig = w.new(f"GSEA — {label}", "Gene Set Enrichment Analysis (MSigDB)")
        _note_page_text(fig, "GSEA not available for this contrast.")
        w.save(fig)
        return
    padj_cut = float(params.get("padj_cutoff", 0.05))
    top_n = int(params.get("top_n_pathways", 15))
    data, msg = engine.filter_top_pathways(tbl, top_n=top_n, padj_cutoff=padj_cut)
    fig = w.new(f"GSEA — {label}", _trunc(msg, 125))

    ax = fig.add_axes([0.40, 0.50, 0.53, 0.38])
    colors = data["Direction"].map({"Up-regulated": C["up"], "Down-regulated": C["down"]})
    names = [_trunc(n, 42) for n in data["pathway_clean"]]
    bars = ax.barh(names, data["NES"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=8)
    ax.tick_params(axis="y", labelsize=6.5)
    for bar, padj in zip(bars, data["padj"]):
        stars = ""
        if pd.notna(padj):
            stars = "***" if padj < 0.001 else "**" if padj < 0.01 else "*" if padj < 0.05 else ""
        if stars:
            wd = bar.get_width()
            ax.text(wd + (0.05 if wd >= 0 else -0.05), bar.get_y() + bar.get_height() / 2,
                    stars, va="center", ha="left" if wd >= 0 else "right",
                    fontsize=8, fontweight="bold")
    _hi, _lo = _dir_names(res)
    ax.legend(handles=[Patch(color=C["up"], label=f"NES > 0: higher in {_trunc(_hi, 28)}"),
                       Patch(color=C["down"], label=f"NES < 0: higher in {_trunc(_lo, 28)}")],
              fontsize=7, frameon=False, loc="lower right")
    _style_axes(ax)

    rows = []
    for _, r in data.sort_values("padj").iterrows():
        rows.append([_trunc(r.get("pathway_clean", ""), 52), _fmt_num(r.get("NES")),
                     _fmt_sci(r.get("padj")), str(r.get("size", ""))])
    _table(fig, rows, ["Pathway", "NES", "padj", "Size"], [0.60, 0.13, 0.15, 0.12], 0.45)
    w.save(fig)


def _page_ora(w: _Writer, res, params: dict):
    label = _contrast_label(res)
    if res.ora_meta is None:
        fig = w.new(f"ORA — {label}", "Over-Representation Analysis (MSigDB)")
        _note_page_text(fig, "ORA was not run for this analysis.")
        w.save(fig)
        return
    tbl = res.ora_tbl
    if tbl is None or tbl.empty:
        fig = w.new(f"ORA — {label}", "Over-Representation Analysis (MSigDB)")
        note = (res.ora_meta.get("note") or "").strip()
        text = "ORA not available for this contrast."
        if note:
            text += "\n\n" + textwrap.fill(note, 80)
        _note_page_text(fig, text)
        w.save(fig)
        return
    padj_cut = float(params.get("padj_cutoff", 0.05))
    top_n = int(params.get("top_n_ora", 15))
    data, msg = engine.filter_top_ora(tbl, top_n=top_n, padj_cutoff=padj_cut)
    n_in = res.ora_meta.get("n_input_genes")
    fig = w.new(f"ORA — {label}",
                _trunc(msg + (f" Genes in input: {n_in}." if n_in else ""), 125))

    dir_colors = {"UP": C["up"], "DOWN": C["down"], "MIXED": C["ns"]}
    ax = fig.add_axes([0.40, 0.50, 0.53, 0.38])
    colors = data["Direction"].map(dir_colors).fillna(C["ns"])
    names = [_trunc(n, 42) for n in data["pathway_clean"]]
    ax.barh(names, data["neglog10"], color=colors)
    ax.axvline(-np.log10(padj_cut), ls="--", color=C["up"], lw=1)
    ax.set_xlabel(r"$-\log_{10}$(padj)", fontsize=8)
    ax.tick_params(axis="y", labelsize=6.5)
    _hi, _lo = _dir_names(res)
    ax.legend(handles=[Patch(color=C["up"], label=f"UP: higher in {_trunc(_hi, 24)} (>= 80% of genes)"),
                       Patch(color=C["down"], label=f"DOWN: higher in {_trunc(_lo, 24)} (>= 80% of genes)"),
                       Patch(color=C["ns"], label="MIXED")],
              fontsize=7, frameon=False, loc="lower right")
    _style_axes(ax)

    rows = []
    for _, r in data.sort_values("p.adjust", kind="stable").iterrows():
        rows.append([
            _trunc(r.get("pathway_clean", ""), 40), f"{int(r['Count'])}/{int(r['SetSize'])}",
            _fmt_num(r.get("FoldEnrichment"), 2), _fmt_sci(r.get("p.adjust")),
            str(r.get("Direction", "")), str(int(r["N_UP"])), str(int(r["N_DOWN"])),
        ])
    _table(fig, rows, ["Pathway", "Count/Size", "FE", "padj", "Direction", "N up", "N down"],
           [0.38, 0.13, 0.08, 0.14, 0.11, 0.08, 0.08], 0.45)
    w.save(fig)


def _pages_log(w: _Writer, log_text: str):
    lines = []
    for raw in (log_text or "").replace("\r", "\n").split("\n"):
        raw = raw.rstrip()
        if not raw:
            continue
        lines.extend(textwrap.wrap(raw, width=125, subsequent_indent="    ") or [""])
    if not lines:
        return
    max_lines = LOG_LINES_PER_PAGE * LOG_MAX_PAGES
    truncated = len(lines) > max_lines
    lines = lines[:max_lines]
    n_pages = (len(lines) + LOG_LINES_PER_PAGE - 1) // LOG_LINES_PER_PAGE
    for p in range(n_pages):
        chunk = lines[p * LOG_LINES_PER_PAGE:(p + 1) * LOG_LINES_PER_PAGE]
        fig = w.new("Appendix — Processing log", f"Part {p + 1} of {n_pages}")
        y = 0.905
        for ln in chunk:
            fig.text(0.05, y, ln, fontsize=6, family="DejaVu Sans Mono",
                     color=C["text"], va="top")
            y -= 0.0098
        if truncated and p == n_pages - 1:
            fig.text(0.05, y - 0.005, "… log truncated …", fontsize=6.5, color=C["muted"])
        w.save(fig)


def build_report(path: str, results, params: dict | None = None, log_text: str = "") -> str:
    if results is None or not results.per_contrast:
        raise ValueError("There are no results to put in the report.")
    params = dict(params or {})
    padj_cut = float(params.get("padj_cutoff", 0.05))
    lfc_cut = float(params.get("lfc_cutoff", 1.0))
    generated = datetime.datetime.now().strftime("%d %B %Y, %H:%M")

    stats = {tag: _contrast_stats(res, padj_cut, lfc_cut)
             for tag, res in results.per_contrast.items()}

    with PdfPages(path) as pdf:
        info = pdf.infodict()
        info["Title"] = "DEA Explorer — Differential Expression Analysis Report"
        info["Author"] = "DEA Explorer"
        info["Subject"] = "Differential expression, GSEA and ORA"
        info["CreationDate"] = datetime.datetime.now()

        w = _Writer(pdf, generated)
        _page_cover(w, results, params, stats, generated)
        _page_guide(w)
        _page_pca(w, results)
        for tag, res in results.per_contrast.items():
            s = stats[tag]
            _page_volcano(w, res, s, params)
            _page_top_degs(w, res, s, params)
            _page_gsea(w, res, params)
            _page_ora(w, res, params)
        _pages_log(w, log_text)
    return path
