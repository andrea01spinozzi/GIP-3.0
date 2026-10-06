"""Gene expression (normalized counts) of the single samples vs a metadata column.

Python counterpart of plot_gene_vs_metadata() of DE_scheletro_FINALE.Rmd (section 8c):
  - numeric metadata (age, ...)     -> dot plot, X = expression, Y = metadata, Spearman correlation
  - categorical metadata (sex, ...) -> horizontal box plot, X = expression, Y = group, Mann-Whitney / Kruskal-Wallis

The POINTS can additionally be coloured by any other metadata column (same rules as the PCA:
numeric -> colour scale + colour bar, categorical -> one colour per group, missing -> grey).

With 'show_trend' on numeric metadata a linear regression is added: fitted line, 95% confidence band
and, in the subtitle, equation, R^2 and p-value of the fit (needs scipy for band and statistics).

Pure pandas/matplotlib (+ scipy, optional, only for the p-values). Used by gui.py and report_pdf.py.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import pca_coloring as pcol

try:                                    # optional: without scipy the plots work, only the p-values are missing
    from scipy import stats as _st
    HAS_SCIPY = True
except Exception:                       # pragma: no cover
    _st = None
    HAS_SCIPY = False

NO_COLOR = "(none)"                     # entry of the 'Color points by' list meaning "single colour"
POINT_COLOR = "#2f74ad"
BOX_COLORS = ["#9ecae1", "#fdae6b", "#a1d99b", "#bcbddc", "#fcbba1", "#c7e9c0", "#dadaeb", "#d9d9d9"]


# --------------------------------------------------------------------------- #
# Gene lookup
# --------------------------------------------------------------------------- #
def find_gene(norm_mat: pd.DataFrame, gene: str) -> str | None:
    """Exact match first, then case-insensitive (mouse/rat/zebrafish symbols are not upper case)."""
    g = str(gene).strip()
    if not g:
        return None
    if g in norm_mat.index:
        return g
    low = norm_mat.index.str.lower()
    hit = np.flatnonzero(low == g.lower())
    return str(norm_mat.index[hit[0]]) if len(hit) else None


def suggest_genes(norm_mat: pd.DataFrame, text: str, n: int = 8) -> list[str]:
    t = str(text).strip().lower()
    if not t:
        return []
    idx = norm_mat.index
    starts = [g for g in idx[idx.str.lower().str.startswith(t)][:n]]
    if len(starts) < n:
        contains = [g for g in idx[idx.str.lower().str.contains(t, regex=False)] if g not in starts]
        starts += contains[: n - len(starts)]
    return [str(g) for g in starts]


# --------------------------------------------------------------------------- #
# Data preparation
# --------------------------------------------------------------------------- #
def build_frame(norm_mat: pd.DataFrame, meta_df: pd.DataFrame, sample_col, gene: str,
                variable: str, kind: str = "Auto") -> dict:
    """Aligns expression of `gene` with the metadata column `variable`.

    Returns {'df': DataFrame(SampleID, expr, val), 'class': 'numeric'|'categorical',
             'levels': [...]|None, 'gene': resolved gene name, 'variable': variable}.
    Raises ValueError with a readable message when the plot cannot be made.
    """
    gene_name = find_gene(norm_mat, gene)
    if gene_name is None:
        sug = suggest_genes(norm_mat, gene)
        hint = f" Similar names: {', '.join(sug)}." if sug else ""
        raise ValueError(f"Gene '{gene}' not found in the normalized matrix (exact symbol? "
                         f"removed by the low-count filter?).{hint}")
    if variable not in meta_df.columns:
        raise ValueError(f"Column '{variable}' not found in the metadata.")

    expr = norm_mat.loc[gene_name]
    df = pd.DataFrame({"SampleID": expr.index.astype(str), "expr": expr.values.astype(float)})

    spec = pcol.build_spec(df, meta_df, sample_col, variable, kind=kind)
    df["val"] = spec["values"].values
    df = df[df["expr"].notna() & df["val"].notna()].reset_index(drop=True)
    if len(df) < 4:
        raise ValueError("Fewer than 4 samples with complete data: plot not meaningful.")

    levels = None
    if spec["mode"] == "categorical":
        df["val"] = df["val"].astype(str)
        levels = [lab for lab, _ in spec["groups"] if (df["val"] == lab).any()]
        if len(levels) < 2:
            raise ValueError("Fewer than 2 groups with data: cannot compare groups.")
    else:
        df["val"] = df["val"].astype(float)
    return {"df": df, "class": spec["mode"], "levels": levels, "gene": gene_name,
            "variable": variable, "n_unmatched": spec.get("n_unmatched", 0)}


def compute_stats(frame: dict, with_regression: bool = False) -> dict:
    """{'test','rho','pvalue','n'} (+ 'slope','intercept','r2','p_lin' when with_regression).
    pvalue is None when scipy is not installed."""
    df = frame["df"]
    out = {"test": None, "rho": None, "pvalue": None, "n": int(len(df)),
           "slope": None, "intercept": None, "r2": None, "p_lin": None}
    if frame["class"] == "numeric":
        out["test"] = "Spearman"
        if HAS_SCIPY and df["expr"].nunique() > 1 and df["val"].nunique() > 1:
            rho, p = _st.spearmanr(df["expr"], df["val"])
            out["rho"], out["pvalue"] = float(rho), float(p)
            if with_regression:
                lr = _st.linregress(df["expr"], df["val"])
                out.update(slope=float(lr.slope), intercept=float(lr.intercept),
                           r2=float(lr.rvalue ** 2), p_lin=float(lr.pvalue))
        return out
    groups = [df.loc[df["val"] == lv, "expr"].values for lv in frame["levels"]]
    if len(groups) == 2:
        out["test"] = "Mann-Whitney U"
        if HAS_SCIPY:
            try:
                p = _st.mannwhitneyu(groups[0], groups[1], alternative="two-sided",
                                     use_continuity=True, method="asymptotic").pvalue
            except TypeError:               # old scipy: no 'method' argument
                p = _st.mannwhitneyu(groups[0], groups[1], alternative="two-sided",
                                     use_continuity=True).pvalue
            out["pvalue"] = float(p)
    else:
        out["test"] = "Kruskal-Wallis"
        if HAS_SCIPY:
            try:
                out["pvalue"] = float(_st.kruskal(*groups).pvalue)
            except ValueError:              # all values identical
                out["pvalue"] = None
    return out


def subtitle_text(stats: dict) -> str:
    if stats["pvalue"] is None:
        extra = "" if HAS_SCIPY else " | p-value unavailable (install scipy)"
        return f"{stats['test']} | n = {stats['n']}{extra}"
    if stats["rho"] is not None:
        txt = f"{stats['test']} rho = {stats['rho']:.2f} | p = {stats['pvalue']:.3g} | n = {stats['n']}"
        if stats.get("slope") is not None:
            sign = "+" if stats["intercept"] >= 0 else "-"
            txt += (f"\nLinear fit: y = {stats['slope']:.3g}·x {sign} {abs(stats['intercept']):.3g}"
                    f" | R² = {stats['r2']:.2f} | p = {stats['p_lin']:.3g}")
        return txt
    return f"{stats['test']} | p = {stats['pvalue']:.3g} | n = {stats['n']}"


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #
def color_spec_for(frame: dict, meta_df: pd.DataFrame, sample_col, color_by: str | None,
                   kind: str = "Auto", colorblind: bool = False, cmap_label: str = "Red → Blue") -> dict:
    """Colour specification for the points of `frame` (aligned to frame['df'] rows)."""
    df = frame["df"]
    if not color_by or color_by == NO_COLOR:
        return pcol.uniform_spec(df, POINT_COLOR, label="Samples")
    return pcol.build_spec(df[["SampleID"]], meta_df, sample_col, color_by, kind, colorblind, cmap_label)


def draw_gene_plot(ax, frame: dict, color_spec: dict, x_label: str, size: float = 60,
                   alpha: float = 0.85, title: str | None = None, show_sample_labels: bool = False,
                   show_trend: bool = False, seed: int = 0):
    """Draws the plot on `ax`. Returns (artists, stats, extra_handles): artists = [(PathCollection, sub_df), ...]
    where sub_df has the columns SampleID, expr, val, _val (value of the colouring column);
    extra_handles = [(label, colour)] legend entry for missing values (numeric colouring only)."""
    df = frame["df"].copy()
    is_num = frame["class"] == "numeric"
    stats = compute_stats(frame, with_regression=show_trend and is_num)
    var = frame["variable"]

    if is_num:
        df["_y"] = df["val"]
        if show_trend and df["expr"].nunique() > 1:
            x = df["expr"].to_numpy(dtype=float)
            y = df["_y"].to_numpy(dtype=float)
            xs = np.linspace(x.min(), x.max(), 100)
            k, b = np.polyfit(x, y, 1)
            ax.plot(xs, k * xs + b, color="black", lw=1.2, ls="--", zorder=2)
            if HAS_SCIPY and len(x) > 2:        # 95% confidence band of the fitted line
                n = len(x)
                s = np.sqrt(np.sum((y - (k * x + b)) ** 2) / (n - 2))
                sxx = np.sum((x - x.mean()) ** 2)
                se = s * np.sqrt(1 / n + (xs - x.mean()) ** 2 / sxx)
                t = _st.t.ppf(0.975, n - 2)
                ax.fill_between(xs, k * xs + b - t * se, k * xs + b + t * se,
                                color="black", alpha=0.10, lw=0, zorder=1)
        ylabel = var
    else:
        levels = frame["levels"]
        pos = {lv: i for i, lv in enumerate(levels)}
        rng = np.random.default_rng(seed)
        df["_y"] = df["val"].map(pos).astype(float) + rng.uniform(-0.12, 0.12, len(df))
        data = [df.loc[df["val"] == lv, "expr"].values for lv in levels]
        n_by = [len(d) for d in data]
        # With a single colour for the points the boxes are coloured per group (like the R plot);
        # when the points carry the colour of another column the boxes stay neutral.
        neutral = not color_spec.get("uniform")
        bp = ax.boxplot(data, positions=list(range(len(levels))), vert=False, widths=0.5,
                        showfliers=False, patch_artist=True, manage_ticks=False,
                        medianprops={"color": "black", "lw": 1.4})
        for i, patch in enumerate(bp["boxes"]):
            patch.set_facecolor("#e3e3e3" if neutral else BOX_COLORS[i % len(BOX_COLORS)])
            patch.set_alpha(0.65)
        ax.set_yticks(range(len(levels)))
        ax.set_yticklabels([f"{lv} (n={n})" for lv, n in zip(levels, n_by)])
        ax.invert_yaxis()                       # first group at the top
        ylabel = var

    artists, handles = pcol.draw_pca_points(ax, df, color_spec, size=size, alpha=alpha, xcol="expr", ycol="_y")
    for sc, _sub in artists:
        sc.set_zorder(3)                        # points always above boxes / trend line
    if show_sample_labels:
        for _, r in df.iterrows():
            ax.annotate(str(r["SampleID"]), (r["expr"], r["_y"]), xytext=(4, 4),
                        textcoords="offset points", fontsize=6)

    ax.set_xlabel(x_label)
    ax.set_ylabel(ylabel)
    two_lines = stats.get("slope") is not None
    ax.set_title(title or frame["gene"], fontweight="bold", loc="left", pad=30 if two_lines else 16)
    ax.text(0.0, 1.01, subtitle_text(stats), transform=ax.transAxes, fontsize=8, va="bottom", color="#444444")
    return artists, stats, handles


def add_legend_or_colorbar(fig, ax, color_spec: dict, extra_handles=None, legend_fn=None):
    """Numeric colouring -> colour bar; categorical -> legend (via legend_fn(title=...) if given).
    Returns the colour bar (or None)."""
    from matplotlib.lines import Line2D
    if color_spec["mode"] == "numeric":
        cb = pcol.add_colorbar(fig, ax, color_spec)
        if extra_handles:
            handles = [Line2D([], [], marker="o", ls="", markerfacecolor=c, markeredgecolor="black", label=l)
                       for l, c in extra_handles]
            (legend_fn or (lambda **k: ax.legend(frameon=False, **k)))(handles=handles)
        return cb
    if not color_spec.get("uniform"):
        (legend_fn or (lambda **k: ax.legend(frameon=False, **k)))(title=color_spec["column"])
    return None
