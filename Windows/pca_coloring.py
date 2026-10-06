"""Colouring of the PCA by any metadata column.

Numeric column (e.g. age)      -> continuous colour scale + colour bar.
Categorical column (e.g. stage) -> one clearly distinct colour per group.
Missing values ("not available", "NA", empty...) are always drawn in grey.

Pure pandas/matplotlib: used both by gui.py and report_pdf.py.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

DEFAULT_COLUMN = "(condition)"          # entry meaning "colour by condition, as before"
MISSING_COLOR = "#bdbdbd"
MISSING_LABEL = "Missing / n.a."

MISSING_TOKENS = {
    "", "na", "n/a", "nan", "none", "null", "nd", "n.a.", "not available",
    "not reported", "not applicable", "unknown", "not specified",
    "not evaluated", "--", "-", "'--",
}

OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
             "#0072B2", "#D55E00", "#CC79A7", "#000000"]

# label shown in the GUI -> matplotlib colormap name
CMAPS = {
    "Red → Blue": "RdYlBu",
    "Blue → Red": "RdYlBu_r",
    "Viridis": "viridis",
    "Plasma": "plasma",
    "Coolwarm": "coolwarm",
}

KIND_CHOICES = ["Auto", "Numeric", "Categorical"]


# --------------------------------------------------------------------------- #
# Metadata reading
# --------------------------------------------------------------------------- #
def read_metadata(path: str) -> pd.DataFrame:
    """Reads the metadata as strings (no automatic type guessing)."""
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        df = pd.read_excel(p, dtype=str)
    elif suffix == ".json":
        df = pd.read_json(p, dtype=False).astype(str)
    else:
        sep = "\t" if suffix in (".tsv", ".txt") else ","
        df = pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def r_clean_id(value) -> str:
    """Python port of clean_identifiers() of the R core (spaces -> '_', make.names(), '.' -> '_').

    The R pipeline renames the samples this way (e.g. 'TCGA-AB-1234' -> 'TCGA_AB_1234'), so the
    sample IDs of the PCA / normalized matrix have to be compared with the metadata IDs after the
    same cleaning, otherwise nothing matches."""
    s = str(value).strip().replace(" ", "_")
    s = re.sub(r"[^\w.]", ".", s)
    if s == "":
        s = "X"
    if re.match(r"^([0-9_]|\.[0-9])", s):
        s = "X" + s
    return s.replace(".", "_")


def resolve_column(df: pd.DataFrame, col) -> str | None:
    """col may be a 1-based index (int) or a column name."""
    cols = list(df.columns)
    if isinstance(col, int):
        return cols[col - 1] if 1 <= col <= len(cols) else None
    col = str(col).strip()
    return col if col in cols else None


# --------------------------------------------------------------------------- #
# Type detection
# --------------------------------------------------------------------------- #
def _clean(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip()
    s = s.where(~s.str.lower().isin(MISSING_TOKENS), other=np.nan)
    return s


def classify(series: pd.Series, override: str = "Auto") -> str:
    """'numeric' or 'categorical'."""
    if override.lower().startswith("num"):
        return "numeric"
    if override.lower().startswith("cat"):
        return "categorical"
    s = _clean(series).dropna()
    if s.empty:
        return "categorical"
    num = pd.to_numeric(s, errors="coerce")
    frac_numeric = num.notna().mean()
    n_unique = num.dropna().nunique()
    # numeric if (almost) everything is a number and there are enough distinct
    # values to make a scale meaningful (0/1/2-coded columns stay categorical)
    if frac_numeric >= 0.9 and n_unique > 6:
        return "numeric"
    return "categorical"


def _natural_key(text: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(text))]


def categorical_palette(n: int, colorblind: bool) -> list:
    base = OKABE_ITO if colorblind else list(plt.cm.tab10.colors)
    if n <= len(base):
        return base[:n]
    if n <= 20:
        return list(plt.cm.tab20.colors)[:n]
    return [plt.cm.gist_rainbow(i / n) for i in range(n)]


# --------------------------------------------------------------------------- #
# Main entry point
# --------------------------------------------------------------------------- #
def build_spec(pca_df: pd.DataFrame, meta_df: pd.DataFrame, sample_col,
               column: str, kind: str = "Auto", colorblind: bool = False,
               cmap_label: str = "Red → Blue") -> dict:
    """Returns a dict describing how to colour each PCA point.

    keys: mode ('categorical'|'numeric'), column, values (Series aligned to pca_df),
          groups (categorical: ordered list of (label, colour)),
          cmap, norm (numeric), n_missing, n_unmatched.
    Raises ValueError with a readable message if the column cannot be used.
    """
    sample_name = resolve_column(meta_df, sample_col)
    if sample_name is None:
        raise ValueError("Sample column not found in the metadata.")
    if column not in meta_df.columns:
        raise ValueError(f"Column '{column}' not found in the metadata.")
    if "SampleID" not in pca_df.columns:
        raise ValueError("The PCA data has no 'SampleID' column: cannot match metadata.")

    # Match on the raw ID first, then on the R-cleaned ID (the PCA/normalized-matrix IDs come from R).
    sid = pca_df["SampleID"].astype(str).str.strip()
    ids = meta_df[sample_name].astype(str).str.strip()
    lookup = meta_df.assign(_sid=ids).drop_duplicates("_sid").set_index("_sid")[column]
    lookup_clean = (meta_df.assign(_sid=ids.map(r_clean_id))
                           .drop_duplicates("_sid").set_index("_sid")[column])
    raw = sid.map(lookup)
    raw = raw.where(sid.isin(lookup.index), sid.map(lookup_clean))
    n_unmatched = int((~(sid.isin(lookup.index) | sid.isin(lookup_clean.index))).sum())
    vals = _clean(raw.astype(object))
    mode = classify(vals, kind)

    spec = {"column": column, "n_unmatched": n_unmatched}
    if mode == "numeric":
        num = pd.to_numeric(vals, errors="coerce")
        if num.notna().sum() == 0:
            raise ValueError(f"Column '{column}' has no numeric values.")
        vmin, vmax = float(num.min()), float(num.max())
        if vmin == vmax:
            vmax = vmin + 1.0
        spec.update(mode="numeric", values=num,
                    cmap=plt.get_cmap(CMAPS.get(cmap_label, "RdYlBu")),
                    norm=Normalize(vmin=vmin, vmax=vmax),
                    n_missing=int(num.isna().sum()))
    else:
        labels = vals.dropna().unique().tolist()
        labels.sort(key=_natural_key)
        colors = categorical_palette(len(labels), colorblind)
        spec.update(mode="categorical", values=vals,
                    groups=list(zip(labels, colors)),
                    n_missing=int(vals.isna().sum()))
    return spec


def draw_pca_points(ax, pca_df: pd.DataFrame, spec: dict, size=80, alpha=0.85,
                    ellipse_fn=None, xcol: str = "PC1", ycol: str = "PC2"):
    """Draws the points on `ax`. Returns (list_of_(artist, sub_df), legend_handles).

    legend_handles is a list of (label, colour) for categorical mode and the
    optional 'missing' entry; for numeric mode the colour bar must be added by
    the caller using spec['cmap'] / spec['norm'] (see add_colorbar).
    """
    artists = []
    handles = []
    df = pca_df.copy()
    df["_val"] = spec["values"].values

    if spec["mode"] == "numeric":
        has = df[df["_val"].notna()]
        miss = df[df["_val"].isna()]
        if not miss.empty:
            sc = ax.scatter(miss[xcol], miss[ycol], color=MISSING_COLOR, s=size,
                            alpha=alpha, edgecolor="black", linewidth=0.5)
            artists.append((sc, miss))
            handles.append((MISSING_LABEL, MISSING_COLOR))
        sc = ax.scatter(has[xcol], has[ycol], c=has["_val"].astype(float),
                        cmap=spec["cmap"], norm=spec["norm"], s=size, alpha=alpha,
                        edgecolor="black", linewidth=0.5)
        artists.append((sc, has))
        spec["_mappable"] = sc
    else:
        for label, color in spec["groups"]:
            sub = df[df["_val"] == label]
            if sub.empty:
                continue
            sc = ax.scatter(sub[xcol], sub[ycol], color=color, s=size, alpha=alpha,
                            edgecolor="black", linewidth=0.5, label=label)
            artists.append((sc, sub))
            if ellipse_fn is not None:
                ellipse_fn(sub[xcol].values, sub[ycol].values, ax, color)
        miss = df[df["_val"].isna()]
        if not miss.empty:
            sc = ax.scatter(miss[xcol], miss[ycol], color=MISSING_COLOR, s=size,
                            alpha=alpha, edgecolor="black", linewidth=0.5,
                            label=MISSING_LABEL)
            artists.append((sc, miss))
    return artists, handles


def add_colorbar(fig, ax, spec: dict, fontsize=None):
    cb = fig.colorbar(spec["_mappable"], ax=ax, pad=0.02, fraction=0.05)
    cb.set_label(spec["column"], fontsize=fontsize)
    return cb


def condition_spec(pca_df: pd.DataFrame, colorblind: bool = False) -> dict:
    """Default colouring: by condition (same behaviour as before)."""
    vals = pca_df["condition"].astype(str)
    labels = sorted(vals.unique())
    colors = categorical_palette(len(labels), colorblind)
    return {"mode": "categorical", "column": "Condition", "values": vals,
            "groups": list(zip(labels, colors)), "n_missing": 0, "n_unmatched": 0}


def uniform_spec(df: pd.DataFrame, color: str = "#2f74ad", label: str = "Samples") -> dict:
    """Spec that paints every point with the same colour (used when no metadata colouring is chosen)."""
    return {"mode": "categorical", "column": label, "values": pd.Series([label] * len(df), index=df.index),
            "groups": [(label, color)], "n_missing": 0, "n_unmatched": 0, "uniform": True}


def spec_from_settings(pca_df: pd.DataFrame, meta_path: str, sample_col, cfg: dict) -> dict:
    """cfg: {'column','kind','colorblind','cmap'} -> spec (falls back to condition)."""
    col = (cfg or {}).get("column", DEFAULT_COLUMN)
    cb = bool((cfg or {}).get("colorblind", False))
    if col == DEFAULT_COLUMN or not meta_path:
        return condition_spec(pca_df, cb)
    meta = read_metadata(meta_path)
    return build_spec(pca_df, meta, sample_col, col, cfg.get("kind", "Auto"), cb,
                      cfg.get("cmap", "Red → Blue"))
