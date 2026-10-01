"""
core_calcolo.py
================
Motore di calcolo lato Python per l'applicazione di analisi di espressione
differenziale (DEA + GSEA).

Questo modulo NON reimplementa la statistica: DESeq2, limma, fgsea e
clusterProfiler (ORA) restano esclusivamente nel file R del core di calcolo
(DE_scheletro_FINALE.Rmd). Qui ci occupiamo solo di:

  - costruire il file di configurazione JSON per la pipeline R;
  - lanciare Rscript come sottoprocesso esterno (subprocess);
  - leggere i risultati grezzi prodotti (CSV/JSON) in DataFrame pandas;
  - fornire funzioni di supporto "leggere" (filtri su q-value/log2FC, pulizia
    nomi pathway, ecc.) usate SOLO per rendere i grafici interattivi lato
    Python, senza toccare i numeri calcolati da R.

Il modulo e' pensato per essere usato dalla GUI (gui.py) ma non dipende da
tkinter: puo' essere testato/usato anche da riga di comando o da altri script.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

# -----------------------------------------------------------------------------
# Percorsi
# -----------------------------------------------------------------------------
def _base_dir() -> Path:
    """Cartella dei file di risorsa: sorgenti in sviluppo, bundle .app (PyInstaller) quando 'frozen'."""
    if getattr(sys, "frozen", False):
        return Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
    return Path(__file__).resolve().parent


THIS_DIR = _base_dir()


def resource_path(name: str) -> str:
    """Percorso assoluto di un file di risorsa (logo.png, ecc.), valido anche nell'app impacchettata."""
    return str(THIS_DIR / name)

# Pacchetti R richiesti dal core di calcolo (usati per il controllo iniziale).
REQUIRED_R_PACKAGES = ["tidyverse", "readxl", "jsonlite", "DESeq2", "msigdbr", "fgsea",
                       "clusterProfiler", "org.Hs.eg.db", "limma", "BiocParallel"]
CORE_R_PATH = THIS_DIR / "DE_scheletro_FINALE.Rmd"
RUNNER_R_PATH = THIS_DIR / "run_pipeline.R"


class RPipelineError(RuntimeError):
    """Sollevata quando run_pipeline.R fallisce (errore statistico, file non
    trovato, formato non supportato, pacchetto R mancante, ecc.)."""
    pass


# -----------------------------------------------------------------------------
# Individuazione di Rscript
# -----------------------------------------------------------------------------
def _subprocess_env() -> dict:
    """Ambiente per Rscript: su macOS l'app da Finder ha PATH/locale minimi e, se impacchettata,
    variabili DYLD_* che potrebbero far caricare a R librerie sbagliate."""
    env = os.environ.copy()
    if sys.platform == "darwin":
        env["PATH"] = ":".join(["/Library/Frameworks/R.framework/Resources/bin", "/opt/homebrew/bin",
                                "/usr/local/bin", env.get("PATH", "")])
        env.setdefault("LANG", "en_US.UTF-8")
        env.setdefault("LC_ALL", "en_US.UTF-8")
        if getattr(sys, "frozen", False):
            for k in ("DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH", "DYLD_INSERT_LIBRARIES"):
                env.pop(k, None)
    return env


def find_rscript() -> Optional[str]:
    """Cerca l'eseguibile Rscript nel PATH di sistema e in alcune posizioni
    tipiche di Windows/macOS, nel caso non sia nel PATH."""
    override = os.environ.get("DEA_RSCRIPT")
    if override and Path(override).exists():
        return override

    exe = shutil.which("Rscript") or shutil.which("Rscript.exe")
    if exe:
        return exe

    candidates = []
    
    # Ricerca robusta tramite Registro di Sistema (Windows)
    if os.name == "nt":
        try:
            import winreg
            for subkey in [r"SOFTWARE\R-core\R", r"SOFTWARE\R-core\R64"]:
                try:
                    with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, subkey) as key:
                        r_path, _ = winreg.QueryValueEx(key, "InstallPath")
                        cand = Path(r_path) / "bin" / "Rscript.exe"
                        if cand.exists():
                            candidates.append(str(cand))
                except EnvironmentError:
                    continue
        except ImportError:
            pass

    # Windows: Fallback su Program Files
    program_files = [os.environ.get("ProgramFiles", r"C:\Program Files"),
                     os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")]
    for pf in program_files:
        r_root = Path(pf) / "R"
        if r_root.is_dir():
            for sub in sorted(r_root.glob("R-*"), reverse=True):
                cand = sub / "bin" / "Rscript.exe"
                if cand.exists():
                    candidates.append(str(cand))

    # macOS: un'app avviata dal Finder ha un PATH minimo, quindi cerchiamo nei
    # percorsi standard (installer CRAN = R.framework, Homebrew Intel/Apple Silicon).
    if sys.platform == "darwin":
        mac_paths = [Path("/Library/Frameworks/R.framework/Resources/bin/Rscript"),
                     Path("/opt/homebrew/bin/Rscript"),
                     Path("/usr/local/bin/Rscript")]
        mac_paths += sorted(Path("/Library/Frameworks/R.framework/Versions").glob("*/Resources/bin/Rscript"),
                            reverse=True)
        for cand in mac_paths:
            if cand.exists():
                candidates.append(str(cand))

    return candidates[0] if candidates else None


def missing_r_packages(rscript: str) -> list[str]:
    """Elenca i pacchetti R richiesti non installati (veloce, non carica i pacchetti)."""
    expr = ('pk <- c(%s); cat(setdiff(pk, rownames(installed.packages())), sep="\\n")'
            % ",".join(f'"{p}"' for p in REQUIRED_R_PACKAGES))
    try:
        out = subprocess.run([rscript, "--vanilla", "-e", expr], capture_output=True, text=True,
                             timeout=60, env=_subprocess_env()).stdout
        return [line.strip() for line in out.splitlines() if line.strip()]
    except Exception:
        return []   # in caso di dubbio non bloccare l'avvio


def check_environment() -> list[str]:
    """Restituisce una lista di problemi di ambiente (vuota se tutto ok)."""
    problems = []
    rscript = find_rscript()
    if rscript is None:
        hint = ("Install R from https://cran.r-project.org/bin/macosx/ (choose the arm64 build "
                "on Apple Silicon Macs)." if sys.platform == "darwin"
                else "Install R and make sure that 'Rscript' is reachable from the terminal.")
        problems.append("Rscript not found. " + hint)
    else:
        missing = missing_r_packages(rscript)
        if missing:
            problems.append("Missing R packages: " + ", ".join(missing) +
                            ". Run once in Terminal: Rscript install_r_packages.R")
    if not CORE_R_PATH.exists():
        problems.append(f"Calculation core file not found at: {CORE_R_PATH}")
    if not RUNNER_R_PATH.exists():
        problems.append(f"run_pipeline.R not found at: {RUNNER_R_PATH}")
    return problems


# -----------------------------------------------------------------------------
# Validazione Preventiva & Diagnostica Errori R
# -----------------------------------------------------------------------------
def parse_r_error(raw_error: str) -> str:
    """Translates the most common R error messages into clear instructions."""
    raw_err_lower = raw_error.lower()

    # 1. Missing package
    if "there is no package called" in raw_err_lower:
        match = re.search(r"no package called ['\"](.*?)['\"]", raw_error)
        pkg = match.group(1) if match else "required"
        return (
            f"❌ R ENVIRONMENT ERROR:\n"
            f"The R package '{pkg}' does not appear to be installed.\n"
            f"Make sure all the pipeline dependencies are installed in R."
        )

    # 2. Design matrix not full rank (collinearity or no replicates)
    if "design matrix not full rank" in raw_err_lower or "not full rank" in raw_err_lower:
        return (
            "❌ EXPERIMENTAL DESIGN ERROR:\n"
            "The design matrix is not full rank (it cannot be computed). This happens if:\n"
            "  - One of your experimental variables is perfectly collinear with another.\n"
            "  - There are conditions/groups with fewer than 2 samples (no biological replicates).\n"
            "  - One of your groups contains exactly the same values as another."
        )

    # 3. Locked or inaccessible files
    if "cannot open connection" in raw_err_lower or "cannot open file" in raw_err_lower:
        return (
            "❌ FILE ACCESS ERROR:\n"
            "Unable to read or write one of the input/output files.\n"
            "Check that the files are not open in Excel or other programs and "
            "that you have write permission on the output folder."
        )

    # 4. Too many zeros for the geometric computation (DESeq2)
    if "every gene has at least one zero" in raw_err_lower:
        return (
            "❌ DESeq2 ERROR (VALID COUNTS):\n"
            "Every gene in your dataset has at least one zero count in some sample.\n"
            "In this situation DESeq2 cannot compute the standard geometric normalization.\n"
            "Check that you have not provided already normalized/log-transformed data, or apply "
            "a stricter gene pre-filter."
        )

    return raw_error  # Return the untranslated error if it matches none of the cases above


def validate_inputs(counts_path: str, metadata_path: str,
                    sample_col: int | str = 1, condition_col: int | str = 2,
                    gene_col: int | str = 1, min_reads: int = 10, min_samples: int = 3) -> None:
    """Runs integrity checks before launching R, to prevent silent or cryptic crashes."""
    c_path = Path(counts_path)
    m_path = Path(metadata_path)

    if not c_path.exists():
        raise ValueError(f"The counts file does not exist: {counts_path}")
    if not m_path.exists():
        raise ValueError(f"The metadata file does not exist: {metadata_path}")

    # Validation of the numeric pre-filtering parameters
    if min_reads < 0:
        raise ValueError("The minimum number of reads (min_reads) cannot be negative.")
    if min_samples < 1:
        raise ValueError("The minimum number of samples (min_samples) must be at least 1.")

    def _read_flexible(path: Path, nrows: int | None = None) -> pd.DataFrame:
        sep = "\t" if path.suffix in [".tsv", ".txt"] else ","
        return pd.read_csv(path, sep=sep, nrows=nrows)

    # Read metadata
    try:
        meta_df = _read_flexible(m_path)
    except Exception as e:
        raise ValueError(f"Unable to read the metadata file: {e}")

    if meta_df.empty:
        raise ValueError("The metadata file is empty.")

    # Resolve metadata columns (handles both numbers and strings)
    cols = list(meta_df.columns)

    if isinstance(sample_col, int):
        if sample_col < 1 or sample_col > len(cols):
            raise ValueError(f"Sample column index ({sample_col}) out of range. The file has {len(cols)} columns.")
        sample_name = cols[sample_col - 1]
    else:
        if sample_col not in cols:
            raise ValueError(f"The sample column '{sample_col}' does not exist in the metadata. Available: {cols}")
        sample_name = sample_col

    if isinstance(condition_col, int):
        if condition_col < 1 or condition_col > len(cols):
            raise ValueError(f"Condition column index ({condition_col}) out of range.")
        cond_name = cols[condition_col - 1]
    else:
        if condition_col not in cols:
            raise ValueError(f"The condition column '{condition_col}' does not exist in the metadata.")
        cond_name = condition_col

    meta_df[cond_name] = (
        meta_df[cond_name]
        .astype(str)
        .str.replace(r"[^a-zA-Z0-9_]", "_", regex=True)  # Replace everything that is not a letter/digit/_
        .str.replace(r"__+", "_", regex=True)            # Avoid adjacent underscores (e.g. "__")
        .str.strip("_"))

    # Check that the minimum-samples filter does not exceed the total number of samples
    samples_in_meta = meta_df[sample_name].dropna().astype(str).tolist()
    total_samples = len(samples_in_meta)
    if min_samples > total_samples:
        raise ValueError(
            f"Invalid filter: you require at least {min_samples} samples, "
            f"but your experimental design has only {total_samples} samples in total."
        )

    # Check for special characters in the factors (not allowed in R condition names)
    invalid_char_pattern = re.compile(r"[^a-zA-Z0-9_]")
    detected_conditions = meta_df[cond_name].dropna().unique()
    for cond in detected_conditions:
        cond_str = str(cond)
        if invalid_char_pattern.search(cond_str):
            raise ValueError(
                f"The condition '{cond_str}' contains invalid characters.\n"
                f"To avoid syntax problems in R, use only letters, digits and underscores '_'. "
                f"Replace spaces, hyphens '-' or math symbols (+, /, *)."
            )

    # Check that samples match between the two files (counts headers only, for speed)
    try:
        counts_headers = list(_read_flexible(c_path, nrows=1).columns)
    except Exception as e:
        raise ValueError(f"Unable to read the header of the counts file: {e}")

    missing_samples = [s for s in samples_in_meta if s not in counts_headers]
    if missing_samples:
        raise ValueError(
            f"The following samples listed in the metadata are not present "
            f"as columns in the counts file:\n{', '.join(missing_samples[:5])}"
            f"{' ...' if len(missing_samples) > 5 else ''}\n\n"
            f"Check that the sample column is correct and that there are no blank spaces."
        )


# -----------------------------------------------------------------------------
# Config / esecuzione pipeline
# -----------------------------------------------------------------------------
@dataclass
class PipelineConfig:
    mode: str  # "list_groups" | "run_dea"
    metadata_path: str
    sample_col: int | str = 1
    condition_col: int | str = 2
    counts_path: Optional[str] = None
    gene_col: int | str = 1
    method: str = "RNAseq"           # "RNAseq" -> DESeq2, "Microarray" -> limma
    contrast: Optional[list] = None  # es. ["Farmaco", "Controllo"]
    pairwise_all: bool = False
    padj_cutoff: float = 0.05
    lfc_cutoff: float = 1.0
    gsea_category: str = "H"
    gsea_subcategory: Optional[str] = None

    # ORA (Over-Representation Analysis): usa la stessa categoria MSigDB della GSEA
    run_ora: bool = True
    ora_direction: str = "all"       # "all" (UP+DOWN) | "up" | "down"

    # Parametri pre-filtering
    min_reads: int = 10
    min_samples: int = 3
    
    # out_dir inizializzato a None per evitare la creazione prematura di cartelle vuote
    out_dir: Optional[str] = None

    def to_json_dict(self) -> dict:
        d = {
            "mode": self.mode,
            "metadata_path": self.metadata_path,
            "sample_col": self.sample_col,
            "condition_col": self.condition_col,
            "counts_path": self.counts_path,
            "gene_col": self.gene_col,
            "method": self.method,
            "contrast": self.contrast,
            "pairwise_all": self.pairwise_all,
            "padj_cutoff": self.padj_cutoff,
            "lfc_cutoff": self.lfc_cutoff,
            "gsea_category": self.gsea_category,
            "gsea_subcategory": self.gsea_subcategory,
            "run_ora": self.run_ora,
            "ora_direction": self.ora_direction,
            "min_reads": self.min_reads,
            "min_samples": self.min_samples,
            "out_dir": self.out_dir,
            "core_path": str(CORE_R_PATH),
        }
        return d


def _run_rscript(config: PipelineConfig, log_callback=None) -> Path:
    """Scrive il config JSON, lancia Rscript e restituisce la out_dir."""
    rscript = find_rscript()
    if rscript is None:
        raise RPipelineError(
            "Rscript not found. Install R (https://cran.r-project.org/) and "
            "make sure it is in the system PATH."
        )
    if not CORE_R_PATH.exists():
        raise RPipelineError(f"Calculation core file missing: {CORE_R_PATH}")
    if not RUNNER_R_PATH.exists():
        raise RPipelineError(f"run_pipeline.R missing: {RUNNER_R_PATH}")

    # Creazione cartella temporanea "on-demand" solo all'effettivo run
    if config.out_dir is None:
        config.out_dir = tempfile.mkdtemp(prefix="dea_run_")
    
    out_dir = Path(config.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    config_path = out_dir / f"config_{uuid.uuid4().hex[:8]}.json"
    with open(config_path, "w", encoding="utf-8") as f:
        json.dump(config.to_json_dict(), f, ensure_ascii=False, indent=2)

    cmd = [rscript, "--vanilla", str(RUNNER_R_PATH), str(config_path)]

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        bufsize=1,
        universal_newlines=True,
        env=_subprocess_env(),
    )

    full_log = []
    assert process.stdout is not None
    for line in process.stdout:
        full_log.append(line)
        if log_callback:
            log_callback(line.rstrip("\n"))
            
    # Timeout di salvaguardia (10 minuti) per non bloccare la GUI indefinitamente
    try:
        process.wait(timeout=600)
    except subprocess.TimeoutExpired:
        process.kill()
        raise RPipelineError("The R computation exceeded the maximum time (10 min) and was stopped.")

    error_json = out_dir / "error.json"
    if process.returncode != 0 or error_json.exists():
        if error_json.exists():
            try:
                with open(error_json, encoding="utf-8") as f:
                    err = json.load(f).get("error", "Unknown error in R.")
            except Exception:
                err = "".join(full_log[-30:])
        else:
            err = "".join(full_log[-30:]) or "Rscript terminated with an unknown error."
        
        # Intercettazione e traduzione umanizzata dell'errore
        friendly_error = parse_r_error(err)
        raise RPipelineError(friendly_error)

    return out_dir


def list_groups(metadata_path: str, sample_col=1, condition_col=2,
                log_callback=None) -> dict:
    """Legge il metadata e restituisce i gruppi disponibili, riusando
    import_metadata()/get_available_groups() definite nel core R (DE_scheletro_FINALE.Rmd)."""
    cfg = PipelineConfig(
        mode="list_groups",
        metadata_path=metadata_path,
        sample_col=sample_col,
        condition_col=condition_col,
    )
    out_dir = _run_rscript(cfg, log_callback=log_callback)
    groups_path = out_dir / "groups.json"
    if not groups_path.exists():
        raise RPipelineError("The R pipeline did not produce groups.json.")
    with open(groups_path, encoding="utf-8") as f:
        data = json.load(f)
    
    # Pulizia preventiva della cartella temporanea creata per list_groups
    cleanup_out_dir(str(out_dir))
    return data


@dataclass
class DEAResult:
    tag: str
    res_tbl: pd.DataFrame                 # tabella grezza DEA (GeneID, log2FoldChange, pvalue, padj, stat)
    ranked_tbl: pd.DataFrame               # tabella dei geni ordinata, gia' filtrata da core_calcolo.r
    gsea_tbl: Optional[pd.DataFrame]       # risultati GSEA grezzi (puo' essere None se fallita)
    group_high: Optional[str] = None       # gruppo "numeratore": log2FC positivo = piu' espresso qui
    group_low: Optional[str] = None        # gruppo "denominatore": log2FC negativo = piu' espresso qui
    ora_tbl: Optional[pd.DataFrame] = None  # risultati ORA grezzi (None se non eseguita / nessun risultato)
    ora_meta: Optional[dict] = None         # {available, direction, n_input_genes, origin, note}


@dataclass
class RunResults:
    pca_data: pd.DataFrame
    pca_x_label: str
    pca_y_label: str
    per_contrast: dict  # tag -> DEAResult


def run_dea(counts_path: str, metadata_path: str, method: str,
            gene_col=1, sample_col=1, condition_col=2,
            contrast: Optional[list] = None, pairwise_all: bool = False,
            padj_cutoff: float = 0.05, lfc_cutoff: float = 1.0,
            gsea_category: str = "H", gsea_subcategory: Optional[str] = None,
            min_reads: int = 10, min_samples: int = 3,  # Nuovi parametri
            out_dir: Optional[str] = None, log_callback=None,
            min_counts=0,
            run_ora: bool = True, ora_direction: str = "all") -> RunResults:
    """Esegue la pipeline completa chiamando run_pipeline.R e legge i
    risultati grezzi prodotti in DataFrame pandas."""

    if ora_direction not in ("all", "up", "down"):
        raise ValueError(f"Invalid ora_direction: '{ora_direction}' (allowed: all, up, down).")

    # 1. Validazione preventiva dei file e dei parametri
    validate_inputs(
        counts_path=counts_path,
        metadata_path=metadata_path,
        sample_col=sample_col,
        condition_col=condition_col,
        gene_col=gene_col,
        min_reads=min_reads,
        min_samples=min_samples
    )

    # 2. Setup configurazione
    cfg = PipelineConfig(
        mode="run_dea",
        metadata_path=metadata_path,
        sample_col=sample_col,
        condition_col=condition_col,
        counts_path=counts_path,
        gene_col=gene_col,
        method=method,
        contrast=contrast,
        pairwise_all=pairwise_all,
        padj_cutoff=padj_cutoff,
        lfc_cutoff=lfc_cutoff,
        gsea_category=gsea_category,
        gsea_subcategory=gsea_subcategory,
        run_ora=run_ora,
        ora_direction=ora_direction,
        min_reads=min_reads,
        min_samples=min_samples,
        out_dir=out_dir
    )
    
    out_path = _run_rscript(cfg, log_callback=log_callback)

    pca_data = pd.read_csv(out_path / "pca_data.csv")
    with open(out_path / "pca_labels.json", encoding="utf-8") as f:
        labels = json.load(f)

    pairs_path = out_path / "pairs.json"
    with open(pairs_path, encoding="utf-8") as f:
        pairs_info = json.load(f)
    tags = pairs_info["pairs"]
    if isinstance(tags, str):
        tags = [tags]

    per_contrast = {}
    for tag in tags:
        res_tbl_path = out_path / f"res_tbl_{tag}.csv"
        ranked_path = out_path / f"ranked_{tag}.csv"
        gsea_path = out_path / f"gsea_{tag}.csv"
        ora_path = out_path / f"ora_{tag}.csv"
        ora_meta_path = out_path / f"ora_meta_{tag}.json"

        if not res_tbl_path.exists():
            continue
        res_tbl = pd.read_csv(res_tbl_path)
        ranked_tbl = pd.read_csv(ranked_path) if ranked_path.exists() else pd.DataFrame()
        gsea_tbl = pd.read_csv(gsea_path) if gsea_path.exists() else None
        ora_tbl = pd.read_csv(ora_path) if ora_path.exists() else None
        ora_meta = None
        if ora_meta_path.exists():
            with open(ora_meta_path, encoding="utf-8") as f:
                ora_meta = json.load(f)

        # Determina quale gruppo e' il "numeratore" (log2FC positivo = piu' espresso li')
        # e quale e' il "denominatore" (log2FC negativo = piu' espresso li').
        # - Confronto singolo (tag == "main"): l'ordine e' esattamente quello scelto in GUI (contrast[0], contrast[1]).
        # - Confronto pairwise: run_pipeline.R costruisce il tag come "<A>_vs_<B>" con lo stesso
        #   ordine (p[1], p[2]) passato a run_differential_expression, quindi possiamo ricavarlo dal tag.
        group_high, group_low = None, None
        if tag == "main" and contrast and len(contrast) == 2:
            group_high, group_low = contrast[0], contrast[1]
        elif "_vs_" in tag:
            group_high, group_low = tag.split("_vs_", 1)

        per_contrast[tag] = DEAResult(
            tag=tag, res_tbl=res_tbl, ranked_tbl=ranked_tbl, gsea_tbl=gsea_tbl,
            group_high=group_high, group_low=group_low,
            ora_tbl=ora_tbl, ora_meta=ora_meta,
        )

    return RunResults(
        pca_data=pca_data,
        pca_x_label=labels.get("x_label", "PC1"),
        pca_y_label=labels.get("y_label", "PC2"),
        per_contrast=per_contrast,
    )


# -----------------------------------------------------------------------------
# Funzioni di supporto per l'interattivita' lato Python
# -----------------------------------------------------------------------------
def categorize_volcano(res_tbl: pd.DataFrame, padj_cutoff: float = 0.05,
                       lfc_cutoff: float = 1.0) -> pd.DataFrame:
    """Aggiunge una colonna 'status' (Upregulated/Downregulated/Not Significant)."""
    df = res_tbl.copy()

    # Logica vettorizzata (veloce, non blocca la GUI anche con decine di migliaia di geni)
    import numpy as np
    
    is_sig = df["padj"] < padj_cutoff
    is_up = df["log2FoldChange"] >= lfc_cutoff
    is_down = df["log2FoldChange"] <= -lfc_cutoff

    condizioni = [
        is_sig & is_up,
        is_sig & is_down
    ]
    scelte = ["Upregulated", "Downregulated"]

    df["status"] = np.select(condizioni, scelte, default="Not Significant")
    return df


def filter_top_pathways(gsea_tbl: pd.DataFrame, top_n: int = 20,
                        padj_cutoff: float = 0.05) -> tuple[pd.DataFrame, str]:
    """Replica la logica di selezione di plot_top_pathways() nel core R (DE_scheletro_FINALE.Rmd)."""
    df = gsea_tbl.dropna(subset=["NES"]).copy()
    n_signif = int((df["padj"] < padj_cutoff).sum())

    if n_signif == 0:
        msg = (f"No pathway has padj < {padj_cutoff:.2f}. Showing the top {top_n} "
               f"pathways by absolute |NES| (not significant, exploratory only).")
        plot_data = df.reindex(df["NES"].abs().sort_values(ascending=False).index).head(top_n)
    elif n_signif < top_n:
        msg = f"Only {n_signif} of the {top_n} requested pathways have padj < {padj_cutoff:.2f}."
        plot_data = df[df["padj"] < padj_cutoff]
    else:
        msg = f"Showing the {top_n} pathways with padj < {padj_cutoff:.2f} and highest |NES|."
        plot_data = df[df["padj"] < padj_cutoff]
        plot_data = plot_data.reindex(plot_data["NES"].abs().sort_values(ascending=False).index).head(top_n)

    plot_data = plot_data.copy()
    plot_data["pathway_clean"] = (
        plot_data["pathway"]
        .str.replace(r"^HALLMARK_", "", regex=True)
        .str.replace(r"^GOBP_", "", regex=True)
        .str.replace("_", " ", regex=False)
    )
    plot_data["Direction"] = plot_data["NES"].apply(lambda v: "Up-regulated" if v > 0 else "Down-regulated")
    plot_data = plot_data.sort_values("NES")
    return plot_data, msg


def filter_top_ora(ora_tbl: pd.DataFrame, top_n: int = 20,
                   padj_cutoff: float = 0.05) -> tuple[pd.DataFrame, str]:
    """Replica la logica di selezione di plot_ora() nel core R: pathway con
    p.adjust < cutoff (o i top N per p-value se nessuno e' significativo).
    Aggiunge le colonne di comodo per il grafico: pathway_clean, SetSize
    (K = geni del pathway nell'universo) e neglog10 (-log10 di p.adjust)."""
    import numpy as np

    df = ora_tbl.sort_values(["p.adjust", "pvalue"], kind="stable").copy()
    n_signif = int((df["p.adjust"] < padj_cutoff).sum())

    if n_signif == 0:
        msg = (f"No pathway has padj < {padj_cutoff:.2f}. Showing the top {top_n} "
               f"pathways by p-value (not significant, exploratory only).")
        plot_data = df.head(top_n)
    elif n_signif < top_n:
        msg = f"Only {n_signif} of the {top_n} requested pathways have padj < {padj_cutoff:.2f}."
        plot_data = df[df["p.adjust"] < padj_cutoff]
    else:
        msg = f"Showing the {top_n} most significant pathways with padj < {padj_cutoff:.2f}."
        plot_data = df[df["p.adjust"] < padj_cutoff].head(top_n)

    plot_data = plot_data.copy()
    plot_data["pathway_clean"] = (
        plot_data["ID"].astype(str)
        .str.replace(r"^HALLMARK_", "", regex=True)
        .str.replace(r"^GOBP_", "", regex=True)
        .str.replace("_", " ", regex=False)
    )
    plot_data["SetSize"] = pd.to_numeric(
        plot_data["BgRatio"].astype(str).str.split("/").str[0], errors="coerce")
    plot_data["neglog10"] = -np.log10(plot_data["p.adjust"].clip(lower=1e-300))
    plot_data = plot_data.sort_values("neglog10", kind="stable")
    return plot_data, msg


def cleanup_out_dir(out_dir: str):
    """Rimuove la cartella temporanea dei risultati di un run."""
    shutil.rmtree(out_dir, ignore_errors=True)
