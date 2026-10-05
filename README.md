# GIP

## What is GIP?

GIP is a free desktop software package that makes transcriptomic analysis accessible to everyone. Its goal is to let researchers, clinicians and students run gene expression analyses without programming skills and without learning complex software. There is nothing to configure, and no code or command line is needed beyond installation.

GIP has two tools that work together: **StandHard** to build RNA-seq datasets from public data, and **GIP** itself to analyze them. Both run locally on your computer.

**StandHard (data acquisition):** this connects to the NCI Genomic Data Commons (GDC) Data Portal and lets you define your cohort with linked dropdown filters: project, primary site, sample type, disease, diagnosis, stage, metastasis, sex, vital status, therapy and outcome. The filters are populated live from the GDC API. Each dropdown shows only the values still compatible with the others, with the number of available files, and you can apply them in any order. Each filter combination becomes a group (condition) with its own number of samples. StandHard then downloads the open-access STAR - Counts files in parallel, extracts the unstranded counts, removes the summary rows (`N_*`) and merges everything into two analysis-ready TSV files: a counts matrix (genes × samples) and a metadata table (sample, barcode, condition, case, sample type, file ID).

**GIP (analysis):** starting from a count matrix and a sample metadata table, such as those produced by StandHard, it runs differential expression analysis with DESeq2 (RNA-seq) or limma (microarray). It then runs enrichment analysis: GSEA with fgsea and over-representation analysis with clusterProfiler, using MSigDB gene sets. Results are shown as interactive, customizable and exportable PCA, volcano and pathway plots.

## Why we did that

Transcriptomic analysis is now central to biology and medicine, but doing it properly still requires a chain of skills and tools: navigating data portals, downloading and merging hundreds of files, writing R scripts for the statistics, and producing figures. Each step is a barrier for people without a bioinformatics background, and even for experts it is repetitive work. Web platforms lower the barrier but often require uploading data to a remote server and offer limited control.

GIP was developed to remove these barriers. The statistics rely on established R/Bioconductor packages, but they are hidden behind a graphical interface, so choosing a cohort, defining a comparison and obtaining publication-ready figures takes a few clicks instead of a script. The output of StandHard is directly the input of GIP, so you can go from public data to results without manual file handling. Your data stay on your machine, and the whole workflow stays transparent and reproducible.

## Requirements

### System
- Windows or macOS (Linux should work but is untested)
- Internet connection: always needed by StandHard (GDC API); needed by GIP only for the one-time R package installation
- A few GB of free disk space; 8 GB of RAM or more is recommended for large cohorts

### Python
Python 3.11 or later (on macOS use the installer from [python.org](https://www.python.org/downloads/), which includes Tk). Tkinter is included in standard Python installers.

| Package | Used by | Required | Purpose |
|---|---|---|---|
| `requests` | StandHard | Yes | Communication with the GDC API |
| `pandas` | StandHard, GIP | Yes | Reading, merging and writing tables |
| `numpy` | GIP | Yes | Numerical operations |
| `matplotlib` | GIP | Yes | PCA, volcano and pathway plots |
| `mplcursors` | GIP | Optional | Interactive tooltips on plots |
| `tkinterdnd2` | GIP | Optional | Drag and drop of input files |

```bash
pip install requests pandas numpy matplotlib mplcursors tkinterdnd2
```

### R (GIP only)
GIP performs the statistics in R, so it needs [R](https://cran.r-project.org/) (4.x recommended; on macOS choose the arm64 build for Apple Silicon and x86_64 for Intel) and `Rscript` reachable by the program. GIP looks for it in the system PATH and in the standard installation folders; you can force a path with the `DEA_RSCRIPT` environment variable.

| R package | Source | Purpose |
|---|---|---|
| `tidyverse` | CRAN | Data handling and plotting |
| `readxl` | CRAN | Reading Excel input files |
| `jsonlite` | CRAN | Exchanging settings and results with the interface |
| `msigdbr` | CRAN | MSigDB gene sets |
| `DESeq2` | Bioconductor | Differential expression for RNA-seq counts |
| `limma` | Bioconductor | Differential expression for microarray data |
| `fgsea` | Bioconductor | Gene set enrichment analysis (GSEA) |
| `clusterProfiler` | Bioconductor | Over-representation analysis (ORA) |
| `org.Hs.eg.db` | Bioconductor | Human gene annotation |
| `BiocParallel` | Bioconductor | Parallel computation |

Install them all once with:

```bash
Rscript install_r_packages.R
```

The first installation can take 10-30 minutes. The current version is set up for **human** data. On macOS you may need the Xcode command line tools (`xcode-select --install`), on Windows Rtools, if a package has to be compiled.

## Quick start

```bash
# 1. Build a dataset from the GDC Data Portal
StandHard.py

# 2. Analyze it
main.py
```

StandHard produces `counts_matrix_unito.tsv` and `metadata_unito.tsv`, which can be loaded directly into GIP.

<img width="918" height="918" alt="Resume pics" src="https://github.com/user-attachments/assets/0a3bba2e-9ba1-4949-aa04-75dc15a5e761" />
```bash
Example results
```
<img width="1917" height="982" alt="SH" src="https://github.com/user-attachments/assets/8ba15438-229a-4d15-b819-2907cd2ed89d" />
```bash
StandHard
```
<img width="1917" height="977" alt="GIP" src="https://github.com/user-attachments/assets/94add01c-13d0-4175-8b48-fb6016f1dddf" />
```bash
GIP
```




