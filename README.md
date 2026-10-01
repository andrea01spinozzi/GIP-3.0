## GIP
# What is GIP?

GIP is a free desktop software package that makes transcriptomic analysis accessible to everyone. Its goal is to let researchers, clinicians and students run gene expression analyses without programming skills and without learning complex software. There is nothing to configure, and no code or command line is needed beyond installation. GIP has two tools that work together: StandHard to build RNA-seq datasets from public data, and GIP itself to analyze them. Both run locally on your computer.

StandHard (data acquisition): this connects to the NCI Genomic Data Commons (GDC) Data Portal and lets you define your cohort with linked dropdown filters: project, primary site, sample type, disease, diagnosis, stage, metastasis, sex, vital status, therapy and outcome. The filters are populated live from the GDC API. Each dropdown shows only the values still compatible with the others, with the number of available files, and you can apply them in any order. Each filter combination becomes a group (condition) with its own number of samples. StandHard then downloads the open-access STAR - Counts files in parallel, extracts the unstranded counts, removes the summary rows (N_*) and merges everything into two analysis-ready TSV files: a counts matrix (genes × samples) and a metadata table (sample, barcode, condition, case, sample type, file ID).

GIP (analysis): starting from a count matrix and a sample metadata table, such as those produced by StandHard, it runs differential expression analysis with DESeq2 (RNA-seq) or limma (microarray). It then runs enrichment analysis: GSEA with fgsea and over-representation analysis with clusterProfiler, using MSigDB gene sets. Results are shown as interactive, customizable and exportable PCA, volcano and pathway plots.

# Why a new one?

Transcriptomic analysis is now central to biology and medicine, but doing it properly still requires a chain of skills and tools: navigating data portals, downloading and merging hundreds of files, writing R scripts for the statistics, and producing figures. Each step is a barrier for people without a bioinformatics background, and even for experts it is repetitive work. Web platforms lower the barrier but often require uploading data to a remote server and offer limited control.

GIP was developed to remove these barriers. The statistics rely on established R/Bioconductor packages, but they are hidden behind a graphical interface, so choosing a cohort, defining a comparison and obtaining publication-ready figures takes a few clicks instead of a script. The output of StandHard is directly the input of GIP, so you can go from public data to results without manual file handling. Your data stay on your machine, and the whole workflow stays transparent and reproducible.

# Requirements
Python 3.11+ and R with the required Bioconductor packages (installed with a single script). Available for Windows and macOS.
