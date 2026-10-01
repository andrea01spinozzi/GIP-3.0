# Installa i pacchetti R richiesti da DEA Explorer (da eseguire UNA volta):
#   Rscript install_r_packages.R
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("tidyverse", "readxl", "jsonlite", "msigdbr",
          "DESeq2", "fgsea", "clusterProfiler", "org.Hs.eg.db", "limma", "BiocParallel")
BiocManager::install(pkgs, ask = FALSE, update = FALSE)
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) stop("Non installati: ", paste(missing, collapse = ", ")) else cat("Tutti i pacchetti R sono installati.\n")
