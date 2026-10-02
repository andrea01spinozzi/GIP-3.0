options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("tidyverse", "readxl", "jsonlite", "msigdbr",
          "DESeq2", "fgsea", "clusterProfiler", "org.Hs.eg.db", "limma", "BiocParallel")
BiocManager::install(pkgs, ask = FALSE, update = FALSE)
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) stop("Not installed: ", paste(missing, collapse = ", ")) else cat("All R packages are installed.\n")
