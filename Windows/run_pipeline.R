suppressWarnings(suppressMessages({
  ok_jsonlite <- requireNamespace("jsonlite", quietly = TRUE)
}))
if (!ok_jsonlite) {
  stop("The 'jsonlite' package is not installed. Run: install.packages('jsonlite')")
}
library(jsonlite)

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# --- Utility: individua la cartella in cui si trova questo script -----------
.script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  hit <- grep(file_flag, cmd_args)
  if (length(hit) > 0) {
    return(dirname(normalizePath(sub(file_flag, "", cmd_args[hit]))))
  }
  return(getwd())
}

# --- Estrae SOLO le funzioni (blocchi 1-6) da DE_scheletro_FINALE.Rmd
.load_core_functions <- function(core_path) {
  if (!file.exists(core_path)) {
    stop("File DE_scheletro_FINALE.Rmd not found at: ", core_path)
  }
  raw_lines <- readLines(core_path, warn = FALSE, encoding = "UTF-8")

  in_chunk <- FALSE
  stop_extraction <- FALSE
  code_lines <- character(0)

  for (ln in raw_lines) {
    trimmed <- trimws(ln)

    # Il blocco 7 ("SIMULAZIONE E TEST COMPLETO") contiene solo una demo con
    # percorsi hardcoded: da qui in poi non estraiamo piu' nulla.
    if (grepl("^#\\s*7\\.", trimmed)) {
      stop_extraction <- TRUE
    }
    if (stop_extraction) next

    if (grepl("^```\\{r", trimmed)) { in_chunk <- TRUE; next }
    if (in_chunk && trimmed == "```") { in_chunk <- FALSE; next }
    if (in_chunk) code_lines <- c(code_lines, ln)
  }

  if (length(code_lines) == 0) {
    stop("Could not extract any function from DE_scheletro_FINALE.Rmd: unexpected format.")
  }

  tmp_core <- tempfile(fileext = ".R")
  writeLines(code_lines, tmp_core, useBytes = TRUE)
  source(tmp_core, echo = FALSE)
}

# --- Scrive un errore in formato JSON e termina con status != 0 -------------
.fail <- function(out_dir, msg) {
  err_path <- file.path(out_dir, "error.json")
  write(toJSON(list(error = msg), auto_unbox = TRUE), err_path)
  cat("ERROR: ", msg, "\n", sep = "")
  quit(status = 1)
}

# MAIN
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_pipeline.R <config.json>")
}
config_path <- args[1]
if (!file.exists(config_path)) {
  stop("Configuration file not found: ", config_path)
}
cfg <- fromJSON(config_path)
out_dir <- cfg$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

core_path <- cfg$core_path %||% file.path(.script_dir(), "DE_scheletro_FINALE.Rmd")

tryCatch({

  suppressWarnings(suppressMessages({
    .load_core_functions(core_path)
  }))

  # MODALITA' 1: list_groups
  # Legge solo il metadata ed estrae i gruppi disponibili, cosi' la GUI puo' farli scegliere all'utente PRIMA di lanciare l'intera pipeline.
  if (cfg$mode == "list_groups") {

    metadata <- import_metadata(cfg$metadata_path,
                                 sample_col = cfg$sample_col,
                                 condition_col = cfg$condition_col)

    groups <- get_available_groups(metadata)
    n_per_group <- as.list(table(as.character(metadata$condition)))

    out <- list(
      groups = groups,
      counts = n_per_group,
      n_samples = nrow(metadata)
    )
    write(toJSON(out, auto_unbox = TRUE, pretty = TRUE),
          file.path(out_dir, "groups.json"))
    cat("OK\n")

  # MODALITA' 2: run_dea
  # Esegue la pipeline completa (import -> DEA -> tabella -> GSEA), riusando ESATTAMENTE le funzioni di DE_scheletro_FINALE.Rmd.
  } else if (cfg$mode == "run_dea") {

    is_micro <- identical(cfg$method, "Microarray")

    counts <- import_expression(cfg$counts_path,
                                 gene_col = cfg$gene_col,
                                 is_microarray = is_micro)
    metadata <- import_metadata(cfg$metadata_path,
                                 sample_col = cfg$sample_col,
                                 condition_col = cfg$condition_col)

    # --- PCA esplorativa: estraiamo solo i dati sottostanti (p$data) e le etichette degli assi per poterla ridisegnare in modo interattivo in Python.
    p_pca <- plot_universal_pca(counts, metadata, title = "PCA dei Dati di Conteggio")
    write.csv(p_pca$data, file.path(out_dir, "pca_data.csv"), row.names = FALSE)
    write(toJSON(list(x_label = p_pca$labels$x, y_label = p_pca$labels$y),
                 auto_unbox = TRUE),
          file.path(out_dir, "pca_labels.json"))

    padj_cutoff <- cfg$padj_cutoff %||% 0.05
    gsea_category <- cfg$gsea_category %||% "H"
    gsea_subcategory <- cfg$gsea_subcategory %||% NULL
    lfc_cutoff <- cfg$lfc_cutoff %||% 1
    run_ora <- isTRUE(cfg$run_ora %||% TRUE)
    ora_direction <- cfg$ora_direction %||% "all"

    run_one <- function(contrast, tag) {
      res_tbl <- run_differential_expression(counts, metadata,
                                              method = cfg$method,
                                              contrast = contrast)
      write.csv(res_tbl, file.path(out_dir, paste0("res_tbl_", tag, ".csv")),
                row.names = FALSE)

      ranked <- generate_ranked_gene_table(res_tbl, padj_cutoff = padj_cutoff)
      write.csv(ranked, file.path(out_dir, paste0("ranked_", tag, ".csv")),
                row.names = FALSE)

      # --- ORA (Over-Representation Analysis) ---
      # Stessi cutoff (padj / log2FC) e stessa categoria MSigDB della GSEA.
      if (run_ora) {
        ora_note <- NULL
        ora_res <- tryCatch({
          withCallingHandlers(
            run_ora_analysis(res_tbl, direction = ora_direction,
                             padj_cutoff = padj_cutoff, lfc_cutoff = lfc_cutoff,
                             category = gsea_category, subcategory = gsea_subcategory),
            warning = function(w) {
              ora_note <<- conditionMessage(w)
              message(">>> [Note] ORA '", tag, "': ", conditionMessage(w))
              invokeRestart("muffleWarning")
            })
        }, error = function(e) {
          ora_note <<- conditionMessage(e)
          message(">>> [Note] ORA failed for '", tag, "': ", conditionMessage(e))
          NULL
        })

        ora_meta <- list(available = FALSE, direction = ora_direction, note = "")
        if (!is.null(ora_res) && nrow(ora_res) > 0) {
          ora_meta$available <- TRUE
          ora_meta$n_input_genes <- attr(ora_res, "n_input_genes")
          ora_meta$origin <- attr(ora_res, "origin")
          write.csv(as.data.frame(ora_res),
                    file.path(out_dir, paste0("ora_", tag, ".csv")), row.names = FALSE)
        } else {
          ora_meta$note <- ora_note %||% "No ORA results."
        }
        write(toJSON(ora_meta, auto_unbox = TRUE),
              file.path(out_dir, paste0("ora_meta_", tag, ".json")))
      }

      gsea_ok <- TRUE
      gsea_res <- tryCatch({
        run_gsea_analysis(res_tbl, category = gsea_category, subcategory = gsea_subcategory)
      }, error = function(e) {
        message(">>> [Note] GSEA failed for '", tag, "': ", conditionMessage(e))
        gsea_ok <<- FALSE
        NULL
      })
      if (gsea_ok && !is.null(gsea_res)) {
        gsea_out <- as.data.frame(gsea_res)
        # leadingEdge e' una list-column: la trasformiamo in stringa per il CSV
        if ("leadingEdge" %in% colnames(gsea_out)) {
          gsea_out$leadingEdge <- vapply(gsea_out$leadingEdge,
                                          function(x) paste(x, collapse = ";"),
                                          character(1))
        }
        write.csv(gsea_out, file.path(out_dir, paste0("gsea_", tag, ".csv")),
                  row.names = FALSE)
      }
      invisible(NULL)
    }

    if (isTRUE(cfg$pairwise_all)) {
      groups <- get_available_groups(metadata)
      if (length(groups) < 2) {
        .fail(out_dir, "At least 2 groups are needed for a pairwise comparison.")
      }
      pairs <- combn(groups, 2, simplify = FALSE)
      pair_tags <- character(0)
      for (p in pairs) {
        tag <- paste(p[1], "vs", p[2], sep = "_")
        pair_tags <- c(pair_tags, tag)
        run_one(p, tag)
      }
      write(toJSON(list(pairs = pair_tags), auto_unbox = TRUE),
            file.path(out_dir, "pairs.json"))
    } else {
      contrast <- as.character(cfg$contrast)
      run_one(contrast, "main")
      write(toJSON(list(pairs = "main"), auto_unbox = TRUE),
            file.path(out_dir, "pairs.json"))
    }

    cat("OK\n")

  } else {
    .fail(out_dir, paste0("Unrecognized mode: ", cfg$mode))
  }

}, error = function(e) {
  .fail(out_dir, conditionMessage(e))
})
