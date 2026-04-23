#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(turboGliph))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: run_gliph2_benchmark.R <input.tsv> <output.tsv> <short_name> <min_cluster_size> [work_dir]")
}

input_path <- args[[1]]
output_path <- args[[2]]
short_name <- args[[3]]
min_cluster_size <- as.integer(args[[4]])
work_dir <- if (length(args) >= 5) args[[5]] else dirname(output_path)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

input_df <- fread(input_path)
epitope <- unique(input_df$antigen.epitope)[1]

search_assignment_table <- function(obj) {
  cdr3_candidates <- c("CDR3b", "cdr3aa", "cdr3", "cdr3_beta_aa")
  cluster_candidates <- c("cluster", "cluster_id", "clusterID", "group", "group_id", "convergence_group", "Cluster")

  if (is.data.frame(obj)) {
    cn <- colnames(obj)
    cdr3_col <- intersect(cdr3_candidates, cn)
    cluster_col <- intersect(cluster_candidates, cn)
    if (length(cdr3_col) > 0 && length(cluster_col) > 0) {
      out <- obj[, c(cdr3_col[1], cluster_col[1]), with = FALSE]
      setnames(out, c("cdr3aa", "cluster"))
      return(out)
    }
  }

  if (is.list(obj)) {
    for (item in obj) {
      res <- tryCatch(search_assignment_table(item), error = function(e) NULL)
      if (!is.null(res)) {
        return(res)
      }
    }
  }

  NULL
}

search_output_files <- function(root_dir) {
  files <- list.files(root_dir, recursive = TRUE, full.names = TRUE)
  files <- files[grepl("\\.(txt|tsv|csv)$", files, ignore.case = TRUE)]

  for (path in files) {
    sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
    tbl <- tryCatch(fread(path, sep = sep), error = function(e) NULL)
    if (is.null(tbl)) {
      next
    }
    res <- tryCatch(search_assignment_table(tbl), error = function(e) NULL)
    if (!is.null(res)) {
      return(res)
    }
  }

  NULL
}

gliph_input <- data.frame(
  CDR3b = input_df$cdr3aa,
  TRBV = input_df$`v.segm`,
  stringsAsFactors = FALSE
)

result <- gliph2(gliph_input, result_folder = work_dir)
assignments <- search_assignment_table(result)
if (is.null(assignments)) {
  assignments <- search_output_files(work_dir)
}
if (is.null(assignments)) {
  saveRDS(result, file.path(work_dir, paste0(short_name, "_gliph2_debug.rds")))
  stop("Could not extract GLIPH2 cluster assignments; saved debug object in work_dir")
}

assignments <- unique(assignments[cluster != "" & !is.na(cluster)])
cluster_sizes <- assignments[, .N, by = cluster][N >= min_cluster_size]
assignments <- assignments[cluster %in% cluster_sizes$cluster]

merged <- merge(input_df, assignments, by = "cdr3aa", allow.cartesian = TRUE)
if (nrow(merged) == 0) {
  fwrite(data.frame(
    gene = character(),
    cdr3aa = character(),
    v.segm = character(),
    j.segm = character(),
    cid = character(),
    antigen.epitope = character()
  ), output_path, sep = "\t")
  quit(save = "no")
}

merged[, cid := paste0("H.B.", epitope, ".", cluster)]
out <- unique(merged[, .(
  gene,
  cdr3aa,
  `v.segm`,
  `j.segm`,
  cid,
  antigen.epitope
)])

fwrite(out, output_path, sep = "\t")
