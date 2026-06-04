#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("PROJECT_ROOT", "/platform_data/p_user/p010/phase0")
task_id <- if (length(args) >= 2) as.integer(args[[2]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

step2_dir <- file.path(project_root, "results/phase0_extension/step2_input_tables")
pkg_dir <- file.path(project_root, "results/phase0_server/nasal4_vs_other_package")
results_dir <- file.path(project_root, "results")
out_dir <- file.path(project_root, "results/phase0_extension/step3_coloc_abf")
per_locus_dir <- file.path(out_dir, "per_locus")
qc_dir <- file.path(out_dir, "qc")
dir.create(per_locus_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

min_snps <- as.integer(Sys.getenv("COLOC_MIN_SNPS", "20"))
pph4_threshold <- as.numeric(Sys.getenv("COLOC_PPH4_THRESHOLD", "0.5"))

write_tsv <- function(dt, path) {
  fwrite(dt, path, sep = "\t", na = "", quote = FALSE)
}

abs_path <- function(path) {
  path <- as.character(path)
  if (!length(path) || is.na(path) || !nzchar(path)) return("")
  if (startsWith(path, "/")) return(path)
  file.path(project_root, path)
}

fread_any <- function(path, ...) {
  path <- abs_path(path)
  if (!file.exists(path)) stop("input_file_missing: ", path)
  if (grepl("\\.gz$", path)) {
    return(fread(cmd = paste("gzip -dc", shQuote(path)), ...))
  }
  fread(path, ...)
}

first_existing <- function(paths) {
  for (path in paths) {
    if (file.exists(path)) return(path)
  }
  ""
}

split_pair <- function(pair_id) {
  parts <- strsplit(pair_id, "__", fixed = TRUE)[[1]]
  if (length(parts) >= 2) return(parts[1:2])
  c(pair_id, "")
}

standard_chr <- function(x) {
  gsub("^chr", "", as.character(x), ignore.case = TRUE)
}

fnum <- function(x) suppressWarnings(as.numeric(x))

read_optional <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(data.table())
  fread(path, showProgress = FALSE)
}

make_targets_from_sources <- function() {
  explicit <- Sys.getenv("PHASE0_COLOC_TARGET_LOCI", "")
  if (nzchar(explicit) && file.exists(explicit)) {
    targets <- fread(explicit, showProgress = FALSE)
    if (!"source" %in% names(targets)) targets[, source := "explicit_target_loci"]
    return(unique(targets, by = c("pair_id", "locus_id")))
  }

  locus_manifest_path <- file.path(pkg_dir, "phase0_v3_coloc_locus_manifest.tsv")
  locus_manifest <- fread(locus_manifest_path, showProgress = FALSE)
  locus_manifest[, `:=`(
    chr_std = standard_chr(chr),
    start_num = as.integer(start),
    end_num = as.integer(end)
  )]

  candidate_rows <- list()

  candidate_snps <- read_optional(file.path(step2_dir, "candidate_snps.tsv"))
  if (nrow(candidate_snps)) {
    candidate_snps[, `:=`(
      chr_std = standard_chr(CHR),
      bp_num = as.integer(BP),
      source_table = "step2_candidate_snps"
    )]
    candidate_rows[[length(candidate_rows) + 1]] <- candidate_snps
  }

  locus_master <- read_optional(file.path(results_dir, "integrated/locus_master_table.tsv"))
  if (nrow(locus_master)) {
    if (!"pair_id" %in% names(locus_master) && "trait_pair" %in% names(locus_master)) {
      setnames(locus_master, "trait_pair", "pair_id")
    }
    if ("chr" %in% names(locus_master) && !"CHR" %in% names(locus_master)) setnames(locus_master, "chr", "CHR")
    if ("pos" %in% names(locus_master) && !"BP" %in% names(locus_master)) setnames(locus_master, "pos", "BP")
    if ("bp" %in% names(locus_master) && !"BP" %in% names(locus_master)) setnames(locus_master, "bp", "BP")
    if ("CHR" %in% names(locus_master) && "BP" %in% names(locus_master)) {
      locus_master[, `:=`(
        chr_std = standard_chr(CHR),
        bp_num = as.integer(BP),
        source_table = "integrated_locus_master"
      )]
      candidate_rows[[length(candidate_rows) + 1]] <- locus_master
    }
  }

  placo_sig <- read_optional(file.path(results_dir, "placo_genomewide/placo_genomewide_significant_snps.tsv"))
  if (nrow(placo_sig)) {
    placo_sig[, `:=`(
      chr_std = standard_chr(CHR),
      bp_num = as.integer(BP),
      source_table = "placo_genomewide"
    )]
    candidate_rows[[length(candidate_rows) + 1]] <- placo_sig
  }

  cpassoc_sig <- read_optional(file.path(results_dir, "cpassoc_genomewide/cpassoc_genomewide_significant_snps.tsv"))
  if (nrow(cpassoc_sig)) {
    cpassoc_sig[, `:=`(
      chr_std = standard_chr(CHR),
      bp_num = as.integer(BP),
      source_table = "cpassoc_genomewide"
    )]
    candidate_rows[[length(candidate_rows) + 1]] <- cpassoc_sig
  }

  if (!length(candidate_rows)) {
    stop("No coloc target sources found. Expected candidate_snps.tsv, integrated locus_master_table.tsv, or genome-wide PLACO/CPASSOC significant SNPs.")
  }

  candidates <- rbindlist(candidate_rows, fill = TRUE)
  candidates <- candidates[nzchar(pair_id)]
  if (!"locus_id" %in% names(candidates)) candidates[, locus_id := ""]

  mapped <- list()
  for (i in seq_len(nrow(candidates))) {
    cand <- candidates[i]
    hit <- data.table()
    if (nzchar(cand$locus_id)) {
      hit <- locus_manifest[pair_id == cand$pair_id & locus_id == cand$locus_id]
    }
    if (!nrow(hit) && nzchar(cand$chr_std) && is.finite(cand$bp_num)) {
      hit <- locus_manifest[
        pair_id == cand$pair_id &
          chr_std == cand$chr_std &
          start_num <= cand$bp_num &
          end_num >= cand$bp_num
      ]
    }
    if (nrow(hit)) {
      hit <- hit[1]
      traits <- split_pair(cand$pair_id)
      mapped[[length(mapped) + 1]] <- data.table(
        pair_id = cand$pair_id,
        trait1 = traits[1],
        trait2 = traits[2],
        locus_id = hit$locus_id,
        CHR = hit$chr,
        START = hit$start,
        END = hit$end,
        region_a_file = hit$region_a_file,
        region_b_file = hit$region_b_file,
        trigger_SNP = if ("SNP" %in% names(cand)) cand$SNP else "",
        trigger_BP = if ("BP" %in% names(cand)) cand$BP else "",
        source = cand$source_table
      )
    }
  }

  if (!length(mapped)) {
    stop("No candidate SNP/locus rows could be mapped to phase0_v3_coloc_locus_manifest.tsv")
  }
  unique(rbindlist(mapped, fill = TRUE), by = c("pair_id", "locus_id"))
}

clean_region <- function(path) {
  dt <- fread_any(path, showProgress = FALSE)
  required <- c("SNP", "A1", "A2", "BETA", "SE", "P", "N", "CHR", "BP")
  missing <- setdiff(required, names(dt))
  if (length(missing)) stop("missing columns: ", paste(missing, collapse = ","))
  dt <- dt[, ..required]
  dt[, `:=`(
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P),
    N = as.numeric(N),
    CHR = as.integer(CHR),
    BP = as.integer(BP)
  )]
  dt <- unique(dt, by = "SNP")
  dt[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & is.finite(N) & N > 0 & nzchar(SNP)]
}

harmonize <- function(a, b) {
  setnames(a, c("A1", "A2", "BETA", "SE", "P", "N", "CHR", "BP"),
           paste0(c("A1", "A2", "BETA", "SE", "P", "N", "CHR", "BP"), "_1"))
  setnames(b, c("A1", "A2", "BETA", "SE", "P", "N", "CHR", "BP"),
           paste0(c("A1", "A2", "BETA", "SE", "P", "N", "CHR", "BP"), "_2"))
  m <- merge(a, b, by = "SNP")
  m <- m[(A1_1 == A1_2 & A2_1 == A2_2) | (A1_1 == A2_2 & A2_1 == A1_2)]
  if (!nrow(m)) return(m)
  flip <- m$A1_1 == m$A2_2 & m$A2_1 == m$A1_2
  if (any(flip)) {
    m$BETA_2[flip] <- -m$BETA_2[flip]
  }
  m
}

empty_result <- function(meta, status, reason, n_trait1 = 0L, n_trait2 = 0L, n_harmonized = 0L) {
  result <- data.table(
    pair_id = meta$pair_id, trait1 = meta$trait1, trait2 = meta$trait2,
    locus_id = meta$locus_id, CHR = meta$CHR, START = meta$START, END = meta$END,
    nsnps = n_harmonized, PP.H0.abf = NA_real_, PP.H1.abf = NA_real_,
    PP.H2.abf = NA_real_, PP.H3.abf = NA_real_, PP.H4.abf = NA_real_,
    lead_SNP = "", lead_min_p = NA_real_, status = status, reason = reason,
    source = meta$source
  )
  write_tsv(result, file.path(per_locus_dir, paste0(meta$locus_id, ".coloc_abf.tsv")))
  qc <- data.table(
    pair_id = meta$pair_id, trait1 = meta$trait1, trait2 = meta$trait2,
    locus_id = meta$locus_id, status = status, reason = reason,
    n_trait1_region = n_trait1, n_trait2_region = n_trait2,
    n_harmonized_snps = n_harmonized
  )
  write_tsv(qc, file.path(qc_dir, paste0(meta$locus_id, ".qc.tsv")))
}

targets <- make_targets_from_sources()
setorder(targets, pair_id, locus_id)
write_tsv(targets, file.path(out_dir, "coloc_target_loci.tsv"))

if (is.na(task_id) || task_id < 1 || task_id > nrow(targets)) {
  stop("SLURM_ARRAY_TASK_ID outside coloc target locus range")
}

meta <- as.list(targets[task_id])

tryCatch({
  a <- clean_region(meta$region_a_file)
  b <- clean_region(meta$region_b_file)
  m <- harmonize(a, b)
  if (nrow(m) < min_snps) {
    empty_result(meta, "failed", "too_few_harmonized_snps", nrow(a), nrow(b), nrow(m))
    quit(status = 0)
  }

  n1 <- median(m$N_1, na.rm = TRUE)
  n2 <- median(m$N_2, na.rm = TRUE)
  d1 <- list(beta = m$BETA_1, varbeta = m$SE_1^2, snp = m$SNP, position = m$BP_1, type = "quant", N = n1)
  d2 <- list(beta = m$BETA_2, varbeta = m$SE_2^2, snp = m$SNP, position = m$BP_2, type = "quant", N = n2)
  res <- NULL
  capture.output({ res <- coloc.abf(dataset1 = d1, dataset2 = d2) })
  summary <- as.data.table(as.data.frame(t(res$summary)))
  for (col in names(summary)) summary[[col]] <- as.numeric(summary[[col]])
  m[, min_p := pmin(P_1, P_2, na.rm = TRUE)]
  lead <- m[which.min(min_p)]
  out <- data.table(
    pair_id = meta$pair_id, trait1 = meta$trait1, trait2 = meta$trait2,
    locus_id = meta$locus_id, CHR = meta$CHR, START = meta$START, END = meta$END,
    nsnps = nrow(m),
    PP.H0.abf = summary$PP.H0.abf,
    PP.H1.abf = summary$PP.H1.abf,
    PP.H2.abf = summary$PP.H2.abf,
    PP.H3.abf = summary$PP.H3.abf,
    PP.H4.abf = summary$PP.H4.abf,
    lead_SNP = lead$SNP,
    lead_min_p = lead$min_p,
    status = "ok",
    reason = ifelse(summary$PP.H4.abf >= pph4_threshold, "coloc_H4_positive", "coloc_H4_below_threshold"),
    source = meta$source
  )
  write_tsv(out, file.path(per_locus_dir, paste0(meta$locus_id, ".coloc_abf.tsv")))
  qc <- data.table(
    pair_id = meta$pair_id, trait1 = meta$trait1, trait2 = meta$trait2,
    locus_id = meta$locus_id, status = "ok", reason = "",
    n_trait1_region = nrow(a), n_trait2_region = nrow(b),
    n_harmonized_snps = nrow(m)
  )
  write_tsv(qc, file.path(qc_dir, paste0(meta$locus_id, ".qc.tsv")))
}, error = function(e) {
  empty_result(meta, "failed", conditionMessage(e))
})
