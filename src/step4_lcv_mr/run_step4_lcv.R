#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("PROJECT_ROOT", "/platform_data/p_user/p010/phase0")

step2_dir <- file.path(project_root, "results/phase0_extension/step2_input_tables")
out_dir <- file.path(project_root, "results/phase0_extension/step4_lcv_mr")
log_dir <- file.path(out_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

lcv_dir <- file.path(project_root, "results/phase0_extension/step0_1/tools/LCV/R")
ldscore_dir <- file.path(project_root, "data/reference/ldsc/1000G_Phase3_ldscores/LDscore")
threads <- as.integer(Sys.getenv("LCV_THREADS", Sys.getenv("SLURM_NTASKS", "1")))
if (!is.finite(threads) || threads < 1) threads <- 1

fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("gzip -dc", shQuote(path)), ...)
  } else {
    fread(path, ...)
  }
}

read_ldscores <- function() {
  files <- file.path(ldscore_dir, paste0("LDscore.", 1:22, ".l2.ldscore.gz"))
  out <- rbindlist(lapply(files, function(f) {
    dt <- fread_any(f, showProgress = FALSE)
    keep <- intersect(c("CHR", "BP", "SNP", "L2"), names(dt))
    dt[, ..keep]
  }), fill = TRUE)
  out[, CHR := as.integer(CHR)]
  out[, BP := as.integer(BP)]
  out <- out[!(CHR == 6 & BP >= 25000000 & BP <= 34000000)]
  out[is.finite(L2) & L2 > 0 & nzchar(SNP)]
}

read_sumstats <- function(path) {
  dt <- fread_any(path, showProgress = FALSE)
  needed <- c("SNP", "A1", "A2", "Z", "N")
  missing <- setdiff(needed, names(dt))
  if (length(missing)) stop("missing sumstats columns: ", paste(missing, collapse = ","))
  dt <- dt[, ..needed]
  dt[, Z := as.numeric(Z)]
  dt[, N := as.numeric(N)]
  dt[is.finite(Z) & is.finite(N) & N > 0 & nzchar(SNP)]
}

classify_lcv <- function(gcp, pval, h2z1, h2z2) {
  if (!is.finite(gcp) || !is.finite(pval) || !is.finite(h2z1) || !is.finite(h2z2)) {
    return(list(direction = "undetermined", interpretation = "insufficient_evidence"))
  }
  if (h2z1 < 4 || h2z2 < 4) {
    return(list(direction = "undetermined", interpretation = "insufficient_evidence"))
  }
  direction <- if (pval < 0.05 && gcp > 0) {
    "trait1_to_trait2"
  } else if (pval < 0.05 && gcp < 0) {
    "trait2_to_trait1"
  } else {
    "undetermined"
  }
  interpretation <- if (pval < 0.05 && abs(gcp) >= 0.6) {
    "strong_partial_causality"
  } else if (pval < 0.05 && abs(gcp) >= 0.3) {
    "moderate_partial_causality"
  } else if (pval >= 0.05) {
    "shared_genetics_only"
  } else {
    "insufficient_evidence"
  }
  list(direction = direction, interpretation = interpretation)
}

pairs <- fread(file.path(step2_dir, "priority_pairs.tsv"), showProgress = FALSE)
pair_manifest <- fread(file.path(step2_dir, "pair_manifest.tsv"), showProgress = FALSE)
pairs <- merge(
  pairs,
  pair_manifest[, .(pair_id, trait1_munged, trait2_munged)],
  by = "pair_id",
  all.x = TRUE
)

oldwd <- getwd()
setwd(lcv_dir)
source("RunLCV.R")

ld <- read_ldscores()
run_one_pair <- function(i) {
  p <- pairs[i]
  row <- data.table(
    pair_id = p$pair_id,
    trait1 = p$trait1,
    trait2 = p$trait2,
    gcp = NA_real_,
    gcp_se = NA_real_,
    p_value = NA_real_,
    direction = "undetermined",
    interpretation = "insufficient_evidence"
  )
  err <- ""
  tryCatch({
    if (!file.exists(p$trait1_munged) || !file.exists(p$trait2_munged)) {
      stop("missing munged sumstats")
    }
    d1 <- read_sumstats(p$trait1_munged)
    d2 <- read_sumstats(p$trait2_munged)
    setnames(d1, c("A1", "A2", "Z", "N"), c("A1_1", "A2_1", "Z_1", "N_1"))
    setnames(d2, c("A1", "A2", "Z", "N"), c("A1_2", "A2_2", "Z_2", "N_2"))
    m <- merge(ld, d1, by = "SNP")
    m <- merge(m, d2, by = "SNP")
    m <- m[(A1_1 == A1_2 & A2_1 == A2_2) | (A1_1 == A2_2 & A2_1 == A1_2)]
    flip <- m$A1_1 == m$A2_2 & m$A2_1 == m$A1_2
    if (any(flip)) m$Z_2[flip] <- -m$Z_2[flip]
    m <- m[order(CHR, BP)]
    if (nrow(m) < 5000) stop("too few SNPs after LD score merge")
    nb <- max(20, min(100, floor(nrow(m) / 5000)))
    n1 <- median(m$N_1, na.rm = TRUE)
    n2 <- median(m$N_2, na.rm = TRUE)
    res <- RunLCV(m$L2, m$Z_1, m$Z_2, no.blocks = nb, n.1 = n1, n.2 = n2)
    cls <- classify_lcv(res$gcp.pm, res$pval.gcpzero.2tailed, res$h2.zscore[1], res$h2.zscore[2])
    row[, `:=`(
      gcp = res$gcp.pm,
      gcp_se = res$gcp.pse,
      p_value = res$pval.gcpzero.2tailed,
      direction = cls$direction,
      interpretation = cls$interpretation
    )]
  }, error = function(e) {
    row[, interpretation := "insufficient_evidence"]
    err <<- paste(p$pair_id, conditionMessage(e))
  })
  list(row = row, error = err)
}

idx <- seq_len(nrow(pairs))
workers <- min(threads, length(idx))
cat("LCV pairs=", length(idx), " workers=", workers, "\n", sep = "")
if (file.exists(file.path(log_dir, "lcv_failures.log"))) {
  file.remove(file.path(log_dir, "lcv_failures.log"))
}
if (workers > 1) {
  results <- mclapply(idx, run_one_pair, mc.cores = workers, mc.preschedule = FALSE)
} else {
  results <- lapply(idx, run_one_pair)
}

rows <- lapply(results, `[[`, "row")
errors <- Filter(nzchar, vapply(results, `[[`, character(1), "error"))
if (length(errors)) {
  writeLines(errors, file.path(log_dir, "lcv_failures.log"))
}
out <- rbindlist(rows, fill = TRUE)
fwrite(out, file.path(out_dir, "lcv_results.tsv"), sep = "\t", quote = FALSE, na = "")
cat("LCV wrote rows=", nrow(out), " failures=", length(errors), "\n", sep = "")
setwd(oldwd)
