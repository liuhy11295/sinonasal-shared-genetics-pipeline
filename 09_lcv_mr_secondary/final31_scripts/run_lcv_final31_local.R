#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

base <- Sys.getenv("BASE")
if (!nzchar(base)) stop("Set BASE for the final 31-pair LCV/MR analysis.")
out_root <- Sys.getenv("OUT", file.path(base, "results_final31"))
lcv_tool_dir <- file.path(base, "lcv/tools/LCV/R")
ldscore_dir <- file.path(base, "lcv/reference_ldscores")
pair_manifest <- file.path(out_root, "00_pair_manifest_final31.tsv")
out_dir <- file.path(out_root, "LCV")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), ...) else fread(path, ...)
}

write_tsv <- function(dt, path) fwrite(dt, path, sep = "\t", quote = FALSE, na = "")

read_ldscores <- function() {
  files <- file.path(ldscore_dir, paste0("LDscore.", 1:22, ".l2.ldscore.gz"))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("missing LD score files: ", paste(missing, collapse = ", "))
  out <- rbindlist(lapply(files, function(f) {
    dt <- fread_any(f, showProgress = FALSE)
    keep <- intersect(c("CHR", "BP", "SNP", "L2"), names(dt))
    dt[, ..keep]
  }), fill = TRUE)
  out[, `:=`(CHR = as.integer(CHR), BP = as.integer(BP), L2 = as.numeric(L2))]
  out <- out[!(CHR == 6 & BP >= 25000000 & BP <= 34000000)]
  out[is.finite(L2) & L2 > 0 & nzchar(SNP)]
}

read_sumstats <- function(path) {
  dt <- fread_any(path, showProgress = FALSE)
  needed <- c("SNP", "A1", "A2", "Z", "N")
  missing <- setdiff(needed, names(dt))
  if (length(missing)) stop("missing sumstats columns: ", paste(missing, collapse = ","))
  dt <- dt[, ..needed]
  dt[, `:=`(A1 = toupper(as.character(A1)), A2 = toupper(as.character(A2)), Z = as.numeric(Z), N = as.numeric(N))]
  dt[is.finite(Z) & is.finite(N) & N > 0 & nzchar(SNP) & nchar(A1) == 1 & nchar(A2) == 1]
}

classify_lcv <- function(gcp, pval, h2z1, h2z2) {
  if (!is.finite(gcp) || !is.finite(pval) || !is.finite(h2z1) || !is.finite(h2z2)) {
    return(list(direction = "undetermined", interpretation = "insufficient_evidence"))
  }
  if (h2z1 < 4 || h2z2 < 4) {
    return(list(direction = "undetermined", interpretation = "insufficient_h2"))
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
    "weak_partial_causality"
  }
  list(direction = direction, interpretation = interpretation)
}

oldwd <- getwd()
setwd(lcv_tool_dir)
source("RunLCV.R")
setwd(oldwd)

log_path <- file.path(out_dir, "lcv_run_log.txt")
sink(log_path, split = TRUE)
cat("LCV final31 local run\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("base:", base, "\n")
cat("LCV tool:", lcv_tool_dir, "\n")
cat("pair manifest:", pair_manifest, "\n")
cat("single_thread: TRUE\n")

pairs <- fread(pair_manifest, showProgress = FALSE)
ld <- read_ldscores()
cat("LD-score SNPs:", nrow(ld), "\n")

rows <- list()
qc <- list()
for (i in seq_len(nrow(pairs))) {
  p <- pairs[i]
  row <- data.table(
    pair_id = p$pair_id, trait1 = p$trait1, trait2 = p$trait2,
    gcp = NA_real_, gcp_se = NA_real_, gcp_p = NA_real_,
    rho_estimate = NA_real_, rho_se = NA_real_,
    h2_trait1 = NA_real_, h2_trait2 = NA_real_,
    zscore = NA_real_,
    status = "failed", notes = ""
  )
  qc_row <- data.table(
    pair_id = p$pair_id, trait1 = p$trait1, trait2 = p$trait2,
    n_snp_trait1 = NA_integer_, n_snp_trait2 = NA_integer_, n_snp_merged_ld = NA_integer_,
    n_blocks = NA_integer_, h2_z_trait1 = NA_real_, h2_z_trait2 = NA_real_,
    status = "failed", notes = ""
  )
  tryCatch({
    if (!file.exists(p$trait1_lcv_file) || !file.exists(p$trait2_lcv_file)) {
      stop("missing LCV input sumstats")
    }
    d1 <- read_sumstats(p$trait1_lcv_file)
    d2 <- read_sumstats(p$trait2_lcv_file)
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
    oldwd2 <- getwd()
    setwd(lcv_tool_dir)
    res <- RunLCV(m$L2, m$Z_1, m$Z_2, no.blocks = nb, n.1 = n1, n.2 = n2)
    setwd(oldwd2)
    cls <- classify_lcv(res$gcp.pm, res$pval.gcpzero.2tailed, res$h2.zscore[1], res$h2.zscore[2])
    row[, `:=`(
      gcp = res$gcp.pm,
      gcp_se = res$gcp.pse,
      gcp_p = res$pval.gcpzero.2tailed,
      rho_estimate = res$rho.est,
      rho_se = res$rho.err,
      h2_trait1 = res$h2.zscore[1],
      h2_trait2 = res$h2.zscore[2],
      zscore = res$zscore,
      status = "success",
      notes = paste(cls$direction, cls$interpretation, sep = ";")
    )]
    qc_row[, `:=`(
      n_snp_trait1 = nrow(d1), n_snp_trait2 = nrow(d2), n_snp_merged_ld = nrow(m),
      n_blocks = nb, h2_z_trait1 = res$h2.zscore[1], h2_z_trait2 = res$h2.zscore[2],
      status = "success", notes = paste("median_n1", n1, "median_n2", n2)
    )]
  }, error = function(e) {
    row[, notes := conditionMessage(e)]
    qc_row[, notes := conditionMessage(e)]
  })
  rows[[i]] <- row
  qc[[i]] <- qc_row
  cat(i, "/", nrow(pairs), " ", p$pair_id, " ", row$status, " ", row$notes, "\n", sep = "")
}

res_dt <- rbindlist(rows, fill = TRUE)
qc_dt <- rbindlist(qc, fill = TRUE)
write_tsv(res_dt, file.path(out_dir, "lcv_pair_results.tsv"))
write_tsv(qc_dt, file.path(out_dir, "lcv_qc_summary.tsv"))
cat("Finished:", as.character(Sys.time()), "\n")
cat("LCV success:", sum(res_dt$status == "success"), " failed:", sum(res_dt$status != "success"), "\n")
sink()
