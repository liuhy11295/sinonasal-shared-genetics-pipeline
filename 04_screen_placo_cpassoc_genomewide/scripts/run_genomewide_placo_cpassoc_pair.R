#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(MASS)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("PROJECT_ROOT", "/platform_data/p_user/p010/phase0")
pair_id_arg <- if (length(args) >= 2) args[[2]] else Sys.getenv("PAIR_ID", "")
task_id <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")))
threads <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("SLURM_NTASKS", "1"))))
if (!is.finite(threads) || threads < 1) threads <- 1L
data.table::setDTthreads(threads)

project_root <- normalizePath(project_root, mustWork = TRUE)
results_dir <- file.path(project_root, "results")
pair_manifest_file <- file.path(results_dir, "phase0_extension/step2_input_tables/pair_manifest.tsv")
evidence_manifest_file <- file.path(results_dir, "step8_mtag/evidence_pair_manifest.tsv")
placo_dir <- file.path(results_dir, "placo_genomewide")
cpassoc_dir <- file.path(results_dir, "cpassoc_genomewide")
scratch_dir <- file.path(results_dir, "scratch", "placo_cpassoc_genomewide", Sys.getenv("SLURM_JOB_ID", "local"))
dir.create(file.path(placo_dir, "significant"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cpassoc_dir, "significant"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(placo_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cpassoc_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)

placo_src <- file.path(results_dir, "phase0_extension/step0_1/tools/PLACO/PLACO_v0.2.0.R")
cpassoc_zip <- file.path(results_dir, "phase0_extension/step0_1/tools/CPASSOC/CPASSOC.zip")
gamma_n <- suppressWarnings(as.integer(Sys.getenv("CPASSOC_GAMMA_N", "200000")))
if (!is.finite(gamma_n) || gamma_n < 10000) gamma_n <- 200000L

suppressMessages(source(placo_src))
suppressMessages(source(unz(cpassoc_zip, "CPASSOC/FunctionSet.R")))

fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), ...)
  else fread(path, ...)
}

write_tsv <- function(dt, path) fwrite(dt, path, sep = "\t", quote = FALSE, na = "")

bh <- function(p) {
  p.adjust(p, method = "BH")
}

clean_p <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p) | p <= 0 | p > 1] <- NA_real_
  p
}

is_ambiguous <- function(a1, a2) {
  pair <- paste0(a1, a2)
  pair %in% c("AT", "TA", "CG", "GC")
}

read_munged <- function(path) {
  dt <- fread_any(path, showProgress = FALSE)
  if ("A1" %in% names(dt) && !"EA" %in% names(dt)) setnames(dt, "A1", "EA")
  if ("A2" %in% names(dt) && !"OA" %in% names(dt)) setnames(dt, "A2", "OA")
  required <- c("SNP", "EA", "OA", "Z", "N")
  missing <- setdiff(required, names(dt))
  if (length(missing)) stop("missing munged GWAS columns in ", path, ": ", paste(missing, collapse = ","))
  keep <- unique(c("SNP", "EA", "OA", "Z", "N", intersect(c("CHR", "BP", "P"), names(dt))))
  dt <- dt[, ..keep]
  dt[, `:=`(
    SNP = as.character(SNP),
    EA = toupper(as.character(EA)),
    OA = toupper(as.character(OA)),
    Z = as.numeric(Z),
    N = as.numeric(N)
  )]
  if (!"P" %in% names(dt)) dt[, P := 2 * pnorm(abs(Z), lower.tail = FALSE)]
  dt[, P := clean_p(P)]
  if (!"CHR" %in% names(dt)) dt[, CHR := NA_integer_]
  if (!"BP" %in% names(dt)) dt[, BP := NA_integer_]
  dt <- dt[nzchar(SNP) & nzchar(EA) & nzchar(OA) & is.finite(Z) & is.finite(N) & N > 0 & !is.na(P)]
  dt <- dt[!is_ambiguous(EA, OA)]
  setorder(dt, SNP)
  unique(dt, by = "SNP")
}

read_annotation <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(data.table(SNP=character(), CHR=integer(), BP=integer()))
  dt <- fread_any(path, select = intersect(c("SNP", "CHR", "BP"), names(fread_any(path, nrows = 0))), showProgress = FALSE)
  if (!all(c("SNP", "CHR", "BP") %in% names(dt))) return(data.table(SNP=character(), CHR=integer(), BP=integer()))
  dt[, `:=`(SNP = as.character(SNP), CHR = as.integer(CHR), BP = as.integer(BP))]
  unique(dt[nzchar(SNP)], by = "SNP")
}

qc_empty <- function(pair_id, trait1, trait2, g1, g2, status, reason) {
  data.table(
    pair_id = pair_id, trait1 = trait1, trait2 = trait2,
    gwas1_file = g1, gwas2_file = g2,
    n_snps_gwas1 = NA_integer_, n_snps_gwas2 = NA_integer_,
    n_snps_overlap = NA_integer_, n_snps_after_qc = NA_integer_,
    n_sig_5e8 = NA_integer_, n_fdr_005 = NA_integer_,
    min_p = NA_real_, median_p = NA_real_, lambda_gc_if_available = NA_real_,
    status = status, failed_reason = reason
  )
}

pair_manifest <- fread(pair_manifest_file, showProgress = FALSE)
evidence <- fread(evidence_manifest_file, showProgress = FALSE)
for (col in intersect(c("pair_id", "trait1", "trait2", "nasal_trait", "partner_trait"), names(evidence))) {
  evidence <- evidence[!grepl("CHRONIC_RHINITIS_PANUKB_J31", get(col), fixed = TRUE)]
}

if (nzchar(pair_id_arg)) {
  target <- evidence[pair_id == pair_id_arg][1]
} else {
  if (!is.finite(task_id) || task_id < 1 || task_id > nrow(evidence)) stop("Invalid SLURM_ARRAY_TASK_ID and no PAIR_ID supplied")
  target <- evidence[task_id]
}
if (nrow(target) != 1 || is.na(target$pair_id)) stop("Unable to resolve target pair")

pair_id <- as.character(target[["pair_id"]][1])
manifest_idx <- which(as.character(pair_manifest[["pair_id"]]) == pair_id)
if (length(manifest_idx) < 1) stop("Unable to resolve pair in pair_manifest.tsv: ", pair_id)
manifest <- pair_manifest[manifest_idx[1], ]
if (nrow(manifest) != 1 || is.na(manifest$pair_id)) stop("Unable to resolve pair in pair_manifest.tsv: ", pair_id)
trait1 <- if ("trait1" %in% names(target)) target$trait1 else manifest$trait1
trait2 <- if ("trait2" %in% names(target)) target$trait2 else manifest$trait2
gwas1 <- manifest$trait1_munged
gwas2 <- manifest$trait2_munged
std1 <- if ("trait1_standard_gwas" %in% names(manifest)) manifest$trait1_standard_gwas else ""
std2 <- if ("trait2_standard_gwas" %in% names(manifest)) manifest$trait2_standard_gwas else ""
safe_pair <- gsub("[^A-Za-z0-9_.-]+", "_", pair_id)

if (identical(Sys.getenv("CHECK_PAIR_ONLY", "0"), "1")) {
  cat("pair_id=", pair_id, "\n", sep = "")
  cat("trait1=", trait1, "\n", sep = "")
  cat("trait2=", trait2, "\n", sep = "")
  cat("trait1_munged=", gwas1, "\n", sep = "")
  cat("trait2_munged=", gwas2, "\n", sep = "")
  quit(status = 0)
}

placo_sig_path <- file.path(placo_dir, "significant", paste0(safe_pair, ".placo.sig.tsv"))
cpassoc_sig_path <- file.path(cpassoc_dir, "significant", paste0(safe_pair, ".cpassoc.sig.tsv"))
qc_placo_path <- file.path(placo_dir, "logs", paste0(safe_pair, ".placo.qc.tsv"))
qc_cpassoc_path <- file.path(cpassoc_dir, "logs", paste0(safe_pair, ".cpassoc.qc.tsv"))

start_time <- proc.time()[["elapsed"]]

empty_placo_sig <- function() {
  data.table(
    pair_id=character(), trait1=character(), trait2=character(), SNP=character(),
    CHR=integer(), BP=integer(), EA=character(), OA=character(),
    Z_trait1=numeric(), Z_trait2=numeric(), P_trait1=numeric(), P_trait2=numeric(),
    PLACO_stat=numeric(), PLACO_p=numeric(), PLACO_bh_q=numeric(),
    significance_class=character()
  )
}

empty_cpassoc_sig <- function() {
  data.table(
    pair_id=character(), trait1=character(), trait2=character(), SNP=character(),
    CHR=integer(), BP=integer(), EA=character(), OA=character(),
    P_trait1=numeric(), P_trait2=numeric(),
    SHet=numeric(), SHet_p=numeric(), SHet_bh_q=numeric(),
    SHom=numeric(), SHom_p=numeric(), SHom_bh_q=numeric(),
    significance_class=character()
  )
}

tryCatch({
  if (is.na(gwas1) || is.na(gwas2) || !file.exists(gwas1) || !file.exists(gwas2)) stop("missing munged GWAS file")
  g1 <- read_munged(gwas1)
  g2 <- read_munged(gwas2)
  n1 <- nrow(g1); n2 <- nrow(g2)
  setnames(g1, c("EA", "OA", "Z", "N", "P", "CHR", "BP"), c("EA1", "OA1", "Z1", "N1", "P1", "CHR1", "BP1"))
  setnames(g2, c("EA", "OA", "Z", "N", "P", "CHR", "BP"), c("EA2", "OA2", "Z2", "N2", "P2", "CHR2", "BP2"))
  m <- merge(g1, g2, by = "SNP")
  n_overlap <- nrow(m)
  if (n_overlap < 10000) stop("SNP overlap below 10000")
  same <- m$EA1 == m$EA2 & m$OA1 == m$OA2
  flip <- m$EA1 == m$OA2 & m$OA1 == m$EA2
  m <- m[same | flip]
  if (nrow(m) < 10000) stop("SNP overlap after allele harmonization below 10000")
  flip2 <- m$EA1 == m$OA2 & m$OA1 == m$EA2
  if (any(flip2)) m$Z2[flip2] <- -m$Z2[flip2]
  m[, `:=`(
    CHR = fifelse(!is.na(CHR1), CHR1, CHR2),
    BP = fifelse(!is.na(BP1), BP1, BP2)
  )]
  if (all(is.na(m$CHR)) || all(is.na(m$BP))) {
    annot <- rbindlist(list(read_annotation(std1), read_annotation(std2)), fill = TRUE)
    annot <- unique(annot[!is.na(CHR) & !is.na(BP)], by = "SNP")
    m <- merge(m, annot, by = "SNP", all.x = TRUE, suffixes = c("", "_ann"))
    m[is.na(CHR), CHR := CHR_ann]
    m[is.na(BP), BP := BP_ann]
    m[, c("CHR_ann", "BP_ann") := NULL]
  }
  m <- m[!is.na(CHR) & !is.na(BP)]
  m <- m[is.finite(Z1) & is.finite(Z2) & !is.na(P1) & !is.na(P2)]
  n_after <- nrow(m)
  if (n_after < 10000) stop("SNP count after QC below 10000")

  zmat <- as.matrix(m[, .(Z1, Z2)])
  pmat <- as.matrix(m[, .(P1, P2)])
  sample_size <- c(median(m$N1, na.rm = TRUE), median(m$N2, na.rm = TRUE))
  cor_scalar <- tryCatch(cor.pearson(zmat, pmat, p.threshold = 1e-4, returnMatrix = FALSE), error = function(e) 0)
  if (!is.finite(cor_scalar)) cor_scalar <- 0
  cor_mat <- matrix(c(1, cor_scalar, cor_scalar, 1), 2, 2)
  varz <- tryCatch(var.placo(zmat, pmat, p.threshold = 1e-4), error = function(e) c(1, 1))

  placo_one <- function(i) {
    res <- tryCatch(placo.plus(Z = zmat[i,], VarZ = varz, CorZ = cor_scalar), error = function(e) list(T.placo.plus = NA_real_, p.placo.plus = NA_real_))
    c(stat = res$T.placo.plus, p = res$p.placo.plus)
  }
  n_workers <- max(1L, min(threads, nrow(m)))
  if (n_workers > 1L) {
    idx <- seq_len(nrow(m))
    chunks <- split(idx, cut(idx, breaks = n_workers, labels = FALSE))
    placo_chunks <- parallel::mclapply(chunks, function(chunk) {
      do.call(rbind, lapply(chunk, placo_one))
    }, mc.cores = n_workers)
    placo_mat <- do.call(rbind, placo_chunks)
  } else {
    placo_mat <- do.call(rbind, lapply(seq_len(nrow(m)), placo_one))
  }
  placo_p <- as.numeric(placo_mat[, "p"])
  placo_q <- bh(placo_p)
  placo_keep <- which(is.finite(placo_p) & placo_p < 5e-8)
  placo_out <- empty_placo_sig()
  if (length(placo_keep)) {
    placo_out <- data.table(
      pair_id = pair_id, trait1 = trait1, trait2 = trait2,
      SNP = m$SNP[placo_keep], CHR = m$CHR[placo_keep], BP = m$BP[placo_keep],
      EA = m$EA1[placo_keep], OA = m$OA1[placo_keep],
      Z_trait1 = m$Z1[placo_keep], Z_trait2 = m$Z2[placo_keep],
      P_trait1 = m$P1[placo_keep], P_trait2 = m$P2[placo_keep],
      PLACO_stat = as.numeric(placo_mat[placo_keep, "stat"]),
      PLACO_p = placo_p[placo_keep], PLACO_bh_q = placo_q[placo_keep],
      significance_class = "genome_wide_5e-8"
    )
  }
  write_tsv(placo_out, placo_sig_path)

  shom <- tryCatch(SHom(zmat, SampleSize = sample_size, CorrMatrix = cor_mat), error = function(e) rep(NA_real_, nrow(m)))
  shet <- tryCatch(SHet(zmat, SampleSize = sample_size, CorrMatrix = cor_mat), error = function(e) rep(NA_real_, nrow(m)))
  gamma <- tryCatch(EstimateGamma(N = gamma_n, SampleSize = sample_size, CorrMatrix = cor_mat), error = function(e) c(NA_real_, NA_real_, NA_real_))
  shet_p <- rep(NA_real_, length(shet))
  if (all(is.finite(gamma))) {
    shet_p <- ifelse(shet > gamma[3], pgamma(shet - gamma[3], shape = gamma[1], scale = gamma[2], lower.tail = FALSE), 1)
  }
  shom_p <- pchisq(shom, df = 1, lower.tail = FALSE)
  shet_q <- bh(shet_p)
  shom_q <- bh(shom_p)
  cpassoc_keep <- which(
    (is.finite(shet_p) & shet_p < 5e-8) | (is.finite(shom_p) & shom_p < 5e-8)
  )
  cpassoc_out <- empty_cpassoc_sig()
  if (length(cpassoc_keep)) {
    cpassoc_out <- data.table(
      pair_id = pair_id, trait1 = trait1, trait2 = trait2,
      SNP = m$SNP[cpassoc_keep], CHR = m$CHR[cpassoc_keep], BP = m$BP[cpassoc_keep],
      EA = m$EA1[cpassoc_keep], OA = m$OA1[cpassoc_keep],
      P_trait1 = m$P1[cpassoc_keep], P_trait2 = m$P2[cpassoc_keep],
      SHet = shet[cpassoc_keep], SHet_p = shet_p[cpassoc_keep], SHet_bh_q = shet_q[cpassoc_keep],
      SHom = shom[cpassoc_keep], SHom_p = shom_p[cpassoc_keep], SHom_bh_q = shom_q[cpassoc_keep],
      significance_class = "genome_wide_5e-8"
    )
  }
  write_tsv(cpassoc_out, cpassoc_sig_path)

  lambda_placo <- median(qchisq(1 - placo_p[is.finite(placo_p)], df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
  qcp <- data.table(
    pair_id = pair_id, trait1 = trait1, trait2 = trait2, gwas1_file = gwas1, gwas2_file = gwas2,
    n_snps_gwas1 = n1, n_snps_gwas2 = n2, n_snps_overlap = n_overlap, n_snps_after_qc = n_after,
    n_sig_5e8 = sum(is.finite(placo_p) & placo_p < 5e-8), n_fdr_005 = sum(is.finite(placo_q) & placo_q < 0.05),
    min_p = suppressWarnings(min(placo_p, na.rm = TRUE)), median_p = suppressWarnings(median(placo_p, na.rm = TRUE)),
    lambda_gc_if_available = lambda_placo, status = "PASS", failed_reason = ""
  )
  qcc <- data.table(
    pair_id = pair_id, trait1 = trait1, trait2 = trait2, gwas1_file = gwas1, gwas2_file = gwas2,
    n_snps_gwas1 = n1, n_snps_gwas2 = n2, n_snps_overlap = n_overlap, n_snps_after_qc = n_after,
    n_sig_5e8 = sum((is.finite(shet_p) & shet_p < 5e-8) | (is.finite(shom_p) & shom_p < 5e-8)),
    n_fdr_005 = sum((is.finite(shet_q) & shet_q < 0.05) | (is.finite(shom_q) & shom_q < 0.05)),
    min_p = suppressWarnings(min(c(shet_p, shom_p), na.rm = TRUE)),
    median_p = suppressWarnings(median(c(shet_p, shom_p), na.rm = TRUE)),
    lambda_gc_if_available = NA_real_, status = "PASS", failed_reason = ""
  )
  write_tsv(qcp, qc_placo_path)
  write_tsv(qcc, qc_cpassoc_path)
}, error = function(e) {
  write_tsv(empty_placo_sig(), placo_sig_path)
  write_tsv(empty_cpassoc_sig(), cpassoc_sig_path)
  write_tsv(qc_empty(pair_id, trait1, trait2, gwas1, gwas2, "FAIL", conditionMessage(e)), qc_placo_path)
  write_tsv(qc_empty(pair_id, trait1, trait2, gwas1, gwas2, "FAIL", conditionMessage(e)), qc_cpassoc_path)
  message("FAILED: ", conditionMessage(e))
})

unlink(scratch_dir, recursive = TRUE, force = TRUE)
elapsed <- proc.time()[["elapsed"]] - start_time
cat("pair_id=", pair_id, " elapsed_seconds=", round(elapsed, 2), " threads=", threads, "\n", sep = "")
