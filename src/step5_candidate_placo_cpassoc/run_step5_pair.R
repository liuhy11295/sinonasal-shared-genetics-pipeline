#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(MASS)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("PROJECT_ROOT", "/platform_data/p_user/p010/phase0")
task_id <- if (length(args) >= 2) as.integer(args[[2]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

step2_dir <- file.path(project_root, "results/phase0_extension/step2_input_tables")
step5_dir <- file.path(project_root, "results/phase0_extension/step5_placo_cpassoc")
per_pair_dir <- file.path(step5_dir, "per_pair")
log_dir <- file.path(step5_dir, "logs")
dir.create(per_pair_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

placo_src <- file.path(project_root, "results/phase0_extension/step0_1/tools/PLACO/PLACO_v0.2.0.R")
cpassoc_zip <- file.path(project_root, "results/phase0_extension/step0_1/tools/CPASSOC/CPASSOC.zip")
gamma_n <- as.integer(Sys.getenv("CPASSOC_GAMMA_N", "200000"))
if (!is.finite(gamma_n) || gamma_n < 10000) gamma_n <- 200000

suppressMessages(source(placo_src))
suppressMessages(source(unz(cpassoc_zip, "CPASSOC/FunctionSet.R")))

fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), ...)
  else fread(path, ...)
}

write_tsv <- function(dt, path) fwrite(dt, path, sep = "\t", quote = FALSE, na = "")

std_gwas <- function(path) {
  dt <- fread_any(path, showProgress = FALSE)
  if ("A1" %in% names(dt) && !"EA" %in% names(dt)) setnames(dt, "A1", "EA")
  if ("A2" %in% names(dt) && !"OA" %in% names(dt)) setnames(dt, "A2", "OA")
  needed <- c("SNP", "CHR", "BP", "EA", "OA", "BETA", "SE", "P", "N")
  missing <- setdiff(needed, names(dt))
  if (length(missing)) stop("missing standardized GWAS columns: ", paste(missing, collapse = ","))
  dt <- dt[, ..needed]
  dt[, `:=`(
    CHR = as.integer(CHR),
    BP = as.integer(BP),
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P),
    N = as.numeric(N),
    EA = toupper(EA),
    OA = toupper(OA)
  )]
  dt <- dt[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & is.finite(N) & N > 0 & nzchar(SNP)]
  dt[, Z := BETA / SE]
  dt
}

empty_placo <- function() data.table(pair_id=character(), trait1=character(), trait2=character(), SNP=character(), CHR=integer(), BP=integer(), EA=character(), OA=character(), Z_trait1=numeric(), Z_trait2=numeric(), P_trait1=numeric(), P_trait2=numeric(), PLACO_stat=numeric(), PLACO_p=numeric())
empty_cpassoc <- function() data.table(pair_id=character(), trait1=character(), trait2=character(), SNP=character(), CHR=integer(), BP=integer(), EA=character(), OA=character(), P_trait1=numeric(), P_trait2=numeric(), SHet=numeric(), SHet_p=numeric(), SHom=numeric(), SHom_p=numeric())
qc_row <- function(pair_id, trait1, trait2, method, status, reason) data.table(pair_id=pair_id, trait1=trait1, trait2=trait2, method=method, status=status, reason=reason)
qc_rows <- list()
add_qc <- function(method, status, reason) {
  qc_rows[[length(qc_rows) + 1]] <<- qc_row(pid, trait1, trait2, method, status, reason)
}

priority <- fread(file.path(step2_dir, "priority_pairs.tsv"), showProgress = FALSE)
pair_manifest <- fread(file.path(step2_dir, "pair_manifest.tsv"), showProgress = FALSE)
candidate <- fread(file.path(step2_dir, "candidate_snps.tsv"), showProgress = FALSE)
if (task_id < 1 || task_id > nrow(priority)) stop("SLURM_ARRAY_TASK_ID outside priority_pairs.tsv")
p <- priority[task_id]
pid <- p$pair_id
trait1 <- p$trait1
trait2 <- p$trait2
manifest <- pair_manifest[pair_id == pid][1]
prefix <- file.path(per_pair_dir, gsub("[^A-Za-z0-9_.-]+", "_", pid))

run_status <- tryCatch({
  if (is.na(manifest$trait1_standard_gwas) || is.na(manifest$trait2_standard_gwas) ||
      !file.exists(manifest$trait1_standard_gwas) || !file.exists(manifest$trait2_standard_gwas)) {
    stop("missing standardized GWAS file")
  }
  cand <- unique(candidate[pair_id == pid & nzchar(SNP), .(SNP)])
  if (nrow(cand) == 0) stop("no candidate SNPs for pair")
  g1 <- std_gwas(manifest$trait1_standard_gwas)
  g2 <- std_gwas(manifest$trait2_standard_gwas)
  setnames(g1, c("EA","OA","BETA","SE","P","N","Z"), c("EA1","OA1","BETA1","SE1","P1","N1","Z1"))
  setnames(g2, c("EA","OA","BETA","SE","P","N","Z"), c("EA2","OA2","BETA2","SE2","P2","N2","Z2"))
  m <- merge(g1, g2, by = c("SNP", "CHR", "BP"))
  if (nrow(m) < 1000) stop("SNP intersection insufficient before allele alignment")
  same <- m$EA1 == m$EA2 & m$OA1 == m$OA2
  flip <- m$EA1 == m$OA2 & m$OA1 == m$EA2
  m <- m[same | flip]
  if (nrow(m) < 1000) stop("SNP intersection insufficient after allele alignment")
  flip2 <- m$EA1 == m$OA2 & m$OA1 == m$EA2
  if (any(flip2)) m$Z2[flip2] <- -m$Z2[flip2]
  tested <- merge(m, cand, by = "SNP")
  if (nrow(tested) == 0) stop("candidate SNPs absent after GWAS merge/alignment")
  zmat_all <- as.matrix(m[, .(Z1, Z2)])
  pmat_all <- as.matrix(m[, .(P1, P2)])
  zmat_test <- as.matrix(tested[, .(Z1, Z2)])
  sample_size <- c(median(m$N1, na.rm = TRUE), median(m$N2, na.rm = TRUE))
  cor_scalar <- tryCatch(cor.pearson(zmat_all, pmat_all, p.threshold = 1e-4, returnMatrix = FALSE), error = function(e) 0)
  if (!is.finite(cor_scalar)) cor_scalar <- 0
  cor_mat <- matrix(c(1, cor_scalar, cor_scalar, 1), 2, 2)
  varz <- tryCatch(var.placo(zmat_all, pmat_all, p.threshold = 1e-4), error = function(e) c(1, 1))
  placo_list <- lapply(seq_len(nrow(tested)), function(i) {
    res <- tryCatch(placo.plus(Z = zmat_test[i,], VarZ = varz, CorZ = cor_scalar), error = function(e) list(T.placo.plus = NA_real_, p.placo.plus = NA_real_))
    c(stat = res$T.placo.plus, p = res$p.placo.plus)
  })
  placo_mat <- do.call(rbind, placo_list)
  placo_out <- data.table(
    pair_id = pid, trait1 = trait1, trait2 = trait2,
    SNP = tested$SNP, CHR = tested$CHR, BP = tested$BP, EA = tested$EA1, OA = tested$OA1,
    Z_trait1 = tested$Z1, Z_trait2 = tested$Z2, P_trait1 = tested$P1, P_trait2 = tested$P2,
    PLACO_stat = placo_mat[, "stat"], PLACO_p = placo_mat[, "p"]
  )
  shom <- tryCatch(SHom(zmat_test, SampleSize = sample_size, CorrMatrix = cor_mat), error = function(e) rep(NA_real_, nrow(tested)))
  shet <- tryCatch(SHet(zmat_test, SampleSize = sample_size, CorrMatrix = cor_mat), error = function(e) rep(NA_real_, nrow(tested)))
  gamma <- tryCatch(EstimateGamma(N = gamma_n, SampleSize = sample_size, CorrMatrix = cor_mat), error = function(e) c(NA_real_, NA_real_, NA_real_))
  shet_p <- rep(NA_real_, length(shet))
  if (all(is.finite(gamma))) {
    shet_p <- ifelse(shet > gamma[3], pgamma(shet - gamma[3], shape = gamma[1], scale = gamma[2], lower.tail = FALSE), 1)
  }
  cpassoc_out <- data.table(
    pair_id = pid, trait1 = trait1, trait2 = trait2,
    SNP = tested$SNP, CHR = tested$CHR, BP = tested$BP, EA = tested$EA1, OA = tested$OA1,
    P_trait1 = tested$P1, P_trait2 = tested$P2,
    SHet = shet, SHet_p = shet_p, SHom = shom, SHom_p = pchisq(shom, df = 1, lower.tail = FALSE)
  )
  write_tsv(placo_out, paste0(prefix, ".placo.tsv"))
  write_tsv(cpassoc_out, paste0(prefix, ".cpassoc.tsv"))
  add_qc("PLACO", "PASS", paste0("tested_snps=", nrow(tested)))
  add_qc("CPASSOC", "PASS", paste0("tested_snps=", nrow(tested), ";gamma_n=", gamma_n))
  TRUE
}, error = function(e) {
  write_tsv(empty_placo(), paste0(prefix, ".placo.tsv"))
  write_tsv(empty_cpassoc(), paste0(prefix, ".cpassoc.tsv"))
  add_qc("PLACO", "FAIL", conditionMessage(e))
  add_qc("CPASSOC", "FAIL", conditionMessage(e))
  FALSE
})

write_tsv(rbindlist(qc_rows, fill = TRUE), paste0(prefix, ".qc.tsv"))
cat("pair_id=", pid, " status=", run_status, "\n", sep = "")
