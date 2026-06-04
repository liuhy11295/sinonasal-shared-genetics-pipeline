#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(MRPRESSO)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("PROJECT_ROOT", "/platform_data/p_user/p010/phase0")
task_id <- if (length(args) >= 2) as.integer(args[[2]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

step2_dir <- file.path(project_root, "results/phase0_extension/step2_input_tables")
out_dir <- file.path(project_root, "results/phase0_extension/step4_lcv_mr")
per_pair_dir <- file.path(out_dir, "per_pair")
tmp_dir <- file.path(out_dir, "tmp", paste0("mr_", task_id))
dir.create(per_pair_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

plink <- Sys.getenv("PLINK_BIN", file.path(Sys.getenv("CONDA_PREFIX", ""), "bin/plink"))
ld_ref_dir <- Sys.getenv("PHASE0_LD_PLINK_DIR", file.path(project_root, "data/reference/lava/1000G_Phase3_plinkfiles/1000G_EUR_Phase3_plink"))
ld_ref_prefix_template <- file.path(ld_ref_dir, "1000G.EUR.QC.%s")
p_threshold <- as.numeric(Sys.getenv("MR_P_THRESHOLD", "5e-8"))
min_iv <- as.integer(Sys.getenv("MR_MIN_IV", "4"))

fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), ...) else fread(path, ...)
}

write_tsv <- function(dt, path) fwrite(dt, path, sep = "\t", quote = FALSE, na = "")

std_gwas <- function(path) {
  dt <- fread_any(path, showProgress = FALSE)
  if ("EA" %in% names(dt) && !"A1" %in% names(dt)) setnames(dt, "EA", "A1")
  if ("OA" %in% names(dt) && !"A2" %in% names(dt)) setnames(dt, "OA", "A2")
  if ("EAF" %in% names(dt) && !"FRQ" %in% names(dt)) setnames(dt, "EAF", "FRQ")
  needed <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N")
  missing <- setdiff(needed, names(dt))
  if (length(missing)) stop("missing GWAS columns: ", paste(missing, collapse = ","))
  if (!"FRQ" %in% names(dt)) dt[, FRQ := NA_real_]
  dt <- dt[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, FRQ, N)]
  dt[, `:=`(CHR = as.integer(CHR), BP = as.integer(BP), BETA = as.numeric(BETA), SE = as.numeric(SE),
            P = as.numeric(P), FRQ = as.numeric(FRQ), N = as.numeric(N))]
  dt[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & is.finite(N) & N > 0 & nzchar(SNP)]
}

clump_instruments <- function(gwas, label) {
  inst <- gwas[P <= p_threshold]
  inst[, FSTAT := (BETA / SE)^2]
  inst <- inst[FSTAT >= 10]
  if (nrow(inst) == 0) return(inst[0])
  out <- list()
  for (chr in sort(unique(inst$CHR))) {
    bfile <- sprintf(ld_ref_prefix_template, chr)
    if (!file.exists(paste0(bfile, ".bed"))) next
    assoc <- inst[CHR == chr, .(SNP, P)]
    assoc_path <- file.path(tmp_dir, paste0(label, ".chr", chr, ".assoc"))
    fwrite(assoc, assoc_path, sep = "\t", quote = FALSE)
    prefix <- file.path(tmp_dir, paste0(label, ".chr", chr, ".clump"))
    cmd <- sprintf(
      "%s --bfile %s --clump %s --clump-snp-field SNP --clump-field P --clump-p1 %g --clump-p2 1 --clump-r2 0.001 --clump-kb 10000 --out %s --allow-no-sex",
      shQuote(plink), shQuote(bfile), shQuote(assoc_path), p_threshold, shQuote(prefix)
    )
    status <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    clumped <- paste0(prefix, ".clumped")
    if (status == 0 && file.exists(clumped) && file.info(clumped)$size > 0) {
      cl <- fread(clumped, fill = TRUE, showProgress = FALSE)
      if ("SNP" %in% names(cl)) out[[as.character(chr)]] <- inst[SNP %in% cl$SNP]
    }
  }
  if (!length(out)) return(inst[0])
  unique(rbindlist(out, fill = TRUE), by = "SNP")
}

is_pal <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")

harmonize <- function(exposure, outcome) {
  e <- copy(exposure)
  o <- outcome[SNP %in% e$SNP]
  setnames(e, c("A1", "A2", "BETA", "SE", "P", "FRQ"), c("A1.exposure", "A2.exposure", "beta.exposure", "se.exposure", "pval.exposure", "eaf.exposure"))
  setnames(o, c("A1", "A2", "BETA", "SE", "P", "FRQ"), c("A1.outcome", "A2.outcome", "beta.outcome", "se.outcome", "pval.outcome", "eaf.outcome"))
  m <- merge(e, o, by = c("SNP", "CHR", "BP"), allow.cartesian = FALSE)
  if (!nrow(m)) return(m)
  keep_same <- m$A1.exposure == m$A1.outcome & m$A2.exposure == m$A2.outcome
  keep_flip <- m$A1.exposure == m$A2.outcome & m$A2.exposure == m$A1.outcome
  m <- m[keep_same | keep_flip]
  if (!nrow(m)) return(m)
  flip <- m$A1.exposure == m$A2.outcome & m$A2.exposure == m$A1.outcome
  if (any(flip)) m$beta.outcome[flip] <- -m$beta.outcome[flip]
  m <- m[!is_pal(A1.exposure, A2.exposure)]
  m[, FSTAT := (beta.exposure / se.exposure)^2]
  m[FSTAT >= 10]
}

p_from_z <- function(z) 2 * pnorm(-abs(z))

ivw <- function(dat) {
  w <- 1 / (dat$se.outcome^2)
  beta <- sum(w * dat$beta.exposure * dat$beta.outcome) / sum(w * dat$beta.exposure^2)
  se <- sqrt(1 / sum(w * dat$beta.exposure^2))
  c(beta = beta, se = se, p = p_from_z(beta / se))
}

egger <- function(dat) {
  w <- 1 / (dat$se.outcome^2)
  fit <- summary(lm(beta.outcome ~ beta.exposure, weights = w, data = dat))
  beta <- coef(fit)["beta.exposure", "Estimate"]
  se <- coef(fit)["beta.exposure", "Std. Error"]
  p <- coef(fit)["beta.exposure", "Pr(>|t|)"]
  intercept <- coef(fit)["(Intercept)", "Estimate"]
  intercept_p <- coef(fit)["(Intercept)", "Pr(>|t|)"]
  c(beta = beta, se = se, p = p, intercept = intercept, intercept_p = intercept_p)
}

weighted_median <- function(dat) {
  ratio <- dat$beta.outcome / dat$beta.exposure
  ratio_se <- abs(dat$se.outcome / dat$beta.exposure)
  w <- 1 / (ratio_se^2)
  ord <- order(ratio)
  ratio <- ratio[ord]; w <- w[ord]
  cs <- cumsum(w) / sum(w)
  beta <- ratio[which(cs >= 0.5)[1]]
  se <- sqrt(1 / sum(w))
  c(beta = beta, se = se, p = p_from_z(beta / se))
}

mr_presso_run <- function(dat) {
  res <- list(beta = NA_real_, se = NA_real_, p = NA_real_, global_p = NA_real_, n_outlier = NA_integer_)
  if (nrow(dat) < 4) return(res)
  tryCatch({
    d <- as.data.frame(dat[, .(beta.outcome, beta.exposure, se.outcome, se.exposure)])
    pr <- mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = FALSE, data = d,
      NbDistribution = 500, seed = 1
    )
    main <- pr$`Main MR results`
    if (!is.null(main) && nrow(main) >= 1) {
      last <- main[nrow(main), ]
      res$beta <- suppressWarnings(as.numeric(last$`Causal Estimate`))
      res$se <- suppressWarnings(as.numeric(last$Sd))
      res$p <- suppressWarnings(as.numeric(last$`P-value`))
    }
    gp <- pr$`MR-PRESSO results`$`Global Test`$Pvalue
    res$global_p <- suppressWarnings(as.numeric(gp))
    outliers <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    res$n_outlier <- if (is.null(outliers)) 0L else length(outliers)
  }, error = function(e) {})
  res
}

run_direction <- function(pair_id, exposure_name, outcome_name, exposure_path, outcome_path) {
  exposure_gwas <- std_gwas(exposure_path)
  outcome_gwas <- std_gwas(outcome_path)
  inst <- clump_instruments(exposure_gwas, paste0(pair_id, "_", exposure_name))
  dat <- harmonize(inst, outcome_gwas)
  if (nrow(dat) < min_iv) {
    return(list(
      feasibility = data.table(direction = paste(exposure_name, outcome_name, sep = "_to_"), n_iv = nrow(dat), status = "insufficient_IV", reason = paste0("independent_IV_lt_", min_iv)),
      mr = data.table(), sens = data.table()
    ))
  }
  methods <- list()
  ivw_res <- ivw(dat)
  methods[[length(methods) + 1]] <- data.table(method = "IVW", beta = ivw_res["beta"], se = ivw_res["se"], p_value = ivw_res["p"])
  egger_res <- egger(dat)
  methods[[length(methods) + 1]] <- data.table(method = "MR-Egger", beta = egger_res["beta"], se = egger_res["se"], p_value = egger_res["p"])
  wm_res <- weighted_median(dat)
  methods[[length(methods) + 1]] <- data.table(method = "Weighted Median", beta = wm_res["beta"], se = wm_res["se"], p_value = wm_res["p"])
  presso <- mr_presso_run(dat)
  methods[[length(methods) + 1]] <- data.table(method = "MR-PRESSO", beta = presso$beta, se = presso$se, p_value = presso$p)
  mr <- rbindlist(methods, fill = TRUE)
  mr[, `:=`(
    pair_id = pair_id, exposure = exposure_name, outcome = outcome_name,
    OR = exp(beta), CI_lower = exp(beta - 1.96 * se), CI_upper = exp(beta + 1.96 * se)
  )]
  setcolorder(mr, c("pair_id", "exposure", "outcome", "method", "beta", "se", "p_value", "OR", "CI_lower", "CI_upper"))
  het_beta <- ivw_res["beta"]
  q <- sum((dat$beta.outcome - het_beta * dat$beta.exposure)^2 / dat$se.outcome^2)
  q_p <- pchisq(q, df = max(1, nrow(dat) - 1), lower.tail = FALSE)
  sens <- data.table(
    pair_id = pair_id, exposure = exposure_name, outcome = outcome_name,
    egger_intercept = egger_res["intercept"], egger_p = egger_res["intercept_p"],
    heterogeneity_q = q, heterogeneity_p = q_p,
    mr_presso_global_p = presso$global_p, n_outlier = presso$n_outlier
  )
  list(
    feasibility = data.table(direction = paste(exposure_name, outcome_name, sep = "_to_"), n_iv = nrow(dat), status = "eligible", reason = "OK"),
    mr = mr, sens = sens
  )
}

candidates <- fread(file.path(out_dir, "mr_candidate_pairs.tsv"), showProgress = FALSE)
manifest <- fread(file.path(step2_dir, "pair_manifest.tsv"), showProgress = FALSE)
if (task_id < 1 || task_id > nrow(candidates)) stop("SLURM_ARRAY_TASK_ID outside MR candidate range")
row <- candidates[task_id]
pm <- manifest[pair_id == row$pair_id][1]

res1 <- run_direction(row$pair_id, row$trait1, row$trait2, pm$trait1_standard_gwas, pm$trait2_standard_gwas)
res2 <- run_direction(row$pair_id, row$trait2, row$trait1, pm$trait2_standard_gwas, pm$trait1_standard_gwas)
feas_dirs <- rbindlist(list(res1$feasibility, res2$feasibility), fill = TRUE)
pair_status <- if (any(feas_dirs$status == "eligible")) "eligible" else "insufficient_IV"
pair_reason <- paste(paste(feas_dirs$direction, feas_dirs$n_iv, feas_dirs$status, sep = ":"), collapse = ";")
feas_pair <- data.table(pair_id = row$pair_id, trait1 = row$trait1, trait2 = row$trait2,
                        n_IV = max(feas_dirs$n_iv), status = pair_status, reason = pair_reason)
write_tsv(feas_pair, file.path(per_pair_dir, paste0(row$pair_id, ".mr_feasibility.tsv")))
write_tsv(rbindlist(list(res1$mr, res2$mr), fill = TRUE), file.path(per_pair_dir, paste0(row$pair_id, ".mr_results.tsv")))
write_tsv(rbindlist(list(res1$sens, res2$sens), fill = TRUE), file.path(per_pair_dir, paste0(row$pair_id, ".mr_sensitivity.tsv")))
