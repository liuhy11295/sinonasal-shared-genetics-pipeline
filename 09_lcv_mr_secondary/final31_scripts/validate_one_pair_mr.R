#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(MRPRESSO)
})

base <- Sys.getenv("BASE")
if (!nzchar(base)) stop("Set BASE for the final 31-pair LCV/MR analysis.")
out_root <- Sys.getenv("OUT", file.path(base, "results_final31"))
pair_id <- Sys.getenv("PAIR_ID", "ALLERGIC_RHINITIS_GCST90038664__ALLERG_ASTHMA")
trait1 <- Sys.getenv("TRAIT1", "ALLERGIC_RHINITIS_GCST90038664")
trait2 <- Sys.getenv("TRAIT2", "ALLERG_ASTHMA")
plink <- Sys.getenv("PLINK_BIN")
if (!nzchar(plink)) stop("Set PLINK_BIN.")
ld_ref_template <- file.path(base, "mr/reference_plink/1000G.EUR.QC.%s")
p_threshold <- as.numeric(Sys.getenv("MR_P_THRESHOLD", "5e-8"))
min_iv <- as.integer(Sys.getenv("MR_MIN_IV", "4"))

val_dir <- file.path(out_root, "MR", "validation_one_pair", pair_id)
std_dir <- file.path(out_root, "MR", "validation_one_pair", "standardized")
tmp_dir <- file.path(val_dir, "tmp")
dir.create(val_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(std_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

download_dir <- file.path(base, "mr/raw_gwas/downloads")
ar_raw <- file.path(download_dir, "33959723-GCST90038664-EFO_0005854.h.tsv.gz")
asthma_raw <- file.path(download_dir, "finngen_R12_ALLERG_ASTHMA.gz")
ar_std <- file.path(std_dir, paste0(trait1, ".standardized.tsv.gz"))
asthma_std <- file.path(std_dir, paste0(trait2, ".standardized.tsv.gz"))

write_tsv <- function(dt, path) fwrite(dt, path, sep = "\t", quote = FALSE, na = "")
fread_gz <- function(path, ...) fread(cmd = paste("gzip -dc", shQuote(path)), showProgress = FALSE, ...)

standardize_gwas_catalog <- function(infile, outfile, n) {
  cols <- c("hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele",
            "hm_beta", "standard_error", "p_value", "hm_effect_allele_frequency")
  dt <- fread_gz(infile, select = cols)
  setnames(dt, cols, c("SNP", "CHR", "BP", "A2", "A1", "BETA", "SE", "P", "FRQ"))
  dt[, `:=`(
    SNP = sub(",.*$", "", as.character(SNP)),
    CHR = as.integer(CHR),
    BP = as.integer(BP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P),
    FRQ = as.numeric(FRQ),
    N = as.numeric(n)
  )]
  dt <- dt[grepl("^rs", SNP) & nchar(A1) == 1 & nchar(A2) == 1 &
             is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1]
  dt <- unique(dt, by = "SNP")
  write_tsv(dt[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, FRQ, N)], outfile)
  invisible(nrow(dt))
}

standardize_finngen <- function(infile, outfile, n) {
  cols <- c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt")
  dt <- fread_gz(infile, select = cols)
  setnames(dt, cols, c("CHR", "BP", "A2", "A1", "SNP", "P", "BETA", "SE", "FRQ"))
  dt[, `:=`(
    SNP = sub(",.*$", "", as.character(SNP)),
    CHR = as.integer(CHR),
    BP = as.integer(BP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P),
    FRQ = as.numeric(FRQ),
    N = as.numeric(n)
  )]
  dt <- dt[grepl("^rs", SNP) & nchar(A1) == 1 & nchar(A2) == 1 &
             is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1]
  dt <- unique(dt, by = "SNP")
  write_tsv(dt[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, FRQ, N)], outfile)
  invisible(nrow(dt))
}

std_gwas <- function(path) {
  dt <- fread(path, showProgress = FALSE)
  dt[, `:=`(CHR = as.integer(CHR), BP = as.integer(BP), BETA = as.numeric(BETA), SE = as.numeric(SE),
            P = as.numeric(P), FRQ = as.numeric(FRQ), N = as.numeric(N))]
  dt[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & is.finite(N) & N > 0 & nzchar(SNP)]
}

clump_instruments <- function(gwas, label) {
  inst_raw <- gwas[P <= p_threshold]
  inst_raw[, FSTAT := (BETA / SE)^2]
  inst <- inst_raw[FSTAT >= 10]
  if (!nrow(inst)) return(list(raw_n = nrow(inst_raw), strong_n = 0L, clumped = inst[0]))
  out <- list()
  for (chr in sort(unique(inst$CHR))) {
    bfile <- sprintf(ld_ref_template, chr)
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
  clumped_dt <- if (length(out)) unique(rbindlist(out, fill = TRUE), by = "SNP") else inst[0]
  list(raw_n = nrow(inst_raw), strong_n = nrow(inst), clumped = clumped_dt)
}

is_pal <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")

harmonize <- function(exposure, outcome) {
  e <- copy(exposure)
  o <- outcome[SNP %in% e$SNP]
  setnames(e, c("A1", "A2", "BETA", "SE", "P", "FRQ"), c("A1.exposure", "A2.exposure", "beta.exposure", "se.exposure", "pval.exposure", "eaf.exposure"))
  setnames(o, c("A1", "A2", "BETA", "SE", "P", "FRQ"), c("A1.outcome", "A2.outcome", "beta.outcome", "se.outcome", "pval.outcome", "eaf.outcome"))
  m0 <- merge(e, o, by = c("SNP", "CHR", "BP"), allow.cartesian = FALSE)
  if (!nrow(m0)) return(list(dat = m0, n_pal = 0L, n_allele_match = 0L))
  keep_same <- m0$A1.exposure == m0$A1.outcome & m0$A2.exposure == m0$A2.outcome
  keep_flip <- m0$A1.exposure == m0$A2.outcome & m0$A2.exposure == m0$A1.outcome
  m <- m0[keep_same | keep_flip]
  if (!nrow(m)) return(list(dat = m, n_pal = 0L, n_allele_match = 0L))
  flip <- m$A1.exposure == m$A2.outcome & m$A2.exposure == m$A1.outcome
  if (any(flip)) m$beta.outcome[flip] <- -m$beta.outcome[flip]
  n_pal <- sum(is_pal(m$A1.exposure, m$A2.exposure))
  m <- m[!is_pal(A1.exposure, A2.exposure)]
  m[, FSTAT := (beta.exposure / se.exposure)^2]
  list(dat = m[FSTAT >= 10], n_pal = n_pal, n_allele_match = nrow(m))
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
  }, error = function(e) {
    res$error <- conditionMessage(e)
  })
  res
}

run_direction <- function(exposure_name, outcome_name, exposure_path, outcome_path) {
  exposure_gwas <- std_gwas(exposure_path)
  outcome_gwas <- std_gwas(outcome_path)
  cl <- clump_instruments(exposure_gwas, paste0(pair_id, "_", exposure_name))
  harm <- harmonize(cl$clumped, outcome_gwas)
  dat <- harm$dat
  direction <- paste(exposure_name, outcome_name, sep = "_to_")
  inst_summary <- data.table(
    pair_id = pair_id, exposure = exposure_name, outcome = outcome_name, direction = direction,
    n_instruments_raw = cl$raw_n,
    n_instruments_after_strength_filter = cl$strong_n,
    n_instruments_after_clump = nrow(cl$clumped),
    n_instruments_after_harmonise = nrow(dat),
    n_palindromic_removed = harm$n_pal,
    status = if (nrow(dat) >= min_iv) "eligible" else "insufficient_instruments"
  )
  if (nrow(dat) < min_iv) return(list(inst = inst_summary, mr = data.table(), sens = data.table()))
  methods <- list()
  ivw_res <- ivw(dat)
  methods[[1]] <- data.table(method = "IVW", beta = ivw_res["beta"], se = ivw_res["se"], p_value = ivw_res["p"])
  egger_res <- egger(dat)
  methods[[2]] <- data.table(method = "MR-Egger", beta = egger_res["beta"], se = egger_res["se"], p_value = egger_res["p"])
  wm_res <- weighted_median(dat)
  methods[[3]] <- data.table(method = "Weighted Median", beta = wm_res["beta"], se = wm_res["se"], p_value = wm_res["p"])
  presso <- mr_presso_run(dat)
  methods[[4]] <- data.table(method = "MR-PRESSO", beta = presso$beta, se = presso$se, p_value = presso$p)
  mr <- rbindlist(methods, fill = TRUE)
  mr[, `:=`(
    pair_id = pair_id, exposure = exposure_name, outcome = outcome_name, direction = direction,
    nsnp = nrow(dat),
    OR = exp(beta), or_lci95 = exp(beta - 1.96 * se), or_uci95 = exp(beta + 1.96 * se),
    status = "eligible", notes = ""
  )]
  setcolorder(mr, c("pair_id", "exposure", "outcome", "direction", "method", "nsnp", "beta", "se", "p_value", "OR", "or_lci95", "or_uci95", "status", "notes"))
  q <- sum((dat$beta.outcome - ivw_res["beta"] * dat$beta.exposure)^2 / dat$se.outcome^2)
  q_p <- pchisq(q, df = max(1, nrow(dat) - 1), lower.tail = FALSE)
  sens <- data.table(
    pair_id = pair_id, exposure = exposure_name, outcome = outcome_name, direction = direction,
    egger_intercept = egger_res["intercept"], egger_p = egger_res["intercept_p"],
    heterogeneity_q = q, heterogeneity_p = q_p,
    mr_presso_global_p = presso$global_p, n_outlier = presso$n_outlier
  )
  list(inst = inst_summary, mr = mr, sens = sens)
}

log_lines <- c(
  paste("Started", Sys.time()),
  paste("pair_id", pair_id),
  paste("plink", plink),
  paste("p_threshold", p_threshold),
  paste("min_iv", min_iv)
)
writeLines(log_lines, file.path(val_dir, "validation_run_log.txt"))

if (!file.exists(ar_std)) standardize_gwas_catalog(ar_raw, ar_std, n = 484598)
if (!file.exists(asthma_std)) standardize_finngen(asthma_raw, asthma_std, n = 283740)

std_inventory <- data.table(
  trait = c(trait1, trait2),
  standardized_file = c(ar_std, asthma_std),
  n_rows = c(nrow(fread(ar_std, select = "SNP", showProgress = FALSE)),
             nrow(fread(asthma_std, select = "SNP", showProgress = FALSE)))
)
write_tsv(std_inventory, file.path(val_dir, "standardized_input_inventory.tsv"))

res1 <- run_direction(trait1, trait2, ar_std, asthma_std)
res2 <- run_direction(trait2, trait1, asthma_std, ar_std)
mr_new <- rbindlist(list(res1$mr, res2$mr), fill = TRUE)
sens_new <- rbindlist(list(res1$sens, res2$sens), fill = TRUE)
inst_new <- rbindlist(list(res1$inst, res2$inst), fill = TRUE)
write_tsv(mr_new, file.path(val_dir, "mr_results_recomputed.tsv"))
write_tsv(sens_new, file.path(val_dir, "mr_sensitivity_recomputed.tsv"))
write_tsv(inst_new, file.path(val_dir, "mr_instrument_summary_recomputed.tsv"))

old_mr_path <- file.path(out_root, "backup_20260607_235212", "mr_results.tsv")
if (!file.exists(old_mr_path)) old_mr_path <- file.path(out_root, "mr_results.tsv")
old <- fread(old_mr_path, showProgress = FALSE)
old <- old[pair_id == !!pair_id]
if ("p_value" %in% names(old) && !"pval" %in% names(old)) setnames(old, "p_value", "p_old")
cmp <- merge(
  old[, .(pair_id, exposure, outcome, method, beta_old = as.numeric(beta), se_old = as.numeric(se), p_old = as.numeric(p_old))],
  mr_new[, .(pair_id, exposure, outcome, method, beta_new = beta, se_new = se, p_new = p_value)],
  by = c("pair_id", "exposure", "outcome", "method"),
  all = TRUE
)
cmp[, `:=`(
  beta_abs_diff = abs(beta_new - beta_old),
  se_abs_diff = abs(se_new - se_old),
  p_log10_abs_diff = abs(-log10(p_new) - -log10(p_old)),
  same_beta_1e_6 = abs(beta_new - beta_old) < 1e-6,
  same_se_1e_6 = abs(se_new - se_old) < 1e-6
)]
write_tsv(cmp, file.path(val_dir, "mr_old_vs_recomputed_comparison.tsv"))

cat("Validation complete\n")
cat("Output:", val_dir, "\n")
print(inst_new)
print(cmp)
