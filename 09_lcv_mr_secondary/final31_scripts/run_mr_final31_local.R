#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(MRPRESSO)
})

thread_limit <- as.integer(Sys.getenv("LCVMR_THREADS", "2"))
if (!is.finite(thread_limit) || thread_limit < 1L) thread_limit <- 2L
setDTthreads(thread_limit)

base <- Sys.getenv("BASE")
if (!nzchar(base)) stop("Set BASE for the final 31-pair LCV/MR analysis.")
out_root <- Sys.getenv("OUT", file.path(base, "results_final31"))
plink <- Sys.getenv("PLINK_BIN")
if (!nzchar(plink)) stop("Set PLINK_BIN.")
ld_ref_template <- file.path(base, "mr/reference_plink/1000G.EUR.QC.%s")
p_threshold <- as.numeric(Sys.getenv("MR_P_THRESHOLD", "5e-8"))
min_iv <- as.integer(Sys.getenv("MR_MIN_IV", "4"))
download_dir <- file.path(base, "mr/raw_gwas/downloads")
std_dir <- file.path(out_root, "MR/standardized")
per_pair_dir <- file.path(out_root, "MR/per_pair")
tmp_root <- file.path(out_root, "MR/tmp")
instrument_dir <- file.path(out_root, "MR/instruments")
lookup_dir <- file.path(out_root, "MR/outcome_lookup")
dir.create(std_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(per_pair_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
dir.create(instrument_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(lookup_dir, recursive = TRUE, showWarnings = FALSE)

pair_manifest <- file.path(out_root, "00_pair_manifest_final31.tsv")
trait_manifest <- file.path(base, "manifests/final_24_traits_input_manifest.tsv")
raw_manifest <- file.path(base, "mr/raw_gwas/MR_raw_GWAS_download_manifest.tsv")
log_path <- file.path(out_root, "MR/mr_run_log.txt")

write_tsv <- function(dt, path) fwrite(dt, path, sep = "\t", quote = FALSE, na = "")
fread_gz <- function(path, ...) fread(cmd = paste("gzip -dc", shQuote(path)), showProgress = FALSE, ...)

clean_url <- function(x) gsub('^"|"$', "", gsub("\r", "", x))

standardize_gwas_catalog <- function(infile, outfile, n) {
  header <- names(fread_gz(infile, nrows = 0))
  harmonized_cols <- c(
    "hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele",
    "hm_beta", "standard_error", "p_value", "hm_effect_allele_frequency"
  )
  standard_cols <- c(
    "rsid", "chromosome", "base_pair_location", "other_allele", "effect_allele",
    "beta", "standard_error", "p_value", "effect_allele_frequency"
  )
  if (all(harmonized_cols %in% header)) {
    cols <- harmonized_cols
  } else if (all(standard_cols %in% header)) {
    cols <- standard_cols
  } else {
    stop("unsupported GWAS Catalog column schema")
  }
  dt <- fread_gz(infile, select = cols)
  setnames(dt, cols, c("SNP", "CHR", "BP", "A2", "A1", "BETA", "SE", "P", "FRQ"))
  dt[, `:=`(
    SNP = sub(",.*$", "", as.character(SNP)),
    CHR = as.integer(CHR), BP = as.integer(BP),
    A1 = toupper(as.character(A1)), A2 = toupper(as.character(A2)),
    BETA = as.numeric(BETA), SE = as.numeric(SE), P = as.numeric(P),
    FRQ = as.numeric(FRQ), N = as.numeric(n)
  )]
  dt <- dt[grepl("^rs", SNP) & nchar(A1) == 1 & nchar(A2) == 1 &
             is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1]
  dt <- unique(dt, by = "SNP")
  write_tsv(dt[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, FRQ, N)], outfile)
  nrow(dt)
}

standardize_finngen <- function(infile, outfile, n) {
  cols <- c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt")
  dt <- fread_gz(infile, select = cols)
  setnames(dt, cols, c("CHR", "BP", "A2", "A1", "SNP", "P", "BETA", "SE", "FRQ"))
  dt[, `:=`(
    SNP = sub(",.*$", "", as.character(SNP)),
    CHR = as.integer(CHR), BP = as.integer(BP),
    A1 = toupper(as.character(A1)), A2 = toupper(as.character(A2)),
    BETA = as.numeric(BETA), SE = as.numeric(SE), P = as.numeric(P),
    FRQ = as.numeric(FRQ), N = as.numeric(n)
  )]
  dt <- dt[grepl("^rs", SNP) & nchar(A1) == 1 & nchar(A2) == 1 &
             is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1]
  dt <- unique(dt, by = "SNP")
  write_tsv(dt[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, FRQ, N)], outfile)
  nrow(dt)
}

std_gwas <- function(path) {
  dt <- fread(path, showProgress = FALSE)
  dt[, `:=`(CHR = as.integer(CHR), BP = as.integer(BP), BETA = as.numeric(BETA), SE = as.numeric(SE),
            P = as.numeric(P), FRQ = as.numeric(FRQ), N = as.numeric(N))]
  dt[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & is.finite(N) & N > 0 & nzchar(SNP)]
}

clump_instruments <- function(gwas, label, tmp_dir) {
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
  c(beta = coef(fit)["beta.exposure", "Estimate"],
    se = coef(fit)["beta.exposure", "Std. Error"],
    p = coef(fit)["beta.exposure", "Pr(>|t|)"],
    intercept = coef(fit)["(Intercept)", "Estimate"],
    intercept_p = coef(fit)["(Intercept)", "Pr(>|t|)"])
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
  parse_presso_p <- function(x) {
    if (is.null(x) || !length(x)) return(NA_real_)
    suppressWarnings(as.numeric(sub("^<", "", as.character(x)[1])))
  }
  res <- list(beta = NA_real_, se = NA_real_, p = NA_real_, global_p = NA_real_,
              n_outlier = NA_integer_, analysis = NA_character_, notes = "")
  if (nrow(dat) < 4) return(res)
  tryCatch({
    d <- as.data.frame(dat[, .(beta.outcome, beta.exposure, se.outcome, se.exposure)])
    pr <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                    SdOutcome = "se.outcome", SdExposure = "se.exposure",
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = d,
                    NbDistribution = 2000, seed = 1)
    main <- pr$`Main MR results`
    if (!is.null(main) && nrow(main) >= 1) {
      valid <- which(is.finite(suppressWarnings(as.numeric(main$`Causal Estimate`))))
      chosen <- if (length(valid)) tail(valid, 1) else integer()
      if (length(chosen)) {
        row <- main[chosen, ]
        res$beta <- suppressWarnings(as.numeric(row$`Causal Estimate`))
        res$se <- suppressWarnings(as.numeric(row$Sd))
        res$p <- parse_presso_p(row$`P-value`)
        res$analysis <- as.character(row$`MR Analysis`)
      }
    }
    gp <- pr$`MR-PRESSO results`$`Global Test`$Pvalue
    res$global_p <- parse_presso_p(gp)
    outliers <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    res$n_outlier <- if (is.null(outliers) ||
      identical(outliers, "No significant outliers")) 0L else length(outliers)
  }, error = function(e) {
    res$notes <- conditionMessage(e)
  })
  res
}

run_direction <- function(pair_id, exposure_name, outcome_name, exposure_path, outcome_path,
                          instrument_prep) {
  exposure <- std_gwas(exposure_path)
  outcome <- std_gwas(outcome_path)
  harm <- harmonize(exposure, outcome)
  dat <- harm$dat
  direction <- paste(exposure_name, outcome_name, sep = "_to_")
  prep <- instrument_prep[trait_id == exposure_name]
  inst <- data.table(pair_id = pair_id, exposure = exposure_name, outcome = outcome_name, direction = direction,
                     n_instruments_raw = prep$n_instruments_raw,
                     n_instruments_after_strength_filter = prep$n_instruments_after_strength_filter,
                     n_instruments_after_clump = prep$n_instruments_after_clump,
                     n_instruments_after_harmonise = nrow(dat),
                     n_palindromic_removed = harm$n_pal,
                     status = if (nrow(dat) >= min_iv) "eligible" else "insufficient_instruments")
  if (nrow(dat) < min_iv) {
    return(list(inst = inst, mr = data.table(pair_id = pair_id, exposure = exposure_name, outcome = outcome_name,
                                             direction = direction, method = NA_character_, nsnp = nrow(dat),
                                             beta = NA_real_, se = NA_real_, pval = NA_real_, or = NA_real_,
                                             or_lci95 = NA_real_, or_uci95 = NA_real_, status = "insufficient_instruments",
                                             notes = paste0("harmonised_IV_lt_", min_iv)),
                sens = data.table(pair_id = pair_id, exposure = exposure_name, outcome = outcome_name,
                                  direction = direction, egger_intercept = NA_real_, egger_p = NA_real_,
                                  heterogeneity_q = NA_real_, heterogeneity_p = NA_real_,
                                  mr_presso_global_p = NA_real_, n_outlier = NA_integer_, status = "insufficient_instruments")))
  }
  ivw_res <- ivw(dat); egger_res <- egger(dat); wm_res <- weighted_median(dat); presso <- mr_presso_run(dat)
  mr <- rbindlist(list(
    data.table(method = "IVW", beta = ivw_res["beta"], se = ivw_res["se"], pval = ivw_res["p"]),
    data.table(method = "MR-Egger", beta = egger_res["beta"], se = egger_res["se"], pval = egger_res["p"]),
    data.table(method = "Weighted Median", beta = wm_res["beta"], se = wm_res["se"], pval = wm_res["p"]),
    data.table(method = "MR-PRESSO", beta = presso$beta, se = presso$se, pval = presso$p)
  ), fill = TRUE)
  mr[, `:=`(pair_id = pair_id, exposure = exposure_name, outcome = outcome_name, direction = direction,
            nsnp = nrow(dat), or = exp(beta), or_lci95 = exp(beta - 1.96 * se),
            or_uci95 = exp(beta + 1.96 * se), status = "eligible", notes = "")]
  setcolorder(mr, c("pair_id", "exposure", "outcome", "direction", "method", "nsnp", "beta", "se", "pval", "or", "or_lci95", "or_uci95", "status", "notes"))
  q <- sum((dat$beta.outcome - ivw_res["beta"] * dat$beta.exposure)^2 / dat$se.outcome^2)
  q_p <- pchisq(q, df = max(1, nrow(dat) - 1), lower.tail = FALSE)
  sens <- data.table(pair_id = pair_id, exposure = exposure_name, outcome = outcome_name, direction = direction,
                     egger_intercept = egger_res["intercept"], egger_p = egger_res["intercept_p"],
                     heterogeneity_q = q, heterogeneity_p = q_p,
                     mr_presso_global_p = presso$global_p, n_outlier = presso$n_outlier, status = "eligible")
  list(inst = inst, mr = mr, sens = sens)
}

sink(log_path, split = TRUE)
cat("MR final31 local run\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("base:", base, "\n")
cat("plink:", plink, "\n")
cat("R_version:", R.version.string, "\n")
cat("data.table_version:", as.character(packageVersion("data.table")), "\n")
cat("MRPRESSO_version:", as.character(packageVersion("MRPRESSO")), "\n")
cat("thread_limit:", thread_limit, "\n")
cat("data.table_threads:", getDTthreads(), "\n")
cat("p_threshold:", p_threshold, " min_iv:", min_iv, "\n")
cat("clump_r2: 0.001 clump_kb: 10000 F_stat_min: 10 palindromic_policy: remove\n")

traits <- fread(trait_manifest, showProgress = FALSE)
raw <- fread(raw_manifest, showProgress = FALSE)
setnames(raw, names(raw), gsub('^"|"$', "", names(raw)))
for (col in names(raw)) raw[[col]] <- clean_url(raw[[col]])
for (col in names(traits)) traits[[col]] <- clean_url(traits[[col]])
trait_info <- merge(raw, traits[, .(trait_id, sample_size)], by = "trait_id", all.x = TRUE)

std_rows <- list()
for (i in seq_len(nrow(trait_info))) {
  tr <- trait_info[i]
  url <- tr$mr_download_url
  infile <- file.path(download_dir, basename(url))
  outfile <- file.path(std_dir, paste0(tr$trait_id, ".standardized.tsv.gz"))
  row <- data.table(trait_id = tr$trait_id, source = tr$source, raw_file = infile,
                    standardized_file = outfile, n_rows = NA_integer_, status = "failed", notes = "")
  tryCatch({
    if (!file.exists(infile) || file.info(infile)$size == 0) stop("missing raw GWAS")
    if (!file.exists(outfile) || file.info(outfile)$size == 0) {
      n <- as.numeric(tr$sample_size)
      nr <- if (grepl("FinnGen", tr$source, ignore.case = TRUE)) {
        standardize_finngen(infile, outfile, n)
      } else {
        standardize_gwas_catalog(infile, outfile, n)
      }
      row[, notes := "created"]
    } else {
      nr <- NA_integer_
      row[, notes := "reused_existing"]
    }
    row[, `:=`(n_rows = nr, status = "ok")]
  }, error = function(e) {
    row[, notes := conditionMessage(e)]
  })
  std_rows[[i]] <- row
  cat("standardize ", i, "/", nrow(trait_info), " ", tr$trait_id, " ", row$status, " ", row$notes, "\n", sep = "")
}
std_inv <- rbindlist(std_rows, fill = TRUE)
write_tsv(std_inv, file.path(out_root, "MR/mr_standardized_input_inventory.tsv"))

instrument_rows <- list()
all_instrument_snps <- character()
for (i in seq_len(nrow(trait_info))) {
  trait_id <- trait_info$trait_id[i]
  std_path <- file.path(std_dir, paste0(trait_id, ".standardized.tsv.gz"))
  instrument_path <- file.path(instrument_dir, paste0(trait_id, ".clumped.tsv"))
  row <- data.table(
    trait_id = trait_id,
    n_instruments_raw = NA_integer_,
    n_instruments_after_strength_filter = NA_integer_,
    n_instruments_after_clump = NA_integer_,
    status = "failed",
    notes = ""
  )
  tryCatch({
    if (!file.exists(std_path) || file.info(std_path)$size == 0) stop("missing standardized GWAS")
    gwas <- std_gwas(std_path)
    trait_tmp <- file.path(tmp_root, "trait_clumping", trait_id)
    dir.create(trait_tmp, recursive = TRUE, showWarnings = FALSE)
    cl <- clump_instruments(gwas, trait_id, trait_tmp)
    write_tsv(cl$clumped, instrument_path)
    row[, `:=`(
      n_instruments_raw = cl$raw_n,
      n_instruments_after_strength_filter = cl$strong_n,
      n_instruments_after_clump = nrow(cl$clumped),
      status = "ok"
    )]
    all_instrument_snps <- union(all_instrument_snps, cl$clumped$SNP)
    rm(gwas, cl)
    gc(verbose = FALSE)
  }, error = function(e) {
    row[, notes := conditionMessage(e)]
  })
  instrument_rows[[i]] <- row
  cat("instrument_prep ", i, "/", nrow(trait_info), " ", trait_id, " ",
      row$status, " ", row$notes, "\n", sep = "")
}
instrument_prep <- rbindlist(instrument_rows, fill = TRUE)
write_tsv(instrument_prep, file.path(out_root, "MR/mr_trait_instrument_preparation.tsv"))

if (!length(all_instrument_snps)) stop("no clumped instruments were available")
for (i in seq_len(nrow(trait_info))) {
  trait_id <- trait_info$trait_id[i]
  std_path <- file.path(std_dir, paste0(trait_id, ".standardized.tsv.gz"))
  lookup_path <- file.path(lookup_dir, paste0(trait_id, ".outcome_lookup.tsv"))
  if (!file.exists(std_path) || file.info(std_path)$size == 0) next
  gwas <- std_gwas(std_path)
  lookup <- gwas[SNP %in% all_instrument_snps]
  write_tsv(lookup, lookup_path)
  cat("outcome_lookup ", i, "/", nrow(trait_info), " ", trait_id,
      " rows=", nrow(lookup), "\n", sep = "")
  rm(gwas, lookup)
  gc(verbose = FALSE)
}

pairs <- fread(pair_manifest, showProgress = FALSE)
all_mr <- list(); all_sens <- list(); all_inst <- list(); qc <- list()
for (i in seq_len(nrow(pairs))) {
  p <- pairs[i]
  t1_instruments <- file.path(instrument_dir, paste0(p$trait1, ".clumped.tsv"))
  t2_instruments <- file.path(instrument_dir, paste0(p$trait2, ".clumped.tsv"))
  t1_lookup <- file.path(lookup_dir, paste0(p$trait1, ".outcome_lookup.tsv"))
  t2_lookup <- file.path(lookup_dir, paste0(p$trait2, ".outcome_lookup.tsv"))
  cat("pair ", i, "/", nrow(pairs), " ", p$pair_id, "\n", sep = "")
  pair_status <- "success"; notes <- ""
  tryCatch({
    required <- c(t1_instruments, t2_instruments, t1_lookup, t2_lookup)
    if (!all(file.exists(required))) stop("missing instrument or outcome lookup")
    r1 <- run_direction(p$pair_id, p$trait1, p$trait2, t1_instruments, t2_lookup,
                        instrument_prep)
    r2 <- run_direction(p$pair_id, p$trait2, p$trait1, t2_instruments, t1_lookup,
                        instrument_prep)
    pair_mr <- rbindlist(list(r1$mr, r2$mr), fill = TRUE)
    pair_sens <- rbindlist(list(r1$sens, r2$sens), fill = TRUE)
    pair_inst <- rbindlist(list(r1$inst, r2$inst), fill = TRUE)
    write_tsv(pair_mr, file.path(per_pair_dir, paste0(p$pair_id, ".mr_results.tsv")))
    write_tsv(pair_sens, file.path(per_pair_dir, paste0(p$pair_id, ".mr_sensitivity.tsv")))
    write_tsv(pair_inst, file.path(per_pair_dir, paste0(p$pair_id, ".mr_instrument_summary.tsv")))
    all_mr[[length(all_mr) + 1]] <- pair_mr
    all_sens[[length(all_sens) + 1]] <- pair_sens
    all_inst[[length(all_inst) + 1]] <- pair_inst
    if (!any(pair_inst$status == "eligible")) pair_status <- "insufficient_instruments"
  }, error = function(e) {
    pair_status <<- "failed"
    notes <<- conditionMessage(e)
  })
  qc[[i]] <- data.table(pair_id = p$pair_id, trait1 = p$trait1, trait2 = p$trait2, status = pair_status, notes = notes)
  cat("pair_status ", p$pair_id, " ", pair_status, " ", notes, "\n", sep = "")
}

mr_all <- if (length(all_mr)) rbindlist(all_mr, fill = TRUE) else data.table()
sens_all <- if (length(all_sens)) rbindlist(all_sens, fill = TRUE) else data.table()
inst_all <- if (length(all_inst)) rbindlist(all_inst, fill = TRUE) else data.table()
qc_all <- rbindlist(qc, fill = TRUE)
write_tsv(mr_all, file.path(out_root, "MR/mr_results_all.tsv"))
write_tsv(sens_all, file.path(out_root, "MR/mr_sensitivity_all.tsv"))
write_tsv(inst_all, file.path(out_root, "MR/mr_instrument_summary.tsv"))
write_tsv(qc_all, file.path(out_root, "MR/mr_qc_summary.tsv"))
cat("Finished:", as.character(Sys.time()), "\n")
cat("pairs:", nrow(qc_all), " success:", sum(qc_all$status == "success"), " insufficient:", sum(qc_all$status == "insufficient_instruments"), " failed:", sum(qc_all$status == "failed"), "\n")
sink()
