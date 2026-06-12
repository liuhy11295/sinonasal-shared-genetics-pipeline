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
mr_dir <- file.path(out_root, "MR")
instrument_dir <- file.path(mr_dir, "instruments")
lookup_dir <- file.path(mr_dir, "outcome_lookup")
checkpoint_dir <- file.path(mr_dir, "presso_checkpoints")
per_pair_dir <- file.path(mr_dir, "per_pair")
log_path <- file.path(mr_dir, "mr_presso_repair_log.txt")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(dt, path) fwrite(dt, path, sep = "\t", quote = FALSE, na = "")

parse_presso_p <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  suppressWarnings(as.numeric(sub("^<", "", as.character(x)[1])))
}

is_pal <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")

harmonize <- function(exposure, outcome) {
  e <- copy(exposure)
  o <- outcome[SNP %in% e$SNP]
  setnames(
    e,
    c("A1", "A2", "BETA", "SE", "P", "FRQ"),
    c("A1.exposure", "A2.exposure", "beta.exposure", "se.exposure",
      "pval.exposure", "eaf.exposure")
  )
  setnames(
    o,
    c("A1", "A2", "BETA", "SE", "P", "FRQ"),
    c("A1.outcome", "A2.outcome", "beta.outcome", "se.outcome",
      "pval.outcome", "eaf.outcome")
  )
  m <- merge(e, o, by = c("SNP", "CHR", "BP"), allow.cartesian = FALSE)
  same <- m$A1.exposure == m$A1.outcome & m$A2.exposure == m$A2.outcome
  flip <- m$A1.exposure == m$A2.outcome & m$A2.exposure == m$A1.outcome
  m <- m[same | flip]
  if (!nrow(m)) return(m)
  reverse <- m$A1.exposure == m$A2.outcome & m$A2.exposure == m$A1.outcome
  m$beta.outcome[reverse] <- -m$beta.outcome[reverse]
  m <- m[!is_pal(A1.exposure, A2.exposure)]
  m[, FSTAT := (beta.exposure / se.exposure)^2]
  m[FSTAT >= 10]
}

run_presso <- function(dat) {
  warnings_seen <- character()
  result <- list(
    beta = NA_real_, se = NA_real_, pval = NA_real_,
    global_p = NA_real_, global_p_raw = NA_character_,
    n_outlier = NA_integer_, analysis = NA_character_,
    status = "failed", notes = ""
  )
  tryCatch({
    d <- as.data.frame(dat[, .(beta.outcome, beta.exposure, se.outcome, se.exposure)])
    pr <- withCallingHandlers(
      mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = d,
        NbDistribution = 2000,
        seed = 1
      ),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    main <- pr$`Main MR results`
    valid <- which(is.finite(suppressWarnings(as.numeric(main$`Causal Estimate`))))
    if (!length(valid)) stop("MR-PRESSO returned no finite causal estimate")
    chosen <- tail(valid, 1)
    selected <- main[chosen, ]
    result$beta <- as.numeric(selected$`Causal Estimate`)
    result$se <- as.numeric(selected$Sd)
    result$pval <- parse_presso_p(selected$`P-value`)
    result$analysis <- as.character(selected$`MR Analysis`)
    global_raw <- pr$`MR-PRESSO results`$`Global Test`$Pvalue
    result$global_p_raw <- as.character(global_raw)[1]
    result$global_p <- parse_presso_p(global_raw)
    outliers <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    if (is.null(outliers) || identical(outliers, "No significant outliers")) {
      result$n_outlier <- 0L
    } else {
      result$n_outlier <- length(outliers)
    }
    result$status <- "success"
    result$notes <- paste(unique(warnings_seen), collapse = " | ")
  }, error = function(e) {
    result$notes <- paste(
      c(conditionMessage(e), unique(warnings_seen)),
      collapse = " | "
    )
  })
  result
}

sink(log_path, split = TRUE)
cat("MR-PRESSO repair for final31\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("thread_limit:", thread_limit, "\n")
cat("NbDistribution: 2000\n")

inst <- fread(file.path(mr_dir, "mr_instrument_summary.tsv"), showProgress = FALSE)
eligible <- inst[status == "eligible"]
for (i in seq_len(nrow(eligible))) {
  row <- eligible[i]
  checkpoint <- file.path(checkpoint_dir, sprintf("%02d.tsv", i))
  if (file.exists(checkpoint) && file.info(checkpoint)$size > 0) {
    existing <- fread(checkpoint, showProgress = FALSE)
    if (nrow(existing) == 1 && existing$status == "success" &&
        existing$direction == row$direction) {
      cat("presso ", i, "/", nrow(eligible), " ", row$direction,
          " reused_checkpoint\n", sep = "")
      next
    }
  }

  exposure_path <- file.path(instrument_dir, paste0(row$exposure, ".clumped.tsv"))
  outcome_path <- file.path(lookup_dir, paste0(row$outcome, ".outcome_lookup.tsv"))
  result <- list(
    beta = NA_real_, se = NA_real_, pval = NA_real_,
    global_p = NA_real_, global_p_raw = NA_character_,
    n_outlier = NA_integer_, analysis = NA_character_,
    status = "failed", notes = ""
  )
  nsnp <- 0L
  tryCatch({
    exposure <- fread(exposure_path, showProgress = FALSE)
    outcome <- fread(outcome_path, showProgress = FALSE)
    dat <- harmonize(exposure, outcome)
    nsnp <- nrow(dat)
    if (nsnp < 4) stop("fewer than four harmonized instruments")
    result <- run_presso(dat)
  }, error = function(e) {
    result$notes <- conditionMessage(e)
  })

  checkpoint_row <- data.table(
    pair_id = row$pair_id,
    exposure = row$exposure,
    outcome = row$outcome,
    direction = row$direction,
    nsnp = nsnp,
    beta = result$beta,
    se = result$se,
    pval = result$pval,
    global_p = result$global_p,
    global_p_raw = result$global_p_raw,
    n_outlier = result$n_outlier,
    analysis = result$analysis,
    status = result$status,
    notes = result$notes
  )
  write_tsv(checkpoint_row, checkpoint)
  cat("presso ", i, "/", nrow(eligible), " ", row$direction,
      " status=", result$status,
      " analysis=", result$analysis,
      " n_outlier=", result$n_outlier,
      "\n", sep = "")
}

checkpoints <- rbindlist(
  lapply(sort(list.files(checkpoint_dir, pattern = "\\.tsv$", full.names = TRUE)), fread),
  fill = TRUE
)
write_tsv(checkpoints, file.path(mr_dir, "mr_presso_qc.tsv"))

mr <- fread(file.path(mr_dir, "mr_results_all.tsv"), showProgress = FALSE)
sens <- fread(file.path(mr_dir, "mr_sensitivity_all.tsv"), showProgress = FALSE)
if (!"mr_presso_analysis" %in% names(sens)) sens[, mr_presso_analysis := NA_character_]
if (!"mr_presso_notes" %in% names(sens)) sens[, mr_presso_notes := NA_character_]
if (!"mr_presso_global_p_raw" %in% names(sens)) sens[, mr_presso_global_p_raw := NA_character_]

for (i in seq_len(nrow(checkpoints))) {
  cp <- checkpoints[i]
  key <- mr$pair_id == cp$pair_id & mr$exposure == cp$exposure &
    mr$outcome == cp$outcome & mr$method == "MR-PRESSO"
  mr[key, `:=`(
    beta = cp$beta,
    se = cp$se,
    pval = cp$pval,
    or = exp(cp$beta),
    or_lci95 = exp(cp$beta - 1.96 * cp$se),
    or_uci95 = exp(cp$beta + 1.96 * cp$se),
    notes = paste0("MR-PRESSO ", cp$analysis, "; ", cp$notes)
  )]
  sens_key <- sens$pair_id == cp$pair_id & sens$exposure == cp$exposure &
    sens$outcome == cp$outcome
  sens[sens_key, `:=`(
    mr_presso_global_p = cp$global_p,
    n_outlier = cp$n_outlier,
    mr_presso_analysis = cp$analysis,
    mr_presso_notes = cp$notes,
    mr_presso_global_p_raw = cp$global_p_raw
  )]
}

write_tsv(mr, file.path(mr_dir, "mr_results_all.tsv"))
write_tsv(sens, file.path(mr_dir, "mr_sensitivity_all.tsv"))
for (current_pair_id in unique(mr$pair_id)) {
  write_tsv(
    mr[pair_id == current_pair_id],
    file.path(per_pair_dir, paste0(current_pair_id, ".mr_results.tsv"))
  )
  write_tsv(
    sens[pair_id == current_pair_id],
    file.path(per_pair_dir, paste0(current_pair_id, ".mr_sensitivity.tsv"))
  )
}

cat("Finished:", as.character(Sys.time()), "\n")
cat("eligible_directions:", nrow(eligible),
    " success:", sum(checkpoints$status == "success"),
    " failed:", sum(checkpoints$status != "success"), "\n")
sink()
