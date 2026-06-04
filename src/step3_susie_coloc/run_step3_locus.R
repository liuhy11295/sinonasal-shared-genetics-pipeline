#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("PROJECT_ROOT", "/platform_data/p_user/p010/phase0")
task_id <- if (length(args) >= 2) as.integer(args[[2]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

step2_dir <- file.path(project_root, "results/phase0_extension/step2_input_tables")
pkg_dir <- file.path(project_root, "results/phase0_server/nasal4_vs_other_package")
out_dir <- file.path(project_root, "results/phase0_extension/step3_susie_coloc")
per_locus_dir <- file.path(out_dir, "per_locus")
tmp_dir <- file.path(out_dir, "tmp", sprintf("task_%s", task_id))
qc_dir <- file.path(out_dir, "qc")
dir.create(per_locus_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

plink <- Sys.getenv("PLINK_BIN", file.path(Sys.getenv("CONDA_PREFIX", ""), "bin/plink"))
if (!nzchar(plink) || !file.exists(plink)) plink <- "plink"
ld_ref_dir <- Sys.getenv(
  "PHASE0_LD_PLINK_DIR",
  file.path(project_root, "data/reference/lava/1000G_Phase3_plinkfiles/1000G_EUR_Phase3_plink")
)
ld_ref_prefix_template <- file.path(ld_ref_dir, "1000G.EUR.QC.%s")
min_snps <- as.integer(Sys.getenv("SUSIE_MIN_SNPS", "20"))
max_l <- as.integer(Sys.getenv("SUSIE_L", "10"))

write_tsv <- function(dt, path) {
  fwrite(dt, path, sep = "\t", na = "", quote = FALSE)
}

append_qc <- function(meta, status, reason, n_snps = NA_integer_, n_ref = NA_integer_,
                      n_shared = NA_integer_, n_cs1 = NA_integer_, n_cs2 = NA_integer_) {
  qc <- data.table(
    pair_id = meta$pair_id,
    trait1 = meta$trait1,
    trait2 = meta$trait2,
    locus_id = meta$locus_id,
    status = status,
    reason = reason,
    n_snps = n_snps,
    n_ref_snps = n_ref,
    n_shared_snps = n_shared,
    n_cs_trait1 = n_cs1,
    n_cs_trait2 = n_cs2
  )
  write_tsv(qc, file.path(per_locus_dir, paste0(meta$locus_id, ".qc.tsv")))
}

abs_path <- function(path) {
  if (!nzchar(path)) return("")
  if (startsWith(path, "/")) return(path)
  file.path(project_root, path)
}

fread_any <- function(path, ...) {
  path <- abs_path(path)
  if (grepl("\\.gz$", path)) {
    return(fread(cmd = paste("gzip -dc", shQuote(path)), ...))
  }
  fread(path, ...)
}

clean_region <- function(path) {
  dt <- fread_any(path, showProgress = FALSE)
  required <- c("SNP", "A1", "A2", "BETA", "SE", "P", "N", "CHR", "BP")
  missing <- setdiff(required, names(dt))
  if (length(missing)) stop("missing columns: ", paste(missing, collapse = ","))
  if (!"Z" %in% names(dt)) {
    dt[["Z"]] <- as.numeric(dt[["BETA"]]) / as.numeric(dt[["SE"]])
  }
  needed <- c("SNP", "A1", "A2", "BETA", "SE", "P", "Z", "N", "CHR", "BP")
  dt <- dt[, ..needed]
  dt[, `:=`(
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P),
    Z = as.numeric(Z),
    N = as.numeric(N),
    CHR = as.integer(CHR),
    BP = as.integer(BP)
  )]
  dt <- dt[!duplicated(SNP)]
  dt[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(Z) & is.finite(N) & N > 0 & nzchar(SNP)]
}

make_susie_rows <- function(fit, meta, trait_label, snp_info) {
  if (is.null(fit$pip) || !length(fit$pip)) {
    return(data.table())
  }
  cs_list <- fit$sets$cs
  cs_membership <- vector("list", length(fit$pip))
  cs_leads <- rep(FALSE, length(fit$pip))
  if (!is.null(cs_list) && length(cs_list)) {
    for (cs_name in names(cs_list)) {
      idx <- cs_list[[cs_name]]
      if (!length(idx)) next
      lead_idx <- idx[which.max(fit$pip[idx])]
      cs_leads[lead_idx] <- TRUE
      for (i in idx) cs_membership[[i]] <- c(cs_membership[[i]], cs_name)
    }
  }
  data.table(
    pair_id = meta$pair_id,
    trait1 = meta$trait1,
    trait2 = meta$trait2,
    locus_id = meta$locus_id,
    CHR = meta$CHR,
    START = meta$START,
    END = meta$END,
    trait = trait_label,
    credible_set_id = vapply(cs_membership, function(x) paste(x, collapse = ";"), character(1)),
    SNP = snp_info$SNP,
    PIP = fit$pip,
    is_cs_lead = cs_leads
  )
}

get_cs_snps <- function(fit, snps, idx) {
  if (is.null(fit$sets$cs) || !length(fit$sets$cs) || is.na(idx)) return(character())
  key <- as.character(idx)
  if (key %in% names(fit$sets$cs)) return(snps[fit$sets$cs[[key]]])
  if (idx <= length(fit$sets$cs)) return(snps[fit$sets$cs[[idx]]])
  character()
}

coloc_susie_safe <- function(s1, s2) {
  cs1 <- s1$sets
  cs2 <- s2$sets
  if (is.null(cs1$cs) || is.null(cs2$cs) || length(cs1$cs) == 0 || length(cs2$cs) == 0) {
    return(list(summary = data.table(nsnps = NA)))
  }
  idx1 <- cs1$cs_index
  idx2 <- cs2$cs_index
  bf1 <- s1$lbf_variable[idx1, , drop = FALSE]
  bf2 <- s2$lbf_variable[idx2, , drop = FALSE]
  ret <- coloc:::coloc.bf_bf(bf1, bf2)
  ret$summary <- as.data.table(ret$summary)
  if ("idx1" %in% names(ret$summary)) ret$summary$idx1 <- idx1[ret$summary$idx1]
  if ("idx2" %in% names(ret$summary)) ret$summary$idx2 <- idx2[ret$summary$idx2]
  ret
}

fail_empty_outputs <- function(meta, status, reason) {
  write_tsv(data.table(
    pair_id = character(), trait1 = character(), trait2 = character(), locus_id = character(),
    CHR = integer(), START = integer(), END = integer(), trait = character(),
    credible_set_id = character(), SNP = character(), PIP = numeric(), is_cs_lead = logical()
  ), file.path(per_locus_dir, paste0(meta$locus_id, ".susie_finemap.tsv")))
  write_tsv(data.table(
    pair_id = character(), trait1 = character(), trait2 = character(), locus_id = character(),
    shared_signal_id = character(), SNP = character(), PIP_trait1 = numeric(),
    PIP_trait2 = numeric(), PP.H4 = numeric()
  ), file.path(per_locus_dir, paste0(meta$locus_id, ".coloc_susie.tsv")))
  append_qc(meta, status, reason)
}

coloc_pos <- fread(file.path(step2_dir, "coloc_positive_loci.tsv"), showProgress = FALSE)
coloc_pos[, PP.H4 := as.numeric(PP.H4)]
coloc_pos <- coloc_pos[PP.H4 >= 0.5]
setorder(coloc_pos, pair_id, locus_id)
if (is.na(task_id) || task_id < 1 || task_id > nrow(coloc_pos)) {
  stop("SLURM_ARRAY_TASK_ID outside coloc-positive locus range")
}

meta <- as.list(coloc_pos[task_id])
locus_manifest <- fread(file.path(pkg_dir, "phase0_v3_coloc_locus_manifest.tsv"), showProgress = FALSE)
loc <- locus_manifest[locus_id == meta$locus_id & pair_id == meta$pair_id]
if (nrow(loc) != 1) {
  fail_empty_outputs(meta, "failed", "locus_manifest_row_missing_or_duplicate")
  quit(status = 0)
}

tryCatch({
  a <- clean_region(loc$region_a_file)
  b <- clean_region(loc$region_b_file)
  setnames(a, c("A1", "A2", "BETA", "SE", "P", "Z", "N", "CHR", "BP"),
           paste0(c("A1", "A2", "BETA", "SE", "P", "Z", "N", "CHR", "BP"), "_1"))
  setnames(b, c("A1", "A2", "BETA", "SE", "P", "Z", "N", "CHR", "BP"),
           paste0(c("A1", "A2", "BETA", "SE", "P", "Z", "N", "CHR", "BP"), "_2"))
  m <- merge(a, b, by = "SNP")
  m <- m[(A1_1 == A1_2 & A2_1 == A2_2) | (A1_1 == A2_2 & A2_1 == A1_2)]
  flip <- m$A1_1 == m$A2_2 & m$A2_1 == m$A1_2
  if (any(flip)) {
    m$BETA_2[flip] <- -m$BETA_2[flip]
    m$Z_2[flip] <- -m$Z_2[flip]
  }
  if (nrow(m) < min_snps) {
    fail_empty_outputs(meta, "failed", "SNP数量不足_after_trait_intersection", n_snps = nrow(m))
    quit(status = 0)
  }

  chr <- as.integer(meta$CHR)
  bfile <- sprintf(ld_ref_prefix_template, chr)
  if (!file.exists(paste0(bfile, ".bed")) || !file.exists(paste0(bfile, ".bim")) || !file.exists(paste0(bfile, ".fam"))) {
    fail_empty_outputs(meta, "failed", "LD缺失_reference_bfile_missing", n_snps = nrow(m))
    quit(status = 0)
  }
  bim <- fread(paste0(bfile, ".bim"), header = FALSE, showProgress = FALSE)
  ref_snps <- unique(bim$V2)
  m <- m[SNP %in% ref_snps]
  if (nrow(m) < min_snps) {
    fail_empty_outputs(meta, "failed", "SNP数量不足_after_LD_reference_intersection",
                       n_snps = nrow(m), n_ref = length(ref_snps), n_shared = nrow(m))
    quit(status = 0)
  }

  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  extract_path <- file.path(tmp_dir, paste0(meta$locus_id, ".snps.txt"))
  fwrite(data.table(SNP = m$SNP), extract_path, col.names = FALSE)
  subset_prefix <- file.path(tmp_dir, paste0(meta$locus_id, ".subset"))
  ld_prefix <- file.path(tmp_dir, paste0(meta$locus_id, ".ld"))
  cmd1 <- sprintf(
    "%s --bfile %s --extract %s --make-bed --out %s --allow-no-sex --keep-allele-order",
    shQuote(plink), shQuote(bfile), shQuote(extract_path), shQuote(subset_prefix)
  )
  status1 <- system(cmd1, ignore.stdout = TRUE, ignore.stderr = FALSE)
  if (status1 != 0 || !file.exists(paste0(subset_prefix, ".bim"))) {
    fail_empty_outputs(meta, "failed", "LD构建失败_plink_make_bed", n_snps = nrow(m), n_ref = length(ref_snps))
    quit(status = 0)
  }
  cmd2 <- sprintf(
    "%s --bfile %s --r square gz --out %s --allow-no-sex",
    shQuote(plink), shQuote(subset_prefix), shQuote(ld_prefix)
  )
  status2 <- system(cmd2, ignore.stdout = TRUE, ignore.stderr = FALSE)
  if (status2 != 0 || !file.exists(paste0(ld_prefix, ".ld.gz"))) {
    fail_empty_outputs(meta, "failed", "LD构建失败_plink_r_square", n_snps = nrow(m), n_ref = length(ref_snps))
    quit(status = 0)
  }
  order_bim <- fread(paste0(subset_prefix, ".bim"), header = FALSE, showProgress = FALSE)
  snp_order <- order_bim$V2
  m <- m[match(snp_order, SNP)]
  if (any(is.na(m$SNP))) {
    fail_empty_outputs(meta, "failed", "LD_order_mismatch", n_snps = length(snp_order), n_ref = length(ref_snps))
    quit(status = 0)
  }
  rmat <- as.matrix(fread_any(paste0(ld_prefix, ".ld.gz"), header = FALSE, showProgress = FALSE))
  storage.mode(rmat) <- "double"
  if (nrow(rmat) != nrow(m) || ncol(rmat) != nrow(m)) {
    fail_empty_outputs(meta, "failed", "LD_matrix_dimension_mismatch", n_snps = nrow(m), n_ref = length(ref_snps))
    quit(status = 0)
  }
  rmat[!is.finite(rmat)] <- 0
  rmat <- (rmat + t(rmat)) / 2
  diag(rmat) <- 1

  n1 <- median(m$N_1, na.rm = TRUE)
  n2 <- median(m$N_2, na.rm = TRUE)
  l_val <- max(1, min(max_l, floor(nrow(m) / 10)))
  rownames(rmat) <- m$SNP
  colnames(rmat) <- m$SNP
  z1 <- m$Z_1
  z2 <- m$Z_2
  names(z1) <- m$SNP
  names(z2) <- m$SNP
  fit1 <- susie_rss(z = z1, R = rmat, n = n1, L = l_val, check_prior = FALSE)
  fit2 <- susie_rss(z = z2, R = rmat, n = n2, L = l_val, check_prior = FALSE)
  cs1_n <- if (is.null(fit1$sets$cs)) 0 else length(fit1$sets$cs)
  cs2_n <- if (is.null(fit2$sets$cs)) 0 else length(fit2$sets$cs)
  susie_rows <- rbindlist(list(
    make_susie_rows(fit1, meta, meta$trait1, data.table(SNP = m$SNP)),
    make_susie_rows(fit2, meta, meta$trait2, data.table(SNP = m$SNP))
  ), fill = TRUE)
  write_tsv(susie_rows, file.path(per_locus_dir, paste0(meta$locus_id, ".susie_finemap.tsv")))

  coloc_rows <- data.table(
    pair_id = character(), trait1 = character(), trait2 = character(), locus_id = character(),
    shared_signal_id = character(), SNP = character(), PIP_trait1 = numeric(),
    PIP_trait2 = numeric(), PP.H4 = numeric()
  )
  if (cs1_n > 0 && cs2_n > 0) {
    cres <- coloc_susie_safe(fit1, fit2)
    if (!is.null(cres$summary) && nrow(cres$summary) > 0 && "PP.H4.abf" %in% names(cres$summary)) {
      summ <- as.data.table(cres$summary)
      for (i in seq_len(nrow(summ))) {
        pp4 <- as.numeric(summ$PP.H4.abf[i])
        idx1 <- as.integer(summ$idx1[i])
        idx2 <- as.integer(summ$idx2[i])
        snps1 <- get_cs_snps(fit1, m$SNP, idx1)
        snps2 <- get_cs_snps(fit2, m$SNP, idx2)
        shared <- intersect(snps1, snps2)
        if (!length(shared)) next
        idx <- match(shared, m$SNP)
        coloc_rows <- rbind(
          coloc_rows,
          data.table(
            pair_id = meta$pair_id,
            trait1 = meta$trait1,
            trait2 = meta$trait2,
            locus_id = meta$locus_id,
            shared_signal_id = paste0("trait1_cs", idx1, "__trait2_cs", idx2),
            SNP = shared,
            PIP_trait1 = fit1$pip[idx],
            PIP_trait2 = fit2$pip[idx],
            PP.H4 = pp4
          ),
          fill = TRUE
        )
      }
    }
  }
  write_tsv(coloc_rows, file.path(per_locus_dir, paste0(meta$locus_id, ".coloc_susie.tsv")))
  status <- if (cs1_n > 0 && cs2_n > 0) "success" else "failed"
  reason <- if (status == "success") "OK" else "credible_set为空"
  append_qc(meta, status, reason, n_snps = nrow(m), n_ref = length(ref_snps),
            n_shared = nrow(m), n_cs1 = cs1_n, n_cs2 = cs2_n)
}, error = function(e) {
  fail_empty_outputs(meta, "failed", paste("区域无法分析", conditionMessage(e), sep = ": "))
})
