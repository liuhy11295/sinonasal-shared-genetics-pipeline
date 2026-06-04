#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("PROJECT_ROOT", "/platform_data/p_user/p010/phase0")
task_id <- if (length(args) >= 2) as.integer(args[[2]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

step2_input_dir <- file.path(project_root, "results/phase0_extension/step2_input_tables")
out_dir <- file.path(project_root, "results/phase0_extension/step2_lava_screen")
per_pair_dir <- file.path(out_dir, "per_pair")
dir.create(per_pair_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(dt, path) fwrite(dt, path, sep = "\t", quote = FALSE, na = "")

targets <- fread(file.path(step2_input_dir, "pair_manifest.tsv"), showProgress = FALSE)
setorder(targets, pair_id)
write_tsv(targets, file.path(out_dir, "lava_input_pairs.tsv"))

if (is.na(task_id) || task_id < 1 || task_id > nrow(targets)) {
  stop("SLURM_ARRAY_TASK_ID outside pair range")
}

target <- targets[task_id]
pair_id <- target$pair_id
safe_pair <- gsub("[^A-Za-z0-9_.-]+", "_", pair_id)
out_path <- file.path(per_pair_dir, paste0(safe_pair, ".lava_local_rg.tsv"))
qc_path <- file.path(per_pair_dir, paste0(safe_pair, ".lava_qc.tsv"))

lava_script <- Sys.getenv("LAVA_LOCAL_RG_SCRIPT", "")
if (!nzchar(lava_script) || !file.exists(lava_script)) {
  write_tsv(data.table(
    pair_id = pair_id,
    trait1 = target$trait1,
    trait2 = target$trait2,
    status = "not_run",
    reason = "Set LAVA_LOCAL_RG_SCRIPT to the project-specific LAVA local-rg runner. LAVA is run for all available pairs, independent of LDSC."
  ), qc_path)
  write_tsv(data.table(), out_path)
  quit(status = 0)
}

cmd <- c(
  lava_script,
  "--trait1", target$trait1,
  "--trait2", target$trait2,
  "--sumstats1", target$trait1_munged,
  "--sumstats2", target$trait2_munged,
  "--out", out_path
)
status <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)
rc <- attr(status, "status")
if (is.null(rc)) rc <- 0
writeLines(status, con = file.path(per_pair_dir, paste0(safe_pair, ".lava.log")))
write_tsv(data.table(
  pair_id = pair_id,
  trait1 = target$trait1,
  trait2 = target$trait2,
  status = ifelse(rc == 0, "PASS", "FAIL"),
  reason = ifelse(rc == 0, "", "LAVA local-rg runner returned non-zero status")
), qc_path)
