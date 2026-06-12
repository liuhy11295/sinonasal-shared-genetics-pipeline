#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

PROJECT_ROOT <- if (dir.exists(file.path(getwd(), "figure_final"))) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
}

run_script <- function(name) {
  path <- file.path(PROJECT_ROOT, "figure_final", "scripts", name)
  message("Running ", path)
  status <- system2("Rscript", shQuote(path))
  if (!identical(status, 0L)) {
    stop("Rendering failed for ", name, call. = FALSE)
  }
}

run_script("render_main_figures_revised.R")
run_script("render_figure3_v2.03_revised.R")
run_script("render_supplementary_figures_revised.R")

message("All revised figures rendered under: ",
        file.path(PROJECT_ROOT, "figure_final"))
