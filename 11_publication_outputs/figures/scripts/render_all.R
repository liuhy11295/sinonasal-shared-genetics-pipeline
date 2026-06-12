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

# Figure 4 currently builds shared main-figure components as an internal
# dependency. Render it first so the dedicated final Figure 1-3 scripts are
# the last writers of their respective artifacts.
run_script("render_figure4.R")
run_script("render_figure1.R")
run_script("render_figure2.R")
run_script("render_figure3.R")
run_script("render_supplementary_figures.R")

message("All final figures rendered under: ",
        file.path(PROJECT_ROOT, "figure_final"))
