options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

find_project_root <- function() {
  configured <- Sys.getenv("PAPER_ROOT")
  if (nzchar(configured)) {
    return(normalizePath(configured, mustWork = TRUE))
  }
  candidates <- c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "..")
  )
  for (candidate in candidates) {
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (dir.exists(file.path(candidate, "figure_final")) &&
        dir.exists(file.path(candidate, "result"))) {
      return(candidate)
    }
  }
  stop("Cannot locate project root containing figure_final and result")
}

ROOT <- find_project_root()
FIG <- file.path(ROOT, "figure_final")
RESULT <- file.path(ROOT, "result")
END <- Sys.getenv("RESULTS_ROOT")
if (!nzchar(END)) stop("Set RESULTS_ROOT to the analysis results directory.")

COL <- c(
  graphite = "#2F3437",
  grey = "#7A858C",
  pale = "#F1F5F6",
  line = "#67B7A8",
  line_dark = "#5E6E9E",
  ar = "#3C5488",
  np = "#E64B35",
  comparator = "#7E6C99",
  chr_light = "#DDE7EE",
  chr_dark = "#8FA5B6",
  teal = "#00A087",
  amber = "#F39B7F",
  blue = "#3C5488",
  rose = "#E64B35",
  green = "#67B7A8",
  violet = "#7E6C99",
  red = "#DC0000",
  absent = "#F6F8F9",
  white = "#FFFFFF"
)
colv <- function(x) unname(COL[x])

pal_anchor <- c("Allergic rhinitis" = colv("ar"),
                "Nasal polyps" = colv("np"))
pal_grade <- c("Gene-A" = colv("teal"),
               "Gene-B" = colv("amber"),
               "A" = colv("teal"),
               "B" = colv("amber"))
pal_status <- c("Bonferroni" = colv("teal"),
                "Nominal" = colv("amber"),
                "Not significant" = colv("grey"),
                "Insufficient" = colv("absent"),
                "Shared" = colv("teal"),
                "Main only" = colv("ar"),
                "Strict only" = colv("amber"))

wrap_axis <- function(x, width = 24) wrap_text(x, width)

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.28, colour = COL["graphite"]),
      axis.ticks = element_line(linewidth = 0.25, colour = COL["graphite"]),
      axis.text = element_text(colour = COL["graphite"]),
      axis.title = element_text(colour = COL["graphite"]),
      plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, colour = COL["grey"]),
      plot.caption = element_text(size = base_size - 1, colour = COL["grey"], hjust = 0),
      plot.margin = margin(6, 8, 6, 8),
      legend.title = element_text(face = "bold", size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.6),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(4.2, "mm"),
      strip.background = element_rect(fill = COL["pale"], colour = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.major.x = element_line(colour = "#EEF2F4", linewidth = 0.18),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}
theme_set(theme_pub())

read_tsv <- function(path, ...) {
  read.delim(path, sep = "\t", quote = "", check.names = FALSE,
             comment.char = "", na.strings = c("", "NA", "NaN"), ...)
}

write_source <- function(x, name) {
  path <- file.path(FIG, "source_data", paste0(name, ".tsv"))
  x_out <- as.data.frame(x, stringsAsFactors = FALSE)
  for (cc in names(x_out)) {
    if (is.factor(x_out[[cc]])) x_out[[cc]] <- as.character(x_out[[cc]])
    if (is.character(x_out[[cc]])) {
      x_out[[cc]] <- gsub("[\r\n]+", " ", x_out[[cc]])
    }
  }
  write.table(x_out, path, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "")
  path
}

short_trait <- function(x) {
  x <- sub("^ALLERGIC_RHINITIS_GCST90038664$", "Allergic rhinitis", x)
  x <- sub("^NASAL_POLYPS_GCST90018883$", "Nasal polyps", x)
  x <- gsub("_GCST[0-9]+", "", x)
  x <- gsub("^H7_", "", x)
  x <- gsub("^H8_", "", x)
  x <- gsub("^J10_", "", x)
  x <- gsub("^K11_", "", x)
  x <- gsub("^L12_", "", x)
  x <- gsub("^M13_", "", x)
  x <- gsub("_EXMORE|_SUGG|_STRICT2|_STRICT|_WIDE", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

pair_partner <- function(x) short_trait(sub("^.*__", "", x))
pair_anchor <- function(x) ifelse(grepl("^ALLERGIC_RHINITIS", x),
                                 "Allergic rhinitis", "Nasal polyps")

wrap_text <- function(x, width = 30) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"),
         character(1))
}

qa_rows <- list()

qa_png <- function(path, artifact, kind, figure, panel) {
  con <- file(path, "rb")
  sig <- readBin(con, what = "raw", n = 8)
  close(con)
  png_signature <- paste(as.integer(sig), collapse = ",") ==
    "137,80,78,71,13,10,26,10"
  pass <- file.info(path)$size > 15000 && png_signature
  qa_rows[[length(qa_rows) + 1]] <<- data.frame(
    artifact = artifact,
    kind = kind,
    figure = figure,
    panel = panel,
    file_size = file.info(path)$size,
    png_signature = png_signature,
    automated_status = ifelse(pass, "pass", "review"),
    stringsAsFactors = FALSE
  )
  invisible(pass)
}

save_artifact <- function(plot, stem, folder, width_mm, height_mm,
                          figure, panel = "", kind = "panel", dpi = 600) {
  out <- file.path(FIG, folder)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(out, stem)
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h,
                       family = "sans")
  print(plot)
  grDevices::dev.off()

  grDevices::png(paste0(base, ".png"), width = width_mm, height = height_mm,
                 units = "mm", res = dpi, type = "cairo")
  print(plot)
  grDevices::dev.off()

  qa_png(paste0(base, ".png"), stem, kind, figure, panel)
  invisible(base)
}

save_panel_set <- function(panels, figure_id, main = TRUE,
                           panel_width = 88, panel_height = 78,
                           composite_width = 183, composite_height = 120,
                           design = NULL, composite_stem = NULL) {
  branch <- if (main) "main" else "supplementary"
  letters <- names(panels)
  for (nm in letters) {
    save_artifact(
      panels[[nm]],
      paste0(figure_id, nm),
      file.path(branch, "panels"),
      panel_width, panel_height,
      figure_id, nm, "panel"
    )
  }
  if (is.null(design)) {
    composite <- patchwork::wrap_plots(panels, nrow = 1) +
      patchwork::plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 10))
  } else {
    composite <- patchwork::wrap_plots(panels, design = design) +
      patchwork::plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 10))
  }
  if (is.null(composite_stem)) {
    composite_stem <- paste0(figure_id, "_", paste0(letters, collapse = ""))
  }
  save_artifact(
    composite,
    composite_stem,
    file.path(branch, "composites"),
    composite_width, composite_height,
    figure_id, paste0(letters, collapse = ""), "composite"
  )
  invisible(composite)
}

copy_composite_to_root <- function(stem, branch = "main") {
  src <- file.path(FIG, branch, "composites", stem)
  for (ext in c("pdf", "png")) {
    source_file <- paste0(src, ".", ext)
    if (file.exists(source_file)) {
      file.copy(source_file, file.path(FIG, basename(source_file)),
                overwrite = TRUE)
    }
  }
}

flow_plot <- function(nodes, edges, title, subtitle = NULL) {
  xr <- range(c(nodes$x, edges$x, edges$xend), na.rm = TRUE)
  yr <- range(c(nodes$y, edges$y, edges$yend), na.rm = TRUE)
  ggplot(nodes) +
    geom_segment(
      data = edges,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.55, colour = COL["grey"],
      arrow = grid::arrow(length = grid::unit(2.3, "mm"), type = "closed")
    ) +
    geom_label(
      aes(x, y, label = label, fill = fill),
      linewidth = 0, colour = COL["graphite"], fontface = "bold",
      size = 3, lineheight = 0.95, label.padding = grid::unit(2.5, "mm")
    ) +
    scale_fill_identity() +
    scale_x_continuous(limits = xr + c(-0.65, 0.65)) +
    scale_y_continuous(limits = yr + c(-0.45, 0.45)) +
    coord_cartesian(clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0),
      plot.subtitle = element_text(size = 7, colour = COL["grey"]),
      plot.margin = margin(8, 10, 8, 10)
    )
}

finish_qa <- function() {
  qa <- if (length(qa_rows)) do.call(rbind, qa_rows) else data.frame()
  out <- file.path(FIG, "qa", "automated_render_QA_other_figures.tsv")
  write.table(qa, out,
              sep = "\t", quote = FALSE, row.names = FALSE)
  qa
}
