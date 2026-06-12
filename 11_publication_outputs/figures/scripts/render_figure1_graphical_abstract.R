options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(grid))

find_figure_root <- function() {
  candidates <- c(getwd(), file.path(getwd(), "figure_final"),
                  file.path(getwd(), ".."))
  for (candidate in candidates) {
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (dir.exists(file.path(candidate, "scripts")) &&
        dir.exists(file.path(candidate, "source_data")) &&
        dir.exists(file.path(candidate, "main"))) {
      return(candidate)
    }
  }
  stop("Cannot locate figure_final project root")
}

FIG <- find_figure_root()
OUT <- file.path(FIG, "main", "composites")
SRC <- file.path(FIG, "source_data", "Figure1_key_numbers.tsv")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(SRC)) stop("Missing source data: ", SRC)

key <- read.delim(SRC, sep = "\t", quote = "", check.names = FALSE,
                  comment.char = "", na.strings = "")
value_of <- function(metric) {
  value <- key$value[key$metric == metric]
  if (length(value) != 1 || is.na(value)) stop("Invalid metric: ", metric)
  as.integer(value)
}

n_nasal <- value_of("nasal_phenotypes")
n_comparators <- value_of("non_nasal_comparators")
n_pairs <- value_of("positive_disease_pairs")
n_ar <- value_of("allergic_rhinitis_pairs")
n_np <- value_of("nasal_polyps_pairs")
stopifnot(n_ar + n_np == n_pairs)

COL <- c(
  ink = "#29343A",
  muted = "#68767D",
  faint = "#98A4AA",
  border = "#D8E0E4",
  arrow = "#AEB9BE",
  panel = "#FBFCFC",
  soft = "#F4F7F8",
  ar = "#3E78A8",
  ar_soft = "#E8F1F7",
  np = "#D8893D",
  np_soft = "#F9EFE4",
  screen = "#557B8F",
  screen_soft = "#EAF1F4",
  locus = "#80699A",
  locus_soft = "#F0ECF4",
  gene = "#3F8A83",
  gene_soft = "#E8F3F1",
  downstream = "#708970",
  downstream_soft = "#EDF2ED",
  white = "#FFFFFF"
)

gp_fill <- function(fill, col = NA, lwd = 1) {
  gpar(fill = fill, col = col, lwd = lwd, linejoin = "round")
}

text_at <- function(label, x, y, size = 10, col = COL["ink"],
                    face = "plain", just = c("left", "centre"),
                    lineheight = 1.05, rot = 0) {
  grid.text(label, x = unit(x, "npc"), y = unit(y, "npc"),
            just = just, rot = rot,
            gp = gpar(fontfamily = "sans", fontsize = size,
                      fontface = face, col = col, lineheight = lineheight))
}

round_rect <- function(x, y, w, h, fill, col = COL["border"],
                       radius = 0.012, lwd = 0.8) {
  grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                 width = unit(w, "npc"), height = unit(h, "npc"),
                 r = unit(radius, "npc"),
                 gp = gp_fill(fill, col, lwd))
}

line_seg <- function(x0, y0, x1, y1, col = COL["arrow"],
                     lwd = 1.2, arrow = FALSE) {
  grid.lines(unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
             gp = gpar(col = col, lwd = lwd, lineend = "round"),
             arrow = if (arrow) {
               grid::arrow(length = unit(2.5, "mm"), type = "closed")
             } else NULL)
}

circle <- function(x, y, r, fill, col = NA, lwd = 0.7) {
  grid.circle(x = unit(x, "npc"), y = unit(y, "npc"),
              r = unit(r, "npc"), gp = gp_fill(fill, col, lwd))
}

chip <- function(label, x, y, w, fill, col = COL["ink"], size = 8.1) {
  round_rect(x, y, w, 0.033, fill, col = NA, radius = 0.009)
  text_at(label, x, y, size = size, col = col,
          face = "bold", just = "centre")
}

bullet <- function(label, x, y, width, col = COL["muted"],
                   dot_col = COL["faint"], size = 8.8, face = "plain") {
  circle(x + 0.004, y + 0.001, 0.0033, dot_col)
  text_at(label, x + 0.013, y, size = size, col = col, face = face,
          just = c("left", "centre"), lineheight = 1.02)
}

draw_gwas_icon <- function(cx, cy) {
  round_rect(cx - 0.031, cy, 0.048, 0.095, COL["white"],
             col = COL["border"], radius = 0.005, lwd = 0.8)
  grid.polygon(
    x = unit(c(cx - 0.007, cx + 0.017, cx + 0.017), "npc"),
    y = unit(c(cy + 0.047, cy + 0.023, cy + 0.047), "npc"),
    gp = gp_fill(COL["soft"], COL["border"], 0.7)
  )
  for (yy in c(cy + 0.018, cy, cy - 0.018)) {
    line_seg(cx - 0.046, yy, cx - 0.018, yy, COL["border"], 1)
  }
  circle(cx + 0.015, cy + 0.018, 0.012, COL["ar"])
  circle(cx + 0.045, cy - 0.012, 0.012, COL["np"])
  circle(cx + 0.017, cy - 0.036, 0.007, "#CBD3D7")
  line_seg(cx - 0.005, cy - 0.005, cx + 0.008, cy + 0.011,
           COL["arrow"], 0.9)
  line_seg(cx - 0.002, cy - 0.020, cx + 0.036, cy - 0.014,
           COL["arrow"], 0.9)
}

draw_network_icon <- function(cx, cy) {
  xs <- c(cx - 0.047, cx - 0.016, cx + 0.023, cx + 0.052, cx + 0.006)
  ys <- c(cy + 0.018, cy + 0.045, cy + 0.025, cy - 0.018, cy - 0.042)
  edges <- rbind(c(1, 2), c(2, 3), c(3, 4), c(3, 5),
                 c(1, 5), c(2, 5))
  for (i in seq_len(nrow(edges))) {
    a <- edges[i, 1]; b <- edges[i, 2]
    line_seg(xs[a], ys[a], xs[b], ys[b], "#C8D3D8", 1)
  }
  fills <- c(COL["ar"], "#A9BAC2", "#7795A4", COL["np"], "#A9BAC2")
  radii <- c(0.013, 0.010, 0.014, 0.013, 0.010)
  for (i in seq_along(xs)) circle(xs[i], ys[i], radii[i], fills[i])
}

draw_locus_icon <- function(cx, cy) {
  x <- seq(cx - 0.058, cx + 0.058, length.out = 9)
  h <- c(0.014, 0.025, 0.017, 0.040, 0.073, 0.035, 0.022, 0.030, 0.015)
  line_seg(cx - 0.064, cy - 0.040, cx + 0.064, cy - 0.040,
           COL["border"], 1)
  for (i in seq_along(x)) {
    grid.rect(x = unit(x[i], "npc"),
              y = unit(cy - 0.040 + h[i] / 2, "npc"),
              width = unit(0.008, "npc"), height = unit(h[i], "npc"),
              gp = gp_fill(if (i == 5) COL["locus"] else "#CFC5DA", NA))
  }
  round_rect(cx, cy - 0.060, 0.112, 0.018, COL["locus_soft"],
             col = COL["locus"], radius = 0.004, lwd = 0.8)
  grid.rect(x = unit(cx + 0.005, "npc"), y = unit(cy - 0.060, "npc"),
            width = unit(0.025, "npc"), height = unit(0.018, "npc"),
            gp = gp_fill(COL["locus"], NA))
}

draw_gene_icon <- function(cx, cy) {
  gene_x <- c(cx - 0.052, cx - 0.020, cx + 0.012)
  gene_y <- c(cy + 0.025, cy - 0.030, cy + 0.032)
  path_x <- c(cx + 0.052, cx + 0.058, cx + 0.035)
  path_y <- c(cy + 0.032, cy - 0.020, cy - 0.050)
  for (i in seq_along(gene_x)) {
    for (j in seq_along(path_x)) {
      if ((i + j) %% 2 == 0 || j == 2) {
        line_seg(gene_x[i], gene_y[i], path_x[j], path_y[j],
                 "#B9D3CF", 0.9)
      }
    }
  }
  for (i in seq_along(gene_x)) circle(gene_x[i], gene_y[i], 0.012,
                                       COL["gene"])
  for (j in seq_along(path_x)) circle(path_x[j], path_y[j],
                                      0.016 - j * 0.0015, "#B9DCD6",
                                      COL["gene"], 0.7)
}

draw_downstream_icon <- function(cx, cy) {
  round_rect(cx - 0.030, cy + 0.015, 0.082, 0.035, COL["white"],
             col = COL["downstream"], radius = 0.018, lwd = 1)
  grid.rect(x = unit(cx - 0.050, "npc"), y = unit(cy + 0.015, "npc"),
            width = unit(0.041, "npc"), height = unit(0.033, "npc"),
            gp = gp_fill(COL["downstream"], NA))
  line_seg(cx - 0.054, cy - 0.045, cx + 0.054, cy - 0.045,
           COL["downstream"], 1.4, arrow = TRUE)
  line_seg(cx + 0.038, cy - 0.061, cx - 0.054, cy - 0.061,
           COL["faint"], 1.1, arrow = TRUE)
}

draw_card <- function(index, cx, title, accent, soft, fig_ref = NULL) {
  card_w <- 0.172
  card_y <- 0.575
  card_h <- 0.585
  round_rect(cx, card_y, card_w, card_h, COL["panel"],
             col = COL["border"], radius = 0.012, lwd = 0.9)
  grid.roundrect(x = unit(cx, "npc"), y = unit(0.853, "npc"),
                 width = unit(card_w, "npc"), height = unit(0.030, "npc"),
                 r = unit(0.010, "npc"), gp = gp_fill(accent, NA))
  circle(cx - card_w / 2 + 0.021, 0.823, 0.0135, accent)
  text_at(as.character(index), cx - card_w / 2 + 0.021, 0.823,
          size = 8.5, col = COL["white"], face = "bold", just = "centre")
  text_at(title, cx - card_w / 2 + 0.040, 0.823, size = 11.0,
          face = "bold", just = c("left", "centre"), lineheight = 0.98)
  if (!is.null(fig_ref)) {
    chip(fig_ref, cx + card_w / 2 - 0.027, 0.332, 0.040,
         soft, col = accent, size = 7.5)
  }
}

draw_flow_arrow <- function(x0, x1) {
  line_seg(x0, 0.575, x1, 0.575, COL["arrow"], 1.25, arrow = TRUE)
}

draw_summary_cell <- function(x0, x1, number = NULL, label,
                              accent = COL["muted"], fill = COL["white"]) {
  cx <- (x0 + x1) / 2
  grid.rect(x = unit(cx, "npc"), y = unit(0.155, "npc"),
            width = unit(x1 - x0, "npc"), height = unit(0.132, "npc"),
            gp = gp_fill(fill, COL["border"], 0.7))
  if (!is.null(number)) {
    text_at(number, cx, 0.176, size = 17, col = accent,
            face = "bold", just = "centre")
    text_at(label, cx, 0.128, size = 8.3, col = COL["muted"],
            face = "bold", just = "centre", lineheight = 1.0)
  } else {
    circle(cx, 0.182, 0.009, accent)
    text_at(label, cx, 0.137, size = 8.6, col = COL["ink"],
            face = "bold", just = "centre", lineheight = 1.0)
  }
}

phase_box <- function(letter, title, x, y, w, h, accent, fill) {
  grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                 width = unit(w, "npc"), height = unit(h, "npc"),
                 r = unit(0.013, "npc"),
                 gp = gpar(fill = fill, col = accent, lwd = 1.1,
                           lty = "22", linejoin = "round"))
  circle(x - w / 2 + 0.020, y + h / 2 - 0.028, 0.014, accent)
  text_at(letter, x - w / 2 + 0.020, y + h / 2 - 0.028,
          10, COL["white"], "bold", "centre")
  text_at(title, x - w / 2 + 0.042, y + h / 2 - 0.028,
          10.3, COL["ink"], "bold", c("left", "centre"))
}

module_title <- function(title, x, y, w, accent, fig_ref = NULL) {
  round_rect(x, y, w, 0.052, COL["white"], accent, 0.008, 0.9)
  grid.rect(x = unit(x - w / 2 + 0.006, "npc"), y = unit(y, "npc"),
            width = unit(0.012, "npc"), height = unit(0.052, "npc"),
            gp = gp_fill(accent, NA))
  text_at(title, x - w / 2 + 0.021, y, 9.1, COL["ink"], "bold",
          c("left", "centre"), 0.98)
  if (!is.null(fig_ref)) {
    text_at(fig_ref, x + w / 2 - 0.010, y, 7.4, accent, "bold",
            c("right", "centre"))
  }
}

draw_correlation_icon <- function(cx, cy) {
  line_seg(cx - 0.050, cy - 0.037, cx - 0.050, cy + 0.040,
           COL["ink"], 0.9)
  line_seg(cx - 0.050, cy - 0.037, cx + 0.052, cy - 0.037,
           COL["ink"], 0.9)
  line_seg(cx - 0.038, cy - 0.025, cx + 0.042, cy + 0.030,
           COL["screen"], 1.2)
  pts <- rbind(
    c(-0.037, -0.020), c(-0.022, -0.004), c(-0.005, 0.003),
    c(0.012, 0.019), c(0.030, 0.020), c(0.043, 0.034)
  )
  for (i in seq_len(nrow(pts))) {
    circle(cx + pts[i, 1], cy + pts[i, 2], 0.0045, COL["screen"])
  }
  line_seg(cx + 0.069, cy - 0.035, cx + 0.069, cy + 0.038,
           COL["border"], 0.8)
  vals <- c(0.012, 0.025, 0.017, 0.040, 0.027)
  for (i in seq_along(vals)) {
    grid.rect(x = unit(cx + 0.082 + (i - 1) * 0.014, "npc"),
              y = unit(cy - 0.035 + vals[i] / 2, "npc"),
              width = unit(0.008, "npc"), height = unit(vals[i], "npc"),
              gp = gp_fill(if (i == 4) COL["np"] else "#AFC2CC", NA))
  }
}

draw_cross_trait_icon <- function(cx, cy) {
  xs <- seq(cx - 0.075, cx + 0.010, length.out = 11)
  hs1 <- c(.010, .018, .012, .030, .015, .052, .020, .014, .026, .012, .018)
  hs2 <- c(.015, .012, .020, .018, .028, .047, .014, .022, .011, .030, .016)
  for (i in seq_along(xs)) {
    grid.rect(x = unit(xs[i], "npc"),
              y = unit(cy + 0.021 + hs1[i] / 2, "npc"),
              width = unit(0.0055, "npc"), height = unit(hs1[i], "npc"),
              gp = gp_fill(if (i == 6) COL["locus"] else "#C9C0D4", NA))
    grid.rect(x = unit(xs[i], "npc"),
              y = unit(cy - 0.022 - hs2[i] / 2, "npc"),
              width = unit(0.0055, "npc"), height = unit(hs2[i], "npc"),
              gp = gp_fill(if (i == 6) COL["locus"] else "#D9D2E0", NA))
  }
  line_seg(cx - 0.080, cy + 0.020, cx + 0.015, cy + 0.020,
           COL["border"], 0.8)
  line_seg(cx - 0.080, cy - 0.020, cx + 0.015, cy - 0.020,
           COL["border"], 0.8)
  diamond_x <- c(cx + 0.046, cx + 0.061, cx + 0.076)
  diamond_fill <- c("#E7DDEB", "#B898C1", COL["locus"])
  for (i in seq_along(diamond_x)) {
    grid.polygon(
      x = unit(c(diamond_x[i], diamond_x[i] + 0.010,
                 diamond_x[i], diamond_x[i] - 0.010), "npc"),
      y = unit(c(cy + 0.012, cy, cy - 0.012, cy), "npc"),
      gp = gp_fill(diamond_fill[i], NA)
    )
  }
}

draw_gene_pathway_icon <- function(cx, cy) {
  genes_x <- c(cx - 0.070, cx - 0.042, cx - 0.014)
  genes_y <- c(cy + 0.027, cy - 0.025, cy + 0.005)
  paths_x <- c(cx + 0.040, cx + 0.068, cx + 0.052)
  paths_y <- c(cy + 0.030, cy - 0.002, cy - 0.038)
  for (i in seq_along(genes_x)) {
    for (j in seq_along(paths_x)) {
      if (i == j || j == 2) {
        line_seg(genes_x[i], genes_y[i], paths_x[j], paths_y[j],
                 "#B7D2CE", 0.9)
      }
    }
  }
  for (i in seq_along(genes_x)) {
    circle(genes_x[i], genes_y[i], 0.010, COL["gene"])
  }
  radii <- c(0.016, 0.013, 0.011)
  for (j in seq_along(paths_x)) {
    circle(paths_x[j], paths_y[j], radii[j], "#B9DCD6",
           COL["gene"], 0.8)
  }
}

draw_downstream_combined_icon <- function(cx, cy) {
  round_rect(cx - 0.047, cy + 0.026, 0.074, 0.032, COL["white"],
             COL["downstream"], 0.016, 1)
  grid.rect(x = unit(cx - 0.065, "npc"), y = unit(cy + 0.026, "npc"),
            width = unit(0.037, "npc"), height = unit(0.030, "npc"),
            gp = gp_fill(COL["downstream"], NA))
  text_at("drug target", cx + 0.018, cy + 0.026, 6.7,
          COL["downstream"], "bold", "centre")
  round_rect(cx - 0.041, cy - 0.039, 0.050, 0.034,
             COL["downstream_soft"], COL["downstream"], 0.005, 0.8)
  round_rect(cx + 0.041, cy - 0.039, 0.050, 0.034,
             COL["screen_soft"], COL["screen"], 0.005, 0.8)
  text_at("nasal", cx - 0.041, cy - 0.039, 6.5,
          COL["downstream"], "bold", "centre")
  text_at("non-nasal", cx + 0.041, cy - 0.039, 6.5,
          COL["screen"], "bold", "centre")
  line_seg(cx - 0.014, cy - 0.034, cx + 0.014, cy - 0.034,
           COL["downstream"], 1.0, TRUE)
  line_seg(cx + 0.014, cy - 0.047, cx - 0.014, cy - 0.047,
           COL["faint"], 0.9, TRUE)
}

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gp_fill(COL["white"], NA))

  text_at("From genetic comorbidity screening to shared biological interpretation",
          0.027, 0.958, size = 16.5, face = "bold",
          just = c("left", "centre"))
  text_at("Study workflow and evidence hierarchy",
          0.027, 0.922, size = 9.3, col = COL["muted"],
          just = c("left", "centre"))
  line_seg(0.027, 0.895, 0.973, 0.895, COL["border"], 0.9)

  phase_y <- 0.590
  phase_h <- 0.540
  phase_box("a", "GWAS input and genetic comorbidity screening",
            0.205, phase_y, 0.356, phase_h, COL["screen"], "#F8FAFB")
  phase_box("b", "Shared genetic evidence and biological interpretation",
            0.602, phase_y, 0.382, phase_h, COL["locus"], "#FCFBFD")
  phase_box("c", "Exploratory downstream\nevidence",
            0.891, phase_y, 0.166, phase_h, COL["downstream"], "#FAFBFA")

  line_seg(0.386, phase_y, 0.407, phase_y, COL["arrow"], 1.6, TRUE)
  line_seg(0.795, phase_y, 0.806, phase_y, COL["arrow"], 1.6, TRUE)

  # Stage a: input and correlation screening.
  module_title("GWAS summary statistics", 0.112, 0.780, 0.145,
               COL["screen"])
  draw_gwas_icon(0.112, 0.695)
  bullet(paste0(n_nasal, " nasal phenotypes screened"),
         0.046, 0.620, 0.13, size = 7.9)
  bullet(paste0(n_comparators, " non-nasal comparator traits"),
         0.046, 0.588, 0.13, size = 7.9)
  bullet("European or European-enriched GWAS",
         0.046, 0.556, 0.13, size = 7.6)
  chip("Allergic rhinitis", 0.087, 0.505, 0.073,
       COL["ar_soft"], COL["ar"], 7.7)
  chip("Nasal polyps", 0.165, 0.505, 0.063,
       COL["np_soft"], COL["np"], 7.7)
  text_at("Other nasal traits screened but not retained",
          0.112, 0.468, 6.7, COL["faint"], "plain", "centre")

  module_title("Genetic comorbidity screening", 0.294, 0.780, 0.145,
               COL["screen"], "Fig. 2")
  draw_correlation_icon(0.277, 0.694)
  text_at("LDSC global  |  LAVA local correlation",
          0.294, 0.621, 7.6, COL["muted"], "bold", "centre")
  round_rect(0.294, 0.552, 0.136, 0.094, COL["screen_soft"],
             NA, 0.008)
  text_at(n_pairs, 0.294, 0.572, 16.0, COL["screen"], "bold", "centre")
  text_at("positive disease pairs", 0.294, 0.533, 8.1,
          COL["ink"], "bold", "centre")
  chip(paste0(n_ar, " AR-related"), 0.259, 0.483, 0.064,
       COL["ar_soft"], COL["ar"], 7.5)
  chip(paste0(n_np, " NP-related"), 0.330, 0.483, 0.064,
       COL["np_soft"], COL["np"], 7.5)

  # Stage b: locus evidence and biological interpretation.
  module_title("Shared SNP and locus evidence", 0.510, 0.780, 0.175,
               COL["locus"])
  draw_cross_trait_icon(0.510, 0.690)
  chip("PLACO", 0.453, 0.622, 0.043,
       COL["locus_soft"], COL["locus"], 7.2)
  chip("CPASSOC", 0.505, 0.622, 0.052,
       COL["locus_soft"], COL["locus"], 7.2)
  chip("coloc", 0.553, 0.622, 0.037,
       COL["locus_soft"], COL["locus"], 7.2)
  chip("SuSiE coloc", 0.481, 0.582, 0.066,
       COL["locus_soft"], COL["locus"], 7.2)
  chip("MTAG", 0.550, 0.582, 0.043,
       COL["locus_soft"], COL["locus"], 7.2)
  text_at("cross-trait SNP and locus evidence",
          0.510, 0.548, 7.1, COL["locus"], "bold", "centre")
  bullet("candidate pleiotropic variants",
         0.437, 0.515, 0.15, size = 7.8, dot_col = COL["locus"])
  bullet("shared and recurrent loci",
         0.437, 0.485, 0.15, size = 7.8, dot_col = COL["locus"])
  bullet("locus-to-gene evidence chains",
         0.437, 0.455, 0.15, size = 7.8, dot_col = COL["locus"])

  line_seg(0.606, 0.590, 0.630, 0.590, COL["arrow"], 1.2, TRUE)

  module_title("Gene and pathway interpretation", 0.704, 0.780, 0.175,
               COL["gene"], "Fig. 3")
  draw_gene_pathway_icon(0.704, 0.694)
  text_at("MAGMA gene analysis  |  FUSION TWAS",
          0.704, 0.626, 7.5, COL["gene"], "bold", "centre")
  text_at("MAGMA + TWAS overlap",
          0.704, 0.596, 7.4, COL["muted"], "bold", "centre")
  text_at("ORA: GO BP / KEGG / Reactome",
          0.704, 0.566, 7.4, COL["muted"], "bold", "centre")
  round_rect(0.704, 0.501, 0.162, 0.088, COL["gene_soft"],
             NA, 0.008)
  text_at("prioritized genes  |  recurrent pathways",
          0.704, 0.523, 7.6, COL["gene"], "bold", "centre")
  text_at("immune, inflammatory and epithelial-barrier biology",
          0.704, 0.483, 7.3, COL["ink"], "bold", "centre")

  # Stage c: deliberately quieter downstream layer.
  module_title("Downstream annotation", 0.890, 0.780, 0.146,
               COL["downstream"], "Fig. 4")
  draw_downstream_combined_icon(0.890, 0.690)
  text_at("approved-drug target annotation",
          0.890, 0.619, 7.6, COL["muted"], "bold", "centre")
  text_at("exploratory drug-target annotation",
          0.890, 0.586, 7.4, COL["downstream"], "bold", "centre")
  text_at("LCV  |  bidirectional MR",
          0.890, 0.538, 7.7, COL["ink"], "bold", "centre")
  text_at("secondary directional evidence",
          0.890, 0.505, 7.4, COL["downstream"], "bold", "centre")
  round_rect(0.890, 0.458, 0.154, 0.052, COL["downstream_soft"],
             NA, 0.007)
  text_at("annotations  |  directional summary",
          0.890, 0.458, 7.2, COL["ink"], "bold", "centre")

  text_at("KEY STUDY OUTPUTS", 0.027, 0.248, 7.6, COL["muted"],
          "bold", c("left", "centre"))
  bounds <- seq(0.027, 0.973, length.out = 8)
  draw_summary_cell(bounds[1], bounds[2], n_pairs,
                    "positive disease\npairs", COL["screen"],
                    COL["screen_soft"])
  draw_summary_cell(bounds[2], bounds[3], n_ar,
                    "AR-related\npairs", COL["ar"], COL["ar_soft"])
  draw_summary_cell(bounds[3], bounds[4], n_np,
                    "NP-related\npairs", COL["np"], COL["np_soft"])
  draw_summary_cell(bounds[4], bounds[5], NULL,
                    "shared\nloci", COL["locus"], COL["locus_soft"])
  draw_summary_cell(bounds[5], bounds[6], NULL,
                    "prioritized\ngenes", COL["gene"], COL["gene_soft"])
  draw_summary_cell(bounds[6], bounds[7], NULL,
                    "recurrent\npathways", COL["gene"], COL["gene_soft"])
  draw_summary_cell(bounds[7], bounds[8], NULL,
                    "downstream\nannotations", COL["downstream"],
                    COL["downstream_soft"])

  text_at("Primary shared-genetic evidence", 0.027, 0.050, 7.5,
          COL["muted"], "bold", c("left", "centre"))
  line_seg(0.150, 0.050, 0.720, 0.050, COL["border"], 1.4)
  text_at("Exploratory / secondary evidence", 0.973, 0.050, 7.5,
          COL["downstream"], "bold", c("right", "centre"))
}

export_figure <- function() {
  width_in <- 16.0
  height_in <- 6.6
  pdf_path <- file.path(OUT, "Figure1.pdf")
  png_path <- file.path(OUT, "Figure1.png")

  grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in,
                       family = "sans", onefile = TRUE)
  draw_figure()
  grDevices::dev.off()

  grDevices::png(png_path, width = width_in, height = height_in,
                 units = "in", res = 600, type = "cairo", bg = "white")
  draw_figure()
  grDevices::dev.off()

  file.copy(pdf_path, file.path(FIG, "Figure1.pdf"), overwrite = TRUE)
  file.copy(png_path, file.path(FIG, "Figure1.png"), overwrite = TRUE)

  message("Wrote ", pdf_path)
  message("Wrote ", png_path)
}

export_figure()
