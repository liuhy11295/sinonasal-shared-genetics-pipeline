PROJECT_ROOT <- if (dir.exists(file.path(getwd(), "figure_final"))) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
}
source(file.path(PROJECT_ROOT, "figure_final", "scripts", "common_revised.R"))

label_norm <- function(x) {
  x <- gsub(" to J10 ", " to ", x)
  x <- gsub(" to H8 ", " to ", x)
  x <- gsub(" to H7 ", " to ", x)
  x
}

format_p <- function(x) {
  ifelse(
    is.na(x), "",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
           sprintf("%.3f", x))
  )
}

format_beta <- function(x) {
  ifelse(is.na(x), "", ifelse(abs(x) >= 1, sprintf("%.2f", x),
                              sprintf("%.3f", x)))
}

pal_dir <- c(
  lcv_bonf = "#00A087",
  mr_bonf = "#46679F",
  mr_nominal = "#D79624",
  mr_ns = "#D5DBDF",
  mr_insufficient = "#F5F7F8",
  strict = "#252A2D",
  outline = "#7A858A",
  text = "#2F3437",
  grid = "#EEF2F4"
)

status_cols <- c(
  "MR Bonferroni" = unname(pal_dir["mr_bonf"]),
  "MR nominal" = unname(pal_dir["mr_nominal"]),
  "Not significant" = unname(pal_dir["mr_ns"]),
  "Insufficient" = unname(pal_dir["mr_insufficient"])
)

build_mr_panel_data <- function(direction_ids, role, mr_direction,
                                presso, strict) {
  out <- mr_direction[mr_direction$direction %in% direction_ids, ]
  out$beta <- suppressWarnings(as.numeric(out$beta))
  out$se <- suppressWarnings(as.numeric(out$se))
  out$pval <- suppressWarnings(as.numeric(out$pval))
  out$n_instruments_after_harmonise <- suppressWarnings(
    as.numeric(out$n_instruments_after_harmonise)
  )
  out <- merge(out, presso, by = "direction", all.x = TRUE)
  out$mr_label <- paste(short_trait(out$exposure), "->",
                        short_trait(out$outcome))
  out$mr_label_norm <- label_norm(gsub(" -> ", " to ", out$mr_label))
  out <- merge(out, strict[, c("label_norm", "strict_robust")],
               by.x = "mr_label_norm", by.y = "label_norm", all.x = TRUE)
  out$strict_robust[is.na(out$strict_robust)] <- FALSE
  out$IVW_status <- ifelse(
    out$status != "eligible" | is.na(out$status), "Insufficient",
    ifelse(out$pval < 0.05 / 58, "MR Bonferroni",
           ifelse(out$pval < 0.05, "MR nominal", "Not significant"))
  )
  out$MR_PRESSO_sig <- ifelse(
    out$status != "eligible" | is.na(out$status) | is.na(out$MR_PRESSO_p),
    "MR-PRESSO not available",
    ifelse(out$MR_PRESSO_p < 0.05, "MR-PRESSO significant",
           "MR-PRESSO not significant")
  )
  out$beta_low <- out$beta - 1.96 * out$se
  out$beta_high <- out$beta + 1.96 * out$se
  out$strict_label <- ifelse(out$strict_robust, "Full strict robust",
                             "Not full strict")
  out$IVW_status <- factor(out$IVW_status,
                           levels = c("MR Bonferroni", "MR nominal",
                                      "Not significant", "Insufficient"))
  out$MR_PRESSO_sig <- factor(out$MR_PRESSO_sig,
                              levels = c("MR-PRESSO significant",
                                         "MR-PRESSO not significant",
                                         "MR-PRESSO not available"))
  out$strict_label <- factor(out$strict_label,
                             levels = c("Not full strict",
                                        "Full strict robust"))
  out$mr_p_label <- ifelse(
    out$IVW_status == "Insufficient",
    "insufficient",
    paste0("P=", format_p(out$pval))
  )
  out$beta_label <- ifelse(
    out$IVW_status == "Insufficient",
    "",
    paste0("beta=", format_beta(out$beta))
  )
  out$mr_direction_role <- role
  out
}

plot_mr_column <- function(df, title, x_limits, x_breaks, x_labels,
                           x_transform = "identity",
                           strict_nudge = 0.025) {
  df$row_label <- factor(df$row_label_chr, levels = row_levels)
  df_obs <- df[!is.na(df$beta), ]
  df_miss <- df[is.na(df$beta), ]
  df_strict <- df_obs[df_obs$strict_robust %in% TRUE, ]
  df_strict$strict_x <- pmin(x_limits[2] * 0.97,
                             df_strict$beta_high + strict_nudge)

  ggplot(df, aes(beta, row_label)) +
    geom_vline(xintercept = 0, colour = "#AEB8BF", linewidth = 0.28) +
    geom_segment(
      data = df_obs,
      aes(x = beta_low, xend = beta_high, yend = row_label,
          colour = IVW_status),
      linewidth = 0.38, alpha = 0.82
    ) +
    geom_point(
      data = df_obs,
      aes(fill = IVW_status, shape = MR_PRESSO_sig,
          colour = strict_label),
      size = 2.25, stroke = 0.62
    ) +
    geom_point(
      data = df_strict,
      aes(x = beta, y = row_label),
      inherit.aes = FALSE,
      shape = 21, fill = NA, colour = "#111111",
      size = 3.25, stroke = 0.72
    ) +
    geom_label(
      data = df_strict,
      aes(x = strict_x, label = "strict"),
      hjust = 0, size = 1.55, linewidth = 0.12,
      label.padding = grid::unit(0.65, "mm"),
      fill = "white", colour = "#111111"
    ) +
    geom_text(
      data = df_miss,
      aes(x = 0.045, label = "insufficient"),
      hjust = 0, size = 1.82, colour = COL["grey"]
    ) +
    scale_fill_manual(values = status_cols, drop = FALSE) +
    scale_colour_manual(
      values = c("MR Bonferroni" = unname(pal_dir["mr_bonf"]),
                 "MR nominal" = unname(pal_dir["mr_nominal"]),
                 "Not significant" = unname(pal_dir["outline"]),
                 "Insufficient" = unname(pal_dir["outline"]),
                 "Not full strict" = unname(pal_dir["outline"]),
                 "Full strict robust" = unname(pal_dir["strict"])),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c("MR-PRESSO significant" = 21,
                 "MR-PRESSO not significant" = 24,
                 "MR-PRESSO not available" = 22),
      drop = FALSE
    ) +
    scale_x_continuous(
      trans = x_transform,
      limits = x_limits,
      breaks = x_breaks,
      labels = x_labels,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(title = title, x = "IVW beta", y = NULL, fill = NULL, shape = NULL) +
    theme_pub(5.8) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      panel.grid.major.x = element_line(colour = pal_dir["grid"],
                                        linewidth = 0.22),
      plot.margin = margin(8, 12, 7, 2)
    )
}

lcv <- read_tsv(file.path(ROOT, "figure_v2_source_data", "lcv.tsv"))
mr_direction <- read_tsv(file.path(ROOT, "figure_v2_source_data",
                                   "mr_directions.tsv"))
mr_all <- read_tsv(file.path(END, "LCVMR", "MR", "mr_results_all.tsv"))
rob <- read_tsv(file.path(ROOT, "figure_v2_source_data",
                          "mr_robustness.tsv"))

lcv$LCV_gcp <- suppressWarnings(as.numeric(lcv$LCV_gcp))
lcv$LCV_gcp_se <- suppressWarnings(as.numeric(lcv$LCV_gcp_se))
lcv$LCV_gcp_p <- suppressWarnings(as.numeric(lcv$LCV_gcp_p))
lcv$LCV_ci_low <- lcv$LCV_gcp - 1.96 * lcv$LCV_gcp_se
lcv$LCV_ci_high <- lcv$LCV_gcp + 1.96 * lcv$LCV_gcp_se
lcv$LCV_sig <- ifelse(
  lcv$LCV_gcp_p < 0.05 / nrow(lcv), "LCV Bonferroni",
  ifelse(lcv$LCV_gcp_p < 0.05, "LCV nominal", "Not significant")
)
lcv$lcv_direction <- ifelse(
  lcv$LCV_gcp >= 0,
  paste(short_trait(lcv$trait1), "->", short_trait(lcv$trait2)),
  paste(short_trait(lcv$trait2), "->", short_trait(lcv$trait1))
)
lcv$mr_direction_id <- ifelse(
  lcv$LCV_gcp >= 0,
  paste(lcv$trait1, lcv$trait2, sep = "_to_"),
  paste(lcv$trait2, lcv$trait1, sep = "_to_")
)
lcv$mr_reverse_direction_id <- ifelse(
  lcv$LCV_gcp >= 0,
  paste(lcv$trait2, lcv$trait1, sep = "_to_"),
  paste(lcv$trait1, lcv$trait2, sep = "_to_")
)
lcv_focus <- lcv[lcv$LCV_gcp_p < 0.05 / nrow(lcv), ]
lcv_focus <- lcv_focus[order(lcv_focus$LCV_gcp_p), ]
lcv_focus$row_order <- seq_len(nrow(lcv_focus))
lcv_focus$disease_pair_label <- paste(pair_anchor(lcv_focus$pair_id),
                                      "-", pair_partner(lcv_focus$pair_id))
lcv_focus$row_label_chr <- lcv_focus$disease_pair_label
row_levels <- rev(lcv_focus$row_label_chr)
lcv_focus$row_label <- factor(lcv_focus$row_label_chr, levels = row_levels)
lcv_focus$LCV_p_label <- paste0("P=", format_p(lcv_focus$LCV_gcp_p))

mr_all$pval <- suppressWarnings(as.numeric(mr_all$pval))
all_directions <- unique(c(lcv_focus$mr_direction_id,
                           lcv_focus$mr_reverse_direction_id))
presso <- mr_all[
  mr_all$method %in% "MR-PRESSO" & mr_all$direction %in% all_directions,
  c("direction", "pval")
]
names(presso)[2] <- "MR_PRESSO_p"

rob$pass_logical <- rob$pass %in% c(TRUE, "True", "TRUE", "true")
strict <- aggregate(pass_logical ~ direction, rob, all)
names(strict) <- c("robust_label", "strict_robust")
strict$label_norm <- label_norm(strict$robust_label)

row_map <- lcv_focus[
  , c("row_order", "pair_id", "trait1", "trait2", "lcv_direction",
      "disease_pair_label", "row_label_chr", "row_label", "mr_direction_id",
      "mr_reverse_direction_id")
]

same_mr <- build_mr_panel_data(
  lcv_focus$mr_direction_id, "LCV direction", mr_direction, presso, strict
)
reverse_mr <- build_mr_panel_data(
  lcv_focus$mr_reverse_direction_id, "Reverse direction", mr_direction,
  presso, strict
)

mr_same <- merge(row_map, same_mr, by.x = "mr_direction_id",
                 by.y = "direction", all.x = TRUE)
mr_reverse <- merge(row_map, reverse_mr, by.x = "mr_reverse_direction_id",
                    by.y = "direction", all.x = TRUE)
mr_same <- mr_same[order(mr_same$row_order), ]
mr_reverse <- mr_reverse[order(mr_reverse$row_order), ]
mr_same$row_label <- factor(mr_same$row_label_chr, levels = row_levels)
mr_reverse$row_label <- factor(mr_reverse$row_label_chr, levels = row_levels)

p_labels <- ggplot(lcv_focus, aes(0, row_label)) +
  geom_text(aes(label = row_label_chr), hjust = 0, size = 2.35,
            colour = unname(pal_dir["text"]), lineheight = 0.86) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Disease pair", x = NULL, y = NULL) +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 6.6,
                              colour = COL["graphite"], hjust = 0),
    plot.margin = margin(8, 10, 7, 2)
  )

p_lcv <- ggplot(lcv_focus, aes(LCV_gcp, row_label)) +
  geom_vline(xintercept = 0, colour = "#AEB8BF", linewidth = 0.28) +
  geom_segment(
    aes(x = LCV_ci_low, xend = LCV_ci_high, yend = row_label),
    colour = unname(pal_dir["lcv_bonf"]),
    linewidth = 0.42, alpha = 0.95
  ) +
  geom_point(fill = unname(pal_dir["lcv_bonf"]), shape = 21, size = 2.3,
             colour = "white", stroke = 0.35) +
  geom_text(aes(x = 1.08, label = LCV_p_label), hjust = 0,
            size = 1.95, colour = unname(pal_dir["text"])) +
  scale_x_continuous(limits = c(-0.62, 1.28),
                     breaks = c(-0.5, 0, 0.5, 1.0),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "LCV GCP", x = "GCP", y = NULL) +
  theme_pub(6.2) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid.major.x = element_line(colour = pal_dir["grid"],
                                      linewidth = 0.22),
    plot.margin = margin(8, 20, 7, 2)
  )

p_mr_same <- plot_mr_column(
  mr_same, "MR, LCV direction",
  x_limits = c(-0.025, 0.33),
  x_breaks = c(0, 0.1, 0.2, 0.3),
  x_labels = c("0", "0.1", "0.2", "0.3"),
  x_transform = "identity",
  strict_nudge = 0.020
)
p_mr_reverse <- plot_mr_column(
  mr_reverse, "Reverse MR",
  x_limits = c(-0.20, 20),
  x_breaks = c(0, 0.5, 1, 5, 10, 15),
  x_labels = c("0", "0.5", "1", "5", "10", "15"),
  x_transform = scales::pseudo_log_trans(sigma = 0.35),
  strict_nudge = 0.40
)

p_key_df <- data.frame(
  x = c(0.025, 0.145, 0.270, 0.405, 0.545, 0.675),
  label = c("LCV Bonf.", "MR Bonf.", "MR nominal",
            "MR not sig.", "insufficient", "Full strict"),
  fill = c(unname(pal_dir["lcv_bonf"]), unname(pal_dir["mr_bonf"]),
           unname(pal_dir["mr_nominal"]), unname(pal_dir["mr_ns"]),
           NA, NA),
  stringsAsFactors = FALSE
)
p_key <- ggplot(p_key_df, aes(x, 1)) +
  geom_segment(
    data = p_key_df[1, ],
    aes(x = x - 0.015, xend = x + 0.015, y = 1, yend = 1),
    inherit.aes = FALSE, colour = unname(pal_dir["lcv_bonf"]),
    linewidth = 0.42
  ) +
  geom_point(
    data = p_key_df[1:4, ],
    aes(fill = fill), shape = 21, size = 2.25,
    colour = "#7A858A", stroke = 0.45
  ) +
  geom_text(
    data = p_key_df[5, ],
    aes(label = "insufficient"),
    hjust = 0.5, size = 1.95, colour = COL["grey"]
  ) +
  geom_point(
    data = p_key_df[6, ],
    shape = 21, fill = NA, colour = "#111111",
    size = 2.8, stroke = 0.75
  ) +
  geom_text(data = p_key_df[-5, ], aes(x = x + 0.018, label = label), hjust = 0,
            size = 2.0, colour = COL["graphite"]) +
  annotate("point", x = 0.805, y = 1, shape = 21, size = 2.25,
           fill = unname(pal_dir["mr_bonf"]), colour = "#7A858A",
           stroke = 0.45) +
  annotate("text", x = 0.825, y = 1, label = "MR-PRESSO sig.",
           hjust = 0, size = 2.0, colour = COL["graphite"]) +
  annotate("point", x = 1.010, y = 1, shape = 24, size = 2.45,
           fill = "#F0F3F5", colour = "#7A858A", stroke = 0.50) +
  annotate("text", x = 1.030, y = 1, label = "MR-PRESSO not sig.",
           hjust = 0, size = 2.0, colour = COL["graphite"]) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1.24), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.88, 1.12), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "sans") +
  theme(plot.margin = margin(0, 18, 3, 4))

directional_body <- (p_labels | p_lcv | p_mr_same | p_mr_reverse) +
  plot_layout(widths = c(1.55, 1.85, 1.85, 1.85)) +
  plot_annotation(
    title = "Directional analyses remain secondary",
    subtitle = paste(
      "Only LCV Bonferroni-positive directions are shown.",
      "MR columns are IVW forest plots in the LCV-consistent and reverse directions;",
      "CI crossing zero marks non-significance, triangle = MR-PRESSO estimate not significant, dark outline = full strict robust."
    )
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 9.2,
                              colour = COL["graphite"]),
    plot.subtitle = element_text(size = 6.0, colour = COL["grey"])
  )

directional_preview <- directional_body / p_key +
  plot_layout(heights = c(1, 0.065))

lcv_source <- lcv_focus[
  , c("pair_id", "trait1", "trait2", "disease_pair_label",
      "lcv_direction", "LCV_gcp", "LCV_gcp_se", "LCV_ci_low",
      "LCV_ci_high", "LCV_gcp_p", "LCV_sig", "mr_direction_id",
      "mr_reverse_direction_id")
]

mr_source_cols <- c(
  "pair_id.y", "pair_id.x", "lcv_direction", "mr_direction_role",
  "disease_pair_label", "exposure", "outcome", "mr_label",
  "n_instruments_after_harmonise", "beta", "se", "beta_low",
  "beta_high", "pval", "IVW_status", "MR_PRESSO_p", "MR_PRESSO_sig",
  "strict_robust", "strict_label"
)
mr_source <- rbind(
  mr_same[, mr_source_cols],
  mr_reverse[, mr_source_cols]
)
names(mr_source)[names(mr_source) == "pair_id.x"] <- "lcv_pair_id"
names(mr_source)[names(mr_source) == "pair_id.y"] <- "mr_pair_id"

write_source(lcv_source, "Figure4C_directional_preview_LCV_source")
write_source(mr_source, "Figure4C_directional_preview_MR_source")
save_artifact(directional_preview, "Figure4C_directional_preview",
              ".", 282, 106, "Figure4C_directional_preview",
              "", "preview", dpi = 600)

cat("LCV Bonferroni directions:", nrow(lcv_source), "\n")
cat("MR rows including reverse:", nrow(mr_source), "\n")
cat("Same-direction strict robust:",
    sum(mr_source$mr_direction_role == "LCV direction" &
          mr_source$strict_robust), "\n")
cat("Reverse-direction strict robust:",
    sum(mr_source$mr_direction_role == "Reverse direction" &
          mr_source$strict_robust), "\n")
