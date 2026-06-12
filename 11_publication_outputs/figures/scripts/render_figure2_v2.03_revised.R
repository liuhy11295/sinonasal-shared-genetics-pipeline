#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(grid)
})

if (!nzchar(Sys.getenv("PAPER_ROOT"))) stop("Set PAPER_ROOT to the manuscript figure project.")
ROOT <- normalizePath(Sys.getenv("PAPER_ROOT"), mustWork = TRUE)
SRC <- file.path(ROOT, "figure_v2_source_data")
GENE_A_PATH <- file.path(
  Sys.getenv("RESULTS_ROOT"), "score_final", "gene_tables",
  "final_geneA_pair_gene_records.tsv"
)
OUT <- file.path(ROOT, "figure_final")
MAIN_COMP <- file.path(OUT, "main", "composites")
MAIN_PAN <- file.path(OUT, "main", "panels")
SUPP_COMP <- file.path(OUT, "supplementary", "composites")
SOURCE_OUT <- file.path(OUT, "source_data")
QA <- file.path(OUT, "qa")

dir.create(MAIN_COMP, recursive = TRUE, showWarnings = FALSE)
dir.create(MAIN_PAN, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_COMP, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(QA, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", quote = "", check.names = FALSE,
             comment.char = "", na.strings = c("", "NA", "NaN"))
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

short_trait <- function(x) {
  x <- sub("^ALLERGIC_RHINITIS_GCST90038664$", "Allergic rhinitis", x)
  x <- sub("^NASAL_POLYPS_GCST90018883$", "Nasal polyps", x)
  x <- gsub("_GCST[0-9]+", "", x)
  x <- gsub("^(H7|H8|J10|K11|L12|M13)_", "", x)
  x <- gsub("_EXMORE|_SUGG|_STRICT2|_STRICT|_WIDE", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

display_partner <- function(x) {
  x <- as.character(x)
  x[x == "Allerg Asthma"] <- "Allergic asthma"
  x[x == "Allergicconjunctivitis"] <- "Allergic conjunctivitis"
  x[x == "Dermatitiseczema"] <- "Dermatitis and eczema"
  x[x == "Nonallerg Asthma"] <- "Non-allergic asthma"
  x[x == "Asthma Eosinophil"] <- "Eosinophilic asthma"
  x[x == "Copd"] <- "COPD"
  x[x == "Atopic"] <- "Atopic dermatitis"
  x[x == "Autoimmune Nonthyroid"] <- "Non-thyroid autoimmune disease"
  x[x == "Chrontonsaden"] <- "Chronic tonsil and adenoid disease"
  x
}

wrap_text <- function(x, width = 30) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

trim_text <- function(x, n = 34) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "..."), x)
}

save_pdf_png <- function(plot, stem, width_in, height_in, dpi = 450) {
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "sans")
  print(plot)
  grDevices::dev.off()
  grDevices::png(paste0(stem, ".png"), width = width_in, height = height_in,
                 units = "in", res = dpi, type = "cairo")
  print(plot)
  grDevices::dev.off()
}

save_panel <- function(plot, stem, width_in, height_in, dpi = 600) {
  save_pdf_png(plot, file.path(MAIN_PAN, stem), width_in, height_in, dpi = dpi)
}

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
  rose = "#A65E8B",
  violet = "#7567A8",
  red = "#DC0000"
)

theme_f2 <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.28, colour = COL["graphite"]),
      axis.ticks = element_line(linewidth = 0.25, colour = COL["graphite"]),
      axis.text = element_text(colour = COL["graphite"]),
      axis.title = element_text(colour = COL["graphite"]),
      panel.grid.major.y = element_line(colour = "#EEF2F4", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, colour = COL["grey"], hjust = 0),
      plot.margin = margin(4, 7, 4, 7),
      legend.title = element_blank(),
      legend.key.height = unit(3.0, "mm"),
      strip.background = element_rect(fill = COL["pale"], colour = NA),
      strip.text = element_text(face = "bold", size = base_size)
    )
}

theme_set(theme_f2())

ldsc <- read_tsv(file.path(SRC, "ldsc.tsv"))
pair_ev <- read_tsv(file.path(SRC, "pair_evidence.tsv"))
locus_ab <- read_tsv(file.path(SRC, "locus_ab.tsv"))
magma_genes <- read_tsv(file.path(SRC, "magma_significant_genes.tsv"))
chr_offsets <- read_tsv(file.path(SRC, "chromosome_offsets.tsv"))
mtag_loci <- read_tsv(file.path(SRC, "mtag_top_loci.tsv"))
gene_a_records <- read_tsv(GENE_A_PATH)
cytoband_hg19 <- read.delim(
  gzfile(file.path(SOURCE_OUT, "cytoBand_hg19.txt.gz")),
  sep = "\t", header = FALSE, quote = "", comment.char = "",
  col.names = c("chrom", "band_start", "band_end", "band", "stain")
)

map_cytoband <- function(chr, pos) {
  vapply(seq_along(chr), function(i) {
    hit <- cytoband_hg19 %>%
      filter(
        chrom == paste0("chr", chr[i]),
        pos[i] >= band_start,
        pos[i] < band_end
      )
    if (nrow(hit)) paste0(chr[i], hit$band[1]) else paste0(chr[i], "q?")
  }, character(1))
}

chr_df <- chr_offsets %>%
  transmute(chr = as.integer(chr), chr_len = as.numeric(chr_len), offset = as.numeric(offset)) %>%
  arrange(chr) %>%
  mutate(mid = offset + chr_len / 2,
         odd = chr %% 2,
         genome_end = offset + chr_len)
genome_total <- max(chr_df$genome_end, na.rm = TRUE)

pair_summary <- pair_ev %>%
  left_join(ldsc %>% select(pair_id, rg = x, rg_se = size, rg_p = pvalue,
                            rg_q = ldsc_bh_q, rg_lo = lo, rg_hi = hi),
            by = "pair_id") %>%
  mutate(
    trait1 = ifelse(anchor == "Allergic rhinitis",
                    "ALLERGIC_RHINITIS_GCST90038664", "NASAL_POLYPS_GCST90018883"),
    trait2 = sub("^.*__", "", pair_id),
    pair_name = paste(anchor_short, display_partner(partner_short), sep = " \u2013 "),
    base_group = anchor,
    lava_sig_regions = as.integer(lava_n_significant_loci),
    rg = as.numeric(rg)
  ) %>%
  arrange(factor(base_group, levels = c("Allergic rhinitis", "Nasal polyps")),
          desc(rg), desc(lava_sig_regions), partner_short) %>%
  slice_head(n = 31) %>%
  mutate(pair_order = row_number(),
         P = paste0("P", pair_order),
         P = factor(P, levels = paste0("P", seq_len(n())))) %>%
  select(P, pair_order, pair_id, trait1, trait2, pair_name, base_group,
         partner_short, phenotype_domain, rg, rg_se, rg_lo, rg_hi, rg_p, rg_q,
         lava_n_loci, lava_sig_regions, magma_sig_gene_count, route)

pair_map <- pair_summary %>%
  select(P, pair_order, pair_id, pair_name, base_group, partner_short, phenotype_domain)

write_tsv(pair_summary, file.path(SOURCE_OUT, "pair_summary.tsv"))

gene_nearest <- magma_genes %>%
  filter(pair_id %in% pair_map$pair_id) %>%
  mutate(gene_pos = as.numeric(genome_pos),
         gene_start = as.numeric(start),
         gene_stop = as.numeric(stop),
         gene_score = as.numeric(score),
         gene_p = as.numeric(p)) %>%
  arrange(pair_id, chr, gene_p) %>%
  select(pair_id, gene_symbol, gene_chr = chr, gene_start, gene_stop,
         gene_pos, gene_score, gene_p) %>%
  distinct()

gene_a_by_locus <- gene_a_records %>%
  filter(gene_grade == "Gene-A", locus_grade == "A") %>%
  arrange(pair_id, locus_id, best_twas_fdr, magma_p, gene_symbol) %>%
  group_by(pair_id, locus_id) %>%
  summarise(
    gene_a_genes = paste(unique(gene_symbol), collapse = ", "),
    n_gene_a = n_distinct(gene_symbol),
    .groups = "drop"
  )

recurrent_loci <- locus_ab %>%
  filter(pair_id %in% pair_map$pair_id) %>%
  left_join(pair_map, by = "pair_id") %>%
  left_join(gene_a_by_locus, by = c("pair_id", "locus_id")) %>%
  mutate(
    chr = as.integer(chr),
    start = as.numeric(start),
    end = as.numeric(end),
    mid = as.numeric(mid),
    genome_pos = as.numeric(genome_pos),
    support_n = as.integer(coloc_positive %in% c(TRUE, "True")) +
      as.integer(susie_positive %in% c(TRUE, "True")) +
      as.integer(has_mtag_or_lava %in% c(TRUE, "True")) +
      ifelse(locus_grade == "A", 2L, 1L),
    support_n = pmax(1L, support_n),
    p_value = NA_real_,
    cytoband = paste0("chr", chr, ":", round(start / 1e6, 1), "-", round(end / 1e6, 1), "Mb"),
    cytoband_name = map_cytoband(chr, mid),
    recurrent_key = paste(chr, round(mid / 1e6), sep = ":")
  ) %>%
  rowwise() %>%
  mutate(
    positive_genes = {
      gd <- gene_nearest[gene_nearest$pair_id == pair_id &
                           gene_nearest$gene_chr == chr, , drop = FALSE]
      overlap <- gd[gd$gene_start <= end & gd$gene_stop >= start, , drop = FALSE]
      if (!nrow(overlap)) {
        NA_character_
      } else {
        overlap <- overlap[order(overlap$gene_p, -overlap$gene_score), , drop = FALSE]
        paste(head(unique(overlap$gene_symbol), 2), collapse = ", ")
      }
    },
    nearest_gene = {
    gd <- gene_nearest[gene_nearest$pair_id == pair_id &
                         gene_nearest$gene_chr == chr, , drop = FALSE]
    if (!nrow(gd)) {
      NA_character_
    } else {
      gd$dist <- abs(gd$gene_pos - genome_pos)
      gd$gene_symbol[which.min(gd$dist)]
    }
  }) %>%
  ungroup()

mtag_locus_support <- mtag_loci %>%
  filter(pair_id %in% pair_map$pair_id) %>%
  mutate(chr = as.integer(CHR),
         bp = as.numeric(BP),
         mtag_score = as.numeric(score),
         recurrent_key = paste(chr, round(bp / 1e6), sep = ":")) %>%
  group_by(pair_id, recurrent_key) %>%
  summarise(best_mtag_score = max(mtag_score, na.rm = TRUE),
            best_mtag_snp = SNP[which.max(mtag_score)][1],
            .groups = "drop")

recurrent_loci <- recurrent_loci %>%
  left_join(mtag_locus_support, by = c("pair_id", "recurrent_key")) %>%
  group_by(chr, recurrent_key) %>%
  mutate(n_pairs_at_locus = n_distinct(pair_id)) %>%
  ungroup() %>%
  mutate(
    support_n = support_n + ifelse(is.finite(best_mtag_score), 1L, 0L),
    priority = n_pairs_at_locus * 100 + support_n * 10 +
      ifelse(locus_grade == "A", 4, 0) +
      ifelse(is.finite(best_mtag_score), pmin(best_mtag_score, 40) / 10, 0),
    nearest_gene = ifelse(is.na(nearest_gene) | nearest_gene == "", lead_snp, nearest_gene),
    positive_genes = ifelse(
      is.na(positive_genes) | positive_genes == "",
      ifelse(is.na(nearest_gene) | nearest_gene == "", cytoband, nearest_gene),
      positive_genes
    ),
    label = paste0(cytoband, " / ", nearest_gene),
    plot_label = paste0(cytoband, "\n", nearest_gene)
  ) %>%
  select(P, pair_order, pair_id, pair_name, base_group, phenotype_domain,
         locus_id, chr, start, end, mid, genome_pos, cytoband, nearest_gene,
         cytoband_name,
         positive_genes, gene_a_genes, n_gene_a, label, plot_label,
         locus_grade, support_n,
         n_pairs_at_locus, best_mtag_score,
         best_mtag_snp, priority, coloc_pph4, susie_pph4_max, n_overlapping_core_snps)

locus_genes <- recurrent_loci %>%
  select(P, pair_id, locus_id, chr, mid, genome_pos, cytoband, nearest_gene,
         cytoband_name,
         positive_genes, gene_a_genes, n_gene_a, locus_grade, support_n,
         n_pairs_at_locus, best_mtag_score)

write_tsv(recurrent_loci %>% mutate(plot_label = gsub("\n", " / ", plot_label, fixed = TRUE)),
          file.path(SOURCE_OUT, "recurrent_loci.tsv"))
write_tsv(locus_genes, file.path(SOURCE_OUT, "locus_genes.tsv"))

pal_anchor <- c("Allergic rhinitis" = unname(COL["ar"]),
                "Nasal polyps" = unname(COL["np"]))

p_rg <- ggplot(pair_summary, aes(P, rg, fill = base_group)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = rg_lo, ymax = rg_hi), width = 0.18,
                linewidth = 0.22, colour = "black", alpha = 0.88) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = COL["graphite"]) +
  scale_fill_manual(values = pal_anchor) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.13))) +
  labs(title = "a  Genetic comorbidity overview across 31 positive disease pairs",
       subtitle = "Numbered pairs link genome-wide rg, membership, local rg regions and recurrent loci",
       x = NULL, y = "LDSC genetic\ncorrelation rg") +
  theme_f2(7.2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#F0F3F5", linewidth = 0.22),
        legend.position = c(0.82, 0.82),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", colour = NA))

membership <- pair_summary %>%
  transmute(P, pair_order, base_group, comparator = partner_short) %>%
  pivot_longer(c(base_group, comparator), names_to = "slot", values_to = "disease") %>%
  mutate(
    dot_group = ifelse(slot == "base_group", disease, "Comparator disease")
  )

row_levels <- unique(c("Allergic rhinitis", "Nasal polyps",
                       rev(unique(membership$disease[membership$slot == "comparator"]))))
membership$disease <- factor(membership$disease, levels = row_levels)

p_matrix <- ggplot(membership, aes(P, disease)) +
  geom_point(aes(colour = dot_group), size = 2.0) +
  scale_colour_manual(values = c(pal_anchor, "Comparator disease" = unname(COL["comparator"]))) +
  labs(x = NULL, y = NULL) +
  theme_f2(6.1) +
  theme(axis.text.x = element_text(size = 5.6, angle = 0),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_line(colour = "#F2F4F5", linewidth = 0.16))

p_lava <- ggplot(pair_summary, aes(P, lava_sig_regions, fill = base_group)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = pal_anchor, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = "Significant\nlocal rg\nregions") +
  theme_f2(7.0) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ideo_chr <- chr_df %>%
  mutate(xmin = offset / genome_total,
         xmax = genome_end / genome_total,
         xmid = mid / genome_total)

ideo_loci <- recurrent_loci %>%
  mutate(x_locus = genome_pos / genome_total,
         x_pair = (pair_order - 0.5) / nrow(pair_summary),
         y_locus = 0.30,
         y_pair = 0.98) %>%
  group_by(P) %>%
  arrange(desc(priority), .by_group = TRUE) %>%
  slice_head(n = 6) %>%
  ungroup()

ideo_points <- recurrent_loci %>%
  mutate(x_locus = genome_pos / genome_total) %>%
  group_by(chr, recurrent_key = paste(chr, round(mid / 1e6), sep = ":")) %>%
  summarise(x_locus = median(x_locus, na.rm = TRUE),
            n_pairs = n_distinct(P),
            max_priority = max(priority, na.rm = TRUE),
            label = paste0(cytoband[which.max(priority)], " / ",
                           nearest_gene[which.max(priority)]),
            .groups = "drop") %>%
  arrange(desc(n_pairs), desc(max_priority)) %>%
  slice_head(n = 65)

ideo_labels <- ideo_points %>%
  arrange(desc(n_pairs), desc(max_priority)) %>%
  slice_head(n = 10) %>%
  arrange(x_locus) %>%
  mutate(label_level = rep(c(1L, 3L, 2L), length.out = n()),
         y_label = c(-0.60, -0.80, -1.00)[label_level],
         short_label = sub(" / ", "\n", label, fixed = TRUE)) %>%
  group_by(label_level) %>%
  mutate(
    label_x = pmin(
      0.97,
      pmax(0.03, cummax(x_locus - row_number() * 0.145) + row_number() * 0.145)
    )
  ) %>%
  ungroup()

p_ideo <- ggplot() +
  geom_segment(data = ideo_chr,
               aes(x = xmin, xend = xmax, y = 0, yend = 0,
                   colour = factor(odd)),
               linewidth = 7, lineend = "butt") +
  geom_segment(data = ideo_loci,
               aes(x = x_locus, xend = x_pair, y = y_locus, yend = y_pair),
               linewidth = 0.18, colour = COL["line"], alpha = 0.34) +
  geom_point(data = ideo_points,
             aes(x = x_locus, y = 0.30, size = n_pairs),
             colour = COL["rose"], alpha = 0.88) +
  geom_segment(data = ideo_labels,
               aes(x = x_locus, xend = label_x, y = 0.23, yend = y_label + 0.06),
               linewidth = 0.20, colour = COL["line_dark"], alpha = 0.68) +
  geom_label(data = ideo_labels,
             aes(x = label_x, y = y_label, label = short_label),
             size = 1.58, linewidth = 0, label.padding = unit(0.65, "mm"),
             fill = "white", colour = COL["graphite"], alpha = 0.92) +
  geom_text(data = chr_df,
            aes(x = mid / genome_total, y = -0.23, label = chr),
            size = 1.85, colour = COL["graphite"]) +
  scale_colour_manual(values = c("0" = unname(COL["chr_light"]),
                                 "1" = unname(COL["chr_dark"])),
                      guide = "none") +
  scale_size_continuous(range = c(1.0, 4.2), guide = "none") +
  coord_cartesian(ylim = c(-1.12, 1.08), clip = "off") +
  labs(title = "Recurrent loci ideogram", x = NULL, y = NULL) +
  theme_void(base_size = 7, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 7.6, hjust = 0),
        plot.margin = margin(0, 7, 2, 7))

fig2a <- p_rg / p_matrix / p_lava / p_ideo +
  plot_layout(heights = c(2.15, 1.95, 1.10, 1.72))

display_gene_count <- recurrent_loci %>%
  filter(!is.na(nearest_gene), nearest_gene != "") %>%
  group_by(pair_id) %>%
  summarise(
    n_display_genes = n_distinct(nearest_gene),
    recurrent_hits = sum(n_pairs_at_locus > 1, na.rm = TRUE),
    .groups = "drop"
  )

pair_summary2 <- pair_summary %>%
  left_join(display_gene_count, by = "pair_id") %>%
  mutate(
    n_display_genes = ifelse(is.na(n_display_genes), 0L, as.integer(n_display_genes)),
    recurrent_hits = ifelse(is.na(recurrent_hits), 0L, as.integer(recurrent_hits))
  )

np_top3 <- pair_summary2 %>%
  filter(
    base_group == "Nasal polyps",
    n_display_genes > 0,
    !grepl("asthma", partner_short, ignore.case = TRUE) |
      partner_short == "Allerg Asthma"
  ) %>%
  arrange(desc(n_display_genes), desc(recurrent_hits), pair_order) %>%
  slice_head(n = 3)

ar_top3 <- pair_summary2 %>%
  filter(
    base_group == "Allergic rhinitis",
    n_display_genes > 0,
    !grepl("asthma", partner_short, ignore.case = TRUE) |
      partner_short == "Allerg Asthma"
  ) %>%
  arrange(desc(n_display_genes), desc(recurrent_hits), pair_order) %>%
  slice_head(n = 3)

selected_pairs <- c(as.character(np_top3$P), as.character(ar_top3$P))
selected_pairs_layout <- as.vector(rbind(as.character(np_top3$P), as.character(ar_top3$P)))

write_tsv(
  bind_rows(
    np_top3 %>% mutate(display_group = "Nasal polyps", display_rank = row_number()),
    ar_top3 %>% mutate(display_group = "Allergic rhinitis", display_rank = row_number())
  ) %>%
    select(display_group, display_rank, P, pair_order, pair_id, pair_name,
           base_group, n_display_genes, recurrent_hits),
  file.path(SOURCE_OUT, "Figure2B_selected_pairs.tsv")
)

adjust_label_y_positions <- function(df, min_gap = 0.12, y_min = -1.28, y_max = 1.28) {
  if (nrow(df) == 0) return(df)

  df <- df[order(df$label_y), , drop = FALSE]

  for (i in seq_len(nrow(df))) {
    if (i > 1) {
      if (df$label_y[i] - df$label_y[i - 1] < min_gap) {
        df$label_y[i] <- df$label_y[i - 1] + min_gap
      }
    }
  }

  overflow_top <- max(df$label_y, na.rm = TRUE) - y_max
  if (is.finite(overflow_top) && overflow_top > 0) {
    df$label_y <- df$label_y - overflow_top
  }

  overflow_bottom <- y_min - min(df$label_y, na.rm = TRUE)
  if (is.finite(overflow_bottom) && overflow_bottom > 0) {
    df$label_y <- df$label_y + overflow_bottom
  }

  df
}

make_circle_data <- function(pairs) {
  panels <- pair_map %>%
    filter(P %in% pairs) %>%
    mutate(P = factor(P, levels = pairs))

  chr_arc <- tidyr::crossing(panels, chr_df) %>%
    mutate(a0 = 2 * pi * offset / genome_total,
           a1 = 2 * pi * genome_end / genome_total,
           amid = 2 * pi * mid / genome_total) %>%
    rowwise() %>%
    mutate(theta = list(seq(a0, a1, length.out = 22))) %>%
    ungroup() %>%
    unnest(theta) %>%
    mutate(x = sin(theta), y = cos(theta))

  chr_lab <- tidyr::crossing(panels, chr_df) %>%
    mutate(theta = 2 * pi * mid / genome_total,
           x = 1.145 * sin(theta),
           y = 1.145 * cos(theta))

  loci_points <- recurrent_loci %>%
    filter(
      P %in% pairs,
      !is.na(nearest_gene),
      nearest_gene != ""
    ) %>%
    mutate(P = factor(P, levels = pairs),
           theta = 2 * pi * genome_pos / genome_total) %>%
    group_by(P) %>%
    arrange(desc(priority), .by_group = TRUE) %>%
    slice_head(n = 14) %>%
    ungroup() %>%
    mutate(x = 1.025 * sin(theta),
           y = 1.025 * cos(theta),
           grade_col = ifelse(locus_grade == "A", "A", "B"))

  loci_labels <- loci_points %>%
    mutate(
      short_label = paste0(nearest_gene, "\n", cytoband_name),
      label_lines = lengths(strsplit(short_label, "\n", fixed = TRUE))
    ) %>%
    group_by(P) %>%
    arrange(desc(priority), .by_group = TRUE) %>%
    distinct(nearest_gene, .keep_all = TRUE) %>%
    slice_head(n = 12) %>%
    ungroup() %>%
    mutate(
      side = ifelse(sin(theta) < 0, "left", "right"),
      label_x = ifelse(side == "left", -1.54, 1.54),
      label_y = 1.30 * cos(theta),
      mid_x = ifelse(side == "left", -1.20, 1.20),
      mid_y = label_y,
      hjust = ifelse(side == "left", 1, 0)
    ) %>%
    group_by(P, side) %>%
    group_modify(~ adjust_label_y_positions(
      .x,
      min_gap = max(0.12, 0.11 * max(.x$label_lines, na.rm = TRUE))
    )) %>%
    ungroup() %>%
    mutate(mid_y = label_y)

  list(
    panels = panels,
    chr_arc = chr_arc,
    chr_lab = chr_lab,
    loci_points = loci_points,
    loci_labels = loci_labels
  )
}

circle_plot <- function(pairs, title_text, ncol = 4, label_size = 1.55, center_size = 2.15) {
  cd <- make_circle_data(pairs)
  ggplot() +
    geom_path(data = cd$chr_arc,
              aes(x, y, group = interaction(P, chr), colour = factor(odd)),
              linewidth = 1.75, lineend = "butt") +
    geom_text(data = cd$chr_lab, aes(x, y, label = chr),
              size = 1.12, colour = COL["graphite"]) +
    geom_point(data = cd$loci_points,
               aes(x, y, size = support_n, fill = grade_col),
               shape = 21, stroke = 0.25, colour = COL["graphite"], alpha = 0.95) +
    geom_segment(
      data = cd$loci_labels,
      aes(x = x, y = y, xend = mid_x, yend = label_y),
      colour = COL["line_dark"], linewidth = 0.16, alpha = 0.58
    ) +
    geom_segment(
      data = cd$loci_labels,
      aes(x = mid_x, y = label_y, xend = label_x, yend = label_y),
      colour = COL["line_dark"], linewidth = 0.16, alpha = 0.58
    ) +
    geom_text(
      data = cd$loci_labels,
      aes(x = label_x, y = label_y, label = short_label, hjust = hjust),
      size = label_size, colour = COL["graphite"], lineheight = 0.82
    ) +
    geom_text(data = cd$panels,
              aes(0, 0, label = wrap_text(pair_name, 16), colour = base_group),
              size = center_size, fontface = "bold", lineheight = 0.82) +
    facet_wrap(~P, ncol = ncol) +
    scale_colour_manual(values = c("0" = unname(COL["chr_light"]),
                                   "1" = unname(COL["chr_dark"]),
                                   pal_anchor),
                        guide = "none") +
    scale_fill_manual(
      values = c("A" = unname(COL["np"]), "B" = unname(COL["amber"])),
      breaks = c("A", "B"),
      labels = c("Locus-A", "Locus-B"),
      guide = guide_legend(
        override.aes = list(size = 3.7, shape = 21,
                            colour = COL["graphite"], stroke = 0.30,
                            alpha = 0.95)
      )
    ) +
    scale_size_continuous(range = c(1.25, 3.25), guide = "none") +
    coord_fixed(xlim = c(-2.28, 2.28), ylim = c(-1.72, 1.72), clip = "off") +
    labs(title = title_text, fill = NULL) +
    theme_void(base_family = "sans", base_size = 7.2) +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0),
          strip.text = element_text(face = "bold", size = 7.5, margin = margin(1, 0, 1, 0)),
          legend.position = "bottom",
          legend.text = element_text(size = 7.0, colour = COL["graphite"]),
          legend.key.height = unit(4.2, "mm"),
          legend.key.width = unit(5.0, "mm"),
          plot.margin = margin(4, 4, 4, 4))
}

fig2b <- circle_plot(selected_pairs_layout,
                     "b  Representative disease-pair circular locus maps",
                     ncol = 2, label_size = 1.78, center_size = 1.95)

fig2 <- fig2a | fig2b
fig2 <- fig2 + plot_layout(widths = c(1.28, 1.05))

save_panel(fig2a, "Figure2A_revised", width_in = 8.9, height_in = 10.6, dpi = 600)
save_panel(fig2b, "Figure2B_revised", width_in = 7.5, height_in = 10.6, dpi = 600)
save_pdf_png(fig2, file.path(MAIN_COMP, "Figure2_revised"),
             width_in = 16.8, height_in = 10.6, dpi = 600)

main_circle_data <- make_circle_data(selected_pairs_layout)
write_tsv(
  main_circle_data$loci_labels %>%
    mutate(short_label = gsub("\n", " / ", short_label, fixed = TRUE)) %>%
    select(P, pair_id, locus_id, priority, side, x, y, mid_x, mid_y,
           label_x, label_y, hjust, short_label),
  file.path(SOURCE_OUT, "Figure2B_label_positions.tsv")
)
write_tsv(
  main_circle_data$loci_points %>%
    select(P, pair_id, locus_id, locus_grade, cytoband, nearest_gene,
           cytoband_name,
           n_pairs_at_locus, priority),
  file.path(SOURCE_OUT, "Figure2B_locus_points.tsv")
)

all_pairs <- as.character(pair_summary$P)
supp_pages <- split(all_pairs, ceiling(seq_along(all_pairs) / 8))
grDevices::cairo_pdf(file.path(SUPP_COMP, "Supplementary_all_pair_circular_loci.pdf"),
                     width = 12.6, height = 9.4, family = "sans", onefile = TRUE)
for (i in seq_along(supp_pages)) {
  print(circle_plot(supp_pages[[i]],
                    paste0("Supplementary circular chromosome locus plots, page ", i),
                    ncol = 4, label_size = 1.18, center_size = 1.72))
}
grDevices::dev.off()

all_circle_labels <- make_circle_data(all_pairs)$loci_labels
label_audit <- all_circle_labels %>%
  arrange(P, side, label_y) %>%
  group_by(P, side) %>%
  summarise(
    n_labels = n(),
    min_vertical_gap = ifelse(n() > 1, min(diff(label_y)), NA_real_),
    .groups = "drop"
  )
write_tsv(label_audit, file.path(QA, "circular_label_audit.tsv"))

invisible(file.copy(file.path(MAIN_COMP, "Figure2_revised.pdf"),
                    file.path(OUT, "Figure2_revised.pdf"), overwrite = TRUE))
invisible(file.copy(file.path(MAIN_COMP, "Figure2_revised.png"),
                    file.path(OUT, "Figure2_revised.png"), overwrite = TRUE))
invisible(file.copy(file.path(MAIN_PAN, "Figure2A_revised.pdf"),
                    file.path(OUT, "Figure2A_revised.pdf"), overwrite = TRUE))
invisible(file.copy(file.path(MAIN_PAN, "Figure2A_revised.png"),
                    file.path(OUT, "Figure2A_revised.png"), overwrite = TRUE))
invisible(file.copy(file.path(MAIN_PAN, "Figure2B_revised.pdf"),
                    file.path(OUT, "Figure2B_revised.pdf"), overwrite = TRUE))
invisible(file.copy(file.path(MAIN_PAN, "Figure2B_revised.png"),
                    file.path(OUT, "Figure2B_revised.png"), overwrite = TRUE))
invisible(file.copy(file.path(SUPP_COMP, "Supplementary_all_pair_circular_loci.pdf"),
                    file.path(OUT, "Supplementary_all_pair_circular_loci.pdf"),
                    overwrite = TRUE))

qa <- data.frame(
  artifact = c("Figure2_revised.pdf", "Figure2_revised.png",
               "Figure2A_revised.pdf", "Figure2A_revised.png",
               "Figure2B_revised.pdf", "Figure2B_revised.png",
               "Supplementary_all_pair_circular_loci.pdf"),
  path = c(file.path(MAIN_COMP, "Figure2_revised.pdf"),
           file.path(MAIN_COMP, "Figure2_revised.png"),
           file.path(MAIN_PAN, "Figure2A_revised.pdf"),
           file.path(MAIN_PAN, "Figure2A_revised.png"),
           file.path(MAIN_PAN, "Figure2B_revised.pdf"),
           file.path(MAIN_PAN, "Figure2B_revised.png"),
           file.path(SUPP_COMP, "Supplementary_all_pair_circular_loci.pdf")),
  stringsAsFactors = FALSE
) %>%
  mutate(file_size = file.info(path)$size,
         exists = file.exists(path),
         automated_status = ifelse(exists & file_size > 15000, "pass", "review"))
write_tsv(qa, file.path(QA, "automated_render_QA.tsv"))

session <- capture.output(sessionInfo())
writeLines(session, file.path(QA, "R_sessionInfo.txt"))

message("Wrote revised Figure 2 outputs to: ", OUT)
