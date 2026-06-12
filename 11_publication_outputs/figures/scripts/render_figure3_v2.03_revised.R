#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(grid)
})

if (!nzchar(Sys.getenv("PAPER_ROOT"))) stop("Set PAPER_ROOT to the manuscript figure project.")
ROOT <- normalizePath(Sys.getenv("PAPER_ROOT"), mustWork = TRUE)
SRC <- file.path(ROOT, "figure_v2_source_data")
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
  x_out <- x
  text_cols <- vapply(x_out, function(z) is.character(z) || is.factor(z), logical(1))
  for (nm in names(x_out)[text_cols]) {
    x_out[[nm]] <- gsub("[\r\n]+", " | ", as.character(x_out[[nm]]))
  }
  write.table(x_out, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

wrap_text <- function(x, width = 28) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

clean_term <- function(x) {
  x <- gsub("^Homo sapiens: ", "", x)
  x <- gsub(" pathway$", "", x, ignore.case = TRUE)
  x <- gsub("signaling pathway", "signalling", x, ignore.case = TRUE)
  x
}

short_term <- function(x, width = 26) {
  x <- clean_term(x)
  x <- gsub("cell surface toll-like receptor", "cell-surface TLR", x, ignore.case = TRUE)
  x <- gsub("toll-like receptor", "TLR", x, ignore.case = TRUE)
  x <- gsub("pattern recognition receptor", "PRR", x, ignore.case = TRUE)
  x <- gsub("cellular response to", "response to", x, ignore.case = TRUE)
  x <- gsub("molecule of bacterial origin", "bacterial molecule", x, ignore.case = TRUE)
  x <- gsub("Diseases associated with the TLR signaling cascade", "TLR disease cascade",
            x, ignore.case = TRUE)
  x <- gsub("Diseases of Immune System", "immune-system disease", x, ignore.case = TRUE)
  wrap_text(x, width)
}

clean_tissue_term <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("cervical c-1", "cervical C1", x, ignore.case = TRUE)
  x
}

format_fuma_group <- function(x) {
  dplyr::case_when(
    x == "global_primary" ~ "Global primary",
    x == "global_strict" ~ "Global strict",
    x == "asthma_lower_airway" ~ "Asthma/lower airway",
    x == "atopic_allergic" ~ "Atopic/allergic",
    x == "autoimmune_IBD" ~ "Autoimmune / IBD",
    x == "ENT_infection" ~ "ENT / infection",
    TRUE ~ x
  )
}

format_fuma_category <- function(x) {
  dplyr::case_when(
    x == "DEG.up" ~ "Up",
    x == "DEG.down" ~ "Down",
    x == "DEG.twoside" ~ "Two-sided",
    TRUE ~ x
  )
}

save_pdf_png <- function(plot, stem, width_in, height_in, dpi = 600) {
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "sans")
  print(plot)
  grDevices::dev.off()
  grDevices::png(paste0(stem, ".png"), width = width_in, height = height_in,
                 units = "in", res = dpi, type = "cairo")
  print(plot)
  grDevices::dev.off()
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
  red = "#DC0000",
  absent = "#F6F8F9"
)

pal_group <- c(
  "Global" = unname(COL["comparator"]),
  "Asthma/lower airway" = unname(COL["ar"]),
  "Atopic/allergic" = unname(COL["np"]),
  "Other domains" = unname(COL["line"])
)
pal_gene <- c("Gene-A" = unname(COL["teal"]),
              "Gene-B" = unname(COL["amber"]))
pal_direction <- c("Up" = unname(COL["np"]),
                   "Down" = unname(COL["ar"]),
                   "Two-sided" = unname(COL["comparator"]))
pal_database <- c("GO_BP" = unname(COL["ar"]),
                  "KEGG" = unname(COL["amber"]),
                  "Reactome" = unname(COL["rose"]))

theme_f3 <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.28, colour = COL["graphite"]),
      axis.ticks = element_line(linewidth = 0.25, colour = COL["graphite"]),
      axis.text = element_text(colour = COL["graphite"]),
      axis.title = element_text(colour = COL["graphite"]),
      panel.grid.major.x = element_line(colour = "#EEF2F4", linewidth = 0.22),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1.4, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.4, colour = COL["grey"], hjust = 0),
      legend.title = element_text(face = "bold", size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = unit(3.0, "mm"),
      strip.background = element_rect(fill = COL["pale"], colour = NA),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(5, 7, 5, 7)
    )
}

theme_set(theme_f3())

pathways <- read_tsv(file.path(SRC, "pathways.tsv"))
pathway_counts <- read_tsv(file.path(
  ROOT, "result", "Table", "Table_R5_pathway_term_count_summary.tsv"
))
twas <- read_tsv(file.path(SRC, "twas_matrix.tsv"))
cell <- read_tsv(file.path(SRC, "cell_audit.tsv"))
RESULTS_END <- normalizePath(file.path(dirname(ROOT), "results_end"), mustWork = TRUE)
TWAS_PAIR_RECORDS <- file.path(
  RESULTS_END,
  "final_gtexv8_twas_total_domain_fuma_enrichment",
  "final_gene_tables",
  "final_twas_supported_pair_gene_records.tsv"
)
twas_pair_records <- read_tsv(TWAS_PAIR_RECORDS)
pair_evidence <- read_tsv(file.path(SRC, "pair_evidence.tsv"))
FUMA_TISSUE <- file.path(RESULTS_END, "integrated_fuma_tissue_cell_ora",
                         "fuma_gene2func_tissue_enrichment")
FUMA_PAIRWISE <- file.path(RESULTS_END, "integrated_fuma_tissue_cell_ora",
                           "fuma_pairwise_tissue_annotation_existing")
fuma_tissue54_fdr <- read_tsv(file.path(
  FUMA_TISSUE, "fuma_gtex_v8_54_tissue_DEG_enrichment_FDR05.tsv"
))
fuma_tissue30_fdr <- read_tsv(file.path(
  FUMA_TISSUE, "fuma_gtex_v8_30_general_tissue_DEG_enrichment_FDR05.tsv"
))
pairwise_fuma_summary <- read_tsv(file.path(
  FUMA_PAIRWISE, "pairwise_fuma_tissue_FDR05_summary.tsv"
))
pairwise_fuma54_top <- read_tsv(file.path(
  FUMA_PAIRWISE, "pairwise_fuma_gtex54_tissue_top_terms.tsv"
))

pathways <- pathways %>%
  filter(!is.na(group), group != "", !is.na(FDR), is.finite(FDR), FDR > 0,
         !is.na(term_name), term_name != "") %>%
  mutate(
    group = factor(group, levels = c("Global", "Asthma/lower airway",
                                     "Atopic/allergic", "Other domains")),
    database = factor(database, levels = c("GO_BP", "KEGG", "Reactome")),
    term_clean = clean_term(term_name),
    term_wrapped = wrap_text(term_clean, 32),
    term_short = short_term(term_clean, 25),
    neglog10_fdr = -log10(FDR),
    overlap_count = as.numeric(overlap_count),
    odds_ratio = as.numeric(odds_ratio)
  ) %>%
  filter(!is.na(group))

panel_a_data <- pathways %>%
  filter(database == "GO_BP", group != "Other domains") %>%
  arrange(FDR, desc(overlap_count)) %>%
  distinct(term_clean, .keep_all = TRUE) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  arrange(FDR) %>%
  mutate(term_short = factor(term_short, levels = rev(unique(term_short))))

p3a <- ggplot(panel_a_data,
              aes(neglog10_fdr, term_short, size = overlap_count,
                  fill = group)) +
  geom_point(shape = 21, colour = "white", alpha = 0.94, stroke = 0.34) +
  geom_vline(xintercept = -log10(0.05), linewidth = 0.25,
             linetype = "dashed", colour = COL["grey"]) +
  scale_fill_manual(values = pal_group) +
  scale_size_continuous(range = c(1.6, 5.4)) +
  labs(title = "a  Top enriched biological-process pathways",
       subtitle = "GO BP terms ranked by FDR across analysis groups",
       x = expression(-log[10](FDR)), y = NULL,
       fill = "Analysis group", size = "Overlap genes") +
  theme_f3(6.8) +
  theme(legend.position = "bottom")

top_terms <- panel_a_data %>%
  arrange(FDR, desc(overlap_count)) %>%
  mutate(term_heatmap = short_term(term_clean, 18))

term_gene <- do.call(rbind, lapply(seq_len(nrow(top_terms)), function(i) {
  genes <- unlist(strsplit(top_terms$overlap_genes[i], ";", fixed = TRUE))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) return(NULL)
  data.frame(term_clean = top_terms$term_clean[i],
             term_short = as.character(top_terms$term_short[i]),
             term_heatmap = top_terms$term_heatmap[i],
             group = as.character(top_terms$group[i]),
             neglog10_fdr = top_terms$neglog10_fdr[i],
             gene = genes,
             row.names = NULL,
             stringsAsFactors = FALSE)
}))

gene_summary <- term_gene %>%
  count(gene, name = "recurrence_count") %>%
  arrange(desc(recurrence_count), gene) %>%
  slice_head(n = 20)

top_genes <- gene_summary %>%
  pull(gene)

term_lookup <- top_terms %>%
  transmute(term_heatmap = as.character(term_heatmap),
            group = as.character(group)) %>%
  distinct()

panel_b_hits <- term_gene %>%
  filter(gene %in% top_genes) %>%
  group_by(term_heatmap, gene) %>%
  summarise(neglog10_fdr = max(neglog10_fdr, na.rm = TRUE), .groups = "drop")

panel_b_data <- tidyr::crossing(
  term_heatmap = unique(as.character(top_terms$term_heatmap)),
  gene = top_genes
) %>%
  left_join(term_lookup, by = "term_heatmap") %>%
  left_join(gene_summary, by = "gene") %>%
  left_join(panel_b_hits, by = c("term_heatmap", "gene")) %>%
  mutate(
    present = !is.na(neglog10_fdr),
    fill_group = ifelse(present, group, "Absent"),
    neglog10_fdr = ifelse(is.na(neglog10_fdr), 0, neglog10_fdr),
    gene = factor(gene, levels = rev(top_genes)),
    term_heatmap = factor(term_heatmap, levels = unique(top_terms$term_heatmap))
  )

p3b <- ggplot(panel_b_data, aes(term_heatmap, gene, fill = fill_group)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  scale_fill_manual(
    values = c("Absent" = unname(COL["absent"]), pal_group),
    breaks = c("Global", "Asthma/lower airway", "Atopic/allergic")
  ) +
  labs(title = "b  Driver genes recurring across top pathways",
       subtitle = "Top recurrent genes mapped to the panel-a pathway set",
       x = NULL, y = NULL, fill = "Pathway group") +
  theme_f3(5.5) +
  theme(axis.text.x = element_text(angle = 58, hjust = 1, vjust = 1),
        legend.position = "bottom",
        panel.grid = element_blank())

tissue_support_all <- twas %>%
  filter(!is.na(tissue), tissue != "") %>%
  count(tissue, grade, name = "n_gene_records") %>%
  group_by(tissue) %>%
  summarise(total = sum(n_gene_records),
            gene_a = sum(n_gene_records[grade == "Gene-A"], na.rm = TRUE),
            gene_b = sum(n_gene_records[grade == "Gene-B"], na.rm = TRUE),
            .groups = "drop") %>%
  mutate(tissue_clean = clean_tissue_term(tissue))

pair_group_map <- pair_evidence %>%
  transmute(pair_id,
            base_group = anchor,
            partner,
            phenotype_domain_pair = phenotype_domain)

twas_support_by_group <- twas_pair_records %>%
  filter(!is.na(significant_tissue_list), significant_tissue_list != "") %>%
  separate_rows(significant_tissue_list, sep = ",") %>%
  mutate(
    tissue = trimws(significant_tissue_list),
    tissue_clean = clean_tissue_term(tissue),
    gene_grade = factor(gene_grade, levels = c("Gene-B", "Gene-A")),
    strict_positive = strict_positive %in% c(TRUE, "True", "TRUE", "true", "1")
  ) %>%
  left_join(pair_group_map, by = "pair_id") %>%
  mutate(
    phenotype_domain = dplyr::coalesce(phenotype_domain, phenotype_domain_pair),
    base_group = factor(base_group, levels = c("Allergic rhinitis", "Nasal polyps"))
  ) %>%
  filter(!is.na(base_group))

panel_c_all <- twas_support_by_group %>%
  group_by(base_group, tissue, tissue_clean, gene_grade) %>%
  summarise(
    n_pair_gene_tissue_records = n(),
    n_unique_genes = n_distinct(gene_symbol),
    n_pairs = n_distinct(pair_id),
    .groups = "drop"
  )

panel_c_tissue_totals <- twas_support_by_group %>%
  group_by(base_group, tissue, tissue_clean) %>%
  summarise(
    total_records = n(),
    total_unique_genes = n_distinct(gene_symbol),
    n_pairs = n_distinct(pair_id),
    gene_a_records = sum(gene_grade == "Gene-A", na.rm = TRUE),
    gene_b_records = sum(gene_grade == "Gene-B", na.rm = TRUE),
    gene_a_genes = n_distinct(gene_symbol[gene_grade == "Gene-A"]),
    gene_b_genes = n_distinct(gene_symbol[gene_grade == "Gene-B"]),
    .groups = "drop"
  )

top_tissues_by_group <- panel_c_tissue_totals %>%
  group_by(base_group) %>%
  arrange(desc(total_records), desc(total_unique_genes), tissue_clean, .by_group = TRUE) %>%
  slice_head(n = 12) %>%
  ungroup() %>%
  mutate(tissue_group = paste(as.character(base_group), tissue_clean, sep = "___"))

tissue_group_levels <- top_tissues_by_group %>%
  arrange(base_group, total_records, total_unique_genes, tissue_clean) %>%
  pull(tissue_group)

panel_c_data <- panel_c_all %>%
  inner_join(
    top_tissues_by_group %>%
      select(base_group, tissue, tissue_clean, total_records, total_unique_genes,
             n_pairs, gene_a_records, gene_b_records, gene_a_genes,
             gene_b_genes, tissue_group),
    by = c("base_group", "tissue", "tissue_clean")
  ) %>%
  mutate(
    tissue_group = factor(tissue_group, levels = tissue_group_levels),
    base_group = factor(base_group, levels = c("Allergic rhinitis", "Nasal polyps"))
  )

panel_c_labels <- top_tissues_by_group %>%
  mutate(
    tissue_group = factor(tissue_group, levels = tissue_group_levels),
    base_group = factor(base_group, levels = c("Allergic rhinitis", "Nasal polyps"))
  )

p3c <- ggplot(panel_c_data,
              aes(n_pair_gene_tissue_records, tissue_group,
                  fill = gene_grade)) +
  geom_col(width = 0.70, alpha = 0.92) +
  geom_text(
    data = panel_c_labels,
    aes(x = total_records, y = tissue_group, label = total_records),
    inherit.aes = FALSE,
    hjust = -0.14,
    size = 2.0,
    colour = COL["graphite"]
  ) +
  facet_wrap(~ base_group, ncol = 2, scales = "free_y") +
  scale_y_discrete(labels = function(z) wrap_text(sub("^.*___", "", z), 23)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.18))) +
  scale_fill_manual(values = pal_gene,
                    breaks = c("Gene-A", "Gene-B"),
                    drop = FALSE) +
  labs(title = "c  TWAS tissue support by disease group",
       subtitle = "Top GTEx tissues supported by Gene-A and Gene-B signals across disease pairs",
       x = "TWAS-supported pair-gene-tissue records", y = NULL,
       fill = "Gene grade") +
  theme_f3(6.2) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 6.8, face = "bold"),
        axis.text.y = element_text(size = 5.6),
        panel.spacing = unit(5, "mm"))

tissue_support <- tissue_support_all %>%
  arrange(desc(total)) %>%
  slice_head(n = 14) %>%
  mutate(tissue_label = tissue_clean,
         tissue_label = wrap_text(tissue_label, 22),
         tissue_label = factor(tissue_label, levels = rev(tissue_label)))

cell_marker <- cell %>%
  filter(run_ORA %in% c(TRUE, "True")) %>%
  group_by(marker_source) %>%
  summarise(n_tested_terms = sum(as.numeric(n_tested_terms), na.rm = TRUE),
            n_FDR05 = sum(as.numeric(n_FDR05), na.rm = TRUE),
            n_nominal_p05 = sum(as.numeric(n_nominal_p05), na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n_nominal_p05)) %>%
  slice_head(n = 8) %>%
  mutate(marker_source = factor(marker_source, levels = rev(marker_source)))

p3d1 <- ggplot(tissue_support, aes(total, tissue_label)) +
  geom_col(fill = COL["ar"], width = 0.70, alpha = 0.90) +
  geom_point(aes(x = gene_a), colour = COL["teal"], size = 1.5) +
  geom_point(aes(x = gene_b), colour = COL["amber"], size = 1.5) +
  labs(title = "a  Tissue support and cell-marker boundary",
       subtitle = "Bars: TWAS tissue records; dots: Gene-A/Gene-B components",
       x = "TWAS gene-tissue records", y = NULL) +
  theme_f3(5.8)

p3d2 <- ggplot(cell_marker, aes(n_nominal_p05, marker_source)) +
  geom_col(fill = COL["pale"], colour = "white", width = 0.70) +
  geom_point(aes(x = n_FDR05), colour = COL["rose"], size = 1.8) +
  geom_text(aes(label = paste0("FDR=", n_FDR05)), hjust = -0.12,
            size = 2.0, colour = COL["graphite"]) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = "Nominal cell-marker terms", y = NULL) +
  theme_f3(5.8)

p3d <- p3d1 / p3d2 + plot_layout(heights = c(1.35, 1.0))

supp_fuma_tissue <- fuma_tissue30_fdr %>%
  mutate(dataset = "GTEx v8 30 general tissue") %>%
  mutate(
    adjP = suppressWarnings(as.numeric(adjP)),
    p = suppressWarnings(as.numeric(p)),
    N_overlap = suppressWarnings(as.numeric(N_overlap)),
    analysis_group = format_fuma_group(analysis_id),
    direction = format_fuma_category(Category),
    direction = factor(direction, levels = c("Up", "Down", "Two-sided")),
    tissue_clean = clean_tissue_term(GeneSet),
    tissue_direction = paste0(tissue_clean, " (", tolower(as.character(direction)), ")"),
    tissue_label = wrap_text(tissue_direction, 28),
    neglog10_adjP = -log10(adjP),
    facet_label = paste(dataset, analysis_group, sep = "\n")
  ) %>%
  filter(analysis_level == "global") %>%
  arrange(dataset, analysis_group, adjP) %>%
  mutate(tissue_label = factor(tissue_label, levels = rev(unique(tissue_label))))

p3s1 <- ggplot(supp_fuma_tissue,
               aes(neglog10_adjP, tissue_label,
                   size = N_overlap, colour = direction)) +
  geom_point(alpha = 0.88, stroke = 0.45) +
  geom_vline(xintercept = -log10(0.05), linewidth = 0.25,
             linetype = "dashed", colour = COL["grey"]) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = pal_direction,
                      drop = FALSE) +
  scale_size_continuous(range = c(1.7, 4.8)) +
  labs(title = "b  GTEx30 general tissue sensitivity",
       subtitle = "FUMA GENE2FUNC FDR-significant general-tissue terms retained outside the main panel",
       x = expression(-log[10](FDR)), y = NULL,
       colour = "GTEx DEG set", size = "Overlap genes") +
  theme_f3(5.8) +
  theme(legend.position = "bottom")

supp_pathways <- pathways %>%
  filter(database %in% c("KEGG", "Reactome")) %>%
  arrange(FDR, desc(overlap_count)) %>%
  distinct(database, term_clean, .keep_all = TRUE) %>%
  group_by(database) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  arrange(database, FDR) %>%
  mutate(
    term_short = short_term(term_clean, 28),
    term_short = factor(term_short, levels = rev(unique(term_short)))
  )

p3s2 <- ggplot(supp_pathways,
               aes(neglog10_fdr, term_short, size = overlap_count, colour = group)) +
  geom_point(alpha = 0.88, stroke = 0.45) +
  geom_vline(xintercept = -log10(0.05), linewidth = 0.25,
             linetype = "dashed", colour = COL["grey"]) +
  facet_wrap(~ database, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = pal_group) +
  scale_size_continuous(range = c(1.4, 4.6)) +
  labs(title = "c  KEGG and Reactome support",
       subtitle = "Top compact pathway terms retained outside the main figure",
       x = expression(-log[10](FDR)), y = NULL,
       colour = "Analysis group", size = "Overlap genes") +
  theme_f3(5.8) +
  theme(legend.position = "bottom")

pathway_count_audit <- pathway_counts %>%
  filter(analysis_type %in% c("primary", "primary_domain"),
         database %in% c("GO_BP", "KEGG", "Reactome")) %>%
  mutate(
    group = case_when(
      analysis_type == "primary" ~ "Global",
      phenotype_domain == "asthma_lower_airway" ~ "Asthma/lower airway",
      phenotype_domain == "atopic_allergic" ~ "Atopic/allergic",
      phenotype_domain == "autoimmune_ibd" ~ "Autoimmune / IBD",
      phenotype_domain == "ent_infection" ~ "ENT / infection",
      TRUE ~ "Other domains"
    ),
    group = factor(group, levels = c("Global", "Asthma/lower airway",
                                     "Atopic/allergic", "Autoimmune / IBD",
                                     "ENT / infection", "Other domains")),
    database = factor(database, levels = c("GO_BP", "KEGG", "Reactome"))
  ) %>%
  group_by(group, database) %>%
  summarise(n_significant_terms = sum(as.numeric(n_FDR_lt_0.05), na.rm = TRUE),
            .groups = "drop")

p3s3 <- ggplot(pathway_count_audit,
               aes(n_significant_terms, group, fill = database)) +
  geom_col(position = position_dodge2(width = 0.78, preserve = "single"),
           width = 0.68, alpha = 0.92) +
  scale_fill_manual(values = pal_database) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(title = "d  Pathway count audit across databases",
       subtitle = "FDR-significant terms in full pathway-count summaries",
       x = "FDR-significant terms", y = NULL, fill = "Database") +
  theme_f3(5.8) +
  theme(legend.position = "bottom")

pairwise_fuma_plot <- pairwise_fuma_summary %>%
  mutate(
    dataset = gsub("_", " ", dataset),
    n_FDR05_terms = suppressWarnings(as.numeric(n_FDR05_terms)),
    n_pairs_with_FDR05 = suppressWarnings(as.numeric(n_pairs_with_FDR05)),
    dataset = factor(dataset, levels = rev(dataset))
  )

p3s4 <- ggplot(pairwise_fuma_plot, aes(n_pairs_with_FDR05, dataset)) +
  geom_col(width = 0.62, fill = COL["pale"], colour = "white") +
  geom_text(aes(label = paste0(n_FDR05_terms, " FDR terms")),
            hjust = -0.12, size = 2.1, colour = COL["graphite"]) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "e  Pairwise FUMA tissue enrichment audit",
       subtitle = "Pair-level tissue enrichment returned no FDR-significant GTEx tissue terms",
       x = "Disease pairs with FDR-significant tissue terms", y = NULL) +
  theme_f3(5.8)

fig3 <- p3a / p3b / p3c +
  plot_layout(heights = c(1.0, 1.15, 1.05))

supp_fig3 <- wrap_plots(
  p3d,
  p3s1,
  p3s2 | p3s3,
  p3s4,
  ncol = 1,
  heights = c(1.0, 0.95, 1.0, 0.45)
)

unlink(file.path(MAIN_PAN, paste0("Figure3D.", c("pdf", "png"))), force = TRUE)
unlink(file.path(SOURCE_OUT, c("Figure3C_pathway_breadth.tsv",
                               "Figure3C_pathway_module_summary.tsv",
                               "Figure3C_pathway_database_summary.tsv",
                               "Figure3C_domain_GTEx54_tissue_enrichment.tsv",
                               "Figure3C_GTEx54_FDR05_tissue_enrichment.tsv",
                               "Figure3C_domain_tissue_support_and_enrichment.tsv",
                               "Figure3C_domain_TWAS_preview_long.tsv",
                               "Figure3C_domain_TWAS_preview_all_tissues.tsv",
                               "Figure3C_domain_TWAS_preview_plot_data.tsv",
                               "Figure3C_TWAS_tissue_support_by_disease_group.tsv",
                               "Figure3C_TWAS_tissue_support_by_disease_group_all_tissues.tsv",
                               "Supplementary_Figure3_global_strict_GTEx_tissue_enrichment.tsv",
                               "Figure3D_tissue_support.tsv",
                               "Figure3D_cell_marker_audit.tsv")),
       force = TRUE)

save_pdf_png(p3a, file.path(MAIN_PAN, "Figure3A"), 8.8, 4.8)
save_pdf_png(p3b, file.path(MAIN_PAN, "Figure3B"), 10.2, 5.8)
save_pdf_png(p3c, file.path(MAIN_PAN, "Figure3C"), 8.8, 4.6)
save_pdf_png(fig3, file.path(MAIN_COMP, "Figure3_revised"), 10.2, 12.8)
file.copy(file.path(MAIN_COMP, "Figure3_revised.pdf"),
          file.path(OUT, "Figure3_revised.pdf"), overwrite = TRUE)
file.copy(file.path(MAIN_COMP, "Figure3_revised.png"),
          file.path(OUT, "Figure3_revised.png"), overwrite = TRUE)
save_pdf_png(supp_fig3, file.path(OUT, "Supplementary_Figure3"), 14.8, 13.0)
save_pdf_png(supp_fig3, file.path(SUPP_COMP, "Supplementary_Figure3"), 14.8, 13.0)

write_tsv(panel_a_data, file.path(SOURCE_OUT, "Figure3A_pathway_dotplot.tsv"))
write_tsv(panel_b_data, file.path(SOURCE_OUT, "Figure3B_pathway_gene_membership.tsv"))
write_tsv(panel_c_data, file.path(SOURCE_OUT, "Figure3C_TWAS_tissue_support_by_disease_group.tsv"))
write_tsv(panel_c_tissue_totals, file.path(SOURCE_OUT, "Figure3C_TWAS_tissue_support_by_disease_group_all_tissues.tsv"))
write_tsv(tissue_support, file.path(SOURCE_OUT, "Supplementary_Figure3_tissue_support.tsv"))
write_tsv(supp_fuma_tissue, file.path(SOURCE_OUT, "Supplementary_Figure3_GTEx30_general_tissue_enrichment.tsv"))
write_tsv(pairwise_fuma_summary, file.path(SOURCE_OUT, "Supplementary_Figure3_pairwise_FUMA_tissue_FDR05_summary.tsv"))
write_tsv(pairwise_fuma54_top, file.path(SOURCE_OUT, "Supplementary_Figure3_pairwise_FUMA_GTEx54_top_terms.tsv"))
write_tsv(cell_marker, file.path(SOURCE_OUT, "Supplementary_Figure3_cell_marker_audit.tsv"))
write_tsv(supp_pathways, file.path(SOURCE_OUT, "Supplementary_Figure3_KEGG_Reactome.tsv"))
write_tsv(pathway_count_audit, file.path(SOURCE_OUT, "Supplementary_Figure3_pathway_count_audit.tsv"))

qa <- data.frame(
  artifact = c("Figure3_revised.pdf", "Figure3_revised.png",
               paste0("Figure3", LETTERS[1:3], ".png"),
               "Supplementary_Figure3.pdf", "Supplementary_Figure3.png"),
  path = c(file.path(MAIN_COMP, "Figure3_revised.pdf"),
           file.path(MAIN_COMP, "Figure3_revised.png"),
           file.path(MAIN_PAN, paste0("Figure3", LETTERS[1:3], ".png")),
           file.path(SUPP_COMP, "Supplementary_Figure3.pdf"),
           file.path(SUPP_COMP, "Supplementary_Figure3.png")),
  stringsAsFactors = FALSE
) %>%
  mutate(file_size = file.info(path)$size,
         exists = file.exists(path),
         automated_status = ifelse(exists & file_size > 15000, "pass", "review"))
write_tsv(qa, file.path(QA, "Figure3_automated_render_QA.tsv"))

writeLines(capture.output(sessionInfo()), file.path(QA, "Figure3_R_sessionInfo.txt"))

message("Wrote Figure 3 v2.03 outputs to: ", OUT)
