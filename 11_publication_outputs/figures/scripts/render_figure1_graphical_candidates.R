PROJECT_ROOT <- if (dir.exists(file.path(getwd(), "figure_final"))) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
}
source(file.path(PROJECT_ROOT, "figure_final", "scripts", "common_revised.R"))

tab <- function(n) file.path(RESULT, "Table", n)
pair_axis_label <- function(pair_id) {
  paste(ifelse(pair_anchor(pair_id) == "Nasal polyps", "NP", "AR"),
        pair_partner(pair_id), sep = " - ")
}

dir.create(file.path(FIG, "main", "candidates"), recursive = TRUE,
           showWarnings = FALSE)

save_candidate <- function(plot, stem, width_in = 16.8, height_in = 10.6,
                           dpi = 600) {
  base <- file.path(FIG, "main", "candidates", stem)
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_in,
                       height = height_in, family = "sans")
  print(plot)
  grDevices::dev.off()
  grDevices::png(paste0(base, ".png"), width = width_in, height = height_in,
                 units = "in", res = dpi, type = "cairo")
  print(plot)
  grDevices::dev.off()
  file.copy(paste0(base, ".pdf"), file.path(FIG, paste0(stem, ".pdf")),
            overwrite = TRUE)
  file.copy(paste0(base, ".png"), file.path(FIG, paste0(stem, ".png")),
            overwrite = TRUE)
  invisible(base)
}

label_count <- function(n, label) paste0(format(n, big.mark = ","), "\n", label)

pair <- read_tsv(tab("Table_R2_pair_screening_and_locus_summary.tsv"))
pair$anchor <- pair_anchor(pair$pair_id)
pair$partner <- pair_partner(pair$pair_id)
pair$pair_axis <- pair_axis_label(pair$pair_id)
pair$ldsc <- tolower(pair$ldsc_sig_bh05) == "true"
pair$lava <- pair$lava_n_significant_loci > 0
pair$route <- ifelse(pair$ldsc & pair$lava, "Both",
                     ifelse(pair$ldsc, "LDSC only", "LAVA only"))
pair$domain <- gsub("_", " ", pair$phenotype_domain)

genes <- read_tsv(file.path(
  END, "final_gtexv8_twas_total_domain_fuma_enrichment",
  "final_gene_tables", "final_twas_supported_pair_gene_records.tsv"
))
recur <- read_tsv(tab("Table_R4_top_recurrent_prioritized_genes.tsv"))
path_counts <- read_tsv(tab("Table_R5_pathway_term_count_summary.tsv"))
path_global <- read_tsv(tab("Table_R6_global_top_GO_terms.tsv"))
path_domain <- read_tsv(tab("Table_R7_domain_top_enriched_terms.tsv"))
tissues <- read_tsv(tab("Table_R8_top_FUSION_TWAS_tissues.tsv"))
lcv <- read_tsv(tab("Table_R14_LCVMR_final31_summary.tsv"))

loci <- read_tsv(file.path(ROOT, "figure_v2_source_data", "locus_ab.tsv"))
locus_counts <- as.data.frame(table(loci$locus_grade))
names(locus_counts) <- c("grade", "n_loci")
gene_counts <- as.data.frame(table(genes$gene_grade))
names(gene_counts) <- c("grade", "n_records")
domain_gene <- read_tsv(tab("Table_R3_TWAS_domain_gene_counts.tsv"))

n_pairs <- nrow(pair)
n_ar <- sum(pair$anchor == "Allergic rhinitis")
n_np <- sum(pair$anchor == "Nasal polyps")
n_loci <- nrow(unique(loci[c("pair_id", "locus_id", "locus_grade")]))
n_locus_a <- sum(locus_counts$n_loci[locus_counts$grade == "A"])
n_locus_b <- sum(locus_counts$n_loci[locus_counts$grade == "B"])
n_magma_records <- sum(pair$magma_sig_gene_count, na.rm = TRUE)
n_magma_genes <- 205
n_twas_records <- nrow(genes)
n_twas_genes <- length(unique(genes$gene_symbol))
n_gene_a <- sum(gene_counts$n_records[gene_counts$grade == "Gene-A"])
n_gene_b <- sum(gene_counts$n_records[gene_counts$grade == "Gene-B"])
n_go_global <- path_counts$n_FDR_lt_0.05[
  path_counts$analysis_type == "primary" & path_counts$scope == "global" &
    path_counts$database == "GO_BP"
][1]
n_lcv_bonf <- sum(lcv$LCV_gcp_p < 0.05 / n_pairs, na.rm = TRUE)

stage_source <- data.frame(
  metric = c("nasal_phenotypes", "comparators", "final_pairs",
             "allergic_rhinitis_pairs", "nasal_polyps_pairs",
             "prioritized_loci", "locus_A", "locus_B",
             "MAGMA_records", "MAGMA_genes", "TWAS_records", "TWAS_genes",
             "Gene_A_records", "Gene_B_records", "global_GO_terms",
             "LCV_bonferroni_pairs"),
  value = c(4, 57, n_pairs, n_ar, n_np, n_loci, n_locus_a, n_locus_b,
            n_magma_records, n_magma_genes, n_twas_records, n_twas_genes,
            n_gene_a, n_gene_b, n_go_global, n_lcv_bonf)
)
write_source(stage_source, "Figure1_candidate_core_counts")
write_source(pair, "Figure1_candidate_pair_summary")

base_theme_void <- theme_void(base_family = "sans") +
  theme(plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(face = "bold", size = 18,
                                  colour = COL["graphite"], hjust = 0),
        plot.subtitle = element_text(size = 10.5, colour = COL["grey"],
                                     hjust = 0))

node_rect <- function(data, label_size = 4.0) {
  list(
    geom_rect(data = data,
              aes(xmin = x - w / 2, xmax = x + w / 2,
                  ymin = y - h / 2, ymax = y + h / 2, fill = fill),
              colour = "white", linewidth = 0.8),
    geom_text(data = data, aes(x, y, label = label),
              size = label_size, lineheight = 0.9,
              colour = COL["graphite"], fontface = "bold")
  )
}

edge_arrow <- function(data, curved = FALSE) {
  if (curved) {
    geom_curve(data = data,
               aes(x = x, y = y, xend = xend, yend = yend),
               curvature = 0.08, linewidth = 0.45,
               colour = "#9AA8B0", alpha = 0.7,
               arrow = grid::arrow(length = grid::unit(2.0, "mm"),
                                   type = "closed"))
  } else {
    geom_segment(data = data,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.45, colour = "#9AA8B0", alpha = 0.7,
                 arrow = grid::arrow(length = grid::unit(2.0, "mm"),
                                     type = "closed"))
  }
}

# Candidate 1: graphical abstract style.
c1_nodes <- data.frame(
  x = c(0.7, 2.0, 3.35, 4.75, 6.35, 7.95, 9.35),
  y = c(3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0),
  w = c(1.05, 1.25, 1.25, 1.35, 1.35, 1.35, 1.25),
  h = c(0.82, 0.82, 0.82, 0.82, 0.82, 0.82, 0.82),
  fill = c(COL["pale"], "#DCE8ED", "#EAF0F4", "#DCE8F2",
           "#CFE7E3", "#E7E0F0", "#F4E4C5"),
  label = c(label_count(61, "GWAS traits"),
            label_count(n_pairs, "genetically\nsupported pairs"),
            label_count(n_loci, "prioritized\nshared loci"),
            label_count(n_magma_records, "MAGMA\npair-gene records"),
            label_count(n_twas_genes, "TWAS-prioritized\ngenes"),
            label_count(n_go_global, "global GO BP\nterms"),
            "annotation\nboundary")
)
c1_edges <- data.frame(
  x = head(c1_nodes$x + c1_nodes$w / 2, -1),
  y = head(c1_nodes$y, -1),
  xend = tail(c1_nodes$x - c1_nodes$w / 2, -1),
  yend = tail(c1_nodes$y, -1)
)
c1_sub <- data.frame(
  x = c(2.0, 2.0, 3.35, 3.35, 6.35, 6.35, 9.35, 9.35),
  y = c(1.85, 1.25, 1.85, 1.25, 1.85, 1.25, 1.85, 1.25),
  w = c(1.05, 1.05, 1.05, 1.05, 1.0, 1.0, 1.05, 1.05),
  h = rep(0.42, 8),
  fill = c(COL["ar"], COL["np"], COL["teal"], COL["amber"],
           COL["teal"], COL["amber"], COL["blue"], COL["grey"]),
  label = c(paste0("AR ", n_ar), paste0("NP ", n_np),
            paste0("Locus-A ", n_locus_a), paste0("Locus-B ", n_locus_b),
            paste0("Gene-A ", n_gene_a), paste0("Gene-B ", n_gene_b),
            paste0("LCV Bonf. ", n_lcv_bonf), "MR secondary")
)
c1_high <- data.frame(
  x = c(5.3, 5.85, 6.4, 6.95, 7.5),
  y = rep(4.25, 5),
  label = head(recur$gene_symbol, 5),
  n = head(recur$n_pair_gene_records, 5)
)
cand1 <- ggplot() +
  edge_arrow(c1_edges) +
  node_rect(c1_nodes, 3.9) +
  node_rect(c1_sub, 3.2) +
  geom_text(data = data.frame(x = 6.4, y = 4.72,
                              label = "recurrent prioritized genes"),
            aes(x, y, label = label), size = 4.1, fontface = "bold",
            colour = COL["graphite"]) +
  geom_point(data = c1_high, aes(x, y, size = n), colour = COL["teal"],
             alpha = 0.88) +
  geom_text(data = c1_high, aes(x, y - 0.36, label = label),
            size = 3.3, fontface = "bold", colour = COL["graphite"]) +
  scale_size_continuous(range = c(5, 9), guide = "none") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.05, 10.1), ylim = c(0.55, 5.0), clip = "off") +
  labs(title = "Candidate 1: graphical summary of nasal shared-genetic prioritization",
       subtitle = "A left-to-right visual abstract: pair selection -> locus grading -> gene prioritization -> biological interpretation -> secondary annotation") +
  base_theme_void

# Candidate 2: evidence landscape style.
pair_ord <- pair[order(pair$anchor, -pair$lava_n_significant_loci,
                       -pair$magma_sig_gene_count), ]
pair_ord$pair_axis <- factor(pair_ord$pair_axis, levels = rev(pair_ord$pair_axis))
pair_long <- rbind(
  data.frame(pair_axis = pair_ord$pair_axis, anchor = pair_ord$anchor,
             domain = pair_ord$domain, layer = "LDSC", value = pair_ord$ldsc + 0),
  data.frame(pair_axis = pair_ord$pair_axis, anchor = pair_ord$anchor,
             domain = pair_ord$domain, layer = "LAVA", value = pmin(pair_ord$lava_n_significant_loci, 25) / 25),
  data.frame(pair_axis = pair_ord$pair_axis, anchor = pair_ord$anchor,
             domain = pair_ord$domain, layer = "MAGMA", value = pmin(pair_ord$magma_sig_gene_count, 90) / 90)
)
pair_long$layer <- factor(pair_long$layer, levels = c("LDSC", "LAVA", "MAGMA"))
c2a <- ggplot(pair_long, aes(layer, pair_axis, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  facet_grid(anchor ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = COL["absent"], high = COL["blue"],
                      name = "Evidence\nintensity") +
  labs(title = "Evidence landscape across 31 disease pairs",
       subtitle = "Rows retain AR/NP identity; columns summarize global sharing, local sharing and gene evidence",
       x = NULL, y = NULL) +
  theme_pub(6.0) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 5.0),
        legend.position = "bottom")

domain_sum <- aggregate(cbind(lava_n_significant_loci, magma_sig_gene_count) ~
                          phenotype_domain, pair, sum, na.rm = TRUE)
domain_sum$domain <- tools::toTitleCase(gsub("_", " ", domain_sum$phenotype_domain))
domain_gene2 <- domain_gene
domain_gene2$domain <- tools::toTitleCase(gsub("_", " ", domain_gene2$phenotype_domain))
domain_sum <- merge(domain_sum, domain_gene2, by = "domain", all.x = TRUE)
domain_sum[is.na(domain_sum)] <- 0
c2b <- ggplot(domain_sum,
              aes(magma_sig_gene_count, reorder(domain, magma_sig_gene_count))) +
  geom_segment(aes(x = 0, xend = magma_sig_gene_count,
                   yend = reorder(domain, magma_sig_gene_count)),
               linewidth = 3.0, colour = COL["pale"]) +
  geom_point(aes(size = n_unique_genes, colour = lava_n_significant_loci),
             alpha = 0.9) +
  scale_colour_gradient(low = COL["amber"], high = COL["teal"],
                        name = "LAVA loci") +
  scale_size_continuous(range = c(4.0, 9.0), name = "TWAS genes") +
  labs(title = "Domain-level concentration",
       subtitle = "Asthma/lower-airway and atopic/allergic domains carry the strongest gene/pathway signal",
       x = "MAGMA-positive pair-gene records", y = NULL) +
  theme_pub(7) + theme(legend.position = "bottom")

route <- as.data.frame(table(pair$route)); names(route) <- c("route", "n")
route$route <- factor(route$route, levels = c("Both", "LDSC only", "LAVA only"))
c2c <- ggplot(route, aes(route, n, fill = route)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), vjust = -0.2, fontface = "bold", size = 5) +
  scale_fill_manual(values = c("Both" = colv("teal"),
                               "LDSC only" = colv("blue"),
                               "LAVA only" = colv("amber"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Pair-selection routes", x = NULL, y = "Disease pairs") +
  theme_pub(7) + theme(legend.position = "none")
cand2 <- (c2a | (c2b / c2c)) +
  plot_layout(widths = c(1.35, 1.0), heights = c(1, 0.8)) +
  plot_annotation(
    title = "Candidate 2: disease-pair evidence map",
    subtitle = "A data-rich Figure 1 that foregrounds which disease pairs enter the shared-locus and gene-prioritization pipeline"
  ) &
  theme(plot.title = element_text(face = "bold", size = 18,
                                  colour = COL["graphite"]),
        plot.subtitle = element_text(size = 10.5, colour = COL["grey"]))

# Candidate 3: biology-first ecosystem style.
ecosystem_nodes <- data.frame(
  name = c("31 genetically\nsupported pairs", "181 Locus-A/B\nregions",
           "123 TWAS-prioritized\ngenes", "Immune and microbial\nresponse pathways",
           "Tissue support\ncontext", "FUMA/drug/MR\nsecondary layers"),
  x = c(0, 2.1, 4.1, 6.1, 4.1, 6.1),
  y = c(2.3, 2.3, 2.3, 2.3, 0.92, 0.92),
  w = c(1.45, 1.45, 1.55, 1.75, 1.55, 1.55),
  h = c(0.78, 0.78, 0.78, 0.78, 0.70, 0.70),
  fill = c("#EAF0F4", "#E7E0F0", "#CFE7E3", "#F4E4C5",
           "#DCE8F2", COL["pale"])
)
ecosystem_nodes$label <- ecosystem_nodes$name
ecosystem_edges <- data.frame(
  x = c(0.75, 2.85, 4.9, 4.1, 4.9),
  y = c(2.3, 2.3, 2.3, 1.92, 1.92),
  xend = c(1.35, 3.35, 5.25, 4.1, 5.35),
  yend = c(2.3, 2.3, 2.3, 1.28, 1.28)
)
top_terms <- head(path_global$term_name, 5)
term_cloud <- data.frame(
  x = c(6.0, 6.68, 7.36, 6.34, 7.05),
  y = c(3.66, 3.66, 3.66, 3.18, 3.18),
  label = wrap_text(top_terms, 18),
  n = head(path_global$overlap_count, 5)
)
tissue_cloud <- head(tissues, 5)
if (!"tissue" %in% names(tissue_cloud)) names(tissue_cloud)[1] <- "tissue"
tissue_cloud$x <- seq(3.15, 5.05, length.out = nrow(tissue_cloud))
tissue_cloud$y <- c(0.15, 0.05, 0.18, 0.05, 0.16)[seq_len(nrow(tissue_cloud))]
tissue_cloud$label <- wrap_text(gsub("_", " ", tissue_cloud[[1]]), 14)
gene_cloud <- head(recur, 10)
gene_cloud$x <- rep(c(2.85, 3.45, 4.05, 4.65, 5.25), 2)
gene_cloud$y <- rep(c(3.55, 3.10), each = 5)
cand3 <- ggplot() +
  edge_arrow(ecosystem_edges) +
  node_rect(ecosystem_nodes, 3.7) +
  geom_point(data = gene_cloud,
             aes(x, y, size = n_pair_gene_records),
             colour = COL["teal"], alpha = 0.86) +
  geom_text(data = gene_cloud, aes(x, y - 0.22, label = gene_symbol),
            size = 2.9, fontface = "bold", colour = COL["graphite"]) +
  geom_label(data = term_cloud, aes(x, y, label = label),
             fill = "#FFF6ED", colour = COL["graphite"],
             linewidth = 0, lineheight = 0.84, size = 2.6) +
  geom_text(data = tissue_cloud, aes(x, y, label = label),
            size = 3.0, colour = COL["blue"], fontface = "bold",
            lineheight = 0.85) +
  geom_text(data = data.frame(
    x = c(0, 2.1, 4.1, 6.1, 4.1, 6.1),
    y = c(1.55, 1.55, 1.55, 1.55, 0.45, 0.45),
    label = c(paste0("AR ", n_ar, " | NP ", n_np),
              paste0("A ", n_locus_a, " | B ", n_locus_b),
              paste0("A ", n_gene_a, " | B ", n_gene_b, " records"),
              paste0(n_go_global, " global GO BP terms"),
              "support, not origin",
              "annotation, not grading")
  ), aes(x, y, label = label), size = 3.35, colour = COL["grey"]) +
  scale_size_continuous(range = c(2.8, 6.4), guide = "none") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-0.85, 8.25), ylim = c(-0.25, 4.0), clip = "off") +
  labs(title = "Candidate 3: biology-first graphical abstract",
       subtitle = "The first read is the biological convergence of prioritized genes; workflow counts remain as supporting anchors") +
  base_theme_void

save_candidate(cand1, "Figure1_candidate1_graphical_summary")
save_candidate(cand2, "Figure1_candidate2_pair_evidence_map")
save_candidate(cand3, "Figure1_candidate3_biology_first")

contact <- (cand1 / cand2 / cand3) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Figure 1 candidate comparison",
    subtitle = "Three alternative graphical-summary directions rendered from the same project source data"
  ) &
  theme(plot.title = element_text(face = "bold", size = 18,
                                  colour = COL["graphite"]),
        plot.subtitle = element_text(size = 10.5, colour = COL["grey"]))
save_candidate(contact, "Figure1_candidate_contact_sheet",
               width_in = 13.5, height_in = 18.0)

candidate_manifest <- data.frame(
  candidate = c("Figure1_candidate1_graphical_summary",
                "Figure1_candidate2_pair_evidence_map",
                "Figure1_candidate3_biology_first"),
  concept = c("Graphical abstract pipeline",
              "Disease-pair evidence landscape",
              "Biology-first graphical abstract"),
  intended_use = c("Best if Figure 1 should replace a technical flowchart with an overview of the whole study.",
                   "Best if readers must first understand which disease-pair evidence drove the analysis.",
                   "Best if the first figure should sell the biological message before technical detail."),
  inspiration = c("Graphical-abstract / sequential overview papers",
                  "Biobank overview and trait-map Figure 1 papers",
                  "Nature-style hero biological summary with workflow as support")
)
write_source(candidate_manifest, "Figure1_candidate_manifest")

message("Wrote Figure 1 candidates to: ", file.path(FIG, "main", "candidates"))
