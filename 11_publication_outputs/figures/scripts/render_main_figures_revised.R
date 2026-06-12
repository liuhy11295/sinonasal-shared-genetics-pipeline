PROJECT_ROOT <- if (dir.exists(file.path(getwd(), "figure_final"))) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
}
source(file.path(PROJECT_ROOT, "figure_final", "scripts", "common_revised.R"))

tab <- function(n) file.path(RESULT, "Table", n)

# Shared source data
pair <- read_tsv(tab("Table_R2_pair_screening_and_locus_summary.tsv"))
pair$anchor <- pair_anchor(pair$pair_id)
pair$partner <- pair_partner(pair$pair_id)
pair$ldsc <- tolower(pair$ldsc_sig_bh05) == "true"
pair$lava <- pair$lava_n_significant_loci > 0
pair$route <- ifelse(pair$ldsc & pair$lava, "Both",
                     ifelse(pair$ldsc, "LDSC only", "LAVA only"))

genes <- read_tsv(file.path(
  END, "final_gtexv8_twas_total_domain_fuma_enrichment",
  "final_gene_tables", "final_twas_supported_pair_gene_records.tsv"
))
recur <- read_tsv(tab("Table_R4_top_recurrent_prioritized_genes.tsv"))
path_counts <- read_tsv(tab("Table_R5_pathway_term_count_summary.tsv"))
path_global <- read_tsv(tab("Table_R6_global_top_GO_terms.tsv"))
path_domain <- read_tsv(tab("Table_R7_domain_top_enriched_terms.tsv"))
fuma_long <- read_tsv(file.path(
  END, "fuma_snp_annotation_concordance", "parsed_fuma_annotation",
  "fuma_mapped_gene_TWAS_concordance.tsv"
))
chain <- read_tsv(file.path(
  END, "fuma_snp_annotation_concordance", "parsed_fuma_annotation",
  "variant_to_TWAS_gene_concordance.tsv"
))

# Figure 1A
n1 <- data.frame(
  x = c(1, 1, 2.5, 4, 4),
  y = c(2.8, 1.2, 2, 2.8, 1.2),
  label = c("4 nasal\nphenotypes", "57 comparator\ntraits",
            "Heritability and\npower audit", "Allergic rhinitis\n13 pairs",
            "Nasal polyps\n18 pairs"),
  fill = c(COL["pale"], COL["pale"], "#DCE8ED", "#F2DCE3", "#DCE8F2")
)
e1 <- data.frame(
  x = c(1.35, 1.35, 2.85, 2.85), y = c(2.8, 1.2, 2.1, 1.9),
  xend = c(2.1, 2.1, 3.55, 3.55), yend = c(2.2, 1.8, 2.8, 1.2)
)
p1a <- flow_plot(n1, e1, "Phenotype screening",
                 "Power assessment retained two nasal anchors and 31 pairs")

# Figure 1B
atlas <- rbind(
  data.frame(pair_id = pair$pair_id, partner = pair$partner,
             anchor = pair$anchor, domain = pair$phenotype_domain,
             metric = "LDSC FDR", value = ifelse(pair$ldsc, 1, 0)),
  data.frame(pair_id = pair$pair_id, partner = pair$partner,
             anchor = pair$anchor, domain = pair$phenotype_domain,
             metric = "LAVA loci", value = pmin(pair$lava_n_significant_loci, 20) / 20)
)
atlas$partner <- factor(atlas$partner, levels = rev(unique(pair$partner)))
p1b <- ggplot(atlas, aes(metric, partner, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(anchor ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = COL["pale"], high = COL["blue"],
                      name = "Evidence") +
  labs(title = "Genome-wide and local sharing",
       subtitle = "21 both | 2 LDSC only | 8 LAVA only",
       x = NULL, y = NULL) +
  theme_pub(6.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")

# Figure 1C
n1c <- data.frame(
  x = c(1, 2.4, 2.4, 4, 4),
  y = c(2, 2.8, 1.2, 2.7, 1.3),
  label = c("Candidate regions\nMTAG + PLACO + CPASSOC + LAVA",
            "Coloc ABF\n182 positive", "SuSiE-coloc\n276 entries\n285 analyses\n76 positive",
            "Locus-A\n68 regions", "Locus-B\n113 regions"),
  fill = c(COL["pale"], "#DCE8F2", "#E5E1F2", "#CFE7E3", "#F4E4C5")
)
e1c <- data.frame(
  x = c(1.45, 1.45, 2.85, 2.85), y = c(2.1, 1.9, 2.8, 1.2),
  xend = c(2.0, 2.0, 3.55, 3.55), yend = c(2.8, 1.2, 2.7, 1.3)
)
p1c <- flow_plot(n1c, e1c, "Regional resolution and grading",
                 "Coloc and SuSiE are parallel summaries, not a numerical funnel")
write_source(pair, "Figure1_pair_evidence")
save_panel_set(list(A = p1a, B = p1b, C = p1c), "Figure1",
               panel_width = 92, panel_height = 98,
               composite_height = 165, design = "AB\nCC",
               composite_stem = "Figure1_revised")
copy_composite_to_root("Figure1_revised")

# Figure 2A
if (FALSE) {
funnel <- data.frame(
  stage = factor(c("Locus-A/B regions", "MAGMA records", "MAGMA genes",
                   "TWAS records", "TWAS genes"),
                 levels = rev(c("Locus-A/B regions", "MAGMA records",
                                "MAGMA genes", "TWAS records", "TWAS genes"))),
  n = c(181, 566, 205, 336, 123),
  type = c("regions", "records", "genes", "records", "genes")
)
p2a <- ggplot(funnel, aes(n, stage, fill = type)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n), hjust = -0.15, fontface = "bold", size = 3) +
  scale_fill_manual(values = c(genes = colv("teal"), records = colv("blue"),
                               regions = colv("amber"))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(title = "Candidate-restricted gene prioritization",
       subtitle = "Records and unique entities use separate denominators",
       x = "Count", y = NULL, fill = NULL) +
  theme_pub(7) + theme(legend.position = "bottom")

# Figure 2B
top_genes <- head(recur$gene_symbol, 24)
gm <- genes[genes$gene_symbol %in% top_genes,
            c("pair_id", "phenotype_domain", "gene_symbol", "gene_grade")]
gm$pair <- pair_partner(gm$pair_id)
gm$gene_symbol <- factor(gm$gene_symbol, levels = rev(top_genes))
p2b <- ggplot(gm, aes(pair, gene_symbol, fill = gene_grade)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("Gene-A" = colv("teal"),
                               "Gene-B" = colv("amber"))) +
  labs(title = "Recurrent support across disease pairs",
       x = NULL, y = NULL, fill = NULL) +
  theme_pub(5.7) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1),
        legend.position = "bottom")

# Figure 2C
tissue_long <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
  z <- genes[i, ]
  tissues <- strsplit(ifelse(is.na(z$significant_tissue_list), "",
                            z$significant_tissue_list), ",", fixed = TRUE)[[1]]
  tissues <- tissues[nzchar(tissues)]
  if (!length(tissues) || !(z$gene_symbol %in% top_genes)) return(NULL)
  data.frame(gene_symbol = z$gene_symbol, tissue = tissues,
             gene_grade = z$gene_grade)
}))
tissue_freq <- sort(table(tissue_long$tissue), decreasing = TRUE)
top_tissues <- names(head(tissue_freq, 15))
tissue_long <- tissue_long[tissue_long$tissue %in% top_tissues, ]
tissue_long$gene_symbol <- factor(tissue_long$gene_symbol, levels = rev(top_genes))
tissue_long$tissue <- factor(tissue_long$tissue, levels = top_tissues)
p2c <- ggplot(tissue_long, aes(tissue, gene_symbol)) +
  geom_point(aes(colour = gene_grade), size = 1.45, alpha = 0.8) +
  scale_colour_manual(values = c("Gene-A" = colv("teal"),
                                 "Gene-B" = colv("amber"))) +
  labs(title = "Distributed TWAS tissue support",
       subtitle = "Top 15 tissues among the displayed recurrent genes",
       x = NULL, y = NULL, colour = NULL) +
  theme_pub(5.7) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        legend.position = "bottom")
write_source(gm, "Figure2_gene_pair_matrix")
write_source(tissue_long, "Figure2_gene_tissue_support")
save_panel_set(list(A = p2a, B = p2b, C = p2c), "Figure2",
               panel_width = 100, panel_height = 105,
               composite_height = 175, design = "AB\nCC")
}

# Figure 3A
if (FALSE) {
pc <- path_counts[path_counts$analysis_type %in%
                    c("primary", "strict", "primary_domain"), ]
pc$analysis_label <- ifelse(is.na(pc$phenotype_domain) | pc$phenotype_domain == "",
                            pc$analysis_name, pc$phenotype_domain)
pc$analysis_label <- gsub("_", " ", pc$analysis_label)
pc$database <- gsub("_", " ", pc$database)
p3a <- ggplot(pc, aes(database, analysis_label, fill = n_FDR_lt_0.05)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n_FDR_lt_0.05), size = 2.5) +
  scale_fill_gradient(low = COL["pale"], high = COL["teal"]) +
  labs(title = "Functional annotation atlas",
       subtitle = "FDR-significant terms; zeros remain visible",
       x = NULL, y = NULL, fill = "Terms") +
  theme_pub(6.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")

# Figure 3B
pg <- data.frame(
  group = "Global", term_name = path_global$term_name,
  overlap_count = path_global$overlap_count, FDR = path_global$FDR,
  database = path_global$database
)
pd <- data.frame(
  group = gsub("_", " ", path_domain$phenotype_domain),
  term_name = path_domain$term_name,
  overlap_count = path_domain$overlap_count, FDR = path_domain$FDR,
  database = path_domain$database
)
pdot <- rbind(pg, pd)
pdot <- pdot[order(pdot$FDR), ]
pdot <- do.call(rbind, lapply(split(pdot, pdot$group), head, 2))
pdot$term <- factor(wrap_text(pdot$term_name, 27),
                    levels = rev(unique(wrap_text(pdot$term_name, 27))))
p3b <- ggplot(pdot, aes(-log10(FDR), term, size = overlap_count,
                        colour = group, shape = database)) +
  geom_point(alpha = 0.85) +
  scale_colour_manual(values = c("Global" = colv("graphite"),
                                 "asthma lower airway" = colv("blue"),
                                 "atopic allergic" = colv("rose"))) +
  labs(title = "Leading non-redundant pathways",
       subtitle = "Global (graphite), asthma/lower-airway (blue), atopic/allergic (rose)",
       x = expression(-log[10](FDR)), y = NULL,
       size = "Overlap", colour = NULL, shape = NULL) +
  theme_pub(6) + theme(legend.position = "none")

# Figure 3C
lane <- data.frame(
  layer = factor(c("FUSION tissue support", "FUMA tissue enrichment",
                   "Global/domain cell type"),
                 levels = rev(c("FUSION tissue support",
                                "FUMA tissue enrichment",
                                "Global/domain cell type"))),
  count = c(123, 24, 0),
  status = c("Supportive", "Post hoc", "No robust signal"),
  note = c("49 GTEx tissues; lung and whole blood each 64 records",
           "17 GTEx54 + 7 GTEx30 significant rows",
           "Only isolated pair-level exploratory findings")
)
p3c <- ggplot(lane, aes(count, layer, colour = status)) +
  geom_segment(aes(x = 0, xend = count, yend = layer), linewidth = 2.5,
               alpha = 0.28) +
  geom_point(size = 4) +
  geom_text(aes(label = note), hjust = -0.03, size = 2.4,
            colour = COL["graphite"]) +
  scale_colour_manual(values = c("Supportive" = colv("blue"),
                                 "Post hoc" = colv("amber"),
                                 "No robust signal" = colv("grey"))) +
  scale_x_continuous(limits = c(0, 360), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Tissue and cell evidence boundary",
       x = "Summary count", y = NULL, colour = NULL) +
  theme_pub(6.5) + theme(legend.position = "bottom")
write_source(pc, "Figure3_pathway_atlas")
write_source(pdot, "Figure3_pathway_dotplot")
save_panel_set(list(A = p3a, B = p3b, C = p3c), "Figure3",
               panel_width = 98, panel_height = 90,
               composite_height = 165, design = "AB\nCC")
}

drug <- read_tsv(file.path(END, "network_pharmacology_res",
                           "drug_recurrence_summary.tsv"))
drug <- drug[drug$approved_any %in% c(TRUE, "True"), ]
split_semicolon <- function(x) {
  x <- unlist(strsplit(paste(x, collapse = ";"), ";", fixed = TRUE))
  sort(unique(x[nzchar(x)]))
}
display_drug <- function(x) tools::toTitleCase(tolower(x))

drug_gene_map <- do.call(rbind, lapply(seq_len(nrow(drug)), function(i) {
  gene_set <- split_semicolon(drug$target_gene_list[i])
  pair_set <- split_semicolon(drug$pair_list[i])
  if (!length(gene_set)) return(NULL)
  data.frame(
    drug_name = drug$drug_name[i],
    display_drug_name = display_drug(drug$drug_name[i]),
    gene_symbol = gene_set,
    pair_list = paste(pair_set, collapse = ";"),
    n_pairs_for_gene_drug = length(pair_set),
    analysis_type = drug$analysis_type[i],
    stringsAsFactors = FALSE
  )
}))
drug_gene_map <- do.call(rbind, lapply(
  split(drug_gene_map, paste(drug_gene_map$drug_name,
                             drug_gene_map$gene_symbol, sep = "__")),
  function(z) {
    pair_set <- split_semicolon(z$pair_list)
    data.frame(
      drug_name = z$drug_name[1],
      display_drug_name = z$display_drug_name[1],
      gene_symbol = z$gene_symbol[1],
      pair_list = paste(pair_set, collapse = ";"),
      n_pairs_for_gene_drug = length(pair_set),
      analysis_type = paste(sort(unique(z$analysis_type)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
))

gene_drug_summary <- do.call(rbind, lapply(
  split(drug_gene_map, drug_gene_map$gene_symbol), function(z) {
    pair_set <- split_semicolon(z$pair_list)
    drug_rank <- aggregate(n_pairs_for_gene_drug ~ drug_name + display_drug_name,
                           z, max)
    drug_rank <- drug_rank[order(-drug_rank$n_pairs_for_gene_drug,
                                 drug_rank$drug_name), ]
    data.frame(
      gene_symbol = z$gene_symbol[1],
      n_approved_drugs = length(unique(z$drug_name)),
      drug_supported_pairs = length(pair_set),
      top_approved_drugs = paste(head(drug_rank$display_drug_name, 2),
                                 collapse = "; "),
      all_approved_drugs = paste(drug_rank$display_drug_name, collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
))

gene_evidence <- do.call(rbind, lapply(split(genes, genes$gene_symbol), function(z) {
  domains <- unique(z$phenotype_domain[!is.na(z$phenotype_domain) &
                                         nzchar(z$phenotype_domain)])
  n_tissues <- suppressWarnings(as.numeric(z$n_significant_tissues))
  twas_fdr <- suppressWarnings(as.numeric(z$best_twas_fdr))
  best_fdr <- suppressWarnings(min(twas_fdr, na.rm = TRUE))
  if (!is.finite(best_fdr)) best_fdr <- NA_real_
  max_tissues <- suppressWarnings(max(n_tissues, na.rm = TRUE))
  if (!is.finite(max_tissues)) max_tissues <- 0
  data.frame(
    gene_symbol = z$gene_symbol[1],
    gene_grade = ifelse(any(z$gene_grade == "Gene-A"), "Gene-A", "Gene-B"),
    n_pair_gene_records = length(unique(z$pair_id)),
    n_geneA_records = sum(z$gene_grade == "Gene-A", na.rm = TRUE),
    n_geneB_records = sum(z$gene_grade == "Gene-B", na.rm = TRUE),
    n_phenotype_domains = length(domains),
    max_twas_tissues = max_tissues,
    best_twas_fdr = best_fdr,
    n_strict_positive = sum(z$strict_positive %in%
                              c(TRUE, "True", "TRUE", "true"), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

target_gene_summary <- merge(gene_drug_summary, gene_evidence,
                             by = "gene_symbol", all.x = TRUE)
target_gene_summary$gene_grade[is.na(target_gene_summary$gene_grade)] <- "Not prioritized"
numeric_cols <- c("n_pair_gene_records", "n_geneA_records", "n_geneB_records",
                  "n_phenotype_domains", "max_twas_tissues",
                  "n_strict_positive")
for (cc in numeric_cols) {
  target_gene_summary[[cc]][is.na(target_gene_summary[[cc]])] <- 0
}
target_gene_summary$gene_grade_rank <- ifelse(target_gene_summary$gene_grade == "Gene-A", 1, 0)
target_gene_summary <- target_gene_summary[order(
  -target_gene_summary$drug_supported_pairs,
  -target_gene_summary$n_pair_gene_records,
  -target_gene_summary$n_approved_drugs,
  -target_gene_summary$gene_grade_rank,
  target_gene_summary$gene_symbol
), ]
network_genes <- head(target_gene_summary, 10)
selected_gene_levels <- network_genes$gene_symbol
edge_pool <- merge(drug_gene_map, network_genes, by = "gene_symbol")
edge_pool <- edge_pool[edge_pool$gene_symbol %in% selected_gene_levels, ]

must_keep_drugs <- unique(unlist(lapply(split(edge_pool, edge_pool$gene_symbol),
                                        function(z) {
  z <- z[order(-z$n_pairs_for_gene_drug, z$drug_name), ]
  head(z$drug_name, 2)
})))
drug_rank <- do.call(rbind, lapply(split(edge_pool, edge_pool$drug_name),
                                   function(z) {
  data.frame(
    drug_name = z$drug_name[1],
    display_drug_name = z$display_drug_name[1],
    max_pair_support = max(z$n_pairs_for_gene_drug, na.rm = TRUE),
    n_connected_genes = length(unique(z$gene_symbol)),
    must_keep = z$drug_name[1] %in% must_keep_drugs,
    stringsAsFactors = FALSE
  )
}))
drug_rank <- drug_rank[order(!drug_rank$must_keep,
                             -drug_rank$max_pair_support,
                             -drug_rank$n_connected_genes,
                             drug_rank$drug_name), ]
selected_drugs <- head(drug_rank$drug_name, 20)
network_edges <- edge_pool[edge_pool$drug_name %in% selected_drugs, ]
network_edges$gene_symbol <- factor(network_edges$gene_symbol,
                                    levels = selected_gene_levels)
network_edges$drug_name <- factor(network_edges$drug_name,
                                  levels = selected_drugs)

gene_nodes <- network_genes[, c("gene_symbol", "gene_grade",
                                "n_pair_gene_records",
                                "n_approved_drugs",
                                "drug_supported_pairs")]
gene_nodes$x <- 0
gene_nodes$y <- seq(from = 9.6, to = 0.6, length.out = nrow(gene_nodes))
gene_nodes$gene_symbol <- factor(gene_nodes$gene_symbol,
                                 levels = selected_gene_levels)

drug_nodes <- drug_rank[match(selected_drugs, drug_rank$drug_name),
                        c("drug_name", "display_drug_name",
                          "max_pair_support", "n_connected_genes")]
drug_nodes$rank <- seq_len(nrow(drug_nodes))
drug_nodes$side <- ifelse(drug_nodes$rank %% 2 == 1, "left", "right")
drug_nodes$x <- ifelse(drug_nodes$side == "left", -2.1, 2.1)
drug_nodes$label_x <- ifelse(drug_nodes$side == "left", -2.34, 2.34)
drug_nodes$hjust <- ifelse(drug_nodes$side == "left", 1, 0)
drug_nodes$y <- NA_real_
for (ss in c("left", "right")) {
  idx <- which(drug_nodes$side == ss)
  drug_nodes$y[idx] <- seq(from = 9.8, to = 0.4, length.out = length(idx))
}

edge_draw <- merge(network_edges, drug_nodes, by = "drug_name")
edge_draw <- merge(edge_draw, gene_nodes[, c("gene_symbol", "x", "y")],
                   by = "gene_symbol", suffixes = c("_drug", "_gene"))
edge_draw$display_drug_name <- edge_draw$display_drug_name.x
edge_draw <- edge_draw[order(edge_draw$n_pairs_for_gene_drug), ]
panel_a_source <- edge_draw[, c("drug_name", "display_drug_name",
                                "gene_symbol", "gene_grade",
                                "n_pairs_for_gene_drug", "pair_list",
                                "n_pair_gene_records", "n_approved_drugs",
                                "drug_supported_pairs", "analysis_type")]

figure4_gene_colours <- c("Gene-A" = "#00A087",
                          "Gene-B" = "#F39B7F")

p4a <- ggplot() +
  geom_curve(
    data = edge_draw[edge_draw$side == "left", ],
    aes(x = x_drug, y = y_drug, xend = x_gene, yend = y_gene,
        linewidth = n_pairs_for_gene_drug),
    curvature = -0.15, colour = "#AEB8BF", alpha = 0.46
  ) +
  geom_curve(
    data = edge_draw[edge_draw$side == "right", ],
    aes(x = x_drug, y = y_drug, xend = x_gene, yend = y_gene,
        linewidth = n_pairs_for_gene_drug),
    curvature = 0.15, colour = "#AEB8BF", alpha = 0.46
  ) +
  geom_point(
    data = drug_nodes,
    aes(x, y, size = max_pair_support),
    shape = 21, fill = "#E6E1F2", colour = COL["violet"],
    stroke = 0.25, alpha = 0.95
  ) +
  geom_text(
    data = drug_nodes,
    aes(label_x, y, label = display_drug_name, hjust = hjust),
    size = 1.82, colour = COL["graphite"], lineheight = 0.84
  ) +
  geom_point(
    data = gene_nodes,
    aes(x, y, size = n_pair_gene_records, fill = gene_grade),
    shape = 21, colour = "white", stroke = 0.45, alpha = 0.98
  ) +
  geom_text(
    data = gene_nodes,
    aes(x = 0.18, y = y, label = gene_symbol),
    inherit.aes = FALSE, hjust = 0, size = 2.15,
    colour = COL["graphite"], fontface = "bold"
  ) +
  annotate("text", x = 0, y = 10.35, label = "Prioritized target genes",
           size = 2.5, fontface = "bold", colour = COL["graphite"]) +
  annotate("text", x = -2.1, y = 10.35, label = "Approved drugs",
           size = 2.35, fontface = "bold", colour = COL["grey"]) +
  annotate("text", x = 2.1, y = 10.35, label = "Approved drugs",
           size = 2.35, fontface = "bold", colour = COL["grey"]) +
  scale_fill_manual(values = figure4_gene_colours) +
  scale_size_continuous(range = c(2.6, 7.6), guide = "none") +
  scale_linewidth_continuous(range = c(0.18, 1.05), guide = "none") +
  coord_cartesian(xlim = c(-2.86, 2.86), ylim = c(0, 10.7),
                  clip = "off") +
  labs(title = "Genetically prioritized drug-target network",
       subtitle = "Edges connect approved drugs to prioritized target genes; line width reflects disease-pair support",
       x = NULL, y = NULL, fill = NULL) +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 8.2,
                              colour = COL["graphite"], hjust = 0),
    plot.subtitle = element_text(size = 5.8, colour = COL["grey"], hjust = 0),
    legend.position = "bottom",
    legend.text = element_text(size = 5.5),
    plot.margin = margin(6, 8, 6, 8)
  )

# Figure 4B
panel_b_source <- network_edges[, c("drug_name", "display_drug_name",
                                    "gene_symbol", "gene_grade",
                                    "n_pairs_for_gene_drug",
                                    "n_pair_gene_records",
                                    "n_approved_drugs",
                                    "drug_supported_pairs",
                                    "analysis_type")]
panel_b_source$gene_symbol <- factor(panel_b_source$gene_symbol,
                                      levels = selected_gene_levels)
panel_b_source$display_drug_name <- factor(
  panel_b_source$display_drug_name,
  levels = rev(drug_nodes$display_drug_name)
)
matrix_bg <- expand.grid(
  gene_symbol = factor(selected_gene_levels, levels = selected_gene_levels),
  display_drug_name = factor(rev(levels(panel_b_source$display_drug_name)),
                             levels = levels(panel_b_source$display_drug_name))
)

p4b <- ggplot(matrix_bg, aes(gene_symbol, display_drug_name)) +
  geom_tile(fill = "#F3F6F7", colour = "white", linewidth = 0.25) +
  geom_point(
    data = panel_b_source,
    aes(size = n_pairs_for_gene_drug, fill = gene_grade),
    shape = 21,
    colour = "white",
    stroke = 0.22,
    alpha = 0.98
  ) +
  scale_fill_manual(values = figure4_gene_colours) +
  scale_size_continuous(range = c(1.5, 4.8),
                        breaks = sort(unique(pmin(panel_b_source$n_pairs_for_gene_drug,
                                                  12)))) +
  guides(
    size = guide_legend(
      override.aes = list(shape = 21, fill = COL["graphite"],
                          colour = COL["graphite"], alpha = 1, stroke = 0.2)
    )
  ) +
  labs(title = "Approved drug-target recurrence matrix",
       subtitle = "Point size shows disease-pair support for each drug-target edge",
       x = NULL, y = NULL, fill = NULL, size = "Pairs") +
  theme_pub(5.4) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# Figure 4C
lcv <- read_tsv(file.path(ROOT, "figure_v2_source_data", "lcv.tsv"))
mr_direction <- read_tsv(file.path(ROOT, "figure_v2_source_data",
                                   "mr_directions.tsv"))
mr_robust <- read_tsv(file.path(ROOT, "figure_v2_source_data",
                                "mr_robustness.tsv"))
robust_direction <- aggregate(
  pass ~ direction, mr_robust,
  function(x) all(x %in% c(TRUE, "True"))
)

direction_summary <- data.frame(
  analysis = factor(
    c("LCV Bonferroni", "LCV nominal", "MR eligible",
      "MR Bonferroni", "MR strict robust", "MR insufficient"),
    levels = rev(c("LCV Bonferroni", "LCV nominal", "MR eligible",
                   "MR Bonferroni", "MR strict robust", "MR insufficient"))
  ),
  n = c(
    sum(lcv$LCV_gcp_p < 0.05 / 31, na.rm = TRUE),
    sum(lcv$LCV_gcp_p < 0.05, na.rm = TRUE),
    sum(mr_direction$status == "eligible", na.rm = TRUE),
    sum(mr_direction$status2 == "Bonferroni", na.rm = TRUE),
    sum(robust_direction$pass, na.rm = TRUE),
    sum(mr_direction$status != "eligible", na.rm = TRUE)
  ),
  family = c("LCV", "LCV", "MR", "MR", "MR", "MR")
)

p4d <- ggplot(direction_summary, aes(n, analysis, fill = family)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n), hjust = -0.16, size = 2.8,
            fontface = "bold", colour = COL["graphite"]) +
  scale_fill_manual(values = c("LCV" = unname(COL["teal"]),
                               "MR" = unname(COL["blue"]))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Directional analyses remain secondary",
       subtitle = "Strict robust MR requires all method and sensitivity checks",
       x = "Number of pairs or bidirectional tests", y = NULL, fill = NULL) +
  theme_pub(6.2) +
  theme(legend.position = "bottom")

unlink(file.path(FIG, "main", "panels", paste0("Figure4D.", c("pdf", "png"))),
       force = TRUE)
unlink(file.path(FIG, "source_data",
                 c("Figure4_FUMA_gene_matrix.tsv",
                   "Figure4_evidence_chains.tsv",
                   "Figure4_approved_drug_recurrence.tsv",
                   "Figure4_directional_summary.tsv")),
       force = TRUE)
write_source(panel_a_source, "Figure4_panelA_source")
write_source(panel_b_source, "Figure4_panelB_source")
write_source(direction_summary, "Figure4_panelC_source")
save_panel_set(list(A = p4a, B = p4b, C = p4d), "Figure4",
               panel_width = 98, panel_height = 92,
               composite_height = 150, design = "AA\nBC",
               composite_stem = "Figure4_revised")
copy_composite_to_root("Figure4_revised")

# Candidate Figure 4 preview 1: all target-gene burden view
candidate_gene_evidence <- gene_evidence[
  gene_evidence$gene_grade %in% c("Gene-A", "Gene-B"), ]
candidate_drug_gene_map <- drug_gene_map[
  drug_gene_map$gene_symbol %in% candidate_gene_evidence$gene_symbol, ]
candidate_gene_drug_summary <- do.call(rbind, lapply(
  split(candidate_drug_gene_map, candidate_drug_gene_map$gene_symbol),
  function(z) {
    pair_set <- split_semicolon(z$pair_list)
    drug_rank <- aggregate(n_pairs_for_gene_drug ~ drug_name + display_drug_name,
                           z, max)
    drug_rank <- drug_rank[order(-drug_rank$n_pairs_for_gene_drug,
                                 drug_rank$drug_name), ]
    data.frame(
      gene_symbol = z$gene_symbol[1],
      n_approved_drugs = length(unique(z$drug_name)),
      drug_supported_pairs = length(pair_set),
      top_approved_drugs = paste(head(drug_rank$display_drug_name, 2),
                                 collapse = "; "),
      all_approved_drugs = paste(drug_rank$display_drug_name, collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
))
candidate_gene_source <- merge(candidate_gene_drug_summary,
                               candidate_gene_evidence,
                               by = "gene_symbol", all.x = TRUE)
candidate_gene_source$gene_grade_rank <- ifelse(
  candidate_gene_source$gene_grade == "Gene-A", 1, 0
)
candidate_gene_source <- candidate_gene_source[order(
  -candidate_gene_source$drug_supported_pairs,
  -candidate_gene_source$n_approved_drugs,
  -candidate_gene_source$n_pair_gene_records,
  -candidate_gene_source$gene_grade_rank,
  candidate_gene_source$gene_symbol
), ]
candidate_gene_source$gene_label <- factor(
  candidate_gene_source$gene_symbol,
  levels = rev(candidate_gene_source$gene_symbol)
)
candidate_gene_source$top_drug_label <- wrap_text(
  candidate_gene_source$top_approved_drugs, 26
)
candidate_gene_source$gene_grade <- factor(
  candidate_gene_source$gene_grade,
  levels = c("Gene-A", "Gene-B")
)
candidate_n_genes <- length(unique(candidate_gene_source$gene_symbol))
candidate_n_edges <- nrow(candidate_drug_gene_map)
gene_burden_xmax <- max(candidate_gene_source$n_approved_drugs, na.rm = TRUE)
journal_pal <- c(
  gene_a = unname(figure4_gene_colours["Gene-A"]),
  gene_b = unname(figure4_gene_colours["Gene-B"]),
  link_a = "#00A087",
  link_b = "#E64B35",
  drug_fill = "#E8EEF7",
  drug_outline = "#3C5488",
  deep_blue = "#3C5488",
  grey_text = "#4C555B"
)

p4_candidate_gene_burden <- ggplot(
  candidate_gene_source,
  aes(n_approved_drugs, gene_label)
) +
  geom_col(aes(fill = gene_grade), width = 0.62, alpha = 0.88) +
  geom_point(aes(size = drug_supported_pairs),
             shape = 21, fill = "white",
             colour = unname(journal_pal["deep_blue"]),
             stroke = 0.25, alpha = 0.92) +
  geom_text(aes(label = n_approved_drugs),
            hjust = -0.25, size = 2.0, colour = COL["graphite"],
            fontface = "bold") +
  geom_text(aes(x = gene_burden_xmax + 4.0, label = top_drug_label),
            hjust = 0, size = 1.85, colour = COL["graphite"],
            lineheight = 0.82) +
  annotate("text", x = gene_burden_xmax + 4.0,
           y = nrow(candidate_gene_source) + 1.15,
           label = "Representative approved drugs",
           hjust = 0, size = 2.6, fontface = "bold",
           colour = COL["graphite"]) +
  scale_fill_manual(values = c("Gene-A" = unname(journal_pal["gene_a"]),
                               "Gene-B" = unname(journal_pal["gene_b"]))) +
  scale_size_continuous(range = c(1.2, 4.4),
                        breaks = pretty(candidate_gene_source$drug_supported_pairs,
                                        n = 4)) +
  scale_x_continuous(limits = c(0, gene_burden_xmax + 22),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Candidate 1: all approved-drug target genes",
       subtitle = paste0("All ", candidate_n_genes,
                         " Gene-A/Gene-B target genes are shown; bars count approved-drug annotations and point size shows disease-pair support"),
       x = "Approved-drug annotations per target gene",
       y = NULL, fill = NULL, size = "Drug-pair\nsupport") +
  theme_pub(6.1) +
  theme(
    axis.text.y = element_text(size = 5.5),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    plot.margin = margin(6, 55, 6, 8)
  )

write_source(candidate_gene_source, "Figure4_candidate_gene_burden_source")
save_artifact(p4_candidate_gene_burden, "Figure4_candidate_gene_burden",
              ".", 183, 190, "Figure4_candidate_gene_burden",
              "", "candidate")

# Candidate Figure 4 preview 2: non-column arc network
arc_edge_source <- merge(candidate_drug_gene_map, candidate_gene_source,
                         by = "gene_symbol", all.x = TRUE)
arc_edge_source <- arc_edge_source[order(
  -arc_edge_source$n_pairs_for_gene_drug,
  arc_edge_source$drug_name,
  arc_edge_source$gene_symbol
), ]

arc_gene_nodes <- candidate_gene_source
arc_gene_nodes$gene_index <- seq_len(nrow(arc_gene_nodes))
arc_gene_nodes$theta <- seq(pi / 2, pi / 2 - 2 * pi,
                            length.out = nrow(arc_gene_nodes) + 1)[
                              seq_len(nrow(arc_gene_nodes))]
arc_gene_nodes$x <- 1.05 * cos(arc_gene_nodes$theta)
arc_gene_nodes$y <- 1.05 * sin(arc_gene_nodes$theta)
arc_gene_nodes$label_x <- 1.25 * cos(arc_gene_nodes$theta)
arc_gene_nodes$label_y <- 1.25 * sin(arc_gene_nodes$theta)
arc_gene_nodes$hjust <- ifelse(arc_gene_nodes$label_x >= 0, 0, 1)

arc_drug_rank <- do.call(rbind, lapply(split(arc_edge_source,
                                             arc_edge_source$drug_name),
                                       function(z) {
  data.frame(
    drug_name = z$drug_name[1],
    display_drug_name = z$display_drug_name[1],
    max_pair_support = max(z$n_pairs_for_gene_drug, na.rm = TRUE),
    n_connected_genes = length(unique(z$gene_symbol)),
    stringsAsFactors = FALSE
  )
}))
arc_must_keep <- unique(unlist(lapply(split(arc_edge_source,
                                            arc_edge_source$gene_symbol),
                                      function(z) {
  z <- z[order(-z$n_pairs_for_gene_drug, z$drug_name), ]
  head(z$drug_name, 1)
})))
arc_drug_rank$must_keep <- arc_drug_rank$drug_name %in% arc_must_keep
arc_drug_rank <- arc_drug_rank[order(!arc_drug_rank$must_keep,
                                     -arc_drug_rank$max_pair_support,
                                     -arc_drug_rank$n_connected_genes,
                                     arc_drug_rank$drug_name), ]
arc_selected_drugs <- unique(c(arc_must_keep,
                               head(arc_drug_rank$drug_name, 38)))
arc_selected_drugs <- arc_selected_drugs[
  arc_selected_drugs %in% arc_drug_rank$drug_name
]
arc_drug_nodes <- arc_drug_rank[
  match(arc_selected_drugs, arc_drug_rank$drug_name), ]
arc_drug_nodes$drug_index <- seq_len(nrow(arc_drug_nodes))
arc_drug_nodes$theta <- seq(pi / 2 + pi / 36,
                            pi / 2 + pi / 36 - 2 * pi,
                            length.out = nrow(arc_drug_nodes) + 1)[
                              seq_len(nrow(arc_drug_nodes))]
arc_drug_nodes$x <- 2.0 * cos(arc_drug_nodes$theta)
arc_drug_nodes$y <- 2.0 * sin(arc_drug_nodes$theta)
arc_drug_nodes$label_x <- 2.24 * cos(arc_drug_nodes$theta)
arc_drug_nodes$label_y <- 2.24 * sin(arc_drug_nodes$theta)
arc_drug_nodes$hjust <- ifelse(arc_drug_nodes$label_x >= 0, 0, 1)
arc_drug_nodes$label_this <- arc_drug_nodes$max_pair_support >= 7 |
  arc_drug_nodes$n_connected_genes > 1
arc_drug_nodes$drug_label <- arc_drug_nodes$display_drug_name

arc_gene_label_nodes <- arc_gene_nodes[order(
  -arc_gene_nodes$n_pair_gene_records,
  -arc_gene_nodes$n_approved_drugs,
  arc_gene_nodes$gene_symbol
), ]
arc_gene_label_nodes <- head(arc_gene_label_nodes, 24)

arc_drug_label_nodes <- arc_drug_nodes[arc_drug_nodes$label_this, ]
arc_drug_label_nodes <- arc_drug_label_nodes[order(
  -arc_drug_label_nodes$max_pair_support,
  -arc_drug_label_nodes$n_connected_genes,
  arc_drug_label_nodes$drug_name
), ]
arc_drug_label_nodes <- head(arc_drug_label_nodes, 16)

arc_edges_plot <- arc_edge_source[
  arc_edge_source$drug_name %in% arc_selected_drugs, ]
arc_edges_plot <- merge(
  arc_edges_plot,
  arc_gene_nodes[, c("gene_symbol", "x", "y")],
  by = "gene_symbol"
)
arc_edges_plot <- merge(
  arc_edges_plot,
  arc_drug_nodes[, c("drug_name", "x", "y")],
  by = "drug_name",
  suffixes = c("_gene", "_drug")
)
arc_edges_plot$shown_in_candidate_arc <- TRUE
arc_edge_source$shown_in_candidate_arc <- paste(
  arc_edge_source$drug_name, arc_edge_source$gene_symbol, sep = "__"
) %in% paste(arc_edges_plot$drug_name, arc_edges_plot$gene_symbol,
             sep = "__")
arc_edge_source$shown_in_candidate_arc[is.na(
  arc_edge_source$shown_in_candidate_arc
)] <- FALSE

p4_candidate_arc_network <- ggplot() +
  geom_segment(
    data = arc_edges_plot,
    aes(x = x_drug, y = y_drug, xend = x_gene, yend = y_gene,
        linewidth = n_pairs_for_gene_drug, colour = gene_grade),
    alpha = 0.26, lineend = "round"
  ) +
  geom_point(
    data = arc_drug_nodes,
    aes(x, y, size = max_pair_support),
    shape = 21, fill = unname(journal_pal["drug_fill"]),
    colour = unname(journal_pal["drug_outline"]),
    stroke = 0.28, alpha = 0.96
  ) +
  geom_point(
    data = arc_gene_nodes,
    aes(x, y, size = n_pair_gene_records, fill = gene_grade),
    shape = 21, colour = "white", stroke = 0.35, alpha = 0.98
  ) +
  geom_text(
    data = arc_gene_label_nodes,
    aes(label_x, label_y, label = gene_symbol, hjust = hjust),
    size = 1.68, colour = COL["graphite"], fontface = "bold",
    check_overlap = TRUE
  ) +
  geom_text(
    data = arc_drug_label_nodes,
    aes(label_x, label_y, label = drug_label, hjust = hjust),
    size = 1.30, colour = unname(journal_pal["grey_text"]),
    lineheight = 0.76, check_overlap = TRUE
  ) +
  annotate("text", x = 0, y = 0,
           label = "Approved-drug\ntarget genes",
           size = 3.0, fontface = "bold",
           colour = unname(journal_pal["deep_blue"]),
           lineheight = 0.85) +
  scale_fill_manual(values = c("Gene-A" = unname(journal_pal["gene_a"]),
                               "Gene-B" = unname(journal_pal["gene_b"]))) +
  scale_colour_manual(values = c("Gene-A" = unname(journal_pal["link_a"]),
                                 "Gene-B" = unname(journal_pal["link_b"])),
                      guide = "none") +
  scale_size_continuous(range = c(1.2, 5.2), guide = "none") +
  scale_linewidth_continuous(range = c(0.12, 0.72), guide = "none") +
  coord_fixed(xlim = c(-2.75, 2.75), ylim = c(-2.55, 2.65),
              clip = "off") +
  labs(title = "Candidate 2: arc network of target genes and representative drugs",
       subtitle = paste0("All ", candidate_n_genes,
                         " target genes are shown; labelled drugs are highest-support representatives, with all ",
                         candidate_n_edges,
                         " drug-target edges retained in source data"),
       x = NULL, y = NULL, fill = NULL) +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 8.2,
                              colour = COL["graphite"], hjust = 0),
    plot.subtitle = element_text(size = 5.7, colour = COL["grey"], hjust = 0),
    legend.position = "bottom",
    legend.text = element_text(size = 5.4),
    plot.margin = margin(6, 18, 6, 18)
  )

write_source(arc_edge_source, "Figure4_candidate_arc_network_source")
save_artifact(p4_candidate_arc_network, "Figure4_candidate_arc_network",
              ".", 183, 170, "Figure4_candidate_arc_network",
              "", "candidate")

finish_qa()
