PROJECT_ROOT <- if (dir.exists(file.path(getwd(), "figure_final"))) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
}
source(file.path(PROJECT_ROOT, "figure_final", "scripts", "common.R"))

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
pair_support_breaks <- sort(unique(as.integer(c(
  min(panel_b_source$n_pairs_for_gene_drug, na.rm = TRUE),
  round(stats::median(panel_b_source$n_pairs_for_gene_drug, na.rm = TRUE)),
  max(panel_b_source$n_pairs_for_gene_drug, na.rm = TRUE)
))))

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
  scale_fill_manual(
    values = figure4_gene_colours,
    name = "Gene class",
    guide = guide_legend(
      order = 1,
      override.aes = list(size = 4.8, shape = 21, colour = "white",
                          alpha = 1, stroke = 0.3)
    )
  ) +
  scale_size_continuous(range = c(1.5, 4.8),
                        breaks = pair_support_breaks,
                        name = "Pairs") +
  guides(
    size = guide_legend(
      order = 2,
      override.aes = list(shape = 21, fill = COL["graphite"],
                          colour = COL["graphite"], alpha = 1, stroke = 0.2)
    )
  ) +
  labs(title = "Approved drug-target recurrence matrix",
       subtitle = "Point size shows disease-pair support for each drug-target edge",
       x = NULL, y = NULL) +
  theme_pub(5.4) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.justification = "center",
    legend.title = element_text(size = 6.3, colour = COL["graphite"]),
    legend.text = element_text(size = 6.0, colour = COL["graphite"]),
    legend.key.size = unit(4.5, "mm"),
    legend.spacing.x = unit(2.0, "mm"),
    legend.spacing.y = unit(0.4, "mm"),
    legend.box.spacing = unit(0.8, "mm"),
    legend.margin = margin(0, 0, 0, 0)
  )
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

write_source(lcv_source, "Figure4C_directional_LCV_source")
write_source(mr_source, "Figure4C_directional_MR_source")
p4a_arc <- p4_candidate_arc_network +
  labs(
    title = "Genetically prioritized drug-target arc network",
    subtitle = paste(
      "Approved drugs connect to prioritized Gene-A/Gene-B target genes;",
      "line width reflects disease-pair support"
    )
  ) +
  theme(
    legend.position = "none"
  )

p4c_wrapped <- wrap_elements(full = directional_preview)

figure4_top <- (p4a_arc | p4b) +
  plot_layout(widths = c(1.28, 1.00))

figure4_assembled <- figure4_top / p4c_wrapped +
  plot_layout(heights = c(1.00, 0.58)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 10,
                            colour = COL["graphite"]),
    plot.margin = margin(4, 4, 4, 4)
  )

save_artifact(p4a_arc, "Figure4A",
              file.path("main", "panels"),
              183, 132, "Figure4", "A", "panel", dpi = 600)
save_artifact(directional_preview, "Figure4C",
              file.path("main", "panels"),
              183, 69, "Figure4", "C", "panel", dpi = 600)
save_artifact(figure4_assembled, "Figure4_revised",
              file.path("main", "composites"),
              183, 198, "Figure4", "ABC", "composite", dpi = 600)
copy_composite_to_root("Figure4_revised")

cat("Assembled Figure 4 panels: top row arc-network A + matrix B, bottom row C\n")

finish_qa()
