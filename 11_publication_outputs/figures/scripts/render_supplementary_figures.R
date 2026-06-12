PROJECT_ROOT <- if (dir.exists(file.path(getwd(), "figure_final"))) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
}
source(file.path(PROJECT_ROOT, "figure_final", "scripts", "common.R"))

tab <- function(n) file.path(RESULT, "Table", n)
save_supp <- function(panels, id, design = NULL, h = 150) {
  save_panel_set(panels, id, main = FALSE, panel_width = 96,
                 panel_height = 88, composite_width = 183,
                 composite_height = h, design = design,
                 composite_stem = paste0(id, "_revised"))
}
save_supp_single <- function(plot, id, h = 95) {
  save_artifact(plot, paste0(id, "A"),
                file.path("supplementary", "panels"),
                183, h, id, "A", "panel")
  save_artifact(plot, paste0(id, "_revised"),
                file.path("supplementary", "composites"),
                183, h, id, "", "composite")
}
cleanup_supplementary_render <- function() {
  unlink(Sys.glob(file.path(FIG, "supplementary", "composites",
                            "SupplementaryFigure*_revised.*")))
  unlink(Sys.glob(file.path(FIG, "supplementary", "panels",
                            "SupplementaryFigure*[A-Z].*")))
  unlink(Sys.glob(file.path(FIG, "source_data",
                            "SupplementaryFigure*.tsv")))
}
cleanup_supplementary_render()

pair_axis_label <- function(pair_id) {
  paste(ifelse(pair_anchor(pair_id) == "Nasal polyps", "NP", "AR"),
        pair_partner(pair_id), sep = " - ")
}

# Supplementary Figure 1: source and power audit
meta <- read_tsv(file.path(ROOT, "method_table",
  "Supplementary_Table_1_GWAS_source_and_phenotype_information.tsv"))
meta$power <- ifelse(grepl("retained|adequate", meta$ldsc_h2_power_class,
                          ignore.case = TRUE), "Adequate/retained",
                     ifelse(grepl("borderline", meta$ldsc_h2_power_class,
                                  ignore.case = TRUE), "Borderline", "Low/not retained"))
power_count <- as.data.frame(table(meta$power))
names(power_count) <- c("power", "n")
n_s1 <- data.frame(
  x = c(1, 1, 2.5, 4, 4), y = c(2.7, 1.3, 2, 2.7, 1.3),
  label = c("4 nasal\nphenotypes", "57 comparator\ntraits",
            "LDSC h2 and\npower audit", "2 retained\nnasal anchors",
            "31 final\npairs"),
  fill = c(colv("pale"), colv("pale"), "#DCE8ED", "#DCE8F2", "#CFE7E3")
)
e_s1 <- data.frame(x = c(1.35, 1.35, 2.85, 2.85),
                   y = c(2.7, 1.3, 2.1, 1.9),
                   xend = c(2.1, 2.1, 3.55, 3.55),
                   yend = c(2.2, 1.8, 2.7, 1.3))
s1a <- flow_plot(n_s1, e_s1, "GWAS universe and retention")
s1b <- ggplot(power_count, aes(n, reorder(power, n), fill = power)) +
  geom_col(width = 0.65) + geom_text(aes(label = n), hjust = -0.15) +
  scale_fill_manual(values = c("Adequate/retained" = colv("teal"),
                               "Borderline" = colv("amber"),
                               "Low/not retained" = colv("grey"))) +
  scale_x_continuous(expand = expansion(mult = c(0, .15))) +
  labs(title = "Trait-level heritability power classes",
       x = "Traits", y = NULL, fill = NULL) +
  theme_pub(7) + theme(legend.position = "none")
ret <- aggregate(cbind(n_final_downstream_pairs,
                       n_ldsc_rg_pairs_bh05,
                       n_lava_significant_pair_regions) ~ trait_label,
                 meta[meta$appears_in_final_31_pairs %in% c(TRUE, "True"), ],
                 sum, na.rm = TRUE)
ret <- head(ret[order(-ret$n_final_downstream_pairs), ], 14)
s1c <- ggplot(ret, aes(n_final_downstream_pairs,
                       reorder(short_trait(trait_label),
                               n_final_downstream_pairs))) +
  geom_col(fill = colv("blue"), width = .65) +
  labs(title = "Traits represented in final pairs",
       x = "Final disease pairs", y = NULL) + theme_pub(6.5)
write_source(power_count, "SupplementaryFigure1_power_classes")
write_source(ret, "SupplementaryFigure1_retained_traits")
save_supp(list(A = s1a, B = s1b, C = s1c), "SupplementaryFigure1",
          design = "AB\nCC", h = 155)

# Supplementary Figure 2: LDSC and LAVA
ldsc <- read_tsv(file.path(END, "ldsc", "ldsc_rg_bubble_input.tsv"))
pair <- read_tsv(tab("Table_R2_pair_screening_and_locus_summary.tsv"))
ldsc <- ldsc[ldsc$pair_id %in% pair$pair_id, ]
ldsc$partner <- pair_partner(ldsc$pair_id)
ldsc$pair_axis <- pair_axis_label(ldsc$pair_id)
ldsc$lo <- ldsc$x - 1.96 * ldsc$size
ldsc$hi <- ldsc$x + 1.96 * ldsc$size
ldsc$sig <- ldsc$ldsc_bh_q < .05
s2a <- ggplot(ldsc, aes(x, reorder(pair_axis, x), colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = colv("grey")) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .4) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("TRUE" = colv("blue"),
                                 "FALSE" = colv("grey")),
                      labels = c("TRUE" = "Significant",
                                 "FALSE" = "Not significant")) +
  labs(title = "LDSC genetic correlation", subtitle = "Point estimate and 95% CI",
       x = expression(r[g]), y = NULL, colour = "BH-FDR < 0.05") +
  theme_pub(6) + theme(legend.position = "bottom")
pair$partner <- pair_partner(pair$pair_id)
pair$pair_axis <- pair_axis_label(pair$pair_id)
s2b <- ggplot(pair, aes(lava_n_significant_loci,
                        reorder(pair_axis, lava_n_significant_loci),
                        fill = pair_anchor(pair_id))) +
  geom_col(width = .65) +
  scale_fill_manual(values = pal_anchor) +
  labs(title = "Significant LAVA regions",
       x = "Regions", y = NULL, fill = NULL) +
  theme_pub(6) + theme(legend.position = "bottom")
pair$route <- ifelse(tolower(pair$ldsc_sig_bh05) == "true" &
                              pair$lava_n_significant_loci > 0, "Both",
                     ifelse(tolower(pair$ldsc_sig_bh05) == "true",
                            "LDSC only", "LAVA only"))
route <- as.data.frame(table(pair$route)); names(route) <- c("route", "n")
s2c <- ggplot(route, aes(route, n, fill = route)) +
  geom_col(width = .62) + geom_text(aes(label = n), vjust = -0.3,
                                    fontface = "bold") +
  scale_fill_manual(values = c("Both" = colv("teal"),
                               "LDSC only" = colv("blue"),
                               "LAVA only" = colv("amber"))) +
  labs(title = "Final-pair inclusion routes", x = NULL, y = "Pairs") +
  theme_pub(7) + theme(legend.position = "none")
write_source(ldsc, "SupplementaryFigure2_LDSC")
write_source(pair, "SupplementaryFigure2_LAVA_pair_summary")
write_source(route, "SupplementaryFigure2_inclusion_routes")
save_supp(list(A = s2a, B = s2b, C = s2c), "SupplementaryFigure2",
          design = "AB\nCC", h = 170)

# Supplementary Figure 3: cross-trait methods and colocalization audit
locus_master <- read_tsv(file.path(END, "integrated", "locus_master_table.tsv"))
mc <- aggregate(SNP ~ pair_id + method, locus_master, length)
names(mc)[3] <- "n"
mc$partner <- pair_partner(mc$pair_id)
mc$pair_axis <- pair_axis_label(mc$pair_id)
mc$method <- gsub("_genomewide", "", mc$method)
s3a <- ggplot(mc, aes(method, pair_axis, fill = log10(n + 1))) +
  geom_tile(colour = "white", linewidth = .2) +
  scale_fill_gradient(low = colv("pale"), high = colv("blue")) +
  labs(title = "Cross-trait SNP evidence by method",
       x = NULL, y = NULL, fill = expression(log[10](n+1))) +
  theme_pub(5.7) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "bottom")
n_s3 <- data.frame(
  x = c(1, 2.5, 2.5, 4.2), y = c(2, 2.7, 1.3, 2),
  label = c("Candidate\npair-loci", "Coloc ABF\n182 positive",
            "SuSiE-coloc\n276 entries\n285 analyses\n76 positive",
            "68 Locus-A\n113 Locus-B"),
  fill = c(colv("pale"), "#DCE8F2", "#E5E1F2", "#D8E9DF")
)
e_s3 <- data.frame(x = c(1.35, 1.35, 2.95, 2.95),
                   y = c(2.1, 1.9, 2.7, 1.3),
                   xend = c(2.05, 2.05, 3.75, 3.75),
                   yend = c(2.7, 1.3, 2.1, 1.9))
s3b <- flow_plot(n_s3, e_s3, "Colocalization audit",
                 "Parallel branches converge only at grading")
write_source(mc, "SupplementaryFigure3_method_counts")
save_supp(list(A = s3a, B = s3b), "SupplementaryFigure3", h = 105)

# Supplementary Figure 4: full Figure 4 evidence archive
bool_true <- function(x) x %in% c(TRUE, "True", "TRUE", "true", "1", 1)
safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}
split_semicolon <- function(x) {
  x <- unlist(strsplit(paste(x, collapse = ";"), ";", fixed = TRUE))
  sort(unique(x[nzchar(x)]))
}

twas_chain <- read_tsv(file.path(
  END, "final_gtexv8_twas_total_domain_fuma_enrichment",
  "final_gene_tables", "final_twas_supported_pair_gene_records.tsv"
))
lstat4 <- read_tsv(file.path(END, "final_coloc_susie_export",
                             "final_all_31_pairs_locus_status.tsv"))
rec4 <- read_tsv(tab("Table_R4_top_recurrent_prioritized_genes.tsv"))
twas_chain <- merge(
  twas_chain,
  lstat4[, c("pair_id", "locus_id", "coloc_pph4", "susie_status")],
  by = c("pair_id", "locus_id"),
  all.x = TRUE
)
twas_chain <- merge(
  twas_chain,
  rec4[, c("gene_symbol", "n_pair_gene_records",
           "n_phenotype_domains", "strict_positive_record_count")],
  by = "gene_symbol",
  all.x = TRUE
)
twas_chain$n_pair_gene_records[is.na(twas_chain$n_pair_gene_records)] <- 1
twas_chain$n_phenotype_domains[is.na(twas_chain$n_phenotype_domains)] <- 1
twas_chain$strict_positive_record_count[
  is.na(twas_chain$strict_positive_record_count)
] <- 0
short_locus <- function(x) paste0("L", sub("^0+", "",
                                           sub("^unified_locus_", "", x)))
evidence_cols <- c("Locus-A/B", "Coloc/SuSiE", "MAGMA",
                   "TWAS FDR", "Multi-tissue TWAS", "Gene-A/B",
                   "Strict positive", "Gene recurrence")
chain_summary <- do.call(rbind, lapply(seq_len(nrow(twas_chain)), function(i) {
  z <- twas_chain[i, ]
  pair_short_label <- paste(
    ifelse(pair_anchor(z$pair_id) == "Allergic rhinitis", "AR", "NP"),
    pair_partner(z$pair_id)
  )
  values <- c(
    "Locus-A/B" = ifelse(z$locus_grade == "A", 2, 1),
    "Coloc/SuSiE" = ifelse(as.numeric(z$coloc_pph4) >= 0.80 &
                             z$susie_status == "ok", 2,
                           ifelse(as.numeric(z$coloc_pph4) >= 0.50, 1, 0)),
    "MAGMA" = ifelse(as.numeric(z$magma_p) <= 5e-8, 2, 1),
    "TWAS FDR" = ifelse(as.numeric(z$best_twas_fdr) <= 0.01, 2, 1),
    "Multi-tissue TWAS" = ifelse(as.numeric(z$n_significant_tissues) >= 3, 2, 1),
    "Gene-A/B" = ifelse(z$gene_grade == "Gene-A", 2, 1),
    "Strict positive" = ifelse(bool_true(z$strict_positive), 2, 0),
    "Gene recurrence" = ifelse(as.numeric(z$n_pair_gene_records) >= 5, 2,
                               ifelse(as.numeric(z$n_pair_gene_records) >= 2, 1, 0))
  )
  data.frame(
    row = paste(pair_short_label, short_locus(z$locus_id),
                z$gene_symbol, sep = " | "),
    pair_id = z$pair_id,
    pair_short = pair_short_label,
    gene = z$gene_symbol,
    locus_id = z$locus_id,
    locus_short = short_locus(z$locus_id),
    locus_grade = z$locus_grade,
    gene_grade = z$gene_grade,
    magma_p = as.numeric(z$magma_p),
    best_TWAS_fdr = as.numeric(z$best_twas_fdr),
    best_TWAS_tissue = z$best_twas_tissue,
    n_significant_tissues = as.numeric(z$n_significant_tissues),
    coloc_pph4 = as.numeric(z$coloc_pph4),
    susie_status = z$susie_status,
    strict_positive = bool_true(z$strict_positive),
    n_pair_gene_records = as.numeric(z$n_pair_gene_records),
    n_phenotype_domains = as.numeric(z$n_phenotype_domains),
    support_score = sum(values),
    n_evidence_types = sum(values > 0),
    t(values),
    check.names = FALSE
  )
}))
chain_summary <- chain_summary[order(-chain_summary$support_score,
                                     -chain_summary$n_pair_gene_records,
                                     -chain_summary$n_significant_tissues,
                                     chain_summary$best_TWAS_fdr,
                                     chain_summary$row), ]
chain_full <- do.call(rbind, lapply(seq_len(nrow(chain_summary)), function(i) {
  z <- chain_summary[i, ]
  data.frame(
    row = z$row,
    pair_id = z$pair_id,
    pair_short = z$pair_short,
    gene = z$gene,
    locus_id = z$locus_id,
    locus_short = z$locus_short,
    locus_grade = z$locus_grade,
    gene_grade = z$gene_grade,
    magma_p = z$magma_p,
    best_TWAS_fdr = z$best_TWAS_fdr,
    best_TWAS_tissue = z$best_TWAS_tissue,
    n_significant_tissues = z$n_significant_tissues,
    coloc_pph4 = z$coloc_pph4,
    susie_status = z$susie_status,
    strict_positive = z$strict_positive,
    n_pair_gene_records = z$n_pair_gene_records,
    n_phenotype_domains = z$n_phenotype_domains,
    support_score = z$support_score,
    evidence = factor(evidence_cols, levels = evidence_cols),
    value = as.numeric(z[, evidence_cols]),
    stringsAsFactors = FALSE
  )
}))
chain_full$support <- factor(ifelse(chain_full$value >= 2, "Strong support",
                             ifelse(chain_full$value == 1, "Moderate support",
                                    "Absent")),
                             levels = c("Absent", "Moderate support",
                                        "Strong support"))
chain_plot <- chain_full[chain_full$row %in% head(chain_summary$row, 32), ]
chain_plot$row <- factor(chain_plot$row, levels = rev(unique(head(chain_summary$row, 32))))
s4a <- ggplot(chain_plot, aes(evidence, row, fill = support)) +
  geom_tile(colour = "white", linewidth = .18) +
  scale_fill_manual(values = c("Absent" = "#F4F7F8",
                               "Moderate support" = colv("amber"),
                               "Strong support" = colv("teal")),
                    drop = FALSE) +
  labs(title = "Non-FUMA locus-to-gene evidence chains",
       subtitle = "Rows ranked by locus, MAGMA, TWAS and recurrence support; source data retains all TWAS-supported pair-gene records",
       x = NULL, y = NULL, fill = NULL) +
  theme_pub(4.8) +
  theme(axis.text.x = element_text(angle = 38, hjust = 1),
        axis.text.y = element_text(size = 3.4),
        legend.position = "bottom",
        panel.grid = element_blank())

drug4 <- read_tsv(file.path(END, "network_pharmacology_res",
                            "drug_recurrence_summary.tsv"))
drug4 <- drug4[drug4$approved_any %in% c(TRUE, "True"), ]
drug_full <- do.call(rbind, lapply(split(drug4, drug4$drug_name), function(z) {
  pair_set <- split_semicolon(z$pair_list)
  gene_set <- split_semicolon(z$target_gene_list)
  data.frame(
    drug_name = z$drug_name[1],
    display_drug_name = tools::toTitleCase(tolower(z$drug_name[1])),
    n_pairs = length(pair_set),
    pair_list = paste(pair_set, collapse = ";"),
    n_target_genes_total = length(gene_set),
    target_gene_list = paste(gene_set, collapse = ";"),
    max_target_genes_per_pair = max(suppressWarnings(
      as.numeric(z$max_target_genes_per_pair)), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
drug_full <- drug_full[order(-drug_full$n_pairs,
                             -drug_full$n_target_genes_total,
                             drug_full$drug_name), ]
drug_plot <- head(drug_full, 30)
drug_plot$display_drug_name <- factor(drug_plot$display_drug_name,
                                      levels = rev(drug_plot$display_drug_name))
s4b <- ggplot(drug_plot, aes(n_pairs, display_drug_name)) +
  geom_segment(aes(x = 0, xend = n_pairs, yend = display_drug_name),
               colour = "#D7DEE2", linewidth = 1.1) +
  geom_point(aes(size = n_target_genes_total), colour = colv("violet"),
             alpha = .88) +
  scale_size_continuous(range = c(1.5, 4.0),
                        breaks = sort(unique(drug_plot$n_target_genes_total))) +
  scale_x_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Full approved-drug target annotation",
       subtitle = "Top recurrent approved drugs displayed; source data retains all approved annotations",
       x = "Number of disease pairs", y = NULL, size = "Target genes") +
  theme_pub(4.9) + theme(legend.position = "bottom",
                         axis.text.y = element_text(size = 3.7))

lcv4 <- read_tsv(tab("Table_R14_LCVMR_final31_summary.tsv"))
lcv4$partner <- pair_partner(lcv4$pair_id)
lcv4$pair_axis <- pair_axis_label(lcv4$pair_id)
lcv4$lo <- lcv4$LCV_gcp - 1.96 * lcv4$LCV_gcp_se
lcv4$hi <- lcv4$LCV_gcp + 1.96 * lcv4$LCV_gcp_se
lcv4$status <- ifelse(lcv4$LCV_gcp_p < .05 / 31, "Bonferroni",
                      ifelse(lcv4$LCV_gcp_p < .05, "Nominal",
                             "Not significant"))
s4c <- ggplot(lcv4, aes(LCV_gcp, reorder(pair_axis, LCV_gcp),
                        colour = status)) +
  geom_vline(xintercept = 0, linetype = 2, colour = colv("grey")) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .28) +
  geom_point(size = 1.35) +
  scale_colour_manual(values = c("Bonferroni" = colv("teal"),
                                 "Nominal" = colv("amber"),
                                 "Not significant" = colv("grey"))) +
  labs(title = "Full LCV directional evidence",
       x = "Genetic causality proportion", y = NULL, colour = NULL) +
  theme_pub(4.8) + theme(legend.position = "bottom",
                         axis.text.y = element_text(size = 3.8))

mr4 <- read_tsv(file.path(END, "LCVMR", "MR", "mr_results_all.tsv"))
ivw4 <- mr4[mr4$method == "IVW", ]
qc4 <- read_tsv(tab("Table_R15_MR_direction_QC_final31.tsv"))
mr_full <- merge(qc4, ivw4[, c("direction", "beta", "se", "pval", "nsnp")],
                 by = "direction", all.x = TRUE)
mr_full$label <- paste(short_trait(mr_full$exposure), "to",
                       short_trait(mr_full$outcome))
mr_full$status2 <- ifelse(mr_full$status != "eligible", "Insufficient",
                          ifelse(mr_full$pval < .05 / 62, "Bonferroni",
                                 ifelse(mr_full$pval < .05, "Nominal",
                                        "Not significant")))
mr_plot <- mr_full[order(mr_full$pval), ]
mr_plot <- head(mr_plot, 32)
mr_plot$label <- factor(mr_plot$label, levels = rev(mr_plot$label))
s4d <- ggplot(mr_plot, aes(-log10(pval), label, colour = status2)) +
  geom_point(size = 1.4, na.rm = TRUE) +
  geom_vline(xintercept = -log10(.05 / 62), linetype = 2,
             colour = colv("red")) +
  scale_colour_manual(values = c("Bonferroni" = colv("blue"),
                                 "Nominal" = colv("amber"),
                                 "Not significant" = colv("grey"),
                                 "Insufficient" = colv("pale"))) +
  labs(title = "Full MR directional evidence",
       subtitle = "Top IVW directions displayed; all directions retained in source data",
       x = expression(-log[10](P)), y = NULL, colour = NULL) +
  theme_pub(4.7) + theme(legend.position = "bottom",
                         axis.text.y = element_text(size = 3.6))

path_breadth <- read_tsv(tab("Table_R5_pathway_term_count_summary.tsv"))
path_breadth$phenotype_domain[is.na(path_breadth$phenotype_domain)] <- ""
path_breadth$group <- ifelse(path_breadth$analysis_type == "primary", "Global primary",
                      ifelse(path_breadth$analysis_type == "strict", "Global strict",
                      ifelse(path_breadth$analysis_type == "conditional", "Global conditional",
                      ifelse(path_breadth$analysis_type == "global_background_sensitivity",
                             "Background sensitivity",
                      ifelse(path_breadth$analysis_type == "conditional_MAGMA",
                             "Conditional MAGMA",
                      ifelse(path_breadth$phenotype_domain == "asthma_lower_airway",
                             "Asthma/lower airway",
                      ifelse(path_breadth$phenotype_domain == "atopic_allergic",
                             "Atopic/allergic",
                      ifelse(path_breadth$phenotype_domain == "autoimmune_IBD",
                             "Autoimmune / IBD",
                      ifelse(path_breadth$phenotype_domain == "ENT_infection",
                             "ENT / infection", "Other domains")))))))))
path_breadth <- path_breadth[path_breadth$database %in%
                               c("GO_BP", "KEGG", "Reactome"), ]
s4e <- ggplot(path_breadth, aes(as.numeric(n_FDR_lt_0.05),
                                reorder(group, as.numeric(n_FDR_lt_0.05)),
                                fill = database)) +
  geom_col(position = position_dodge2(width = .78, preserve = "single"),
           width = .65, alpha = .9) +
  scale_fill_manual(values = c("GO_BP" = colv("teal"),
                               "KEGG" = colv("amber"),
                               "Reactome" = colv("violet"))) +
  scale_x_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Pathway database breadth",
       x = "FDR-significant pathway terms", y = NULL, fill = "Database") +
  theme_pub(5.4) + theme(legend.position = "bottom")

tissue4 <- read_tsv(file.path(
  END, "final_gtexv8_twas_total_domain_fuma_enrichment",
  "tissue_summary", "tissue_level_support_counts.tsv"
))
tissue4$tissue_clean <- gsub("_", " ", tissue4$tissue)
tissue4 <- tissue4[order(-tissue4$n_supported_pair_gene_records,
                         tissue4$tissue), ]
tissue_plot <- head(tissue4, 22)
tissue_long <- rbind(
  data.frame(tissue_clean = tissue_plot$tissue_clean,
             grade = "Gene-A", n = tissue_plot$n_geneA_records),
  data.frame(tissue_clean = tissue_plot$tissue_clean,
             grade = "Gene-B", n = tissue_plot$n_geneB_records)
)
tissue_long$tissue_clean <- factor(
  tissue_long$tissue_clean,
  levels = rev(tissue_plot$tissue_clean)
)
s4f <- ggplot(tissue_long, aes(n, tissue_clean, fill = grade)) +
  geom_col(width = .68, alpha = .92) +
  scale_fill_manual(values = c("Gene-A" = colv("teal"),
                               "Gene-B" = colv("amber"))) +
  labs(title = "Non-FUMA tissue support",
       subtitle = "GTEx TWAS-supported pair-gene-tissue records",
       x = "Supported records", y = NULL, fill = NULL) +
  theme_pub(4.9) + theme(legend.position = "bottom",
                         axis.text.y = element_text(size = 3.8))

cell4 <- read_tsv(file.path(ROOT, "figure_v2_source_data", "cell_audit.tsv"))
cell4 <- cell4[cell4$run_ORA %in% c(TRUE, "True"), ]
cell_summary <- aggregate(cbind(n_tested_terms = as.numeric(n_tested_terms),
                                n_FDR05 = as.numeric(n_FDR05),
                                n_nominal_p05 = as.numeric(n_nominal_p05)) ~
                            marker_source, cell4, sum, na.rm = TRUE)
cell_summary <- cell_summary[order(-cell_summary$n_nominal_p05,
                                   -cell_summary$n_FDR05,
                                   cell_summary$marker_source), ]
s4g <- ggplot(cell_summary, aes(n_nominal_p05,
                                reorder(marker_source, n_nominal_p05))) +
  geom_col(fill = colv("pale"), colour = "white", width = .65) +
  geom_point(aes(x = n_FDR05), colour = colv("red"), size = 2.0) +
  geom_text(aes(label = paste0("FDR=", n_FDR05)), hjust = -0.12,
            size = 2.0, colour = colv("graphite")) +
  scale_x_continuous(expand = expansion(mult = c(0, .22))) +
  labs(title = "Non-FUMA cell-marker support",
       subtitle = "Bars: nominal terms; red points: FDR-significant terms",
       x = "Cell-marker terms", y = NULL) +
  theme_pub(5.6)

unlink(file.path(FIG, "source_data",
                 c("SupplementaryFigure4_regional_profiles.tsv",
                   paste0("Supplementary_Figure4_panel", LETTERS[1:7],
                          "_source.tsv"))),
       force = TRUE)
write_source(chain_full, "Supplementary_Figure4_panelA_source")
write_source(drug_full, "Supplementary_Figure4_panelB_source")
write_source(lcv4, "Supplementary_Figure4_panelC_source")
write_source(mr_full, "Supplementary_Figure4_panelD_source")
write_source(path_breadth, "Supplementary_Figure4_panelE_source")
write_source(tissue4, "Supplementary_Figure4_panelF_source")
write_source(cell_summary, "Supplementary_Figure4_panelG_source")
save_supp(list(A = s4a, B = s4b, C = s4c, D = s4d,
               E = s4e, F = s4f, G = s4g),
          "SupplementaryFigure4", design = "AA\nBB\nCD\nEF\nGG",
          h = 285)
for (ext in c("pdf", "png")) {
  src <- file.path(FIG, "supplementary", "composites",
                   paste0("SupplementaryFigure4_revised.", ext))
  file.copy(src, file.path(FIG, paste0("Supplementary_Figure4.", ext)),
            overwrite = TRUE)
  file.copy(src, file.path(FIG, "supplementary", "composites",
                           paste0("Supplementary_Figure4.", ext)),
            overwrite = TRUE)
}

# Supplementary Figure 9: MAGMA candidates
pair <- read_tsv(tab("Table_R2_pair_screening_and_locus_summary.tsv"))
pair$partner <- pair_partner(pair$pair_id)
pair$pair_axis <- pair_axis_label(pair$pair_id)
s5a <- ggplot(pair, aes(magma_sig_gene_count,
                        reorder(pair_axis, magma_sig_gene_count),
                        fill = pair_anchor(pair_id))) +
  geom_col(width = .65) +
  scale_fill_manual(values = pal_anchor) +
  labs(title = "MAGMA-positive candidate records by pair",
       x = "Bonferroni-positive genes", y = NULL, fill = NULL) +
  theme_pub(6) + theme(legend.position = "bottom")
mag <- read_tsv(file.path(END, "magma", "significant_genes.tsv"))
mr <- sort(table(mag$gene_symbol), decreasing = TRUE)
mr <- data.frame(gene = names(head(mr, 25)), n = as.integer(head(mr, 25)))
s5b <- ggplot(mr, aes(n, reorder(gene, n))) +
  geom_segment(aes(x = 0, xend = n, yend = reorder(gene, n)),
               colour = colv("pale"), linewidth = 2) +
  geom_point(colour = colv("teal"), size = 2.4) +
  labs(title = "Recurrent MAGMA candidates", x = "Disease pairs", y = NULL) +
  theme_pub(6.2)
write_source(mr, "SupplementaryFigure9_MAGMA_recurrence")
write_source(pair, "SupplementaryFigure9_MAGMA_pair_counts")
save_supp(list(A = s5a, B = s5b), "SupplementaryFigure9", h = 115)

# Supplementary Figure 11: complete TWAS tissue-by-gene support
genes <- read_tsv(file.path(END,
  "final_gtexv8_twas_total_domain_fuma_enrichment", "final_gene_tables",
  "final_twas_supported_pair_gene_records.tsv"))
tl <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
  z <- genes[i, ]
  ts <- trimws(strsplit(z$significant_tissue_list, ",", fixed = TRUE)[[1]])
  ts <- ts[nzchar(ts)]
  data.frame(gene = z$gene_symbol, tissue = ts,
             grade = z$gene_grade, strict = z$strict_positive,
             pair_id = z$pair_id, stringsAsFactors = FALSE)
}))
tl <- unique(tl)
tl$tissue_clean <- gsub("_", " ", tl$tissue)
cell_key <- aggregate(pair_id ~ gene + tissue + tissue_clean, tl, length)
names(cell_key)[4] <- "n_records"
cell_grade <- aggregate(grade ~ gene + tissue + tissue_clean, tl,
                        function(x) paste(sort(unique(x)), collapse = " + "))
twas_cell <- merge(cell_key, cell_grade,
                   by = c("gene", "tissue", "tissue_clean"), all.x = TRUE)
twas_cell$grade <- ifelse(twas_cell$grade == "Gene-A + Gene-B",
                          "Both grades", twas_cell$grade)
gf <- sort(tapply(twas_cell$n_records, twas_cell$gene, sum), decreasing = TRUE)
tf <- sort(tapply(twas_cell$n_records, twas_cell$tissue_clean, sum),
           decreasing = TRUE)
twas_cell$gene <- factor(twas_cell$gene, levels = rev(names(gf)))
twas_cell$tissue_clean <- factor(twas_cell$tissue_clean, levels = names(tf))
tissue_bar <- data.frame(tissue_clean = factor(names(tf), levels = names(tf)),
                         n = as.numeric(tf))
gene_bar <- data.frame(gene = factor(names(gf), levels = rev(names(gf))),
                       n = as.numeric(gf))
s6_top <- ggplot(tissue_bar, aes(tissue_clean, n)) +
  geom_col(fill = colv("ar"), width = .72, alpha = .82) +
  labs(x = NULL, y = "Records") +
  theme_pub(4.2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 4.2),
        plot.margin = margin(2, 8, 0, 8))
s6_heat <- ggplot(twas_cell, aes(tissue_clean, gene, fill = grade)) +
  geom_tile(colour = "white", linewidth = .05, alpha = .94) +
  scale_fill_manual(values = c("Gene-A" = colv("teal"),
                               "Gene-B" = colv("amber"),
                               "Both grades" = colv("violet")),
                    drop = FALSE) +
  labs(title = "Complete TWAS tissue-by-gene support",
       subtitle = "Body: 123 prioritized genes x GTEx tissues; margins show supported record counts",
       x = NULL, y = NULL, fill = NULL) +
  theme_pub(3.8) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1,
                                   size = 3.3),
        axis.text.y = element_text(size = 3.0),
        legend.position = "bottom",
        panel.grid = element_blank())
s6_right <- ggplot(gene_bar, aes(n, gene)) +
  geom_col(fill = colv("comparator"), width = .72, alpha = .72) +
  labs(x = "Records", y = NULL) +
  theme_pub(4.2) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 4.2),
        plot.margin = margin(18, 6, 6, 0))
s6a <- (s6_top + plot_spacer()) / (s6_heat + s6_right) +
  plot_layout(widths = c(8.5, 1.25), heights = c(1.1, 8.5))
write_source(tl, "SupplementaryFigure11_TWAS_matrix")
write_source(twas_cell, "SupplementaryFigure11_TWAS_matrix_cells")
save_supp(list(A = s6a), "SupplementaryFigure11", h = 165)

# Supplementary Figure 12: strict TWAS sensitivity
ov <- read_tsv(file.path(END,
  "final_gtexv8_twas_total_domain_fuma_enrichment", "audit",
  "strict_vs_main_overlap.tsv"))
ovs <- read_tsv(file.path(END,
  "final_gtexv8_twas_total_domain_fuma_enrichment", "audit",
  "strict_vs_main_overlap_summary.tsv"))
status_label <- c(
  shared_pair_gene_record = "Shared",
  main_only_pair_gene_record = "Main only",
  strict_only_pair_gene_record = "Strict only",
  shared_unique_gene = "Shared",
  main_only_unique_gene = "Main only",
  strict_only_unique_gene = "Strict only"
)
ovs$status_label <- unname(status_label[ovs$status])
s7a <- ggplot(ovs, aes(status_label, n, fill = status_label)) +
  geom_col(width = .65) + geom_text(aes(label = n), vjust = -.25) +
  facet_wrap(~ level, scales = "free_x") +
  scale_fill_manual(values = c("Shared" = colv("teal"),
                               "Main only" = colv("blue"),
                               "Strict only" = colv("amber"))) +
  labs(title = "Main versus strict TWAS overlap",
       x = NULL, y = "Count") +
  theme_pub(6) + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                       legend.position = "none")
pair_status <- ov[ov$level == "pair_gene", ]
pair_count <- aggregate(gene_symbol ~ pair_id + status, pair_status, length)
names(pair_count)[3] <- "n"
pair_count$partner <- pair_partner(pair_count$pair_id)
pair_count$pair_axis <- pair_axis_label(pair_count$pair_id)
pair_count$status_label <- unname(status_label[pair_count$status])
pair_count$total <- ave(pair_count$n, pair_count$pair_axis, FUN = sum)
s7b <- ggplot(pair_count, aes(n, reorder(pair_axis, total),
                             fill = status_label)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("Shared" = colv("teal"),
                               "Main only" = colv("blue"),
                               "Strict only" = colv("amber"))) +
  labs(title = "Sensitivity overlap by disease pair",
       x = "Pair-gene records", y = NULL, fill = NULL) +
  theme_pub(5.5) +
  theme(axis.text.y = element_text(size = 4.2),
        legend.position = "bottom")
write_source(ovs, "SupplementaryFigure12_overlap_summary")
write_source(pair_count, "SupplementaryFigure12_pair_overlap")
save_supp(list(A = s7a, B = s7b), "SupplementaryFigure12", h = 145)

# Supplementary Figure 13: pathway detail
pt <- read_tsv(file.path(END,
  "final_gtexv8_twas_total_domain_fuma_enrichment", "enrichment_results",
  "domain", "enrichment_global_and_domain_summary_top_terms.tsv"))
pt <- pt[is.finite(pt$FDR) & pt$FDR > 0, ]
pt <- pt[order(pt$FDR), ]
pt$group <- ifelse(grepl("global", pt$analysis_name), "Global",
                   ifelse(grepl("asthma", pt$analysis_name),
                          "Asthma/lower airway",
                          ifelse(grepl("atopic", pt$analysis_name),
                                 "Atopic/allergic", "Other domains")))
pt <- do.call(rbind, lapply(split(pt, interaction(pt$group, pt$database)),
                            head, 4))
pt$term_clean <- gsub("^Homo sapiens: ", "", pt$term_name)
pt$term_clean <- gsub(" pathway$", "", pt$term_clean, ignore.case = TRUE)
pt$term <- wrap_text(pt$term_clean, 38)
pt$neglog10_fdr <- -log10(pt$FDR)
pt$x_key <- paste(pt$group, pt$database, sep = "\n")
x_levels <- unique(pt$x_key[order(match(pt$group,
                                        c("Global", "Asthma/lower airway",
                                          "Atopic/allergic", "Other domains")),
                                 match(pt$database,
                                       c("GO_BP", "KEGG", "Reactome")))])
term_score <- tapply(pt$neglog10_fdr, pt$term, max, na.rm = TRUE)
term_levels <- names(sort(term_score, decreasing = FALSE))
pt$x_key <- factor(pt$x_key, levels = x_levels)
pt$term <- factor(pt$term, levels = term_levels)
pt_top <- aggregate(term_id ~ x_key, pt, function(x) length(unique(x)))
names(pt_top)[2] <- "n_terms"
pt_side <- aggregate(x_key ~ term, pt, function(x) length(unique(x)))
names(pt_side)[2] <- "n_panels"
s8_top <- ggplot(pt_top, aes(x_key, n_terms)) +
  geom_col(fill = colv("np"), width = .68, alpha = .82) +
  geom_text(aes(label = n_terms), vjust = -0.15, size = 1.8,
            colour = colv("graphite")) +
  scale_y_continuous(expand = expansion(mult = c(0, .22))) +
  labs(x = NULL, y = "Terms") +
  theme_pub(4.5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 4.4),
        plot.margin = margin(2, 8, 0, 8))
s8_main <- ggplot(pt, aes(x_key, term)) +
  geom_point(aes(size = overlap_count, fill = database,
                 alpha = pmin(neglog10_fdr, 5)),
             shape = 21, colour = "white", stroke = .18) +
  scale_fill_manual(values = c("GO_BP" = colv("ar"),
                               "KEGG" = colv("amber"),
                               "Reactome" = colv("violet"))) +
  scale_alpha_continuous(range = c(.55, 1), guide = "none") +
  scale_size_continuous(range = c(1.2, 4.2)) +
  labs(title = "Leading pathway patterns by domain and database",
       subtitle = "Central matrix displays top terms per analysis group/database; margins count repeated support",
       x = NULL, y = NULL, fill = "Database", size = "Overlap") +
  theme_pub(4.8) +
  theme(axis.text.x = element_text(angle = 42, hjust = 1, size = 4.0,
                                   lineheight = .86),
        axis.text.y = element_text(size = 3.35, lineheight = .76),
        legend.position = "bottom",
        panel.grid.major.x = element_line(colour = "#EEF2F4", linewidth = .16),
        panel.grid.major.y = element_line(colour = "#F5F7F8", linewidth = .12))
s8_side <- ggplot(pt_side, aes(n_panels, term)) +
  geom_col(fill = colv("comparator"), width = .68, alpha = .72) +
  geom_text(aes(label = n_panels), hjust = -0.2, size = 1.8,
            colour = colv("graphite")) +
  scale_x_continuous(expand = expansion(mult = c(0, .25))) +
  labs(x = "Panels", y = NULL) +
  theme_pub(4.5) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 4.3),
        plot.margin = margin(18, 6, 6, 0))
s8a <- (s8_top + plot_spacer()) / (s8_main + s8_side) +
  plot_layout(widths = c(8.6, 1.2), heights = c(1.0, 8.0))
write_source(pt, "SupplementaryFigure13_pathways")
save_supp(list(A = s8a), "SupplementaryFigure13", h = 260)

# Supplementary Figure 14: FUMA GENE2FUNC tissues
f54 <- read_tsv(tab("Table_R9_FUMA_GTEx54_tissue_FDR05.tsv"))
f30 <- read_tsv(tab("Table_R10_FUMA_GTEx30_general_tissue_FDR05.tsv"))
fuma_plot <- function(x, title) {
  x$analysis <- gsub("_", " ", x$analysis_id)
  x$GeneSet <- gsub("_", " ", x$GeneSet)
  category_colours <- c("DEG.down" = colv("blue"),
                        "DEG.twoside" = colv("violet"),
                        "DEG.up" = colv("amber"))
  ggplot(x, aes(-log10(adjP), reorder(GeneSet, -log10(adjP)),
                size = N_overlap, colour = Category)) +
    geom_point(alpha = .8) + facet_wrap(~ analysis, scales = "free_y") +
    scale_colour_manual(values = category_colours, drop = FALSE) +
    scale_size_continuous(range = c(1.4, 4.0)) +
    labs(title = title,
         subtitle = "Colour denotes DEG category; point size denotes overlap",
         x = expression(-log[10](FDR)), y = NULL,
         colour = "DEG category", size = "Overlap") +
    guides(colour = guide_legend(order = 1, nrow = 1),
           size = guide_legend(order = 2, nrow = 1)) +
    theme_pub(5.2) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 4.2),
          legend.title = element_text(size = 4.4))
}
s9a <- fuma_plot(f54, "GTEx v8 54-tissue enrichment")
s9b <- fuma_plot(f30, "GTEx v8 30-general-tissue enrichment") +
  theme(legend.position = "none")
write_source(f54, "SupplementaryFigure14_GTEx54_tissue_enrichment")
write_source(f30, "SupplementaryFigure14_GTEx30_tissue_enrichment")
save_supp(list(A = s9a, B = s9b), "SupplementaryFigure14", h = 125)

# Supplementary Figure 15: cell-type sensitivity
ca <- read_tsv(file.path(END, "integrated_fuma_tissue_cell_ora",
  "cell_marker_ora", "gene_grade_background_ora",
  "gene_grade_background_cell_marker_ORA_audit.tsv"))
ca <- ca[ca$background_name == "AB1456" |
           grepl("_AB$", ca$background_name), ]
ca <- ca[ca$comparison %in%
           c("GeneA_or_B", "remaining_AB_background_not_TWAS_primary"), ]
ca$marker_source <- gsub("_", " ", ca$marker_source)
pair_cell <- read_tsv(tab("Table_R11_pair_level_cell_marker_ORA_FDR05.tsv"))
web <- read_tsv(tab("Table_R12_WebCSEA_Bonferroni_result.tsv"))
focus <- rbind(
  data.frame(source = "Marker ORA", label = pair_cell$term_name,
             score = -log10(pair_cell$FDR)),
  data.frame(source = "WebCSEA", label = web$Tissue_cell_type_name,
             score = -log10(web$input_list_combined_p))
)
focus$label <- trimws(gsub("_+", " ", focus$label))
s10b <- ggplot(focus, aes(score, reorder(label, score), colour = source)) +
  geom_segment(aes(x = 0, xend = score, yend = reorder(label, score)),
               colour = colv("pale"), linewidth = 2) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("Marker ORA" = colv("amber"),
                                 "WebCSEA" = colv("violet"))) +
  labs(title = "Exploratory cell-type findings",
       x = expression(-log[10](adjusted~P)), y = NULL, colour = NULL) +
  theme_pub(6.2) + theme(legend.position = "bottom")
write_source(ca, "SupplementaryFigure15_cell_audit")
write_source(focus, "SupplementaryFigure15_exploratory_findings")
save_supp_single(s10b, "SupplementaryFigure15", h = 95)

# Supplementary Figure 16: FUMA annotation and concordance
fa <- read_tsv(file.path(END, "fuma_snp_annotation_concordance",
                         "parsed_fuma_annotation",
                         "fuma_AB_snp_functional_annotation.tsv"))
fa_sum <- aggregate(cbind(CADD_gt12 = as.integer(fa$CADD_gt12 %in% c(TRUE, "True")),
                          RegulomeDB = as.integer(fa$strong_RegulomeDB_evidence %in%
                                                   c(TRUE, "True"))) ~
                      locus_grade, fa, sum)
fa_long <- rbind(
  data.frame(grade = fa_sum$locus_grade, annotation = "CADD > 12",
             n = fa_sum$CADD_gt12),
  data.frame(grade = fa_sum$locus_grade, annotation = "Strong RegulomeDB",
             n = fa_sum$RegulomeDB)
)
s11a <- ggplot(fa_long, aes(annotation, n, fill = grade)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(A = colv("teal"), B = colv("amber"))) +
  labs(title = "FUMA SNP functional annotation burden",
       x = NULL, y = "Annotated SNPs", fill = "Locus") +
  theme_pub(6.5) + theme(legend.position = "bottom")
chain <- read_tsv(file.path(END, "fuma_snp_annotation_concordance",
                            "parsed_fuma_annotation",
                            "variant_to_TWAS_gene_concordance.tsv"))
chain <- chain[order(chain$best_TWAS_fdr), ]
chain <- chain[!duplicated(chain$mapped_gene), ]
chain <- head(chain, 20)
cl <- do.call(rbind, lapply(seq_len(nrow(chain)), function(i) {
  z <- chain[i, ]
  data.frame(gene = z$mapped_gene,
    evidence = c("Gene-A", "Gene-B", "Strict TWAS", "CADD>12", "RegulomeDB"),
    present = c(z$in_GeneA, z$in_GeneB, z$in_strict_positive,
                z$CADD_gt12, z$strong_RegulomeDB_evidence) %in% c(TRUE, "True"))
}))
s11b <- ggplot(cl, aes(evidence, factor(gene, levels = rev(unique(gene))),
                       fill = present)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c("TRUE" = colv("teal"),
                               "FALSE" = colv("pale")), guide = "none") +
  labs(title = "Variant-to-TWAS-gene concordance", x = NULL, y = NULL) +
  theme_pub(6) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
write_source(fa_long, "SupplementaryFigure16_FUMA_burden")
write_source(cl, "SupplementaryFigure16_variant_TWAS_concordance")
save_supp(list(A = s11a, B = s11b), "SupplementaryFigure16", h = 115)

# Supplementary Figure 17: LCV and bidirectional MR
lcv <- read_tsv(tab("Table_R14_LCVMR_final31_summary.tsv"))
lcv$partner <- pair_partner(lcv$pair_id)
lcv$pair_axis <- pair_axis_label(lcv$pair_id)
lcv$lo <- lcv$LCV_gcp - 1.96 * lcv$LCV_gcp_se
lcv$hi <- lcv$LCV_gcp + 1.96 * lcv$LCV_gcp_se
lcv$status <- ifelse(lcv$LCV_gcp_p < .05 / 31, "Bonferroni",
                     ifelse(lcv$LCV_gcp_p < .05, "Nominal", "Not significant"))
s12a <- ggplot(lcv, aes(LCV_gcp, reorder(pair_axis, LCV_gcp),
                        colour = status)) +
  geom_vline(xintercept = 0, linetype = 2, colour = colv("grey")) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .35) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c("Bonferroni" = colv("teal"),
                                 "Nominal" = colv("amber"),
                                 "Not significant" = colv("grey"))) +
  labs(title = "LCV directionality across 31 pairs",
       x = "Genetic causality proportion", y = NULL, colour = NULL) +
  theme_pub(5.8) + theme(legend.position = "bottom")
mr <- read_tsv(file.path(ROOT, "LCVMR_data", "results_final31",
                         "MR", "mr_results_all.tsv"))
ivw <- mr[mr$method == "IVW", ]
qc <- read_tsv(tab("Table_R15_MR_direction_QC_final31.tsv"))
all_dir <- merge(qc, ivw[, c("direction", "beta", "se", "pval")],
                 by = "direction", all.x = TRUE)
all_dir$label <- paste(short_trait(all_dir$exposure), "->",
                       short_trait(all_dir$outcome))
all_dir$status2 <- ifelse(all_dir$status != "eligible", "Insufficient",
                          ifelse(all_dir$pval < .05 / 62, "Bonferroni",
                                 ifelse(all_dir$pval < .05, "Nominal",
                                        "Not significant")))
all_dir$lo <- all_dir$beta - 1.96 * all_dir$se
all_dir$hi <- all_dir$beta + 1.96 * all_dir$se
all_dir$exposure_anchor <- ifelse(grepl("^ALLERGIC_RHINITIS", all_dir$exposure),
                                  "AR as exposure",
                                  ifelse(grepl("^NASAL_POLYPS", all_dir$exposure),
                                         "NP as exposure", "Comparator as exposure"))
finite_ci <- c(all_dir$lo, all_dir$hi)
finite_ci <- finite_ci[is.finite(finite_ci)]
clip_rng <- as.numeric(stats::quantile(finite_ci, c(.03, .97), na.rm = TRUE))
all_dir$lo_plot <- pmax(all_dir$lo, clip_rng[1], na.rm = TRUE)
all_dir$hi_plot <- pmin(all_dir$hi, clip_rng[2], na.rm = TRUE)
all_dir$beta_plot <- pmin(pmax(all_dir$beta, clip_rng[1]), clip_rng[2])
all_dir$label_f <- factor(all_dir$label,
                          levels = all_dir$label[order(all_dir$beta_plot,
                                                       all_dir$pval,
                                                       na.last = TRUE)])
s12b <- ggplot(all_dir, aes(beta_plot, label_f, colour = status2)) +
  geom_vline(xintercept = 0, linetype = 2, colour = colv("grey")) +
  geom_errorbarh(aes(xmin = lo_plot, xmax = hi_plot),
                 height = 0, linewidth = .25, na.rm = TRUE) +
  geom_point(aes(shape = exposure_anchor), size = 1.35, na.rm = TRUE) +
  scale_colour_manual(values = pal_status, drop = FALSE) +
  scale_shape_manual(values = c("AR as exposure" = 16,
                                "NP as exposure" = 17,
                                "Comparator as exposure" = 15)) +
  labs(title = "IVW forest across 62 prespecified directions",
       subtitle = "CIs are clipped at the display limits; complete estimates are retained in source data",
       x = "IVW beta (95% CI, display-clipped)", y = NULL,
       colour = NULL, shape = NULL) +
  theme_pub(4.5) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        axis.text.y = element_text(size = 3.3))
sens <- read_tsv(tab("Table_R16_MR_sensitivity_final31.tsv"))
methods <- mr[mr$method %in% c("IVW", "MR-Egger", "Weighted Median", "MR-PRESSO"), ]
mm <- reshape(methods[, c("direction", "method", "pval")],
              idvar = "direction", timevar = "method", direction = "wide")
rob <- merge(sens, mm, by = "direction", all.x = TRUE)
rob_long <- do.call(rbind, lapply(seq_len(nrow(rob)), function(i) {
  z <- rob[i, ]
  data.frame(direction = paste(short_trait(z$exposure), "->",
                               short_trait(z$outcome)),
    check = c("IVW", "MR-Egger", "Weighted median", "MR-PRESSO",
              "No heterogeneity", "Egger intercept OK", "PRESSO global OK"),
    pass = c(is.finite(z[["pval.IVW"]]) && z[["pval.IVW"]] < .05,
             is.finite(z[["pval.MR-Egger"]]) && z[["pval.MR-Egger"]] < .05,
             is.finite(z[["pval.Weighted Median"]]) &&
               z[["pval.Weighted Median"]] < .05,
             is.finite(z[["pval.MR-PRESSO"]]) &&
               z[["pval.MR-PRESSO"]] < .05,
             is.finite(z$heterogeneity_p) && z$heterogeneity_p >= .05,
             is.finite(z$egger_p) && z$egger_p >= .05,
             is.finite(z$mr_presso_global_p) && z$mr_presso_global_p >= .05)
  )
}))
s12c <- ggplot(rob_long, aes(check, factor(direction,
                                           levels = rev(unique(direction))),
                               fill = pass)) +
  geom_tile(colour = "white", linewidth = .15) +
  scale_fill_manual(values = c("TRUE" = colv("teal"),
                               "FALSE" = colv("absent")), guide = "none") +
  labs(title = "Method concordance and sensitivity flags",
       x = NULL, y = NULL) +
  theme_pub(4.3) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
write_source(all_dir, "SupplementaryFigure17_MR_directions")
write_source(rob_long, "SupplementaryFigure17_MR_robustness")
save_supp(list(A = s12a, B = s12b, C = s12c),
          "SupplementaryFigure17", design = "AB\nCC", h = 185)

# Supplementary Figure 5: genomic positions of prioritized loci
chr_df <- read_tsv(file.path(ROOT, "figure_v2_source_data",
                             "chromosome_offsets.tsv"))
chr_df$chr <- as.integer(chr_df$chr)
chr_df$chr_start <- chr_df$offset
chr_df$chr_end <- chr_df$offset + chr_df$chr_len
chr_df$chr_mid <- (chr_df$chr_start + chr_df$chr_end) / 2

loci14 <- read_tsv(file.path(ROOT, "figure_v2_source_data", "locus_ab.tsv"))
loci14$anchor <- pair_anchor(loci14$pair_id)
loci14$pair_axis <- pair_axis_label(loci14$pair_id)
loci14$locus_grade <- factor(loci14$locus_grade, levels = c("A", "B"))
loci14$anchor <- factor(loci14$anchor,
                         levels = c("Nasal polyps", "Allergic rhinitis"))
s14a <- ggplot() +
  geom_segment(data = chr_df,
               aes(x = chr_start, xend = chr_end, y = 0, yend = 0,
                   colour = factor(chr %% 2)),
               linewidth = 3.2, lineend = "butt") +
  geom_point(data = loci14,
             aes(genome_pos, as.numeric(anchor), fill = locus_grade,
                 size = n_overlapping_core_snps),
             shape = 21, colour = "white", stroke = .16, alpha = .90) +
  geom_text(data = chr_df, aes(chr_mid, 2.55, label = chr),
            size = 1.8, colour = colv("grey")) +
  scale_colour_manual(values = c("0" = colv("chr_light"),
                                 "1" = colv("chr_dark")),
                      guide = "none") +
  scale_fill_manual(values = pal_grade, drop = FALSE) +
  scale_size_continuous(range = c(1.0, 4.2)) +
  scale_y_continuous(breaks = c(1, 2),
                     labels = c("Nasal polyps", "Allergic rhinitis"),
                     limits = c(-.25, 2.8)) +
  labs(title = "Genomic positions of prioritized Locus-A/B regions",
       subtitle = "Points show prioritized loci; chromosome labels mark genome-wide position",
       x = NULL, y = NULL, fill = "Locus grade", size = "Core SNPs") +
  theme_pub(6.0) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank())

locus_chr <- aggregate(locus_id ~ pair_axis + chr + locus_grade, loci14, length)
names(locus_chr)[4] <- "n_loci"
locus_chr$pair_axis <- factor(
  locus_chr$pair_axis,
  levels = rev(names(sort(table(loci14$pair_axis))))
)
s14b <- ggplot(locus_chr, aes(factor(chr), pair_axis, fill = n_loci)) +
  geom_tile(colour = "white", linewidth = .16) +
  facet_wrap(~ locus_grade, nrow = 1) +
  scale_fill_gradient(low = colv("absent"), high = colv("comparator")) +
  labs(title = "Pair-by-chromosome locus burden",
       x = "Chromosome", y = NULL, fill = "Loci") +
  theme_pub(5.2) +
  theme(axis.text.x = element_text(size = 4.2, angle = 90,
                                   hjust = 1, vjust = .5),
        axis.text.y = element_text(size = 4.2),
        legend.position = "bottom",
        panel.grid = element_blank())

locus_chr_total <- aggregate(locus_id ~ chr + locus_grade, loci14, length)
names(locus_chr_total)[3] <- "n_loci"
s14c <- ggplot(locus_chr_total, aes(factor(chr), n_loci, fill = locus_grade)) +
  geom_col(width = .72, position = "stack") +
  scale_fill_manual(values = pal_grade, drop = FALSE) +
  labs(title = "Locus-A/B counts by chromosome",
       x = "Chromosome", y = "Prioritized loci", fill = "Locus") +
  theme_pub(6.2) +
  theme(legend.position = "bottom")

recurrent14 <- read_tsv(file.path(FIG, "source_data", "recurrent_loci.tsv"))
recurrent14$anchor <- pair_anchor(recurrent14$pair_id)
recurrent14$pair_axis <- pair_axis_label(recurrent14$pair_id)
recurrent14$region_key <- paste(recurrent14$chr,
                                round(as.numeric(recurrent14$mid) / 1e6),
                                sep = ":")
region_rank14 <- aggregate(pair_id ~ region_key, recurrent14,
                           function(x) length(unique(x)))
names(region_rank14)[2] <- "n_pairs"
region_support14 <- aggregate(n_overlapping_core_snps ~ region_key,
                              recurrent14, sum, na.rm = TRUE)
names(region_support14)[2] <- "support_n"
region_rank14 <- merge(region_rank14, region_support14,
                       by = "region_key", all.x = TRUE)
region_rank14 <- region_rank14[order(-region_rank14$n_pairs,
                                     -region_rank14$support_n,
                                     region_rank14$region_key), ]
top_region14 <- head(region_rank14$region_key, 12)
locus_edges <- recurrent14[recurrent14$region_key %in% top_region14, ]
locus_edges <- locus_edges[order(locus_edges$region_key,
                                 locus_edges$locus_grade != "A",
                                 -locus_edges$support_n), ]
locus_edges <- locus_edges[!duplicated(locus_edges[, c("region_key",
                                                        "pair_id")]), ]
region_nodes14 <- do.call(rbind, lapply(top_region14, function(k) {
  z <- locus_edges[locus_edges$region_key == k, ]
  grade <- if (any(z$locus_grade == "A")) "A" else "B"
  gene <- z$nearest_gene[!is.na(z$nearest_gene) & nzchar(z$nearest_gene)][1]
  band <- z$cytoband_name[!is.na(z$cytoband_name) &
                            nzchar(z$cytoband_name)][1]
  if (!length(gene) || is.na(gene)) gene <- ""
  if (!length(band) || is.na(band)) band <- paste0("chr", z$chr[1])
  data.frame(region_key = k,
             node_label = ifelse(nzchar(gene), paste(band, gene), band),
             node_class = grade,
             support_n = sum(z$n_overlapping_core_snps, na.rm = TRUE))
}))
pair_nodes14 <- unique(locus_edges[, c("pair_axis", "anchor", "P")])
pair_nodes14$pair_number <- as.integer(sub("^P", "", pair_nodes14$P))
pair_nodes14 <- pair_nodes14[order(pair_nodes14$pair_number), ]
region_nodes14$angle <- seq(115, 245, length.out = nrow(region_nodes14))
pair_nodes14$angle <- seq(-82, 82, length.out = nrow(pair_nodes14))
node_xy14 <- rbind(
  data.frame(node_id = region_nodes14$region_key,
             node_label = region_nodes14$node_label,
             node_class = region_nodes14$node_class,
             angle = region_nodes14$angle),
  data.frame(node_id = pair_nodes14$pair_axis,
             node_label = pair_nodes14$P,
             node_class = as.character(pair_nodes14$anchor),
             angle = pair_nodes14$angle)
)
node_xy14$x <- cos(node_xy14$angle * pi / 180)
node_xy14$y <- sin(node_xy14$angle * pi / 180)
node_xy14$label_x <- 1.17 * node_xy14$x
node_xy14$label_y <- 1.17 * node_xy14$y
node_xy14$hjust <- ifelse(node_xy14$x < 0, 1, 0)
node_xy14$label_angle <- 0
pair_label_rows14 <- node_xy14$node_class %in%
  c("Allergic rhinitis", "Nasal polyps")
node_xy14$label_x[pair_label_rows14] <- 1.24 * node_xy14$x[pair_label_rows14]
node_xy14$label_y[pair_label_rows14] <- 1.24 * node_xy14$y[pair_label_rows14]
node_xy14$label_angle[pair_label_rows14] <-
  node_xy14$angle[pair_label_rows14]
locus_edges <- merge(locus_edges,
                     node_xy14[, c("node_id", "x", "y")],
                     by.x = "region_key", by.y = "node_id", all.x = TRUE)
names(locus_edges)[names(locus_edges) %in% c("x", "y")] <- c("x0", "y0")
locus_edges <- merge(locus_edges,
                     node_xy14[, c("node_id", "x", "y")],
                     by.x = "pair_axis", by.y = "node_id", all.x = TRUE)
names(locus_edges)[names(locus_edges) %in% c("x", "y")] <- c("x1", "y1")
chord14 <- do.call(rbind, lapply(seq_len(nrow(locus_edges)), function(i) {
  z <- locus_edges[i, ]
  tt <- seq(0, 1, length.out = 45)
  data.frame(edge_id = i, t = tt,
             x = (1 - tt)^2 * z$x0 + tt^2 * z$x1,
             y = (1 - tt)^2 * z$y0 + tt^2 * z$y1,
             anchor = z$anchor,
             support = z$n_overlapping_core_snps,
             region_key = z$region_key, pair_id = z$pair_id)
}))
s14d <- ggplot() +
  annotate("path", x = cos(seq(0, 2 * pi, length.out = 240)),
           y = sin(seq(0, 2 * pi, length.out = 240)),
           colour = colv("pale"), linewidth = 2.2) +
  geom_path(data = chord14,
            aes(x, y, group = edge_id, colour = anchor,
                linewidth = pmax(1, support)),
            alpha = .20, lineend = "round") +
  geom_point(data = node_xy14,
             aes(x, y, fill = node_class),
             shape = 21, size = 2.7, colour = "white", stroke = .22) +
  geom_text(data = node_xy14,
            aes(label_x, label_y, label = node_label, hjust = hjust,
                angle = label_angle),
            size = 1.85, colour = colv("graphite")) +
  scale_colour_manual(values = pal_anchor, guide = "none") +
  scale_fill_manual(values = c(pal_grade, pal_anchor), guide = "none") +
  scale_linewidth_continuous(range = c(.18, .72), guide = "none") +
  coord_equal(xlim = c(-1.48, 1.48), ylim = c(-1.30, 1.30),
              clip = "off") +
  labs(title = "Recurrent loci shared across disease pairs",
       subtitle = paste0(
         "Top 12 recurrent regions linked to ",
         length(unique(locus_edges$pair_id)),
         " disease pairs (P1-P31); regions are grouped by chromosome and rounded megabase position"
       )) +
  theme_void(base_size = 6) +
  theme(plot.title = element_text(face = "bold", colour = colv("graphite"),
                                  size = 7.3),
        plot.subtitle = element_text(colour = colv("grey"), size = 5.2),
        plot.margin = margin(8, 22, 8, 22))
write_source(loci14, "SupplementaryFigure5_prioritized_loci")
write_source(locus_edges, "SupplementaryFigure5_circular_locus_pair_map")
save_supp(list(A = s14a, B = s14b, C = s14c, D = s14d),
          "SupplementaryFigure5", design = "AA\nBC\nDD", h = 225)

# Supplementary Figure 6: cross-trait SNP association atlas
mtag15 <- read_tsv(file.path(ROOT, "figure_v2_source_data",
                             "mtag_top_loci.tsv"))
mtag15$score <- pmin(-log10(as.numeric(mtag15$p_mtag)), 60)
mtag15$display_class <- factor(mtag15$display_class,
                               levels = c("suggestive", "candidate"))
mtag15$region_mb <- floor(as.numeric(mtag15$BP) / 1e6)
mtag15$region_bin <- paste0("chr", mtag15$chr, ":",
                            mtag15$region_mb, "-",
                            mtag15$region_mb + 1, "Mb")
s15a <- ggplot(mtag15, aes(genome_pos, score, colour = display_class)) +
  geom_point(size = .34, alpha = .62) +
  geom_hline(yintercept = -log10(5e-8), linetype = 2,
             colour = colv("grey"), linewidth = .25) +
  geom_text(data = chr_df, aes(chr_mid, 61, label = chr),
            inherit.aes = FALSE, size = 1.7, colour = colv("grey")) +
  scale_colour_manual(values = c("candidate" = colv("np"),
                                 "suggestive" = colv("line_dark")),
                      drop = FALSE) +
  scale_y_continuous(limits = c(0, 62), expand = expansion(mult = c(0, .02))) +
  labs(title = "Cross-trait MTAG association atlas",
       subtitle = "Top SNP associations retained for downstream cross-trait locus evaluation",
       x = NULL, y = expression(-log[10](MTAG~P)), colour = NULL) +
  theme_pub(5.7) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid.major.x = element_blank())

candidate15 <- mtag15[mtag15$display_class == "candidate", ]
pair_cand <- aggregate(SNP ~ pair_id, candidate15,
                       function(x) length(unique(x)))
names(pair_cand)[2] <- "n_candidate_snps"
pair_cand$anchor <- pair_anchor(pair_cand$pair_id)
pair_cand$partner <- pair_partner(pair_cand$pair_id)
pair_cand$pair_axis <- pair_axis_label(pair_cand$pair_id)
pair_cand <- pair_cand[order(-pair_cand$n_candidate_snps,
                             pair_cand$anchor, pair_cand$partner), ]
pair_plot <- head(pair_cand, 24)
pair_plot$pair_axis <- factor(pair_plot$pair_axis,
                              levels = rev(pair_plot$pair_axis))
s15b <- ggplot(pair_plot, aes(n_candidate_snps, pair_axis, fill = anchor)) +
  geom_col(width = .68) +
  geom_text(aes(label = n_candidate_snps), hjust = -0.15,
            size = 1.9, colour = colv("graphite")) +
  scale_fill_manual(values = pal_anchor) +
  scale_x_continuous(expand = expansion(mult = c(0, .16))) +
  labs(title = "Candidate SNP burden by pair",
       x = "Candidate SNPs", y = NULL, fill = NULL) +
  theme_pub(5.4) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 3.8))

region_cand <- aggregate(pair_id ~ region_bin + chr, candidate15,
                         function(x) length(unique(x)))
names(region_cand)[3] <- "n_pairs"
region_counts <- as.data.frame(table(candidate15$region_bin))
names(region_counts) <- c("region_bin", "total_snps")
region_cand <- merge(region_cand, region_counts, by = "region_bin", all.x = TRUE)
region_cand <- region_cand[order(-region_cand$n_pairs,
                                 -region_cand$total_snps,
                                 region_cand$region_bin), ]
region_plot <- head(region_cand, 28)
region_plot$region_bin <- factor(region_plot$region_bin,
                                 levels = rev(region_plot$region_bin))
s15c <- ggplot(region_plot, aes(n_pairs, region_bin)) +
  geom_segment(aes(x = 0, xend = n_pairs, yend = region_bin),
               colour = colv("pale"), linewidth = 1.5) +
  geom_point(aes(size = total_snps), colour = colv("teal"), alpha = .88) +
  scale_size_continuous(range = c(1.2, 3.8)) +
  labs(title = "Recurrent candidate regions",
       x = "Disease pairs", y = NULL, size = "SNPs") +
  theme_pub(5.2) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 3.8))
write_source(mtag15, "SupplementaryFigure6_MTAG_top_loci")

# Supplementary Figure 6D: radial recurrent-region summary
reg_pair <- aggregate(SNP ~ region_bin + pair_id, candidate15,
                      function(x) length(unique(x)))
names(reg_pair)[3] <- "n_snps"
reg_pair$anchor <- pair_anchor(reg_pair$pair_id)
reg_pair$partner <- pair_partner(reg_pair$pair_id)
reg_pair$pair_axis <- pair_axis_label(reg_pair$pair_id)
region_anchor16 <- aggregate(pair_id ~ region_bin + anchor, reg_pair,
                             function(x) length(unique(x)))
names(region_anchor16)[3] <- "n_pairs"
region_snp16 <- aggregate(n_snps ~ region_bin, reg_pair, sum)
region_total16 <- aggregate(pair_id ~ region_bin, reg_pair,
                            function(x) length(unique(x)))
names(region_total16)[2] <- "total_pairs"
region_rank16 <- merge(region_total16, region_snp16,
                       by = "region_bin", all = TRUE)
region_rank16 <- region_rank16[order(-region_rank16$total_pairs,
                                     -region_rank16$n_snps,
                                     region_rank16$region_bin), ]
top_regions16 <- head(region_rank16$region_bin, 18)
radial16 <- region_anchor16[region_anchor16$region_bin %in% top_regions16, ]
radial16 <- merge(radial16, region_rank16, by = "region_bin", all.x = TRUE)
radial16$region_bin <- factor(radial16$region_bin,
                              levels = rev(top_regions16))
radial16$x_id <- as.numeric(radial16$region_bin) + 2
region_points16 <- unique(radial16[, c("region_bin", "x_id",
                                       "total_pairs", "n_snps")])
region_points16$label_y <- region_points16$total_pairs + 1.65
region_points16$angle <- 90 -
  360 * (region_points16$x_id - .5) / 22
region_points16$hjust <- ifelse(region_points16$angle < -90, 1, 0)
region_points16$angle <- ifelse(region_points16$angle < -90,
                                region_points16$angle + 180,
                                region_points16$angle)
s16a <- ggplot(radial16,
               aes(x_id, n_pairs, fill = anchor)) +
  geom_col(width = .78, colour = "white", linewidth = .16) +
  geom_point(data = region_points16,
             aes(x_id, total_pairs + .55, size = n_snps),
             inherit.aes = FALSE, shape = 21, fill = colv("teal"),
             colour = "white", stroke = .18) +
  geom_text(data = region_points16,
            aes(x_id, label_y, label = region_bin,
                angle = angle, hjust = hjust),
            inherit.aes = FALSE, size = 1.90,
            colour = colv("graphite")) +
  scale_fill_manual(values = pal_anchor, name = NULL) +
  scale_size_continuous(range = c(1.4, 4.2), name = "Candidate SNPs") +
  scale_x_continuous(limits = c(.5, 22.5)) +
  scale_y_continuous(limits = c(-3.2, max(region_points16$label_y) + 2.5),
                     breaks = NULL) +
  coord_polar(start = 0, clip = "off") +
  labs(title = "Candidate regions recurring across disease pairs",
       subtitle = "Radial height denotes pair recurrence; stacked colours separate AR and NP pairs",
       x = NULL, y = NULL) +
  theme_void(base_size = 6) +
  theme(plot.title = element_text(face = "bold", colour = colv("graphite"),
                                  size = 7.3),
        plot.subtitle = element_text(colour = colv("grey"), size = 5.2),
        plot.margin = margin(8, 24, 8, 24),
        legend.position = "bottom",
        legend.box = "horizontal")
write_source(radial16, "SupplementaryFigure6_radial_region_summary")
save_supp(list(A = s15a, B = s15b, C = s15c, D = s16a),
          "SupplementaryFigure6", design = "AA\nBC\nDD", h = 225)

# Supplementary Figure 7: candidate variant evidence summary
cand7 <- read_tsv(file.path(END, "fuma_snp_annotation_concordance",
                            "fuma_upload_inputs",
                            "AB_loci_candidate_snps_for_FUMA.tsv"))
cand7$anchor <- pair_anchor(cand7$pair_id)
cand7$partner <- pair_partner(cand7$pair_id)
cand7$pair_axis <- paste(ifelse(cand7$anchor == "Nasal polyps", "NP", "AR"),
                         cand7$partner, sep = " - ")
cand7$phenotype_domain <- gsub("_", " ", cand7$phenotype_domain)
cand7$locus_grade <- factor(cand7$locus_grade, levels = c("A", "B"))
cand7$snp_grade <- factor(cand7$snp_grade, levels = c("A", "B", "C", "D"))
cand7$method_group <- cand7$evidence_methods
method_totals <- sort(table(cand7$method_group), decreasing = TRUE)
top_method_groups <- names(head(method_totals, 8))
cand7$method_group <- ifelse(cand7$method_group %in% top_method_groups,
                             cand7$method_group, "Other combinations")
cand7$method_group <- factor(cand7$method_group,
                             levels = c(top_method_groups,
                                        "Other combinations"))

grade_domain7 <- aggregate(SNP ~ phenotype_domain + locus_grade + snp_grade,
                           cand7, function(x) length(unique(x)))
names(grade_domain7)[4] <- "n_snps"
s7a <- ggplot(grade_domain7, aes(phenotype_domain, n_snps, fill = snp_grade)) +
  geom_col(width = .72) +
  facet_wrap(~ locus_grade, nrow = 1) +
  scale_fill_manual(values = c(A = colv("teal"), B = colv("amber"),
                               C = colv("ar"), D = colv("grey")),
                    drop = FALSE) +
  labs(title = "Candidate variant distribution by locus and SNP grade",
       x = NULL, y = "Unique candidate SNPs", fill = "SNP grade") +
  theme_pub(5.6) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")

method7 <- aggregate(SNP ~ method_group, cand7,
                     function(x) length(unique(x)))
names(method7)[2] <- "n_snps"
method7 <- method7[order(-method7$n_snps, method7$method_group), ]
method7$method_group <- factor(method7$method_group,
                               levels = rev(as.character(method7$method_group)))
method7$label_x <- method7$n_snps + max(method7$n_snps) * .035
s7b <- ggplot(method7, aes(n_snps, method_group)) +
  geom_segment(aes(x = 0, xend = n_snps, yend = method_group),
               colour = colv("pale"), linewidth = 1.8) +
  geom_point(colour = colv("comparator"), size = 2.2) +
  geom_text(aes(x = label_x, label = n_snps), hjust = 0, size = 1.9,
            colour = colv("graphite")) +
  scale_x_continuous(expand = expansion(mult = c(0, .14))) +
  labs(title = "Evidence-method support",
       x = "Unique candidate SNPs", y = NULL) +
  theme_pub(5.4)

pair_variant7 <- aggregate(SNP ~ pair_axis + pair_id + anchor, cand7,
                           function(x) length(unique(x)))
names(pair_variant7)[4] <- "n_snps"
pair_variant7 <- pair_variant7[order(-pair_variant7$n_snps,
                                     pair_variant7$anchor,
                                     pair_variant7$pair_axis), ]
pair_variant7_plot <- head(pair_variant7, 26)
pair_variant7_plot$pair_axis <- factor(pair_variant7_plot$pair_axis,
                                       levels = rev(pair_variant7_plot$pair_axis))
s7c <- ggplot(pair_variant7_plot, aes(n_snps, pair_axis, fill = anchor)) +
  geom_col(width = .68) +
  geom_text(aes(label = n_snps), hjust = -0.15, size = 1.85,
            colour = colv("graphite")) +
  scale_fill_manual(values = pal_anchor) +
  scale_x_continuous(expand = expansion(mult = c(0, .16))) +
  labs(title = "Disease-pair burden of candidate variants",
       x = "Unique candidate SNPs", y = NULL, fill = NULL) +
  theme_pub(5.2) +
  theme(axis.text.y = element_text(size = 3.9),
        legend.position = "bottom")
write_source(cand7, "SupplementaryFigure7_candidate_variant_evidence")
save_supp(list(A = s7a, B = s7b, C = s7c),
          "SupplementaryFigure7", design = "AB\nCC", h = 170)

# Supplementary Figure 8: representative regional association profiles
profile_n <- aggregate(SNP ~ region_bin + chr, mtag15,
                       function(x) length(unique(x)))
names(profile_n)[3] <- "total_snps"
profile_pairs <- aggregate(pair_id ~ region_bin, candidate15,
                           function(x) length(unique(x)))
names(profile_pairs)[2] <- "n_candidate_pairs"
profile_score <- aggregate(score ~ region_bin, mtag15,
                           function(x) max(x, na.rm = TRUE))
names(profile_score)[2] <- "max_score"
profile_regions <- merge(profile_n, profile_pairs, by = "region_bin",
                         all.x = TRUE)
profile_regions <- merge(profile_regions, profile_score, by = "region_bin",
                         all.x = TRUE)
profile_regions$n_candidate_pairs[is.na(profile_regions$n_candidate_pairs)] <- 0
profile_regions <- profile_regions[order(-profile_regions$n_candidate_pairs,
                                         -profile_regions$max_score,
                                         -profile_regions$total_snps,
                                         profile_regions$region_bin), ]
profile_regions <- head(profile_regions, 8)
prof17 <- mtag15[mtag15$region_bin %in% profile_regions$region_bin, ]
prof17$anchor <- pair_anchor(prof17$pair_id)
prof17$bp_mb <- prof17$BP / 1e6
prof17$facet <- factor(prof17$region_bin, levels = profile_regions$region_bin)
lead17 <- prof17[order(prof17$facet,
                       prof17$display_class != "candidate",
                       -prof17$score,
                       prof17$SNP), ]
lead17 <- lead17[!duplicated(lead17$facet), ]
lead17$label_y <- lead17$score + 2.4
facet_range17 <- aggregate(bp_mb ~ facet, prof17,
                           function(x) c(min = min(x), max = max(x)))
facet_range17 <- data.frame(
  facet = facet_range17$facet,
  facet_min = facet_range17$bp_mb[, "min"],
  facet_max = facet_range17$bp_mb[, "max"]
)
lead17 <- merge(lead17, facet_range17, by = "facet", all.x = TRUE,
                sort = FALSE)
lead17$label_hjust <- ifelse(
  lead17$bp_mb <= lead17$facet_min +
    .12 * (lead17$facet_max - lead17$facet_min), 0,
  ifelse(lead17$bp_mb >= lead17$facet_max -
           .12 * (lead17$facet_max - lead17$facet_min), 1, .5)
)
s17a <- ggplot(prof17, aes(bp_mb, score, colour = display_class)) +
  geom_point(aes(shape = anchor), size = .40, alpha = .70) +
  geom_hline(yintercept = -log10(5e-8), linetype = 2,
             colour = colv("grey"), linewidth = .22) +
  geom_point(data = lead17,
             aes(bp_mb, score),
             inherit.aes = FALSE, shape = 23, size = 2.0,
             fill = colv("violet"), colour = "white", stroke = .22) +
  geom_text(data = lead17,
            aes(bp_mb, label_y, label = SNP, hjust = label_hjust),
            inherit.aes = FALSE, size = 1.75,
            colour = colv("violet"), fontface = "bold",
            check_overlap = TRUE) +
  facet_wrap(~ facet, scales = "free_x", ncol = 2) +
  scale_colour_manual(values = c("candidate" = colv("np"),
                                 "suggestive" = colv("line_dark")),
                      drop = FALSE) +
  scale_shape_manual(values = c("Allergic rhinitis" = 16,
                                "Nasal polyps" = 17)) +
  labs(title = "Representative regional cross-trait association profiles",
       subtitle = "Purple diamond marks the strongest retained SNP per region; LD colouring was not available locally",
       x = "Position (Mb)", y = expression(-log[10](MTAG~P)),
       colour = NULL, shape = NULL) +
  theme_pub(5.4) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 5.0),
        panel.grid.major.y = element_line(colour = "#EEF2F4", linewidth = .15))
write_source(prof17, "SupplementaryFigure8_regional_association_profiles")
write_source(lead17, "SupplementaryFigure8_regional_lead_snps")
save_supp_single(s17a, "SupplementaryFigure8", h = 170)

# Supplementary Figure 10: MAGMA gene Manhattan and recurrent genes
mg_full10 <- read_tsv(file.path(END, "visualization_inputs",
                                "10_magma_gene_manhattan.tsv"))
mg_full10$chr <- as.integer(mg_full10$chr)
mg_score10 <- aggregate(minus_log10_p ~ gene_symbol + chr + start + stop,
                        mg_full10, function(x) max(x, na.rm = TRUE))
mg_bonf10 <- aggregate(gene_bonf_p ~ gene_symbol + chr + start + stop,
                       mg_full10, function(x) min(x, na.rm = TRUE))
mg_gene10 <- merge(mg_score10, mg_bonf10,
                   by = c("gene_symbol", "chr", "start", "stop"),
                   all = TRUE)
mg_gene10 <- merge(mg_gene10, chr_df[, c("chr", "offset")],
                   by = "chr", all.x = TRUE)
mg_gene10$mid <- (as.numeric(mg_gene10$start) + as.numeric(mg_gene10$stop)) / 2
mg_gene10$genome_pos <- mg_gene10$offset + mg_gene10$mid
mg_gene10$score <- pmin(-log10(as.numeric(mg_gene10$gene_bonf_p)), 40)
mg_gene10$significant_gene <- mg_gene10$gene_bonf_p < .05
mg_gene10$chr_group <- factor(mg_gene10$chr %% 2)
label10 <- mg_gene10[mg_gene10$significant_gene, ]
label10 <- label10[order(-label10$score, label10$gene_symbol), ]
label10 <- head(label10, 26)
s18a <- ggplot(mg_gene10, aes(genome_pos, score)) +
  geom_point(aes(colour = chr_group), size = .20, alpha = .38) +
  geom_point(data = mg_gene10[mg_gene10$significant_gene, ],
             colour = colv("np"), size = .48, alpha = .90) +
  geom_hline(yintercept = -log10(.05), linetype = 2,
             colour = colv("grey"), linewidth = .22) +
  geom_text(data = chr_df, aes(chr_mid, 41, label = chr),
            inherit.aes = FALSE, size = 1.7, colour = colv("grey")) +
  geom_text(data = label10,
            aes(label = gene_symbol), size = 1.65,
            colour = colv("graphite"), vjust = -0.4,
            check_overlap = TRUE) +
  scale_colour_manual(values = c("0" = colv("chr_light"),
                                 "1" = colv("chr_dark")),
                      guide = "none") +
  scale_y_continuous(limits = c(0, 42), expand = expansion(mult = c(0, .02))) +
  labs(title = "MAGMA gene-level Manhattan plot",
       subtitle = "Each point is a gene represented by its strongest Bonferroni-corrected MAGMA signal across final pairs",
       x = NULL, y = expression(-log[10](MAGMA~Bonferroni~P))) +
  theme_pub(5.8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())

mg18 <- read_tsv(file.path(ROOT, "figure_v2_source_data",
                           "magma_significant_genes.tsv"))
mg_summary18 <- read_tsv(file.path(ROOT, "figure_v2_source_data",
                                   "magma_gene_summary.tsv"))
mg18 <- merge(mg18, mg_summary18[, c("gene_symbol", "twas_supported")],
              by = "gene_symbol", all.x = TRUE)
mg18$twas_supported <- mg18$twas_supported %in% c(TRUE, "True", "TRUE", "1")
mg18$pair_axis <- pair_axis_label(mg18$pair_id)
mg18$anchor <- pair_anchor(mg18$pair_id)
gene_rec18 <- aggregate(pair_id ~ gene_symbol + twas_supported, mg18,
                        function(x) length(unique(x)))
names(gene_rec18)[3] <- "n_pairs"
gene_rec18 <- gene_rec18[order(-gene_rec18$n_pairs,
                               gene_rec18$twas_supported,
                               gene_rec18$gene_symbol), ]
gene_plot18 <- head(gene_rec18, 30)
gene_plot18$gene_symbol <- factor(gene_plot18$gene_symbol,
                                  levels = rev(gene_plot18$gene_symbol))
s18b <- ggplot(gene_plot18, aes(n_pairs, gene_symbol,
                                colour = twas_supported)) +
  geom_segment(aes(x = 0, xend = n_pairs, yend = gene_symbol),
               colour = colv("pale"), linewidth = 1.4) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c("TRUE" = colv("teal"),
                                 "FALSE" = colv("grey")),
                      labels = c("FALSE" = "MAGMA only",
                                 "TRUE" = "TWAS-supported")) +
  labs(title = "Recurrent MAGMA-positive genes",
       x = "Disease pairs", y = NULL, colour = NULL) +
  theme_pub(5.4) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 4.0))

gene_chr18 <- aggregate(gene_symbol ~ chr + twas_supported, mg18,
                        function(x) length(unique(x)))
names(gene_chr18)[3] <- "n_genes"
s18c <- ggplot(gene_chr18, aes(factor(chr), n_genes,
                               fill = twas_supported)) +
  geom_col(width = .72) +
  scale_fill_manual(values = c("TRUE" = colv("teal"),
                               "FALSE" = colv("grey")),
                    labels = c("FALSE" = "MAGMA only",
                               "TRUE" = "TWAS-supported")) +
  labs(title = "MAGMA-positive genes by chromosome",
       x = "Chromosome", y = "Unique genes", fill = NULL) +
  theme_pub(6.0) +
  theme(legend.position = "bottom")

gene_map_nodes <- gene_rec18[order(-gene_rec18$n_pairs,
                                   -as.numeric(gene_rec18$twas_supported),
                                   gene_rec18$gene_symbol), ]
top_gene_map <- head(gene_map_nodes$gene_symbol, 15)
pair_nodes18 <- unique(mg18[, c("pair_id", "pair_axis", "anchor")])
pair_code18 <- read_tsv(file.path(FIG, "source_data", "pair_summary.tsv"))
pair_code18 <- unique(pair_code18[, c("pair_id", "P")])
pair_nodes18 <- merge(pair_nodes18, pair_code18,
                      by = "pair_id", all.x = TRUE)
pair_nodes18$pair_number <- as.integer(sub("^P", "", pair_nodes18$P))
pair_nodes18$anchor <- factor(
  pair_nodes18$anchor,
  levels = c("Allergic rhinitis", "Nasal polyps")
)
pair_nodes18 <- pair_nodes18[order(pair_nodes18$anchor,
                                   pair_nodes18$pair_number), ]
n_ar18 <- sum(pair_nodes18$anchor == "Allergic rhinitis")
pair_nodes18$x_pos <- c(seq_len(n_ar18),
                        n_ar18 + 2 + seq_len(nrow(pair_nodes18) - n_ar18))
pair_nodes18$pair_label <- pair_nodes18$P
max_x18 <- max(pair_nodes18$x_pos)
gene_nodes18 <- gene_map_nodes[gene_map_nodes$gene_symbol %in%
                                 top_gene_map, ]
gene_nodes18 <- gene_nodes18[match(top_gene_map, gene_nodes18$gene_symbol), ]
gene_nodes18$y_pos <- rev(seq_len(nrow(gene_nodes18)))
gene_grid18 <- expand.grid(gene_symbol = top_gene_map,
                           pair_id = pair_nodes18$pair_id,
                           stringsAsFactors = FALSE)
gene_grid18 <- merge(gene_grid18,
                     pair_nodes18[, c("pair_id", "pair_axis", "anchor",
                                      "P", "pair_number", "x_pos")],
                     by = "pair_id", all.x = TRUE)
gene_grid18 <- merge(gene_grid18,
                     gene_nodes18[, c("gene_symbol", "y_pos",
                                      "twas_supported", "n_pairs")],
                     by = "gene_symbol", all.x = TRUE)
gene_signal18 <- aggregate(gene_bonf_p ~ gene_symbol + pair_id, mg18,
                           function(x) min(as.numeric(x), na.rm = TRUE))
gene_signal18$signal_score <- pmin(-log10(gene_signal18$gene_bonf_p), 12)
gene_grid18 <- merge(gene_grid18, gene_signal18,
                     by = c("gene_symbol", "pair_id"), all.x = TRUE)
gene_grid18$present <- is.finite(gene_grid18$signal_score)
anchor_strip18 <- pair_nodes18
anchor_strip18$y_pos <- nrow(gene_nodes18) + .75
anchor_label18 <- aggregate(x_pos ~ anchor, pair_nodes18, mean)
anchor_label18$label <- c("AR pairs", "NP pairs")
anchor_label18$y_pos <- nrow(gene_nodes18) + 1.45
s18d <- ggplot(gene_grid18, aes(x_pos, y_pos)) +
  geom_tile(data = gene_grid18,
            width = .84, height = .76,
            fill = colv("absent"), colour = "white", linewidth = .12) +
  geom_point(data = gene_grid18[gene_grid18$present, ],
             aes(size = signal_score, colour = twas_supported),
             alpha = .92) +
  geom_tile(data = anchor_strip18,
            aes(x_pos, y_pos, fill = anchor),
            width = .88, height = .34, colour = "white", linewidth = .10,
            inherit.aes = FALSE) +
  geom_text(data = anchor_label18,
            aes(x_pos, y_pos, label = label, colour = anchor),
            size = 2.35, fontface = "bold", inherit.aes = FALSE) +
  scale_colour_manual(values = c("TRUE" = colv("teal"),
                                 "FALSE" = colv("grey"),
                                 "Allergic rhinitis" = colv("ar"),
                                 "Nasal polyps" = colv("np")),
                      labels = c("FALSE" = "MAGMA only",
                                 "TRUE" = "TWAS-supported"),
                      breaks = c("FALSE", "TRUE"),
                      name = NULL) +
  scale_fill_manual(values = pal_anchor, guide = "none") +
  scale_size_continuous(range = c(.55, 2.5),
                        name = expression(-log[10](Bonferroni~P))) +
  scale_x_continuous(
    limits = c(.3, max_x18 + .7),
    breaks = pair_nodes18$x_pos,
    labels = pair_nodes18$pair_label,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(.4, nrow(gene_nodes18) + 1.8),
    breaks = gene_nodes18$y_pos,
    labels = gene_nodes18$gene_symbol,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(clip = "off") +
  labs(title = "Recurrent MAGMA genes across disease pairs",
       subtitle = "Top 15 recurrent genes; point size shows gene-level significance and colour denotes TWAS support",
       x = NULL, y = NULL) +
  theme_pub(5.8) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 4.7, angle = 90,
                                   hjust = 1, vjust = .5),
        axis.text.y = element_text(size = 5.2, face = "italic"),
        panel.grid = element_blank(),
        plot.margin = margin(8, 10, 6, 10),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.x = grid::unit(2, "mm"))
write_source(mg_gene10, "SupplementaryFigure10_MAGMA_gene_manhattan")
write_source(mg18, "SupplementaryFigure10_MAGMA_significant_genes")
write_source(gene_grid18, "SupplementaryFigure10_gene_pair_dot_matrix")
save_supp(list(A = s18a, B = s18b, C = s18c, D = s18d),
          "SupplementaryFigure10", design = "AA\nBC\nDD", h = 225)

# Supplementary Figure 18: Methods/Results citation map
cite_map <- data.frame(
  step = factor(c("Trait and pair screening",
                  "Local/shared-locus resolution",
                  "Genome-wide SNP association atlas",
                  "Candidate variant recurrence",
                  "Regional association profiles",
                  "MAGMA gene prioritization",
                  "TWAS sensitivity",
                  "Pathway interpretation",
                  "Tissue support",
                  "Cell-marker sensitivity",
                  "FUMA concordance",
                  "Drug-target annotation",
                  "Directional analyses"),
                levels = rev(c("Trait and pair screening",
                               "Local/shared-locus resolution",
                               "Genome-wide SNP association atlas",
                               "Candidate variant recurrence",
                               "Regional association profiles",
                               "MAGMA gene prioritization",
                               "TWAS sensitivity",
                               "Pathway interpretation",
                               "Tissue support",
                               "Cell-marker sensitivity",
                               "FUMA concordance",
                               "Drug-target annotation",
                               "Directional analyses"))),
  methods_results_position = c("Results para. 1-2; GWAS QC and LDSC/LAVA Methods",
                               "Results para. 3-4; coloc/SuSiE Methods",
                               "Results para. 3-4; MTAG/PLACO/CPASSOC Methods",
                               "Results para. 3-4; candidate variant Methods",
                               "Results para. 3-4; regional association Methods",
                               "Results para. 5-7; MAGMA Methods",
                               "Results para. 5; strict TWAS sensitivity",
                               "Results para. 8-9; enrichment Methods",
                               "Results para. 10-11; tissue-support Methods",
                               "Results para. 12; exploratory cell Methods",
                               "Results para. 13; post hoc FUMA Methods",
                               "Results para. 14; drug annotation Methods",
                               "Results para. 15-16; LCV/MR Methods"),
  main_figure = c("Fig. 1A-B", "Fig. 1C",
                  "Fig. 1C", "Fig. 1C", "Fig. 1C", "Fig. 2A-B",
                  "Fig. 2", "Fig. 3A-B", "Fig. 3C", "None",
                  "None", "Fig. 4A-B", "Fig. 4C"),
  supplementary_figure = c("Supplementary Figs. 1-2",
                           "Supplementary Figs. 3, 5",
                           "Supplementary Fig. 6",
                           "Supplementary Fig. 7",
                           "Supplementary Fig. 8",
                           "Supplementary Figs. 9-10",
                           "Supplementary Figs. 11-12",
                           "Supplementary Fig. 13",
                           "Supplementary Fig. 14",
                           "Supplementary Fig. 15",
                           "Supplementary Fig. 16",
                           "Supplementary Fig. 4",
                           "Supplementary Fig. 17"),
  supporting_tables = c("Tables S1, S2, S5, S6, S8",
                        "Tables S6, S7",
                        "MTAG top-loci source table",
                        "candidate SNP source table",
                        "Regional association source table",
                        "MAGMA significant-gene source table",
                        "TWAS audit source tables",
                        "Tables R5-R7",
                        "Tables R8-R10",
                        "Tables R11-R12",
                        "Table R13",
                        "drug-target source table",
                        "Tables R14-R16"),
  evidence_role = c("Select disease pairs",
                    "Prioritize shared regions",
                    "Show SNP-level association breadth",
                    "Show candidate variant grade and method support",
                    "Show representative regional tracks",
                    "Prioritize and localize genes",
                    "Check robustness",
                    "Interpret biology",
                    "Contextualize tissues",
                    "Bound exploratory cell results",
                    "Post hoc variant-to-gene concordance",
                    "Secondary annotation only",
                    "Secondary directionality only"),
  stringsAsFactors = FALSE
)
cite_long <- rbind(
  data.frame(step = cite_map$step, track = "Methods/Results position",
             label = cite_map$methods_results_position, stringsAsFactors = FALSE),
  data.frame(step = cite_map$step, track = "Main figure",
             label = cite_map$main_figure, stringsAsFactors = FALSE),
  data.frame(step = cite_map$step, track = "Supplementary figure",
             label = cite_map$supplementary_figure, stringsAsFactors = FALSE),
  data.frame(step = cite_map$step, track = "Supporting tables",
             label = cite_map$supporting_tables, stringsAsFactors = FALSE),
  data.frame(step = cite_map$step, track = "Evidence role",
             label = cite_map$evidence_role, stringsAsFactors = FALSE)
)
cite_long$track <- factor(cite_long$track,
                          levels = c("Methods/Results position",
                                     "Main figure", "Supplementary figure",
                                     "Supporting tables", "Evidence role"))
s13a <- ggplot(cite_long, aes(track, step, fill = track)) +
  geom_tile(colour = "white", linewidth = .35, alpha = .92) +
  geom_text(aes(label = wrap_axis(label, 24)), size = 2.05,
            lineheight = .86, colour = colv("graphite")) +
  scale_fill_manual(values = c("Methods/Results position" = colv("pale"),
                               "Main figure" = colv("ar"),
                               "Supplementary figure" = colv("np"),
                               "Supporting tables" = colv("comparator"),
                               "Evidence role" = colv("line")),
                    guide = "none") +
  labs(title = "Methods/Results citation map for supplementary figures",
       subtitle = "Each row marks where the supplementary evidence should be cited in the manuscript",
       x = NULL, y = NULL) +
  theme_pub(6.2) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank())
write_source(cite_map, "SupplementaryFigure18_citation_map")
save_supp(list(A = s13a), "SupplementaryFigure18", h = 175)

finish_qa()
