# Publication outputs

The `figures/scripts` directory archives the current final R scripts used to
assemble the main and supplementary figures. Stable entry points are
`render_figure1.R` through `render_figure4.R`,
`render_supplementary_figures.R`, and `render_all.R`. Set `PAPER_ROOT` and
`RESULTS_ROOT` before rendering; the scripts expect the manuscript project's
derived figure source data and the analysis result tree.

The `tables/scripts` directory archives the five-table builder and QA scripts.
The `tables/source_data` directory contains the nine final TSV sheets:

| Supplementary table | Source-data sheets |
|---|---|
| ST1 | `ST1_GWAS_and_power_QC.tsv` |
| ST2 | `ST2_Analysis_parameters.tsv`, `ST2_Harmonization_QC.tsv` |
| ST3 | `ST3_Final_31_pairs.tsv` |
| ST4 | `ST4_Prioritized_loci.tsv`, `ST4_Prioritized_genes.tsv` |
| ST5 | `ST5_Significant_pathways.tsv`, `ST5_Approved_drug_targets.tsv`, `ST5_Directional_evidence.tsv` |

The nine sheets are components of five supplementary tables, not 9 or 17
separate supplementary tables.
