# Downstream annotation

This directory contains post-prioritization analyses used for interpretation:

- `run_ora_enrichment.R`: GO BP, KEGG, and Reactome over-representation analysis.
- `build_fuma_snp_annotation_concordance.py`: concordance with externally generated FUMA SNP2GENE annotations.
- `run_gene_grade_background_cell_marker_ora.py`: exploratory cell-marker ORA using explicit gene-grade backgrounds.
- `run_pairwise_network_pharmacology.py`: exploratory DGIdb/STRING target annotation.

These analyses annotate prioritized evidence. They do not upgrade locus or gene
evidence grades, establish treatment efficacy, or provide primary causal proof.
Required directories must be passed by command-line arguments or the
environment variables documented in `.env.example`.
