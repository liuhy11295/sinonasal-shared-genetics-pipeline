# Downstream annotation

This directory contains post-prioritization analyses used for interpretation:

- `run_pathway_ora.R`: GO BP, KEGG, and Reactome over-representation analysis.
- `annotate_fuma_concordance.py`: concordance with externally generated FUMA SNP2GENE annotations.
- `run_cell_marker_ora.py`: exploratory cell-marker ORA using explicit gene-grade backgrounds.
- `annotate_drug_targets.py`: exploratory DGIdb/STRING target annotation.

These analyses annotate prioritized evidence. They do not upgrade locus or gene
evidence grades, establish treatment efficacy, or provide primary causal proof.
Required directories must be passed by command-line arguments or the
environment variables documented in `.env.example`.
