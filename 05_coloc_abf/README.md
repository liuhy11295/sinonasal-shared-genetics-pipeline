# 05 Coloc ABF

This step fills the explicit coloc stage between PLACO/CPASSOC SNP/locus screening and coloc-SuSiE fine-mapping.

## Inputs

- `results/phase0_extension/step2_input_tables/candidate_snps.tsv`
- `results/integrated/locus_master_table.tsv`
- `results/placo_genomewide/placo_genomewide_significant_snps.tsv`
- `results/cpassoc_genomewide/cpassoc_genomewide_significant_snps.tsv`
- `results/phase0_server/nasal4_vs_other_package/phase0_v3_coloc_locus_manifest.tsv`
- Regional GWAS files referenced by `region_a_file` and `region_b_file`

`PHASE0_COLOC_TARGET_LOCI` can be set to a custom target-locus table with at least `pair_id` and `locus_id`.

## Scripts

- `scripts/run_locus.R`: maps screened candidate SNPs/loci to regional GWAS files, harmonizes alleles, and runs `coloc::coloc.abf`.
- `scripts/aggregate_results.py`: aggregates per-locus coloc output and writes coloc-positive handoff tables.

## Outputs

- `results/phase0_extension/step3_coloc_abf/coloc_target_loci.tsv`
- `results/phase0_extension/step3_coloc_abf/per_locus/*.coloc_abf.tsv`
- `results/phase0_extension/step3_coloc_abf/coloc_abf_results.tsv`
- `results/phase0_extension/step3_coloc_abf/coloc_abf_qc.tsv`
- `results/phase0_extension/step3_coloc_abf/coloc_positive_loci.tsv`
- `results/phase0_extension/step3_coloc_abf/candidate_snps_with_coloc_abf.tsv`

For compatibility with coloc-SuSiE, `aggregate_results.py` also syncs:

- `results/phase0_extension/step2_input_tables/coloc_positive_loci.tsv`
- `results/phase0_extension/step2_input_tables/candidate_snps.tsv`
