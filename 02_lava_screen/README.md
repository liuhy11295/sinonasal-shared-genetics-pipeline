# 02 LAVA Screen

LAVA is a parallel screen to LDSC, not a downstream filter after LDSC.

## Scripts

- `scripts/run_lava_local_rg.R`: runs a project-specific LAVA local-rg command for every pair in `step2_input_tables/pair_manifest.tsv`. Set `LAVA_LOCAL_RG_SCRIPT` to the actual local-rg runner.
- `scripts/aggregate_lava_local_rg.py`: aggregates per-pair LAVA output and applies Bonferroni `0.05 / n_lava_tests`.
- `scripts/select_lava_all_pairs.py`: applies the same Bonferroni rule to an already existing all-pair LAVA result table.

## Handoff

MTAG selection uses:

- LDSC Bonferroni-positive pairs from `step1_ldsc_screen/ldsc_rg_significant.tsv`
- LAVA Bonferroni-positive pairs from `step2_lava_screen/lava_local_rg_significant.tsv`

The rule is `LDSC_positive OR LAVA_positive`.
