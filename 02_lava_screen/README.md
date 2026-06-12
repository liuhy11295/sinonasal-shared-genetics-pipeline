# 02 LAVA Screen

LAVA is a parallel screen to LDSC, not a downstream filter after LDSC.

## Scripts

- `scripts/run_lava_local_rg.R`: runs a project-specific LAVA local-rg command for every pair in `step2_input_tables/pair_manifest.tsv`. Set `LAVA_LOCAL_RG_SCRIPT` to the actual local-rg runner.
- `scripts/aggregate_lava_local_rg.py`: aggregates per-pair LAVA output and retains positive local correlations with LAVA-adjusted `p_adj <= 0.05`.
- `scripts/select_lava_all_pairs.py`: applies the same adjusted-P rule to an existing all-pair LAVA result table.

## Handoff

MTAG selection uses:

- LDSC BH-FDR-positive pairs from `step1_ldsc_screen/ldsc_rg_significant.tsv`
- LAVA positive-local-rg, adjusted-P-positive pairs from `step2_lava_screen/lava_local_rg_significant.tsv`

The rule is `LDSC_positive OR LAVA_positive`.
