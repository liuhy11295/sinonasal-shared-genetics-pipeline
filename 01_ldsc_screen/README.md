# 01 LDSC Screen

This step runs genetic-correlation screening with LDSC and applies
Benjamini-Hochberg false-discovery-rate correction across tested disease pairs.

## Script

- `scripts/run_ldsc.py`

## Inputs

- `results/phase0_extension/step2_input_tables/pair_manifest.tsv`
- Pair-specific LDSC munged sumstats from `trait1_munged` and `trait2_munged`
- LDSC reference:
  - `data/reference/ldsc/1000G_Phase3_ldscores/LDscore`
  - `data/reference/ldsc/1000G_Phase3_weights_hm3_no_MHC/1000G_Phase3_weights_hm3_no_MHC`

## Outputs

- `results/phase0_extension/step1_ldsc_screen/ldsc_rg_summary.tsv`
- `results/phase0_extension/step1_ldsc_screen/ldsc_rg_significant.tsv`
- `results/phase0_extension/step1_ldsc_screen/ldsc_qc.tsv`

## Threshold

LDSC-positive pairs are selected at:

`BH-FDR q < 0.05` and `rg > 0`

LDSC runs in parallel with LAVA. MTAG receives pairs positive in either LDSC or LAVA.
