# 03 MTAG After LDSC/LAVA

This step selects trait pairs for MTAG after the parallel LDSC and LAVA screens.

## Selection Rule

A pair enters MTAG if:

`LDSC Bonferroni-positive OR LAVA Bonferroni-positive`

## Scripts

- `scripts/select_and_run_mtag_pairs.py`: writes selected pairs and MTAG commands; use `--run` to execute MTAG.
- `scripts/build_step2_input_tables.py`: retained input-table construction script from the previous bundle.

## MTAG SNP Threshold

Downstream MTAG SNP evidence should be retained at:

`P < 5e-8`
