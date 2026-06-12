# 03 MTAG After LDSC/LAVA

This step selects trait pairs for MTAG after the parallel LDSC and LAVA screens.

## Selection Rule

A pair enters MTAG if:

`LDSC BH-FDR-positive OR LAVA adjusted-P-positive`

## Scripts

- `scripts/prepare_inputs.py`: constructs the pair and candidate input tables.
- `scripts/run_mtag.py`: writes selected pairs and MTAG commands; use `--run` to execute MTAG.

## MTAG SNP Threshold

Downstream MTAG SNP evidence should be retained at:

`P < 5e-8`
