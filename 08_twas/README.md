# 08 TWAS

The main TWAS branch is MAGMA-candidate restricted:

1. Select MAGMA Bonferroni-positive genes.
2. Run FUSION TWAS only for those candidate genes.
3. Use GTEx v8 EUR 49-tissue weights.
4. Call TWAS support by BH-FDR within each `pair_id x tissue`.

## Final scripts

- `scripts/select_candidate_genes.py`: writes MAGMA Bonferroni-positive pair-gene inputs.
- `scripts/prepare_sumstats.py`: audits references and normalizes MTAG summary statistics.
- `scripts/download_weights.sh`: downloads GTEx v8 EUR FUSION weights.
- `scripts/prepare_weights.py`: restricts tissue weights to candidate genes.
- `scripts/run_twas.sh`: runs pair-by-tissue-by-chromosome jobs.
- `scripts/run_twas_job.sh`: executes and monitors one FUSION job.
- `scripts/aggregate_results.py`: computes FDR and collapses output to gene-level support.
