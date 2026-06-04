# 08 TWAS

The main TWAS branch is MAGMA-candidate restricted:

1. Select MAGMA Bonferroni-positive genes.
2. Run FUSION TWAS only for those candidate genes.
3. Use GTEx v8 EUR 49-tissue weights.
4. Call TWAS support by BH-FDR within each `pair_id x tissue`.

## Main 49-Tissue Branch

Scripts are under `scripts/gtexv8_magma_candidates/`.

- `prepare_magma_positive_inputs.py`: writes MAGMA Bonferroni-positive pair-gene inputs.
- `prepare_branch_inputs.py`: audits FUSION/LDREF/weights and can normalize MTAG full sumstats.
- `download_gtexv8_eur_weights.sh`: downloads GTEx v8 EUR FUSION weights.
- `prepare_restricted_weights.py`: restricts each tissue weight `.pos` to MAGMA candidate genes.
- `run_all_fusion_candidate_twas.sh`: runs pair x tissue x chromosome FUSION jobs.
- `merge_fusion_outputs.py`: merges FUSION output, computes FDR, and collapses to gene-level support.

The older `step9_twas_pipeline.py` is retained as historical six-tissue full-MTAG TWAS code, not the main workflow.
