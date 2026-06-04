# Final Code Bundle

This directory collects the final reproducible scripts retained after cleanup. Original result files remain in `results/` and `results_end/`. Script contents are preserved; only the bundle folders are ordered and named to reflect the intended analysis flow.

## Folder order

1. `00_environment_setup/`: environment setup, tool repair, and symlink manifest scripts.
2. `01_ldsc_screen/`: LDSC rg runner and Bonferroni screening.
3. `02_lava_screen/`: LAVA local-rg runner/aggregation for all pairs, parallel to LDSC, with Bonferroni screening.
4. `03_mtag_after_ldsc_lava/`: pair/input tables plus MTAG selection/run scripts for pairs positive in LDSC or LAVA.
5. `04_screen_placo_cpassoc_genomewide/`: genome-wide PLACO/CPASSOC scripts retaining P < 5e-8 SNPs.
6. `05_coloc_abf/`: ABF coloc on screened candidate loci; writes coloc-positive loci for SuSiE.
7. `06_coloc_susie/`: coloc-positive locus SuSiE fine-mapping and coloc-SuSiE aggregation.
8. `07_magma_gene_pathway/`: MAGMA gene/pathway analysis after cross-trait and coloc-SuSiE support.
9. `08_twas/`: MAGMA-positive candidate-gene TWAS, including the GTEx v8 49-tissue FUSION branch.
10. `09_lcv_mr_secondary/`: secondary LCV/MR evidence scripts.

The intended interpretation is: LDSC and LAVA run in parallel; pairs positive in either screen enter MTAG; MTAG, PLACO, and CPASSOC retain genome-wide SNP evidence at P < 5e-8; screened loci proceed to ABF coloc and coloc-SuSiE; MAGMA uses Bonferroni-positive gene evidence; MAGMA-positive genes enter 49-tissue GTEx v8 FUSION TWAS with TWAS support called by FDR.

The bundle intentionally excludes backup scripts, failed test scripts, Python bytecode caches, and pilot TWAS-GSEA variants.
