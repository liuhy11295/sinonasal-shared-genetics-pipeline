# Nasal inflammatory disease shared genetics pipeline

This repository contains the code used to prepare and analyze pairwise shared
genetic evidence between nasal inflammatory diseases and immune/respiratory
partner traits.

The repository contains analysis code and the derived source-data tables needed
to audit the five final supplementary tables. It does not include raw GWAS
summary statistics, LD reference panels, credentials, logs, or private paths.

## Final Workflow

The numbered folders follow the intended final analysis order:

1. `00_environment_setup/`: environment setup, method installation, and input manifest scripts.
2. `01_ldsc_screen/`: LDSC genetic-correlation runner and BH-FDR screening.
3. `02_lava_screen/`: LAVA local-rg runner and aggregation across all pairs; this is parallel to LDSC, not downstream of LDSC.
4. `03_mtag_after_ldsc_lava/`: pair/input table construction and MTAG selection for pairs positive in LDSC or LAVA.
5. `04_screen_placo_cpassoc_genomewide/`: genome-wide PLACO and CPASSOC screening; positive SNP evidence is retained at `P < 5e-8`.
6. `05_coloc_abf/`: ABF coloc on screened candidate loci before SuSiE.
7. `06_coloc_susie/`: coloc-positive loci enter coloc-SuSiE fine-mapping and aggregation.
8. `07_magma_gene_pathway/`: MAGMA gene and pathway analysis after cross-trait and coloc-SuSiE support.
9. `08_twas/`: MAGMA-positive candidate genes enter GTEx v8 49-tissue FUSION TWAS; TWAS support is called by FDR.
10. `09_lcv_mr_secondary/`: secondary LCV and Mendelian-randomization evidence scripts.
11. `10_downstream_annotation/`: ORA, FUMA-concordance, cell-marker, and exploratory drug-target annotation.
12. `11_publication_outputs/`: final figure/table scripts and the nine source-data sheets underlying five supplementary tables.

In short:

```text
LDSC runner  ┐
             ├─> LDSC-positive or LAVA-positive pairs -> MTAG
LAVA runner  ┘

MTAG / PLACO / CPASSOC positive SNP evidence
-> coloc ABF
-> coloc-SuSiE
-> MAGMA Bonferroni-positive genes
-> GTEx v8 49-tissue TWAS
```

## Thresholds

- LDSC: Benjamini-Hochberg FDR across tested disease pairs, `q < 0.05`, with `rg > 0`.
- LAVA: positive local genetic correlation with the supplied adjusted probability `p_adj <= 0.05`.
- MTAG: genome-wide SNP evidence at `P < 5e-8`.
- PLACO: genome-wide SNP evidence at `PLACO_p < 5e-8`.
- CPASSOC: genome-wide SNP evidence at `SHet_p < 5e-8` or `SHom_p < 5e-8`.
- coloc ABF and coloc-SuSiE: high-confidence evidence requires `PP.H4 >= 0.80` and, where `PP.H3` is available, `PP.H4 > PP.H3`.
- MAGMA: Bonferroni-positive genes are used for downstream TWAS candidate selection.
- TWAS: FDR-based tissue-level support across the GTEx v8 49-tissue branch.

## Directory Layout

```text
00_environment_setup/                 Environment and manifest setup
01_ldsc_screen/                       LDSC rg screening
02_lava_screen/                       LAVA local-rg screening
03_mtag_after_ldsc_lava/              MTAG after LDSC/LAVA positive pair selection
04_screen_placo_cpassoc_genomewide/   Genome-wide PLACO and CPASSOC screening
05_coloc_abf/                         coloc ABF before SuSiE
06_coloc_susie/                       coloc-SuSiE
07_magma_gene_pathway/                MAGMA gene/pathway analysis
08_twas/                              TWAS, including GTEx v8 49-tissue candidate branch
09_lcv_mr_secondary/                  Secondary LCV/MR scripts
10_downstream_annotation/             Exploratory downstream annotations
11_publication_outputs/               Figure/table code and final table source data
config/                               Example path configuration
tests/                                Threshold and repository-safety tests
```

## Environment

Python 3.10+ is recommended. R 4.3+ is recommended for the R scripts. The
publication build was audited with Python 3.13 and R 4.5.3; method-specific
legacy tools may require their own isolated environments.

Python packages used by different steps include:

- pandas
- numpy
- scipy
- statsmodels
- pyarrow
- polars
- tqdm
- pyyaml
- joblib

R packages used by different steps include:

- data.table
- dplyr
- readr
- tidyr
- stringr
- ggplot2
- optparse
- argparse
- coloc
- susieR
- Matrix
- foreach
- doParallel

Third-party methods such as LDSC, MTAG, LAVA, PLACO, CPASSOC, MAGMA, LCV, and
TwoSampleMR are not vendored here. Install them separately in your analysis
environment and point the configuration or wrappers to their locations.

## Minimal Setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Configuration

Copy `config/example_paths.yaml` to a private location and edit paths for your
server or workstation. Do not commit real private paths, raw data locations, or
credentials.

The main runtime variables are documented in `.env.example`. Scripts fail
early when required roots are missing instead of falling back to author-local
paths.

## Script organization

Each numbered workflow directory contains only its final reproducible scripts.
Names follow a small stable vocabulary: `prepare_*`, `run_*`,
`aggregate_results`, and `finalize_outputs`. Historical benchmarks, repair
utilities, previews, and superseded branches are intentionally excluded; Git
history remains the archive.

## Publication Outputs

`11_publication_outputs/tables/source_data/` contains nine machine-readable
sheets grouped into five final supplementary tables:

- ST1: GWAS metadata and power QC.
- ST2: analysis parameters and harmonization QC.
- ST3: the final 31 disease pairs.
- ST4: prioritized loci and prioritized genes.
- ST5: significant pathways, approved-drug target annotations, and secondary directional evidence.

Drug-target results are exploratory annotations and do not imply therapeutic
efficacy. LCV and bidirectional MR are secondary directional analyses rather
than definitive causal proof.

## Verification

```bash
python -m unittest discover -s tests -v
python -m compileall -q .
find . -type f \( -name '*.sh' -o -name '*.slurm' \) -exec bash -n {} \;
```

## Excluded Data

The repository excludes raw and generated data:

- raw GWAS files
- 1000 Genomes or LD score reference files
- MAGMA/MTAG/LDSC third-party source trees
- large intermediate result tables (the compact final supplementary-table source data are included)
- server logs
- private keys or SSH configuration
- compressed archives and large binary reference files

Scheduler scripts that are part of the reproducible code bundle are tracked in
the numbered folders. Runtime scheduler outputs such as `.out`, `.err`, `.log`,
and `slurm-*.out` remain excluded.
