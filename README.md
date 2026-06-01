# Nasal inflammatory disease shared genetics pipeline

This repository contains the code used to prepare and analyze pairwise shared
genetic evidence between nasal inflammatory diseases and immune/respiratory
partner traits.

The repository is intentionally code-only. It does not include GWAS summary
statistics, LD reference panels, private server paths, credentials, logs, or
large generated result files.

## Analysis scope

The final analysis workflow supports:

- standardized GWAS and LDSC munged summary statistics indexing
- LDSC and LAVA result table integration
- MTAG result filtering and full MTAG downstream inputs
- coloc and SuSiE-coloc result aggregation
- LCV and Mendelian randomization summaries
- candidate-based PLACO/CPASSOC
- genome-wide PLACO/CPASSOC sensitivity analysis from paired munged GWAS
- MAGMA gene-based analysis and GO/KEGG/Reactome pathway summaries
- optional TWAS post-processing
- final result-folder and visualization-input table construction

## What is excluded

No scheduler templates are included. In particular, this export intentionally
omits all `.slurm` files. The scripts are scheduler-agnostic and can be called
from any cluster scheduler, local wrapper, or workflow manager as long as large
computations are submitted according to the user's compute policy.

Also excluded:

- raw GWAS files
- 1000 Genomes or LD score reference files
- MAGMA/MTAG/LDSC third-party source trees
- generated result tables
- server logs
- private keys or SSH configuration
- compressed archives

## Directory layout

```text
src/
  phase0_pipeline/                 Core package builders and validators
  pipeline_builders/               Nasal4-vs-other package construction utilities
  step3_susie_coloc/               Fine-mapping and coloc.SuSiE helper scripts
  step4_lcv_mr/                    LCV and MR helper scripts
  step5_candidate_placo_cpassoc/   Candidate-based PLACO/CPASSOC scripts
  genomewide_placo_cpassoc/        Genome-wide PLACO/CPASSOC sensitivity scripts
  step8_mtag_magma_gsea/           Full MTAG -> MAGMA -> GSEA post-processing
  step9_twas_optional/             Optional TWAS post-processing script
tests/                             Unit tests for package-generation code
tools/                             Publication safety audit
config/                            Example path configuration
docs/                              Reproducibility and result-contract notes
```

## Environment

Python 3.10+ is recommended. R 4.3+ is recommended for the R scripts.

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
environment and point the configuration to their locations.

## Minimal setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m pytest
```

## Configuration

Copy `config/example_paths.yaml` to a private location and edit paths for your
server or workstation. Do not commit real private paths, raw data locations, or
credentials.

## Typical workflow

The exact command sequence depends on your local manifests and compute
environment, but the final pipeline is organized as:

1. Build or validate standardized GWAS and munged-GWAS manifests.
2. Build pair-level input tables.
3. Run coloc and SuSiE-coloc per locus.
4. Run LCV and MR summaries on priority pairs.
5. Run candidate-based PLACO/CPASSOC where needed.
6. Run genome-wide PLACO/CPASSOC sensitivity from paired munged GWAS.
7. Run full MTAG, MAGMA gene-based analysis, and GSEA/ORA summaries.
8. Build final result folders and visualization input tables.

For the genome-wide PLACO/CPASSOC sensitivity step, use:

```bash
Rscript src/genomewide_placo_cpassoc/run_genomewide_placo_cpassoc_pair.R \
  /path/to/project_root PAIR_ID

python src/genomewide_placo_cpassoc/aggregate_genomewide_placo_cpassoc.py
```

Only significant rows are retained by the genome-wide PLACO/CPASSOC scripts:
`P < 5e-8` or `BH-FDR < 0.05`. Full non-significant genome-wide outputs are not
kept.

## Publication safety audit

Before pushing to GitHub, run:

```bash
python tools/audit_github_code_export.py .
```

The audit blocks common data files, archives, logs, secrets, and unexpectedly
large files.

