# Reproducibility notes

This repository stores code and table contracts only. Data and references must
be acquired by the user and kept outside git.

## Compute policy

Large analyses should be submitted through an HPC scheduler or workflow manager.
Scheduler-specific templates are intentionally not included. Wrap the scripts
with the scheduler used by your institution.

## Main safeguards

- Do not run large GWAS-wide R/Python jobs on login nodes.
- Do not commit raw GWAS, LD reference panels, or generated full genome-wide
  result files.
- Do not commit SSH keys, tokens, passwords, or private server paths.
- Keep non-significant genome-wide PLACO/CPASSOC temporary rows in scratch and
  delete them after extracting positive rows.

## Positive-result thresholds

Genome-wide PLACO/CPASSOC sensitivity keeps rows meeting either:

- `P < 5e-8`
- `BH-FDR < 0.05`

These two classes should be reported separately in QC summaries when needed.

