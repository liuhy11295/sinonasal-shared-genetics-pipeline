# Repository Script Cleanup Implementation Plan

**Goal:** Keep only the final reproducible analysis and publication scripts, with stable names and no development, benchmark, repair, preview, or duplicate entry points.

**Architecture:** Each numbered analysis directory retains one runner per distinct final method plus one aggregator when required. Final TWAS and LCV/MR implementations replace superseded branches. Publication figures retain four figure renderers, one supplementary renderer, one shared helper, and one all-figures runner.

**Tech Stack:** Python, R, Bash/SLURM, unittest, Git.

---

### Task 1: Consolidate publication figures

- Merge the Figure 4 component and directional-panel code into `render_figure4.R`.
- Delete exploratory and standalone QA scripts.
- Retain exactly seven R files in the figure script directory.
- Update the public copy and inventory tests.

### Task 2: Remove duplicate analysis branches

- Merge the alternate LAVA table-selection behavior into the retained aggregator.
- Remove one-time PLACO/CPASSOC repair, trim, and benchmark scripts.
- Remove the superseded monolithic MAGMA runner.
- Retain the final GTEx v8 candidate TWAS branch and remove old step9/phase0 wrappers.
- Move the final 31-pair LCV/MR scripts into the primary scripts directory and remove the old step4 branch.

### Task 3: Simplify publication tables and documentation

- Retain the final table builder as the formal table script.
- Remove standalone table audit wrappers now covered by repository tests.
- Update directory READMEs and the root workflow inventory.

### Task 4: Verify

- Assert exact script inventories for every workflow directory.
- Parse all R scripts, compile all Python scripts, and check all shell scripts.
- Run repository tests and full figure rendering/QA.
- Confirm no internal reference points to a removed filename.
- Commit and push to `main`.
