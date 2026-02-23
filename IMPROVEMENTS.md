# Improvements and Consistency Fixes

This repository copy was improved without changing the original source repository.

## Inconsistencies identified
- Project folder names were inconsistent (`TWASTutorial`, `TWAStutorial`, `TWAS_tutorial`).
- Several commands assumed folders existed but did not create them with `mkdir -p`.
- `kallisto` module command used `ml load` instead of `ml` (inconsistent with other sections).
- `samples.trimmed.txt` creation was shown as a commented line but later required by Slurm.
- R plotting section used `fread()` without loading `data.table`.
- Multiple typos and wording issues reduced clarity (for example: "cromosomes", "specie", "copy pase").
- `downloadSamples.sh` lacked a shebang, error handling, and output-directory creation.

## Improvements applied
- Rewrote `README.md` for consistency, clearer sequencing, and reproducibility.
- Standardized working directory naming (`TWAS_tutorial`) and path usage.
- Added explicit `mkdir -p` and setup/verification commands.
- Corrected module-loading and script snippets for FastQC, Trimmomatic, and Kallisto.
- Added missing `library(data.table)` in the Manhattan plot section.
- Modernized `downloadSamples.sh` with:
  - `#!/usr/bin/env bash`
  - `set -euo pipefail`
  - automatic output directory creation
  - loop-based URL download for easier maintenance

## Scope
- No biological assumptions or model choices were changed.
- No original data files were modified.
