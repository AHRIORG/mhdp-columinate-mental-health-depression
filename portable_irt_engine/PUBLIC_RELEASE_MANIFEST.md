# Public Release Manifest

This document defines what should and should not be copied from
`portable_irt_engine/` into the public
`mhdp-columinate-mental-health-depression` repository.

## Safe To Publish

These files are intended to go public:

- `DESCRIPTION`
- `NAMESPACE`
- `LICENSE`
- `.Rbuildignore`
- `.gitignore`
- `README.md`
- `R/`
- `shiny-app/app.R`
- `tests/testthat.R`
- `tests/testthat/`
- `inst/engines/engine_catalog.csv`
- `inst/engines/*.rds`
  Note: only after metadata sanitization
- `inst/extdata/example_batch.csv`
  Note: keep only if it remains synthetic/demo-only
- `data-raw/sanitize_engine_metadata.R`
- `data-raw/stage_public_repo_bundle.R`

## Must Be Sanitized Before Public Release

These artifacts require sanitization:

- `inst/engines/*_meta.rds`
  Reason: the copied Objective 3 metadata originally contained absolute
  internal file paths in `source_file`.

The provided helper script sanitizes these metadata objects by:

- replacing absolute `source_file` paths with the source artifact basename
- adding `source_scope = "public_release_sanitized"`
- adding `source_note` to document that the path was intentionally removed

## Must Not Be Published

Do not push the following into the public repo as part of the app bundle:

- any raw participant data
- `statistical_analysis/`
- private output caches from Objective 3
- `.Rhistory`
- `.RData`
- `.Rproj.user/`
- `.quarto/`
- any future logs, temporary uploads, or ad hoc scored outputs

## Recommended Public Repo Target

Recommended public destination:

- `public/mhdp-columinate-mental-health-depression/portable_irt_engine/`

That keeps the package core and Shiny app self-contained and separate from the
website asset bundle.

## Current Engine Status

- `ssq10_sens_at_spec_1f_overall`
  Status: functional default for the public scaffold
- `ssq10_sens_at_spec_2f_overall`
  Status: retained for audit; appears to have flat SSQ slopes in the exported
  artifact and should not be treated as production-ready until repaired
