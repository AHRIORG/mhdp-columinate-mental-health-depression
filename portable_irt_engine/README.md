# portableIRTEngine

`portableIRTEngine` is the Objective 3 portable scoring scaffold for the
SSQ-10 IRT harmonization engine developed in the CO-LUMINATE psychometric
workspace.

## What It Does

- loads a selected packaged SSQ-10 harmonization engine
- validates incoming SSQ-10 data
- reconstructs the fixed-parameter IRT scoring model
- returns `Theta_Harmonized`
- maps theta back to PHQ-oriented quantities
- applies stored depression cutoffs
- exposes a Shiny interface for manual and batch scoring

## Privacy Design

The intended deployment model is privacy-first:

- uploaded data are processed in memory
- raw uploads are not written to persistent storage
- no raw scoring inputs are logged
- only user-requested downloads are materialized during the active session

## Current Engine Bundle

The scaffold now packages a wider Objective 3 engine family:

- `sens_at_spec` 1-factor engines for `overall`, `sex`, `age group`, and `sex + age group`
- `spec_at_sens` 1-factor engines for `overall`, `sex`, `age group`, and `sex + age group`
- `youden` 1-factor engines for `overall`, `sex`, `age group`, and `sex + age group`
- `model_tcc` 1-factor overall engine

The default production engine remains:

- `ssq10_sens_at_spec_1f_overall`

The bundle also retains audit-only multidimensional exports:

- `sens_at_spec` 2-factor engines
- `youden` 4-factor overall engine

These multidimensional exports currently carry flat SSQ item slopes in the
saved artifacts, so they are kept for provenance only and are hidden from live
scoring until the upstream Objective 3 exports are repaired.

## Public Release Readiness

The scaffold is designed to be mirrored into the public
`mhdp-columinate-mental-health-depression` repository without participant-level
data. The packaged engine objects do not contain raw records. The metadata
artifacts included in `inst/engines/` should be sanitized before public
release, and this repo includes a helper script for that purpose:

- `data-raw/sanitize_engine_metadata.R`

See `PUBLIC_RELEASE_MANIFEST.md` for the exact file-level release guidance.

This engine was copied from the Objective 3 winners directory so the package is
anchored to the current reporting workspace rather than to an external repo.

## Project Layout

- `R/`: package scoring functions
- `inst/engines/`: packaged engine objects and engine catalog
- `inst/extdata/`: example batch input
- `shiny-app/`: browser-based scoring interface
- `deploy/`: hosted deployment entrypoint and environment examples
- `tests/testthat/`: starter unit tests

## Hosted Deployment

The package now includes a container-ready deployment scaffold:

- `Dockerfile`
- `deploy/run-shiny.R`
- `deploy/.env.example`
- `DEPLOYMENT.md`

This supports hosting the Shiny workspace behind a real HTTPS URL rather than a
localhost-only runtime.

## Next Steps

1. add formal package documentation and examples
2. harden input validation for deployment edge cases
3. repair and regenerate the multidimensional Objective 3 exports upstream
4. deploy the Shiny app to a managed host and register the final HTTPS URL in the website config
