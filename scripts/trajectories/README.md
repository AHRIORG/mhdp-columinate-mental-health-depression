# Public Trajectory Bundle

This folder contains the public-safe script bundle for the CO-LUMINATE social-determinant trajectory workflow.

## Scope

- `programs/analysis/run_lcga.R`: orchestrates univariate and joint LCGA fitting.
- `programs/analysis/report_univariate_lcga.R`: builds first-pass univariate reporting outputs.
- `programs/analysis/report_joint_lcga.R`: builds first-pass joint reporting outputs.
- `programs/utils/`: bundle-local helpers sourced by the analysis scripts.
- `data/metadata/class_registry.csv`: materialized class-label registry used during reporting.

## Repo-Local Inputs

By default, the bundle expects these repo-local inputs:

- `data/inputs/dt_childhood_exposure.RData`
- `data/inputs/dt_multistudies.RData`

You can override those paths at runtime with:

- `COLU_CHILDHOOD_RDATA`
- `COLU_MULTISTUDIES_RDATA`

## Dependencies

The bundle expects the following R packages to be available:

- `dplyr`
- `ggplot2`
- `glue`
- `gt`
- `here`
- `MplusAutomation`
- `purrr`
- `readr`
- `rlang`
- `scales`
- `stringr`
- `tibble`
- `tidyr`

`run_lcga.R` also assumes a working local Mplus installation.

## Public-Safety Notes

- The bundle no longer resolves data from private internal handoff paths.
- Prior class-label curation is materialized inside `data/metadata/class_registry.csv`, so the public bundle does not point back to unpublished internal project names or source scripts.
- The targeted joint baseline script was intentionally excluded because it depended on a non-public baseline results flow.
