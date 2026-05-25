# Public Causal Mediation Bundle

This folder contains the public-safe script bundle for the CO-LUMINATE causal
mediation workflow. It is staged from the private analysis workspace but uses
repo-local paths and does not point back to private project folders.

## Scope

- `programs/analysis/build_strict_joint_history_baseline_imputations.R`: builds
  baseline-covariate imputations from the restricted subject-level analysis
  input and transition-history bundle.
- `programs/analysis/fit_joint_history_mediation.R`: fits total-effect,
  single-mediator, and joint-mediator models for a supplied branch.
- `programs/analysis/pool_joint_history_mediation_results.R`: pools
  imputation-specific mediation results.
- `programs/analysis/run_strict_joint_history_mediation_pipeline.R`: orchestrates
  the strict joint-history mediation pipeline.
- `config/`: public-safe contrast, adjustment-set, pathway-variable, and
  retained-solution metadata.

## Restricted Inputs

The scripts expect controlled-access inputs under `input/` inside this bundle:

- `input/jmaes_strict_subject_level_analysis.rds`
- `input/jmaes_strict_transition_bundle.rds`
- `input/imputations/strict_jmaes_final/jmaes_strict_subject_level_final_imp_*.rds`
  when running only the fitting or pooling stages.

These files are not included in the public repository. They may be supplied
through an approved data-access process.

## Public Outputs

Disclosure-reviewed aggregate outputs are staged separately under:

- `results/release/causal_mediation/` in the private handoff repo;
- `results/causal_mediation/` after approved public promotion.

The release outputs are aggregate tables only. Subject-level records,
imputation objects, model objects, and draw-level files are not public handoff
materials.

## Dependencies

The bundle expects the following R packages:

- `dplyr`
- `MASS`
- `mice`
- `purrr`
- `readr`
- `sandwich`
- `tibble`
- `tidyr`
