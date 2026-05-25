# Public Pathway Screening Bundle

This folder contains the public-safe script bundle for the CO-LUMINATE pathway
screening and behavioral mediator prioritisation workflow.

## Scope

- `programs/analysis/build_contextual_inventory_tables.R`: converts a
  Markdown pathway inventory into release-ready long and summary tables.
- `programs/analysis/sanitize_pathway_screening_outputs.R`: stages aggregate
  pathway-screening outputs from a controlled analysis run into a public-safe
  release folder.
- `programs/utils/pathway_screening_helpers.R`: bundle-local path and table
  helpers.

## Restricted Inputs

The scripts can run against local controlled-access inputs placed under:

- `data/inputs/pathway_inventory.md`
- `data/inputs/pathway_screening_raw/`

These inputs are not included in the public repository unless explicitly
approved. The public release should contain aggregate tables and figures only.

## Public Outputs

Disclosure-reviewed aggregate outputs are staged separately under:

- `results/release/pathway_screening/` in the private handoff repo;
- `results/pathway_screening/` after approved public promotion.

Do not publish subject-level data, fitted model objects, draw-level files, or
files containing local/private paths.

## Dependencies

The bundle expects the following R packages:

- `dplyr`
- `glue`
- `purrr`
- `readr`
- `stringr`
- `tibble`
- `tidyr`
