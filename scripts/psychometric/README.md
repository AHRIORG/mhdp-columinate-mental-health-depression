# Public Psychometric Script Bundle

This bundle contains public-safe psychometric analysis scripts copied from the private development workspace.

## Scope

- Includes active scripts only.
- Excludes archived/retired scripts.
- Uses repo-relative defaults for input data and environment-variable overrides.

## Data Inputs

By default, scripts expect:

- `data/inputs/dt_psychometric.RData`
- `data/inputs/dt_dreams_multi_ssq.RData` (only for optional sensitivity scoring)

You can override with environment variables:

- `COLU_PSYCHOMETRIC_RDATA`
- `COLU_DREAMS_RDATA`

## Notes

- This bundle is for reproducible method transfer and review.
- Do not add private data files to this public repository.
