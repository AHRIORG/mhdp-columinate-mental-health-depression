# OBJ00 Public Rebuild Bundle

This bundle allows approved users to rebuild OBJ00 analysis datasets locally by supplying their own private data directory.

No private paths are hardcoded. You must provide:

- `--private-data-dir=/path/to/co-luminate`  
  or set environment variable `CO_LUMINATE_PRIVATE_DATA_DIR`.

## Main Entry Point

Run from the public repository root:

```sh
Rscript scripts/obj00_rebuild_pipeline.R \
  --private-data-dir=/path/to/co-luminate \
  --include-final-steps=false \
  --build-examples=true \
  --self-repro-check=true
```

## Outputs

Default outputs are repo-local:

- synthetic examples: `data/data_examples/`
- validation report: `results/validation/obj00_pipeline_validation_report.md`
- optional ADaM save path: `results/intermediate/adam/`

## Final-Step Mode

If you enable `--include-final-steps=true`, the default scoring engine is:

- `scripts/obj00/data_management/data_examples/1_staging_snippets/derived_models/irt_joint_models.rds`

You can still override it with:

- `--score-engines-rds=/path/to/irt_joint_models.rds`

Optional:

- `--scoring-script=/path/to/02b_irt_model.R`

If the default scoring-engine file is missing, provide `--score-engines-rds` explicitly.
