# Public Scripts

Place public-safe scripts here. Any script that depends on private data locations, internal credentials, or unpublished methods must remain in the private repository.

## Release Helpers

- `release_preflight.sh`: local release scan for staged files (or provided paths) before pushing to public.

## OBJ00 Rebuild Pipeline

- `obj00_rebuild_pipeline.R`: public-safe rebuild runner for OBJ00 analysis datasets. Users must plug in their own approved private data path via `--private-data-dir` (or `CO_LUMINATE_PRIVATE_DATA_DIR`).
- `obj00/`: minimal script bundle used by the rebuild runner (`01b` builder, metadata helper, IRT scoring helper, and approved scoring engine `.rds` needed for final-step rebuild mode).
