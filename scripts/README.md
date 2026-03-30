# Public Scripts

Place public-safe scripts here. Any script that depends on private data locations, internal credentials, or unpublished methods must remain in the private repository.

## Release Helpers

- `release_preflight.sh`: local release scan for staged files (or provided paths) before pushing to public.

## Psychometric Scripts

- `psychometric/`: public-safe psychometric script bundle.
- `psychometric/programs/`: active analysis and helper scripts copied from private development (retired scripts excluded).
- Scripts use repo-local defaults (`data/inputs/*.RData`) and support overrides via `COLU_PSYCHOMETRIC_RDATA` and `COLU_DREAMS_RDATA`.

## Data Management Assets

- `data_management/`: retained data-management attributes from the previous OBJ00 structure, without the `obj00/` folder wrapper.
