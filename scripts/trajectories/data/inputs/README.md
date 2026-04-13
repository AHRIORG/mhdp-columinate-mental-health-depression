# Trajectory Input Files

Place the public-safe trajectory input `.RData` files in this folder before running the trajectory bundle.

Expected defaults:

- `dt_childhood_exposure.RData`
- `dt_multistudies.RData`

Expected objects inside those files:

- `dt_childhood_exposure`
- `dt_multistudies`

If you store the inputs elsewhere, set:

- `COLU_CHILDHOOD_RDATA`
- `COLU_MULTISTUDIES_RDATA`

These repo-local inputs replace the old private internal lookup pattern.
