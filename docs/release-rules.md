# Public Release Rules

This repository is the public release surface for the CO-LUMINATE depression workstream. Development happens in the private `CO-LUMINATE` repository, and promotion into this repository should happen only after a deliberate release review.

The goal is to support open science by sharing reusable code, metadata, documentation, and release-ready outputs without exposing sensitive data or internal working details.

## Core Principles

- Keep anything private if there is uncertainty about whether it is safe to publish.
- Share methods, code, metadata, and approved outputs in ways that others can inspect and reuse.
- Never publish material that reveals participant-level data, internal storage locations, or machine-specific environments.
- Prefer small, reviewable releases over large directory copies.

## Preferred Source Locations In The Private Repository

Promotion into this public repository should normally start from private locations already intended for release handoff:

- `results/release/`
- `metadata/`
- selected files from `scripts/`
- approved rendered website output from `website/`

Do not copy directly from the following without an explicit release review:

- `data/raw/`
- `data/processed/`
- `data/access/`
- `_quarantine/`
- exploratory notebooks, scratch folders, local caches, or objective-specific work areas that have not been cleared for public release

## Never Publish

The following material must not be transferred into this public repository:

- raw, processed, or participant-level data
- row-level household or clinical extracts, even when they look small
- small-cell summaries or other outputs that could support re-identification
- serialized or binary objects that can carry data or model state, including `.RData`, `.rda`, `.rds`, `.qs`, `.fst`, `.feather`, `.parquet`, `.pkl`, `.joblib`, and `.sqlite`
- internal access notes, signed URLs, secrets, API keys, tokens, passwords, or environment files
- absolute or machine-specific paths outside this repository, including patterns such as `/Users/...`, `C:\...`, `OneDrive...`, network shares, or `_private_use/...`
- rendered notebooks, preview files, logs, comments, or metadata that expose local usernames, absolute paths, or internal filenames
- drafts, internal review artefacts, or unpublished materials that have not been approved for external sharing

## Allowed After Review

The following are normally appropriate for the public repository once they have been checked:

- scripts that use repo-relative paths and document how restricted inputs can be requested
- synthetic, fully anonymised, or aggregated example data approved for public release
- codebooks, dictionaries, variable inventories, and release notes
- final tables, figures, and reports approved for dissemination
- rendered website files after checking that the rendered output does not leak local paths or hidden source metadata

## Required Sanitisation Before Transfer

Before copying any file from the private repository into this one:

- replace absolute paths with repo-relative paths or documented placeholders
- remove or rewrite code that depends on private objects stored outside the public repository
- strip execution artefacts that can reveal local state, private directories, or hidden comments
- prefer text-based, reviewable formats over opaque binary artefacts
- add or update a `README.md` when public users need provenance, scope, or access-request guidance

## Validation Checklist

Every transfer should pass both automated scans and manual review. The scans below are heuristics only; passing them does not replace judgement.

1. Define the candidate files and their destination paths in the public repository.
2. Scan candidate text files for leaked paths and private markers.
3. Scan candidate files for banned object types.
4. Manually review any tabular files to confirm they are synthetic, anonymised, or appropriately aggregated.
5. Inspect rendered website and documentation output for embedded comments, search indexes, or metadata that reveal local paths.
6. Copy only the approved files into the matching public location.
7. Review the final `git diff` in the public repository before committing or pushing.

Suggested preflight commands, run from the repository root:

```sh
scripts/release_preflight.sh <candidate-paths>
rg -n '/Users/|OneDrive|_private_use|C:\\|file://' <candidate-paths>
rg --files <candidate-paths> | rg '\.(RData|rda|rds|qs|fst|feather|parquet|pkl|joblib|sqlite)$'
rg --files <candidate-paths> | rg '\.(csv|tsv|xlsx|xls|sav|dta)$'
rg -n 'readRDS|saveRDS|load\\(|read_csv|fread\\(|read_excel|read_dta|read_sav' <candidate-paths>
```

If the release includes rendered website assets, run an additional scan over the rendered output:

```sh
rg -n '/Users/|OneDrive|_private_use|file://|nb-[0-9]+' website
```

## Promotion Workflow

Use the following workflow when moving approved content from `private/CO-LUMINATE` into this repository:

1. Prepare and review the candidate files in the private repository.
2. Confirm that each file belongs in one of the public folders such as `docs/`, `scripts/`, `results/`, `metadata/`, `data/data_examples/`, or `website/`.
3. Run the validation checklist and fix any issues before copying files across.
4. Copy the approved files into the public repository.
5. Review the public diff as if it will be read by an external collaborator with no access to the private repo.
6. Commit and push only after the release review is complete.

## Rule Of Thumb

Every public commit should be explainable to an external reader and safe to share without private context. If a file requires private data, private paths, or internal explanation to be safe, it is not ready for this repository.
