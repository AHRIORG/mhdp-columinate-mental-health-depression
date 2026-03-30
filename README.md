# mhdp-columinate-mental-health-depression

Public publication repository for the CO-LUMINATE depression workstream within MHDP.

This repository is the approved release surface. Development happens in the private `CO-LUMINATE` repository, and only reviewed content should be promoted here.

Before promoting any content from the private repository, follow the release rules in `docs/release-rules.md`.

## Public Website Handoff

Publication rule: do not develop the public website directly in this public repository. Publish only approved rendered output.

Required handoff flow:

1. Render the website in private `website/`.
2. Review the output in private `website/`.
3. Copy the rendered files into public `website/`.
4. Run `scripts/release_preflight.sh --staged`.
5. Commit and push the public repository.
6. GitHub Pages deploys automatically via GitHub Actions.

## Public Website

GitHub Pages URL: https://ahriorg.github.io/mhdp-columinate-mental-health-depression/

## Intended Structure

- `docs/`: public-facing documentation and study notes.
- `data/raw/`: instructions for obtaining raw data through approved AHRI channels. No raw data is stored here.
- `data/data_examples/`: approved public-safe example data only.
- `scripts/`: public-safe analysis and publication scripts.
- `results/`: release-ready reports, figures, and tables.
- `metadata/`: codebooks, dictionaries, and metadata approved for release.
- `website/`: rendered static site payload that GitHub Pages deploys from this repository.
