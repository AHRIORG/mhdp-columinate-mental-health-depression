# Publishing And Release Flow

Follow the rules in [release-rules.md](release-rules.md) before transferring any content from the private `CO-LUMINATE` repository into this public repository.

## General Release Steps

1. Select only approved files from release-safe locations in the private repository.
2. Run the path, object-type, and manual review checks described in [release-rules.md](release-rules.md).
3. Copy the approved files into the matching location in this public repository.
4. Review the final `git diff` in this repository before committing.
5. Commit and push only after the release review is complete.

## Website Publishing Flow

The public website is deployed from the committed contents of `website/`.

## Source of Truth

- Website development happens in the private `CO-LUMINATE` repository.
- The public repository stores only approved rendered output and other release-safe assets.

## Public Website Handoff (Required)

Publication rule: do not develop the public website directly in the public repository. Publish only approved rendered output.

Automated consistency guard: the `Release Guard` GitHub Actions workflow runs on `push` and `pull_request` to `main` and fails when blocked file types, `_quarantine` paths, `.DS_Store`, secret-like tokens, or private local path markers are detected.

Before commit/push, run a local preflight scan in this public repo:

```sh
scripts/release_preflight.sh --staged
```

## Publication Steps

1. Render the website in private `website/`.
2. Review the output in private `website/`.
3. Run the release checks in [release-rules.md](release-rules.md), paying special attention to rendered HTML, metadata comments, and `search.json`.
4. Copy the rendered files into public `website/`.
5. Commit and push the public repository.
6. GitHub Actions deploys the site to GitHub Pages automatically.
