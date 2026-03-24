# Website Publishing Flow

The public website is deployed from the committed contents of `website/`.

## Source of Truth

- Website development happens in the private `CO-LUMINATE` repository.
- The public repository stores only approved rendered output and other release-safe assets.

## Publication Steps

1. Render the approved website in the private repository.
2. Review the rendered output internally.
3. Copy the rendered files from private `website/_site/` into public `website/`.
4. Commit and push the changes to `main`.
5. GitHub Actions deploys the site to GitHub Pages.
