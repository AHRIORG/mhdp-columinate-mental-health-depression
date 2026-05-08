# Hosted Deployment Guide

`portableIRTEngine` is designed to be deployed as a hosted Shiny application
behind a real HTTPS URL. GitHub can host the code, but it cannot run the Shiny
runtime itself.

## Deployment Shape

Use two layers:

1. this repository as the versioned code and engine bundle source
2. a managed container or Shiny host to run the live app

The website should then point to the hosted HTTPS URL rather than to
`127.0.0.1`.

## Files Included For Hosting

- `Dockerfile`
- `deploy/run-shiny.R`
- `shiny-app/`
- `R/`
- `inst/engines/`
- `inst/extdata/`

## Environment Variables

The container honours the following environment variables:

- `PORT`
  Hosting platforms commonly set this automatically.
- `PORTABLE_IRT_HOST`
  Default: `0.0.0.0`
- `PORTABLE_IRT_PORT`
  Default: `8080`

In most hosted environments, `PORT` is sufficient.

## Build Locally

From the `portable_irt_engine/` directory:

```bash
docker build -t portable-irt-engine .
```

## Run Locally As A Container

```bash
docker run --rm -p 8080:8080 portable-irt-engine
```

Then open:

```text
http://127.0.0.1:8080/
```

## Recommended Hosting Pattern

Deploy the container to a platform that provides:

- HTTPS
- a stable public URL
- container support
- environment-variable based port binding

Examples include container hosts such as Render, Fly.io, Cloud Run, or similar
services.

## Website Connection

After the app is hosted, update the website configuration file:

`private/CO-LUMINATE/website/api-config.json`

Set:

```json
{
  "api_base_url": "",
  "irt_scoring_url": "https://your-hosted-irt-app.example.org"
}
```

The website tool page will then use that hosted URL for:

- the launch button
- the embedded workspace iframe

## Privacy Reminder

The app is designed to process uploads in memory and return outputs directly to
the user. Hosted deployment still requires careful platform settings so that raw
uploads are not logged or persisted outside the user-requested download flow.
