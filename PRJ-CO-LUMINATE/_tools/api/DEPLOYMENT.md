# CO-LUMINATE API Deployment

This API can be deployed on any container platform that provides HTTPS.

## 1) Prerequisites

1. Ensure `irt_joint_models.rds` is available at:
   - `PRJ-CO-LUMINATE/_tools/_objects/irt_joint_models.rds`
2. Confirm local API works:
   - `Rscript _tools/api/run_api.R`

## 2) Build container image

From repository root:

```bash
docker build -f PRJ-CO-LUMINATE/_tools/api/Dockerfile -t columinate-irt-api .
```

## 3) Run locally (container)

```bash
docker run --rm -p 8080:8080 \
  --env-file PRJ-CO-LUMINATE/_tools/api/.env.example \
  columinate-irt-api
```

Health check:

```bash
curl -sS http://127.0.0.1:8080/health
```

## 4) Deploy to hosted HTTPS

Use the same image in your preferred host. Set these environment variables:

1. `COLUMINATE_API_HOST=0.0.0.0`
2. `COLUMINATE_API_PORT=8080` (or host-provided port)
3. `COLUMINATE_API_CORS_ORIGIN=https://mapel-lab.github.io` (or your site URL)
4. `COLUMINATE_API_KEY=<strong-secret-token>`
5. `COLUMINATE_SCORING_SCRIPT=/app/PRJ-CO-LUMINATE/_tools/api/scoring_runtime.R`

If the platform provides a dynamic `PORT`, map it to `COLUMINATE_API_PORT`.

## 5) Connect website Tools tab

1. Open Tools tab.
2. Set `API Base URL` to your HTTPS endpoint.
3. Enter API key in `API Key (Bearer token)` if enabled.
4. Click `Load Engines`.
