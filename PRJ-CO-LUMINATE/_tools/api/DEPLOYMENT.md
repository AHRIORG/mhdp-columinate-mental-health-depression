# CO-LUMINATE API Deployment

This API can be deployed on any HTTPS container host. The simplest path for this repository is a Render web service backed by the repo-root `render.yaml`.

## 1) Prerequisites

1. Ensure `irt_joint_models.rds` is available at:
   - `PRJ-CO-LUMINATE/_tools/_objects/irt_joint_models.rds`
2. Confirm the local API works:
   - `cd PRJ-CO-LUMINATE`
   - `Rscript _tools/api/run_api.R`
3. Confirm the health endpoint:
   - `curl -sS http://127.0.0.1:8088/health`

## 2) Build container image locally

From repository root:

```bash
docker build -f PRJ-CO-LUMINATE/_tools/api/Dockerfile -t columinate-irt-api .
```

## 3) Run locally in Docker

```bash
docker run --rm -p 8080:8080 \
  --env-file PRJ-CO-LUMINATE/_tools/api/.env.example \
  columinate-irt-api
```

Health check:

```bash
curl -sS http://127.0.0.1:8080/health
```

## 4) Deploy on Render

1. Push the repository with the following deployment files:
   - `render.yaml`
   - `.dockerignore`
2. In Render, create a new Blueprint from the GitHub repository.
3. Render will build `PRJ-CO-LUMINATE/_tools/api/Dockerfile` using the repo root as Docker context.
4. Confirm the health check path is `/health`.
5. Recommended runtime environment variables:
   - `COLUMINATE_API_HOST=0.0.0.0`
   - `COLUMINATE_API_CORS_ORIGIN=https://mapel-lab.github.io`
   - `COLUMINATE_SCORING_SCRIPT=/app/PRJ-CO-LUMINATE/_tools/api/scoring_runtime.R`
6. Optional:
   - `COLUMINATE_API_KEY=<strong-secret-token>` if you want Bearer auth

Note:

- Render sets `PORT` automatically for web services, and `run_api.R` now falls back to that value.
- If `COLUMINATE_API_KEY` is omitted, the API stays public and the website's API key field can be left blank.

## 5) Wire the hosted API into the website

After Render creates the service URL:

1. Open `PRJ-CO-LUMINATE/api-config.json`
2. Set `api_base_url` to your HTTPS endpoint
3. Re-render the site and redeploy GitHub Pages

Example:

```json
{
  "api_base_url": "https://your-service.onrender.com"
}
```

## 6) Connect from the Tools tab

The Tools tab resolves the API URL in this order:

1. `?api=https://...` query string
2. Browser local storage
3. `api-config.json`
4. Local default `http://127.0.0.1:8088` when running on `localhost` / `127.0.0.1`

If the live site is loaded over HTTPS, the API must also use HTTPS.
