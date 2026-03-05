# CO-LUMINATE IRT API Contract

Version: `0.1.0`

## Purpose

Expose calibrated IRT scoring engines as an HTTP API for:

1. Engine discovery and metadata (`/engines`).
2. Single respondent scoring (`/score/single`).
3. Batch scoring (`/score/batch`).

The API uses the bundled `irt_joint_models.rds` scoring engines and returns:

1. `Theta_Harmonized`
2. `Depression_Binary`
3. `PHQ_Expected_Sum`
4. `PHQ_Prob_GE10`

## Authentication

Authentication is optional and controlled by environment variable:

1. If `COLUMINATE_API_KEY` is empty, scoring endpoints are open.
2. If `COLUMINATE_API_KEY` is set, send `Authorization: Bearer <token>` for:
   - `POST /score/single`
   - `POST /score/batch`

## Endpoints

### `GET /health`

Returns service health and loaded engine source metadata.

### `GET /engines`

Returns available engines and key metadata:

1. `engine_id`
2. `items`
3. `n_factors`
4. `cutoff_apply_mode`
5. `overall_cutoff_theta`

### `GET /contract`

Returns a machine-readable version of this contract.

### `POST /score/single`

Scores one respondent.

Request JSON:

```json
{
  "engine_id": "TCC 1-Factor Model (Overall)",
  "responses": {
    "SSQ01": 1, "SSQ02": 0, "SSQ03": 1, "SSQ08": 0, "SSQ09": 1,
    "SSQ10": 0, "SSQ11": 1, "SSQ12": 0, "SSQ13": 1, "SSQ14": 0
  },
  "sex": "Female",
  "age": 19
}
```

Response JSON:

```json
{
  "engine_id": "TCC 1-Factor Model (Overall)",
  "output": {
    "Theta_Harmonized": 1.23,
    "Depression_Binary": 0,
    "PHQ_Expected_Sum": 6.10,
    "PHQ_Prob_GE10": 0.22
  }
}
```

### `POST /score/batch`

Scores multiple respondents.

Request JSON:

```json
{
  "engine_id": "TCC 1-Factor Model (Overall)",
  "records": [
    {
      "SSQ01": 1, "SSQ02": 0, "SSQ03": 1, "SSQ08": 0, "SSQ09": 1,
      "SSQ10": 0, "SSQ11": 1, "SSQ12": 0, "SSQ13": 1, "SSQ14": 0,
      "SEX": "Female", "AGE": 19
    },
    {
      "SSQ01": 0, "SSQ02": 0, "SSQ03": 0, "SSQ08": 0, "SSQ09": 0,
      "SSQ10": 0, "SSQ11": 0, "SSQ12": 0, "SSQ13": 0, "SSQ14": 0,
      "SEX": "Male", "AGE": 22
    }
  ]
}
```

Response JSON:

1. `engine_id`
2. `n_records`
3. `output` (row-wise scored data frame)

## Validation Rules

1. All 10 SSQ items are required per record.
2. `engine_id` must exist in `/engines`.
3. Batch size is limited by `COLUMINATE_API_MAX_BATCH` (default `5000`).
4. Rows with all SSQ items missing produce missing outputs.

## CORS

`Access-Control-Allow-Origin` is set from:

1. `COLUMINATE_API_CORS_ORIGIN` (recommended in production).
2. `*` by default.

## Local Run

From `PRJ-CO-LUMINATE`:

```bash
Rscript _tools/_scripts/assess_irt_api_model.R
Rscript _tools/api/run_api.R
```

Default URL: `http://127.0.0.1:8088`
