# Browser-Only IRT Scoring Architecture

This document describes the GitHub Pages path for the CO-LUMINATE portable
SSQ-10 IRT engine.

## Why A Browser Build Is Feasible

The current production engine family is well suited to client-side execution:

- production scoring is based on **1-factor** SSQ-10 engines
- each scored item carries a compact `a1` and `d1` parameter pair
- PHQ-oriented outputs are already stored as a theta lookup table
- threshold application logic is deterministic and lightweight

That means the browser does not need to run `mirt` itself. Instead, it can:

1. load a JSON export of the production-safe engines
2. estimate `Theta_Harmonized` on a fixed theta grid using EAP scoring
3. interpolate `PHQ_Expected_Sum` and `PHQ_Prob_GE10`
4. apply the stored threshold rule

## Scope For The First Browser Version

The first static build should support:

- production-safe **1-factor** engines only
- manual single-record scoring
- browser-local batch CSV scoring
- downloadable scored CSV
- no server round trips for scoring

The browser build should **exclude**:

- audit-only 2-factor and 4-factor exports
- server-side persistence
- direct write-back of uploaded raw data

## File Structure

```text
portable_irt_engine/
  BROWSER_JS_ARCHITECTURE.md
  data-raw/
    export_browser_engine_bundle.R

website/
  _tools/
    irt-scoring-static/
      README.md
      index.html
      styles.css
      app.js
      scorer.js
      engines/
        production_engines.json
        example_batch.csv
```

## Runtime Flow

### 1. Load The Engine Bundle

The app loads `engines/production_engines.json`, which contains:

- engine metadata
- item parameter pairs
- cutoff tables
- PHQ theta maps
- display labels and notes

### 2. Standardize Optional Group Fields

Before thresholding:

- `SEX` is standardized to `Male` / `Female`
- `AGEGRP` is standardized to `Age_17_19` / `Age_20_24`
- `SEX_AGEGRP` is derived when possible

### 3. Estimate Theta

For each response pattern:

- compute the response likelihood across a fixed theta grid
- combine with a standard normal prior
- normalize posterior weights
- calculate the EAP theta estimate

The browser version uses the exported `Theta` grid from the PHQ map to avoid a
second grid definition.

### 4. Map Theta To PHQ-Oriented Outputs

Use linear interpolation over the stored theta map to derive:

- `PHQ_Expected_Sum`
- `PHQ_Prob_GE10`

### 5. Apply The Stored Threshold Rule

Apply one of:

- overall threshold
- sex-specific threshold
- age-group-specific threshold
- sex-by-age-group threshold

If subgroup fields are missing or unmatched, fall back to the stored overall
threshold.

### 6. Return Reader-Friendly Outputs

Per scored row:

- `SSQ10_Total_Score`
- `Theta_Harmonized`
- `Depression_Binary`
- `PHQ_Expected_Sum`
- `PHQ_Prob_GE10`
- `Cutoff_Theta_Applied`
- `Cutoff_Grouping_Used`
- `Cutoff_Group_Level_Used`
- `Cutoff_Apply_Mode`
- `Cutoff_Method`
- `Engine_ID`
- `Engine_Label`
- `Items_Answered`

## Validation Plan

The browser build should be validated against the R engine with a fixed test set:

1. export a known example batch
2. score it in R
3. score it in the browser build
4. compare row-by-row outputs

Acceptance target for the 1-factor engines:

- exact match for threshold outputs and metadata
- near-identical theta estimates within a very small tolerance
- near-identical interpolated PHQ outputs

## Current Prototype

A first static prototype is placed at:

`website/_tools/irt-scoring-static/index.html`

It is intended as a proof-of-feasibility for GitHub Pages hosting rather than a
fully validated production release.
