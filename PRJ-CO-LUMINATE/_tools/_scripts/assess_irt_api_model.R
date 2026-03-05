#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

project_root <- "."
objects_dir <- file.path(project_root, "_tools", "_objects")
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

resolve_engine_bundle <- function(root) {
  env_path <- Sys.getenv("COLUMINATE_ENGINES_PATH", unset = "")
  candidates <- c(
    env_path,
    file.path(root, "_tools", "_objects", "irt_joint_models.rds"),
    file.path(
      root, "..", "OBJ00-Datasets Preperation", "data_management",
      "data_examples", "1_staging_snippets", "derived_models", "irt_joint_models.rds"
    )
  )
  candidates <- unique(candidates[nzchar(candidates)])
  idx <- which(file.exists(candidates))[1]
  if (is.na(idx)) {
    stop(
      "Unable to locate irt_joint_models.rds.\n",
      "Set COLUMINATE_ENGINES_PATH or place file in _tools/_objects."
    )
  }
  candidates[idx]
}

bundle_path <- resolve_engine_bundle(project_root)
engines <- readRDS(bundle_path)

if (!is.list(engines) || length(engines) == 0) {
  stop("Engine bundle is empty or invalid: ", bundle_path)
}

summarise_engine <- function(engine_id, engine_obj) {
  ct <- engine_obj$cutoff_table %||% data.frame()
  overall <- if (nrow(ct) > 0) {
    ct[ct$Grouping == "none" & ct$Group_Level == "Overall", , drop = FALSE]
  } else {
    data.frame()
  }

  data.frame(
    engine_id = engine_id,
    n_items = length(engine_obj$items %||% character()),
    n_factors = as.integer(engine_obj$n_factors %||% NA_integer_),
    model_type = as.character(engine_obj$model_type %||% NA_character_),
    cutoff_apply_mode = as.character(engine_obj$cutoff_apply_mode %||% NA_character_),
    cutoff_groupings = if (nrow(ct) > 0) {
      paste(sort(unique(as.character(ct$Grouping))), collapse = "|")
    } else {
      NA_character_
    },
    overall_cutoff_theta = if (nrow(overall) > 0) as.numeric(overall$Cutoff_Theta[1]) else NA_real_,
    overall_accuracy = if (nrow(overall) > 0) as.numeric(overall$Accuracy[1]) else NA_real_,
    overall_sensitivity = if (nrow(overall) > 0) as.numeric(overall$Sensitivity[1]) else NA_real_,
    overall_specificity = if (nrow(overall) > 0) as.numeric(overall$Specificity[1]) else NA_real_,
    overall_ppv = if (nrow(overall) > 0) as.numeric(overall$PPV[1]) else NA_real_,
    overall_npv = if (nrow(overall) > 0) as.numeric(overall$NPV[1]) else NA_real_,
    overall_kappa = if (nrow(overall) > 0) as.numeric(overall$Kappa[1]) else NA_real_,
    has_theta_map = is.data.frame(engine_obj$phq_theta_map) && nrow(engine_obj$phq_theta_map) > 0,
    stringsAsFactors = FALSE
  )
}

engine_table <- do.call(
  rbind,
  Map(summarise_engine, names(engines), engines)
)

all_items <- unique(unlist(lapply(engines, function(x) x$items)))

assessment <- list(
  generated_at = as.character(Sys.time()),
  engine_bundle_path = normalizePath(bundle_path, mustWork = FALSE),
  n_engines = length(engines),
  input_items = sort(all_items),
  optional_covariates = c("SEX", "AGE", "AGEGRP", "SEX_AGEGRP"),
  outputs = c("Theta_Harmonized", "Depression_Binary", "PHQ_Expected_Sum", "PHQ_Prob_GE10"),
  engines = engine_table
)

saveRDS(assessment, file.path(objects_dir, "irt_api_assessment.rds"))
jsonlite::write_json(
  assessment,
  path = file.path(objects_dir, "irt_api_assessment.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

cat("Created IRT API assessment objects in:", objects_dir, "\n")
cat("Engine bundle source:", normalizePath(bundle_path, mustWork = FALSE), "\n")
cat("Engines:", length(engines), "| Unique SSQ items:", length(all_items), "\n")
