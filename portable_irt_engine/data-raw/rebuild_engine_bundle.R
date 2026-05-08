`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

engine_diagnostics <- function(engine) {
  ssq_a <- subset(
    engine$parameters,
    grepl("^SSQ", item) & grepl("^a", name)
  )

  total_a <- nrow(ssq_a)
  nonzero_a <- sum(abs(ssq_a$value) > 1e-10, na.rm = TRUE)

  map_var_sum <- if (!is.null(engine$phq_theta_map)) {
    stats::var(engine$phq_theta_map$PHQ_Expected_Sum, na.rm = TRUE)
  } else {
    NA_real_
  }

  map_var_prob <- if (!is.null(engine$phq_theta_map)) {
    stats::var(engine$phq_theta_map$PHQ_Prob_GE10, na.rm = TRUE)
  } else {
    NA_real_
  }

  scoring_ready <- isTRUE(
    total_a > 0 &&
      nonzero_a > 0 &&
      is.finite(map_var_sum) &&
      map_var_sum > 0 &&
      is.finite(map_var_prob) &&
      map_var_prob > 0
  )

  list(
    total_a = total_a,
    nonzero_a = nonzero_a,
    map_var_sum = map_var_sum,
    map_var_prob = map_var_prob,
    scoring_ready = scoring_ready
  )
}

sanitize_meta <- function(meta) {
  meta$source_file <- if (!is.null(meta$source_file) && length(meta$source_file) == 1) {
    basename(meta$source_file)
  } else {
    meta$source_file
  }

  meta$source_scope <- "portable_irt_engine_bundle"
  meta$source_note <- paste(
    "Metadata were rebuilt for the portable engine bundle.",
    "Absolute local source paths were removed from the packaged artifact."
  )

  meta
}

load_source_artifact <- function(spec, winners_dir, batch_dir) {
  if (identical(spec$source_kind, "winner")) {
    engine <- readRDS(file.path(winners_dir, spec$source_engine))
    meta <- readRDS(file.path(winners_dir, spec$source_meta))
    meta$source_kind <- "winner"
    return(list(engine = engine, meta = sanitize_meta(meta)))
  }

  if (identical(spec$source_kind, "batch")) {
    batch_file <- file.path(batch_dir, paste0("IRT_batch__", spec$run_id, ".rds"))
    res <- readRDS(batch_file)
    meta <- list(
      run_id = spec$run_id,
      source_file = basename(batch_file),
      source_kind = "batch",
      input_structure = res$input_structure %||% NULL,
      threshold_applied = res$threshold %||% NULL,
      metrics_applied = res$metrics %||% NULL,
      batch_exported_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    return(list(engine = res$scoring_engine, meta = sanitize_meta(meta)))
  }

  stop("Unknown source_kind: ", spec$source_kind, call. = FALSE)
}

build_engine_specs <- function() {
  data.frame(
    engine_id = c(
      "ssq10_sens_at_spec_1f_overall",
      "ssq10_sens_at_spec_1f_sex",
      "ssq10_sens_at_spec_1f_agegrp",
      "ssq10_sens_at_spec_1f_sex_agegrp",
      "ssq10_spec_at_sens_1f_overall",
      "ssq10_spec_at_sens_1f_sex",
      "ssq10_spec_at_sens_1f_agegrp",
      "ssq10_spec_at_sens_1f_sex_agegrp",
      "ssq10_youden_1f_overall",
      "ssq10_youden_1f_sex",
      "ssq10_youden_1f_agegrp",
      "ssq10_youden_1f_sex_agegrp",
      "ssq10_tcc_1f_overall",
      "ssq10_sens_at_spec_2f_overall",
      "ssq10_sens_at_spec_2f_sex",
      "ssq10_sens_at_spec_2f_agegrp",
      "ssq10_sens_at_spec_2f_sex_agegrp",
      "ssq10_youden_4f_overall"
    ),
    label = c(
      "SSQ-10 Sens@Spec 1-Factor Overall",
      "SSQ-10 Sens@Spec 1-Factor by Sex",
      "SSQ-10 Sens@Spec 1-Factor by Age Group",
      "SSQ-10 Sens@Spec 1-Factor by Sex and Age Group",
      "SSQ-10 Spec@Sens 1-Factor Overall",
      "SSQ-10 Spec@Sens 1-Factor by Sex",
      "SSQ-10 Spec@Sens 1-Factor by Age Group",
      "SSQ-10 Spec@Sens 1-Factor by Sex and Age Group",
      "SSQ-10 Youden 1-Factor Overall",
      "SSQ-10 Youden 1-Factor by Sex",
      "SSQ-10 Youden 1-Factor by Age Group",
      "SSQ-10 Youden 1-Factor by Sex and Age Group",
      "SSQ-10 Model-TCC 1-Factor Overall",
      "SSQ-10 Sens@Spec 2-Factor Overall",
      "SSQ-10 Sens@Spec 2-Factor by Sex",
      "SSQ-10 Sens@Spec 2-Factor by Age Group",
      "SSQ-10 Sens@Spec 2-Factor by Sex and Age Group",
      "SSQ-10 Youden 4-Factor Overall"
    ),
    source_kind = c(
      "winner", "winner", "winner", "winner",
      "batch", "batch", "batch", "batch",
      "batch", "batch", "batch", "batch",
      "winner",
      "winner", "winner", "winner", "winner", "winner"
    ),
    source_engine = c(
      "SSQ10_Harmonized_Scoring_Engine__sens_at_spec_1_factor_model_overall.rds",
      "SSQ10_Harmonized_Scoring_Engine__sens_at_spec_1_factor_model_gender.rds",
      "SSQ10_Harmonized_Scoring_Engine__sens_at_spec_1_factor_model_age_group.rds",
      "SSQ10_Harmonized_Scoring_Engine__sens_at_spec_1_factor_model_gender_age_group.rds",
      NA, NA, NA, NA,
      NA, NA, NA, NA,
      "SSQ10_Harmonized_Scoring_Engine__tcc_1_factor_model_overall.rds",
      "SSQ10_Harmonized_Scoring_Engine__sens_at_spec_2_factor_model_overall.rds",
      "SSQ10_Harmonized_Scoring_Engine__sens_at_spec_2_factor_model_gender.rds",
      "SSQ10_Harmonized_Scoring_Engine__sens_at_spec_2_factor_model_age_group.rds",
      "SSQ10_Harmonized_Scoring_Engine__sens_at_spec_2_factor_model_gender_age_group.rds",
      "SSQ10_Harmonized_Scoring_Engine__youden_4_factor_model_overall.rds"
    ),
    source_meta = c(
      "SELECTED_SCENARIO_META__sens_at_spec_1_factor_model_overall.rds",
      "SELECTED_SCENARIO_META__sens_at_spec_1_factor_model_gender.rds",
      "SELECTED_SCENARIO_META__sens_at_spec_1_factor_model_age_group.rds",
      "SELECTED_SCENARIO_META__sens_at_spec_1_factor_model_gender_age_group.rds",
      NA, NA, NA, NA,
      NA, NA, NA, NA,
      "SELECTED_SCENARIO_META__tcc_1_factor_model_overall.rds",
      "SELECTED_SCENARIO_META__sens_at_spec_2_factor_model_overall.rds",
      "SELECTED_SCENARIO_META__sens_at_spec_2_factor_model_gender.rds",
      "SELECTED_SCENARIO_META__sens_at_spec_2_factor_model_age_group.rds",
      "SELECTED_SCENARIO_META__sens_at_spec_2_factor_model_gender_age_group.rds",
      "SELECTED_SCENARIO_META__youden_4_factor_model_overall.rds"
    ),
    run_id = c(
      NA, NA, NA, NA,
      "23737ef9c596f9ee",
      "baee49342bab13d7",
      "78cac47ab6a33b18",
      "88d1f0c1206bf332",
      "f40efd1d11d6399d",
      "b7b7eacc43b2612b",
      "74a8da946f18a514",
      "2adf849547c2e8ae",
      NA,
      NA, NA, NA, NA, NA
    ),
    default = c(
      TRUE, rep(FALSE, 17)
    ),
    available_in_app = c(
      rep(TRUE, 13),
      rep(FALSE, 5)
    ),
    description = c(
      "Working 1-factor overall SSQ-10 harmonisation engine using Sens@Spec cutoffs.",
      "Working 1-factor SSQ-10 harmonisation engine using Sens@Spec cutoffs applied by sex.",
      "Working 1-factor SSQ-10 harmonisation engine using Sens@Spec cutoffs applied by age group.",
      "Working 1-factor SSQ-10 harmonisation engine using Sens@Spec cutoffs applied by sex and age group.",
      "Working 1-factor overall SSQ-10 harmonisation engine using Spec@Sens cutoffs.",
      "Working 1-factor SSQ-10 harmonisation engine using Spec@Sens cutoffs applied by sex.",
      "Working 1-factor SSQ-10 harmonisation engine using Spec@Sens cutoffs applied by age group.",
      "Working 1-factor SSQ-10 harmonisation engine using Spec@Sens cutoffs applied by sex and age group.",
      "Working 1-factor overall SSQ-10 harmonisation engine using Youden cutoffs.",
      "Working 1-factor SSQ-10 harmonisation engine using Youden cutoffs applied by sex.",
      "Working 1-factor SSQ-10 harmonisation engine using Youden cutoffs applied by age group.",
      "Working 1-factor SSQ-10 harmonisation engine using Youden cutoffs applied by sex and age group.",
      "Working 1-factor overall SSQ-10 harmonisation engine using the model-TCC cutoff.",
      "Audit-only 2-factor overall export retained because the saved slopes are flat.",
      "Audit-only 2-factor sex export retained because the saved slopes are flat.",
      "Audit-only 2-factor age-group export retained because the saved slopes are flat.",
      "Audit-only 2-factor sex-and-age-group export retained because the saved slopes are flat.",
      "Audit-only 4-factor overall export retained because the saved slopes are flat."
    ),
    notes = c(
      "Default production engine for the app.",
      "Group-aware engine using sex-specific thresholds.",
      "Group-aware engine using age-group-specific thresholds.",
      "Group-aware engine using sex-by-age-group-specific thresholds.",
      "Alternative 1-factor cutoff family sourced from the batch scenario registry.",
      "Alternative 1-factor cutoff family sourced from the batch scenario registry.",
      "Alternative 1-factor cutoff family sourced from the batch scenario registry.",
      "Alternative 1-factor cutoff family sourced from the batch scenario registry.",
      "Alternative 1-factor cutoff family sourced from the batch scenario registry.",
      "Alternative 1-factor cutoff family sourced from the batch scenario registry.",
      "Alternative 1-factor cutoff family sourced from the batch scenario registry.",
      "Alternative 1-factor cutoff family sourced from the batch scenario registry.",
      "Legacy/overall TCC engine retained for comparison with empirical cutoffs.",
      "Keep hidden from live scoring until the multidimensional export is repaired upstream.",
      "Keep hidden from live scoring until the multidimensional export is repaired upstream.",
      "Keep hidden from live scoring until the multidimensional export is repaired upstream.",
      "Keep hidden from live scoring until the multidimensional export is repaired upstream.",
      "Keep hidden from live scoring until the multidimensional export is repaired upstream."
    ),
    stringsAsFactors = FALSE
  )
}

rebuild_engine_bundle <- function(
    package_root = ".",
    winners_dir = "../statistical_analysis/output/objects/irt/winners",
    batch_dir = "../statistical_analysis/output/objects/irt/batch",
    engine_dir = file.path("inst", "engines")) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  package_root <- normalizePath(package_root, mustWork = TRUE)
  setwd(package_root)

  winners_dir <- normalizePath(winners_dir, mustWork = TRUE)
  batch_dir <- normalizePath(batch_dir, mustWork = TRUE)
  engine_dir <- file.path(package_root, engine_dir)
  dir.create(engine_dir, recursive = TRUE, showWarnings = FALSE)

  specs <- build_engine_specs()
  catalog_rows <- vector("list", nrow(specs))

  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, , drop = FALSE]
    artifact <- load_source_artifact(spec, winners_dir = winners_dir, batch_dir = batch_dir)
    engine <- artifact$engine
    meta <- artifact$meta
    diag <- engine_diagnostics(engine)

    if (isTRUE(spec$available_in_app[[1]]) && !isTRUE(diag$scoring_ready)) {
      stop(
        "Engine marked available_in_app failed diagnostics: ",
        spec$engine_id[[1]],
        call. = FALSE
      )
    }

    engine_file_name <- paste0(spec$engine_id[[1]], ".rds")
    meta_file_name <- paste0(spec$engine_id[[1]], "_meta.rds")

    saveRDS(engine, file.path(engine_dir, engine_file_name))
    saveRDS(meta, file.path(engine_dir, meta_file_name))

    cutoff_method <- NA_character_
    if (!is.null(engine$cutoff_table) && nrow(engine$cutoff_table) > 0) {
      cutoff_method <- as.character(engine$cutoff_table$Cutoff_Method[[1]])
    }

    catalog_rows[[i]] <- data.frame(
      engine_id = spec$engine_id[[1]],
      label = spec$label[[1]],
      file_name = engine_file_name,
      meta_file_name = meta_file_name,
      default = spec$default[[1]],
      available_in_app = spec$available_in_app[[1]],
      scoring_ready = diag$scoring_ready,
      n_factors = engine$n_factors %||% NA_integer_,
      cutoff_apply_mode = engine$cutoff_apply_mode %||% NA_character_,
      cutoff_method = cutoff_method,
      source_kind = spec$source_kind[[1]],
      source_ref = if (identical(spec$source_kind[[1]], "batch")) spec$run_id[[1]] else spec$source_engine[[1]],
      description = spec$description[[1]],
      notes = spec$notes[[1]],
      total_a = diag$total_a,
      nonzero_a = diag$nonzero_a,
      map_var_sum = diag$map_var_sum,
      map_var_prob = diag$map_var_prob,
      stringsAsFactors = FALSE
    )
  }

  catalog <- do.call(rbind, catalog_rows)
  utils::write.csv(catalog, file.path(engine_dir, "engine_catalog.csv"), row.names = FALSE)

  invisible(catalog)
}

if (sys.nframe() == 0) {
  rebuild_engine_bundle()
}
