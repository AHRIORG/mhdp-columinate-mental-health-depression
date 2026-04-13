# ==============================================================================
# UTILITY: Reporting Standards & Sensitivity Diagnostics
# PURPOSE: Build LCGA sensitivity diagnostics and GRoLTS/SMART-oriented
#          reporting checklists for the public trajectory bundle.
# ==============================================================================

library(dplyr)
library(glue)
library(gt)
library(purrr)
library(rlang)
library(stringr)
library(tibble)
library(tidyr)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

format_pct_string <- function(x, decimals = 1) {
  out <- sprintf(paste0("%.", decimals, "f%%"), x * 100)
  out[is.na(x)] <- "NA"
  out
}

format_num_string <- function(x, decimals = 3) {
  out <- sprintf(paste0("%.", decimals, "f"), x)
  out[is.na(x)] <- "NA"
  out
}

style_status_column <- function(gt_tbl, status_col = "Status") {
  pass_expr <- rlang::parse_expr(glue("{status_col} %in% c('Pass', 'Reported')"))
  watch_expr <- rlang::parse_expr(glue("{status_col} %in% c('Watch', 'Partial', 'Not applicable', 'Info', 'Unavailable')"))
  fail_expr <- rlang::parse_expr(glue("{status_col} %in% c('Fail', 'Not reported')"))

  gt_tbl %>%
    tab_style(
      style = cell_text(weight = "bold", color = "#0B6E4F"),
      locations = cells_body(columns = all_of(status_col), rows = !!pass_expr)
    ) %>%
    tab_style(
      style = cell_text(weight = "bold", color = "#A15C00"),
      locations = cells_body(columns = all_of(status_col), rows = !!watch_expr)
    ) %>%
    tab_style(
      style = cell_text(weight = "bold", color = "#B22222"),
      locations = cells_body(columns = all_of(status_col), rows = !!fail_expr)
    )
}

compute_lcga_class_weights <- function(model_obj) {
  if (is.null(model_obj$savedata)) {
    return(tibble())
  }

  savedata <- tibble::as_tibble(model_obj$savedata)
  class_col <- grep("^C$|^MLCC$", names(savedata), value = TRUE)[1]
  if (is.na(class_col)) {
    return(tibble())
  }

  savedata %>%
    transmute(LatentClass = as.integer(.data[[class_col]])) %>%
    filter(!is.na(LatentClass)) %>%
    count(LatentClass, name = "ClassCount") %>%
    mutate(ClassProportion = ClassCount / sum(ClassCount))
}

resolve_expected_levels <- function(meta) {
  levels <- meta$transform_levels %||% meta$levels
  if (is.null(levels)) character(0) else as.character(levels)
}

map_value_to_category <- function(value, expected_levels) {
  if (length(expected_levels) == 0) {
    return(as.character(value))
  }
  idx <- suppressWarnings(as.integer(value)) + 1L
  ifelse(!is.na(idx) & idx >= 1 & idx <= length(expected_levels), expected_levels[idx], as.character(value))
}

collect_empirical_category_probs <- function(df_mplus, item_names_list, outcome_meta, ages, missing_code = -9999) {
  if (is.null(df_mplus)) {
    return(tibble())
  }

  map_dfr(outcome_meta, function(meta) {
    if (!tolower(meta$type %||% "") %in% c("binary", "ordinal")) {
      return(tibble())
    }

    items <- item_names_list[[meta$wide_prefix]]
    if (is.null(items) || length(items) == 0) {
      return(tibble())
    }

    expected_levels <- resolve_expected_levels(meta)
    long_df <- df_mplus %>%
      select(all_of(items)) %>%
      pivot_longer(cols = everything(), names_to = "item", values_to = "value") %>%
      mutate(Age = as.integer(str_extract(item, "\\d+$"))) %>%
      filter(!is.na(value), value != missing_code)

    if (nrow(long_df) == 0) {
      return(tidyr::expand_grid(
        prefix = meta$wide_prefix,
        variable = meta$label,
        Age = ages,
        Category = factor(expected_levels, levels = expected_levels)
      ) %>%
        mutate(Observed_N = 0, Observed_Prob = 0))
    }

    long_df %>%
      mutate(
        Category = map_value_to_category(value, expected_levels),
        Category = factor(Category, levels = expected_levels)
      ) %>%
      count(Age, Category, name = "Observed_N") %>%
      complete(Age = ages, Category = factor(expected_levels, levels = expected_levels), fill = list(Observed_N = 0)) %>%
      group_by(Age) %>%
      mutate(Observed_Prob = ifelse(sum(Observed_N) > 0, Observed_N / sum(Observed_N), 0)) %>%
      ungroup() %>%
      mutate(
        prefix = meta$wide_prefix,
        variable = meta$label
      ) %>%
      select(prefix, variable, Age, Category, Observed_N, Observed_Prob)
  })
}

compute_model_implied_category_probs <- function(model_obj, item_names_list, outcome_meta, ages, time_scores, is_joint) {
  if (is.null(model_obj$parameters$unstandardized)) {
    return(tibble())
  }

  class_weights <- compute_lcga_class_weights(model_obj)
  if (nrow(class_weights) == 0) {
    return(tibble())
  }

  unstd <- model_obj$parameters$unstandardized

  map_dfr(seq_along(outcome_meta), function(j) {
    meta <- outcome_meta[[j]]
    if (!tolower(meta$type %||% "") %in% c("binary", "ordinal")) {
      return(tibble())
    }

    expected_levels <- resolve_expected_levels(meta)
    items <- item_names_list[[meta$wide_prefix]]
    if (is.null(items) || length(items) == 0) {
      return(tibble())
    }

    gf_means <- if (is_joint) extract_joint_gf_means(unstd, j = j) else extract_univariate_gf_means(unstd)
    if (nrow(gf_means) == 0) {
      return(tibble())
    }

    thresholds <- extract_all_thresholds(unstd, items, ages)
    if (is.null(thresholds) || nrow(thresholds) == 0) {
      return(tibble())
    }

    model_probs <- compute_probs_all_cats(gf_means, ages, time_scores, thresholds, expected_levels) %>%
      mutate(
        LatentClass = suppressWarnings(as.integer(as.character(LatentClass))),
        Category = factor(as.character(Category), levels = expected_levels)
      ) %>%
      left_join(class_weights, by = "LatentClass") %>%
      group_by(Age, Category) %>%
      summarise(Model_Prob = sum(Prob * ClassProportion, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        prefix = meta$wide_prefix,
        variable = meta$label
      ) %>%
      select(prefix, variable, Age, Category, Model_Prob)

    model_probs
  })
}

summarise_lcga_sensitivity <- function(empirical_probs, model_probs, outcome_meta) {
  compare_df <- full_join(
    empirical_probs,
    model_probs,
    by = c("prefix", "variable", "Age", "Category")
  ) %>%
    mutate(
      AbsError = abs(Observed_Prob - Model_Prob),
      SqError = (Observed_Prob - Model_Prob)^2
    )

  summary_df <- map_dfr(outcome_meta, function(meta) {
    if (!tolower(meta$type %||% "") %in% c("binary", "ordinal")) {
      return(tibble())
    }

    expected_levels <- resolve_expected_levels(meta)
    compare_sub <- compare_df %>% filter(prefix == meta$wide_prefix)
    observed_present <- empirical_probs %>%
      filter(prefix == meta$wide_prefix, Observed_N > 0) %>%
      distinct(Category) %>%
      pull(Category) %>%
      as.character()

    eval_categories <- if (tolower(meta$type %||% "") == "binary" && length(expected_levels) > 0) {
      tail(expected_levels, 1)
    } else {
      expected_levels
    }

    compare_eval <- compare_sub %>%
      filter(as.character(Category) %in% eval_categories) %>%
      filter(!is.na(Observed_Prob), !is.na(Model_Prob))

    missing_categories <- setdiff(expected_levels, observed_present)
    mae <- if (nrow(compare_eval) > 0) mean(compare_eval$AbsError, na.rm = TRUE) else NA_real_
    rmse <- if (nrow(compare_eval) > 0) sqrt(mean(compare_eval$SqError, na.rm = TRUE)) else NA_real_
    max_abs <- if (nrow(compare_eval) > 0) max(compare_eval$AbsError, na.rm = TRUE) else NA_real_

    status <- case_when(
      length(missing_categories) > 0 ~ "Fail",
      is.na(mae) ~ "Unavailable",
      mae <= 0.05 ~ "Pass",
      mae <= 0.10 ~ "Watch",
      TRUE ~ "Fail"
    )

    tibble(
      prefix = meta$wide_prefix,
      variable = meta$label,
      evaluated_categories = paste(eval_categories, collapse = ", "),
      missing_categories = if (length(missing_categories) == 0) "None" else paste(missing_categories, collapse = ", "),
      NComparisons = nrow(compare_eval),
      MAE = mae,
      RMSE = rmse,
      MaxAbs = max_abs,
      Status = status
    )
  })

  list(summary = summary_df, detail = compare_df)
}

generate_lcga_sensitivity_table <- function(
    fit_df,
    best_k,
    class_assign_data,
    sensitivity_summary,
    analysis_label = "Latent Class Growth Analysis"
) {
  avepp_summary <- summarise_lcga_avepp(class_assign_data)

  selected_row <- fit_df %>% filter(Classes == best_k) %>% slice(1)
  local_rows <- fit_df %>%
    filter(Classes %in% c(best_k - 1, best_k, best_k + 1)) %>%
    arrange(Classes)

  local_evidence <- if (nrow(local_rows) == 0) {
    "Local K comparison unavailable."
  } else {
    paste(
      apply(local_rows, 1, function(row) {
        glue(
          "K={row[['Classes']]} | BIC={format_num_string(as.numeric(row[['BIC']]), 1)} | Entropy={format_num_string(as.numeric(row[['Entropy']]), 3)} | MinClass={format_pct_string(as.numeric(row[['MinClassPct']]), 1)}"
        )
      }),
      collapse = " || "
    )
  }

  min_class_pct <- if (nrow(selected_row) > 0 && "MinClassPct" %in% names(selected_row)) as.numeric(selected_row$MinClassPct[[1]]) else NA_real_
  avepp_status <- case_when(
    is.na(avepp_summary$min_diag) ~ "Unavailable",
    avepp_summary$min_diag >= 0.80 ~ "Pass",
    avepp_summary$min_diag >= 0.70 ~ "Watch",
    TRUE ~ "Fail"
  )

  model_rows <- tibble(
    Domain = "Model-level",
    Diagnostic = c("Local K comparison", "Selected-model class size", "Selected-model classification precision"),
    Evidence = c(
      local_evidence,
      glue("Min class proportion at K={best_k}: {format_pct_string(min_class_pct, 1)}"),
      glue("Mean diagonal AvePP = {format_num_string(avepp_summary$mean_diag, 3)}; Min diagonal AvePP = {format_num_string(avepp_summary$min_diag, 3)}")
    ),
    Status = c(
      "Info",
      if (is.na(min_class_pct)) "Unavailable" else if (min_class_pct >= 0.05) "Pass" else "Fail",
      avepp_status
    )
  )

  outcome_rows <- sensitivity_summary %>%
    transmute(
      Domain = "Outcome-level",
      Diagnostic = variable,
      Evidence = paste0(
        "Checked: ", evaluated_categories,
        "; Missing observed categories: ", missing_categories,
        "; MAE = ", format_num_string(MAE, 3),
        "; Max abs error = ", format_num_string(MaxAbs, 3)
      ),
      Status = Status
    )

  bind_rows(model_rows, outcome_rows) %>%
    gt(groupname_col = "Domain") %>%
    tab_header(
      title = md(glue("**Table 7: Misspecification & Sensitivity Diagnostics - {analysis_label}**")),
      subtitle = "Observed-versus-model-implied marginal checks plus local class-enumeration context"
    ) %>%
    cols_label(
      Diagnostic = "Diagnostic",
      Evidence = "Evidence",
      Status = "Status"
    ) %>%
    tab_source_note("Note: Outcome-level diagnostics compare observed marginal category probabilities against model-implied marginals at the fitted ages. MAE thresholds are operational heuristics: <=0.05 pass, 0.05-0.10 watch, >0.10 fail.") %>%
    style_status_column("Status")
}

build_lcga_final_solution_summary <- function(
    best_model,
    class_assign_data,
    outcome_meta,
    analysis_label,
    is_joint = FALSE
) {
  if (is.null(best_model$parameters$unstandardized)) {
    return(tibble())
  }

  unstd <- best_model$parameters$unstandardized
  class_weights <- class_assign_data %>%
    count(Assigned_Class, Assigned_Class_Label, name = "ClassN") %>%
    mutate(ClassProportion = ClassN / sum(ClassN))

  if (is_joint) {
    est_df <- purrr::imap_dfr(outcome_meta, function(meta, idx) {
      unstd %>%
        filter(paramHeader == "Means", !is.na(LatentClass)) %>%
        mutate(
          paramU = toupper(trimws(param)),
          outcome_index = suppressWarnings(as.integer(stringr::str_extract(paramU, "\\d+$"))),
          parameter_letter = substr(paramU, 1, 1)
        ) %>%
        filter(outcome_index == idx, parameter_letter %in% c("I", "S", "Q")) %>%
        transmute(
          Outcome = meta$label %||% analysis_label,
          Assigned_Class = as.integer(LatentClass),
          Parameter = parameter_letter,
          Estimate = as.numeric(est),
          SE = as.numeric(se)
        )
    })
  } else {
    est_df <- unstd %>%
      filter(paramHeader == "Means", !is.na(LatentClass)) %>%
      mutate(
        paramU = toupper(trimws(param)),
        parameter_letter = substr(paramU, 1, 1)
      ) %>%
      filter(parameter_letter %in% c("I", "S", "Q")) %>%
      group_by(LatentClass, parameter_letter) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(
        Outcome = outcome_meta[[1]]$label %||% analysis_label,
        Assigned_Class = as.integer(LatentClass),
        Parameter = parameter_letter,
        Estimate = as.numeric(est),
        SE = as.numeric(se)
      )
  }

  if (nrow(est_df) == 0) {
    return(tibble())
  }

  est_df %>%
    left_join(class_weights, by = "Assigned_Class") %>%
    mutate(
      Parameter = recode(
        Parameter,
        I = "Intercept (I)",
        S = "Linear Slope (S)",
        Q = "Quadratic Slope (Q)",
        .default = Parameter
      ),
      Lower95CI = Estimate - 1.96 * SE,
      Upper95CI = Estimate + 1.96 * SE
    ) %>%
    select(
      Outcome,
      Assigned_Class,
      Assigned_Class_Label,
      ClassN,
      ClassProportion,
      Parameter,
      Estimate,
      SE,
      Lower95CI,
      Upper95CI
    ) %>%
    arrange(Outcome, Assigned_Class, Parameter)
}

generate_lcga_final_solution_summary_table <- function(summary_df, analysis_label) {
  if (nrow(summary_df) == 0) {
    return(NULL)
  }

  summary_df %>%
    mutate(
      Assigned_Class_Label = as.character(Assigned_Class_Label)
    ) %>%
    gt(groupname_col = "Outcome") %>%
    tab_header(
      title = md(glue("**Table 6A: Final Solution Numerical Summary - {analysis_label}**")),
      subtitle = "Selected-model class sizes and growth-factor estimates with uncertainty intervals"
    ) %>%
    cols_label(
      Assigned_Class_Label = "Assigned Class",
      ClassN = "N",
      ClassProportion = "Class Proportion",
      Parameter = "Parameter",
      Estimate = "Estimate",
      SE = "SE",
      Lower95CI = "Lower 95% CI",
      Upper95CI = "Upper 95% CI"
    ) %>%
    cols_hide(columns = "Assigned_Class") %>%
    fmt_percent(columns = "ClassProportion", decimals = 1) %>%
    fmt_number(columns = c("Estimate", "SE", "Lower95CI", "Upper95CI"), decimals = 3) %>%
    tab_source_note("Note: Values are unstandardized selected-model growth-factor estimates with Wald-style 95% confidence intervals and selected-model class sizes.")
}

build_missingness_correlates_summary <- function(
    df_mplus,
    item_names_list,
    outcome_meta,
    id_var = "USUBJID",
    missing_code = -9999,
    dt_multi = NULL,
    candidate_vars = NULL
) {
  modeled_items <- unique(unlist(item_names_list, use.names = FALSE))
  if (is.null(df_mplus) || length(modeled_items) == 0 || !id_var %in% names(df_mplus)) {
    return(tibble())
  }

  subject_missing <- df_mplus %>%
    select(any_of(c(id_var, modeled_items))) %>%
    mutate(
      TotalModeled = length(modeled_items),
      ObservedCount = rowSums(across(all_of(modeled_items), ~ !is.na(.x) & .x != missing_code)),
      ObservedPct = ObservedCount / TotalModeled,
      AnyMissing = ObservedCount < TotalModeled
    ) %>%
    select(any_of(id_var), TotalModeled, ObservedCount, ObservedPct, AnyMissing)

  earliest_item_map <- purrr::imap(item_names_list, function(items, prefix) {
    item_ages <- suppressWarnings(as.integer(stringr::str_extract(items, "\\d+$")))
    items[which.min(item_ages)]
  }) %>%
    unlist(use.names = TRUE)

  modeled_candidate_df <- df_mplus %>%
    select(any_of(c(id_var, unname(earliest_item_map))))

  modeled_lookup <- purrr::imap_dfr(outcome_meta, function(meta, idx) {
    source_var <- earliest_item_map[[meta$wide_prefix]]
    lvls <- meta$transform_levels %||% meta$levels
    if (is.null(source_var) || is.null(lvls) || length(lvls) == 0) {
      return(tibble())
    }
    tibble(
      Variable = source_var,
      RawValue = as.character(seq_along(lvls) - 1),
      Level = lvls
    )
  })

  modeled_meta <- tibble(
    Variable = unname(earliest_item_map),
    VariableLabel = purrr::map_chr(outcome_meta, ~ .x$label %||% .x$wide_prefix),
    Source = "Modeled indicator (earliest wave)"
  )

  baseline_candidate_df <- NULL
  baseline_meta <- tibble()
  if (!is.null(dt_multi) && nrow(dt_multi) > 0) {
    present_candidates <- intersect(candidate_vars %||% character(), names(dt_multi))
    if (length(present_candidates) > 0 && id_var %in% names(dt_multi)) {
      baseline_candidate_df <- dt_multi %>%
        arrange(.data[[id_var]]) %>%
        group_by(.data[[id_var]]) %>%
        summarise(
          across(
            all_of(present_candidates),
            ~ {
              vals <- .x[!is.na(.x)]
              if (length(vals) == 0) NA else vals[[1]]
            }
          ),
          .groups = "drop"
        )

      baseline_meta <- tibble(
        Variable = present_candidates,
        VariableLabel = present_candidates,
        Source = "Baseline covariate"
      )
    }
  }

  candidate_df <- subject_missing %>%
    left_join(modeled_candidate_df, by = id_var)
  if (!is.null(baseline_candidate_df)) {
    candidate_df <- candidate_df %>% left_join(baseline_candidate_df, by = id_var)
  }

  candidate_meta <- bind_rows(modeled_meta, baseline_meta) %>% distinct(Variable, .keep_all = TRUE)
  candidate_vars_all <- intersect(candidate_meta$Variable, names(candidate_df))
  if (length(candidate_vars_all) == 0) {
    return(tibble())
  }

  candidate_df <- candidate_df %>%
    mutate(across(all_of(candidate_vars_all), as.character))

  long_df <- candidate_df %>%
    pivot_longer(cols = all_of(candidate_vars_all), names_to = "Variable", values_to = "RawValue") %>%
    filter(!is.na(RawValue), RawValue != missing_code) %>%
    left_join(modeled_lookup, by = c("Variable", "RawValue")) %>%
    left_join(candidate_meta, by = "Variable") %>%
    mutate(Level = coalesce(Level, RawValue))

  if (nrow(long_df) == 0) {
    return(tibble())
  }

  long_df %>%
    group_by(Source, VariableLabel, Level) %>%
    summarise(
      N = n(),
      AnyMissingPct = mean(AnyMissing, na.rm = TRUE),
      MeanObservedPct = mean(ObservedPct, na.rm = TRUE),
      MeanObservedCount = mean(ObservedCount, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Source, VariableLabel, Level)
}

generate_missingness_correlates_table <- function(missingness_df, analysis_label) {
  if (nrow(missingness_df) == 0) {
    return(NULL)
  }

  missingness_df %>%
    gt(groupname_col = "Source") %>%
    tab_header(
      title = md(glue("**Table 6B: Missingness Correlates - {analysis_label}**")),
      subtitle = "Observed missingness by baseline-like candidate variables"
    ) %>%
    cols_label(
      VariableLabel = "Variable",
      Level = "Level",
      N = "N",
      AnyMissingPct = "Any Missing",
      MeanObservedPct = "Mean Observed Proportion",
      MeanObservedCount = "Mean Observed Count"
    ) %>%
    fmt_percent(columns = c("AnyMissingPct", "MeanObservedPct"), decimals = 1) %>%
    fmt_number(columns = "MeanObservedCount", decimals = 1) %>%
    tab_source_note("Note: Candidate variables are drawn from earliest-wave modeled indicators and available baseline covariates. The table is descriptive and is intended to document variables associated with incomplete modeled data.")
}

summarise_structural_sensitivity_reporting <- function(structural_sensitivity_audit) {
  summary_df <- structural_sensitivity_audit$summary %||% tibble()
  if (is.null(summary_df) || nrow(summary_df) == 0) {
    return(list(
      within_class_status = "Not reported",
      covariance_status = "Not reported",
      within_class_evidence = "No structural sensitivity audit was captured for the selected solution.",
      covariance_evidence = "No structural sensitivity audit was captured for the selected solution."
    ))
  }

  summary_df <- as_tibble(summary_df) %>%
    mutate(
      Baseline = as.logical(Baseline),
      AlternativeWithinClassHeterogeneity = as.logical(AlternativeWithinClassHeterogeneity),
      AlternativeCovarianceStructure = as.logical(AlternativeCovarianceStructure),
      estimated = Status %in% c("Estimated")
    )

  build_evidence <- function(rows, success_rows, fallback) {
    if (nrow(rows) == 0) {
      return(fallback)
    }

    attempted <- paste(unique(rows$ProfileLabel), collapse = "; ")
    estimated <- if (nrow(success_rows) > 0) {
      paste(unique(success_rows$ProfileLabel), collapse = "; ")
    } else {
      "none"
    }

    glue("Profiles attempted: {attempted}. Profiles with parsed output: {estimated}.")
  }

  within_rows <- summary_df %>% filter(!Baseline, AlternativeWithinClassHeterogeneity)
  within_success <- within_rows %>% filter(estimated)
  covariance_rows <- summary_df %>% filter(!Baseline, AlternativeCovarianceStructure)
  covariance_success <- covariance_rows %>% filter(estimated)

  list(
    within_class_status = case_when(
      nrow(within_success) > 0 ~ "Reported",
      nrow(within_rows) > 0 ~ "Partial",
      TRUE ~ "Not reported"
    ),
    covariance_status = case_when(
      nrow(covariance_success) > 0 ~ "Reported",
      nrow(covariance_rows) > 0 ~ "Partial",
      TRUE ~ "Not reported"
    ),
    within_class_evidence = build_evidence(
      within_rows,
      within_success,
      "No alternative within-class heterogeneity profiles were attempted."
    ),
    covariance_evidence = build_evidence(
      covariance_rows,
      covariance_success,
      "No alternative variance-covariance profiles were attempted."
    )
  )
}

generate_lcga_structural_sensitivity_table <- function(structural_sensitivity_audit, analysis_label) {
  summary_df <- structural_sensitivity_audit$summary %||% tibble()
  if (is.null(summary_df) || nrow(summary_df) == 0) {
    return(
      tibble(
        ProfileLabel = "Structural sensitivity audit",
        Baseline = "No",
        Description = "No structural sensitivity audit object was captured for this variant branch.",
        Status = "Unavailable",
        SelectedK = NA_integer_,
        BIC = NA_real_,
        Entropy = NA_real_,
        MinClassPct = NA_real_,
        Evidence = "This report variant was rendered without a structural sensitivity audit. The checklist is still shown, but alternative within-class and covariance-structure evidence is unavailable for this branch unless a dedicated structural audit was run."
      ) %>%
        gt() %>%
        tab_header(
          title = md(glue("**Table 7A: Structural Sensitivity Audit - {analysis_label}**")),
          subtitle = "Selected-solution alternative growth-factor variance and covariance structures"
        ) %>%
        cols_label(
          ProfileLabel = "Profile",
          Baseline = "Baseline",
          Description = "Structure",
          Status = "Run Status",
          SelectedK = "K",
          BIC = "BIC",
          Entropy = "Entropy",
          MinClassPct = "Min Class",
          Evidence = "Evidence"
        ) %>%
        fmt_number(columns = c("BIC"), decimals = 1) %>%
        fmt_number(columns = c("Entropy"), decimals = 3) %>%
        fmt_percent(columns = c("MinClassPct"), decimals = 1) %>%
        tab_source_note("Note: This placeholder table is emitted when the variant branch was reported without a captured structural sensitivity audit object.") %>%
        style_status_column("Status")
    )
  }

  summary_df %>%
    mutate(Baseline = ifelse(Baseline, "Yes", "No")) %>%
    select(
      ProfileLabel,
      Baseline,
      Description,
      Status,
      SelectedK,
      BIC,
      Entropy,
      MinClassPct,
      Evidence
    ) %>%
    gt() %>%
    tab_header(
      title = md(glue("**Table 7A: Structural Sensitivity Audit - {analysis_label}**")),
      subtitle = "Selected-solution alternative growth-factor variance and covariance structures"
    ) %>%
    cols_label(
      ProfileLabel = "Profile",
      Baseline = "Baseline",
      Description = "Structure",
      Status = "Run Status",
      SelectedK = "K",
      BIC = "BIC",
      Entropy = "Entropy",
      MinClassPct = "Min Class",
      Evidence = "Evidence"
    ) %>%
    fmt_number(columns = c("BIC"), decimals = 1) %>%
    fmt_number(columns = c("Entropy"), decimals = 3) %>%
    fmt_percent(columns = c("MinClassPct"), decimals = 1) %>%
    tab_source_note("Note: Structural sensitivity models are run at the selected K and are intended to document whether alternative within-class growth-factor variance and covariance assumptions materially change the selected-solution diagnostics.")
}

generate_lcga_grolts_checklist_table <- function(
    analysis_label,
    time_label,
    ages,
    time_var = "Time",
    outcome_meta,
    fit_df,
    best_k,
    class_assign_data,
    run_profile = NULL,
    syntax_dir = NULL,
    availability_available = FALSE,
    descriptive_available = FALSE,
    sensitivity_available = FALSE,
    estimates_available = FALSE,
    weighted_available = FALSE,
    final_plot_available = FALSE,
    per_model_plots_available = FALSE,
    combined_observed_plot_available = FALSE,
    missingness_correlates_available = FALSE,
    numerical_summary_available = FALSE,
    structural_sensitivity_audit = NULL,
    custom_labels_applied = FALSE,
    manual_override_key = NULL
) {
  has_fit_columns <- all(c("BIC", "Entropy", "MinClassPct") %in% names(fit_df))
  has_test_col <- any(c("BLRT_p", "VLMR_p") %in% names(fit_df))
  syntax_available <- !is.null(syntax_dir) && dir.exists(syntax_dir) && length(list.files(syntax_dir, pattern = "\\.inp$", full.names = TRUE)) > 0
  k_range <- if (nrow(fit_df) > 0) paste(range(fit_df$Classes, na.rm = TRUE), collapse = " to ") else "Unavailable"
  avepp_summary <- summarise_lcga_avepp(class_assign_data)
  fitted_model_count <- if (nrow(fit_df) > 0) nrow(fit_df) else 0
  age_mean <- mean(ages)
  age_var <- stats::var(ages)
  selected_min_class <- if ("MinClassPct" %in% names(fit_df)) (fit_df %>% filter(Classes == best_k) %>% slice(1) %>% pull(MinClassPct) %>% as.numeric())[1] else NA_real_
  entropy_value <- if ("Entropy" %in% names(fit_df)) (fit_df %>% filter(Classes == best_k) %>% slice(1) %>% pull(Entropy) %>% as.numeric())[1] else NA_real_
  comparison_tools <- intersect(c("BIC", "aBIC", "Entropy", "BLRT_p", "VLMR_p", "MinClassPct"), names(fit_df))
  structural_reporting <- summarise_structural_sensitivity_reporting(structural_sensitivity_audit)
  starts_evidence <- if (!is.null(run_profile) && !is.null(run_profile$starts)) {
    glue("STARTS = {run_profile$starts}; processors = {run_profile$processors %||% 'NA'}; stage = {run_profile$stage %||% 'NA'}")
  } else {
    "Execution settings not captured in report object"
  }
  class_size_evidence <- if (!is.null(class_assign_data)) {
    class_assign_data %>%
      count(Assigned_Class_Label, name = "n") %>%
      mutate(pct = n / sum(n)) %>%
      transmute(text = paste0(as.character(Assigned_Class_Label), ": ", n, " (", format_pct_string(pct, 1), ")")) %>%
      pull(text) %>%
      paste(collapse = "; ")
  } else {
    "Class counts unavailable"
  }

  checklist_df <- tibble(
    Domain = c(
      "Time metric", "Time metric",
      "Missing data", "Missing data", "Missing data",
      "Observed data", "Software",
      "Model specification", "Model specification", "Model specification", "Model specification",
      "Enumeration", "Enumeration", "Enumeration", "Enumeration", "Classification",
      "Visualization", "Visualization", "Visualization",
      "Reporting", "Reproducibility"
    ),
    Item = c(
      "Metric of time used in the model",
      "Mean and variance of time within the wave",
      "3a. Missing data mechanism reported",
      "3b. Variables related to attrition/missingness described",
      "3c. How missing data were handled in the analysis",
      "Distribution of observed variables",
      "Software used",
      "6a. Alternative specifications of within-class heterogeneity considered and justified",
      "6b. Alternative specifications of between-class variance-covariance structure considered and justified",
      "Alternative trajectory shapes / functional forms described",
      "If covariates were used, whether the analysis is still replicable",
      "Number of random starts and final iterations",
      "Model comparison / selection tools described statistically",
      "Total number of fitted models, including the 1-class solution",
      "Number or proportion of cases per class for each model",
      "Entropy, if classification is a goal",
      "14a. Plot of estimated mean trajectories for the final solution",
      "14b. Plots of estimated mean trajectories for each model",
      "14c. Plot combining estimated means with observed individual trajectories by class",
      "Numerical description of the final class solution, such as means, SE/SD, n, and CIs",
      "Availability of syntax files"
    ),
    Status = c(
      "Reported",
      "Reported",
      "Partial",
      if (isTRUE(missingness_correlates_available)) "Reported" else "Not reported",
      if (isTRUE(availability_available)) "Reported" else "Partial",
      if (isTRUE(descriptive_available)) "Reported" else "Not reported",
      "Reported",
      structural_reporting$within_class_status,
      structural_reporting$covariance_status,
      if (isTRUE(estimates_available)) "Partial" else "Not reported",
      "Not applicable",
      if (!is.null(run_profile) && !is.null(run_profile$starts)) "Partial" else "Not reported",
      if (has_fit_columns && has_test_col) "Reported" else "Partial",
      if (fitted_model_count > 0) "Reported" else "Not reported",
      if (!is.null(class_assign_data)) "Partial" else "Not reported",
      if (!is.na(entropy_value)) "Reported" else "Not reported",
      if (isTRUE(final_plot_available)) "Reported" else "Not reported",
      if (isTRUE(per_model_plots_available)) "Reported" else "Not reported",
      if (isTRUE(combined_observed_plot_available)) "Reported" else "Not reported",
      if (isTRUE(numerical_summary_available)) "Reported" else if (isTRUE(estimates_available)) "Partial" else "Not reported",
      if (isTRUE(syntax_available)) "Reported" else "Not reported"
    ),
    Evidence = c(
      glue("{time_var}; fitted ages = {paste(ages, collapse = ', ')}"),
      glue("Mean age = {format_num_string(age_mean, 2)}; variance = {format_num_string(age_var, 2)}"),
      if (isTRUE(sensitivity_available)) "Availability and sensitivity outputs document the extent of missingness, but the missing-data mechanism is not formally tested in the report." else "Missingness extent is visible from the available tables, but the missing-data mechanism is not formally tested in the report.",
      if (isTRUE(missingness_correlates_available)) "A missingness-correlates table summarises incomplete modeled data by candidate baseline variables and earliest-wave indicators." else "No explicit attrition or missingness-predictor table is included in the current LCGA report.",
      "Missing values are coded with the configured missing code and handled through Mplus estimation; the report documents availability rather than complete-case deletion.",
      if (isTRUE(descriptive_available) && length(outcome_meta) > 0) paste(purrr::map_chr(outcome_meta, "label"), collapse = "; ") else "Observed distribution table not currently available",
      "R, MplusAutomation, and Mplus are used in the current workflow.",
      structural_reporting$within_class_evidence,
      structural_reporting$covariance_evidence,
      if (isTRUE(estimates_available)) "Intercept, slope, and quadratic terms are reported for the selected model; alternative functional forms are not systematically compared in the report." else "No trajectory-shape evidence available",
      "Covariates are not used in the LCGA fitting stage, so this item is not applicable for the current report.",
      starts_evidence,
      if (length(comparison_tools) > 0) paste(comparison_tools, collapse = ", ") else "No model-comparison statistics available",
      glue("Estimated K values: {k_range}; selected K = {best_k}; fitted models counted = {fitted_model_count}"),
      if (!is.null(class_assign_data)) glue("Selected-model class counts are available: {class_size_evidence}. Non-selected models are summarized via fit statistics and minimum class size only.") else "Class counts unavailable",
      if (!is.na(entropy_value)) glue("Selected-model entropy = {format_num_string(entropy_value, 3)}; mean diagonal AvePP = {format_num_string(avepp_summary$mean_diag, 3)}") else "Entropy not available",
      if (isTRUE(final_plot_available)) "Final selected-solution trajectory plot generated." else "Final trajectory plot missing",
      if (isTRUE(per_model_plots_available)) "Per-K estimated trajectory plots generated." else "The current report focuses on the selected model and does not generate trajectory plots for every fitted K.",
      if (isTRUE(combined_observed_plot_available)) "Combined estimated-plus-observed individual trajectory plots by assigned class generated." else "A true estimated-plus-individual-trajectories-by-class plot is not currently generated.",
      if (isTRUE(numerical_summary_available)) glue("Selected-model class sizes, estimates, SEs, and 95% CIs are tabulated; selected-model min class = {format_pct_string(selected_min_class, 1)}; class labels customized = {ifelse(custom_labels_applied, 'yes', 'no')}.") else if (isTRUE(estimates_available)) glue("Selected-model estimates and SEs are reported; selected-model min class = {format_pct_string(selected_min_class, 1)}; class labels customized = {ifelse(custom_labels_applied, 'yes', 'no')}; CIs are not currently tabulated.") else "No selected-model numerical summary table generated",
      if (isTRUE(syntax_available)) syntax_dir else "No syntax directory detected"
    )
  )

  checklist_df %>%
    gt(groupname_col = "Domain") %>%
    tab_header(
      title = md(glue("**Table 8: GRoLTS-Oriented Reporting Checklist - {analysis_label}**")),
      subtitle = "Operational reporting checklist for the selected latent trajectory model"
    ) %>%
    cols_label(
      Item = "Checklist Item",
      Status = "Status",
      Evidence = "Evidence"
    ) %>%
    tab_source_note("Note: This is an operational project checklist aligned to GRoLTS-style reporting domains. It is intended to document whether the selected-model report includes the expected elements, not to certify formal compliance.") %>%
    style_status_column("Status")
}

generate_lca_smart_checklist_table <- function(
    analysis_label,
    indicator_vars,
    fit_summary_aug,
    best_k,
    report_tables,
    starts = NULL,
    syntax_dir = NULL,
    profile_plot_available = FALSE,
    class_labels = NULL,
    manual_override_key = NULL
) {
  has_fit_columns <- all(c("BIC", "Entropy", "MinClassPct") %in% names(fit_summary_aug))
  has_test_col <- any(c("BLRT_p", "VLMR_p") %in% names(fit_summary_aug))
  custom_labels_applied <- !is.null(class_labels) && any(!grepl("^Class\\s+[0-9]+$", class_labels))
  syntax_available <- !is.null(syntax_dir) && dir.exists(syntax_dir) && length(list.files(syntax_dir, pattern = "\\.inp$", full.names = TRUE)) > 0
  k_range <- if (nrow(fit_summary_aug) > 0) paste(range(fit_summary_aug$Classes, na.rm = TRUE), collapse = " to ") else "Unavailable"
  starts_evidence <- if (!is.null(starts)) glue("STARTS = {starts}") else "Execution settings not captured"
  selected_min_class <- if ("MinClassPct" %in% names(fit_summary_aug)) (fit_summary_aug %>% filter(Classes == best_k) %>% slice(1) %>% pull(MinClassPct) %>% as.numeric())[1] else NA_real_
  comparison_tools <- intersect(c("BIC", "aBIC", "Entropy", "BLRT_p", "VLMR_p", "MinClassPct"), names(fit_summary_aug))

  checklist_df <- tibble(
    Domain = c(
      "Study setup", "Observed data", "Preprocessing", "Missing data",
      "Model specification", "Class enumeration", "Class enumeration", "Class enumeration",
      "Transparency", "Estimation", "Reporting", "Inference", "Visualization", "Follow-up analyses"
    ),
    Item = c(
      "Pre-analysis",
      "Examining observed data",
      "Data preprocessing",
      "Missing data",
      "Model specification",
      "Number of classes",
      "Criteria for class enumeration",
      "Criteria for eliminating models from consideration",
      "Transparency and reproducibility",
      "Estimation and convergence",
      "Reporting results",
      "Inference",
      "Visualization",
      "Follow-up analyses, including accounting for classification inaccuracy"
    ),
    Status = c(
      if (length(indicator_vars) > 0) "Partial" else "Not reported",
      if (!is.null(report_tables$availability) && !is.null(report_tables$descriptive)) "Reported" else "Partial",
      "Reported",
      "Reported",
      if (length(indicator_vars) > 0) "Reported" else "Not reported",
      if (nrow(fit_summary_aug) > 0) "Reported" else "Not reported",
      if (has_fit_columns && has_test_col) "Reported" else "Partial",
      if ("MinClassPct" %in% names(fit_summary_aug)) "Reported" else "Not reported",
      if (isTRUE(syntax_available)) "Reported" else "Partial",
      if (!is.null(starts)) "Reported" else "Partial",
      if (!is.null(report_tables$fit) && !is.null(report_tables$avepp) && !is.null(report_tables$estimates)) "Reported" else "Partial",
      if (!is.null(report_tables$weighted)) "Partial" else "Not reported",
      if (isTRUE(profile_plot_available)) "Reported" else "Not reported",
      if (!is.null(report_tables$weighted)) "Partial" else "Not reported"
    ),
    Evidence = c(
      glue("Indicators prespecified for {analysis_label}: {paste(indicator_vars, collapse = ', ')}"),
      paste(c(if (!is.null(report_tables$availability)) "availability", if (!is.null(report_tables$descriptive)) "descriptive"), collapse = ", "),
      "Indicators are recoded through prep_lca_data() and written as categorical Mplus inputs.",
      "Configured missing code is passed to Mplus and data availability is reported explicitly.",
      glue("Selected model family = {analysis_label}; indicators = {length(indicator_vars)}"),
      glue("Estimated K values: {k_range}; selected K = {best_k}"),
      if (length(comparison_tools) > 0) paste(comparison_tools, collapse = ", ") else "No model-comparison statistics available",
      if ("MinClassPct" %in% names(fit_summary_aug)) glue("Selected K min class = {format_pct_string(selected_min_class, 1)}; models failing size thresholds are excluded from consideration.") else "MinClassPct not available",
      if (isTRUE(syntax_available)) glue("Syntax files available under {syntax_dir}; custom labels applied = {ifelse(custom_labels_applied, 'yes', 'no')}; manual override = {manual_override_key %||% 'none'}") else "Syntax files not detected",
      glue("{starts_evidence}; convergence is evaluated through the fitted output summaries."),
      paste(c(if (!is.null(report_tables$fit)) "fit", if (!is.null(report_tables$avepp)) "AvePP", if (!is.null(report_tables$estimates)) "item-response estimates"), collapse = ", "),
      if (!is.null(report_tables$weighted)) "Weighted class-characteristics table generated for follow-up description; formal model-based distal inference is not yet included in the report." else "Follow-up inference table unavailable",
      if (isTRUE(profile_plot_available)) "Selected-model profile plot generated." else "Profile plot unavailable",
      if (!is.null(report_tables$weighted)) "Weighted summaries use posterior-probability information for descriptive follow-up, but a full dedicated follow-up module is not yet included." else "Classification-uncertainty follow-up output unavailable"
    )
  )

  checklist_df %>%
    gt(groupname_col = "Domain") %>%
    tab_header(
      title = md(glue("**Table 9: SMART-LCA Reporting Checklist - {analysis_label}**")),
      subtitle = "Operational reporting checklist for the selected latent class model"
    ) %>%
    cols_label(
      Item = "Checklist Item",
      Status = "Status",
      Evidence = "Evidence"
    ) %>%
    tab_source_note("Note: This is an operational project checklist aligned to SMART-LCA-style reporting domains. It documents whether the selected-model report includes the expected elements, not whether a formal external checklist has been audited.") %>%
    style_status_column("Status")
}
