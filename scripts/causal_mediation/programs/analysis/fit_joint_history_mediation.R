#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(paste0("^", flag, "="), "", hit[[1]])
}

n_sim <- as.integer(arg_value("--sim", "300"))
block_sim <- as.integer(arg_value("--block-sim", as.character(min(n_sim, 120L))))
expanded_block_sim <- as.integer(arg_value("--expanded-block-sim", "0"))
seed <- as.integer(arg_value("--seed", "20260413"))
weight_floor <- as.numeric(arg_value("--weight-floor", "1e-10"))
branch <- arg_value("--branch", "strict")
subject_input_arg <- arg_value("--subject-input", NULL)
adjustment_config_arg <- arg_value("--adjustment-config", NULL)
output_prefix_arg <- arg_value("--output-prefix", NULL)

branch_specs <- list(
  strict = list(
    branch = "strict",
    input_prefix = "jmaes_strict",
    output_prefix = "strict_jmaes",
    scenario_label = "Stricter restricted (>=3 complete occasions)",
    exposure_solution = "Stricter restricted Joint Maternal Absence + Household Socioeconomic Status LCGA"
  ),
  restricted = list(
    branch = "restricted",
    input_prefix = "jmaes_restricted",
    output_prefix = "restricted_jmaes",
    scenario_label = "Restricted (>=3 occasions)",
    exposure_solution = "Restricted Joint Maternal Absence + Household Socioeconomic Status LCGA"
  ),
  unrestricted = list(
    branch = "unrestricted",
    input_prefix = "jmaes_unrestricted",
    output_prefix = "unrestricted_jmaes",
    scenario_label = "Unrestricted",
    exposure_solution = "Unrestricted Joint Maternal Absence + Household Socioeconomic Status LCGA"
  )
)

if (!branch %in% names(branch_specs)) {
  stop("Use --branch=strict, --branch=restricted, or --branch=unrestricted")
}

branch_spec <- branch_specs[[branch]]

set.seed(seed)

script_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][[1]])
script_arg <- gsub("~\\+~", " ", script_arg)
script_path <- normalizePath(script_arg, mustWork = TRUE)
workspace_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

input_dir <- file.path(workspace_root, "input")
config_dir <- file.path(workspace_root, "config")
listing_dir <- file.path(workspace_root, "output", "listings")
dir.create(listing_dir, recursive = TRUE, showWarnings = FALSE)
subject_input_path <- if (!is.null(subject_input_arg)) {
  normalizePath(subject_input_arg, mustWork = TRUE)
} else {
  file.path(input_dir, paste0(branch_spec$input_prefix, "_subject_level_analysis.rds"))
}
adjustment_config_path <- if (!is.null(adjustment_config_arg)) {
  normalizePath(adjustment_config_arg, mustWork = TRUE)
} else {
  file.path(config_dir, "adjustment_sets.csv")
}
output_prefix <- if (!is.null(output_prefix_arg)) output_prefix_arg else branch_spec$output_prefix

format_rd_ci <- function(rd, low, high) {
  if (!is.finite(rd)) {
    return("")
  }
  if (!is.finite(low) || !is.finite(high)) {
    return(sprintf("%.3f", rd))
  }
  sprintf("%.3f (%.3f, %.3f)", rd, low, high)
}

format_percent <- function(x) {
  if (!is.finite(x)) {
    return("")
  }
  sprintf("%.1f%%", 100 * x)
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  sum(w[ok] * x[ok]) / sum(w[ok])
}

safe_vcov <- function(fit, data) {
  vc <- tryCatch(
    sandwich::vcovCL(fit, cluster = data$USUBJID, type = "HC0"),
    error = function(e) vcov(fit)
  )
  vc <- as.matrix(vc)
  vc[!is.finite(vc)] <- 0
  vc <- (vc + t(vc)) / 2
  eig <- eigen(vc, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- suppressWarnings(min(eig, na.rm = TRUE))
  if (!is.finite(min_eig)) {
    vc <- diag(1e-8, nrow(vc))
    rownames(vc) <- colnames(vc) <- names(coef(fit))
  } else if (min_eig < 1e-10) {
    vc <- vc + diag(abs(min_eig) + 1e-8, nrow(vc))
  }
  vc
}

model_draws <- function(fit, data, n) {
  beta <- coef(fit)
  keep <- is.finite(beta)
  beta <- beta[keep]
  vc <- safe_vcov(fit, data)
  vc <- vc[names(beta), names(beta), drop = FALSE]
  if (n <= 0) {
    draws <- matrix(numeric(0), nrow = 0, ncol = length(beta))
    colnames(draws) <- names(beta)
    return(draws)
  }
  draws <- MASS::mvrnorm(n = n, mu = beta, Sigma = vc)
  draws <- as.matrix(draws)
  if (n == 1) {
    draws <- matrix(draws, nrow = 1)
  }
  colnames(draws) <- names(beta)
  draws
}

model_matrix_for_draws <- function(fit, data, draws) {
  terms_obj <- delete.response(terms(fit))
  mm <- model.matrix(terms_obj, data = data)
  draw_names <- colnames(draws)
  missing_cols <- setdiff(draw_names, colnames(mm))
  if (length(missing_cols) > 0) {
    for (missing_col in missing_cols) {
      mm <- cbind(mm, 0)
      colnames(mm)[ncol(mm)] <- missing_col
    }
  }
  mm[, draw_names, drop = FALSE]
}

predict_mean_draws <- function(fit, data, draws, weights) {
  if (nrow(draws) == 0) {
    return(numeric(0))
  }
  mm <- model_matrix_for_draws(fit, data, draws)
  vapply(
    seq_len(nrow(draws)),
    function(i) weighted_mean(plogis(drop(mm %*% draws[i, ])), weights),
    numeric(1)
  )
}

predict_prob_draw_matrix <- function(fit, data, draws) {
  if (nrow(draws) == 0) {
    return(matrix(numeric(0), nrow = nrow(data), ncol = 0))
  }
  mm <- model_matrix_for_draws(fit, data, draws)
  eta <- mm %*% t(draws)
  plogis(eta)
}

set_exposure <- function(data, exposure_var, value) {
  out <- data
  out[[exposure_var]] <- factor(value, levels = levels(data[[exposure_var]]))
  out
}

set_mediator <- function(data, mediator_var, value) {
  out <- data
  out[[mediator_var]] <- value
  out
}

set_mediator_pattern <- function(data, mediator_vars, values) {
  out <- data
  for (i in seq_along(mediator_vars)) {
    out[[mediator_vars[[i]]]] <- values[[i]]
  }
  out
}

weighted_col_mean <- function(mat, w) {
  ok <- is.finite(w) & w > 0
  if (!any(ok)) {
    return(rep(NA_real_, ncol(mat)))
  }
  colSums(sweep(mat[ok, , drop = FALSE], 1, w[ok], `*`), na.rm = TRUE) / sum(w[ok])
}

fit_glm_safe <- function(formula, data) {
  warnings <- character()
  fit <- withCallingHandlers(
    glm(
      formula = formula,
      data = data,
      weights = posterior_weight,
      family = quasibinomial(),
      control = glm.control(maxit = 100)
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  attr(fit, "fit_warnings") <- unique(warnings)
  fit
}

make_formula <- function(lhs, rhs_terms) {
  as.formula(paste(lhs, "~", paste(rhs_terms, collapse = " + ")))
}

complete_analysis_data <- function(data, exposure_var, covariates, extra_vars = character()) {
  needed <- unique(c("USUBJID", "posterior_weight", "DPBN", exposure_var, covariates, extra_vars))
  data %>%
    filter(if_all(all_of(needed), ~ !is.na(.x))) %>%
    filter(is.finite(posterior_weight), posterior_weight > weight_floor) %>%
    droplevels()
}

estimate_total_effect <- function(fit, model_data, exposure_var, reference, comparator, draws) {
  ref_data <- set_exposure(model_data, exposure_var, reference)
  comp_data <- set_exposure(model_data, exposure_var, comparator)
  weights <- model_data$posterior_weight

  p_ref <- predict(fit, newdata = ref_data, type = "response")
  p_comp <- predict(fit, newdata = comp_data, type = "response")
  mean_p_ref <- weighted_mean(p_ref, weights)
  mean_p_comp <- weighted_mean(p_comp, weights)
  rd <- weighted_mean(p_comp - p_ref, weights)

  ref_draw <- predict_mean_draws(fit, ref_data, draws, weights)
  comp_draw <- predict_mean_draws(fit, comp_data, draws, weights)
  rd_draw <- comp_draw - ref_draw
  ci <- if (length(rd_draw) > 1) {
    unname(quantile(rd_draw, c(0.025, 0.975), na.rm = TRUE))
  } else {
    c(NA_real_, NA_real_)
  }

  tibble(
    p_ref = mean_p_ref,
    p_comp = mean_p_comp,
    rd = rd,
    ci_low = ci[[1]],
    ci_high = ci[[2]],
    within_variance = if (length(rd_draw) > 1) stats::var(rd_draw, na.rm = TRUE) else NA_real_,
    rd_ci = format_rd_ci(rd, ci[[1]], ci[[2]])
  )
}

estimate_single_mediator <- function(mediator_fit,
                                     outcome_fit,
                                     model_data,
                                     exposure_var,
                                     mediator_var,
                                     reference,
                                     comparator,
                                     mediator_draws,
                                     outcome_draws) {
  weights <- model_data$posterior_weight

  ref_data <- set_exposure(model_data, exposure_var, reference)
  comp_data <- set_exposure(model_data, exposure_var, comparator)

  mediator_ref <- predict(mediator_fit, newdata = ref_data, type = "response")
  mediator_comp <- predict(mediator_fit, newdata = comp_data, type = "response")

  y_ref_m0 <- predict(outcome_fit, newdata = set_mediator(ref_data, mediator_var, 0), type = "response")
  y_ref_m1 <- predict(outcome_fit, newdata = set_mediator(ref_data, mediator_var, 1), type = "response")
  y_comp_m0 <- predict(outcome_fit, newdata = set_mediator(comp_data, mediator_var, 0), type = "response")
  y_comp_m1 <- predict(outcome_fit, newdata = set_mediator(comp_data, mediator_var, 1), type = "response")

  psi_ref_ref <- weighted_mean((1 - mediator_ref) * y_ref_m0 + mediator_ref * y_ref_m1, weights)
  psi_comp_ref <- weighted_mean((1 - mediator_ref) * y_comp_m0 + mediator_ref * y_comp_m1, weights)
  psi_comp_comp <- weighted_mean((1 - mediator_comp) * y_comp_m0 + mediator_comp * y_comp_m1, weights)

  point <- c(
    `Total effect` = psi_comp_comp - psi_ref_ref,
    `Direct effect` = psi_comp_ref - psi_ref_ref,
    `Indirect effect` = psi_comp_comp - psi_comp_ref
  )

  if (nrow(mediator_draws) > 1 && nrow(outcome_draws) > 1) {
    mediator_ref_draw <- predict_prob_draw_matrix(mediator_fit, ref_data, mediator_draws)
    mediator_comp_draw <- predict_prob_draw_matrix(mediator_fit, comp_data, mediator_draws)
    y_ref_m0_draw <- predict_prob_draw_matrix(outcome_fit, set_mediator(ref_data, mediator_var, 0), outcome_draws)
    y_ref_m1_draw <- predict_prob_draw_matrix(outcome_fit, set_mediator(ref_data, mediator_var, 1), outcome_draws)
    y_comp_m0_draw <- predict_prob_draw_matrix(outcome_fit, set_mediator(comp_data, mediator_var, 0), outcome_draws)
    y_comp_m1_draw <- predict_prob_draw_matrix(outcome_fit, set_mediator(comp_data, mediator_var, 1), outcome_draws)

    sim_count <- min(ncol(mediator_ref_draw), ncol(y_ref_m0_draw))
    draw_effects <- vapply(
      seq_len(sim_count),
      function(i) {
        pm_ref <- mediator_ref_draw[, i]
        pm_comp <- mediator_comp_draw[, i]
        psi00 <- weighted_mean((1 - pm_ref) * y_ref_m0_draw[, i] + pm_ref * y_ref_m1_draw[, i], weights)
        psi10 <- weighted_mean((1 - pm_ref) * y_comp_m0_draw[, i] + pm_ref * y_comp_m1_draw[, i], weights)
        psi11 <- weighted_mean((1 - pm_comp) * y_comp_m0_draw[, i] + pm_comp * y_comp_m1_draw[, i], weights)
        c(
          `Total effect` = psi11 - psi00,
          `Direct effect` = psi10 - psi00,
          `Indirect effect` = psi11 - psi10
        )
      },
      numeric(3)
    )
  } else {
    draw_effects <- matrix(NA_real_, nrow = 3, ncol = 0, dimnames = list(names(point), NULL))
  }

  total_rd <- unname(point[["Total effect"]])

  tibble(
    effect = names(point),
    rd = unname(point),
    ci_low = map_dbl(
      names(point),
      ~ if (ncol(draw_effects) > 1) unname(quantile(draw_effects[.x, ], 0.025, na.rm = TRUE)) else NA_real_
    ),
    ci_high = map_dbl(
      names(point),
      ~ if (ncol(draw_effects) > 1) unname(quantile(draw_effects[.x, ], 0.975, na.rm = TRUE)) else NA_real_
    ),
    within_variance = map_dbl(
      names(point),
      ~ if (ncol(draw_effects) > 1) stats::var(draw_effects[.x, ], na.rm = TRUE) else NA_real_
    ),
    percent_mediated = if_else(
      effect == "Indirect effect" &
        is.finite(total_rd) &
        abs(total_rd) >= 0.01 &
        sign(rd) == sign(total_rd),
      rd / total_rd,
      NA_real_
    ),
    rd_ci = pmap_chr(list(rd, ci_low, ci_high), format_rd_ci),
    percent_mediated_display = map_chr(percent_mediated, format_percent)
  )
}

estimate_joint_mediator_block <- function(mediator_fits,
                                          outcome_fit,
                                          model_data,
                                          exposure_var,
                                          mediator_vars,
                                          reference,
                                          comparator,
                                          mediator_draws,
                                          outcome_draws) {
  weights <- model_data$posterior_weight
  patterns <- expand.grid(
    rep(list(c(0, 1)), length(mediator_vars)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  names(patterns) <- mediator_vars

  psi_point <- function(exposure_value, mediator_distribution_value) {
    psi <- rep(0, nrow(model_data))
    dist_base <- set_exposure(model_data, exposure_var, mediator_distribution_value)
    outcome_base <- set_exposure(model_data, exposure_var, exposure_value)

    for (pattern_i in seq_len(nrow(patterns))) {
      values <- as.numeric(patterns[pattern_i, mediator_vars, drop = TRUE])
      pattern_prob <- rep(1, nrow(model_data))

      for (k in seq_along(mediator_vars)) {
        previous_vars <- mediator_vars[seq_len(k - 1)]
        pred_data <- if (length(previous_vars) > 0) {
          set_mediator_pattern(dist_base, previous_vars, values[seq_len(k - 1)])
        } else {
          dist_base
        }
        pm <- predict(mediator_fits[[k]], newdata = pred_data, type = "response")
        pattern_prob <- pattern_prob * if (values[[k]] == 1) pm else (1 - pm)
      }

      outcome_data <- set_mediator_pattern(outcome_base, mediator_vars, values)
      mu <- predict(outcome_fit, newdata = outcome_data, type = "response")
      psi <- psi + pattern_prob * mu
    }

    weighted_mean(psi, weights)
  }

  psi_draw <- function(exposure_value, mediator_distribution_value) {
    sim_count <- min(c(nrow(outcome_draws), vapply(mediator_draws, nrow, integer(1))))
    if (!is.finite(sim_count) || sim_count <= 1) {
      return(numeric(0))
    }

    mediator_draws_i <- lapply(mediator_draws, function(x) x[seq_len(sim_count), , drop = FALSE])
    outcome_draws_i <- outcome_draws[seq_len(sim_count), , drop = FALSE]
    psi <- rep(0, sim_count)
    dist_base <- set_exposure(model_data, exposure_var, mediator_distribution_value)
    outcome_base <- set_exposure(model_data, exposure_var, exposure_value)

    for (pattern_i in seq_len(nrow(patterns))) {
      values <- as.numeric(patterns[pattern_i, mediator_vars, drop = TRUE])
      pattern_prob <- matrix(1, nrow = nrow(model_data), ncol = sim_count)

      for (k in seq_along(mediator_vars)) {
        previous_vars <- mediator_vars[seq_len(k - 1)]
        pred_data <- if (length(previous_vars) > 0) {
          set_mediator_pattern(dist_base, previous_vars, values[seq_len(k - 1)])
        } else {
          dist_base
        }
        pm <- predict_prob_draw_matrix(mediator_fits[[k]], pred_data, mediator_draws_i[[k]])
        pattern_prob <- pattern_prob * if (values[[k]] == 1) pm else (1 - pm)
      }

      outcome_data <- set_mediator_pattern(outcome_base, mediator_vars, values)
      mu <- predict_prob_draw_matrix(outcome_fit, outcome_data, outcome_draws_i)
      psi <- psi + weighted_col_mean(pattern_prob * mu, weights)
    }

    psi
  }

  psi_ref_ref <- psi_point(reference, reference)
  psi_comp_ref <- psi_point(comparator, reference)
  psi_comp_comp <- psi_point(comparator, comparator)

  point <- c(
    `Total effect` = psi_comp_comp - psi_ref_ref,
    `Block direct effect` = psi_comp_ref - psi_ref_ref,
    `Total indirect effect` = psi_comp_comp - psi_comp_ref
  )

  psi00_draw <- psi_draw(reference, reference)
  psi10_draw <- psi_draw(comparator, reference)
  psi11_draw <- psi_draw(comparator, comparator)

  draw_effects <- if (length(psi00_draw) > 1) {
    rbind(
      `Total effect` = psi11_draw - psi00_draw,
      `Block direct effect` = psi10_draw - psi00_draw,
      `Total indirect effect` = psi11_draw - psi10_draw
    )
  } else {
    matrix(NA_real_, nrow = 3, ncol = 0, dimnames = list(names(point), NULL))
  }

  total_rd <- unname(point[["Total effect"]])

  tibble(
    effect = names(point),
    rd = unname(point),
    ci_low = map_dbl(
      names(point),
      ~ if (ncol(draw_effects) > 1) unname(quantile(draw_effects[.x, ], 0.025, na.rm = TRUE)) else NA_real_
    ),
    ci_high = map_dbl(
      names(point),
      ~ if (ncol(draw_effects) > 1) unname(quantile(draw_effects[.x, ], 0.975, na.rm = TRUE)) else NA_real_
    ),
    within_variance = map_dbl(
      names(point),
      ~ if (ncol(draw_effects) > 1) stats::var(draw_effects[.x, ], na.rm = TRUE) else NA_real_
    ),
    percent_mediated = if_else(
      effect == "Total indirect effect" &
        is.finite(total_rd) &
        abs(total_rd) >= 0.01 &
        sign(rd) == sign(total_rd),
      rd / total_rd,
      NA_real_
    ),
    rd_ci = pmap_chr(list(rd, ci_low, ci_high), format_rd_ci),
    percent_mediated_display = map_chr(percent_mediated, format_percent)
  )
}

effect_flag <- function(rd, low, high) {
  if (!is.finite(rd)) {
    return("Not estimable")
  }
  direction <- if (rd > 0) "positive" else if (rd < 0) "negative" else "null"
  interval <- if (is.finite(low) && is.finite(high) && (low > 0 || high < 0)) {
    "interval excludes 0"
  } else if (is.finite(low) && is.finite(high)) {
    "interval includes 0"
  } else {
    "interval unavailable"
  }
  paste(direction, interval, sep = "; ")
}

subject <- readRDS(subject_input_path) %>%
  mutate(
    USUBJID = as.character(USUBJID),
    SEX = factor(SEX),
    RURAL = factor(RURAL),
    across(c(DNKA, ESXC, FDSC, GOVG, VLNC, SCSP, SCHL, DPBN), as.numeric)
  )

transition_bundle <- readRDS(file.path(input_dir, paste0(branch_spec$input_prefix, "_transition_bundle.rds")))
class_labels <- read_csv(file.path(config_dir, paste0(branch_spec$input_prefix, "_class_labels.csv")), show_col_types = FALSE)
strict_class_labels <- read_csv(file.path(config_dir, "jmaes_strict_class_labels.csv"), show_col_types = FALSE)
strict_total_effect_hypotheses <- read_csv(file.path(config_dir, "total_effect_hypotheses.csv"), show_col_types = FALSE) %>%
  filter(model_family == "JMAES") %>%
  separate_rows(comparator_histories, sep = ";\\s*")
pathway_variables <- read_csv(file.path(config_dir, "pathway_variables.csv"), show_col_types = FALSE)
adjustment_sets <- read_csv(adjustment_config_path, show_col_types = FALSE) %>%
  mutate(covariate_list = strsplit(covariates, ";\\s*"))

strict_family_lookup <- strict_class_labels %>%
  transmute(
    window = .data$window,
    strict_reporting_label = .data$reporting_label,
    harmonized_family_id = .data$harmonized_family_id
  )

branch_family_lookup <- class_labels %>%
  transmute(
    window = .data$window,
    harmonized_family_id = .data$harmonized_family_id,
    branch_reporting_label = .data$reporting_label
  )

branch_total_effect_hypotheses <- strict_total_effect_hypotheses %>%
  left_join(
    strict_family_lookup %>%
      rename(reference_history = strict_reporting_label, reference_family_id = harmonized_family_id),
    by = c("window", "reference_history")
  ) %>%
  left_join(
    strict_family_lookup %>%
      rename(comparator_histories = strict_reporting_label, comparator_family_id = harmonized_family_id),
    by = c("window", "comparator_histories")
  ) %>%
  left_join(
    branch_family_lookup %>%
      rename(reference_family_id = harmonized_family_id, branch_reference = branch_reporting_label),
    by = c("window", "reference_family_id")
  ) %>%
  left_join(
    branch_family_lookup %>%
      rename(comparator_family_id = harmonized_family_id, branch_comparator = branch_reporting_label),
    by = c("window", "comparator_family_id")
  ) %>%
  filter(!is.na(.data$reference_family_id), !is.na(.data$comparator_family_id)) %>%
  filter(!is.na(.data$branch_reference), !is.na(.data$branch_comparator)) %>%
  mutate(
    scenario = branch_spec$scenario_label,
    model_label = branch_spec$exposure_solution,
    analysis_role = paste0(branch_spec$exposure_solution, " total-effect contrast"),
    reference_history = .data$branch_reference,
    comparator_histories = .data$branch_comparator
  ) %>%
  select(names(strict_total_effect_hypotheses), reference_family_id, comparator_family_id) %>%
  distinct()

if (nrow(branch_total_effect_hypotheses) == 0) {
  stop("No branch-valid JMAES total-effect contrasts could be derived for ", branch, ".")
}

block_specs <- list(
  list(
    block_id = "base_pathway_block",
    block_label = "Base pathway block",
    mediator_vars = pathway_variables %>% filter(base_model) %>% pull(variable),
    mediator_labels = pathway_variables %>% filter(base_model) %>% pull(label)
  ),
  list(
    block_id = "expanded_pathway_block",
    block_label = "Expanded pathway block",
    mediator_vars = pathway_variables %>% filter(extension_model) %>% pull(variable),
    mediator_labels = pathway_variables %>% filter(extension_model) %>% pull(label)
  )
)

early_levels <- class_labels %>%
  filter(window == "0-5 years") %>%
  arrange(adversity_order) %>%
  pull(reporting_label)
mid_levels <- class_labels %>%
  filter(window == "6-12 years") %>%
  arrange(adversity_order) %>%
  pull(reporting_label)

prob_cols <- grep("^CPROB", names(transition_bundle$saved_probs), value = TRUE)

pseudo_data <- transition_bundle$saved_probs %>%
  transmute(USUBJID = as.character(USUBJID), across(all_of(prob_cols), as.numeric)) %>%
  pivot_longer(all_of(prob_cols), names_to = "cprob_col", values_to = "posterior_weight_raw") %>%
  filter(is.finite(posterior_weight_raw), posterior_weight_raw > weight_floor) %>%
  group_by(USUBJID) %>%
  mutate(
    posterior_row_total_raw = sum(posterior_weight_raw, na.rm = TRUE),
    posterior_weight = posterior_weight_raw / posterior_row_total_raw
  ) %>%
  ungroup() %>%
  left_join(transition_bundle$joint_map, by = "cprob_col") %>%
  mutate(
    early_exposure = factor(early_reporting_label, levels = early_levels),
    mid_exposure = factor(mid_reporting_label, levels = mid_levels)
  ) %>%
  inner_join(subject, by = "USUBJID")

validate_labels <- branch_total_effect_hypotheses %>%
  mutate(
    exposure_var = if_else(window == "0-5 years", "early_exposure", "mid_exposure"),
    valid_reference = if_else(
      exposure_var == "early_exposure",
      reference_history %in% early_levels,
      reference_history %in% mid_levels
    ),
    valid_comparator = if_else(
      exposure_var == "early_exposure",
      comparator_histories %in% early_levels,
      comparator_histories %in% mid_levels
    )
  )

if (any(!validate_labels$valid_reference | !validate_labels$valid_comparator)) {
  bad <- validate_labels %>%
    filter(!valid_reference | !valid_comparator) %>%
    select(hypothesis_id, window, reference_history, comparator_histories, valid_reference, valid_comparator)
  stop(
    "Some ", branch, " JMAES contrast labels are not present in the transition bundle:\n",
    paste(capture.output(print(bad, n = Inf)), collapse = "\n")
  )
}

total_results <- list()
mediator_results <- list()
joint_block_results <- list()
model_diagnostics <- list()

for (adjustment_i in seq_len(nrow(adjustment_sets))) {
  adjustment_set <- adjustment_sets$adjustment_set[[adjustment_i]]
  adjustment_label <- adjustment_sets$label[[adjustment_i]]
  covariates <- adjustment_sets$covariate_list[[adjustment_i]]

  for (window_value in unique(branch_total_effect_hypotheses$window)) {
    exposure_var <- if (window_value == "0-5 years") "early_exposure" else "mid_exposure"
    window_hypotheses <- branch_total_effect_hypotheses %>%
      filter(window == window_value)

    total_model_data <- complete_analysis_data(pseudo_data, exposure_var, covariates)
    total_fit <- fit_glm_safe(make_formula("DPBN", c(exposure_var, covariates)), total_model_data)
    total_draws <- model_draws(total_fit, total_model_data, n_sim)

    model_diagnostics[[length(model_diagnostics) + 1]] <- tibble(
      adjustment_set = adjustment_set,
      window = window_value,
      mediator = NA_character_,
      model_type = "Total-effect outcome model",
      formula = deparse(formula(total_fit)),
      n_subjects = n_distinct(total_model_data$USUBJID),
      n_pseudo_rows = nrow(total_model_data),
      weighted_outcome_events = sum(total_model_data$posterior_weight * total_model_data$DPBN, na.rm = TRUE),
      weighted_mediator_events = NA_real_,
      converged = isTRUE(total_fit$converged),
      warnings = paste(attr(total_fit, "fit_warnings"), collapse = " | ")
    )

    for (contrast_i in seq_len(nrow(window_hypotheses))) {
      hypothesis <- window_hypotheses[contrast_i, ]
      te <- estimate_total_effect(
        fit = total_fit,
        model_data = total_model_data,
        exposure_var = exposure_var,
        reference = hypothesis$reference_history,
        comparator = hypothesis$comparator_histories,
        draws = total_draws
      )

      total_results[[length(total_results) + 1]] <- te %>%
        mutate(
          adjustment_set = adjustment_set,
          adjustment_label = adjustment_label,
          tier = hypothesis$tier,
          hypothesis_id = hypothesis$hypothesis_id,
          window = hypothesis$window,
          scenario = hypothesis$scenario,
          reference = hypothesis$reference_history,
          comparator = hypothesis$comparator_histories,
          contrast = paste(comparator, "vs", reference),
          contrast_direction = "Comparator minus reference",
          expected_direction = hypothesis$expected_total_effect_direction,
          n_subjects = n_distinct(total_model_data$USUBJID),
          n_pseudo_rows = nrow(total_model_data),
          interval_method = paste0("cluster-robust coefficient simulation, ", n_sim, " draws"),
          proceed_to_decomposition = if_else(
            is.finite(ci_low) & is.finite(ci_high) & (ci_low > 0 | ci_high < 0),
            "Prespecified decomposition; total-effect interval excludes 0",
            "Prespecified decomposition; interpret descriptively"
          ),
          interpretation_flag = pmap_chr(list(rd, ci_low, ci_high), effect_flag)
        ) %>%
        select(
          adjustment_set, adjustment_label, tier, hypothesis_id, window, scenario,
          reference, comparator, contrast, contrast_direction, expected_direction,
          n_subjects, n_pseudo_rows, p_ref, p_comp, rd, ci_low, ci_high, rd_ci,
          within_variance,
          proceed_to_decomposition, interpretation_flag, interval_method
        )
    }

    for (mediator_i in seq_len(nrow(pathway_variables))) {
      mediator_var <- pathway_variables$variable[[mediator_i]]
      mediator_label <- pathway_variables$label[[mediator_i]]
      mediator_scope <- if (isTRUE(pathway_variables$base_model[[mediator_i]])) "Base mediator" else "Extension mediator"

      mediator_model_data <- complete_analysis_data(
        pseudo_data,
        exposure_var,
        covariates,
        extra_vars = mediator_var
      )

      mediator_fit <- fit_glm_safe(
        make_formula(mediator_var, c(exposure_var, covariates)),
        mediator_model_data
      )
      outcome_fit <- fit_glm_safe(
        make_formula("DPBN", c(exposure_var, mediator_var, covariates)),
        mediator_model_data
      )
      mediator_draws <- model_draws(mediator_fit, mediator_model_data, n_sim)
      outcome_draws <- model_draws(outcome_fit, mediator_model_data, n_sim)

      model_diagnostics[[length(model_diagnostics) + 1]] <- tibble(
        adjustment_set = adjustment_set,
        window = window_value,
        mediator = mediator_label,
        model_type = "Mediator model",
        formula = deparse(formula(mediator_fit)),
        n_subjects = n_distinct(mediator_model_data$USUBJID),
        n_pseudo_rows = nrow(mediator_model_data),
        weighted_outcome_events = sum(mediator_model_data$posterior_weight * mediator_model_data$DPBN, na.rm = TRUE),
        weighted_mediator_events = sum(mediator_model_data$posterior_weight * mediator_model_data[[mediator_var]], na.rm = TRUE),
        converged = isTRUE(mediator_fit$converged),
        warnings = paste(attr(mediator_fit, "fit_warnings"), collapse = " | ")
      )
      model_diagnostics[[length(model_diagnostics) + 1]] <- tibble(
        adjustment_set = adjustment_set,
        window = window_value,
        mediator = mediator_label,
        model_type = "Mediator-adjusted outcome model",
        formula = deparse(formula(outcome_fit)),
        n_subjects = n_distinct(mediator_model_data$USUBJID),
        n_pseudo_rows = nrow(mediator_model_data),
        weighted_outcome_events = sum(mediator_model_data$posterior_weight * mediator_model_data$DPBN, na.rm = TRUE),
        weighted_mediator_events = sum(mediator_model_data$posterior_weight * mediator_model_data[[mediator_var]], na.rm = TRUE),
        converged = isTRUE(outcome_fit$converged),
        warnings = paste(attr(outcome_fit, "fit_warnings"), collapse = " | ")
      )

      for (contrast_i in seq_len(nrow(window_hypotheses))) {
        hypothesis <- window_hypotheses[contrast_i, ]
        decomp <- estimate_single_mediator(
          mediator_fit = mediator_fit,
          outcome_fit = outcome_fit,
          model_data = mediator_model_data,
          exposure_var = exposure_var,
          mediator_var = mediator_var,
          reference = hypothesis$reference_history,
          comparator = hypothesis$comparator_histories,
          mediator_draws = mediator_draws,
          outcome_draws = outcome_draws
        )

        mediator_results[[length(mediator_results) + 1]] <- decomp %>%
          mutate(
            adjustment_set = adjustment_set,
            adjustment_label = adjustment_label,
            tier = hypothesis$tier,
            hypothesis_id = hypothesis$hypothesis_id,
            window = hypothesis$window,
            scenario = hypothesis$scenario,
            reference = hypothesis$reference_history,
            comparator = hypothesis$comparator_histories,
            contrast = paste(comparator, "vs", reference),
            contrast_direction = "Comparator minus reference",
            expected_direction = hypothesis$expected_total_effect_direction,
            mediator_variable = mediator_var,
            mediator = mediator_label,
            mediator_scope = mediator_scope,
            n_subjects = n_distinct(mediator_model_data$USUBJID),
            n_pseudo_rows = nrow(mediator_model_data),
            interval_method = paste0("cluster-robust coefficient simulation, ", n_sim, " independent mediator/outcome draws"),
            interpretation_flag = pmap_chr(list(rd, ci_low, ci_high), effect_flag)
          ) %>%
          select(
            adjustment_set, adjustment_label, tier, hypothesis_id, window, scenario,
            reference, comparator, contrast, contrast_direction, expected_direction,
            mediator_scope, mediator_variable, mediator, effect,
            n_subjects, n_pseudo_rows, rd, ci_low, ci_high, rd_ci,
            within_variance,
            percent_mediated, percent_mediated_display,
            interpretation_flag, interval_method
          )
      }
    }

    for (block_spec in block_specs) {
      mediator_vars <- block_spec$mediator_vars
      mediator_labels <- block_spec$mediator_labels
      block_label <- block_spec$block_label
      block_id <- block_spec$block_id
      current_block_sim <- if (identical(block_id, "expanded_pathway_block")) {
        expanded_block_sim
      } else {
        block_sim
      }

      block_model_data <- complete_analysis_data(
        pseudo_data,
        exposure_var,
        covariates,
        extra_vars = mediator_vars
      )

      block_mediator_fits <- vector("list", length(mediator_vars))
      block_mediator_draws <- vector("list", length(mediator_vars))

      for (k in seq_along(mediator_vars)) {
        previous_vars <- mediator_vars[seq_len(k - 1)]
        rhs_terms <- c(exposure_var, covariates, previous_vars)
        block_mediator_fits[[k]] <- fit_glm_safe(
          make_formula(mediator_vars[[k]], rhs_terms),
          block_model_data
        )
        block_mediator_draws[[k]] <- model_draws(block_mediator_fits[[k]], block_model_data, current_block_sim)

        model_diagnostics[[length(model_diagnostics) + 1]] <- tibble(
          adjustment_set = adjustment_set,
          window = window_value,
          mediator = paste0(block_label, ": ", mediator_labels[[k]]),
          model_type = "Joint-block sequential mediator model",
          formula = deparse(formula(block_mediator_fits[[k]])),
          n_subjects = n_distinct(block_model_data$USUBJID),
          n_pseudo_rows = nrow(block_model_data),
          weighted_outcome_events = sum(block_model_data$posterior_weight * block_model_data$DPBN, na.rm = TRUE),
          weighted_mediator_events = sum(block_model_data$posterior_weight * block_model_data[[mediator_vars[[k]]]], na.rm = TRUE),
          converged = isTRUE(block_mediator_fits[[k]]$converged),
          warnings = paste(attr(block_mediator_fits[[k]], "fit_warnings"), collapse = " | ")
        )
      }

      block_outcome_fit <- fit_glm_safe(
        make_formula("DPBN", c(exposure_var, mediator_vars, covariates)),
        block_model_data
      )
      block_outcome_draws <- model_draws(block_outcome_fit, block_model_data, current_block_sim)

      model_diagnostics[[length(model_diagnostics) + 1]] <- tibble(
        adjustment_set = adjustment_set,
        window = window_value,
        mediator = block_label,
        model_type = "Joint-block mediator-adjusted outcome model",
        formula = deparse(formula(block_outcome_fit)),
        n_subjects = n_distinct(block_model_data$USUBJID),
        n_pseudo_rows = nrow(block_model_data),
        weighted_outcome_events = sum(block_model_data$posterior_weight * block_model_data$DPBN, na.rm = TRUE),
        weighted_mediator_events = NA_real_,
        converged = isTRUE(block_outcome_fit$converged),
        warnings = paste(attr(block_outcome_fit, "fit_warnings"), collapse = " | ")
      )

      for (contrast_i in seq_len(nrow(window_hypotheses))) {
        hypothesis <- window_hypotheses[contrast_i, ]
        block_decomp <- estimate_joint_mediator_block(
          mediator_fits = block_mediator_fits,
          outcome_fit = block_outcome_fit,
          model_data = block_model_data,
          exposure_var = exposure_var,
          mediator_vars = mediator_vars,
          reference = hypothesis$reference_history,
          comparator = hypothesis$comparator_histories,
          mediator_draws = block_mediator_draws,
          outcome_draws = block_outcome_draws
        )

        joint_block_results[[length(joint_block_results) + 1]] <- block_decomp %>%
          mutate(
            adjustment_set = adjustment_set,
            adjustment_label = adjustment_label,
            tier = hypothesis$tier,
            hypothesis_id = hypothesis$hypothesis_id,
            window = hypothesis$window,
            scenario = hypothesis$scenario,
            reference = hypothesis$reference_history,
            comparator = hypothesis$comparator_histories,
            contrast = paste(comparator, "vs", reference),
            contrast_direction = "Comparator minus reference",
            expected_direction = hypothesis$expected_total_effect_direction,
            mediator_block_id = block_id,
            mediator_block = block_label,
            mediators_in_block = paste(mediator_labels, collapse = "; "),
            n_mediators = length(mediator_vars),
            n_patterns = 2^length(mediator_vars),
            n_subjects = n_distinct(block_model_data$USUBJID),
            n_pseudo_rows = nrow(block_model_data),
            interval_method = if (current_block_sim > 1) {
              paste0("cluster-robust coefficient simulation, ", current_block_sim, " sequential joint-block draws")
            } else {
              "Point estimate only; interval simulation not run for this block"
            },
            interpretation_flag = pmap_chr(list(rd, ci_low, ci_high), effect_flag)
          ) %>%
          select(
            adjustment_set, adjustment_label, tier, hypothesis_id, window, scenario,
            reference, comparator, contrast, contrast_direction, expected_direction,
            mediator_block_id, mediator_block, mediators_in_block, n_mediators, n_patterns,
            effect, n_subjects, n_pseudo_rows, rd, ci_low, ci_high, rd_ci,
            within_variance,
            percent_mediated, percent_mediated_display,
            interpretation_flag, interval_method
          )
      }
    }
  }
}

total_results <- bind_rows(total_results)
mediator_results <- bind_rows(mediator_results)
joint_block_results <- bind_rows(joint_block_results)
model_diagnostics <- bind_rows(model_diagnostics)

fixed_context_comparison <- mediator_results %>%
  filter(effect == "Indirect effect", mediator_scope == "Base mediator") %>%
  select(
    adjustment_set, hypothesis_id, window, contrast, mediator, rd, ci_low, ci_high, rd_ci,
    percent_mediated_display, n_subjects
  ) %>%
  pivot_wider(
    id_cols = c(hypothesis_id, window, contrast, mediator),
    names_from = adjustment_set,
    values_from = c(rd, ci_low, ci_high, rd_ci, percent_mediated_display, n_subjects),
    names_glue = "{adjustment_set}_{.value}"
  ) %>%
  mutate(
    absolute_change = fixed_context_rd - base_rd,
    relative_attenuation = if_else(
      is.finite(base_rd) &
        is.finite(fixed_context_rd) &
        abs(base_rd) >= 0.002 &
        sign(base_rd) == sign(fixed_context_rd),
      1 - (fixed_context_rd / base_rd),
      NA_real_
    ),
    base_indirect_rd_ci = base_rd_ci,
    fixed_context_indirect_rd_ci = fixed_context_rd_ci,
    absolute_change_display = if_else(is.finite(absolute_change), sprintf("%.3f", absolute_change), ""),
    relative_attenuation_display = if_else(is.finite(relative_attenuation), sprintf("%.1f%%", 100 * relative_attenuation), ""),
    interpretation = case_when(
      !is.finite(base_rd) | !is.finite(fixed_context_rd) ~ "Not estimable",
      abs(base_rd) < 0.002 ~ "Base indirect effect near null; attenuation not summarised",
      sign(base_rd) != sign(fixed_context_rd) ~ "Indirect-effect direction changed after fixed-context adjustment",
      abs(fixed_context_rd) < abs(base_rd) ~ "Attenuated after fixed-context adjustment",
      abs(fixed_context_rd) > abs(base_rd) ~ "Larger after fixed-context adjustment",
      TRUE ~ "Unchanged"
    )
  ) %>%
  select(
    hypothesis_id, window, contrast, mediator,
    base_indirect_rd_ci, fixed_context_indirect_rd_ci,
    absolute_change_display, relative_attenuation_display,
    interpretation,
    base_n_subjects, fixed_context_n_subjects
  )

run_manifest <- tibble(
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  script = script_path,
  workspace_root = workspace_root,
  branch = branch_spec$branch,
  subject_input_path = subject_input_path,
  adjustment_config_path = adjustment_config_path,
  n_simulation_draws = n_sim,
  seed = seed,
  n_block_simulation_draws = block_sim,
  n_expanded_block_simulation_draws = expanded_block_sim,
  transition_subjects = n_distinct(pseudo_data$USUBJID),
  transition_pseudo_rows = nrow(pseudo_data),
  posterior_weight_floor = weight_floor,
  posterior_weights_normalised_within_subject = TRUE,
  exposure_solution = branch_spec$exposure_solution,
  n_total_effect_contrasts = nrow(branch_total_effect_hypotheses),
  interval_method = "Cluster-robust coefficient simulation; bootstrap remains a later manuscript-facing uncertainty check."
)

total_results_path <- file.path(listing_dir, paste0(output_prefix, "_total_effect_results.csv"))
single_results_path <- file.path(listing_dir, paste0(output_prefix, "_single_mediator_results.csv"))
joint_results_path <- file.path(listing_dir, paste0(output_prefix, "_joint_mediator_results.csv"))
fixed_context_path <- file.path(listing_dir, paste0(output_prefix, "_fixed_context_indirect_comparison.csv"))
diagnostics_path <- file.path(listing_dir, paste0(output_prefix, "_model_diagnostics.csv"))
manifest_path <- file.path(listing_dir, paste0(output_prefix, "_mediation_run_manifest.csv"))

write_csv(total_results, total_results_path)
write_csv(mediator_results, single_results_path)
write_csv(joint_block_results, joint_results_path)
write_csv(fixed_context_comparison, fixed_context_path)
write_csv(model_diagnostics, diagnostics_path)
write_csv(run_manifest, manifest_path)

message("Wrote total-effect results: ", total_results_path)
message("Wrote single-mediator results: ", single_results_path)
message("Wrote joint-mediator block results: ", joint_results_path)
message("Wrote fixed-context comparison: ", fixed_context_path)
message("Wrote diagnostics: ", diagnostics_path)
message("Wrote run manifest: ", manifest_path)
