#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(paste0("^", flag, "="), "", hit[[1]])
}

input_prefix <- arg_value("--input-prefix", "strict_jmaes_final_imp_")
output_prefix <- arg_value("--output-prefix", "strict_jmaes_final")
m <- as.integer(arg_value("--m", "20"))

script_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][[1]])
script_arg <- gsub("~\\+~", " ", script_arg)
script_path <- normalizePath(script_arg, mustWork = TRUE)
workspace_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
listing_dir <- file.path(workspace_root, "output", "listings")

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

pool_interval <- function(q, u) {
  keep <- is.finite(q)
  q <- q[keep]
  u <- u[keep]
  if (!length(q)) {
    return(list(
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      within_variance = NA_real_,
      between_variance = NA_real_,
      total_variance = NA_real_,
      df = NA_real_,
      m_used = 0L
    ))
  }

  u <- ifelse(is.finite(u), u, NA_real_)
  estimate <- mean(q)
  ubar <- if (all(is.na(u))) NA_real_ else mean(u, na.rm = TRUE)
  b <- if (length(q) > 1) stats::var(q) else 0
  total_variance <- if (is.finite(ubar)) ubar + (1 + 1 / length(q)) * b else NA_real_

  if (!is.finite(total_variance) || total_variance < 0) {
    return(list(
      estimate = estimate,
      ci_low = NA_real_,
      ci_high = NA_real_,
      within_variance = ubar,
      between_variance = b,
      total_variance = total_variance,
      df = NA_real_,
      m_used = length(q)
    ))
  }

  if (isTRUE(all.equal(b, 0)) || !is.finite(ubar) || isTRUE(all.equal(ubar, 0))) {
    df <- Inf
  } else {
    r <- ((1 + 1 / length(q)) * b) / ubar
    df <- if (is.finite(r) && r > 0) (length(q) - 1) * (1 + 1 / r)^2 else Inf
  }

  crit <- if (is.finite(df)) stats::qt(0.975, df = df) else stats::qnorm(0.975)
  se <- sqrt(total_variance)

  list(
    estimate = estimate,
    ci_low = estimate - crit * se,
    ci_high = estimate + crit * se,
    within_variance = ubar,
    between_variance = b,
    total_variance = total_variance,
    df = df,
    m_used = length(q)
  )
}

pool_results <- function(data, group_cols, carry_cols = character(), mean_cols = character()) {
  grouped <- group_split(group_by(data, across(all_of(group_cols))))

  bind_rows(lapply(grouped, function(df) {
    pooled <- pool_interval(df$rd, df$within_variance)
    out <- df[1, unique(c(group_cols, carry_cols)), drop = FALSE]
    for (col in intersect(mean_cols, names(df))) {
      values <- df[[col]]
      out[[col]] <- if (any(is.finite(values))) mean(values[is.finite(values)]) else NA_real_
    }
    out$rd <- pooled$estimate
    out$ci_low <- pooled$ci_low
    out$ci_high <- pooled$ci_high
    out$rd_ci <- format_rd_ci(pooled$estimate, pooled$ci_low, pooled$ci_high)
    out$within_variance <- pooled$within_variance
    out$between_variance <- pooled$between_variance
    out$total_variance <- pooled$total_variance
    out$rubin_df <- pooled$df
    out$m_imputations <- pooled$m_used
    out
  }))
}

derive_percent_mediated_single <- function(data) {
  totals <- data %>%
    filter(effect == "Total effect") %>%
    transmute(
      adjustment_set, adjustment_label, tier, hypothesis_id, window, scenario,
      reference, comparator, contrast, contrast_direction, expected_direction,
      mediator_scope, mediator_variable, mediator,
      total_rd = rd
    )

  data %>%
    left_join(
      totals,
      by = c(
        "adjustment_set", "adjustment_label", "tier", "hypothesis_id", "window", "scenario",
        "reference", "comparator", "contrast", "contrast_direction", "expected_direction",
        "mediator_scope", "mediator_variable", "mediator"
      )
    ) %>%
    mutate(
      percent_mediated = if_else(
        effect == "Indirect effect" &
          is.finite(.data$total_rd) &
          abs(.data$total_rd) >= 0.01 &
          sign(.data$rd) == sign(.data$total_rd),
        .data$rd / .data$total_rd,
        NA_real_
      ),
      percent_mediated_display = map_chr(.data$percent_mediated, format_percent),
      interpretation_flag = pmap_chr(list(.data$rd, .data$ci_low, .data$ci_high), effect_flag)
    ) %>%
    select(-total_rd)
}

derive_percent_mediated_joint <- function(data) {
  totals <- data %>%
    filter(effect == "Total effect") %>%
    transmute(
      adjustment_set, adjustment_label, tier, hypothesis_id, window, scenario,
      reference, comparator, contrast, contrast_direction, expected_direction,
      mediator_block_id, mediator_block, mediators_in_block, n_mediators, n_patterns,
      total_rd = rd
    )

  data %>%
    left_join(
      totals,
      by = c(
        "adjustment_set", "adjustment_label", "tier", "hypothesis_id", "window", "scenario",
        "reference", "comparator", "contrast", "contrast_direction", "expected_direction",
        "mediator_block_id", "mediator_block", "mediators_in_block", "n_mediators", "n_patterns"
      )
    ) %>%
    mutate(
      percent_mediated = if_else(
        effect == "Total indirect effect" &
          is.finite(.data$total_rd) &
          abs(.data$total_rd) >= 0.01 &
          sign(.data$rd) == sign(.data$total_rd),
        .data$rd / .data$total_rd,
        NA_real_
      ),
      percent_mediated_display = map_chr(.data$percent_mediated, format_percent),
      interpretation_flag = pmap_chr(list(.data$rd, .data$ci_low, .data$ci_high), effect_flag)
    ) %>%
    select(-total_rd)
}

result_paths <- function(suffix) {
  file.path(listing_dir, sprintf("%s%02d_%s.csv", input_prefix, seq_len(m), suffix))
}

read_bound <- function(paths) {
  keep <- file.exists(paths)
  if (!all(keep)) {
    stop("Missing expected imputation result files:\n", paste(paths[!keep], collapse = "\n"))
  }
  map2_dfr(paths, seq_along(paths), function(path, idx) {
    read_csv(path, show_col_types = FALSE) %>%
      mutate(imputation_id = idx)
  })
}

total_raw <- read_bound(result_paths("total_effect_results"))
single_raw <- read_bound(result_paths("single_mediator_results"))
joint_raw <- read_bound(result_paths("joint_mediator_results"))
diagnostics_raw <- read_bound(result_paths("model_diagnostics"))
manifest_raw <- read_bound(result_paths("mediation_run_manifest"))

total_group_cols <- c(
  "adjustment_set", "adjustment_label", "tier", "hypothesis_id", "window", "scenario",
  "reference", "comparator", "contrast", "contrast_direction", "expected_direction",
  "n_subjects", "n_pseudo_rows"
)
single_group_cols <- c(
  "adjustment_set", "adjustment_label", "tier", "hypothesis_id", "window", "scenario",
  "reference", "comparator", "contrast", "contrast_direction", "expected_direction",
  "mediator_scope", "mediator_variable", "mediator", "effect",
  "n_subjects", "n_pseudo_rows"
)
joint_group_cols <- c(
  "adjustment_set", "adjustment_label", "tier", "hypothesis_id", "window", "scenario",
  "reference", "comparator", "contrast", "contrast_direction", "expected_direction",
  "mediator_block_id", "mediator_block", "mediators_in_block", "n_mediators", "n_patterns",
  "effect", "n_subjects", "n_pseudo_rows"
)

total_results <- pool_results(
  total_raw,
  total_group_cols,
  carry_cols = c("interval_method"),
  mean_cols = c("p_ref", "p_comp")
) %>%
  mutate(
    proceed_to_decomposition = if_else(
      is.finite(.data$ci_low) & is.finite(.data$ci_high) & (.data$ci_low > 0 | .data$ci_high < 0),
      "Prespecified decomposition; total-effect interval excludes 0",
      "Prespecified decomposition; interpret descriptively"
    ),
    interpretation_flag = pmap_chr(list(.data$rd, .data$ci_low, .data$ci_high), effect_flag),
    interval_method = paste0("Rubin pooled over ", .data$m_imputations, " imputations; ", .data$interval_method)
  ) %>%
  select(all_of(total_raw %>% select(-rd, -ci_low, -ci_high, -rd_ci, -within_variance) %>% names() %>% setdiff(c("imputation_id", "interpretation_flag"))),
         rd, ci_low, ci_high, rd_ci, within_variance, between_variance, total_variance, rubin_df, m_imputations, interpretation_flag)

single_results <- pool_results(single_raw, single_group_cols, carry_cols = c("interval_method")) %>%
  mutate(
    interval_method = paste0("Rubin pooled over ", .data$m_imputations, " imputations; ", .data$interval_method)
  ) %>%
  derive_percent_mediated_single()

joint_results <- pool_results(joint_raw, joint_group_cols, carry_cols = c("interval_method")) %>%
  mutate(
    interval_method = paste0("Rubin pooled over ", .data$m_imputations, " imputations; ", .data$interval_method)
  ) %>%
  derive_percent_mediated_joint()

fixed_context_comparison <- single_results %>%
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

diagnostics_summary <- diagnostics_raw %>%
  group_by(.data$adjustment_set, .data$window, .data$mediator, .data$model_type, .data$formula) %>%
  summarise(
    n_subjects = round(mean(.data$n_subjects, na.rm = TRUE)),
    n_pseudo_rows = round(mean(.data$n_pseudo_rows, na.rm = TRUE)),
    weighted_outcome_events = mean(.data$weighted_outcome_events, na.rm = TRUE),
    weighted_mediator_events = mean(.data$weighted_mediator_events, na.rm = TRUE),
    converged_all = all(.data$converged),
    n_warning_fits = sum(nzchar(.data$warnings)),
    warnings = paste(unique(.data$warnings[nzchar(.data$warnings)]), collapse = " | "),
    m_imputations = n(),
    .groups = "drop"
  )

run_manifest <- tibble(
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  script = script_path,
  workspace_root = workspace_root,
  branch = "strict",
  subject_input_path = "Rubin-pooled across imputed strict subject-level datasets",
  adjustment_config_path = file.path(workspace_root, "config", "adjustment_sets_final_inferential.csv"),
  n_simulation_draws = suppressWarnings(as.integer(manifest_raw$n_simulation_draws[[1]])),
  seed = NA_integer_,
  n_block_simulation_draws = suppressWarnings(as.integer(manifest_raw$n_block_simulation_draws[[1]])),
  n_expanded_block_simulation_draws = suppressWarnings(as.integer(manifest_raw$n_expanded_block_simulation_draws[[1]])),
  transition_subjects = round(mean(manifest_raw$transition_subjects, na.rm = TRUE)),
  transition_pseudo_rows = round(mean(manifest_raw$transition_pseudo_rows, na.rm = TRUE)),
  posterior_weight_floor = mean(manifest_raw$posterior_weight_floor, na.rm = TRUE),
  posterior_weights_normalised_within_subject = all(manifest_raw$posterior_weights_normalised_within_subject),
  exposure_solution = manifest_raw$exposure_solution[[1]],
  n_total_effect_contrasts = round(mean(manifest_raw$n_total_effect_contrasts, na.rm = TRUE)),
  interval_method = paste0("Rubin pooled over ", m, " imputations; within-imputation intervals from cluster-robust coefficient simulation."),
  m_imputations = m
)

write_csv(total_results, file.path(listing_dir, paste0(output_prefix, "_total_effect_results.csv")))
write_csv(single_results, file.path(listing_dir, paste0(output_prefix, "_single_mediator_results.csv")))
write_csv(joint_results, file.path(listing_dir, paste0(output_prefix, "_joint_mediator_results.csv")))
write_csv(fixed_context_comparison, file.path(listing_dir, paste0(output_prefix, "_fixed_context_indirect_comparison.csv")))
write_csv(diagnostics_summary, file.path(listing_dir, paste0(output_prefix, "_model_diagnostics.csv")))
write_csv(run_manifest, file.path(listing_dir, paste0(output_prefix, "_mediation_run_manifest.csv")))

message("Wrote pooled MI results with prefix ", output_prefix)
