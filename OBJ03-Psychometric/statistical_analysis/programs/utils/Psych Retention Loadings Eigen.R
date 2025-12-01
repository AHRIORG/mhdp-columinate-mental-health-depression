# Psychometric Bootstrap EFA – Retention, Loadings & Eigenvalue Helpers

# Dependencies: dplyr, tidyr, ggplot2, scales

# ------------------------------------------------------------------------------
# 1. Factor Retention Analysis
# ------------------------------------------------------------------------------

analyze_factor_retention <- function(df, custom_subtitle = "") {
  clean_df <- df %>% dplyr::filter(Is_Suitable == TRUE)
  n_total  <- nrow(clean_df)
  if (n_total == 0) return(NULL)
  
  long_df <- clean_df %>%
    dplyr::select(Iteration, k_factor_f, p_factor_f) %>%
    tidyr::pivot_longer(
      cols      = c(k_factor_f, p_factor_f),
      names_to  = "Method_Code",
      values_to = "Factors_Suggested"
    ) %>%
    dplyr::filter(!is.na(Factors_Suggested)) %>%
    dplyr::mutate(
      Method = dplyr::case_when(
        Method_Code == "k_factor_f" ~ "Kaiser Criterion",
        Method_Code == "p_factor_f" ~ "Parallel Analysis"
      )
    )
  
  stats_df <- long_df %>%
    dplyr::group_by(Method, Factors_Suggested) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(Method) %>%
    dplyr::mutate(
      Total_Valid_For_Method = sum(Count),
      Probability = Count / Total_Valid_For_Method,
      SE = sqrt((Probability * (1 - Probability)) / Total_Valid_For_Method),
      Lower_CI = pmax(0, Probability - 1.96 * SE),
      Upper_CI = pmin(1, Probability + 1.96 * SE)
    ) %>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(
    stats_df,
    ggplot2::aes(
      x    = factor(Factors_Suggested),
      y    = Probability,
      fill = Method
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      width    = 0.7,
      alpha    = 0.8
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Lower_CI, ymax = Upper_CI),
      position = ggplot2::position_dodge(width = 0.8),
      width    = 0.25,
      alpha    = 0.7
    ) +
    ggplot2::geom_point(
      position   = ggplot2::position_dodge(width = 0.8),
      size       = 2,
      color      = "black",
      show.legend = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels  = scales::percent_format(),
      limits  = c(0, 1.1),
      sec.axis = ggplot2::sec_axis(~ . * n_total, name = "Count")
    ) +
    ggplot2::labs(
      title    = NULL,
      subtitle = custom_subtitle,
      caption  = NULL,
      x        = "Number of Factors",
      y        = "Probability"
    ) +
    ggplot2::theme_minimal(base_family = "Arial", base_size = 14) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(family = "Arial", size = 12)
    )
  
  list(stats = stats_df, plot = p)
}

# ------------------------------------------------------------------------------
# 2. Loadings Extraction Helper
# ------------------------------------------------------------------------------

#' Extract long-format loadings with CI and default plot
#'
#' @param results_collection The list returned by run_psychometric_pipeline().
#' @param set_name Character, name of the item set (e.g., "SSQ-14 (Original)").
#' @param method "Kaiser" or "Parallel".
#' @param factor_filter Optional integer or vector of factors to keep.
#' @param ci_level Confidence level for CI (default 0.95).
#'
#' @return A list with elements: data (long df), summary (CI table), plot (ggplot).
extract_loadings_long <- function(results_collection,
                                  set_name,
                                  method = c("Kaiser", "Parallel"),
                                  factor_filter = NULL,
                                  ci_level = 0.95) {
  method <- match.arg(method)
  if (!set_name %in% names(results_collection)) {
    stop("set_name not found in results_collection: ", set_name)
  }
  raw_df <- results_collection[[set_name]]$raw_data
  if (is.null(raw_df)) {
    stop("No raw_data found for set_name: ", set_name)
  }
  
  col_name <- if (method == "Kaiser") "k_loadings" else "p_loadings"
  if (!col_name %in% names(raw_df)) {
    stop("Column ", col_name, " not found in raw_data.")
  }
  
  tmp <- raw_df[, c("Iteration", col_name)]
  names(tmp)[2] <- "loadings"
  long_df <- tidyr::unnest(tmp, loadings)
  
  if (!is.null(factor_filter)) {
    long_df <- long_df %>% dplyr::filter(Factor %in% factor_filter)
  }
  
  long_df <- long_df %>%
    dplyr::mutate(Method = method, Set = set_name)
  
  alpha <- 1 - ci_level
  summary_df <- long_df %>%
    dplyr::group_by(Item, Factor) %>%
    dplyr::summarise(
      n            = dplyr::n(),
      mean_loading = mean(Loading, na.rm = TRUE),
      sd_loading   = sd(Loading, na.rm = TRUE),
      se_loading   = sd_loading / sqrt(n),
      ci_lower     = mean_loading + stats::qnorm(alpha / 2) * se_loading,
      ci_upper     = mean_loading + stats::qnorm(1 - alpha / 2) * se_loading,
      .groups      = "drop"
    ) %>%
    dplyr::mutate(Method = method, Set = set_name)
  
  p <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = Loading, fill = factor(Factor))
  ) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    ggplot2::facet_wrap(~ Item, scales = "free_y") +
    ggplot2::theme_minimal(base_family = "Arial", base_size = 14) +
    ggplot2::labs(
      title = paste0("Bootstrap factor loadings (", method, "): ", set_name),
      x     = "Loading",
      y     = "Count",
      fill  = "Factor"
    )
  
  list(data = long_df, summary = summary_df, plot = p)
}

# ------------------------------------------------------------------------------
# 3. Eigenvalue Extraction Helper
# ------------------------------------------------------------------------------

#' Extract eigenvalues (actual vs simulated) with CI and scree plot
#'
#' @param results_collection The list returned by run_psychometric_pipeline().
#' @param set_name Character, name of the item set.
#' @param ci_level Confidence level for CI (default 0.95).
#'
#' @return A list with elements: data (long df), summary (CI table), plot (ggplot).
extract_eigen_scree <- function(results_collection,
                                set_name,
                                ci_level = 0.95) {
  if (!set_name %in% names(results_collection)) {
    stop("set_name not found in results_collection: ", set_name)
  }
  raw_df <- results_collection[[set_name]]$raw_data
  if (is.null(raw_df)) {
    stop("No raw_data found for set_name: ", set_name)
  }
  
  tmp_a <- raw_df[, c("Iteration", "eig_actual")]
  names(tmp_a)[2] <- "eig"
  actual_long <- tidyr::unnest(tmp_a, eig)
  if (nrow(actual_long) > 0) {
    actual_long <- actual_long %>%
      dplyr::group_by(Iteration) %>%
      dplyr::mutate(Component = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::rename(Value = eig) %>%
      dplyr::mutate(Type = "Actual", Set = set_name)
  }
  
  tmp_s <- raw_df[, c("Iteration", "eig_simulated")]
  names(tmp_s)[2] <- "eig"
  sim_long <- tidyr::unnest(tmp_s, eig)
  if (nrow(sim_long) > 0) {
    sim_long <- sim_long %>%
      dplyr::group_by(Iteration) %>%
      dplyr::mutate(Component = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::rename(Value = eig) %>%
      dplyr::mutate(Type = "Simulated", Set = set_name)
  }
  
  long_df <- dplyr::bind_rows(actual_long, sim_long)
  if (nrow(long_df) == 0) {
    return(list(data = long_df, summary = NULL, plot = NULL))
  }
  
  alpha <- 1 - ci_level
  summary_df <- long_df %>%
    dplyr::group_by(Type, Component) %>%
    dplyr::summarise(
      n          = dplyr::n(),
      mean_value = mean(Value, na.rm = TRUE),
      sd_value   = sd(Value, na.rm = TRUE),
      se_value   = sd_value / sqrt(n),
      ci_lower   = mean_value + stats::qnorm(alpha / 2) * se_value,
      ci_upper   = mean_value + stats::qnorm(1 - alpha / 2) * se_value,
      .groups    = "drop"
    )
  
  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = Component, y = mean_value, color = Type, fill = Type)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      alpha = 0.2,
      color = NA
    ) +
    ggplot2::theme_minimal(base_family = "Arial", base_size = 14) +
    ggplot2::labs(
      title = paste0("Bootstrap Scree with ", round(ci_level * 100), "% CI: ", set_name),
      x     = "Component",
      y     = "Eigenvalue"
    )
  
  list(data = long_df, summary = summary_df, plot = p)
}
