############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 04_summary.R - Summaries & gt tables
############################################################

# Helper: format estimate with 95% CI ---------------------------------
# Vectorised: est, lower, upper can be vectors of same length
format_est_ci <- function(est, lower, upper, digits = 3) {
  n <- length(est)
  out <- rep(NA_character_, n)
  valid <- !is.na(est) & !is.na(lower) & !is.na(upper)
  if (any(valid)) {
    out[valid] <- paste0(
      round(est[valid],   digits),
      " (", round(lower[valid], digits), ", ", round(upper[valid], digits), ")"
    )
  }
  out
}

# 1. Model fit summary gt ----------------------------------------------
summarise_fit_gt <- function(res_fit) {
  if (is.null(res_fit) || nrow(res_fit) == 0) return(NULL)
  
  summary_fit <- res_fit %>%
    summarise(
      N_Valid         = n(),
      Mean_CFI        = mean(cfi,   na.rm = TRUE),
      Lower_CFI       = quantile(cfi,   0.025, na.rm = TRUE),
      Upper_CFI       = quantile(cfi,   0.975, na.rm = TRUE),
      Mean_RMSEA      = mean(rmsea, na.rm = TRUE),
      Lower_RMSEA     = quantile(rmsea, 0.025, na.rm = TRUE),
      Upper_RMSEA     = quantile(rmsea, 0.975, na.rm = TRUE),
      Prop_Good_CFI   = mean(cfi   > 0.90, na.rm = TRUE),
      Prop_Good_RMSEA = mean(rmsea < 0.08, na.rm = TRUE)
    ) %>%
    mutate(
      `CFI (95% CI)`   = format_est_ci(Mean_CFI,   Lower_CFI,   Upper_CFI),
      `RMSEA (95% CI)` = format_est_ci(Mean_RMSEA, Lower_RMSEA, Upper_RMSEA)
    )
  
  tab_fit <- summary_fit %>%
    select(
      N_Valid,
      `CFI (95% CI)`,
      Prop_Good_CFI,
      `RMSEA (95% CI)`,
      Prop_Good_RMSEA
    )
  
  summary_gt <- gt(tab_fit) %>%
    tab_header(
      title = "Global Fit Indices Across Monte Carlo Samples"
    ) %>%
    cols_label(
      N_Valid         = "Number of Valid Runs",
      `CFI (95% CI)`   = "CFI (95% CI)",
      Prop_Good_CFI   = "Proportion CFI > 0.90",
      `RMSEA (95% CI)` = "RMSEA (95% CI)",
      Prop_Good_RMSEA = "Proportion RMSEA < 0.08"
    ) %>%
    fmt_number(
      columns  = c(Prop_Good_CFI, Prop_Good_RMSEA),
      decimals = 3
    ) %>%
    tab_source_note(
      source_note = paste(
        "Interpretation guide: CFI ≥ 0.90 indicates acceptable fit;",
        "RMSEA ≤ 0.08 indicates acceptable fit."
      )
    )
  
  list(
    summary_source = summary_fit,
    summary_gt     = summary_gt
  )
}

# 2. Factor loadings summary gt ----------------------------------------
summarise_loadings_gt <- function(res_loadings) {
  if (is.null(res_loadings) || nrow(res_loadings) == 0) return(NULL)
  
  summary_loadings <- res_loadings %>%
    group_by(lhs, rhs) %>%
    summarise(
      Mean_Std_Est = mean(est.std, na.rm = TRUE),
      Lower_95     = quantile(est.std, 0.025, na.rm = TRUE),
      Upper_95     = quantile(est.std, 0.975, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    mutate(
      `Loading (95% CI)` = format_est_ci(Mean_Std_Est, Lower_95, Upper_95)
    )
  
  tab_load <- summary_loadings %>%
    select(
      Factor = lhs,
      Item   = rhs,
      `Loading (95% CI)`
    )
  
  summary_gt <- gt(tab_load) %>%
    tab_header(
      title = "Standardized Factor Loadings Across Monte Carlo Samples"
    ) %>%
    cols_label(
      Factor             = "Factor",
      Item               = "Item",
      `Loading (95% CI)` = "Standardized Loading (95% CI)"
    ) %>%
    tab_source_note(
      source_note = paste(
        "Interpretation guide: loadings ≥ 0.40 are often considered practically important,",
        "with higher values indicating stronger association between item and factor."
      )
    )
  
  list(
    summary_source = summary_loadings,
    summary_gt     = summary_gt
  )
}

# 3. Measurement invariance summary gt ---------------------------------
summarise_invariance_gt <- function(res_invariance) {
  if (is.null(res_invariance) || nrow(res_invariance) == 0) return(NULL)
  
  summary_inv <- res_invariance %>%
    group_by(Group) %>%
    summarise(
      N_Runs           = n(),
      Mean_Delta_CFI   = round(mean(Delta_CFI,   na.rm = TRUE), 4),
      Max_Delta_CFI    = max(Delta_CFI,   na.rm = TRUE),
      Prop_Inv_CFI     = mean(abs(Delta_CFI)   <= 0.010, na.rm = TRUE),
      Mean_Delta_RMSEA = mean(Delta_RMSEA, na.rm = TRUE),
      Prop_Inv_RMSEA   = mean(abs(Delta_RMSEA) <= 0.015, na.rm = TRUE),
      .groups          = "drop"
    )
  
  summary_gt <- gt(summary_inv) %>%
    tab_header(
      title = "Measurement Invariance: Scalar vs Configural"
    ) %>%
    cols_label(
      Group            = "Grouping Variable",
      N_Runs           = "Number of Runs",
      Mean_Delta_CFI   = "Mean ΔCFI (Scalar - Configural)",
      Max_Delta_CFI    = "Max ΔCFI",
      Prop_Inv_CFI     = "Proportion |ΔCFI| ≤ 0.010",
      Mean_Delta_RMSEA = "Mean ΔRMSEA (Scalar - Configural)",
      Prop_Inv_RMSEA   = "Proportion |ΔRMSEA| ≤ 0.015"
    ) %>%
    fmt_number(
      columns  = c(
        Mean_Delta_CFI, Max_Delta_CFI, Prop_Inv_CFI,
        Mean_Delta_RMSEA, Prop_Inv_RMSEA
      ),
      decimals = 3
    ) %>%
    tab_source_note(
      source_note = paste(
        "Interpretation guide: |ΔCFI| ≤ 0.010 and |ΔRMSEA| ≤ 0.015",
        "are commonly used thresholds for supporting scalar invariance."
      )
    )
  
  list(
    summary_source = summary_inv,
    summary_gt     = summary_gt
  )
}

# 4. ROC summary gt (all groups) ---------------------------------------
summarise_roc_gt <- function(res_roc) {
  if (is.null(res_roc) || nrow(res_roc) == 0) return(NULL)
  
  summary_roc <- res_roc %>%
    group_by(Group) %>%
    summarise(
      N_Runs          = n(),
      Mean_AUC        = mean(AUC, na.rm = TRUE),
      Lower_AUC       = quantile(AUC, 0.025, na.rm = TRUE),
      Upper_AUC       = quantile(AUC, 0.975, na.rm = TRUE),
      Mean_Sens       = mean(Sensitivity, na.rm = TRUE),
      Lower_Sens      = quantile(Sensitivity, 0.025, na.rm = TRUE),
      Upper_Sens      = quantile(Sensitivity, 0.975, na.rm = TRUE),
      Mean_Spec       = mean(Specificity, na.rm = TRUE),
      Lower_Spec      = quantile(Specificity, 0.025, na.rm = TRUE),
      Upper_Spec      = quantile(Specificity, 0.975, na.rm = TRUE),
      Mean_Kappa      = mean(Kappa, na.rm = TRUE),
      Lower_Kappa     = quantile(Kappa, 0.025, na.rm = TRUE),
      Upper_Kappa     = quantile(Kappa, 0.975, na.rm = TRUE),
      Mean_Cutoff     = if ("Optimal_Cutoff" %in% colnames(cur_data_all())) mean(Optimal_Cutoff, na.rm = TRUE) else NA_real_,
      Lower_Cutoff    = if ("Optimal_Cutoff" %in% colnames(cur_data_all())) quantile(Optimal_Cutoff, 0.025, na.rm = TRUE) else NA_real_,
      Upper_Cutoff    = if ("Optimal_Cutoff" %in% colnames(cur_data_all())) quantile(Optimal_Cutoff, 0.975, na.rm = TRUE) else NA_real_,
      Mean_Prev_PHQ   = if ("Prevalence_PHQ" %in% colnames(cur_data_all())) mean(Prevalence_PHQ, na.rm = TRUE) else NA_real_,
      Lower_Prev_PHQ  = if ("Prevalence_PHQ" %in% colnames(cur_data_all())) quantile(Prevalence_PHQ, 0.025, na.rm = TRUE) else NA_real_,
      Upper_Prev_PHQ  = if ("Prevalence_PHQ" %in% colnames(cur_data_all())) quantile(Prevalence_PHQ, 0.975, na.rm = TRUE) else NA_real_,
      Mean_Prev_SSQ   = if ("Prevalence_SSQ" %in% colnames(cur_data_all())) mean(Prevalence_SSQ, na.rm = TRUE) else NA_real_,
      Lower_Prev_SSQ  = if ("Prevalence_SSQ" %in% colnames(cur_data_all())) quantile(Prevalence_SSQ, 0.025, na.rm = TRUE) else NA_real_,
      Upper_Prev_SSQ  = if ("Prevalence_SSQ" %in% colnames(cur_data_all())) quantile(Prevalence_SSQ, 0.975, na.rm = TRUE) else NA_real_,
      .groups         = "drop"
    ) %>%
    mutate(
      `AUC (95% CI)`        = format_est_ci(Mean_AUC,      Lower_AUC,      Upper_AUC),
      `Sens (95% CI)`       = format_est_ci(Mean_Sens,     Lower_Sens,     Upper_Sens),
      `Spec (95% CI)`       = format_est_ci(Mean_Spec,     Lower_Spec,     Upper_Spec),
      `Kappa (95% CI)`      = format_est_ci(Mean_Kappa,    Lower_Kappa,    Upper_Kappa),
      `Cutoff (95% CI)`     = format_est_ci(Mean_Cutoff,   Lower_Cutoff,   Upper_Cutoff),
      `Prev PHQ (95% CI)`   = format_est_ci(Mean_Prev_PHQ * 100, Lower_Prev_PHQ * 100, Upper_Prev_PHQ * 100, digits = 1),
      `Prev SSQ (95% CI)`   = format_est_ci(Mean_Prev_SSQ * 100, Lower_Prev_SSQ * 100, Upper_Prev_SSQ * 100, digits = 1)
    )
  
  tab_roc <- summary_roc %>%
    select(
      Group,
      `AUC (95% CI)`,
      `Sens (95% CI)`,
      `Spec (95% CI)`,
      `Kappa (95% CI)`,
      `Cutoff (95% CI)`,
      `Prev PHQ (95% CI)`,
      `Prev SSQ (95% CI)`,
      N_Runs
    )
  
  summary_gt <- gt(tab_roc) %>%
    tab_header(
      title = "Diagnostic Accuracy: ROC-based Metrics (PHQ-9 ≥ 10 as Benchmark)"
    ) %>%
    cols_label(
      Group               = "Group",
      `AUC (95% CI)`      = "AUC (95% CI)",
      `Sens (95% CI)`     = "Sensitivity (95% CI)",
      `Spec (95% CI)`     = "Specificity (95% CI)",
      `Kappa (95% CI)`    = "Cohen's Kappa (95% CI)",
      `Cutoff (95% CI)`   = "Optimal SSQ-10 Cutoff (95% CI)",
      `Prev PHQ (95% CI)` = "Prevalence PHQ-09 ≥ 10 (95% CI)",
      `Prev SSQ (95% CI)` = "Prevalence SSQ-10 above cutoff (95% CI)",
      N_Runs              = "Number of Runs"
    ) %>%
    fmt_number(
      columns  = c(N_Runs),
      decimals = 0
    ) %>%
    tab_source_note(
      source_note = paste(
        "Interpretation guide: AUC 0.70–0.80 acceptable, 0.80–0.90 excellent, >0.90 outstanding.",
        "Kappa > 0.60 indicates substantial agreement; >0.80 indicates almost perfect agreement.",
        "Prevalence is calculated separately using PHQ-09 (gold standard) and SSQ-10 at the chosen cutoff.")
    )
  
  list(
    summary_source = summary_roc,
    summary_gt     = summary_gt
  )
}

# 5. Latent factor correlations summary gt -----------------------------
summarise_factor_corr_gt <- function(res_factor_corr) {
  if (is.null(res_factor_corr) || nrow(res_factor_corr) == 0) return(NULL)
  
  summary_corr <- res_factor_corr %>%
    group_by(F1, F2) %>%
    summarise(
      N_Runs     = n(),
      Mean_Corr  = mean(Corr, na.rm = TRUE),
      Lower_Corr = quantile(Corr, 0.025, na.rm = TRUE),
      Upper_Corr = quantile(Corr, 0.975, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      `Corr (95% CI)` = format_est_ci(Mean_Corr, Lower_Corr, Upper_Corr)
    )
  
  tab_corr <- summary_corr %>%
    select(
      F1,
      F2,
      `Corr (95% CI)`,
      N_Runs
    )
  
  summary_gt <- gt(tab_corr) %>%
    tab_header(
      title = "Latent Factor Correlations Across Monte Carlo Samples"
    ) %>%
    cols_label(
      F1              = "Factor 1",
      F2              = "Factor 2",
      `Corr (95% CI)` = "Correlation (95% CI)",
      N_Runs          = "Number of Runs"
    ) %>%
    tab_source_note(
      source_note = paste(
        "Latent factor correlations estimated from CFA models across Monte Carlo splits,",
        "reported as mean and 95% Monte Carlo quantile interval."
      )
    )
  
  list(
    summary_source = summary_corr,
    summary_gt     = summary_gt
  )
}

# 6. Focused ROC comparison: Overall vs Overall (Excl. 0-0) ------------
summarise_overall_restricted_roc_gt <- function(res_roc) {
  if (is.null(res_roc) || nrow(res_roc) == 0) return(NULL)
  
  comp_df <- res_roc %>%
    filter(Group %in% c("Overall", "Overall (Excl. 0-0)")) %>%
    group_by(Group) %>%
    summarise(
      N_Runs         = n(),
      Mean_AUC       = mean(AUC, na.rm = TRUE),
      Lower_AUC      = quantile(AUC, 0.025, na.rm = TRUE),
      Upper_AUC      = quantile(AUC, 0.975, na.rm = TRUE),
      Mean_Sens      = mean(Sensitivity, na.rm = TRUE),
      Lower_Sens     = quantile(Sensitivity, 0.025, na.rm = TRUE),
      Upper_Sens     = quantile(Sensitivity, 0.975, na.rm = TRUE),
      Mean_Spec      = mean(Specificity, na.rm = TRUE),
      Lower_Spec     = quantile(Specificity, 0.025, na.rm = TRUE),
      Upper_Spec     = quantile(Specificity, 0.975, na.rm = TRUE),
      Mean_Kappa     = mean(Kappa, na.rm = TRUE),
      Lower_Kappa    = quantile(Kappa, 0.025, na.rm = TRUE),
      Upper_Kappa    = quantile(Kappa, 0.975, na.rm = TRUE),
      Mean_Cutoff    = if ("Optimal_Cutoff" %in% colnames(cur_data_all())) mean(Optimal_Cutoff, na.rm = TRUE) else NA_real_,
      Lower_Cutoff   = if ("Optimal_Cutoff" %in% colnames(cur_data_all())) quantile(Optimal_Cutoff, 0.025, na.rm = TRUE) else NA_real_,
      Upper_Cutoff   = if ("Optimal_Cutoff" %in% colnames(cur_data_all())) quantile(Optimal_Cutoff, 0.975, na.rm = TRUE) else NA_real_,
      Mean_Prev_PHQ  = if ("Prevalence_PHQ" %in% colnames(cur_data_all())) mean(Prevalence_PHQ, na.rm = TRUE) else NA_real_,
      Lower_Prev_PHQ = if ("Prevalence_PHQ" %in% colnames(cur_data_all())) quantile(Prevalence_PHQ, 0.025, na.rm = TRUE) else NA_real_,
      Upper_Prev_PHQ = if ("Prevalence_PHQ" %in% colnames(cur_data_all())) quantile(Prevalence_PHQ, 0.975, na.rm = TRUE) else NA_real_,
      Mean_Prev_SSQ  = if ("Prevalence_SSQ" %in% colnames(cur_data_all())) mean(Prevalence_SSQ, na.rm = TRUE) else NA_real_,
      Lower_Prev_SSQ = if ("Prevalence_SSQ" %in% colnames(cur_data_all())) quantile(Prevalence_SSQ, 0.025, na.rm = TRUE) else NA_real_,
      Upper_Prev_SSQ = if ("Prevalence_SSQ" %in% colnames(cur_data_all())) quantile(Prevalence_SSQ, 0.975, na.rm = TRUE) else NA_real_,
      .groups        = "drop"
    ) %>%
    mutate(
      `AUC (95% CI)`        = format_est_ci(Mean_AUC,      Lower_AUC,      Upper_AUC),
      `Sens (95% CI)`       = format_est_ci(Mean_Sens,     Lower_Sens,     Upper_Sens),
      `Spec (95% CI)`       = format_est_ci(Mean_Spec,     Lower_Spec,     Upper_Spec),
      `Kappa (95% CI)`      = format_est_ci(Mean_Kappa,    Lower_Kappa,    Upper_Kappa),
      `Cutoff (95% CI)`     = format_est_ci(Mean_Cutoff,   Lower_Cutoff,   Upper_Cutoff),
      `Prev PHQ (95% CI)`   = format_est_ci(Mean_Prev_PHQ * 100, Lower_Prev_PHQ * 100, Upper_Prev_PHQ * 100, digits = 1),
      `Prev SSQ (95% CI)`   = format_est_ci(Mean_Prev_SSQ * 100, Lower_Prev_SSQ * 100, Upper_Prev_SSQ * 100, digits = 1)
    )
  
  if (nrow(comp_df) == 0) return(NULL)
  
  tab_comp <- comp_df %>%
    select(
      Group,
      `AUC (95% CI)`,
      `Sens (95% CI)`,
      `Spec (95% CI)`,
      `Kappa (95% CI)`,
      `Cutoff (95% CI)`,
      `Prev PHQ (95% CI)`,
      `Prev SSQ (95% CI)`,
      N_Runs
    )
  
  summary_gt <- gt(tab_comp) %>%
    tab_header(
      title = "ROC Comparison: Full Sample vs Restricted (Excl. 0-0)"
    ) %>%
    cols_label(
      Group               = "Group",
      `AUC (95% CI)`      = "AUC (95% CI)",
      `Sens (95% CI)`     = "Sensitivity (95% CI)",
      `Spec (95% CI)`     = "Specificity (95% CI)",
      `Kappa (95% CI)`    = "Cohen's Kappa (95% CI)",
      `Cutoff (95% CI)`   = "Optimal SSQ-10 Cutoff (95% CI)",
      `Prev PHQ (95% CI)` = "Prevalence PHQ-09 ≥ 10 (95% CI)",
      `Prev SSQ (95% CI)` = "Prevalence SSQ-10 above cutoff (95% CI)",
      N_Runs              = "Number of Runs"
    ) %>%
    fmt_number(
      columns  = c(N_Runs),
      decimals = 0
    ) %>%
    tab_source_note(
      source_note = paste(
        "Comparison of ROC metrics between the full sample (\"Overall\") and the restricted sample",
        "(\"Overall (Excl. 0-0)\") where participants with both SSQ-10 = 0 and PHQ-09 = 0",
        "were excluded from the ROC calculation. Prevalence is reported for both PHQ-09 and SSQ-10.")
    )
  
  list(
    summary_source = comp_df,
    summary_gt     = summary_gt
  )
}

# 7. Reliability summary gt (Cronbach's alpha, McDonald's omega) ------- (Cronbach's alpha, McDonald's omega) -------
# Expect res_reliability to contain rows from semTools::reliability(fit)
# with columns including alpha and omega (or similar), plus Scale and Seed.
summarise_reliability_gt <- function(res_reliability) {
  if (is.null(res_reliability) || nrow(res_reliability) == 0) return(NULL)
  
  # Try to detect alpha and omega column names flexibly
  alpha_col <- intersect(c("alpha", "alpha.total", "Alpha"), colnames(res_reliability))
  omega_col <- intersect(c("omega", "omega.tot", "omega.total"), colnames(res_reliability))
  
  if (length(alpha_col) == 0 && length(omega_col) == 0) {
    # Nothing usable
    return(NULL)
  }
  
  summary_rel <- res_reliability %>%
    group_by(Scale) %>%
    summarise(
      N_Runs       = n(),
      Mean_alpha   = if (length(alpha_col) > 0) mean(.data[[alpha_col[1]]], na.rm = TRUE) else NA_real_,
      Lower_alpha  = if (length(alpha_col) > 0) quantile(.data[[alpha_col[1]]], 0.025, na.rm = TRUE) else NA_real_,
      Upper_alpha  = if (length(alpha_col) > 0) quantile(.data[[alpha_col[1]]], 0.975, na.rm = TRUE) else NA_real_,
      Mean_omega   = if (length(omega_col) > 0) mean(.data[[omega_col[1]]], na.rm = TRUE) else NA_real_,
      Lower_omega  = if (length(omega_col) > 0) quantile(.data[[omega_col[1]]], 0.025, na.rm = TRUE) else NA_real_,
      Upper_omega  = if (length(omega_col) > 0) quantile(.data[[omega_col[1]]], 0.975, na.rm = TRUE) else NA_real_,
      .groups      = "drop"
    ) %>%
    mutate(
      `Alpha (95% CI)` = if (length(alpha_col) > 0) format_est_ci(Mean_alpha, Lower_alpha, Upper_alpha) else NA_character_,
      `Omega (95% CI)` = if (length(omega_col) > 0) format_est_ci(Mean_omega, Lower_omega, Upper_omega) else NA_character_
    )
  
  tab_rel <- summary_rel %>%
    select(
      Scale,
      `Alpha (95% CI)`,
      `Omega (95% CI)`,
      N_Runs
    )
  
  summary_gt <- gt(tab_rel) %>%
    tab_header(
      title = "Internal Consistency Reliability Across Monte Carlo Samples"
    ) %>%
    cols_label(
      Scale            = "Scale / Factor",
      `Alpha (95% CI)` = "Cronbach's α (95% CI)",
      `Omega (95% CI)` = "McDonald's ω (95% CI)",
      N_Runs           = "Number of Runs"
    ) %>%
    tab_source_note(
      source_note = paste(
        "Interpretation guide: α, ω ≥ 0.70 acceptable, ≥ 0.80 good, ≥ 0.90 excellent.",
        "Estimates are averaged across Monte Carlo splits with 95% quantile intervals."
      )
    )
  
  list(
    summary_source = summary_rel,
    summary_gt     = summary_gt
  )
}

# Master wrapper: summarise_results ------------------------------------
summarise_results <- function(res_fit,
                              res_loadings,
                              res_invariance,
                              res_roc,
                              res_factor_corr = NULL,
                              res_reliability = NULL) {
  gt_fit                    <- summarise_fit_gt(res_fit)
  gt_loadings               <- summarise_loadings_gt(res_loadings)
  gt_invariance             <- summarise_invariance_gt(res_invariance)
  gt_roc                    <- summarise_roc_gt(res_roc)
  gt_roc_overall_restricted <- summarise_overall_restricted_roc_gt(res_roc)
  gt_factor_corr            <- summarise_factor_corr_gt(res_factor_corr)
  gt_reliability            <- summarise_reliability_gt(res_reliability)
  
  list(
    gt_fit_source                    = if (!is.null(gt_fit)) gt_fit[["summary_source"]] else NULL,
    fit                              = if (!is.null(gt_fit)) gt_fit[["summary_gt"]] else NULL,
    
    gt_loadings_source               = if (!is.null(gt_loadings)) gt_loadings[["summary_source"]] else NULL,
    loadings                         = if (!is.null(gt_loadings)) gt_loadings[["summary_gt"]] else NULL,
    
    gt_invariance_source             = if (!is.null(gt_invariance)) gt_invariance[["summary_source"]] else NULL,
    invariance                       = if (!is.null(gt_invariance)) gt_invariance[["summary_gt"]] else NULL,
    
    gt_roc_source                    = if (!is.null(gt_roc)) gt_roc[["summary_source"]] else NULL,
    roc                              = if (!is.null(gt_roc)) gt_roc[["summary_gt"]] else NULL,
    
    gt_roc_overall_restricted_source = if (!is.null(gt_roc_overall_restricted)) gt_roc_overall_restricted[["summary_source"]] else NULL,
    roc_overall_restricted           = if (!is.null(gt_roc_overall_restricted)) gt_roc_overall_restricted[["summary_gt"]] else NULL,
    
    gt_factor_correlations_source    = if (!is.null(gt_factor_corr)) gt_factor_corr[["summary_source"]] else NULL,
    factor_correlations              = if (!is.null(gt_factor_corr)) gt_factor_corr[["summary_gt"]] else NULL,
    
    gt_reliability_source            = if (!is.null(gt_reliability)) gt_reliability[["summary_source"]] else NULL,
    reliability                      = if (!is.null(gt_reliability)) gt_reliability[["summary_gt"]] else NULL
  )
}
