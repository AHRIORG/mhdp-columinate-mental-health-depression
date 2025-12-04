############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 04_summary.R - Summary tables & saving results
############################################################

summarise_results <- function(res_fit, res_loadings, res_invariance, res_roc) {
  message("

========================================================")
  message("  MONTE CARLO CROSS-VALIDATION RESULTS")
  message("========================================================")
  
  gt_fit       <- NULL
  gt_loadings  <- NULL
  gt_inv       <- NULL
  gt_roc       <- NULL
  
  # 1. Global Fit Indices ----------------------------------------------
  if (nrow(res_fit) > 0) {
    summary_fit <- res_fit %>%
      summarise(
        N_Valid         = n(),
        Mean_CFI        = mean(cfi,  na.rm = TRUE),
        Lower_CFI       = quantile(cfi,  0.025, na.rm = TRUE),
        Upper_CFI       = quantile(cfi,  0.975, na.rm = TRUE),
        Mean_RMSEA      = mean(rmsea, na.rm = TRUE),
        Lower_RMSEA     = quantile(rmsea, 0.025, na.rm = TRUE),
        Upper_RMSEA     = quantile(rmsea, 0.975, na.rm = TRUE),
        Prop_Good_CFI   = mean(cfi   > 0.90, na.rm = TRUE),
        Prop_Good_RMSEA = mean(rmsea < 0.08, na.rm = TRUE)
      ) %>%
      mutate(
        CFI_Est_CI   = sprintf("%.3f (%.3f, %.3f)", Mean_CFI,   Lower_CFI,   Upper_CFI),
        RMSEA_Est_CI = sprintf("%.3f (%.3f, %.3f)", Mean_RMSEA, Lower_RMSEA, Upper_RMSEA)
      )
    
    gt_fit <- summary_fit %>%
      gt() %>%
      cols_label(
        N_Valid         = "Valid runs (N)",
        CFI_Est_CI      = "CFI (95% CI)",
        RMSEA_Est_CI    = "RMSEA (95% CI)",
        Prop_Good_CFI   = "Prop. CFI ≥ 0.90",
        Prop_Good_RMSEA = "Prop. RMSEA ≤ 0.08"
      ) %>%
      cols_hide(c(Mean_CFI, Lower_CFI, Upper_CFI,
                  Mean_RMSEA, Lower_RMSEA, Upper_RMSEA,
                  Prop_Good_CFI, Prop_Good_RMSEA)) %>%
      fmt_number(
        columns = c(N_Valid),
        decimals = 0
      ) %>%
      tab_source_note(
        md(
          "**Interpretation:** CFI/TLI ≥ 0.90 (often ≥ 0.95 for *good* fit) and RMSEA ≤ 0.08 (≤ 0.06 for *good* fit) indicate acceptable model fit. Proportions (not shown) correspond to the fraction of Monte Carlo samples meeting these benchmarks."
        )
      )
  }
  
  # 2. Factor Loadings --------------------------------------------------
  if (nrow(res_loadings) > 0) {
    summary_loadings <- res_loadings %>%
      group_by(lhs, rhs) %>%
      summarise(
        Mean_Std_Est = mean(est.std),
        SD_Std_Est   = sd(est.std),
        Lower_95     = quantile(est.std, 0.025),
        Upper_95     = quantile(est.std, 0.975),
        .groups      = "drop"
      ) %>%
      mutate(
        Loading_Est_CI = sprintf("%.3f (%.3f, %.3f)",
                                 Mean_Std_Est, Lower_95, Upper_95)
      )
    
    gt_loadings <- summary_loadings %>%
      arrange(lhs, rhs) %>%
      gt() %>%
      cols_label(
        lhs            = "Factor",
        rhs            = "Item",
        Loading_Est_CI = "Std. loading (95% CI)"
      ) %>%
      cols_hide(c(Mean_Std_Est, SD_Std_Est, Lower_95, Upper_95)) %>%
      tab_source_note(
        md(
          "**Interpretation:** Standardized loadings ≥ 0.40 are typically considered *salient*; ≥ 0.70 are considered *strong*."
        )
      )
  }
  
  # 3. Measurement Invariance -------------------------------------------
  if (nrow(res_invariance) > 0) {
    summary_inv <- res_invariance %>%
      group_by(Group) %>%
      summarise(
        N_Runs            = n(),
        Mean_Delta_CFI    = mean(Delta_CFI),
        Lower_Delta_CFI   = quantile(Delta_CFI, 0.025),
        Upper_Delta_CFI   = quantile(Delta_CFI, 0.975),
        Prop_Inv_CFI      = mean(abs(Delta_CFI) <= 0.010),
        Mean_Delta_RMSEA  = mean(Delta_RMSEA),
        Lower_Delta_RMSEA = quantile(Delta_RMSEA, 0.025),
        Upper_Delta_RMSEA = quantile(Delta_RMSEA, 0.975),
        Prop_Inv_RMSEA    = mean(abs(Delta_RMSEA) <= 0.015)
      ) %>%
      mutate(
        Delta_CFI_Est_CI   = sprintf("%.3f (%.3f, %.3f)",
                                     Mean_Delta_CFI, Lower_Delta_CFI, Upper_Delta_CFI),
        Delta_RMSEA_Est_CI = sprintf("%.3f (%.3f, %.3f)",
                                     Mean_Delta_RMSEA, Lower_Delta_RMSEA, Upper_Delta_RMSEA)
      )
    
    gt_inv <- summary_inv %>%
      arrange(Group) %>%
      gt() %>%
      cols_label(
        Group             = "Grouping variable",
        N_Runs            = "Valid runs (N)",
        Delta_CFI_Est_CI  = "ΔCFI (95% CI)",
        Prop_Inv_CFI      = "Prop. |ΔCFI| ≤ 0.010",
        Delta_RMSEA_Est_CI = "ΔRMSEA (95% CI)",
        Prop_Inv_RMSEA    = "Prop. |ΔRMSEA| ≤ 0.015"
      ) %>%
      cols_hide(c(
        Mean_Delta_CFI, Lower_Delta_CFI, Upper_Delta_CFI,
        Mean_Delta_RMSEA, Lower_Delta_RMSEA, Upper_Delta_RMSEA,
        Prop_Inv_CFI, Prop_Inv_RMSEA
      )) %>%
      fmt_number(columns = c(N_Runs), decimals = 0) %>%
      tab_source_note(
        md(
          "**Interpretation:** Scalar invariance is typically supported when ΔCFI ≥ -0.01 and ΔRMSEA ≤ 0.015 between configural and scalar models. Proportions (not shown) correspond to the fraction of Monte Carlo samples satisfying these criteria across runs."
        )
      )
  }
  
  # 4. Diagnostic Accuracy (ROC) ----------------------------------------
  if (nrow(res_roc) > 0) {
    summary_roc <- res_roc %>%
      group_by(Group) %>%
      summarise(
        Mean_AUC      = mean(AUC),
        AUC_Lower     = quantile(AUC, 0.025),
        AUC_Upper     = quantile(AUC, 0.975),
        Mean_Kappa    = mean(Kappa, na.rm = TRUE),
        Median_Cutoff = median(Optimal_Cutoff),
        Cutoff_Lower  = quantile(Optimal_Cutoff, 0.025),
        Cutoff_Upper  = quantile(Optimal_Cutoff, 0.975),
        Mean_Sens     = mean(Sensitivity),
        Sens_Lower    = quantile(Sensitivity, 0.025),
        Sens_Upper    = quantile(Sensitivity, 0.975),
        Mean_Spec     = mean(Specificity),
        Spec_Lower    = quantile(Specificity, 0.025),
        Spec_Upper    = quantile(Specificity, 0.975)
      ) %>%
      mutate(
        AUC_Est_CI    = sprintf("%.3f (%.3f, %.3f)", Mean_AUC,      AUC_Lower,     AUC_Upper),
        Cutoff_Est_CI = sprintf("%.3f (%.3f, %.3f)", Median_Cutoff, Cutoff_Lower,  Cutoff_Upper),
        Sens_Est_CI   = sprintf("%.3f (%.3f, %.3f)", Mean_Sens,     Sens_Lower,    Sens_Upper),
        Spec_Est_CI   = sprintf("%.3f (%.3f, %.3f)", Mean_Spec,     Spec_Lower,    Spec_Upper)
      )
    
    gt_roc <- summary_roc %>%
      arrange(Group) %>%
      gt() %>%
      cols_label(
        Group         = "Group",
        AUC_Est_CI    = "AUC (95% CI)",
        Mean_Kappa    = "Mean κ",
        Cutoff_Est_CI = "Cut-off (95% CI)",
        Sens_Est_CI   = "Sensitivity (95% CI)",
        Spec_Est_CI   = "Specificity (95% CI)"
      ) %>%
      cols_hide(c(
        Mean_AUC, AUC_Lower, AUC_Upper,
        Median_Cutoff, Cutoff_Lower, Cutoff_Upper,
        Mean_Sens, Sens_Lower, Sens_Upper,
        Mean_Spec, Spec_Lower, Spec_Upper
      )) %>%
      fmt_number(columns = c(Mean_Kappa), decimals = 3) %>%
      tab_source_note(
        md(
          "**Interpretation:** AUC 0.70–0.80 = *acceptable*, 0.80–0.90 = *excellent*, >0.90 = *outstanding* discrimination. Cohen's κ < 0.40 = *poor–fair*, 0.40–0.60 = *moderate*, 0.60–0.80 = *substantial*, >0.80 = *almost perfect* agreement. Higher sensitivity/specificity indicate better screening performance relative to PHQ-9 ≥ 10 as the reference."
        )
      )
  }
  
  # Return a list of gt tables (entries may be NULL if corresponding data is empty)
  list(
    fit        = gt_fit,
    loadings   = gt_loadings,
    invariance = gt_inv,
    roc        = gt_roc
  )
}

save_mc_results <- function(res_fit, res_loadings, res_invariance,
                            res_roc, res_roc_curves, safe_set_name) {
  output_file <- paste0("mc_cfa_validation_", safe_set_name, ".rds")
  saveRDS(
    list(
      fit        = res_fit,
      loadings   = res_loadings,
      invariance = res_invariance,
      roc        = res_roc,
      roc_curves = res_roc_curves
    ),
    file = here("statistical_analysis/output/objects", output_file)
  )
  message(paste("\nFull Monte Carlo results saved to:", output_file))
}
