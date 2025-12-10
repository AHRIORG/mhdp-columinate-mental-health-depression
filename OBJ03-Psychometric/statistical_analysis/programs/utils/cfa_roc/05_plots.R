############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 05_plots.R - Visualisation helpers
############################################################

generate_plots <- function(res_fit, res_loadings, res_invariance,
                           res_roc, res_roc_curves,
                           safe_set_name, target_set_name, n_mc_samples) {
  message("\nGenerating distribution plots...")
  
  figure_dir <- here("statistical_analysis/output/figures")
  if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)
  
  # Initialise plot objects
  p_fit <- p_load <- p_inv <- p_roc <- NULL
  p_roc_overall <- p_roc_sex <- p_roc_age <- NULL
  
  # 1. Fit Indices ------------------------------------------------------
  if (nrow(res_fit) > 0) {
    p_fit_data <- res_fit %>%
      select(Seed, cfi, tli, rmsea, srmr) %>%
      pivot_longer(-Seed, names_to = "Metric", values_to = "Value") %>%
      mutate(Metric = toupper(Metric))
    
    cutoffs <- data.frame(
      Metric    = c("CFI", "TLI", "RMSEA", "SRMR"),
      Intercept = c(0.90, 0.90, 0.08, 0.08)
    )
    
    p_fit <- ggplot(p_fit_data, aes(x = Value, fill = Metric)) +
      geom_histogram(bins = 20, color = "black", alpha = 0.7) +
      facet_wrap(~Metric, scales = "free", ncol = 4) +
      geom_vline(data = cutoffs, aes(xintercept = Intercept),
                 linetype = "dashed", color = "red", linewidth = 0.8) +
      scale_fill_npg() +
      theme_minimal() +
      labs(
        title   = paste("Distribution of Fit Indices (", n_mc_samples, " Samples)", sep = ""),
        x       = "Index value",
        y       = "Frequency",
        caption = paste(
          "Set:", target_set_name,
          "| Red lines: CFI/TLI ≥ 0.90, RMSEA/SRMR ≤ 0.08"
        )
      ) +
      theme(legend.position = "none")
    
    ggsave(
      filename = file.path(figure_dir, paste0("dist_fit_indices_", safe_set_name, ".pdf")),
      plot     = p_fit,
      width    = 8,
      height   = 6
    )
  }
  
  # 2. Factor Loadings --------------------------------------------------
  if (nrow(res_loadings) > 0) {
    p_load <- ggplot(res_loadings, aes(x = rhs, y = est.std, fill = lhs)) +
      geom_boxplot(outlier.shape = 21, outlier.fill = "white") +
      scale_fill_npg() +
      theme_minimal() +
      labs(
        title   = "Distribution of Standardized Factor Loadings",
        x       = "Item",
        y       = "Standardized loading",
        fill    = "Factor",
        caption = paste("Set:", target_set_name)
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      filename = file.path(figure_dir, paste0("dist_loadings_", safe_set_name, ".pdf")),
      plot     = p_load,
      width    = 10,
      height   = 6
    )
  }
  
  # 3. Measurement Invariance ------------------------------------------
  if (nrow(res_invariance) > 0) {
    p_inv_data <- res_invariance %>%
      select(Group, Delta_CFI, Delta_RMSEA) %>%
      pivot_longer(cols = c("Delta_CFI", "Delta_RMSEA"),
                   names_to = "Metric", values_to = "Value")
    
    p_inv <- ggplot(p_inv_data, aes(x = Value, fill = Group)) +
      geom_histogram(bins = 20, color = "black", alpha = 0.6, position = "identity") +
      facet_grid(Group ~ Metric, scales = "free") +
      scale_fill_npg() +
      theme_minimal() +
      labs(
        title   = "Distribution of Invariance Deltas (Scalar - Configural)",
        x       = "Delta value",
        y       = "Frequency",
        caption = paste(
          "Cutoffs: ΔCFI ≥ -0.01, ΔRMSEA ≤ 0.015 | Groups include AGEGRP, SEX,",
          "PHQBIN (Depressed vs Non-Depressed)"
        )
      )
    
    ggsave(
      filename = file.path(figure_dir, paste0("dist_invariance_", safe_set_name, ".pdf")),
      plot     = p_inv,
      width    = 9,
      height   = 7
    )
  }
  
  # 4. Diagnostic Accuracy (ROC + Kappa) -------------------------------
  if (nrow(res_roc) > 0) {
    p_roc_data <- res_roc %>%
      select(Seed, Group, AUC, Sensitivity, Specificity, Kappa) %>%
      pivot_longer(cols = c("AUC", "Sensitivity", "Specificity", "Kappa"),
                   names_to = "Metric", values_to = "Value")
    
    p_roc <- ggplot(p_roc_data, aes(x = Group, y = Value, fill = Group)) +
      geom_violin(trim = FALSE, alpha = 0.5) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
      facet_wrap(~Metric, scales = "free_y", ncol = 2) +
      scale_fill_npg() +
      theme_minimal() +
      labs(
        title   = "Distribution of Diagnostic Accuracy Metrics",
        x       = "Subgroup",
        y       = "Value",
        caption = paste(
          "Benchmark: PHQ-9 ≥ 10 | Set:", target_set_name,
          "| AUC: 0.70–0.80 acceptable, 0.80–0.90 excellent, >0.90 outstanding"
        )
      ) +
      theme(
        legend.position = "none",
        axis.text.x     = element_text(angle = 45, hjust = 1)
      )
    
    ggsave(
      filename = file.path(figure_dir, paste0("dist_roc_metrics_", safe_set_name, ".pdf")),
      plot     = p_roc,
      width    = 10,
      height   = 8
    )
  }
  
  # 5. ROC Curves with Ribbons -----------------------------------------
  if (nrow(res_roc_curves) > 0) {
    roc_agg <- res_roc_curves %>%
      group_by(Group, Specificity) %>%
      summarise(
        Mean_Sens  = mean(Sensitivity, na.rm = TRUE),
        Lower_Sens = quantile(Sensitivity, 0.025, na.rm = TRUE),
        Upper_Sens = quantile(Sensitivity, 0.975, na.rm = TRUE),
        .groups    = "drop"
      ) %>%
      mutate(InvSpec = 1 - Specificity)
    
    roc_pts <- res_roc %>%
      group_by(Group) %>%
      summarise(
        Mean_Sens = mean(Sensitivity, na.rm = TRUE),
        Mean_Spec = mean(Specificity, na.rm = TRUE),
        .groups   = "drop"
      ) %>%
      mutate(InvSpec = 1 - Mean_Spec)
    
    plot_roc_ribbon <- function(data_curve, data_pt, title_suffix) {
      ggplot(data_curve, aes(x = InvSpec, y = Mean_Sens)) +
        geom_line(linewidth = 1) +
        geom_ribbon(aes(ymin = Lower_Sens, ymax = Upper_Sens), alpha = 0.2) +
        geom_point(data = data_pt, aes(x = InvSpec, y = Mean_Sens),
                   color = "red", size = 3) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
        theme_minimal() +
        labs(
          title   = paste("ROC Curve with 95% CI:", title_suffix),
          x       = "1 - Specificity",
          y       = "Sensitivity",
          caption = paste(
            "Set:", target_set_name,
            "| Red point: mean optimal cutoff"
          )
        )
    }
    
    # Overall
    if ("Overall" %in% unique(roc_agg$Group)) {
      p_roc_overall <- plot_roc_ribbon(
        filter(roc_agg, Group == "Overall"),
        filter(roc_pts, Group == "Overall"),
        "Overall"
      )
      
      ggsave(
        filename = file.path(figure_dir, paste0("roc_curve_overall_", safe_set_name, ".pdf")),
        plot     = p_roc_overall,
        width    = 6,
        height   = 6
      )
    }
    
    # Sex groups
    if (any(grepl("^Sex:", roc_agg$Group))) {
      p_roc_sex <- ggplot(filter(roc_agg, grepl("^Sex:", Group)),
                          aes(x = InvSpec, y = Mean_Sens, color = Group, fill = Group)) +
        geom_line(linewidth = 1) +
        geom_ribbon(aes(ymin = Lower_Sens, ymax = Upper_Sens), alpha = 0.2, color = NA) +
        geom_point(data = filter(roc_pts, grepl("^Sex:", Group)),
                   aes(x = InvSpec, y = Mean_Sens), size = 3, show.legend = FALSE) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
        scale_color_npg() +
        scale_fill_npg() +
        theme_minimal() +
        labs(
          title   = "ROC Curves by Sex with 95% CI",
          x       = "1 - Specificity",
          y       = "Sensitivity",
          caption = paste("Set:", target_set_name)
        )
      
      ggsave(
        filename = file.path(figure_dir, paste0("roc_curve_sex_", safe_set_name, ".pdf")),
        plot     = p_roc_sex,
        width    = 7,
        height   = 6
      )
    }
    
    # Age groups
    if (any(grepl("^Age:", roc_agg$Group))) {
      p_roc_age <- ggplot(filter(roc_agg, grepl("^Age:", Group)),
                          aes(x = InvSpec, y = Mean_Sens, color = Group, fill = Group)) +
        geom_line(linewidth = 1) +
        geom_ribbon(aes(ymin = Lower_Sens, ymax = Upper_Sens), alpha = 0.2, color = NA) +
        geom_point(data = filter(roc_pts, grepl("^Age:", Group)),
                   aes(x = InvSpec, y = Mean_Sens), size = 3, show.legend = FALSE) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
        scale_color_npg() +
        scale_fill_npg() +
        theme_minimal() +
        labs(
          title   = "ROC Curves by Age Group with 95% CI",
          x       = "1 - Specificity",
          y       = "Sensitivity",
          caption = paste("Set:", target_set_name)
        )
      
      ggsave(
        filename = file.path(figure_dir, paste0("roc_curve_age_", safe_set_name, ".pdf")),
        plot     = p_roc_age,
        width    = 7,
        height   = 6
      )
    }
  }
  
  # Collect all plot objects into a list and return --------------------
  plot_list <- list(
    fit_indices  = p_fit,
    loadings     = p_load,
    invariance   = p_inv,
    roc_metrics  = p_roc,
    roc_overall  = p_roc_overall,
    roc_sex      = p_roc_sex,
    roc_age      = p_roc_age
  )
  
  plot_list
}

