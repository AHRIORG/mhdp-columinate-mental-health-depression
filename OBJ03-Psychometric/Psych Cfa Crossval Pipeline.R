# Load necessary libraries
if (!require("here")) install.packages("here")
if (!require("psych")) install.packages("psych")
if (!require("lavaan")) install.packages("lavaan")
if (!require("semTools")) install.packages("semTools")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("pROC")) install.packages("pROC")
if (!require("ggplot2")) install.packages("ggplot2")

library(here)
library(psych)
library(lavaan)
library(semTools)
library(dplyr)
library(tidyr)
library(pROC)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Configuration & Data Loading
# -----------------------------------------------------------------------------

# Load Main Data
data_path <- here("../../data_management/data/co-luminate/adam/dt_psychometric.RData")
if (file.exists(data_path)) {
  load(file = data_path)
  # Ensure the object is named 'dt_main'
  if (!exists("dt_main")) dt_main <- get(ls()[1]) 
  message("Data loaded successfully.")
} else {
  stop("Data file not found.")
}

# Define Item Sets
item_sets_list <- list(
  "SSQ-14 (Original)" = paste0("SSQ", sprintf("%02d", 1:14)),
  "SSQ-10 (Theoretical Depressive Symptoms)" = paste0("SSQ", sprintf("%02d", c(1:3, 8:14))),
  "SSQ-08 (Dixon Chibanda et.al)" = paste0("SSQ", sprintf("%02d", c(1, 8:14))),
  "PHQ-9 (Validation)" = paste0("PHQ9", sprintf("%02d", 1:9))
)

# --- CONFIGURATION: Select which set to validate ---
target_set_name <- "SSQ-10 (Theoretical Depressive Symptoms)"
target_items    <- item_sets_list[[target_set_name]]
dataset_label   <- "Isisikelo Sempilo"

# --- CONFIGURATION: Consensus Model Override ---
# Leave NULL to programmatically select the most frequent model across methods.
manual_method    <- NULL 
manual_n_factors <- NULL

# --- CONFIGURATION: Monte Carlo Settings ---
n_mc_samples <- 100 # Number of cross-validation splits to run

message(paste("Target Set:", target_set_name))
message(paste("Items:", paste(target_items, collapse = ", ")))

# -----------------------------------------------------------------------------
# 2. Identify Consensus Structure & Generate Syntax
# -----------------------------------------------------------------------------

# Construct filename based on target set name
safe_set_name <- gsub("[\n ]+", "_", target_set_name)
safe_ds_label <- gsub("[\n ]+", "_", dataset_label)
# Note: Adjust 'boot_1000' if you ran a different number of EFA boots
boot_file <- paste0("boot_1000_", safe_ds_label, "_", safe_set_name, "_results.rds")
boot_path <- here("statistical_analysis/output/objects", boot_file)

if (!file.exists(boot_path)) stop("EFA Bootstrap results not found. Run EFA script first.")

boot_results <- readRDS(boot_path)

# A. Determine Winning Factor Count/Method
if (is.null(manual_method) || is.null(manual_n_factors)) {
  k_counts <- table(boot_results$k_factor_f)
  p_counts <- table(boot_results$p_factor_f)
  
  k_win_n <- as.numeric(names(which.max(k_counts)))
  p_win_n <- as.numeric(names(which.max(p_counts)))
  
  if (max(k_counts) >= max(p_counts)) {
    target_method <- "Kaiser"
    target_n <- k_win_n
  } else {
    target_method <- "Parallel"
    target_n <- p_win_n
  }
  message(sprintf("Consensus: %s suggested %d Factors", target_method, target_n))
} else {
  target_method <- manual_method
  target_n <- manual_n_factors
  message(sprintf("Manual Override: %s with %d Factors", target_method, target_n))
}

# B. Identify Winning Item Structure (Signature)
target_col_prefix <- if (target_method == "Kaiser") "k_factor_" else "p_factor_"
cols_to_check <- paste0(target_col_prefix, 1:target_n)

# Filter to valid runs of the winning type
valid_runs <- boot_results %>%
  filter(
    Is_Suitable == TRUE,
    if(target_method=="Kaiser") k_factor_f == target_n else p_factor_f == target_n
  )

# Create structure signatures
signatures <- apply(valid_runs[, cols_to_check, drop=FALSE], 1, paste, collapse="|")
winning_sig <- names(which.max(table(signatures)))

# C. Build Lavaan Syntax
factor_defs <- unlist(strsplit(winning_sig, "\\|"))
cfa_syntax_parts <- c()

for (i in seq_along(factor_defs)) {
  f_items <- trimws(unlist(strsplit(factor_defs[i], ",")))
  cfa_syntax_parts <- c(cfa_syntax_parts, paste0("F", i, " =~ ", paste(f_items, collapse = " + ")))
}
cfa_model_syntax <- paste(cfa_syntax_parts, collapse = "\n")

message("\n--- Consensus CFA Syntax ---")
cat(cfa_model_syntax)
message("\n----------------------------")

# Select seeds for Monte Carlo
target_seeds <- valid_runs$Seed[signatures == winning_sig]

if (length(target_seeds) > n_mc_samples) {
  set.seed(123)
  run_seeds <- sample(target_seeds, n_mc_samples)
} else {
  run_seeds <- target_seeds
}

# -----------------------------------------------------------------------------
# 3. Monte Carlo Cross-Validation Loop
# -----------------------------------------------------------------------------

# --- Helper Functions ---

# Invariance Helper
get_invariance_deltas <- function(data, model, group_var, items) {
  if (!group_var %in% names(data) || length(unique(na.omit(data[[group_var]]))) < 2) return(NULL)
  
  # Check subgroup sizes
  min_n <- min(table(data[[group_var]]))
  if (min_n < 50) return(NULL) 
  
  tryCatch({
    # Configural
    fit_conf <- cfa(model, data = data, group = group_var, estimator = "WLSMV", ordered = items)
    # Scalar (Strong)
    fit_scal <- cfa(model, data = data, group = group_var, group.equal = c("loadings", "thresholds"), estimator = "WLSMV", ordered = items)
    
    fm_conf <- fitMeasures(fit_conf, c("cfi", "rmsea"))
    fm_scal <- fitMeasures(fit_scal, c("cfi", "rmsea"))
    
    data.frame(
      Group = group_var,
      Delta_CFI = fm_conf["cfi"] - fm_scal["cfi"],
      Delta_RMSEA = fm_scal["rmsea"] - fm_conf["rmsea"]
    )
  }, error = function(e) NULL)
}

# ROC Helper
get_roc_metrics <- function(data, seed, group_label, items) {
  if (nrow(data) < 20 || length(unique(data$Depressed)) < 2) return(NULL)
  tryCatch({
    # Calculate Sum Score
    # Note: Data passed here must be the numeric subset for calculation
    # But we need 'Depressed' column. Assumes 'data' has both metadata and numeric items.
    # Since items might be factor in original df, we calculate sum on the fly
    
    # Convert items to numeric 0/1 just for sum
    item_data <- data[, items]
    item_data[] <- lapply(item_data, function(x) if(is.factor(x)) as.numeric(x)-1 else x)
    
    sum_score <- rowSums(item_data, na.rm = TRUE)
    
    roc_obj <- roc(data$Depressed, sum_score, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    coords_res <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
    
    if(nrow(coords_res) > 1) coords_res <- coords_res[1, ]
    
    # Calculate Kappa
    thresh <- coords_res$threshold
    pred <- ifelse(sum_score >= thresh, 1, 0)
    obs <- data$Depressed
    
    # Create confusion matrix
    tbl <- table(factor(pred, levels = c(0, 1)), factor(obs, levels = c(0, 1)))
    
    # Calculate Cohen's Kappa manually or via psych
    # Po = (TP + TN) / N
    # Pe = ( (TP+FP)*(TP+FN) + (TN+FN)*(TN+FP) ) / N^2
    # Kappa = (Po - Pe) / (1 - Pe)
    
    po <- sum(diag(tbl)) / sum(tbl)
    pe <- ( (sum(tbl[2,]) * sum(tbl[,2])) + (sum(tbl[1,]) * sum(tbl[,1])) ) / sum(tbl)^2
    kappa_val <- (po - pe) / (1 - pe)
    
    data.frame(Seed = seed, Group = group_label, AUC = auc_val,
               Optimal_Cutoff = coords_res$threshold, 
               Sensitivity = coords_res$sensitivity, 
               Specificity = coords_res$specificity,
               Kappa = kappa_val)
  }, error = function(e) NULL)
}

# ROC Curve Data Helper
get_roc_curve_data <- function(data, seed, group_label, items) {
  if (nrow(data) < 20 || length(unique(data$Depressed)) < 2) return(NULL)
  tryCatch({
    # Calculate Sum Score
    item_data <- data[, items]
    item_data[] <- lapply(item_data, function(x) if(is.factor(x)) as.numeric(x)-1 else x)
    sum_score <- rowSums(item_data, na.rm = TRUE)
    
    roc_obj <- roc(data$Depressed, sum_score, quiet = TRUE)
    
    # Interpolate sensitivities at fixed specificities for averaging
    # We use a dense grid of specificities (0 to 1)
    specs <- seq(0, 1, 0.01)
    # Note: coords with input="specificity" performs linear interpolation
    senss <- coords(roc_obj, x = specs, input = "specificity", ret = "sensitivity", transpose = FALSE)
    
    data.frame(Seed = seed, Group = group_label, Specificity = specs, Sensitivity = as.numeric(unlist(senss)))
  }, error = function(e) {
    message(paste("Error in get_roc_curve_data for Group:", group_label, "-", e$message))
    NULL
  })
}

# --- Execution Loop ---

# Storage
res_fit <- data.frame()
res_loadings <- data.frame()
res_invariance <- data.frame()
res_roc <- data.frame()
res_roc_curves <- data.frame()

message(sprintf("\nRunning Monte Carlo CFA on %d splits...", length(run_seeds)))
pb <- txtProgressBar(min = 0, max = length(run_seeds), style = 3)

for (i in seq_along(run_seeds)) {
  curr_seed <- run_seeds[i]
  
  # A. Recreate Split
  set.seed(curr_seed)
  n_rows <- nrow(dt_main)
  idx_explore <- sample(seq_len(n_rows), floor(0.5 * n_rows), replace = FALSE)
  idx_confirm <- setdiff(seq_len(n_rows), idx_explore)
  dt_cfa <- dt_main[idx_confirm, ]
  
  # B. Data Prep (Numeric for Lavaan)
  dt_cfa_num <- dt_cfa
  dt_cfa_num[target_items] <- lapply(dt_cfa_num[target_items], function(x) {
    if (is.factor(x)) as.numeric(x) - 1 else x
  })
  
  # Add Grouping Vars
  if ("PHQSCR" %in% names(dt_cfa)) dt_cfa_num$PHQBIN <- ifelse(dt_cfa$PHQSCR >= 10, "Depressed", "Non-Depressed")
  if ("AGE" %in% names(dt_cfa)) dt_cfa_num$AGEGRP <- ifelse(dt_cfa$AGE < 20, "17-19", "20-24")
  if ("SEX" %in% names(dt_cfa)) dt_cfa_num$SEX <- dt_cfa$SEX
  
  # C. Run CFA
  fit <- tryCatch({
    cfa(cfa_model_syntax, data = dt_cfa_num, estimator = "WLSMV", ordered = target_items)
  }, error = function(e) NULL)
  
  if (!is.null(fit) && inspect(fit, "converged")) {
    
    # 1. Fit Indices
    fits <- fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr", "chisq", "df"))
    res_fit <- rbind(res_fit, data.frame(Seed = curr_seed, as.list(fits)))
    
    # 2. Loadings
    est <- standardizedSolution(fit) %>% filter(op == "=~") %>% 
      select(lhs, rhs, est.std, se, pvalue) %>% mutate(Seed = curr_seed)
    res_loadings <- rbind(res_loadings, est)
    
    # 3. Invariance (MGCFA)
    if (target_set_name != "PHQ-9 (Validation)") {
      inv_phq <- get_invariance_deltas(dt_cfa_num, cfa_model_syntax, "PHQBIN", target_items)
      if(!is.null(inv_phq)) res_invariance <- rbind(res_invariance, cbind(Seed=curr_seed, inv_phq))
    }
    inv_sex <- get_invariance_deltas(dt_cfa_num, cfa_model_syntax, "SEX", target_items)
    if(!is.null(inv_sex)) res_invariance <- rbind(res_invariance, cbind(Seed=curr_seed, inv_sex))
    
    inv_age <- get_invariance_deltas(dt_cfa_num, cfa_model_syntax, "AGEGRP", target_items)
    if(!is.null(inv_age)) res_invariance <- rbind(res_invariance, cbind(Seed=curr_seed, inv_age))
    
    # 4. ROC Analysis
    if ("PHQSCR" %in% names(dt_cfa)) {
      roc_dat <- dt_cfa
      roc_dat$Depressed <- ifelse(roc_dat$PHQSCR >= 10, 1, 0)
      
      # Overall
      res_roc <- rbind(res_roc, get_roc_metrics(roc_dat, curr_seed, "Overall", target_items))
      res_roc_curves <- rbind(res_roc_curves, get_roc_curve_data(roc_dat, curr_seed, "Overall", target_items))
      
      # Subgroups
      if ("SEX" %in% names(roc_dat)) {
        for(g in unique(na.omit(roc_dat$SEX))) {
          res_roc <- rbind(res_roc, get_roc_metrics(roc_dat[roc_dat$SEX == g,], curr_seed, paste0("Sex: ", g), target_items))
          res_roc_curves <- rbind(res_roc_curves, get_roc_curve_data(roc_dat[roc_dat$SEX == g,], curr_seed, paste0("Sex: ", g), target_items))
        }
      }
      if ("AGE" %in% names(roc_dat)) {
        roc_dat$AGEGRP <- ifelse(roc_dat$AGE < 20, "17-19", "20-24")
        for(g in unique(na.omit(roc_dat$AGEGRP))) {
          res_roc <- rbind(res_roc, get_roc_metrics(roc_dat[roc_dat$AGEGRP == g,], curr_seed, paste0("Age: ", g), target_items))
          res_roc_curves <- rbind(res_roc_curves, get_roc_curve_data(roc_dat[roc_dat$AGEGRP == g,], curr_seed, paste0("Age: ", g), target_items))
        }
      }
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# -----------------------------------------------------------------------------
# 4. Summary & Output
# -----------------------------------------------------------------------------

message("\n\n========================================================")
message("  MONTE CARLO CROSS-VALIDATION RESULTS")
message("========================================================")

# --- A. Model Fit Summary ---
if(nrow(res_fit) > 0) {
  summary_fit <- res_fit %>%
    summarise(
      N_Valid = n(),
      Mean_CFI = mean(cfi, na.rm=T), Lower_CFI = quantile(cfi, 0.025, na.rm=T), Upper_CFI = quantile(cfi, 0.975, na.rm=T),
      Mean_RMSEA = mean(rmsea, na.rm=T), Lower_RMSEA = quantile(rmsea, 0.025, na.rm=T), Upper_RMSEA = quantile(rmsea, 0.975, na.rm=T),
      Prop_Good_CFI = mean(cfi > 0.90), 
      Prop_Good_RMSEA = mean(rmsea < 0.08)
    )
  message("\n1. Global Fit Indices (Distribution across held-out samples):")
  print(t(summary_fit))
}

# --- B. Factor Loadings Summary ---
if(nrow(res_loadings) > 0) {
  summary_loadings <- res_loadings %>%
    group_by(lhs, rhs) %>%
    summarise(
      Mean_Std_Est = mean(est.std),
      SD_Std_Est = sd(est.std),
      Lower_95 = quantile(est.std, 0.025),
      Upper_95 = quantile(est.std, 0.975),
      .groups = 'drop'
    )
  message("\n2. Standardized Factor Loadings (Aggregated):")
  print(as.data.frame(summary_loadings))
}

# --- C. Measurement Invariance Summary ---
if (nrow(res_invariance) > 0) {
  summary_inv <- res_invariance %>%
    group_by(Group) %>%
    summarise(
      N_Runs = n(),
      Mean_Delta_CFI = round(mean(Delta_CFI), 4),
      Max_Delta_CFI = max(Delta_CFI),
      Prop_Inv_CFI = mean(abs(Delta_CFI) <= 0.010),
      Mean_Delta_RMSEA = mean(Delta_RMSEA),
      Prop_Inv_RMSEA = mean(abs(Delta_RMSEA) <= 0.015)
    )
  
  message("\n3. Measurement Invariance (Scalar vs Configural):")
  print(as.data.frame(summary_inv))
}

# --- D. Diagnostic Accuracy (ROC) Summary ---
if (nrow(res_roc) > 0) {
  summary_roc <- res_roc %>%
    group_by(Group) %>%
    summarise(
      Mean_AUC = mean(AUC), 
      AUC_Lower = quantile(AUC, 0.025), AUC_Upper = quantile(AUC, 0.975),
      Mean_Kappa = mean(Kappa, na.rm=TRUE),
      Median_Cutoff = median(Optimal_Cutoff),
      Mean_Sens = mean(Sensitivity),
      Mean_Spec = mean(Specificity)
    )
  
  message("\n4. Diagnostic Accuracy (Benchmark: PHQ-9 >= 10):")
  print(as.data.frame(summary_roc))
}

# -----------------------------------------------------------------------------
# 5. Save Monte Carlo Results
# -----------------------------------------------------------------------------
output_file <- paste0("mc_cfa_validation_", safe_set_name, ".rds")
saveRDS(list(fit = res_fit, loadings = res_loadings, invariance = res_invariance, roc = res_roc, roc_curves = res_roc_curves), 
        file = here("statistical_analysis/output/objects", output_file))

message(paste("\nFull Monte Carlo results saved to:", output_file))

# -----------------------------------------------------------------------------
# 6. Visualization
# -----------------------------------------------------------------------------

message("\nGenerating distribution plots...")
plot_dir <- here("statistical_analysis/output/plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 1. Fit Indices
if (nrow(res_fit) > 0) {
  p_fit_data <- res_fit %>%
    select(Seed, cfi, tli, rmsea, srmr) %>%
    pivot_longer(-Seed, names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = toupper(Metric))
  
  # Define cutoffs
  cutoffs <- data.frame(
    Metric = c("CFI", "TLI", "RMSEA", "SRMR"),
    Intercept = c(0.90, 0.90, 0.08, 0.08)
  )
  
  p_fit <- ggplot(p_fit_data, aes(x = Value, fill = Metric)) +
    geom_histogram(bins = 20, color = "black", alpha = 0.7) +
    facet_wrap(~Metric, scales = "free", ncol = 2) +
    geom_vline(data = cutoffs, aes(xintercept = Intercept), 
               linetype = "dashed", color = "red", linewidth = 0.8) +
    theme_minimal() +
    labs(title = paste("Distribution of Fit Indices (", n_mc_samples, " Samples)"),
         subtitle = paste("Set:", target_set_name, "| Red lines indicate recommended cutoffs (CFI/TLI >= 0.90, RMSEA/SRMR <= 0.08)"),
         y = "Frequency", x = "Index Value") +
    theme(legend.position = "none")
  
  ggsave(filename = file.path(plot_dir, paste0("dist_fit_indices_", safe_set_name, ".pdf")), 
         plot = p_fit, width = 8, height = 6)
}

# 2. Factor Loadings
if (nrow(res_loadings) > 0) {
  p_load <- ggplot(res_loadings, aes(x = rhs, y = est.std, fill = lhs)) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "white") +
    theme_minimal() +
    labs(title = "Distribution of Standardized Factor Loadings",
         subtitle = paste("Set:", target_set_name),
         x = "Item", y = "Standardized Loading", fill = "Factor") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(filename = file.path(plot_dir, paste0("dist_loadings_", safe_set_name, ".pdf")), 
         plot = p_load, width = 10, height = 6)
}

# 3. Measurement Invariance
if (nrow(res_invariance) > 0) {
  p_inv <- res_invariance %>%
    select(Group, Delta_CFI, Delta_RMSEA) %>%
    pivot_longer(cols = c("Delta_CFI", "Delta_RMSEA"), names_to = "Metric", values_to = "Value") %>%
    ggplot(aes(x = Value, fill = Group)) +
    geom_histogram(bins = 20, color = "black", alpha = 0.6, position = "identity") +
    facet_grid(Group ~ Metric, scales = "free") +
    geom_vline(data = data.frame(Metric = "Delta_CFI", x = -0.01), aes(xintercept = x), linetype = "dashed", color = "red") +
    geom_vline(data = data.frame(Metric = "Delta_RMSEA", x = 0.015), aes(xintercept = x), linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "Distribution of Invariance Deltas (Scalar - Configural)",
         subtitle = "Red dashed lines indicate recommended cutoffs (-0.01 CFI, +0.015 RMSEA)",
         y = "Frequency")
  
  ggsave(filename = file.path(plot_dir, paste0("dist_invariance_", safe_set_name, ".pdf")), 
         plot = p_inv, width = 9, height = 7)
}

# 4. Diagnostic Accuracy (ROC + Kappa)
if (nrow(res_roc) > 0) {
  # Pivot longer to plot AUC, Sensitivity, Specificity, Kappa
  p_roc_data <- res_roc %>%
    select(Seed, Group, AUC, Sensitivity, Specificity, Kappa) %>%
    pivot_longer(cols = c("AUC", "Sensitivity", "Specificity", "Kappa"), 
                 names_to = "Metric", values_to = "Value")
  
  p_roc <- ggplot(p_roc_data, aes(x = Group, y = Value, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    theme_minimal() +
    labs(title = "Distribution of Diagnostic Accuracy Metrics",
         subtitle = paste("Benchmark: PHQ-9 >= 10 | Set:", target_set_name),
         y = "Value", x = "Subgroup") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(filename = file.path(plot_dir, paste0("dist_roc_metrics_", safe_set_name, ".pdf")), 
         plot = p_roc, width = 10, height = 8)
}

# 5. ROC Curves with Ribbons
if (nrow(res_roc_curves) > 0) {
  # Aggregate Curve Data
  roc_agg <- res_roc_curves %>%
    group_by(Group, Specificity) %>%
    summarise(
      Mean_Sens = mean(Sensitivity, na.rm = TRUE),
      Lower_Sens = quantile(Sensitivity, 0.025, na.rm = TRUE),
      Upper_Sens = quantile(Sensitivity, 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(InvSpec = 1 - Specificity)
  
  # Aggregate Optimal Cutoff Points
  roc_pts <- res_roc %>%
    group_by(Group) %>%
    summarise(
      Mean_Sens = mean(Sensitivity, na.rm = TRUE),
      Mean_Spec = mean(Specificity, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(InvSpec = 1 - Mean_Spec)
  
  # Plot Function
  plot_roc_ribbon <- function(data_curve, data_pt, title_suffix) {
    ggplot(data_curve, aes(x = InvSpec, y = Mean_Sens)) +
      geom_line(color = "blue", size = 1) +
      geom_ribbon(aes(ymin = Lower_Sens, ymax = Upper_Sens), alpha = 0.2, fill = "blue") +
      geom_point(data = data_pt, aes(x = InvSpec, y = Mean_Sens), color = "red", size = 3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      theme_minimal() +
      labs(title = paste("ROC Curve with 95% CI:", title_suffix),
           subtitle = paste("Set:", target_set_name, "| Red Point: Mean Optimal Cutoff"),
           x = "1 - Specificity", y = "Sensitivity")
  }
  
  # Overall
  p_roc_overall <- plot_roc_ribbon(
    filter(roc_agg, Group == "Overall"), 
    filter(roc_pts, Group == "Overall"), 
    "Overall"
  )
  ggsave(filename = file.path(plot_dir, paste0("roc_curve_overall_", safe_set_name, ".pdf")), 
         plot = p_roc_overall, width = 6, height = 6)
  
  # Sex
  p_roc_sex <- ggplot(filter(roc_agg, grepl("Sex:", Group)), aes(x = InvSpec, y = Mean_Sens, color = Group, fill = Group)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = Lower_Sens, ymax = Upper_Sens), alpha = 0.2, color = NA) +
    geom_point(data = filter(roc_pts, grepl("Sex:", Group)), aes(x = InvSpec, y = Mean_Sens), size = 3, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_minimal() +
    labs(title = "ROC Curves by Sex with 95% CI",
         subtitle = paste("Set:", target_set_name),
         x = "1 - Specificity", y = "Sensitivity")
  
  ggsave(filename = file.path(plot_dir, paste0("roc_curve_sex_", safe_set_name, ".pdf")), 
         plot = p_roc_sex, width = 7, height = 6)
  
  # Age
  p_roc_age <- ggplot(filter(roc_agg, grepl("Age:", Group)), aes(x = InvSpec, y = Mean_Sens, color = Group, fill = Group)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = Lower_Sens, ymax = Upper_Sens), alpha = 0.2, color = NA) +
    geom_point(data = filter(roc_pts, grepl("Age:", Group)), aes(x = InvSpec, y = Mean_Sens), size = 3, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_minimal() +
    labs(title = "ROC Curves by Age Group with 95% CI",
         subtitle = paste("Set:", target_set_name),
         x = "1 - Specificity", y = "Sensitivity")
  
  ggsave(filename = file.path(plot_dir, paste0("roc_curve_age_", safe_set_name, ".pdf")), 
         plot = p_roc_age, width = 7, height = 6)
}

message(paste("Plots saved to:", plot_dir))