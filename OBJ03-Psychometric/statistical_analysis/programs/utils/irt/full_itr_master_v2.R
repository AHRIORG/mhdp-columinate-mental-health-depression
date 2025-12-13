# Script 03: Latent Harmonization via Item Response Theory (IRT)
# ==============================================================================
# PURPOSE:
# 1. Fit Graded Response Models (GRM) to SSQ-10 and PHQ-9.
#    (Supports both Exploratory and Confirmatory/Syntax-based models)
# 2. Check for Differential Item Functioning (DIF) / Group Differences (Age/Sex).
# 3. Estimate latent trait scores (Theta) for all participants.
# 4. Establish a latent cut-off on the Theta scale equivalent to PHQ-9 >= 10.
# 5. Generate harmonized binary classification based on Theta.
# 6. EXPORT SCORING ENGINE: Create tools to score NEW datasets (SSQ only).
# ==============================================================================

# Load necessary libraries
if (!require("here")) install.packages("here")
if (!require("mirt")) install.packages("mirt")     # Core IRT package
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("tidyr")) install.packages("tidyr")
if (!require("psych")) install.packages("psych")   # For Kappa
if (!require("gt")) install.packages("gt")         # For Tables

library(here)
library(mirt)
library(dplyr)
library(ggplot2)
library(tidyr)
library(psych)
library(gt)

# ==============================================================================
# A. HELPER FUNCTIONS
# ==============================================================================

# -----------------------------------------------------------------------------
# 0. Syntax Converter (Lavaan -> Mirt)
#    FIX: Proper pairwise COV statement for >2 factors
# -----------------------------------------------------------------------------
lavaan_to_mirt <- function(syntax) {
  # If numeric (e.g. 1 or 2 factors exploratory), return as is
  if (is.numeric(syntax) || is.null(syntax)) return(syntax)
  
  # Basic cleaning
  s <- gsub("=~", "=", syntax)
  s <- gsub("\\+", ",", s)
  
  # Identify factor names
  lines <- trimws(unlist(strsplit(s, "\n")))
  lines <- lines[lines != ""]
  factor_lines <- grep("=", lines)
  factor_names <- unique(sapply(strsplit(lines[factor_lines], "="), function(x) trimws(x[1])))
  
  # If > 1 factor, append pairwise COV statement
  if (length(factor_names) > 1) {
    pairs <- combn(factor_names, 2, simplify = FALSE)
    cov_terms <- vapply(pairs, function(p) paste0(p[1], "*", p[2]), character(1))
    cov_stmt <- paste0("COV = ", paste(cov_terms, collapse = ", "))
    s <- paste(s, cov_stmt, sep = "\n")
  }
  
  return(mirt.model(s))
}

# -----------------------------------------------------------------------------
# 1. Data Loading & Preparation
# -----------------------------------------------------------------------------
load_and_clean_data <- function(ssq_idx, phq_idx) {
  # Locate Data
  possible_paths <- c(
    here("../../_private_use/data_management/data/co-luminate/adam/dt_psychometric.RData")
  )
  
  data_path <- NULL
  for (p in possible_paths) {
    if (file.exists(p)) {
      data_path <- p
      break
    }
  }
  
  if (is.null(data_path)) stop("Data file not found.")
  
  message(paste("Loading data from:", data_path))
  loaded_objects <- load(file = data_path)
  
  if ("dt_psychometric" %in% loaded_objects) {
    dt_main <- dt_psychometric |> filter(!(SSQSCR==0 & PHQSCR>0))
  } else if ("dt_main" %in% loaded_objects) {
    # dt_main already loaded
  } else {
    dt_main <- get(loaded_objects[1])
  }
  
  # Define Items
  ssq_items <- paste0("SSQ", sprintf("%02d", ssq_idx))
  phq_items <- paste0("PHQ9", sprintf("%02d", phq_idx))
  
  # Subset and Clean
  df_irt <- dt_main[, c("USUBJID", "SEX", "AGE", ssq_items, phq_items)]
  
  if ("AGE" %in% names(df_irt)) {
    df_irt$AGEGRP <- ifelse(df_irt$AGE < 20, "Age_17_19", "Age_20_24")
  }
  
  # Convert Factors to Numeric (0-based)
  df_irt[ssq_items] <- lapply(df_irt[ssq_items], function(x) if (is.factor(x)) as.numeric(x) - 1 else x)
  df_irt[phq_items] <- lapply(df_irt[phq_items], function(x) if (is.factor(x)) as.numeric(x) - 1 else x)
  
  return(list(data = df_irt, ssq_items = ssq_items, phq_items = phq_items))
}

# -----------------------------------------------------------------------------
# 2. Fit Statistics Table Generator
# -----------------------------------------------------------------------------
create_fit_table <- function(fit_obj, title_str) {
  target_cols <- c("M2", "df", "p", "RMSEA", "RMSEA_5", "RMSEA_95", "SRMSR", "TLI", "CFI")
  valid_cols <- intersect(target_cols, names(fit_obj))
  
  if (length(valid_cols) == 0) return(NULL)
  
  fit_df <- as.data.frame(fit_obj)[, valid_cols, drop = FALSE]
  
  gt(fit_df) %>%
    tab_header(title = title_str) %>%
    fmt_number(columns = any_of(c("M2", "RMSEA", "RMSEA_5", "RMSEA_95", "SRMSR", "TLI", "CFI")), decimals = 3) %>%
    fmt_number(columns = any_of(c("p")), decimals = 4) %>%
    fmt_number(columns = any_of("df"), decimals = 0) %>%
    cols_label(M2 = "M2 Statistic", df = "df", p = "p-value", RMSEA_5 = "RMSEA (5%)", RMSEA_95 = "RMSEA (95%)") %>%
    tab_source_note(source_note = "Good Fit Criteria: CFI > 0.95, TLI > 0.95, RMSEA < 0.06, SRMSR < 0.08.")
}

# -----------------------------------------------------------------------------
# 3. Model Fitting Wrapper
# -----------------------------------------------------------------------------
fit_grm_models <- function(df_irt, ssq_items, phq_items, model_phq = 1, model_ssq = 1, model_joint = 1) {
  
  # Convert inputs if they are strings (Lavaan syntax)
  mirt_phq <- lavaan_to_mirt(model_phq)
  mirt_ssq <- lavaan_to_mirt(model_ssq)
  mirt_joint <- lavaan_to_mirt(model_joint)
  
  # A. PHQ-9
  message("Fitting GRM to PHQ-9...")
  phq_dat <- df_irt[, phq_items]
  phq_dat <- phq_dat[rowSums(!is.na(phq_dat)) > 0, ]
  mod_phq <- mirt(phq_dat, mirt_phq, itemtype = "graded", verbose = FALSE)
  m2_phq <- M2(mod_phq, type = "C2", calcNull = FALSE)
  
  # B. SSQ-10
  message("Fitting GRM to SSQ-10...")
  ssq_dat <- df_irt[, ssq_items]
  ssq_dat <- ssq_dat[rowSums(!is.na(ssq_dat)) > 0, ]
  mod_ssq <- mirt(ssq_dat, mirt_ssq, itemtype = "graded", verbose = FALSE)
  m2_ssq <- M2(mod_ssq, type = "C2", calcNull = FALSE)
  
  # C. Joint Model
  message("Fitting Joint GRM (Harmonization)...")
  joint_dat <- df_irt[, c(phq_items, ssq_items)]
  joint_dat <- joint_dat[rowSums(!is.na(joint_dat)) > 0, ]
  mod_joint <- mirt(joint_dat, mirt_joint, itemtype = "graded", verbose = FALSE)
  
  # Generate Tables
  tbl_phq <- create_fit_table(m2_phq, "PHQ-9 Model Fit (GRM)")
  tbl_ssq <- create_fit_table(m2_ssq, "SSQ-10 Model Fit (GRM)")
  
  return(list(
    models = list(phq = mod_phq, ssq = mod_ssq, joint = mod_joint),
    tables = list(phq = tbl_phq, ssq = tbl_ssq)
  ))
}

# --- CACHING LOGIC FOR FITTING ---
run_fitting_workflow <- function(df_irt, ssq_items, phq_items, phq_model, ssq_model, joint_model, output_dir) {
  
  if (is.numeric(joint_model)) {
    model_label <- paste0(joint_model, "factor")
  } else {
    model_label <- "custom_structure"
  }
  
  fit_file_name <- paste0("03_IRT_fitted_models_", model_label, ".rds")
  fit_file_path <- file.path(output_dir, fit_file_name)
  run_checks <- TRUE
  
  current_model_defs <- list(phq = phq_model, ssq = ssq_model, joint = joint_model)
  
  if (file.exists(fit_file_path)) {
    saved_obj <- readRDS(fit_file_path)
    if (!is.null(saved_obj$model_defs) && identical(saved_obj$model_defs, current_model_defs)) {
      message(paste0("\nExisting Fitted Models found for ", model_label, ". (Matches current structure)"))
      if (interactive()) {
        if (tolower(substring(readline(prompt = "  Rerun Model Fitting (y/n)? "), 1, 1)) != "y") {
          run_checks <- FALSE
          fit_res <- saved_obj$results
        }
      } else {
        message("  Loading existing models (non-interactive mode).")
        run_checks <- FALSE
        fit_res <- saved_obj$results
      }
    } else {
      message("\nExisting Fitted Models found but model structure has CHANGED.")
      message("  Forcing Re-run.")
    }
  }
  
  if (run_checks) {
    fit_res <- fit_grm_models(df_irt, ssq_items, phq_items,
                              model_phq = phq_model, model_ssq = ssq_model, model_joint = joint_model)
    
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(results = fit_res, model_defs = current_model_defs), file = fit_file_path)
    message(paste("Fitted models saved to:", fit_file_path))
  } else {
    message("Fitted models loaded.")
  }
  
  return(fit_res)
}

# -----------------------------------------------------------------------------
# 4. Invariance Analysis Helpers
#    FIX: Avoid na.omit() on items; drop only missing group labels.
# -----------------------------------------------------------------------------
check_group_invariance <- function(item_data, group_vec, group_name, model_def) {
  if (length(unique(na.omit(group_vec))) < 2) return(NULL)
  
  message(paste0("Checking Invariance by ", group_name, "..."))
  
  keep <- !is.na(group_vec)
  dat_grp <- data.frame(item_data[keep, , drop = FALSE], Group = as.factor(group_vec[keep]))
  
  # Convert syntax if necessary
  mirt_model <- lavaan_to_mirt(model_def)
  
  mod_mg <- multipleGroup(
    dat_grp[, -ncol(dat_grp), drop = FALSE],
    model = mirt_model,
    group = dat_grp$Group,
    invariance = c("slopes", "intercepts", "free_means", "free_var"),
    verbose = FALSE
  )
  
  coefs <- coef(mod_mg, simplify = TRUE)
  dif_res <- DIF(mod_mg, which.par = c("a1", "d"), scheme = "drop", p.adjust = "fdr")
  
  return(list(model = mod_mg, dif = dif_res, latent_pars = coefs))
}

summarize_invariance <- function(res_obj, group_label) {
  if (is.null(res_obj)) return(NULL)
  if (!is.list(res_obj$latent_pars) || length(res_obj$latent_pars) < 2) return(NULL)
  
  grp_names <- names(res_obj$latent_pars)
  ref_grp <- if (!is.null(grp_names)) grp_names[1] else "Group 1"
  foc_grp <- if (!is.null(grp_names)) grp_names[2] else "Group 2"
  
  focal_pars <- res_obj$latent_pars[[2]]$GroupPars
  latent_mean <- NA
  latent_var <- NA
  
  if (!is.null(focal_pars) && length(focal_pars) > 0) {
    if ("MEAN_1" %in% names(focal_pars)) latent_mean <- focal_pars["MEAN_1"] else if (length(focal_pars) >= 1) latent_mean <- focal_pars[1]
    if ("COV_11" %in% names(focal_pars)) latent_var <- focal_pars["COV_11"] else if (length(focal_pars) >= 2) latent_var <- focal_pars[2]
  }
  
  dif_df <- res_obj$dif
  pcol <- if ("adj_p" %in% names(dif_df)) "adj_p" else "p"
  sig_items <- rownames(dif_df)[dif_df[[pcol]] < 0.05]
  if (is.null(sig_items)) sig_items <- character(0)
  
  data.frame(
    Grouping = group_label,
    Reference_Group = ref_grp,
    Focal_Group = foc_grp,
    Latent_Mean_Diff = round(as.numeric(latent_mean), 3),
    Latent_Var_Ratio = round(as.numeric(latent_var), 3),
    DIF_Items = if (length(sig_items) > 0) paste(sig_items, collapse = ", ") else "None",
    DIF_Count = length(sig_items),
    stringsAsFactors = FALSE
  )
}

create_dif_table <- function(res_obj, group_label) {
  if (is.null(res_obj)) return(NULL)
  dif_df <- res_obj$dif
  dif_df$Item <- rownames(dif_df)
  dif_df$Grouping_Variable <- group_label
  
  cols_to_keep <- intersect(c("Grouping_Variable", "Item", "AIC", "SABIC", "HQ", "BIC", "X2", "df", "p", "adj_p"), names(dif_df))
  
  gt(dif_df[, cols_to_keep]) %>%
    tab_header(title = paste("Detailed DIF Statistics:", group_label)) %>%
    fmt_number(columns = any_of(c("AIC", "SABIC", "HQ", "BIC", "X2")), decimals = 2) %>%
    fmt_number(columns = any_of(c("p", "adj_p")), decimals = 4) %>%
    fmt_number(columns = any_of("df"), decimals = 0) %>%
    cols_label(Item = "Item", X2 = "Chi-Sq", adj_p = "Adj P-Val") %>%
    tab_style(
      style = list(cell_fill(color = "#F9E3E3"), cell_text(weight = "bold")),
      locations = if ("adj_p" %in% names(dif_df)) cells_body(rows = adj_p < 0.05) else cells_body(rows = p < 0.05)
    )
}

run_invariance_workflow <- function(df_irt, ssq_items, phq_items, joint_model_def, output_dir) {
  
  if (is.numeric(joint_model_def)) {
    model_label <- paste0(joint_model_def, "factor")
  } else {
    model_label <- "custom_structure"
  }
  
  inv_file_name <- paste0("03_IRT_invariance_models_", model_label, ".rds")
  inv_file_path <- file.path(output_dir, inv_file_name)
  run_checks <- TRUE
  
  if (file.exists(inv_file_path)) {
    saved_obj <- readRDS(inv_file_path)
    if (!is.null(saved_obj$model_def) && identical(saved_obj$model_def, joint_model_def)) {
      message(paste0("\nExisting Invariance Models found for ", model_label, ". (Matches structure)"))
      if (interactive()) {
        if (tolower(substring(readline(prompt = "  Rerun Invariance checks (y/n)? "), 1, 1)) != "y") {
          run_checks <- FALSE
          inv_results <- saved_obj
        }
      } else {
        message("  Loading existing models (non-interactive mode).")
        run_checks <- FALSE
        inv_results <- saved_obj
      }
    } else {
      message("\nExisting Invariance Models found but model structure has CHANGED.")
      message("  Forcing Re-run.")
    }
  }
  
  if (run_checks) {
    joint_data <- df_irt[, c(phq_items, ssq_items)]
    valid_joint <- rowSums(!is.na(joint_data)) > 0
    clean_dat <- joint_data[valid_joint, , drop = FALSE]
    
    res_sex <- check_group_invariance(clean_dat, df_irt$SEX[valid_joint], "SEX", joint_model_def)
    res_age <- check_group_invariance(clean_dat, df_irt$AGEGRP[valid_joint], "AGEGRP", joint_model_def)
    
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    inv_results <- list(sex = res_sex, age = res_age, model_def = joint_model_def)
    saveRDS(inv_results, file = inv_file_path)
    message(paste("Invariance models saved to:", inv_file_path))
  } else {
    message("Invariance models loaded.")
  }
  
  res_sex <- inv_results$sex
  res_age <- inv_results$age
  summ_sex <- summarize_invariance(res_sex, "Biological Sex")
  summ_age <- summarize_invariance(res_age, "Age Group")
  
  tbl_summary <- gt(rbind(summ_sex, summ_age)) %>%
    tab_header(title = paste("Measurement Invariance (DIF) & Latent Impact -", model_label, "Model")) %>%
    fmt_number(columns = c("Latent_Mean_Diff", "Latent_Var_Ratio"), decimals = 3) %>%
    tab_source_note(source_note = "Impact: Latent Mean Diff > 0 indicates Focal group has higher depression. DIF: p < 0.05 (FDR if available).")
  
  tbl_dif_sex <- create_dif_table(res_sex, "Biological Sex")
  tbl_dif_age <- create_dif_table(res_age, "Age Group")
  
  return(list(sex = res_sex, age = res_age, tables = list(summary = tbl_summary, sex_dif = tbl_dif_sex, age_dif = tbl_dif_age)))
}

# -----------------------------------------------------------------------------
# 5a. Cutoff calibration helpers (Empirical calibration vs PHQ>=10)
# -----------------------------------------------------------------------------
calc_binary_metrics <- function(theta, y, cutoff) {
  keep <- !is.na(theta) & !is.na(y)
  theta <- theta[keep]
  y <- as.integer(y[keep])
  yhat <- as.integer(theta >= cutoff)
  
  TP <- sum(yhat == 1 & y == 1)
  TN <- sum(yhat == 0 & y == 0)
  FP <- sum(yhat == 1 & y == 0)
  FN <- sum(yhat == 0 & y == 1)
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  acc  <- (TP + TN) / (TP + TN + FP + FN)
  ppv  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  npv  <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
  f1   <- if (!is.na(ppv) && !is.na(sens) && (ppv + sens) > 0) 2 * ppv * sens / (ppv + sens) else NA_real_
  
  list(TP = TP, TN = TN, FP = FP, FN = FN, sensitivity = sens, specificity = spec,
       accuracy = acc, ppv = ppv, npv = npv, f1 = f1)
}

calibrate_theta_cutoff <- function(theta, y,
                                   method = c("youden", "sens_at_spec", "spec_at_sens", "prevalence_match"),
                                   target_sensitivity = 0.80,
                                   target_specificity = 0.95,
                                   grid_n = 2000) {
  method <- match.arg(method)
  
  keep <- !is.na(theta) & !is.na(y)
  theta <- theta[keep]
  y <- as.integer(y[keep])
  
  if (length(unique(y)) < 2) {
    return(list(cutoff = NA_real_, roc = NULL, note = "Outcome has <2 classes; cannot calibrate."))
  }
  
  # Candidate cutoffs across observed theta range
  rng <- range(theta)
  grid <- seq(rng[1], rng[2], length.out = grid_n)
  
  roc <- lapply(grid, function(cut) {
    m <- calc_binary_metrics(theta, y, cut)
    data.frame(
      cutoff = cut,
      sensitivity = m$sensitivity,
      specificity = m$specificity,
      accuracy = m$accuracy,
      ppv = m$ppv,
      npv = m$npv,
      f1 = m$f1,
      stringsAsFactors = FALSE
    )
  })
  roc_df <- do.call(rbind, roc)
  
  # Choose cutoff
  if (method == "youden") {
    roc_df$youden <- roc_df$sensitivity + roc_df$specificity - 1
    best <- roc_df[which.max(roc_df$youden), ]
  } else if (method == "sens_at_spec") {
    ok <- roc_df[!is.na(roc_df$specificity) & roc_df$specificity >= target_specificity, ]
    if (nrow(ok) == 0) {
      # fallback: closest specificity
      roc_df$spec_gap <- abs(roc_df$specificity - target_specificity)
      best <- roc_df[which.min(roc_df$spec_gap), ]
    } else {
      best <- ok[which.max(ok$sensitivity), ]
    }
  } else if (method == "spec_at_sens") {
    ok <- roc_df[!is.na(roc_df$sensitivity) & roc_df$sensitivity >= target_sensitivity, ]
    if (nrow(ok) == 0) {
      # fallback: closest sensitivity
      roc_df$sens_gap <- abs(roc_df$sensitivity - target_sensitivity)
      best <- roc_df[which.min(roc_df$sens_gap), ]
    } else {
      best <- ok[which.max(ok$specificity), ]
    }
  } else if (method == "prevalence_match") {
    prev <- mean(y == 1)
    roc_df$prev_hat <- 1 - roc_df$specificity  # not prevalence; compute directly
    # better: compute predicted positive rate at each cutoff
    roc_df$pos_rate <- sapply(grid, function(cut) mean(theta >= cut))
    roc_df$gap <- abs(roc_df$pos_rate - prev)
    best <- roc_df[which.min(roc_df$gap), ]
  }
  
  list(cutoff = as.numeric(best$cutoff), roc = roc_df, note = NULL)
}

# -----------------------------------------------------------------------------
# 5. Linking, Scoring & Visualization Wrapper
#    FIX: Multi-dim equivalent SSQ sum should hold other dims at 0.
#    ADD: Optional empirical cutoff calibration to improve sensitivity.
# -----------------------------------------------------------------------------
perform_linking_and_scoring <- function(mod_joint, df_irt, ssq_items, phq_items,
                                        cutoff_method = c("model_tcc", "youden", "sens_at_spec", "spec_at_sens"),
                                        target_sensitivity = 0.80,
                                        target_specificity = 0.95) {
  
  cutoff_method <- match.arg(cutoff_method)
  
  joint_data <- df_irt[, c(phq_items, ssq_items)]
  valid_joint <- rowSums(!is.na(joint_data)) > 0
  
  # 1. Estimate Theta
  theta_joint <- fscores(mod_joint, method = "EAP", full.scores = TRUE)
  df_irt$Theta_Harmonized <- NA_real_
  df_irt$Theta_Harmonized[valid_joint] <- theta_joint[, 1]
  
  # 2. PHQ Raw binary (reference)
  df_irt$Sum_PHQ <- rowSums(df_irt[, phq_items, drop = FALSE], na.rm = TRUE)
  df_irt$Sum_SSQ <- rowSums(df_irt[, ssq_items, drop = FALSE], na.rm = TRUE)
  df_irt$Depression_PHQ_Raw <- ifelse(df_irt$Sum_PHQ >= 10, 1, 0)
  
  # 3. Model-based linking cutoff: Theta where expected PHQ sum == 10
  phq_idx <- 1:length(phq_items)
  n_dims <- mod_joint@Model$nfact
  t_seq <- seq(-4, 4, 0.01)
  
  if (n_dims == 1) {
    theta_grid_link <- matrix(t_seq)
  } else {
    theta_grid_link <- matrix(0, nrow = length(t_seq), ncol = n_dims)
    theta_grid_link[, 1] <- t_seq
  }
  
  tcc_phq <- expected.test(mod_joint, Theta = theta_grid_link, which.items = phq_idx)
  min_idx <- which.min(abs(tcc_phq - 10))
  latent_cutoff_tcc <- t_seq[min_idx]
  
  # 4. Optional empirical calibration against PHQ>=10
  roc_obj <- NULL
  latent_cutoff <- latent_cutoff_tcc
  
  if (cutoff_method != "model_tcc") {
    roc_obj <- calibrate_theta_cutoff(
      theta = df_irt$Theta_Harmonized[valid_joint],
      y = df_irt$Depression_PHQ_Raw[valid_joint],
      method = cutoff_method,
      target_sensitivity = target_sensitivity,
      target_specificity = target_specificity
    )
    if (!is.na(roc_obj$cutoff)) latent_cutoff <- roc_obj$cutoff
  }
  
  message(sprintf("Model-based cutoff (TCC PHQ=10): %.3f", latent_cutoff_tcc))
  message(sprintf("Chosen cutoff (%s): %.3f", cutoff_method, latent_cutoff))
  
  # 5. Classification
  df_irt$Depression_Harmonized <- ifelse(df_irt$Theta_Harmonized >= latent_cutoff, 1, 0)
  
  # 6. Metrics + confusion
  tbl <- table(
    Harmonized = factor(df_irt$Depression_Harmonized, levels = c(0, 1)),
    PHQ_Raw = factor(df_irt$Depression_PHQ_Raw, levels = c(0, 1))
  )
  
  TN <- tbl[1, 1]
  FN <- tbl[1, 2]
  FP <- tbl[2, 1]
  TP <- tbl[2, 2]
  total <- sum(tbl)
  
  # Equivalent SSQ Sum Score lookup (hold other dims = 0)
  ssq_idx_cols <- (length(phq_items) + 1):(length(phq_items) + length(ssq_items))
  
  theta_at_cut <- rep(0, n_dims)
  theta_at_cut[1] <- latent_cutoff
  
  equiv_sum <- expected.test(mod_joint, Theta = matrix(theta_at_cut, nrow = 1), which.items = ssq_idx_cols)
  
  metrics_df <- data.frame(
    Metric = c("Accuracy", "Sensitivity (Recall)", "Specificity", "Cohen's Kappa", "Equivalent SSQ-10 Sum"),
    Value = c((TP + TN) / total, TP / (TP + FN), TN / (TN + FP), cohen.kappa(tbl)$kappa, equiv_sum)
  )
  
  tbl_class <- gt(metrics_df) %>%
    tab_header(title = "Classification Agreement: Harmonized IRT vs PHQ-9 Raw Cutoff") %>%
    fmt_number(columns = "Value", decimals = 3) %>%
    fmt_percent(columns = "Value", rows = 1:3, decimals = 1)
  
  tbl_conf <- gt(as.data.frame(tbl)) %>%
    tab_header(title = "Confusion Matrix") %>%
    cols_label(Freq = "Count")
  
  # 7. Visualization (Test Info)
  p_seq <- seq(-4, 4, 0.05)
  if (n_dims == 1) {
    p_grid <- matrix(p_seq)
  } else {
    p_grid <- matrix(0, nrow = length(p_seq), ncol = n_dims)
    p_grid[, 1] <- p_seq
  }
  
  info_phq <- testinfo(mod_joint, p_grid, which.items = phq_idx)
  info_ssq <- testinfo(mod_joint, p_grid, which.items = ssq_idx_cols)
  
  p_info <- ggplot(
    data.frame(
      Theta = rep(p_seq, 2),
      Information = c(info_phq, info_ssq),
      Scale = rep(c("PHQ-9", "SSQ-10"), each = length(p_seq))
    ),
    aes(x = Theta, y = Information, color = Scale)
  ) +
    geom_line(size = 1) +
    geom_vline(xintercept = latent_cutoff, linetype = "dashed", color = "black") +
    annotate(
      "text",
      x = latent_cutoff + 0.5,
      y = max(c(info_phq, info_ssq)),
      label = paste0("Cut-off
(Theta = ", round(latent_cutoff, 2), ")"),
      family = "Arial",
      size = 3.5
    ) +
    labs(
      title = NULL,
      subtitle = "Test Information Functions",
      x = "Latent Depression (Theta)",
      y = "Information"
    ) +
    theme_minimal(base_family = "Arial", base_size = 14) +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = "bottom")
  
  # 8. Dist Plots
  p_sex <- if (!all(is.na(df_irt$SEX))) {
    ggplot(df_irt, aes(x = Theta_Harmonized, fill = SEX)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(title = "Latent Dist by Sex")
  } else NULL
  
  p_age <- if (!all(is.na(df_irt$AGEGRP))) {
    ggplot(df_irt, aes(x = Theta_Harmonized, fill = AGEGRP)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(title = "Latent Dist by Age")
  } else NULL
  
  return(list(
    scores = df_irt,
    threshold = latent_cutoff,
    threshold_model_tcc = latent_cutoff_tcc,
    cutoff_method = cutoff_method,
    roc = roc_obj,
    metrics = metrics_df,
    tables = list(classification = tbl_class, confusion = tbl_conf),
    plots = list(info = p_info, dist_sex = p_sex, dist_age = p_age)
  ))
}

# -----------------------------------------------------------------------------
# 6. Scoring Engine Export
#    FIX: Correct SSQ parameter extraction using item indices (not names)
#    NOTE: Best practice is to export the joint model itself for perfect scoring.
# -----------------------------------------------------------------------------
extract_scoring_engine <- function(mod_joint, ssq_items, phq_items, cutoff_theta, joint_model_def = NULL) {
  message("\nExtracting parameters for scoring engine...")
  
  vals <- mod2values(mod_joint)
  
  # We fitted joint model with items ordered as: c(phq_items, ssq_items)
  joint_names <- c(phq_items, ssq_items)
  ssq_idx <- match(ssq_items, joint_names)
  if (anyNA(ssq_idx)) stop("Could not map SSQ items to joint item order.")
  
  vals_ssq <- vals[vals$item %in% ssq_idx, , drop = FALSE]
  
  scoring_engine <- list(
    values_ssq = vals_ssq,
    ssq_items = ssq_items,
    phq_items = phq_items,
    threshold = cutoff_theta,
    itemtype = "graded",
    n_factors = mod_joint@Model$nfact,
    joint_model_def = joint_model_def
  )
  
  return(scoring_engine)
}

# -----------------------------------------------------------------------------
# 7. Score NEW datasets (SSQ only)
#    FIX: Score on the *original joint model* using response.pattern with PHQ set to NA
# -----------------------------------------------------------------------------
score_new_cohort <- function(new_data, mod_joint, ssq_items, phq_items, cutoff_theta) {
  
  missing_ssq <- setdiff(ssq_items, names(new_data))
  if (length(missing_ssq) > 0) stop(paste("New data missing SSQ items:", paste(missing_ssq, collapse = ", ")))
  
  # Build response.pattern with full joint item set
  rp <- as.data.frame(matrix(NA, nrow = nrow(new_data), ncol = length(c(phq_items, ssq_items))))
  colnames(rp) <- c(phq_items, ssq_items)
  
  rp[, ssq_items] <- new_data[, ssq_items]
  
  # Factor -> numeric (0-based)
  rp[, ssq_items] <- lapply(rp[, ssq_items, drop = FALSE], function(x) if (is.factor(x)) as.numeric(x) - 1 else x)
  
  keep <- rowSums(!is.na(rp[, ssq_items, drop = FALSE])) > 0
  if (any(!keep)) warning(paste("Removed", sum(!keep), "rows with all SSQ missing."))
  
  out_theta <- rep(NA_real_, nrow(new_data))
  out_bin <- rep(NA_integer_, nrow(new_data))
  
  th <- fscores(mod_joint, method = "EAP", response.pattern = rp[keep, , drop = FALSE], full.scores = TRUE)[, 1]
  out_theta[keep] <- th
  out_bin[keep] <- as.integer(th >= cutoff_theta)
  
  return(data.frame(Theta_Harmonized = out_theta, Depression_Binary = out_bin))
}

# ==============================================================================
# B. MAIN EXECUTION PIPELINE
# ==============================================================================

run_irt_harmonization_pipeline <- function(phq_model = 1, ssq_model = 1, joint_model = 1,
                                           cutoff_method = c("model_tcc", "youden", "sens_at_spec", "spec_at_sens"),
                                           target_sensitivity = 0.80,
                                           target_specificity = 0.95,
                                           run_sensitivity = FALSE) {
  
  message("=== STARTING IRT HARMONIZATION PIPELINE ===")
  
  # 1. Load Data
  prep <- load_and_clean_data(ssq_idx = c(1:3, 8:14), phq_idx = 1:9)
  out_dir <- here("statistical_analysis/output/objects")
  
  # 2. Fit Models (With flexible structure & Caching)
  fit_res <- run_fitting_workflow(prep$data, prep$ssq_items, prep$phq_items,
                                  phq_model, ssq_model, joint_model, out_dir)
  
  # 3. Invariance Checks (With flexible structure & Caching)
  inv_res <- run_invariance_workflow(prep$data, prep$ssq_items, prep$phq_items, joint_model, out_dir)
  
  # 4. Linking & Scoring
  link_res <- perform_linking_and_scoring(
    mod_joint = fit_res$models$joint,
    df_irt = prep$data,
    ssq_items = prep$ssq_items,
    phq_items = prep$phq_items,
    cutoff_method = cutoff_method,
    target_sensitivity = target_sensitivity,
    target_specificity = target_specificity
  )
  
  # 5. Extract & Save Scoring Engine
  # NOTE: For robust future scoring, also save the joint model object.
  scoring_engine <- extract_scoring_engine(
    mod_joint = fit_res$models$joint,
    ssq_items = prep$ssq_items,
    phq_items = prep$phq_items,
    cutoff_theta = link_res$threshold,
    joint_model_def = joint_model
  )
  
  saveRDS(scoring_engine, file = file.path(out_dir, "SSQ10_Harmonized_Scoring_Engine.rds"))
  saveRDS(fit_res$models$joint, file = file.path(out_dir, "Joint_IRT_Model_For_Scoring.rds"))
  message("Scoring Engine + Joint Model saved.")
  
  # 6. Sensitivity Analysis (Optional Scoring of External Cohort)
  sensitivity_scores <- NULL
  
  if (run_sensitivity) {
    message("\n--- Running Sensitivity Analysis (External Cohort Scoring) ---")
    dreams_path <- here("../../_private_use/data_management/data/co-luminate/adam/dt_dreams_multi_ssq.RData")
    
    if (file.exists(dreams_path)) {
      message(paste("Loading Dreams Cohort from:", dreams_path))
      load(dreams_path)
      
      if (exists("dt_dreams_multi_ssq")) {
        message("Object 'dt_dreams_multi_ssq' found. Scoring...")
        
        scored_dreams <- tryCatch({
          score_new_cohort(
            new_data = dt_dreams_multi_ssq,
            mod_joint = fit_res$models$joint,
            ssq_items = prep$ssq_items,
            phq_items = prep$phq_items,
            cutoff_theta = link_res$threshold
          )
        }, error = function(e) {
          message(paste("Error scoring Dreams cohort:", e$message))
          return(NULL)
        })
        
        if (!is.null(scored_dreams)) {
          if ("USUBJID" %in% names(dt_dreams_multi_ssq)) {
            sensitivity_scores <- cbind(dt_dreams_multi_ssq[, "USUBJID", drop = FALSE], scored_dreams)
          } else {
            sensitivity_scores <- scored_dreams
          }
          
          dreams_out_file <- file.path(out_dir, "Scored_Dreams_Cohort.rds")
          saveRDS(sensitivity_scores, file = dreams_out_file)
          message(paste("Dreams Cohort scored successfully. Results saved to:", dreams_out_file))
        }
      } else {
        warning("File loaded, but object 'dt_dreams_multi_ssq' was not found.")
      }
    } else {
      warning("Sensitivity Data file not found (Skipping step).")
    }
  }
  
  # 7. Save Final Object
  final_output <- list(
    input_structure = list(ssq_factors = ssq_model, phq_factors = phq_model, joint_factors = joint_model),
    models = fit_res$models,
    invariance_models = list(sex = inv_res$sex, age = inv_res$age),
    tables = c(fit_res$tables, inv_res$tables, link_res$tables),
    plots = link_res$plots,
    scores = link_res$scores,
    metrics = link_res$metrics,
    threshold = link_res$threshold,
    scoring_engine = scoring_engine,
    sensitivity_scores = sensitivity_scores
  )
  
  save_path <- file.path(out_dir, "03_IRT_Harmonization_Results.rds")
  saveRDS(final_output, file = save_path)
  
  message(paste("\nPipeline Complete. Results saved to:", save_path))
  return(final_output)
}

# Execute
if (sys.nframe() == 0) {
  # Example: Run with default 1-factor models and Sensitivity Check enabled
  res <- run_irt_harmonization_pipeline(
    cutoff_method = "model_tcc",
    run_sensitivity = TRUE
  )
}

res_2 <- run_irt_harmonization_pipeline(
  phq_model = 1,
  ssq_model = 1,
  joint_model = 1,
  cutoff_method = "sens_at_spec",
  target_specificity = 0.95,
  run_sensitivity = TRUE
)

res_cal <- run_irt_harmonization_pipeline(
  phq_model = 1,
  ssq_model = 1,
  joint_model = 1,
  cutoff_method = "youden",
  run_sensitivity = TRUE
)

# 
# # 1. Define the specific syntax for each scale
# phq_structure <- "F1 =~ PHQ906 + PHQ907 + PHQ908 + PHQ909\nF2 =~ PHQ901 + PHQ902 + PHQ903 + PHQ904 + PHQ905"
# ssq_structure <- "F1 =~ SSQ09 + SSQ10 + SSQ11 + SSQ12 + SSQ13 + SSQ14 + SSQ08\nF2 =~ SSQ01 + SSQ02 + SSQ03"
# joint_split_structure <- "
#   F1 = PHQ901, PHQ902, PHQ903, PHQ904, PHQ905, PHQ906, PHQ907, PHQ908, PHQ909
#   F2 = SSQ01, SSQ02, SSQ03, SSQ08, SSQ09, SSQ10, SSQ11, SSQ12, SSQ13, SSQ14
#   COV = F1*F2
# "
# # 2. Run the pipeline
# # We pass the custom structures for the individual model checks.
# # We leave 'joint_model = 1' because the goal is to link them onto a SINGLE severity scale.
# res_2 <- run_irt_harmonization_pipeline(
#   phq_model = 1,#phq_structure,
#   ssq_model = 1,#ssq_structure,
#   joint_model = 1,#joint_split_structure,
#   run_sensitivity = TRUE
# )

SUBJEXC_<-res_2[["scores"]] |> filter(Depression_PHQ_Raw==1 & Depression_Harmonized==0) |> pull(USUBJID)
View(
  res_2[["scores"]] |> 
    select(USUBJID,SEX,AGEGRP,
             Sum_PHQ,Sum_SSQ,
             Theta_Harmonized,Depression_PHQ_Raw,Depression_Harmonized) |>  
    arrange(desc(Depression_Harmonized),desc(Depression_PHQ_Raw))
  )
View(
  res_cal[["scores"]] |> 
    select(USUBJID,SEX,AGEGRP,
           Sum_PHQ,Sum_SSQ,
           Theta_Harmonized,Depression_PHQ_Raw,Depression_Harmonized) |>  
    arrange(desc(Depression_Harmonized),desc(Depression_PHQ_Raw))
)
res_2[["metrics"]]
res_cal[["metrics"]]
