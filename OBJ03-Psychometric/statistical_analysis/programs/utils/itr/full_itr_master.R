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
# -----------------------------------------------------------------------------
lavaan_to_mirt <- function(syntax) {
  # If numeric (e.g. 1 or 2 factors exploratory), return as is
  if (is.numeric(syntax) || is.null(syntax)) return(syntax)
  
  # Basic cleaning
  # Replace "=~" with "="
  s <- gsub("=~", "=", syntax)
  # Replace "+" with ","
  s <- gsub("\\+", ",", s)
  
  # Logic to ensure factors are correlated (Lavaan default)
  # Mirt requires explicit COV statement for multidimensional models defined by syntax
  
  # 1. Split into lines to identify factors
  lines <- unlist(strsplit(s, "\n"))
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  
  # 2. Extract factor names (Left hand side of =)
  factor_lines <- grep("=", lines)
  factor_names <- sapply(strsplit(lines[factor_lines], "="), function(x) trimws(x[1]))
  factor_names <- unique(factor_names)
  
  # 3. If > 1 factor, append COV statement
  if (length(factor_names) > 1) {
    cov_stmt <- paste0("COV = ", paste(factor_names, collapse = "*"))
    s <- paste(s, cov_stmt, sep = "\n")
  }
  
  # Convert to mirt.model object
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
    dt_main <- dt_psychometric
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
  
  if(length(valid_cols) == 0) return(NULL)
  
  fit_df <- as.data.frame(fit_obj)[, valid_cols, drop=FALSE]
  
  gt(fit_df) %>%
    tab_header(title = title_str) %>%
    fmt_number(columns = any_of(c("M2", "RMSEA", "RMSEA_5", "RMSEA_95", "SRMSR", "TLI", "CFI")), decimals = 3) %>%
    fmt_number(columns = any_of("p"), decimals = 4) %>%
    fmt_number(columns = any_of("df"), decimals = 0) %>%
    cols_label(M2 = "M2 Statistic", df = "df", p = "p-value", RMSEA_5 = "RMSEA (5%)", RMSEA_95 = "RMSEA (95%)") %>%
    tab_source_note(source_note = "Good Fit Criteria: CFI > 0.95, TLI > 0.95, RMSEA < 0.06, SRMSR < 0.08.")
}

# -----------------------------------------------------------------------------
# 3. Model Fitting Wrapper (Updated)
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
  mod_phq <- mirt(phq_dat, mirt_phq, itemtype = 'graded', verbose = FALSE)
  m2_phq <- M2(mod_phq, type = "C2", calcNull = FALSE)
  
  # B. SSQ-10
  message("Fitting GRM to SSQ-10...")
  ssq_dat <- df_irt[, ssq_items]
  ssq_dat <- ssq_dat[rowSums(!is.na(ssq_dat)) > 0, ]
  mod_ssq <- mirt(ssq_dat, mirt_ssq, itemtype = 'graded', verbose = FALSE)
  m2_ssq <- M2(mod_ssq, type = "C2", calcNull = FALSE)
  
  # C. Joint Model
  message("Fitting Joint GRM (Harmonization)...")
  joint_dat <- df_irt[, c(phq_items, ssq_items)]
  joint_dat <- joint_dat[rowSums(!is.na(joint_dat)) > 0, ]
  mod_joint <- mirt(joint_dat, mirt_joint, itemtype = 'graded', verbose = FALSE)
  
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
  
  # Determine label for filename
  # Since models can be large strings, we create a hash or simplified label
  if (is.numeric(joint_model)) {
    model_label <- paste0(joint_model, "factor")
  } else {
    model_label <- "custom_structure" 
  }
  
  fit_file_name <- paste0("03_IRT_fitted_models_", model_label, ".rds")
  fit_file_path <- file.path(output_dir, fit_file_name)
  run_checks <- TRUE
  
  # Define current structure signature
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
# -----------------------------------------------------------------------------
check_group_invariance <- function(item_data, group_vec, group_name, model_def) {
  if (length(unique(na.omit(group_vec))) < 2) return(NULL)
  
  message(paste0("Checking Invariance by ", group_name, "..."))
  dat_grp <- data.frame(item_data, Group = as.factor(group_vec)) 
  dat_grp <- na.omit(dat_grp) 
  
  # Convert syntax if necessary
  mirt_model <- lavaan_to_mirt(model_def)
  
  mod_mg <- multipleGroup(dat_grp[, -ncol(dat_grp)], 
                          model = mirt_model,
                          group = dat_grp$Group, 
                          invariance = c('slopes', 'intercepts', 'free_means', 'free_var'),
                          verbose = FALSE)
  
  coefs <- coef(mod_mg, simplify = TRUE)
  dif_res <- DIF(mod_mg, which.par = c('a1', 'd'), scheme = 'drop', p.adjust = 'fdr')
  
  return(list(model = mod_mg, dif = dif_res, latent_pars = coefs))
}

summarize_invariance <- function(res_obj, group_label) {
  if(is.null(res_obj)) return(NULL)
  if (!is.list(res_obj$latent_pars) || length(res_obj$latent_pars) < 2) return(NULL)
  
  grp_names <- names(res_obj$latent_pars)
  ref_grp <- if(!is.null(grp_names)) grp_names[1] else "Group 1"
  foc_grp <- if(!is.null(grp_names)) grp_names[2] else "Group 2"
  
  focal_pars <- res_obj$latent_pars[[2]]$GroupPars
  latent_mean <- NA; latent_var <- NA
  
  if (!is.null(focal_pars) && length(focal_pars) > 0) {
    if ("MEAN_1" %in% names(focal_pars)) latent_mean <- focal_pars["MEAN_1"] else if (length(focal_pars) >= 1) latent_mean <- focal_pars[1]
    if ("COV_11" %in% names(focal_pars)) latent_var <- focal_pars["COV_11"] else if (length(focal_pars) >= 2) latent_var <- focal_pars[2]
  }
  
  sig_items <- rownames(res_obj$dif)[res_obj$dif$p < 0.05]
  if(is.null(sig_items)) sig_items <- character(0)
  
  data.frame(
    Grouping = group_label,
    Reference_Group = ref_grp,
    Focal_Group = foc_grp,
    Latent_Mean_Diff = round(as.numeric(latent_mean), 3),
    Latent_Var_Ratio = round(as.numeric(latent_var), 3),
    DIF_Items = if(length(sig_items) > 0) paste(sig_items, collapse=", ") else "None",
    DIF_Count = length(sig_items),
    stringsAsFactors = FALSE
  )
}

create_dif_table <- function(res_obj, group_label) {
  if(is.null(res_obj)) return(NULL)
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
    tab_style(style = list(cell_fill(color = "#F9E3E3"), cell_text(weight = "bold")), locations = cells_body(rows = adj_p < 0.05))
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
    clean_dat <- joint_data[valid_joint, ]
    
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
    tab_source_note(source_note = "Impact: Latent Mean Diff > 0 indicates Focal group has higher depression. DIF: p < 0.05 (FDR).")
  
  tbl_dif_sex <- create_dif_table(res_sex, "Biological Sex")
  tbl_dif_age <- create_dif_table(res_age, "Age Group")
  
  return(list(sex = res_sex, age = res_age, tables = list(summary = tbl_summary, sex_dif = tbl_dif_sex, age_dif = tbl_dif_age)))
}

# -----------------------------------------------------------------------------
# 5. Linking, Scoring & Visualization Wrapper
# -----------------------------------------------------------------------------
perform_linking_and_scoring <- function(mod_joint, df_irt, ssq_items, phq_items) {
  
  joint_data <- df_irt[, c(phq_items, ssq_items)]
  valid_joint <- rowSums(!is.na(joint_data)) > 0
  
  # 1. Estimate Theta
  theta_joint <- fscores(mod_joint, method = "EAP", full.scores = TRUE)
  df_irt$Theta_Harmonized <- NA
  df_irt$Theta_Harmonized[valid_joint] <- theta_joint[, 1]
  
  # 2. Link Theta to PHQ Cutoff (10)
  phq_idx <- 1:length(phq_items)
  
  # Detect dimensionality from model object
  n_dims <- mod_joint@Model$nfact
  t_seq <- seq(-4, 4, 0.01)
  
  # Build Theta Grid for Linking (varying primary dim, holding others at 0)
  if (n_dims == 1) {
    theta_grid_link <- matrix(t_seq)
  } else {
    theta_grid_link <- matrix(0, nrow=length(t_seq), ncol=n_dims)
    theta_grid_link[,1] <- t_seq
  }
  
  tcc_phq <- expected.test(mod_joint, Theta = theta_grid_link, which.items = phq_idx)
  min_idx <- which.min(abs(tcc_phq - 10))
  latent_cutoff <- t_seq[min_idx]
  message(sprintf("Harmonized Latent Threshold (Theta for PHQ=10): %.3f", latent_cutoff))
  
  # 3. Classification
  df_irt$Sum_PHQ <- rowSums(df_irt[, phq_items], na.rm = TRUE)
  df_irt$Sum_SSQ <- rowSums(df_irt[, ssq_items], na.rm = TRUE)
  df_irt$Depression_Harmonized <- ifelse(df_irt$Theta_Harmonized >= latent_cutoff, 1, 0)
  df_irt$Depression_PHQ_Raw <- ifelse(df_irt$Sum_PHQ >= 10, 1, 0)
  
  # 4. Metrics
  tbl <- table(Harmonized = factor(df_irt$Depression_Harmonized, levels=c(0,1)), 
               PHQ_Raw = factor(df_irt$Depression_PHQ_Raw, levels=c(0,1)))
  
  TN <- tbl[1,1]; FN <- tbl[1,2]; FP <- tbl[2,1]; TP <- tbl[2,2]; total <- sum(tbl)
  
  # Equivalent SSQ Sum Score lookup
  ssq_idx_cols <- (length(phq_items)+1):(length(phq_items)+length(ssq_items))
  equiv_sum <- expected.test(mod_joint, Theta = matrix(latent_cutoff, ncol=n_dims), which.items = ssq_idx_cols)
  
  metrics_df <- data.frame(
    Metric = c("Accuracy", "Sensitivity (Recall)", "Specificity", "Cohen's Kappa", "Equivalent SSQ-10 Sum"),
    Value = c((TP + TN)/total, TP/(TP + FN), TN/(TN + FP), cohen.kappa(tbl)$kappa, equiv_sum)
  )
  
  tbl_class <- gt(metrics_df) %>%
    tab_header(title = "Classification Agreement: Harmonized IRT vs PHQ-9 Raw Cutoff") %>%
    fmt_number(columns = "Value", decimals = 3) %>%
    fmt_percent(columns = "Value", rows = 1:3, decimals = 1)
  
  tbl_conf <- gt(as.data.frame(tbl)) %>%
    tab_header(title = "Confusion Matrix") %>% cols_label(Freq = "Count")
  
  # 5. Visualization (Test Info)
  # Plotting Grid (coarser for speed)
  p_seq <- seq(-4, 4, 0.05)
  if (n_dims == 1) {
    p_grid <- matrix(p_seq)
  } else {
    p_grid <- matrix(0, nrow=length(p_seq), ncol=n_dims)
    p_grid[,1] <- p_seq
  }
  
  info_phq <- testinfo(mod_joint, p_grid, which.items = phq_idx)
  info_ssq <- testinfo(mod_joint, p_grid, which.items = ssq_idx_cols)
  
  p_info <- ggplot(data.frame(Theta = rep(p_seq, 2), Information = c(info_phq, info_ssq), Scale = rep(c("PHQ-9", "SSQ-10"), each = length(p_seq))), 
                   aes(x = Theta, y = Information, color = Scale)) +
    geom_line(size = 1) +
    geom_vline(xintercept = latent_cutoff, linetype = "dashed", color = "black") +
    annotate("text", x = latent_cutoff + 0.5, y = max(c(info_phq, info_ssq)), 
             label = paste0("Cut-off\n(Theta = ", round(latent_cutoff, 2), ")"), family = "Arial", size = 3.5) +
    labs(title = NULL, subtitle = "Test Information Functions", x = "Latent Depression (Theta)", y = "Information") +
    theme_minimal(base_family = "Arial", base_size = 14) + scale_color_brewer(palette = "Set1") + theme(legend.position = "bottom")
  
  # 6. Dist Plots
  p_sex <- if(!all(is.na(df_irt$SEX))) ggplot(df_irt, aes(x=Theta_Harmonized, fill=SEX)) + geom_density(alpha=0.5) + theme_minimal() + labs(title="Latent Dist by Sex") else NULL
  p_age <- if(!all(is.na(df_irt$AGEGRP))) ggplot(df_irt, aes(x=Theta_Harmonized, fill=AGEGRP)) + geom_density(alpha=0.5) + theme_minimal() + labs(title="Latent Dist by Age") else NULL
  
  return(list(
    scores = df_irt,
    threshold = latent_cutoff,
    metrics = metrics_df,
    tables = list(classification = tbl_class, confusion = tbl_conf),
    plots = list(info = p_info, dist_sex = p_sex, dist_age = p_age)
  ))
}

# -----------------------------------------------------------------------------
# 6. Scoring Engine Export
# -----------------------------------------------------------------------------
extract_scoring_engine <- function(mod_joint, ssq_items, cutoff_theta) {
  message("\nExtracting SSQ-10 Parameters for Scoring Engine...")
  pars <- mod2values(mod_joint)
  pars_ssq <- pars[pars$item %in% ssq_items, ]
  nfact <- mod_joint@Model$nfact
  
  # Also store the model syntax/structure if possible
  # For now, storing parameters + nfact is sufficient for reconstruction
  
  scoring_engine <- list(
    parameters = pars_ssq,
    items = ssq_items,
    threshold = cutoff_theta,
    model_type = 'graded', 
    n_factors = nfact      
  )
  return(scoring_engine)
}

score_new_cohort <- function(new_data, scoring_engine) {
  missing <- setdiff(scoring_engine$items, names(new_data))
  if (length(missing) > 0) stop(paste("New data missing items:", paste(missing, collapse=", ")))
  
  dat_score <- new_data[, scoring_engine$items]
  dat_score[] <- lapply(dat_score, function(x) if(is.factor(x)) as.numeric(x)-1 else x)
  
  keep_rows <- rowSums(!is.na(dat_score)) > 0
  if (sum(!keep_rows) > 0) warning(paste("Removed", sum(!keep_rows), "rows with all missing data."))
  dat_score_clean <- dat_score[keep_rows, , drop = FALSE]
  
  message("Estimating Theta for new cohort using fixed harmonization parameters...")
  
  n_items <- ncol(dat_score_clean)
  item_types_vec <- rep(scoring_engine$model_type, n_items)
  n_factors <- if(!is.null(scoring_engine$n_factors)) scoring_engine$n_factors else 1
  
  template_pars <- mirt(dat_score_clean, n_factors, itemtype = item_types_vec, pars = 'values', verbose = FALSE)
  saved_pars <- scoring_engine$parameters
  
  for(i in 1:nrow(template_pars)) {
    itm <- template_pars$item[i]
    pnm <- template_pars$name[i]
    if (itm == "GROUP") next
    match_row <- saved_pars[saved_pars$item == itm & saved_pars$name == pnm, ]
    if (nrow(match_row) == 1) {
      template_pars$value[i] <- match_row$value
      template_pars$est[i]   <- FALSE 
    }
  }
  
  mod_fixed <- mirt(dat_score_clean, n_factors, itemtype = item_types_vec, 
                    pars = template_pars, verbose = FALSE, calcNull = FALSE)
  
  scores_clean <- fscores(mod_fixed, method = "EAP", full.scores = TRUE)
  classifications_clean <- ifelse(scores_clean[,1] >= scoring_engine$threshold, 1, 0)
  
  final_scores <- rep(NA, nrow(new_data))
  final_class  <- rep(NA, nrow(new_data))
  final_scores[keep_rows] <- scores_clean[,1]
  final_class[keep_rows]  <- classifications_clean
  
  return(data.frame(Theta_Harmonized = final_scores, Depression_Binary = final_class))
}


# ==============================================================================
# B. MAIN EXECUTION PIPELINE
# ==============================================================================

run_irt_harmonization_pipeline <- function(phq_model = 1, ssq_model = 1, joint_model = 1) {
  
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
  link_res <- perform_linking_and_scoring(fit_res$models$joint, prep$data, prep$ssq_items, prep$phq_items)
  
  # 5. Extract & Save Scoring Engine
  scoring_engine <- extract_scoring_engine(fit_res$models$joint, prep$ssq_items, link_res$threshold)
  saveRDS(scoring_engine, file = file.path(out_dir, "SSQ10_Harmonized_Scoring_Engine.rds"))
  message("Scoring Engine saved.")
  
  # 6. Save Final Object
  final_output <- list(
    input_structure = list(ssq_factors = ssq_model, phq_factors = phq_model, joint_factors = joint_model),
    models = fit_res$models,
    invariance_models = list(sex = inv_res$sex, age = inv_res$age),
    tables = c(fit_res$tables, inv_res$tables, link_res$tables),
    plots = link_res$plots,
    scores = link_res$scores,
    metrics = link_res$metrics,
    threshold = link_res$threshold,
    scoring_engine = scoring_engine 
  )
  
  save_path <- file.path(out_dir, "03_IRT_Harmonization_Results.rds")
  saveRDS(final_output, file = save_path)
  
  message(paste("\nPipeline Complete. Results saved to:", save_path))
  return(final_output) 
}

# Execute
if (sys.nframe() == 0) { 
  # --- EXAMPLE USAGE WITH SYNTAX ---
  # If you have consensus syntax from Script 02, paste it here.
  # Otherwise, defaults to '1' factor.
  # res <- run_irt_harmonization_pipeline(joint_model = "F1 = 1-10 \n F2 = 11-19")
  
  res <- run_irt_harmonization_pipeline()
  
  # Scoring Example (Commented out to prevent auto-run on source)
  # dreams_path <- here("../../_private_use/data_management/data/co-luminate/adam/dt_dreams_multi_ssq.RData")
  # if (file.exists(dreams_path)) { load(dreams_path); score_new_cohort(dt_dreams_multi_ssq, res$scoring_engine) }
}