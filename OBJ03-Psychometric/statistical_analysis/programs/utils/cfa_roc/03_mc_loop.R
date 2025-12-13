############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 03_mc_loop.R - Monte Carlo CFA loop
############################################################

# Helper: compute reliability (alpha, omega) from a lavaan fit --------
# Uses standardized loadings (est.std) and assumes a simple structure
# with residuals uncorrelated. For each latent factor, we compute:
#   Var(T_factor) = (sum(lambda))^2 + sum(theta)
#   alpha         = (k/(k-1)) * (1 - k / Var(T_factor))
#   omega_total   = (sum(lambda))^2 / Var(T_factor)
# where theta = 1 - lambda^2 for standardized solution.
compute_reliability_from_fit <- function(fit) {
  est_std <- standardizedSolution(fit)
  est_std <- est_std[est_std$op == "=~", , drop = FALSE]
  if (nrow(est_std) == 0) return(NULL)
  
  out_list <- list()
  factors <- unique(est_std$lhs)
  
  for (f in factors) {
    lambda <- est_std$est.std[est_std$lhs == f]
    k <- length(lambda)
    if (k < 2) next
    
    theta <- 1 - lambda^2
    sum_lambda <- sum(lambda)
    sum_theta  <- sum(theta)
    var_T      <- sum_lambda^2 + sum_theta
    
    if (!is.finite(var_T) || var_T <= 0) next
    
    alpha <- (k / (k - 1)) * (1 - k / var_T)
    omega <- (sum_lambda^2) / var_T
    
    out_list[[f]] <- data.frame(
      Scale = f,
      alpha = alpha,
      omega = omega,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(out_list) == 0) return(NULL)
  do.call(rbind, out_list)
}


get_roc_metrics <- function(data, seed, group_label, items) {
  if (nrow(data) < 20 || length(unique(data$Depressed)) < 2) return(NULL)
  tryCatch({
    # Calculate Sum Score from SSQ/PHQ items
    item_data <- data[, items, drop = FALSE]
    item_data[] <- lapply(item_data, function(x) if (is.factor(x)) as.numeric(x) - 1 else x)
    sum_score <- rowSums(item_data, na.rm = TRUE)
    
    roc_obj <- roc(data$Depressed, sum_score, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    coords_res <- coords(
      roc_obj,
      x = "best",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    )
    
    if (nrow(coords_res) > 1) coords_res <- coords_res[1, ]
    
    # Calculate Cohen's Kappa
    thresh <- coords_res$threshold
    pred   <- ifelse(sum_score >= thresh, 1, 0)
    obs    <- data$Depressed
    
    tbl <- table(factor(pred, levels = c(0, 1)), factor(obs, levels = c(0, 1)))
    
    po <- sum(diag(tbl)) / sum(tbl)
    pe <- ((sum(tbl[2, ]) * sum(tbl[, 2])) + (sum(tbl[1, ]) * sum(tbl[, 1]))) / sum(tbl)^2
    kappa_val <- (po - pe) / (1 - pe)
    
    # Prevalence according to PHQ-09 (gold standard) and SSQ-10 (test) at chosen cutoff
    prevalence_phq <- mean(obs == 1, na.rm = TRUE)
    prevalence_ssq <- mean(pred == 1, na.rm = TRUE)
    
    data.frame(
      Seed           = seed,
      Group          = group_label,
      AUC            = auc_val,
      Optimal_Cutoff = coords_res$threshold,
      Sensitivity    = coords_res$sensitivity,
      Specificity    = coords_res$specificity,
      Kappa          = kappa_val,
      Prevalence_PHQ = prevalence_phq,
      Prevalence_SSQ = prevalence_ssq
    )
  }, error = function(e) NULL)
}

get_roc_curve_data <- function(data, seed, group_label, items) {
  if (nrow(data) < 20 || length(unique(data$Depressed)) < 2) return(NULL)
  tryCatch({
    item_data <- data[, items, drop = FALSE]
    item_data[] <- lapply(item_data, function(x) if (is.factor(x)) as.numeric(x) - 1 else x)
    sum_score <- rowSums(item_data, na.rm = TRUE)
    
    roc_obj <- roc(data$Depressed, sum_score, quiet = TRUE)
    
    # Interpolate sensitivities at fixed specificities for averaging
    specs <- seq(0, 1, 0.01)
    senss <- coords(
      roc_obj,
      x    = specs,
      input = "specificity",
      ret   = "sensitivity",
      transpose = FALSE
    )
    
    data.frame(
      Seed        = seed,
      Group       = group_label,
      Specificity = specs,
      Sensitivity = as.numeric(unlist(senss))
    )
  }, error = function(e) NULL)
}

run_monte_carlo_cfa <- function(dt_main, cfa_model_syntax, target_items,
                                target_set_name, run_seeds, n_mc_samples) {
  res_fit         <- data.frame()
  res_loadings    <- data.frame()
  res_invariance  <- data.frame()
  res_roc         <- data.frame()
  res_roc_curves  <- data.frame()
  res_factor_corr <- data.frame()
  res_reliability <- data.frame()
  
  # Diagnostics log for each seed
  diagnostics <- data.frame(
    Seed   = integer(),
    Status = character(),
    stringsAsFactors = FALSE
  )
  
  valid_counter <- 0
  
  for (i in seq_along(run_seeds)) {
    if (valid_counter >= n_mc_samples) break
    
    curr_seed <- run_seeds[i]
    status <- "initial"
    
    set.seed(curr_seed)
    n_rows <- nrow(dt_main)
    idx_explore <- sample(seq_len(n_rows), floor(0.5 * n_rows), replace = FALSE)
    idx_confirm <- setdiff(seq_len(n_rows), idx_explore)
    dt_cfa <- dt_main[idx_confirm, ]
    
    if ("PHQSCR" %in% names(dt_cfa)) {
      dt_cfa$PHQBIN <- ifelse(dt_cfa$PHQSCR >= 10, "Depressed", "Non-Depressed")
      dt_cfa$PHQBIN <- factor(dt_cfa$PHQBIN,
                              levels = c("Non-Depressed", "Depressed"))
    }
    
    if ("AGE" %in% names(dt_cfa)) {
      dt_cfa$AGEGRP <- ifelse(dt_cfa$AGE < 20, "17-19", "20-24")
    }
    
    dt_cfa_num <- dt_cfa
    
    fit <- tryCatch({
      cfa(cfa_model_syntax, data = dt_cfa_num,
          estimator = "WLSMV", ordered = target_items)
    }, error = function(e) {
      status <<- "fit_error"
      NULL
    })
    
    is_valid_model <- FALSE
    if (!is.null(fit)) {
      if (!inspect(fit, "converged")) {
        status <- "no_convergence"
      } else if (!lavInspect(fit, "post.check")) {
        status <- "failed_post_check"
      } else {
        is_valid_model <- TRUE
      }
    } else if (status == "initial") {
      status <- "fit_null"
    }
    
    if (is_valid_model) {
      status <- "valid"
      valid_counter <- valid_counter + 1
      
      fits <- fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr", "chisq", "df"))
      res_fit <- rbind(res_fit, data.frame(Seed = curr_seed, as.list(fits)))
      
      est <- standardizedSolution(fit) %>%
        filter(op == "=~") %>%
        select(lhs, rhs, est.std, se, pvalue) %>%
        mutate(Seed = curr_seed)
      res_loadings <- rbind(res_loadings, est)
      
      # Latent factor correlations (only if 2+ latent variables)
      lv_names <- lavNames(fit, type = "lv")
      if (length(lv_names) >= 2) {
        phi <- lavInspect(fit, "cor.lv")
        phi_df <- as.data.frame(phi)
        # Store rownames in a separate column that does not clash with latent names
        phi_df$LV1 <- rownames(phi_df)
        
        phi_long <- phi_df %>%
          pivot_longer(
            cols      = all_of(lv_names),
            names_to  = "LV2",
            values_to = "Corr"
          ) %>%
          rename(F1 = LV1, F2 = LV2) %>%
          filter(F1 < F2) %>%
          mutate(Seed = curr_seed)
        
        res_factor_corr <- rbind(res_factor_corr, phi_long)
      }
      
      # Internal consistency reliability (Cronbach's alpha, McDonald's omega)
      rel_df <- compute_reliability_from_fit(fit)
      if (!is.null(rel_df) && nrow(rel_df) > 0) {
        rel_df$Seed <- curr_seed
        res_reliability <- rbind(res_reliability, rel_df)
      }
      
      
      if ("PHQBIN" %in% names(dt_cfa_num)) {
        inv_phq <- get_invariance_deltas(dt_cfa_num, cfa_model_syntax, "PHQBIN", target_items)
        if (!is.null(inv_phq))
          res_invariance <- rbind(res_invariance, cbind(Seed = curr_seed, inv_phq))
      }
      
      inv_sex <- get_invariance_deltas(dt_cfa_num, cfa_model_syntax, "SEX", target_items)
      if (!is.null(inv_sex))
        res_invariance <- rbind(res_invariance, cbind(Seed = curr_seed, inv_sex))
      
      inv_age <- get_invariance_deltas(dt_cfa_num, cfa_model_syntax, "AGEGRP", target_items)
      if (!is.null(inv_age))
        res_invariance <- rbind(res_invariance, cbind(Seed = curr_seed, inv_age))
      
      if ("PHQSCR" %in% names(dt_cfa)) {
        roc_dat <- dt_cfa
        roc_dat$Depressed <- ifelse(roc_dat$PHQSCR >= 10, 1, 0)
        
        if ("AGE" %in% names(roc_dat)) {
          roc_dat$AGEGRP <- ifelse(roc_dat$AGE < 20, "17-19", "20-24")
        }
        
        # Secondary restricted ROC: exclude cases with all target items = 0 AND PHQSCR = 0
        item_numeric <- roc_dat[, target_items, drop = FALSE]
        item_numeric[] <- lapply(item_numeric, function(x) {
          if (is.factor(x)) as.numeric(x) - 1 else x
        })
        sum_items <- rowSums(item_numeric, na.rm = TRUE)
        zero_items <- sum_items == 0
        zero_phq   <- roc_dat$PHQSCR == 0
        keep_restricted <- !(zero_items & zero_phq)
        roc_dat_restricted <- roc_dat[keep_restricted, , drop = FALSE]
        
        # Primary ROC on full sample
        res_roc <- rbind(res_roc,
                         get_roc_metrics(roc_dat, curr_seed, "Overall", target_items))
        res_roc_curves <- rbind(res_roc_curves,
                                get_roc_curve_data(roc_dat, curr_seed, "Overall", target_items))
        
        # Secondary ROC: restricted sample (exclude 0-0)
        res_roc <- rbind(res_roc,
                         get_roc_metrics(roc_dat_restricted, curr_seed,
                                         "Overall (Excl. 0-0)", target_items))
        res_roc_curves <- rbind(res_roc_curves,
                                get_roc_curve_data(roc_dat_restricted, curr_seed,
                                                   "Overall (Excl. 0-0)", target_items))
        
        if ("SEX" %in% names(roc_dat)) {
          for (g in unique(na.omit(roc_dat$SEX))) {
            sub_dat <- roc_dat[roc_dat$SEX == g, ]
            res_roc <- rbind(res_roc,
                             get_roc_metrics(sub_dat, curr_seed, paste0("Sex: ", g), target_items))
            res_roc_curves <- rbind(res_roc_curves,
                                    get_roc_curve_data(sub_dat, curr_seed, paste0("Sex: ", g), target_items))
            
            # Restricted (exclude 0-0)
            sub_dat_restricted <- roc_dat_restricted[roc_dat_restricted$SEX == g, , drop = FALSE]
            res_roc <- rbind(res_roc,
                             get_roc_metrics(sub_dat_restricted, curr_seed,
                                             paste0("Sex: ", g, " (Excl. 0-0)"), target_items))
            res_roc_curves <- rbind(res_roc_curves,
                                    get_roc_curve_data(sub_dat_restricted, curr_seed,
                                                       paste0("Sex: ", g, " (Excl. 0-0)"), target_items))
          }
        }
        
        if ("AGEGRP" %in% names(roc_dat)) {
          for (g in unique(na.omit(roc_dat$AGEGRP))) {
            sub_dat <- roc_dat[roc_dat$AGEGRP == g, ]
            res_roc <- rbind(res_roc,
                             get_roc_metrics(sub_dat, curr_seed, paste0("Age: ", g), target_items))
            res_roc_curves <- rbind(res_roc_curves,
                                    get_roc_curve_data(sub_dat, curr_seed, paste0("Age: ", g), target_items))
            
            # Restricted (exclude 0-0)
            sub_dat_restricted <- roc_dat_restricted[roc_dat_restricted$AGEGRP == g, , drop = FALSE]
            res_roc <- rbind(res_roc,
                             get_roc_metrics(sub_dat_restricted, curr_seed,
                                             paste0("Age: ", g, " (Excl. 0-0)"), target_items))
            res_roc_curves <- rbind(res_roc_curves,
                                    get_roc_curve_data(sub_dat_restricted, curr_seed,
                                                       paste0("Age: ", g, " (Excl. 0-0)"), target_items))
          }
        }
      }
    }
    
    diagnostics <- rbind(
      diagnostics,
      data.frame(Seed = curr_seed, Status = status, stringsAsFactors = FALSE)
    )
    
  }
  
  
  list(
    res_fit         = res_fit,
    res_loadings    = res_loadings,
    res_invariance  = res_invariance,
    res_roc         = res_roc,
    res_roc_curves  = res_roc_curves,
    res_factor_corr = res_factor_corr,
    res_reliability = res_reliability,
    diagnostics     = diagnostics,
    n_valid         = valid_counter,
    n_tried         = length(run_seeds)
  )
}

# Alternative: full-seed Monte Carlo loop (no early stopping) ---------
# Evaluates *all* seeds in run_seeds, computing validity once and for all.
# This is useful if you want a fixed pool of evaluated seeds and then
# later choose how many valid seeds (n_mc_samples) to use in summaries
# without refitting the models.
run_monte_carlo_cfa_full <- function(dt_main, cfa_model_syntax, target_items,
                                     target_set_name, run_seeds) {
  res_fit         <- data.frame()
  res_loadings    <- data.frame()
  res_invariance  <- data.frame()
  res_roc         <- data.frame()
  res_roc_curves  <- data.frame()
  res_factor_corr <- data.frame()
  res_reliability <- data.frame()
  
  diagnostics <- data.frame(
    Seed   = integer(),
    Status = character(),
    stringsAsFactors = FALSE
  )
  
  valid_counter <- 0
  
  for (i in seq_along(run_seeds)) {
    curr_seed <- run_seeds[i]
    status <- "initial"
    
    set.seed(curr_seed)
    n_rows <- nrow(dt_main)
    idx_explore <- sample(seq_len(n_rows), floor(0.5 * n_rows), replace = FALSE)
    idx_confirm <- setdiff(seq_len(n_rows), idx_explore)
    dt_cfa <- dt_main[idx_confirm, ]
    
    if ("PHQSCR" %in% names(dt_cfa)) {
      dt_cfa$PHQBIN <- ifelse(dt_cfa$PHQSCR >= 10, "Depressed", "Non-Depressed")
      dt_cfa$PHQBIN <- factor(dt_cfa$PHQBIN,
                              levels = c("Non-Depressed", "Depressed"))
    }
    
    if ("AGE" %in% names(dt_cfa)) {
      dt_cfa$AGEGRP <- ifelse(dt_cfa$AGE < 20, "17-19", "20-24")
    }
    
    dt_cfa_num <- dt_cfa
    
    fit <- tryCatch({
      cfa(cfa_model_syntax, data = dt_cfa_num,
          estimator = "WLSMV", ordered = target_items)
    }, error = function(e) {
      status <<- "fit_error"
      NULL
    })
    
    is_valid_model <- FALSE
    if (!is.null(fit)) {
      if (!inspect(fit, "converged")) {
        status <- "no_convergence"
      } else if (!lavInspect(fit, "post.check")) {
        status <- "failed_post_check"
      } else {
        is_valid_model <- TRUE
      }
    } else if (status == "initial") {
      status <- "fit_null"
    }
    
    if (is_valid_model) {
      status <- "valid"
      valid_counter <- valid_counter + 1
      
      fits <- fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr", "chisq", "df"))
      res_fit <- rbind(res_fit, data.frame(Seed = curr_seed, as.list(fits)))
      
      est <- standardizedSolution(fit) %>%
        filter(op == "=~") %>%
        select(lhs, rhs, est.std, se, pvalue) %>%
        mutate(Seed = curr_seed)
      res_loadings <- rbind(res_loadings, est)
      
      # Latent factor correlations (only if 2+ latent variables)
      lv_names <- lavNames(fit, type = "lv")
      if (length(lv_names) >= 2) {
        phi <- lavInspect(fit, "cor.lv")
        phi_df <- as.data.frame(phi)
        phi_df$LV1 <- rownames(phi_df)
        
        phi_long <- phi_df %>%
          pivot_longer(
            cols      = all_of(lv_names),
            names_to  = "LV2",
            values_to = "Corr"
          ) %>%
          rename(F1 = LV1, F2 = LV2) %>%
          filter(F1 < F2) %>%
          mutate(Seed = curr_seed)
        
        res_factor_corr <- rbind(res_factor_corr, phi_long)
      }
      
      # Internal consistency reliability (Cronbach's alpha, McDonald's omega)
      rel_df <- compute_reliability_from_fit(fit)
      if (!is.null(rel_df) && nrow(rel_df) > 0) {
        rel_df$Seed <- curr_seed
        res_reliability <- rbind(res_reliability, rel_df)
      }
      
      if ("PHQBIN" %in% names(dt_cfa_num)) {
        inv_phq <- get_invariance_deltas(dt_cfa_num, cfa_model_syntax, "PHQBIN", target_items)
        if (!is.null(inv_phq))
          res_invariance <- rbind(res_invariance, cbind(Seed = curr_seed, inv_phq))
      }
      
      inv_sex <- get_invariance_deltas(dt_cfa_num, cfa_model_syntax, "SEX", target_items)
      if (!is.null(inv_sex))
        res_invariance <- rbind(res_invariance, cbind(Seed = curr_seed, inv_sex))
      
      inv_age <- get_invariance_deltas(dt_cfa_num, cfa_model_syntax, "AGEGRP", target_items)
      if (!is.null(inv_age))
        res_invariance <- rbind(res_invariance, cbind(Seed = curr_seed, inv_age))
      
      if ("PHQSCR" %in% names(dt_cfa)) {
        roc_dat <- dt_cfa
        roc_dat$Depressed <- ifelse(roc_dat$PHQSCR >= 10, 1, 0)
        
        if ("AGE" %in% names(roc_dat)) {
          roc_dat$AGEGRP <- ifelse(roc_dat$AGE < 20, "17-19", "20-24")
        }
        
        # Secondary restricted ROC: exclude cases with all target items = 0 AND PHQSCR = 0
        item_numeric <- roc_dat[, target_items, drop = FALSE]
        item_numeric[] <- lapply(item_numeric, function(x) {
          if (is.factor(x)) as.numeric(x) - 1 else x
        })
        sum_items <- rowSums(item_numeric, na.rm = TRUE)
        zero_items <- sum_items == 0
        zero_phq   <- roc_dat$PHQSCR == 0
        keep_restricted <- !(zero_items & zero_phq)
        roc_dat_restricted <- roc_dat[keep_restricted, , drop = FALSE]
        
        # Primary ROC on full sample
        res_roc <- rbind(res_roc,
                         get_roc_metrics(roc_dat, curr_seed, "Overall", target_items))
        res_roc_curves <- rbind(res_roc_curves,
                                get_roc_curve_data(roc_dat, curr_seed, "Overall", target_items))
        
        # Secondary ROC: restricted sample (exclude 0-0)
        res_roc <- rbind(res_roc,
                         get_roc_metrics(roc_dat_restricted, curr_seed,
                                         "Overall (Excl. 0-0)", target_items))
        res_roc_curves <- rbind(res_roc_curves,
                                get_roc_curve_data(roc_dat_restricted, curr_seed,
                                                   "Overall (Excl. 0-0)", target_items))
        
        if ("SEX" %in% names(roc_dat)) {
          for (g in unique(na.omit(roc_dat$SEX))) {
            sub_dat <- roc_dat[roc_dat$SEX == g, ]
            res_roc <- rbind(res_roc,
                             get_roc_metrics(sub_dat, curr_seed, paste0("Sex: ", g), target_items))
            res_roc_curves <- rbind(res_roc_curves,
                                    get_roc_curve_data(sub_dat, curr_seed, paste0("Sex: ", g), target_items))
            
            # Restricted (exclude 0-0)
            sub_dat_restricted <- roc_dat_restricted[roc_dat_restricted$SEX == g, , drop = FALSE]
            res_roc <- rbind(res_roc,
                             get_roc_metrics(sub_dat_restricted, curr_seed,
                                             paste0("Sex: ", g, " (Excl. 0-0)"), target_items))
            res_roc_curves <- rbind(res_roc_curves,
                                    get_roc_curve_data(sub_dat_restricted, curr_seed,
                                                       paste0("Sex: ", g, " (Excl. 0-0)"), target_items))
          }
        }
        
        if ("AGEGRP" %in% names(roc_dat)) {
          for (g in unique(na.omit(roc_dat$AGEGRP))) {
            sub_dat <- roc_dat[roc_dat$AGEGRP == g, ]
            res_roc <- rbind(res_roc,
                             get_roc_metrics(sub_dat, curr_seed, paste0("Age: ", g), target_items))
            res_roc_curves <- rbind(res_roc_curves,
                                    get_roc_curve_data(sub_dat, curr_seed, paste0("Age: ", g), target_items))
            
            # Restricted (exclude 0-0)
            sub_dat_restricted <- roc_dat_restricted[roc_dat_restricted$AGEGRP == g, , drop = FALSE]
            res_roc <- rbind(res_roc,
                             get_roc_metrics(sub_dat_restricted, curr_seed,
                                             paste0("Age: ", g, " (Excl. 0-0)"), target_items))
            res_roc_curves <- rbind(res_roc_curves,
                                    get_roc_curve_data(sub_dat_restricted, curr_seed,
                                                       paste0("Age: ", g, " (Excl. 0-0)"), target_items))
          }
        }
      }
    }
    
    diagnostics <- rbind(
      diagnostics,
      data.frame(Seed = curr_seed, Status = status, stringsAsFactors = FALSE)
    )
  }
  
  list(
    res_fit         = res_fit,
    res_loadings    = res_loadings,
    res_invariance  = res_invariance,
    res_roc         = res_roc,
    res_roc_curves  = res_roc_curves,
    res_factor_corr = res_factor_corr,
    res_reliability = res_reliability,
    diagnostics     = diagnostics,
    n_valid         = valid_counter,
    n_tried         = length(run_seeds)
  )
}
