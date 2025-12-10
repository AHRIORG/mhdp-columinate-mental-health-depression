############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 03_mc_loop.R - Monte Carlo CFA loop
############################################################

run_monte_carlo_cfa <- function(dt_main, cfa_model_syntax, target_items,
                                target_set_name, run_seeds, n_mc_samples) {
  res_fit         <- data.frame()
  res_loadings    <- data.frame()
  res_invariance  <- data.frame()
  res_roc         <- data.frame()
  res_roc_curves  <- data.frame()
  res_factor_corr <- data.frame()
  
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
    diagnostics     = diagnostics,
    n_valid         = valid_counter,
    n_tried         = length(run_seeds)
  )
}
