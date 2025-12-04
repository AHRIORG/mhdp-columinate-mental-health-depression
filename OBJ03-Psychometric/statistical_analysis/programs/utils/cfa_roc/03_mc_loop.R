############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 03_mc_loop.R - Monte Carlo CFA loop
############################################################

run_monte_carlo_cfa <- function(dt_main, cfa_model_syntax, target_items,
                                target_set_name, run_seeds, n_mc_samples) {
  res_fit        <- data.frame()
  res_loadings   <- data.frame()
  res_invariance <- data.frame()
  res_roc        <- data.frame()
  res_roc_curves <- data.frame()
  
  valid_counter <- 0
  message(sprintf(
    "\nAttempting to collect %d valid CFA runs (from up to %d seeds)...",
    n_mc_samples, length(run_seeds)
  ))
  pb <- txtProgressBar(min = 0, max = length(run_seeds), style = 3)
  
  for (i in seq_along(run_seeds)) {
    if (valid_counter >= n_mc_samples) break
    
    curr_seed <- run_seeds[i]
    
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
    }, error = function(e) NULL)
    
    is_valid_model <- FALSE
    if (!is.null(fit)) {
      if (inspect(fit, "converged") && lavInspect(fit, "post.check")) {
        is_valid_model <- TRUE
      }
    }
    
    if (is_valid_model) {
      valid_counter <- valid_counter + 1
      
      fits <- fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr", "chisq", "df"))
      res_fit <- rbind(res_fit, data.frame(Seed = curr_seed, as.list(fits)))
      
      est <- standardizedSolution(fit) %>%
        filter(op == "=~") %>%
        select(lhs, rhs, est.std, se, pvalue) %>%
        mutate(Seed = curr_seed)
      res_loadings <- rbind(res_loadings, est)
      
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
        
        res_roc <- rbind(res_roc,
                         get_roc_metrics(roc_dat, curr_seed, "Overall", target_items))
        res_roc_curves <- rbind(res_roc_curves,
                                get_roc_curve_data(roc_dat, curr_seed, "Overall", target_items))
        
        if ("SEX" %in% names(roc_dat)) {
          for (g in unique(na.omit(roc_dat$SEX))) {
            sub_dat <- roc_dat[roc_dat$SEX == g, ]
            res_roc <- rbind(res_roc,
                             get_roc_metrics(sub_dat, curr_seed, paste0("Sex: ", g), target_items))
            res_roc_curves <- rbind(res_roc_curves,
                                    get_roc_curve_data(sub_dat, curr_seed, paste0("Sex: ", g), target_items))
          }
        }
        
        if ("AGEGRP" %in% names(roc_dat)) {
          for (g in unique(na.omit(roc_dat$AGEGRP))) {
            sub_dat <- roc_dat[roc_dat$AGEGRP == g, ]
            res_roc <- rbind(res_roc,
                             get_roc_metrics(sub_dat, curr_seed, paste0("Age: ", g), target_items))
            res_roc_curves <- rbind(res_roc_curves,
                                    get_roc_curve_data(sub_dat, curr_seed, paste0("Age: ", g), target_items))
          }
        }
      }
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  message(sprintf("Collected %d valid CFA runs (out of %d seeds tried).",
                  valid_counter, length(run_seeds)))
  
  list(
    res_fit        = res_fit,
    res_loadings   = res_loadings,
    res_invariance = res_invariance,
    res_roc        = res_roc,
    res_roc_curves = res_roc_curves
  )
}
