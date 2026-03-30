############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 06_master.R - Orchestrator
############################################################

# 1. Diagnostics summary table -----------------------------------------

summarise_diagnostics_gt <- function(diagnostics) {
  if (is.null(diagnostics) || nrow(diagnostics) == 0) return(NULL)
  
  tab <- diagnostics %>%
    count(Status, name = "N") %>%
    mutate(Proportion = N / sum(N))
  
  gt(tab) %>%
    fmt_number(columns = "Proportion", decimals = 3) %>%
    cols_label(
      Status     = "Status",
      N          = "Number of Seeds",
      Proportion = "Proportion"
    ) %>%
    tab_header(
      title = "Monte Carlo Seed Diagnostics"
    ) %>%
    tab_source_note(
      source_note = paste(
        "Status interpretation: 'valid' = retained CFA run;",
        "'no_convergence' = model did not converge;",
        "'failed_post_check' = failed lavaan post.check;",
        "'fit_error'/'fit_null' = lavaan failed to provide a usable solution."
      )
    )
}

# 2. ROC refresh from stored diagnostics (no CFA re-fit) ---------------

refresh_roc_from_valid_seeds <- function(dt_main,
                                         target_items,
                                         diagnostics,
                                         target_set_name) {
  # If PHQ-09 is the target set, ROC analysis is redundant
  if (identical(target_set_name, "PHQ-09 (Original)")) {
    return(list(
      res_roc        = data.frame(),
      res_roc_curves = data.frame()
    ))
  }
  
  if (is.null(diagnostics) || nrow(diagnostics) == 0) return(NULL)
  if (!all(c("Seed", "Status") %in% names(diagnostics))) return(NULL)
  
  seeds_valid <- unique(diagnostics$Seed[diagnostics$Status == "valid"])
  if (length(seeds_valid) == 0) return(NULL)
  
  res_roc        <- data.frame()
  res_roc_curves <- data.frame()
  
  for (curr_seed in seeds_valid) {
    set.seed(curr_seed)
    n_rows <- nrow(dt_main)
    idx_explore <- sample(seq_len(n_rows), floor(0.5 * n_rows), replace = FALSE)
    idx_confirm <- setdiff(seq_len(n_rows), idx_explore)
    dt_cfa <- dt_main[idx_confirm, , drop = FALSE]
    
    # Need PHQSCR for ROC benchmark
    if (!("PHQSCR" %in% names(dt_cfa))) next
    
    roc_dat <- dt_cfa
    roc_dat$Depressed <- ifelse(roc_dat$PHQSCR >= 10, 1, 0)
    
    # Age grouping
    if ("AGE" %in% names(roc_dat)) {
      roc_dat$AGEGRP <- ifelse(roc_dat$AGE < 20, "17-19", "20-24")
    }
    
    # ---- Full-sample ROC ----
    tmp  <- get_roc_metrics(roc_dat, curr_seed, "Overall", target_items)
    if (!is.null(tmp))  res_roc        <- rbind(res_roc, tmp)
    tmpc <- get_roc_curve_data(roc_dat, curr_seed, "Overall", target_items)
    if (!is.null(tmpc)) res_roc_curves <- rbind(res_roc_curves, tmpc)
    
    # ---- Build restricted sample: exclude SSQ=0 & PHQ=0 ----
    item_data <- roc_dat[, target_items, drop = FALSE]
    item_data[] <- lapply(item_data, function(x) {
      if (is.factor(x)) as.numeric(x) - 1 else x
    })
    ssq_sum <- rowSums(item_data, na.rm = TRUE)
    zero_items <- ssq_sum == 0
    zero_phq   <- roc_dat$PHQSCR == 0
    
    keep_restricted <- !(zero_items & zero_phq)
    roc_restricted  <- roc_dat[keep_restricted, , drop = FALSE]
    
    # Overall restricted
    tmp  <- get_roc_metrics(roc_restricted, curr_seed,
                            "Overall (Excl. 0-0)", target_items)
    if (!is.null(tmp))  res_roc        <- rbind(res_roc, tmp)
    tmpc <- get_roc_curve_data(roc_restricted, curr_seed,
                               "Overall (Excl. 0-0)", target_items)
    if (!is.null(tmpc)) res_roc_curves <- rbind(res_roc_curves, tmpc)
    
    # ---- Sex subgroups: full + restricted ----
    if ("SEX" %in% names(roc_dat)) {
      for (g in unique(na.omit(roc_dat$SEX))) {
        # Full-sample subgroup
        sub_dat <- roc_dat[roc_dat$SEX == g, , drop = FALSE]
        tmp  <- get_roc_metrics(sub_dat, curr_seed,
                                paste0("Sex: ", g), target_items)
        if (!is.null(tmp))  res_roc        <- rbind(res_roc, tmp)
        tmpc <- get_roc_curve_data(sub_dat, curr_seed,
                                   paste0("Sex: ", g), target_items)
        if (!is.null(tmpc)) res_roc_curves <- rbind(res_roc_curves, tmpc)
        
        # Restricted subgroup
        sub_dat_re <- roc_restricted[roc_restricted$SEX == g, , drop = FALSE]
        tmp  <- get_roc_metrics(sub_dat_re, curr_seed,
                                paste0("Sex: ", g, " (Excl. 0-0)"), target_items)
        if (!is.null(tmp))  res_roc        <- rbind(res_roc, tmp)
        tmpc <- get_roc_curve_data(sub_dat_re, curr_seed,
                                   paste0("Sex: ", g, " (Excl. 0-0)"), target_items)
        if (!is.null(tmpc)) res_roc_curves <- rbind(res_roc_curves, tmpc)
      }
    }
    
    # ---- Age subgroups: full + restricted ----
    if ("AGEGRP" %in% names(roc_dat)) {
      for (g in unique(na.omit(roc_dat$AGEGRP))) {
        # Full-sample subgroup
        sub_dat <- roc_dat[roc_dat$AGEGRP == g, , drop = FALSE]
        tmp  <- get_roc_metrics(sub_dat, curr_seed,
                                paste0("Age: ", g), target_items)
        if (!is.null(tmp))  res_roc        <- rbind(res_roc, tmp)
        tmpc <- get_roc_curve_data(sub_dat, curr_seed,
                                   paste0("Age: ", g), target_items)
        if (!is.null(tmpc)) res_roc_curves <- rbind(res_roc_curves, tmpc)
        
        # Restricted subgroup
        sub_dat_re <- roc_restricted[roc_restricted$AGEGRP == g, , drop = FALSE]
        tmp  <- get_roc_metrics(sub_dat_re, curr_seed,
                                paste0("Age: ", g, " (Excl. 0-0)"), target_items)
        if (!is.null(tmp))  res_roc        <- rbind(res_roc, tmp)
        tmpc <- get_roc_curve_data(sub_dat_re, curr_seed,
                                   paste0("Age: ", g, " (Excl. 0-0)"), target_items)
        if (!is.null(tmpc)) res_roc_curves <- rbind(res_roc_curves, tmpc)
      }
    }
  }
  
  list(
    res_roc        = res_roc,
    res_roc_curves = res_roc_curves
  )
}

# 3. Master pipeline ---------------------------------------------------

.default_psychometric_data_path <- function() {
  override <- Sys.getenv("COLU_PSYCHOMETRIC_RDATA", unset = "")
  if (nzchar(override)) return(override)
  here::here("data/inputs/dt_psychometric.RData")
}

run_ssq_phq_mc_pipeline <- function(
    data_path            = .default_psychometric_data_path(),
    data_object_name     = NULL,
    target_set_name      = "SSQ-10 (Theoretical Depressive Symptoms)",
    dataset_label        = "Isisekelo Sempilo",
    manual_method        = "Parallel",
    manual_n_factors     = 3,
    n_mc_samples         = 100,
    min_items_per_factor = 3,
    manual_merge_pairs   = list(c(3, 1))
) {
  # Load all required packages (01_packages.R)
  load_required_packages()
  
  # Flexible data loading: either from an existing object or from file (02_config.R)
  dt_main <- load_main_data(path = data_path, object_name = data_object_name)
  
  # Item sets and analysis configuration (02_config.R)
  item_sets_list <- create_item_sets()
  cfg <- analysis_config(
    item_sets_list       = item_sets_list,
    target_set_name      = target_set_name,
    dataset_label        = dataset_label,
    manual_method        = manual_method,
    manual_n_factors     = manual_n_factors,
    n_mc_samples         = n_mc_samples,
    min_items_per_factor = min_items_per_factor,
    manual_merge_pairs   = manual_merge_pairs
  )
  
  # Flag: PHQ-09 runs use PHQ-09 as the gold standard; ROC analysis is redundant
  is_phq <- identical(cfg$target_set_name, "PHQ-09 (Original)")
  
  # Load bootstrap EFA results and safe name (02_config.R)
  boot_info     <- load_bootstrap_results(cfg$dataset_label, cfg$target_set_name)
  boot_results  <- boot_info$boot_results
  safe_set_name <- boot_info$safe_set_name
  
  # Select factor solution (Kaiser vs Parallel, #factors) (02_config.R)
  fac_sel <- select_factor_solution(
    boot_results       = boot_results,
    manual_method      = cfg$manual_method,
    manual_n_factors   = cfg$manual_n_factors
  )
  
  # Short run identifier for file naming (method, k, n)
  method_code <- switch(
    fac_sel$target_method,
    "Parallel" = "pa",
    "Kaiser"   = "ka",
    tolower(fac_sel$target_method)
  )
  
  run_tag <- paste0(
    method_code,
    "_k", fac_sel$target_n,
    "_n", cfg$n_mc_samples
  )
  
  # Extract winning structure (02_config.R)
  struct_info <- extract_winning_structure(
    boot_results  = boot_results,
    target_method = fac_sel$target_method,
    target_n      = fac_sel$target_n
  )
  
  factor_list  <- struct_info$factor_list
  valid_runs   <- struct_info$valid_runs
  winning_sig  <- struct_info$winning_sig
  
  # Apply any manual factor merges (if specified) (02_config.R)
  if (!is.null(cfg$manual_merge_pairs)) {
    factor_list <- merge_factors_manual(factor_list, cfg$manual_merge_pairs)
  }
  
  # Build CFA syntax from final factor list (02_config.R)
  cfa_model_syntax <- generate_cfa_syntax(factor_list)
  
  # Monte Carlo seeds (03_mc_loop.R) – use valid runs matching winning signature
  run_seeds <- select_mc_seeds(
    valid_runs    = valid_runs,
    winning_sig   = winning_sig,
    n_mc_samples  = cfg$n_mc_samples
  )
  
  # Monte Carlo results path (with method/factor/n tag)
  mc_file <- paste0(
    "mc_cfa_validation_",
    safe_set_name, "_",
    run_tag,
    ".rds"
  )
  mc_path <- here("statistical_analysis/output/objects/mc_cfa", mc_file)
  
  # 3a. Run or reuse Monte Carlo CFA ----------------------------------
  
  if (!file.exists(mc_path)) {
    message("No existing Monte Carlo results found. Running Monte Carlo validation...")
    
    mc_res <- suppressWarnings(
      run_monte_carlo_cfa(
        dt_main          = dt_main,
        cfa_model_syntax = cfa_model_syntax,
        target_items     = cfg$target_items,
        target_set_name  = cfg$target_set_name,
        run_seeds        = run_seeds,
        n_mc_samples     = cfg$n_mc_samples
      )
    )
    
    # For PHQ-09 runs, ensure ROC objects are empty data.frames before saving
    if (is_phq) {
      mc_res$res_roc        <- data.frame()
      mc_res$res_roc_curves <- data.frame()
    }
    
    saveRDS(
      list(
        fit         = mc_res$res_fit,
        loadings    = mc_res$res_loadings,
        invariance  = mc_res$res_invariance,
        roc         = mc_res$res_roc,
        roc_curves  = mc_res$res_roc_curves,
        factor_corr = mc_res$res_factor_corr,
        reliability = mc_res$res_reliability,
        diagnostics = mc_res$diagnostics,
        n_valid     = mc_res$n_valid,
        n_tried     = mc_res$n_tried
      ),
      file = mc_path
    )
  } else {
    message("Re-using existing Monte Carlo results from disk.")
    
    stored <- readRDS(mc_path)
    mc_res <- list(
      res_fit         = stored$fit,
      res_loadings    = stored$loadings,
      res_invariance  = stored$invariance,
      res_roc         = if ("roc"        %in% names(stored)) stored$roc        else NULL,
      res_roc_curves  = if ("roc_curves" %in% names(stored)) stored$roc_curves else NULL,
      res_factor_corr = if ("factor_corr"%in% names(stored)) stored$factor_corr else NULL,
      res_reliability = if ("reliability"%in% names(stored)) stored$reliability else NULL,
      diagnostics     = if ("diagnostics"%in% names(stored)) stored$diagnostics else NULL,
      n_valid         = if ("n_valid"    %in% names(stored)) stored$n_valid    else NA_integer_,
      n_tried         = if ("n_tried"    %in% names(stored)) stored$n_tried    else NA_integer_
    )
    
    if (is_phq) {
      # For PHQ-09 runs, ROC is redundant; drop any stored ROC objects
      mc_res$res_roc        <- data.frame()
      mc_res$res_roc_curves <- data.frame()
      
      stored$roc        <- NULL
      stored$roc_curves <- NULL
      saveRDS(stored, file = mc_path)
    } else {
      # If older results are missing prevalence information, refresh ROC
      needs_roc_refresh <- is.null(mc_res$res_roc) ||
        !all(c("Prevalence_PHQ", "Prevalence_SSQ") %in% names(mc_res$res_roc))
      
      if (needs_roc_refresh && !is.null(mc_res$diagnostics) && nrow(mc_res$diagnostics) > 0) {
        message("Refreshing ROC metrics and prevalence using valid seeds only (no CFA re-fitting)...")
        roc_ref <- refresh_roc_from_valid_seeds(
          dt_main         = dt_main,
          target_items    = cfg$target_items,
          diagnostics     = mc_res$diagnostics,
          target_set_name = cfg$target_set_name
        )
        if (!is.null(roc_ref)) {
          mc_res$res_roc        <- roc_ref$res_roc
          mc_res$res_roc_curves <- roc_ref$res_roc_curves
          
          stored$roc        <- mc_res$res_roc
          stored$roc_curves <- mc_res$res_roc_curves
          
          saveRDS(stored, file = mc_path)
        }
      }
    }
  }
  
  # For PHQ-09 runs, enforce empty ROC components before summarising/plotting
  if (is_phq) {
    mc_res$res_roc        <- data.frame()
    mc_res$res_roc_curves <- data.frame()
  }
  
  # 3b. Summaries & tables (04_summary.R) ------------------------------
  
  gt_tables <- summarise_results(
    res_fit         = mc_res$res_fit,
    res_loadings    = mc_res$res_loadings,
    res_invariance  = mc_res$res_invariance,
    res_roc         = mc_res$res_roc,
    res_factor_corr = mc_res$res_factor_corr,
    res_reliability = mc_res$res_reliability
  )
  
  diag_table <- summarise_diagnostics_gt(mc_res$diagnostics)
  gt_tables$diagnostics <- diag_table
  
  summary_tables_file <- paste0("mc_cfa_summary_tables_", safe_set_name, "_", run_tag, ".rds")
  saveRDS(
    gt_tables,
    file = here("statistical_analysis/output/objects/mc_cfa", summary_tables_file)
  )
  message(paste("Summary tables saved to:", summary_tables_file))
  
  # 3c. Plots (05_plots.R) --------------------------------------------
  
  plot_list <- generate_plots(
    res_fit         = mc_res$res_fit,
    res_loadings    = mc_res$res_loadings,
    res_invariance  = mc_res$res_invariance,
    res_roc         = mc_res$res_roc,
    res_roc_curves  = mc_res$res_roc_curves,
    safe_set_name   = safe_set_name,
    target_set_name = cfg$target_set_name,
    n_mc_samples    = cfg$n_mc_samples
  )
  
  plots_file <- paste0("mc_cfa_plot_objects_", safe_set_name, "_", run_tag, ".rds")
  saveRDS(
    plot_list,
    file = here("statistical_analysis/output/objects/mc_cfa", plots_file)
  )
  message(paste("Plot objects saved to:", plots_file))
  
  invisible(list(
    mc_results       = mc_res,
    gt_tables        = gt_tables,
    plot_objects     = plot_list,
    cfa_model_syntax = cfa_model_syntax,
    seeds            = run_seeds,
    factor_list      = factor_list
  ))
}

# To run the full pipeline after sourcing all helper scripts:
# run_ssq_phq_mc_pipeline()
