############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 06_master.R - Orchestrator
############################################################

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

run_ssq_phq_mc_pipeline <- function(
    data_path            = here("../../data_management/data/co-luminate/adam/dt_psychometric.RData"),
    data_object_name     = NULL,
    target_set_name      = "SSQ-10 (Theoretical Depressive Symptoms)",
    dataset_label        = "Isisekelo Sempilo",
    manual_method        = "Parallel",
    manual_n_factors     = 3,
    n_mc_samples         = 100,
    min_items_per_factor = 3,
    manual_merge_pairs   = list(c(3, 1))
) {
  load_required_packages()
  
  # Flexible data loading: either from an existing object or from file
  dt_main <- load_main_data(path = data_path, object_name = data_object_name)
  
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
  
  boot_info    <- load_bootstrap_results(cfg$dataset_label, cfg$target_set_name)
  boot_results <- boot_info$boot_results
  safe_set_name <- boot_info$safe_set_name
  
  fac_sel <- select_factor_solution(
    boot_results,
    cfg$manual_method,
    cfg$manual_n_factors
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
  
  struct_info <- extract_winning_structure(
    boot_results,
    fac_sel$target_method,
    fac_sel$target_n
  )
  
  factor_list <- struct_info$factor_list
  
  if (!is.null(cfg$manual_merge_pairs)) {
    factor_list <- merge_factors_manual(factor_list, cfg$manual_merge_pairs)
  } else {
    factor_list <- merge_factors_programmatic(
      factor_list          = factor_list,
      valid_runs           = struct_info$valid_runs,
      target_col_prefix    = struct_info$target_col_prefix,
      target_items         = cfg$target_items,
      min_items_per_factor = cfg$min_items_per_factor
    )
  }
  
  cfa_model_syntax <- build_cfa_syntax(factor_list)
  
  run_seeds <- select_mc_seeds(
    valid_runs   = struct_info$valid_runs,
    signatures   = struct_info$signatures,
    winning_sig  = struct_info$winning_sig,
    n_mc_samples = cfg$n_mc_samples * 2  # oversample seeds in case of invalid fits
  )
  
  # -------------------------------------------------------------------
  # Check for existing Monte Carlo results on disk to avoid re-running
  # heavy computations unnecessarily.
  # -------------------------------------------------------------------
  mc_file <- paste0("mc_cfa_validation_", safe_set_name, "_", run_tag, ".rds")
  mc_path <- here("statistical_analysis/output/objects/mc_cfa", mc_file)
  
  if (file.exists(mc_path)) {
    message(paste("Existing Monte Carlo results found at:", mc_path))
    resp <- readline(prompt = "Do you want to re-run the Monte Carlo validation? [y/N]: ")
    
    if (tolower(trimws(resp)) %in% c("y", "yes")) {
      message("Re-running Monte Carlo validation...")
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
      
      saveRDS(
        list(
          fit        = mc_res$res_fit,
          loadings   = mc_res$res_loadings,
          invariance = mc_res$res_invariance,
          roc        = mc_res$res_roc,
          roc_curves = mc_res$res_roc_curves
        ),
        file = mc_path
      )
    } else {
      message("Re-using existing Monte Carlo results from disk.")
      stored <- readRDS(mc_path)
      mc_res <- list(
        res_fit        = stored$fit,
        res_loadings   = stored$loadings,
        res_invariance = stored$invariance,
        res_roc        = stored$roc,
        res_roc_curves = stored$roc_curves,
        diagnostics    = if ("diagnostics" %in% names(stored)) stored$diagnostics else NULL,
        n_valid        = if ("n_valid" %in% names(stored)) stored$n_valid else NA_integer_,
        n_tried        = if ("n_tried" %in% names(stored)) stored$n_tried else NA_integer_
      )
    }
  } else {
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
    
    saveRDS(
      list(
        fit         = mc_res$res_fit,
        loadings    = mc_res$res_loadings,
        invariance  = mc_res$res_invariance,
        roc         = mc_res$res_roc,
        roc_curves  = mc_res$res_roc_curves,
        diagnostics = mc_res$diagnostics,
        n_valid     = mc_res$n_valid,
        n_tried     = mc_res$n_tried
      ),
      file = mc_path
    )
  }
  
  gt_tables <- summarise_results(
    res_fit        = mc_res$res_fit,
    res_loadings   = mc_res$res_loadings,
    res_invariance = mc_res$res_invariance,
    res_roc        = mc_res$res_roc
  )
  
  diag_table <- summarise_diagnostics_gt(mc_res$diagnostics)
  gt_tables$diagnostics <- diag_table
  
  summary_tables_file <- paste0("mc_cfa_summary_tables_", safe_set_name, "_", run_tag, ".rds")
  saveRDS(
    gt_tables,
    file = here("statistical_analysis/output/objects/mc_cfa", summary_tables_file)
  )
  message(paste("Summary tables saved to:", summary_tables_file))
  
  save_mc_results(
    res_fit        = mc_res$res_fit,
    res_loadings   = mc_res$res_loadings,
    res_invariance = mc_res$res_invariance,
    res_roc        = mc_res$res_roc,
    res_roc_curves = mc_res$res_roc_curves,
    safe_set_name  = safe_set_name
  )
  
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