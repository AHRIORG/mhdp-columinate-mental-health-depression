############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 06_master.R - Orchestrator
############################################################

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
  
  mc_res <- run_monte_carlo_cfa(
    dt_main          = dt_main,
    cfa_model_syntax = cfa_model_syntax,
    target_items     = cfg$target_items,
    target_set_name  = cfg$target_set_name,
    run_seeds        = run_seeds,
    n_mc_samples     = cfg$n_mc_samples
  )
  
  summarise_results(
    res_fit        = mc_res$res_fit,
    res_loadings   = mc_res$res_loadings,
    res_invariance = mc_res$res_invariance,
    res_roc        = mc_res$res_roc
  )
  
  save_mc_results(
    res_fit        = mc_res$res_fit,
    res_loadings   = mc_res$res_loadings,
    res_invariance = mc_res$res_invariance,
    res_roc        = mc_res$res_roc,
    res_roc_curves = mc_res$res_roc_curves,
    safe_set_name  = safe_set_name
  )
  
  generate_plots(
    res_fit         = mc_res$res_fit,
    res_loadings    = mc_res$res_loadings,
    res_invariance  = mc_res$res_invariance,
    res_roc         = mc_res$res_roc,
    res_roc_curves  = mc_res$res_roc_curves,
    safe_set_name   = safe_set_name,
    target_set_name = cfg$target_set_name,
    n_mc_samples    = cfg$n_mc_samples
  )
}