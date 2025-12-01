# Psychometric Bootstrap EFA – Pipeline Wrapper

# Dependencies:
# library(psych)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(here)
# library(scales)
# library(parallel)
# library(doParallel)
# library(foreach)
#
# Source helper dependance functions:
source("statistical_analysis/programs/utils/Psych Split & EFA Helpers.R")
source("statistical_analysis/programs/utils/Psych Evaluate All Splits.R")
source("statistical_analysis/programs/utils/Psych Retention Loadings Eigen.R")
source("statistical_analysis/programs/utils/Psych Cooccurrence Consensus.R")

# ------------------------------------------------------------------------------
# 1. Pipeline Wrapper
# ------------------------------------------------------------------------------

run_psychometric_pipeline <- function(df_data,
                                      dataset_label,
                                      item_sets_list,
                                      n_boot,
                                      n_pa,
                                      loading_cutoff = 0.4,
                                      use_parallel   = TRUE,
                                      n_cores        = parallel::detectCores() - 4
                                      ) {
  
  message("\n========================================================")
  message(paste0("  STARTING ANALYSIS PIPELINE: ", dataset_label))
  message("========================================================")
  
  splits_list <- generate_ec_splits(DATA = df_data, NBOOT = n_boot)
  message(paste(
    "Generated", length(splits_list),
    "sets of split indices for", dataset_label, "(includes buffer)"
  ))
  
  results_collection <- list()
  
  all_required_items <- unique(unlist(item_sets_list))
  missing_items <- setdiff(all_required_items, names(df_data))
  if (length(missing_items) > 0) {
    message("Error: Some required items not found in dataset columns.")
    message("Missing: ", paste(missing_items, collapse = ", "))
    return(NULL)
  }
  
  output_dir <- here::here("statistical_analysis/output/objects")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (set_name in names(item_sets_list)) {
    current_items <- item_sets_list[[set_name]]
    message(paste0("\n  Processing Item Set: ", set_name, "..."))
    
    safe_set_name <- gsub("[\n ]+", "_", set_name)
    safe_ds_label <- gsub("[\n ]+", "_", dataset_label)
    file_name     <- paste0(
      "boot_", n_boot, "_", safe_ds_label, "_", safe_set_name, "_results.rds"
    )
    file_path <- file.path(output_dir, file_name)
    
    run_analysis <- TRUE
    if (file.exists(file_path)) {
      message("    Found existing saved results.")
      if (interactive()) {
        user_resp <- readline(prompt = "    Do you want to RERUN the analysis? (y/n): ")
        if (tolower(substr(user_resp, 1, 1)) != "y") run_analysis <- FALSE
      } else {
        message("    Loading existing file (non-interactive mode).")
        run_analysis <- FALSE
      }
    }
    
    if (run_analysis) {
      message(paste0("    Running EFA bootstrapping (Target N=", n_boot, ")..."))
      
      assessment_results <- evaluate_all_splits(
        splits_list    = splits_list,
        original_data  = df_data,
        items          = current_items,
        n_iter         = n_pa,
        target_n       = n_boot,
        loading_cutoff = loading_cutoff,
        use_parallel   = use_parallel,
        n_cores        = n_cores
      )
      
      if (!is.null(assessment_results)) {
        saveRDS(assessment_results, file = file_path)
        message(paste0("    Saved results to: ", file_path))
      } else {
        message("    Warning: No results to save.")
      }
    } else {
      assessment_results <- readRDS(file_path)
      message(paste0("    Loaded results from: ", file_path))
    }
    
    if (!is.null(assessment_results)) {
      current_results <- list()
      current_results$raw_data <- assessment_results
      
      plot_sub <- set_name
      current_results$retention_analysis <-
        analyze_factor_retention(assessment_results, custom_subtitle = plot_sub)
      
      analyze_winning_structure <- function(method_name, prefix, count_col) {
        valid_counts <- stats::na.omit(assessment_results[[count_col]])
        if (length(valid_counts) == 0) return(NULL)
        
        counts_table  <- table(valid_counts)
        winning_count <- as.numeric(names(which.max(counts_table)))
        winning_freq  <- max(counts_table)
        percent_win   <- (winning_freq / sum(counts_table)) * 100
        
        sub_df <- assessment_results %>%
          dplyr::filter(.data[[count_col]] == winning_count, Is_Suitable == TRUE)
        if (nrow(sub_df) == 0) return(NULL)
        
        consensus <- identify_consensus_model(
          subset_df     = sub_df,
          method_prefix = prefix,
          winning_count = winning_count
        )
        
        sort_res <- get_sorted_items_by_stability(
          subset_df     = sub_df,
          method_prefix = prefix,
          winning_count = winning_count,
          all_items     = current_items
        )
        final_item_order <- sort_res$sorted_items
        probs_mat        <- sort_res$probs_mat
        
        mat <- calculate_cooccurrence_matrix(sub_df, prefix, current_items)
        
        p_heat <- plot_cooccurrence_heatmap(
          mat,
          subtitle   = paste0(set_name, "\n", paste0(method_name, " (A ", winning_count, "-Factor Solution)")),
          caption    = NULL,
          item_order = final_item_order
        )
        
        p_stab <- plot_item_factor_stability(
          probs_mat      = probs_mat,
          item_order     = final_item_order,
          winning_count  = winning_count,
          consensus_info = consensus,
          subtitle       = paste0(set_name, "\n", paste0(method_name, " (A ", winning_count, "-Factor Solution)"))
        )
        
        list(
          winning_count   = winning_count,
          percent_win     = percent_win,
          matrix          = mat,
          plot_heatmap    = p_heat,
          plot_stability  = p_stab,
          consensus_model = consensus
        )
      }
      
      current_results$cooccurrence <- list(
        Kaiser   = analyze_winning_structure("Kaiser",   "k_factor_", "k_factor_f"),
        Parallel = analyze_winning_structure("Parallel", "p_factor_", "p_factor_f")
      )
      
      results_collection[[set_name]] <- current_results
    }
  }
  
  results_collection
}

