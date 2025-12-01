# Psychometric Bootstrap EFA – Bootstrap Evaluation (Sequential & Parallel)

# ------------------------------------------------------------------------------
# Dependencies assumed:
# library(parallel)
# library(doParallel)
# library(foreach)
# library(dplyr)
# library(tidyr)
#
# check_efa_suitability() must be sourced before this file.
# ------------------------------------------------------------------------------

#' Assess EFA Suitability for All Splits (Sequential or Parallel)
#'
#' Iterates through the provided splits list. If a split is unsuitable,
#' it is discarded. Stops when 'target_n' suitable splits are found.
#' Optionally parallelised using foreach + doParallel.
#'
#' @param splits_list List output from generate_ec_splits (includes buffer).
#' @param original_data The original dataframe.
#' @param items The vector of item names to analyze.
#' @param n_iter Number of parallel analysis iterations.
#' @param target_n The required number of valid/suitable splits (usually NBOOT).
#' @param loading_cutoff Cutoff value for factor loadings.
#' @param use_parallel Logical. If TRUE, evaluate splits in parallel.
#' @param n_cores Integer. Number of cores to use when use_parallel = TRUE.
#' 
#' @return A dataframe summary of the 'target_n' suitable iterations.
evaluate_all_splits <- function(splits_list,
                                original_data,
                                items,
                                n_iter = 20,
                                target_n = NULL,
                                loading_cutoff = 0.4,
                                use_parallel = FALSE,
                                n_cores = parallel::detectCores() - 1) {
  
  if (is.null(target_n)) target_n <- length(splits_list)
  
  standardize_and_bind <- function(results_list) {
    if (length(results_list) == 0) return(NULL)
    
    all_cols <- unique(unlist(lapply(results_list, names)))
    
    standardize_df <- function(df, all_cols) {
      missing <- setdiff(all_cols, names(df))
      if (length(missing) > 0) df[missing] <- NA
      df[, all_cols]
    }
    
    final_df <- do.call(rbind, lapply(results_list, standardize_df, all_cols = all_cols))
    
    base_cols   <- c("Iteration", "Seed", "KMO_Overall", "Is_Suitable", "k_factor_f", "p_factor_f")
    factor_cols <- setdiff(names(final_df), base_cols)
    factor_cols <- factor_cols[order(factor_cols)]
    
    final_df[, c(base_cols, factor_cols)]
  }
  
  # ---------------- SEQUENTIAL VERSION -----------------
  if (!use_parallel) {
    
    results_list    <- list()
    suitable_count  <- 0L
    
    for (i in seq_along(splits_list)) {
      if (suitable_count >= target_n) break
      
      split_info <- splits_list[[i]]
      dt_explore <- original_data[split_info$exploratory, , drop = FALSE]
      
      row_result <- tryCatch({
        res  <- check_efa_suitability(
          dt_explore,
          items          = items,
          n_iter         = n_iter,
          loading_cutoff = loading_cutoff
        )
        summ <- res$Suitability_Summary
        
        row_df <- data.frame(
          Iteration   = split_info$iteration,
          Seed        = split_info$seed,
          KMO_Overall = summ$KMO_Overall,
          Is_Suitable = summ$Is_Suitable,
          k_factor_f  = summ$k_factor_f,
          p_factor_f  = summ$p_factor_f,
          stringsAsFactors = FALSE
        )
        
        if (length(summ$factor_items_k) > 0) {
          for (f_idx in names(summ$factor_items_k)) {
            col_name <- paste0("k_factor_", f_idx)
            row_df[[col_name]] <- summ$factor_items_k[[f_idx]]
          }
        }
        
        if (length(summ$factor_items_p) > 0) {
          for (f_idx in names(summ$factor_items_p)) {
            col_name <- paste0("p_factor_", f_idx)
            row_df[[col_name]] <- summ$factor_items_p[[f_idx]]
          }
        }
        
        row_df$k_loadings    <- list(summ$factor_loadings_k)
        row_df$p_loadings    <- list(summ$factor_loadings_p)
        row_df$eig_actual    <- list(summ$eigenvalues_actual)
        row_df$eig_simulated <- list(summ$eigenvalues_simulated)
        
        row_df
      }, error = function(e) {
        row_df <- data.frame(
          Iteration   = split_info$iteration,
          Seed        = split_info$seed,
          KMO_Overall = NA_real_,
          Is_Suitable = FALSE,
          k_factor_f  = NA_integer_,
          p_factor_f  = NA_integer_,
          stringsAsFactors = FALSE
        )
        row_df$k_loadings    <- list(NULL)
        row_df$p_loadings    <- list(NULL)
        row_df$eig_actual    <- list(NULL)
        row_df$eig_simulated <- list(NULL)
        row_df
      })
      
      if (isTRUE(row_result$Is_Suitable)) {
        suitable_count <- suitable_count + 1L
        results_list[[length(results_list) + 1L]] <- row_result
      }
    }
    
    if (suitable_count < target_n) {
      warning(sprintf(
        "Warning: Could not find %d suitable splits. Only found %d out of %d generated splits. Increase buffer size?",
        target_n, suitable_count, length(splits_list)
      ))
    }
    
    return(standardize_and_bind(results_list))
  }
  
  # ---------------- PARALLEL VERSION -----------------
  
  n_cores <- max(1L, n_cores)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }, add = TRUE)
  
  results_list <- foreach::foreach(
    i = seq_along(splits_list),
    .packages = c("psych"),
    .export   = c("check_efa_suitability")
  ) %dopar% {
    split_info <- splits_list[[i]]
    dt_explore <- original_data[split_info$exploratory, , drop = FALSE]
    
    tryCatch({
      res  <- check_efa_suitability(
        dt_explore,
        items          = items,
        n_iter         = n_iter,
        loading_cutoff = loading_cutoff
      )
      summ <- res$Suitability_Summary
      
      row_df <- data.frame(
        Iteration   = split_info$iteration,
        Seed        = split_info$seed,
        KMO_Overall = summ$KMO_Overall,
        Is_Suitable = summ$Is_Suitable,
        k_factor_f  = summ$k_factor_f,
        p_factor_f  = summ$p_factor_f,
        stringsAsFactors = FALSE
      )
      
      if (length(summ$factor_items_k) > 0) {
        for (f_idx in names(summ$factor_items_k)) {
          col_name <- paste0("k_factor_", f_idx)
          row_df[[col_name]] <- summ$factor_items_k[[f_idx]]
        }
      }
      
      if (length(summ$factor_items_p) > 0) {
        for (f_idx in names(summ$factor_items_p)) {
          col_name <- paste0("p_factor_", f_idx)
          row_df[[col_name]] <- summ$factor_items_p[[f_idx]]
        }
      }
      
      row_df$k_loadings    <- list(summ$factor_loadings_k)
      row_df$p_loadings    <- list(summ$factor_loadings_p)
      row_df$eig_actual    <- list(summ$eigenvalues_actual)
      row_df$eig_simulated <- list(summ$eigenvalues_simulated)
      
      row_df
    }, error = function(e) {
      row_df <- data.frame(
        Iteration   = split_info$iteration,
        Seed        = split_info$seed,
        KMO_Overall = NA_real_,
        Is_Suitable = FALSE,
        k_factor_f  = NA_integer_,
        p_factor_f  = NA_integer_,
        stringsAsFactors = FALSE
      )
      row_df$k_loadings    <- list(NULL)
      row_df$p_loadings    <- list(NULL)
      row_df$eig_actual    <- list(NULL)
      row_df$eig_simulated <- list(NULL)
      row_df
    })
  }
  
  final_df <- standardize_and_bind(results_list)
  if (is.null(final_df)) return(NULL)
  
  final_df <- final_df[final_df$Is_Suitable, , drop = FALSE]
  n_suitable <- nrow(final_df)
  
  if (n_suitable == 0L) {
    warning("Warning: No suitable splits found in parallel evaluation.")
    return(NULL)
  }
  
  if (n_suitable < target_n) {
    warning(sprintf(
      "Warning: Could not find %d suitable splits. Only found %d out of %d generated splits. Increase buffer size?",
      target_n, n_suitable, length(splits_list)
    ))
  } else if (n_suitable > target_n) {
    final_df <- final_df[1:target_n, , drop = FALSE]
  }
  
  final_df
}
