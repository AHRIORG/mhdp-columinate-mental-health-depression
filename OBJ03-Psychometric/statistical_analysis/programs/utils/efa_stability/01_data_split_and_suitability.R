# Psychometric Bootstrap EFA – Splitting & EFA Suitability Helpers

# ------------------------------------------------------------------------------
# Generate Exploratory/Confirmatory Split Indices
# ------------------------------------------------------------------------------

#' Generate Exploratory/Confirmatory Split Indices
#'
#' @param DATA A data frame or matrix to be split.
#' @param NBOOT Integer. The number of bootstrap/split iterations to generate.
#' @param seed Integer. Master seed to control the generation of iteration-specific seeds.
#'
#' @return A list of length NBOOT + buffer. Each element contains a list with:
#'         - iteration: The iteration number
#'         - seed: The specific seed used for this split
#'         - exploratory: Row indices for the exploratory set (50%)
#'         - confirmatory: Row indices for the confirmatory set (50%)

generate_ec_splits <- function(DATA, NBOOT, seed = 123) {
  set.seed(seed)
  
  n_buffer <- NBOOT * 2
  n_total  <- NBOOT + n_buffer
  
  iteration_seeds <- sample(1:9999999, size = n_total, replace = FALSE)
  
  n_rows <- nrow(DATA)
  split_size <- floor(0.5 * n_rows)
  
  split_results <- vector("list", n_total)
  
  for (i in seq_len(n_total)) {
    
    current_seed <- iteration_seeds[i]
    set.seed(current_seed)
    
    idx_explore <- sample(seq_len(n_rows), size = split_size, replace = FALSE)
    idx_confirm <- setdiff(seq_len(n_rows), idx_explore)
    
    split_results[[i]] <- list(
      iteration    = i,
      seed         = current_seed,
      exploratory  = idx_explore,
      confirmatory = idx_confirm
    )
  }
  
  split_results
}

# ------------------------------------------------------------------------------
# Check EFA Suitability & Extract Factor Assignments / Loadings
# ------------------------------------------------------------------------------

#' Assess Data Suitability and Determine Factors for EFA
#'
#' Calculates KMO, Bartlett's Test, and determines optimal factors.
#' Runs FA to assign items to factors for both Kaiser and Parallel suggestions,
#' and returns per-item loadings and eigenvalues.
#'
#' @param data_subset The dataframe containing the exploratory data.
#' @param items Vector of item names to assess.
#' @param n_iter Number of iterations for parallel analysis simulation.
#' @param loading_cutoff Minimum loading required to assign an item to a factor.
#'
#' @return A list with a Suitability_Summary element.
check_efa_suitability <- function(data_subset, 
                                  items = paste0("SSQ", sprintf("%02d", 1:14)),
                                  n_iter = 20,
                                  loading_cutoff = 0.4) {
  
  missing_items <- setdiff(items, names(data_subset))
  if (length(missing_items) > 0) {
    stop(
      "The following items are missing from the dataset: ",
      paste(missing_items, collapse = ", ")
    )
  }
  
  dt_items <- data_subset[, items]
  
  dt_items[] <- lapply(dt_items, function(x) {
    if (is.factor(x)) {
      as.numeric(x) - 1
    } else if (is.character(x)) {
      if (all(na.omit(x) %in% c("Yes", "No"))) {
        as.numeric(x == "Yes")
      } else {
        as.numeric(x)
      }
    } else {
      x
    }
  })
  
  col_vars <- sapply(dt_items, var, na.rm = TRUE)
  zero_var_cols <- names(col_vars)[!is.na(col_vars) & col_vars == 0]
  
  if (length(zero_var_cols) > 0) {
    return(list(
      Suitability_Summary = list(
        KMO_Overall           = NA_real_,
        Bartlett_p_value      = NA_real_,
        Is_Suitable           = FALSE,
        Zero_Var_Cols         = paste(zero_var_cols, collapse = ", "),
        k_factor_f            = NA_integer_,
        p_factor_f            = NA_integer_,
        factor_items_k        = list(),
        factor_items_p        = list(),
        factor_loadings_k     = NULL,
        factor_loadings_p     = NULL,
        eigenvalues_actual    = NULL,
        eigenvalues_simulated = NULL
      )
    ))
  }
  
  is_binary  <- all(sapply(dt_items, max, na.rm = TRUE) == 1)
  cor_method <- if (is_binary) "tet" else "poly"
  
  kmo_res      <- KMO(dt_items)
  bartlett_res <- cortest.bartlett(dt_items)
  
  cor_mat <- tryCatch({
    suppressMessages(suppressWarnings({
      if (is_binary) {
        psych::tetrachoric(dt_items)$rho
      } else {
        psych::polychoric(dt_items)$rho
      }
    }))
  }, error = function(e) NULL)
  
  if (!is.null(cor_mat)) {
    eigen_values <- eigen(cor_mat)$values
    n_kaiser     <- sum(eigen_values > 1.0)
  } else {
    eigen_values <- NULL
    n_kaiser     <- NA_integer_
  }
  
  pa_res <- tryCatch({
    capture.output(
      res <- suppressMessages(suppressWarnings(
        psych::fa.parallel(
          dt_items,
          fm          = "pa",
          fa          = "fa",
          cor         = cor_method,
          n.iter      = n_iter,
          plot        = FALSE,
          show.legend = FALSE
        )
      )),
      file = nullfile()
    )
    res
  }, error = function(e) NULL)
  
  n_parallel <- if (!is.null(pa_res)) pa_res$nfact else NA_integer_
  
  eigenvalues_actual    <- eigen_values
  eigenvalues_simulated <- if (!is.null(pa_res) && !is.null(pa_res$fa.sim)) pa_res$fa.sim else NULL
  
  is_phq_set <- all(grepl("^PHQ", items))
  if (is_phq_set) {
    if (!is.na(n_kaiser))   n_kaiser   <- min(n_kaiser, 3L)
    if (!is.na(n_parallel)) n_parallel <- min(n_parallel, 3L)
  }
  
  get_factor_solution <- function(n_factors, r_mat, n_obs, cutoff) {
    if (is.na(n_factors) || n_factors < 1L || is.null(r_mat)) {
      return(list(
        factor_items  = list(),
        item_loadings = NULL
      ))
    }
    
    model <- tryCatch({
      capture.output(
        res <- suppressMessages(suppressWarnings(
          psych::fa(
            r        = r_mat,
            nfactors = n_factors,
            n.obs    = n_obs,
            rotate   = "promax",
            fm       = "pa",
            warnings = FALSE
          )
        )),
        file = nullfile()
      )
      res
    }, error = function(e) NULL)
    
    if (is.null(model)) {
      return(list(
        factor_items  = list(),
        item_loadings = NULL
      ))
    }
    
    loads_mat <- as.matrix(model$loadings)
    
    max_loads   <- apply(abs(loads_mat), 1, max)
    assignments <- apply(abs(loads_mat), 1, which.max)
    assignments[max_loads < cutoff] <- NA_integer_
    
    raw_groups <- vector("list", n_factors)
    for (f in seq_len(n_factors)) {
      raw_groups[[f]] <- names(assignments)[!is.na(assignments) & assignments == f]
    }
    
    anchor_1_item <- "SSQ14"
    anchor_2_item <- "SSQ01"
    
    all_assigned <- unlist(raw_groups)
    if (length(all_assigned) > 0 && any(grepl("PHQ", all_assigned))) {
      anchor_1_item <- "PHQ909"
      anchor_2_item <- "PHQ904"
    }
    
    idx_1 <- 0L
    idx_2 <- 0L
    for (f in seq_len(n_factors)) {
      if (anchor_1_item %in% raw_groups[[f]]) idx_1 <- f
      if (anchor_2_item %in% raw_groups[[f]]) idx_2 <- f
    }
    
    ordered_indices <- integer(0)
    if (idx_1 > 0L) ordered_indices <- c(ordered_indices, idx_1)
    if (idx_2 > 0L && !(idx_2 %in% ordered_indices)) ordered_indices <- c(ordered_indices, idx_2)
    
    all_indices <- seq_len(n_factors)
    remaining   <- setdiff(all_indices, ordered_indices)
    final_order <- c(ordered_indices, remaining)
    
    factor_items_list <- list()
    for (i in seq_along(final_order)) {
      old_f <- final_order[i]
      factor_items_list[[as.character(i)]] <- paste(raw_groups[[old_f]], collapse = ", ")
    }
    
    loads_reordered <- loads_mat[, final_order, drop = FALSE]
    
    new_assignments <- assignments
    valid_idx <- !is.na(assignments)
    new_assignments[valid_idx] <- match(assignments[valid_idx], final_order)
    
    item_names <- rownames(loads_reordered)
    df_loadings <- data.frame(
      Item    = item_names,
      Factor  = new_assignments,
      Loading = NA_real_,
      stringsAsFactors = FALSE
    )
    idx_assigned <- which(!is.na(df_loadings$Factor))
    if (length(idx_assigned) > 0) {
      df_loadings$Loading[idx_assigned] <-
        loads_reordered[cbind(idx_assigned, df_loadings$Factor[idx_assigned])]
    }
    df_loadings <- df_loadings[!is.na(df_loadings$Factor), , drop = FALSE]
    
    list(
      factor_items  = factor_items_list,
      item_loadings = df_loadings
    )
  }
  
  k_sol <- get_factor_solution(n_kaiser,   cor_mat, nrow(dt_items), loading_cutoff)
  p_sol <- get_factor_solution(n_parallel, cor_mat, nrow(dt_items), loading_cutoff)
  
  list(
    Suitability_Summary = list(
      KMO_Overall           = kmo_res$MSA,
      Bartlett_p_value      = bartlett_res$p.value,
      Is_Suitable           = (kmo_res$MSA > 0.6 &&
                                 bartlett_res$p.value < 0.05 &&
                                 !is.na(n_kaiser) && !is.na(n_parallel)),
      Zero_Var_Cols         = NA_character_,
      k_factor_f            = n_kaiser,
      p_factor_f            = n_parallel,
      factor_items_k        = k_sol$factor_items,
      factor_items_p        = p_sol$factor_items,
      factor_loadings_k     = k_sol$item_loadings,
      factor_loadings_p     = p_sol$item_loadings,
      eigenvalues_actual    = eigenvalues_actual,
      eigenvalues_simulated = eigenvalues_simulated
    )
  )
}
