############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 01_efa_consensus.R - EFA bootstrap consensus & CFA syntax
############################################################

calculate_cooccurrence_matrix <- function(subset_df, method_prefix, all_items) {
  n_iter <- nrow(subset_df)
  if (n_iter == 0) return(NULL)
  
  mat <- matrix(0, nrow = length(all_items), ncol = length(all_items))
  rownames(mat) <- all_items
  colnames(mat) <- all_items
  
  factor_cols <- grep(paste0("^", method_prefix, "[0-9]+$"),
                      names(subset_df), value = TRUE)
  
  for (i in seq_len(n_iter)) {
    row_data <- subset_df[i, ]
    for (col in factor_cols) {
      item_string <- row_data[[col]]
      if (is.na(item_string) || item_string == "") next
      
      items_in_factor <- trimws(unlist(strsplit(item_string, ",")))
      valid_items     <- intersect(items_in_factor, all_items)
      
      if (length(valid_items) >= 2) {
        pairs <- t(combn(valid_items, 2))
        for (p in seq_len(nrow(pairs))) {
          itm1 <- pairs[p, 1]
          itm2 <- pairs[p, 2]
          mat[itm1, itm2] <- mat[itm1, itm2] + 1
          mat[itm2, itm1] <- mat[itm2, itm1] + 1
        }
      }
      for (itm in valid_items) {
        mat[itm, itm] <- mat[itm, itm] + 1
      }
    }
  }
  mat / n_iter
}

load_bootstrap_results <- function(dataset_label, target_set_name) {
  safe_set_name <- gsub("[\n ]+", "_", target_set_name)
  safe_ds_label <- gsub("[\n ]+", "_", dataset_label)
  boot_file <- paste0("boot_1000_", safe_ds_label, "_", safe_set_name, "_results.rds")
  boot_path <- here("statistical_analysis/output/objects", boot_file)
  
  if (!file.exists(boot_path)) {
    stop(paste("EFA Bootstrap results not found at:", boot_path))
  }
  
  list(
    boot_results   = readRDS(boot_path),
    safe_set_name  = safe_set_name,
    safe_ds_label  = safe_ds_label,
    boot_path      = boot_path
  )
}

select_factor_solution <- function(boot_results, manual_method, manual_n_factors) {
  if (is.null(manual_method) || is.null(manual_n_factors)) {
    k_counts <- table(boot_results$k_factor_f)
    p_counts <- table(boot_results$p_factor_f)
    
    k_win_n <- as.numeric(names(which.max(k_counts)))
    p_win_n <- as.numeric(names(which.max(p_counts)))
    
    if (max(k_counts) >= max(p_counts)) {
      target_method <- "Kaiser"
      target_n      <- k_win_n
    } else {
      target_method <- "Parallel"
      target_n      <- p_win_n
    }
    message(sprintf("Initial Consensus: %s suggested %d Factors", target_method, target_n))
  } else {
    target_method <- manual_method
    target_n      <- manual_n_factors
    message(sprintf("Manual Override: %s with %d Factors", target_method, target_n))
  }
  
  list(target_method = target_method, target_n = target_n)
}

extract_winning_structure <- function(boot_results, target_method, target_n) {
  target_col_prefix <- if (target_method == "Kaiser") "k_factor_" else "p_factor_"
  cols_to_check     <- paste0(target_col_prefix, 1:target_n)
  
  valid_runs <- boot_results %>%
    filter(
      Is_Suitable == TRUE,
      if (target_method == "Kaiser") k_factor_f == target_n else p_factor_f == target_n
    )
  
  signatures <- apply(valid_runs[, cols_to_check, drop = FALSE], 1, paste, collapse = "|")
  winning_sig <- names(which.max(table(signatures)))
  
  factor_defs_raw <- unlist(strsplit(winning_sig, "\\|"))
  factor_list <- lapply(factor_defs_raw, function(x) trimws(unlist(strsplit(x, ","))))
  
  list(
    factor_list       = factor_list,
    valid_runs        = valid_runs,
    signatures        = signatures,
    winning_sig       = winning_sig,
    target_col_prefix = target_col_prefix
  )
}

merge_factors_manual <- function(factor_list, manual_merge_pairs) {
  if (is.null(manual_merge_pairs)) return(factor_list)
  
  final_factor_list <- factor_list
  message("\nExecuting Manual Merge(s)...")
  
  for (pair in manual_merge_pairs) {
    source_idx <- pair[1]
    target_idx <- pair[2]
    
    if (source_idx > length(final_factor_list) || target_idx > length(final_factor_list)) {
      warning(sprintf("Skipping invalid merge pair: %d -> %d (Indices out of bounds)",
                      source_idx, target_idx))
      next
    }
    
    if (!is.null(final_factor_list[[source_idx]])) {
      message(sprintf(
        "  Merging Factor %d (Items: %s) into Factor %d",
        source_idx, paste(final_factor_list[[source_idx]], collapse = ","), target_idx
      ))
      final_factor_list[[target_idx]] <- c(final_factor_list[[target_idx]],
                                           final_factor_list[[source_idx]])
      final_factor_list[[source_idx]] <- NULL
    }
  }
  
  final_factor_list <- final_factor_list[!sapply(final_factor_list, is.null)]
  message(sprintf("  Manual Merge Complete. New Factor Count: %d", length(final_factor_list)))
  final_factor_list
}

merge_factors_programmatic <- function(factor_list, valid_runs,
                                       target_col_prefix, target_items,
                                       min_items_per_factor) {
  needs_merge <- any(sapply(factor_list, length) < min_items_per_factor)
  if (!needs_merge) {
    message(sprintf("\nAll factors have >= %d items. No merge required.", min_items_per_factor))
    return(factor_list)
  }
  
  message(sprintf(
    "\nDetected factors with <%d items. Initiating programmatic merge...",
    min_items_per_factor
  ))
  
  cooc_mat <- calculate_cooccurrence_matrix(valid_runs, target_col_prefix, target_items)
  final_factor_list <- factor_list
  small_indices <- which(sapply(final_factor_list, length) < min_items_per_factor)
  
  for (idx in small_indices) {
    small_items <- final_factor_list[[idx]]
    other_indices <- setdiff(seq_along(final_factor_list), small_indices)
    if (length(other_indices) == 0) next
    
    best_target_idx <- NA
    max_avg_cooc    <- -1
    
    for (target_idx in other_indices) {
      target_items_cand <- final_factor_list[[target_idx]]
      sub_mat <- cooc_mat[small_items, target_items_cand, drop = FALSE]
      avg_cooc <- mean(sub_mat)
      
      if (avg_cooc > max_avg_cooc) {
        max_avg_cooc    <- avg_cooc
        best_target_idx <- target_idx
      }
    }
    
    if (!is.na(best_target_idx)) {
      message(sprintf(
        "  Merging Small Factor %d (Items: %s) into Factor %d (Avg Co-occurrence: %.2f)",
        idx, paste(small_items, collapse = ","), best_target_idx, max_avg_cooc
      ))
      final_factor_list[[best_target_idx]] <- c(final_factor_list[[best_target_idx]],
                                                small_items)
      final_factor_list[[idx]] <- NULL
    }
  }
  
  final_factor_list <- final_factor_list[!sapply(final_factor_list, is.null)]
  message(sprintf("  Merge Complete. New Factor Count: %d", length(final_factor_list)))
  final_factor_list
}

build_cfa_syntax <- function(factor_list) {
  parts <- vapply(seq_along(factor_list), function(i) {
    paste0("F", i, " =~ ", paste(factor_list[[i]], collapse = " + "))
  }, character(1))
  
  syntax <- paste(parts, collapse = "\n")
  message("\n--- Final CFA Syntax ---")
  cat(syntax)
  message("\n------------------------")
  syntax
}

select_mc_seeds <- function(valid_runs, signatures, winning_sig, n_mc_samples) {
  target_seeds <- valid_runs$Seed[signatures == winning_sig]
  if (length(target_seeds) > n_mc_samples) {
    set.seed(123)
    sample(target_seeds, n_mc_samples)
  } else {
    target_seeds
  }
}
