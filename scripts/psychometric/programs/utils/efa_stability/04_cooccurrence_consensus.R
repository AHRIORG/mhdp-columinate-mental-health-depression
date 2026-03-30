# Psychometric EFA – Co-occurrence & Consensus Model Helpers

# Dependencies: dplyr, tidyr, ggplot2, scales, grid

# ------------------------------------------------------------------------------
# 1. Co-occurrence Matrix
# ------------------------------------------------------------------------------

#' Calculate co-occurrence matrix for items across bootstrap iterations
#'
#' @param subset_df Data frame of bootstrap results filtered to a single factor count.
#' @param method_prefix String, e.g., "k_factor_" or "p_factor_".
#' @param all_items Character vector of all item names.
#'
#' @return Matrix of co-occurrence probabilities.
calculate_cooccurrence_matrix <- function(subset_df, method_prefix, all_items) {
  n_iter <- nrow(subset_df)
  if (n_iter == 0) return(NULL)
  
  mat <- matrix(0, nrow = length(all_items), ncol = length(all_items))
  rownames(mat) <- all_items
  colnames(mat) <- all_items
  
  factor_cols <- grep(paste0("^", method_prefix, "[0-9]+$"), names(subset_df), value = TRUE)
  
  for (i in seq_len(n_iter)) {
    row_data <- subset_df[i, , drop = FALSE]
    for (col in factor_cols) {
      item_string <- row_data[[col]]
      if (is.na(item_string) || item_string == "") next
      items_in_factor <- trimws(unlist(strsplit(item_string, ",")))
      valid_items <- intersect(items_in_factor, all_items)
      
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

# ------------------------------------------------------------------------------
# 2. Item-Factor Stability & Ordering
# ------------------------------------------------------------------------------

#' Calculate item-factor probabilities and item ordering by stability
#'
#' @param subset_df Data frame of bootstrap results (single factor count).
#' @param method_prefix "k_factor_" or "p_factor_".
#' @param winning_count Integer number of factors.
#' @param all_items Character vector of all item names.
#'
#' @return List with probs_mat (items x factors) and sorted_items.
get_sorted_items_by_stability <- function(subset_df, method_prefix, winning_count, all_items) {
  counts_mat <- matrix(0, nrow = length(all_items), ncol = winning_count)
  rownames(counts_mat) <- all_items
  colnames(counts_mat) <- paste("Factor", seq_len(winning_count))
  
  factor_cols <- paste0(method_prefix, seq_len(winning_count))
  n_total <- nrow(subset_df)
  
  for (i in seq_len(n_total)) {
    row_data <- subset_df[i, , drop = FALSE]
    for (f in seq_len(winning_count)) {
      col_name <- factor_cols[f]
      if (!col_name %in% names(row_data)) next
      item_str <- row_data[[col_name]]
      if (is.na(item_str) || item_str == "") next
      items_in_f <- trimws(unlist(strsplit(item_str, ",")))
      valid_items <- intersect(items_in_f, all_items)
      if (length(valid_items) > 0) {
        counts_mat[valid_items, f] <- counts_mat[valid_items, f] + 1
      }
    }
  }
  
  probs_mat <- counts_mat / n_total
  
  primary_factor <- apply(probs_mat, 1, which.max)
  max_probs      <- apply(probs_mat, 1, max)
  
  sort_df <- data.frame(Item = all_items, Primary = primary_factor, MaxProb = max_probs)
  sort_df <- sort_df[order(sort_df$Primary, -sort_df$MaxProb), , drop = FALSE]
  
  list(probs_mat = probs_mat, sorted_items = as.character(sort_df$Item))
}

# ------------------------------------------------------------------------------
# 3. Item-Factor Stability Plot
# ------------------------------------------------------------------------------

#' Plot item-factor stability (stacked membership probabilities)
#'
#' @param probs_mat Matrix of membership probabilities (items x factors).
#' @param item_order Vector of item names in desired plotting order.
#' @param winning_count Number of factors.
#' @param consensus_info List from identify_consensus_model().
#' @param subtitle Plot subtitle.
#'
#' @return ggplot object.
plot_item_factor_stability <- function(probs_mat, item_order, winning_count, consensus_info, subtitle) {
  plot_df <- as.data.frame(probs_mat)
  plot_df$Item <- rownames(plot_df)
  plot_long <- tidyr::pivot_longer(
    plot_df,
    cols      = dplyr::starts_with("Factor"),
    names_to  = "Factor",
    values_to = "Probability"
  )
  
  plot_long$Item <- factor(plot_long$Item, levels = item_order)
  
  n_success <- consensus_info$frequency_n
  n_trials  <- consensus_info$total_runs_in_subset
  prop      <- consensus_info$frequency_pct / 100
  
  se <- sqrt(prop * (1 - prop) / n_trials)
  ci_lower <- max(0, prop - 1.96 * se)
  ci_upper <- min(1, prop + 1.96 * se)
  
  stats_str <- sprintf(
    "Frequency: %d/%d (%.1f%%); 95%% CI: (%.1f%%, %.1f%%)",
    n_success, n_trials, prop * 100, ci_lower * 100, ci_upper * 100
  )
  
  model_str_vec <- character()
  struct <- consensus_info$structure
  factor_names <- paste("Factor", seq_len(winning_count))
  for (f_name in factor_names) {
    if (f_name %in% names(struct)) {
      f_num    <- as.numeric(sub("Factor ", "", f_name))
      items_str <- struct[[f_name]]
      model_str_vec <- c(model_str_vec, sprintf("FACTOR %d: %s", f_num, items_str))
    }
  }
  caption_str <- paste(c(stats_str, "", paste(model_str_vec, collapse = "; ")), collapse = "\n")
  
  ggplot2::ggplot(
    plot_long,
    ggplot2::aes(x = Item, y = Probability, fill = Factor)
  ) +
    ggplot2::geom_col(position = "fill", width = 0.7, alpha = 0.9) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::theme_minimal(base_family = "Arial", base_size = 14) +
    ggplot2::labs(
      title    = NULL,
      subtitle = subtitle,
      x        = NULL,
      y        = "Membership Probability",
      caption  = caption_str
    ) +
    ggplot2::theme(
      legend.position    = "bottom",
      plot.caption       = ggplot2::element_text(hjust = 0, size = 12, family = "Arial"),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 14, family = "Arial")
    )
}

# ------------------------------------------------------------------------------
# 4. Co-occurrence Heatmap
# ------------------------------------------------------------------------------

#' Plot co-occurrence heatmap
#'
#' @param mat Co-occurrence probability matrix.
#' @param subtitle Subtitle for the plot.
#' @param caption Optional caption.
#' @param item_order Optional item order.
#'
#' @return ggplot object.
plot_cooccurrence_heatmap <- function(mat, subtitle, caption = NULL, item_order = NULL) {
  long_dat <- as.data.frame(as.table(mat))
  colnames(long_dat) <- c("Item1", "Item2", "Probability")
  
  sort_levels <- if (!is.null(item_order)) item_order else colnames(mat)
  long_dat$Item1 <- factor(long_dat$Item1, levels = sort_levels)
  long_dat$Item2 <- factor(long_dat$Item2, levels = sort_levels)
  
  long_dat$x_idx <- match(long_dat$Item1, sort_levels)
  long_dat$y_idx <- match(long_dat$Item2, sort_levels)
  
  long_dat <- dplyr::mutate(
    long_dat,
    Type = dplyr::case_when(
      x_idx == y_idx ~ "Diagonal",
      x_idx >  y_idx ~ "Lower",
      x_idx <  y_idx ~ "Upper"
    )
  )
  
  p <- ggplot2::ggplot(long_dat, ggplot2::aes(x = Item1, y = Item2)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_family = "Arial", base_size = 14) +
    ggplot2::labs(
      title    = NULL,
      subtitle = subtitle,
      caption  = caption,
      x        = NULL,
      y        = NULL
    ) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, family = "Arial", size = 14),
      axis.text.y  = ggplot2::element_text(family = "Arial", size = 14),
      plot.caption = ggplot2::element_text(family = "Arial", size = 12),
      legend.position = "top"
    )
  
  p <- p +
    ggplot2::geom_tile(
      data = subset(long_dat, Type %in% c("Upper", "Lower")),
      ggplot2::aes(fill = Probability),
      color = "white"
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("white", "orange", "#A93226"),
      limits = c(0, 1),
      name   = "Co-occurrence\nProbability",
      guide  = ggplot2::guide_colorbar(
        title          = "Co-occurrence\nProbability",
        title.position = "left",
        barwidth       = grid::unit(12, "lines"),
        barheight      = grid::unit(1.25, "lines")
      )
    ) +
    ggplot2::geom_tile(
      data = subset(long_dat, Type == "Diagonal"),
      fill  = "grey90",
      color = "white"
    ) +
    ggplot2::geom_text(
      data = subset(long_dat, Type %in% c("Upper", "Lower")),
      ggplot2::aes(
        label = sprintf("%.1f", Probability * 100),
        color = ifelse(Probability < 0.4, "black", "white")
      ),
      size   = 4,
      family = "Arial"
    ) +
    ggplot2::scale_color_identity()
  
  p
}

# ------------------------------------------------------------------------------
# 5. Consensus Model Identification
# ------------------------------------------------------------------------------

#' Identify consensus model from bootstrap results
#'
#' @param subset_df Data frame of bootstrap results filtered to a single factor count.
#' @param method_prefix "k_factor_" or "p_factor_".
#' @param winning_count Integer, number of factors.
#'
#' @return List with structure and frequency statistics.
identify_consensus_model <- function(subset_df, method_prefix, winning_count) {
  factor_cols  <- paste0(method_prefix, seq_len(winning_count))
  missing_cols <- setdiff(factor_cols, names(subset_df))
  if (length(missing_cols) > 0) {
    return(list(error = "Missing factor columns in data"))
  }
  
  model_signatures <- apply(
    subset_df[, factor_cols, drop = FALSE],
    1,
    function(row) paste(row, collapse = " | ")
  )
  
  sig_counts  <- table(model_signatures)
  winning_sig <- names(which.max(sig_counts))
  freq_n      <- max(sig_counts)
  freq_pct    <- (freq_n / sum(sig_counts)) * 100
  
  factor_contents <- unlist(strsplit(winning_sig, " \\| "))
  consensus_structure <- list()
  for (i in seq_along(factor_contents)) {
    consensus_structure[[paste("Factor", i)]] <- factor_contents[i]
  }
  
  list(
    structure            = consensus_structure,
    frequency_n          = as.integer(freq_n),
    frequency_pct        = as.numeric(freq_pct),
    total_runs_in_subset = as.integer(sum(sig_counts))
  )
}
