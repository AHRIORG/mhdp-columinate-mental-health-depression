if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it before running this script.")
}
if (!requireNamespace("mirt", quietly = TRUE)) {
  stop("Package 'mirt' is required. Install it before running this script.")
}

library(here)
library(mirt)
score_new_cohort <- function(new_data, scoring_engine) {
  missing <- setdiff(scoring_engine$items, names(new_data))
  if (length(missing) > 0) stop("New data missing items: ", paste(missing, collapse = ", "))
  
  dat_score <- new_data[, scoring_engine$items, drop = FALSE]
  dat_score[] <- lapply(dat_score, function(x) if (is.factor(x)) as.numeric(x) - 1 else x)
  
  keep_rows <- rowSums(!is.na(dat_score)) > 0
  if (sum(!keep_rows) > 0) warning("Removed ", sum(!keep_rows), " rows with all-missing SSQ items.")
  
  dat_score_clean <- dat_score[keep_rows, , drop = FALSE]
  
  n_items <- ncol(dat_score_clean)
  item_types_vec <- rep(scoring_engine$model_type, n_items)
  n_factors <- if (!is.null(scoring_engine$n_factors)) scoring_engine$n_factors else 1
  
  template_pars <- mirt(dat_score_clean, n_factors, itemtype = item_types_vec, pars = "values", verbose = FALSE)
  saved_pars <- scoring_engine$parameters
  
  for (i in seq_len(nrow(template_pars))) {
    itm <- template_pars$item[i]
    pnm <- template_pars$name[i]
    if (itm == "GROUP") next
    match_row <- saved_pars[saved_pars$item == itm & saved_pars$name == pnm, , drop = FALSE]
    if (nrow(match_row) == 1) {
      template_pars$value[i] <- match_row$value
      template_pars$est[i] <- FALSE
    }
  }
  
  mod_fixed <- mirt(dat_score_clean, n_factors, itemtype = item_types_vec,
                    pars = template_pars, verbose = FALSE, calcNull = FALSE)
  
  # Use QMC for high-dimensional scoring models
  nfact_fixed <- tryCatch(mod_fixed@Model$nfact, error = function(e) NA_integer_)
  use_qmc <- !is.na(nfact_fixed) && nfact_fixed >= 3
  theta <- fscores(mod_fixed, method = "EAP", full.scores = TRUE, QMC = use_qmc)[, 1]
  
  out <- data.frame(
    Theta_Harmonized = NA_real_,
    Depression_Binary = NA_integer_,
    PHQ_Expected_Sum = NA_real_,
    PHQ_Prob_GE10 = NA_real_
  )
  out <- out[rep(1, nrow(new_data)), , drop = FALSE]
  out$Theta_Harmonized[keep_rows] <- theta
  
  # Map Theta -> PHQ expected sum and P(PHQ>=10) using stored theta map
  if (!is.null(scoring_engine$phq_theta_map)) {
    tm <- scoring_engine$phq_theta_map
    out$PHQ_Expected_Sum[keep_rows] <- approx(tm$Theta, tm$PHQ_Expected_Sum, xout = theta, rule = 2)$y
    out$PHQ_Prob_GE10[keep_rows] <- approx(tm$Theta, tm$PHQ_Prob_GE10, xout = theta, rule = 2)$y
  }
  
  df_tmp <- new_data
  df_tmp$Theta_Harmonized <- out$Theta_Harmonized
  df_tmp <- .ensure_groups(df_tmp)
  
  cutoff_table <- scoring_engine$cutoff_table
  apply_mode <- scoring_engine$cutoff_apply_mode
  hierarchy <- scoring_engine$cutoff_apply_hierarchy
  
  df_tmp <- .apply_cutoffs(df_tmp, cutoff_table, apply_mode = apply_mode, hierarchy = hierarchy)
  out$Depression_Binary <- df_tmp$Depression_Harmonized
  
  out
}

.ensure_groups <- function(df) {
  if (!"AGEGRP" %in% names(df) && "AGE" %in% names(df)) {
    df$AGEGRP <- ifelse(df$AGE < 20, "Age_17_19", "Age_20_24")
  }
  if (!"SEX_AGEGRP" %in% names(df)) {
    sx <- if ("SEX" %in% names(df)) as.character(df$SEX) else NA_character_
    ag <- if ("AGEGRP" %in% names(df)) as.character(df$AGEGRP) else NA_character_
    df$SEX_AGEGRP <- ifelse(is.na(sx) | is.na(ag), NA_character_, paste0(sx, "__", ag))
  }
  df
}

.apply_cutoffs <- function(df, cutoff_table,
                           apply_mode = c("none", "SEX", "AGEGRP", "SEX_AGEGRP", "hierarchical"),
                           hierarchy = c("SEX_AGEGRP", "SEX", "AGEGRP", "none")) {
  
  apply_mode <- match.arg(apply_mode)
  
  df$Depression_Harmonized <- NA_integer_
  
  overall_cut <- cutoff_table %>%
    filter(Grouping == "none", Group_Level == "Overall") %>%
    slice(1) %>%
    pull(Cutoff_Theta)
  if (length(overall_cut) == 0 || is.na(overall_cut)) overall_cut <- NA_real_
  
  get_cut <- function(g, lv) {
    row <- cutoff_table %>% filter(Grouping == g, Group_Level == lv) %>% slice(1)
    if (nrow(row) == 0) return(NA_real_)
    as.numeric(row$Cutoff_Theta[1])
  }
  
  if (apply_mode == "none") {
    if (!is.na(overall_cut)) {
      df$Depression_Harmonized <- ifelse(df$Theta_Harmonized >= overall_cut, 1L, 0L)
    }
    return(df)
  }
  
  if (apply_mode == "hierarchical") {
    hierarchy <- as.character(hierarchy)
    
    for (i in seq_len(nrow(df))) {
      th <- df$Theta_Harmonized[i]
      if (is.na(th)) next
      
      cut <- NA_real_
      for (g in hierarchy) {
        if (g == "none") {
          cut <- overall_cut
        } else {
          if (!g %in% names(df)) next
          lv <- as.character(df[[g]][i])
          if (!is.na(lv)) cut <- get_cut(g, lv)
        }
        if (!is.na(cut)) break
      }
      
      if (!is.na(cut)) df$Depression_Harmonized[i] <- ifelse(th >= cut, 1L, 0L)
    }
    
    return(df)
  }
  
  g <- apply_mode
  if (!g %in% names(df)) {
    if (!is.na(overall_cut)) {
      df$Depression_Harmonized <- ifelse(df$Theta_Harmonized >= overall_cut, 1L, 0L)
    }
    return(df)
  }
  
  for (lv in unique(na.omit(df[[g]]))) {
    cut <- get_cut(g, as.character(lv))
    if (is.na(cut)) cut <- overall_cut
    idx <- which(df[[g]] == lv & !is.na(df$Theta_Harmonized))
    if (!is.na(cut) && length(idx) > 0) {
      df$Depression_Harmonized[idx] <- ifelse(df$Theta_Harmonized[idx] >= cut, 1L, 0L)
    }
  }
  
  if (!is.na(overall_cut)) {
    idx <- which(is.na(df$Depression_Harmonized) & !is.na(df$Theta_Harmonized))
    if (length(idx) > 0) df$Depression_Harmonized[idx] <- ifelse(df$Theta_Harmonized[idx] >= overall_cut, 1L, 0L)
  }
  
  df
}
