############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 02_invariance_roc_helpers.R - Invariance & ROC helpers
############################################################

get_invariance_deltas <- function(data, model, group_var, items) {
  if (!group_var %in% names(data) ||
      length(unique(na.omit(data[[group_var]]))) < 2) return(NULL)
  
  min_n <- min(table(data[[group_var]]))
  if (min_n < 50) return(NULL)
  
  tryCatch({
    fit_conf <- cfa(model, data = data, group = group_var,
                    estimator = "WLSMV", ordered = items)
    fit_scal <- cfa(model, data = data, group = group_var,
                    group.equal = c("loadings", "thresholds"),
                    estimator = "WLSMV", ordered = items)
    
    fm_conf <- fitMeasures(fit_conf, c("cfi", "rmsea"))
    fm_scal <- fitMeasures(fit_scal, c("cfi", "rmsea"))
    
    data.frame(
      Group       = group_var,
      Delta_CFI   = fm_conf["cfi"] - fm_scal["cfi"],
      Delta_RMSEA = fm_scal["rmsea"] - fm_conf["rmsea"]
    )
  }, error = function(e) NULL)
}

get_roc_metrics <- function(data, seed, group_label, items) {
  if (nrow(data) < 20 || length(unique(data$Depressed)) < 2) return(NULL)
  
  tryCatch({
    item_data <- data[, items]
    item_data[] <- lapply(item_data, function(x)
      if (is.factor(x)) as.numeric(x) - 1 else x)
    sum_score <- rowSums(item_data, na.rm = TRUE)
    
    roc_obj <- roc(data$Depressed, sum_score, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    coords_res <- coords(
      roc_obj, "best",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    )
    if (nrow(coords_res) > 1) coords_res <- coords_res[1, ]
    
    thresh <- coords_res$threshold
    pred   <- ifelse(sum_score >= thresh, 1, 0)
    obs    <- data$Depressed
    
    tbl <- table(factor(pred, levels = c(0, 1)),
                 factor(obs, levels = c(0, 1)))
    
    po <- sum(diag(tbl)) / sum(tbl)
    pe <- ((sum(tbl[2, ]) * sum(tbl[, 2])) +
             (sum(tbl[1, ]) * sum(tbl[, 1]))) / sum(tbl)^2
    kappa_val <- (po - pe) / (1 - pe)
    
    data.frame(
      Seed           = seed,
      Group          = group_label,
      AUC            = auc_val,
      Optimal_Cutoff = coords_res$threshold,
      Sensitivity    = coords_res$sensitivity,
      Specificity    = coords_res$specificity,
      Kappa          = kappa_val
    )
  }, error = function(e) NULL)
}

get_roc_curve_data <- function(data, seed, group_label, items) {
  if (nrow(data) < 20 || length(unique(data$Depressed)) < 2) return(NULL)
  
  tryCatch({
    item_data <- data[, items]
    item_data[] <- lapply(item_data, function(x)
      if (is.factor(x)) as.numeric(x) - 1 else x)
    sum_score <- rowSums(item_data, na.rm = TRUE)
    
    roc_obj <- roc(data$Depressed, sum_score, quiet = TRUE)
    
    specs <- seq(0, 1, 0.01)
    senss <- coords(
      roc_obj, x = specs, input = "specificity",
      ret = "sensitivity", transpose = FALSE
    )
    
    data.frame(
      Seed        = seed,
      Group       = group_label,
      Specificity = specs,
      Sensitivity = as.numeric(unlist(senss))
    )
  }, error = function(e) {
    message(paste("Error in get_roc_curve_data for Group:", group_label, "-", e$message))
    NULL
  })
}
