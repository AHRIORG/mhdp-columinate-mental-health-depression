# Batch Scoring + Metrics Across Scenarios (Run IDs)
# -----------------------------------------------------------------------------
# PURPOSE
#   Given a vector of scenario run_ids (from your batch outputs), this script:
#   1) Loads each scenario .rds and extracts its scoring engine.
#   2) Scores a NEW dataset in TWO ways (when possible):
#        a) JOINT theta scoring (PHQ+SSQ; uses scenario$obj$models$joint)  -> matches run_id metrics on dt_psychometric
#        b) SSQ-only theta scoring (fixed-parameter engine; SSQ only)      -> portable to SSQ-only cohorts
#   3) Computes classification metrics + 95% CI at the scenario-applied cutoff for BOTH theta sources.
#   4) Computes traditional ROC summaries (AUC + CI) + Youden operating point metrics for BOTH theta sources.
#   5) Computes traditional ROC summaries (AUC + CI) + Youden operating point metrics for SSQ-10 sum score.
#   6) Returns a single tidy dt_metrics_scenario data.frame.
#
# NOTES
#   - If your new dataset does not contain PHQ items (or an outcome variable),
#     ROC and supervised metrics that require truth will be returned as NA.
#   - Joint-theta requires the scenario object to contain models$joint and the new_data to contain PHQ+SSQ items.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
  if (!requireNamespace("mirt", quietly = TRUE)) install.packages("mirt")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
})

library(here)
library(mirt)
library(dplyr)
library(pROC)
library(ggplot2)

# -----------------------------------------------------------------------------
# A. Small utilities
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

.default_psychometric_data_path <- function() {
  override <- Sys.getenv("COLU_PSYCHOMETRIC_RDATA", unset = "")
  if (nzchar(override)) return(override)
  here::here("data/inputs/dt_psychometric.RData")
}

# Coerce ordinal item responses to numeric 0-based values
# - factor: as.numeric(level) - 1
# - character: try as.numeric; if fails, treat as categorical (factor) -> as.numeric-1
# - haven_labelled: best-effort numeric
.coerce_item_0based <- function(x) {
  if (is.factor(x)) return(as.numeric(x) - 1)
  if (is.character(x)) {
    nx <- suppressWarnings(as.numeric(x))
    if (!all(is.na(nx))) return(nx)
    return(as.numeric(factor(x, exclude = NULL)) - 1)
  }
  if (inherits(x, "haven_labelled")) {
    nx <- suppressWarnings(as.numeric(x))
    if (!all(is.na(nx))) return(nx)
  }
  x
}

.wilson_ci <- function(x, n, conf.level = 0.95) {
  if (is.na(x) || is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  denom <- 1 + (z^2 / n)
  center <- (p + (z^2 / (2 * n))) / denom
  half <- (z * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2)))) / denom
  c(max(0, center - half), min(1, center + half))
}

.compute_kappa <- function(y_true, y_pred) {
  ok <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]
  if (length(y_true) == 0) return(NA_real_)
  
  tab <- table(factor(y_pred, levels = c(0, 1)), factor(y_true, levels = c(0, 1)))
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  
  po <- (tab[1, 1] + tab[2, 2]) / n
  pe <- ((sum(tab[1, ]) * sum(tab[, 1])) + (sum(tab[2, ]) * sum(tab[, 2]))) / (n^2)
  if (isTRUE(all.equal(1, pe))) return(NA_real_)
  (po - pe) / (1 - pe)
}

.bootstrap_ci_kappa <- function(y_true, y_pred, conf.level = 0.95, n_boot = 500, seed = 1) {
  ok <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]
  n <- length(y_true)
  if (n < 20) return(c(NA_real_, NA_real_))
  
  set.seed(seed)
  ks <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    ks[b] <- .compute_kappa(y_true[idx], y_pred[idx])
  }
  ks <- ks[!is.na(ks)]
  if (length(ks) < 50) return(c(NA_real_, NA_real_))
  
  alpha <- (1 - conf.level) / 2
  stats::quantile(ks, probs = c(alpha, 1 - alpha), names = FALSE)
}

.compute_metrics_with_ci <- function(y_true, y_pred, conf.level = 0.95, n_boot_kappa = 500) {
  ok <- !is.na(y_true) & !is.na(y_pred)
  y_true <- as.integer(y_true[ok])
  y_pred <- as.integer(y_pred[ok])
  
  if (length(y_true) == 0 || length(unique(y_true)) < 2) {
    return(list(
      Accuracy = NA_real_, Accuracy_LCL = NA_real_, Accuracy_UCL = NA_real_,
      Sensitivity = NA_real_, Sensitivity_LCL = NA_real_, Sensitivity_UCL = NA_real_,
      Specificity = NA_real_, Specificity_LCL = NA_real_, Specificity_UCL = NA_real_,
      PPV = NA_real_, PPV_LCL = NA_real_, PPV_UCL = NA_real_,
      NPV = NA_real_, NPV_LCL = NA_real_, NPV_UCL = NA_real_,
      Cohen_Kappa = NA_real_, Cohen_Kappa_LCL = NA_real_, Cohen_Kappa_UCL = NA_real_,
      Observed_Prevalence = NA_real_, Predicted_Prevalence = NA_real_, Predicted_Positive_N = NA_integer_
    ))
  }
  
  tab <- table(factor(y_pred, levels = c(0, 1)), factor(y_true, levels = c(0, 1)))
  TN <- tab[1, 1]; FN <- tab[1, 2]; FP <- tab[2, 1]; TP <- tab[2, 2]
  n <- sum(tab)
  
  acc <- (TP + TN) / n
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  ppv <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  npv <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
  
  acc_ci <- .wilson_ci(TP + TN, n, conf.level)
  sens_ci <- .wilson_ci(TP, TP + FN, conf.level)
  spec_ci <- .wilson_ci(TN, TN + FP, conf.level)
  ppv_ci <- .wilson_ci(TP, TP + FP, conf.level)
  npv_ci <- .wilson_ci(TN, TN + FN, conf.level)
  
  kappa <- .compute_kappa(y_true, y_pred)
  k_ci <- .bootstrap_ci_kappa(y_true, y_pred, conf.level = conf.level, n_boot = n_boot_kappa)
  
  list(
    Accuracy = acc, Accuracy_LCL = acc_ci[1], Accuracy_UCL = acc_ci[2],
    Sensitivity = sens, Sensitivity_LCL = sens_ci[1], Sensitivity_UCL = sens_ci[2],
    Specificity = spec, Specificity_LCL = spec_ci[1], Specificity_UCL = spec_ci[2],
    PPV = ppv, PPV_LCL = ppv_ci[1], PPV_UCL = ppv_ci[2],
    NPV = npv, NPV_LCL = npv_ci[1], NPV_UCL = npv_ci[2],
    Cohen_Kappa = kappa, Cohen_Kappa_LCL = k_ci[1], Cohen_Kappa_UCL = k_ci[2],
    Observed_Prevalence = mean(y_true == 1, na.rm = TRUE),
    Predicted_Prevalence = mean(y_pred == 1, na.rm = TRUE),
    Predicted_Positive_N = as.integer(sum(y_pred == 1, na.rm = TRUE))
  )
}

# -----------------------------------------------------------------------------
# B. Load scenario objects by run_id
# -----------------------------------------------------------------------------

.find_scenario_file <- function(out_dir, run_id) {
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  patt <- paste0("IRT_batch__", run_id, "\\.rds$")
  cand <- list.files(out_dir, pattern = patt, recursive = TRUE, full.names = TRUE)
  if (length(cand) == 0) {
    patt2 <- paste0(run_id, ".*\\.rds$")
    cand <- list.files(out_dir, pattern = patt2, recursive = TRUE, full.names = TRUE)
  }
  if (length(cand) == 0) return(NA_character_)
  cand[[1]]
}

.load_scenario <- function(out_dir, run_id) {
  f <- .find_scenario_file(out_dir, run_id)
  if (is.na(f)) stop("Scenario file not found for run_id: ", run_id)
  obj <- readRDS(f)
  list(run_id = run_id, file = f, obj = obj)
}

# -----------------------------------------------------------------------------
# C. Scoring helpers (SSQ-only theta OR JOINT theta)
# -----------------------------------------------------------------------------

.add_agegrp_if_needed <- function(df, age_var = "AGE") {
  if (!"AGEGRP" %in% names(df) && age_var %in% names(df)) {
    df$AGEGRP <- ifelse(df[[age_var]] < 20, "Age_17_19", "Age_20_24")
  }
  df
}

.get_row_cutoff <- function(df, scoring_engine, sex_var = "SEX") {
  ct <- scoring_engine$cutoff_table
  if (is.null(ct) || !all(c("Grouping", "Group_Level", "Cutoff_Theta") %in% names(ct))) {
    thr <- scoring_engine$threshold %||% scoring_engine$cutoff_theta %||% NA_real_
    return(rep(as.numeric(thr), nrow(df)))
  }
  
  ct$Cutoff_Theta <- suppressWarnings(as.numeric(ct$Cutoff_Theta))
  ct <- ct[!is.na(ct$Cutoff_Theta), , drop = FALSE]
  
  pull <- function(g, lv) {
    r <- ct[ct$Grouping == g & ct$Group_Level == lv, , drop = FALSE]
    if (nrow(r) == 0) return(NA_real_)
    as.numeric(r$Cutoff_Theta[1])
  }
  
  overall <- pull("none", "Overall")
  mode <- scoring_engine$cutoff_apply_mode %||% "single"
  hier <- scoring_engine$cutoff_apply_hierarchy %||% c("SEX_AGEGRP", "SEX", "AGEGRP", "none")
  
  if (!"SEX_AGEGRP" %in% names(df)) {
    if (sex_var %in% names(df) && "AGEGRP" %in% names(df)) {
      df$SEX_AGEGRP <- paste0(df[[sex_var]], "__", df$AGEGRP)
    }
  }
  
  n <- nrow(df)
  
  if (mode %in% c("none", "single")) {
    return(rep(overall, n))
  }
  
  if (mode %in% c("SEX", "AGEGRP", "SEX_AGEGRP")) {
    g <- mode
    if (!g %in% names(df)) return(rep(overall, n))
    lv <- as.character(df[[g]])
    out <- vapply(lv, function(z) {
      v <- pull(g, z)
      if (is.na(v)) overall else v
    }, numeric(1))
    return(out)
  }
  
  if (mode == "hierarchical") {
    out <- rep(NA_real_, n)
    for (g in hier) {
      if (g == "none") break
      if (!g %in% names(df)) next
      lv <- as.character(df[[g]])
      m <- is.na(out)
      if (any(m)) {
        out[m] <- vapply(lv[m], function(z) {
          v <- pull(g, z)
          if (is.na(v)) NA_real_ else v
        }, numeric(1))
      }
    }
    out[is.na(out)] <- overall
    return(out)
  }
  
  rep(overall, n)
}

.score_theta_ssq <- function(new_data, scoring_engine, method = "EAP") {
  items <- scoring_engine$items
  missing <- setdiff(items, names(new_data))
  if (length(missing) > 0) stop("New data missing SSQ items: ", paste(missing, collapse = ", "))
  
  dat <- new_data[, items, drop = FALSE]
  dat[] <- lapply(dat, .coerce_item_0based)
  
  keep <- rowSums(!is.na(dat)) > 0
  if (!any(keep)) {
    out <- rep(NA_real_, nrow(new_data))
    return(list(theta = out, keep = keep))
  }
  
  dat2 <- dat[keep, , drop = FALSE]
  
  n_items <- ncol(dat2)
  item_types <- rep(scoring_engine$model_type %||% "graded", n_items)
  n_factors <- scoring_engine$n_factors %||% 1
  
  template <- mirt(dat2, n_factors, itemtype = item_types, pars = "values", verbose = FALSE)
  saved <- scoring_engine$parameters
  
  for (i in seq_len(nrow(template))) {
    itm <- template$item[i]
    nm <- template$name[i]
    if (itm == "GROUP") next
    m <- saved[saved$item == itm & saved$name == nm, , drop = FALSE]
    if (nrow(m) == 1) {
      template$value[i] <- m$value
      template$est[i] <- FALSE
    }
  }
  
  mod_fixed <- mirt(dat2, n_factors, itemtype = item_types, pars = template, verbose = FALSE, calcNull = FALSE)
  
  QMC_flag <- isTRUE(n_factors > 1)
  sc <- fscores(mod_fixed, method = method, full.scores = TRUE, QMC = QMC_flag)
  
  theta_out <- rep(NA_real_, nrow(new_data))
  theta_out[keep] <- sc[, 1]
  
  list(theta = theta_out, keep = keep)
}

.score_theta_joint <- function(new_data,
                               mod_joint,
                               ssq_items,
                               phq_items,
                               method = "EAP") {
  if (is.null(mod_joint)) stop("Scenario does not contain a joint model.")
  
  missing <- setdiff(c(phq_items, ssq_items), names(new_data))
  if (length(missing) > 0) stop("New data missing joint items: ", paste(missing, collapse = ", "))
  
  dat <- new_data[, c(phq_items, ssq_items), drop = FALSE]
  dat[] <- lapply(dat, .coerce_item_0based)
  
  valid_joint <- rowSums(!is.na(dat)) > 0
  if (!any(valid_joint)) {
    out <- rep(NA_real_, nrow(new_data))
    return(list(theta = out, keep = valid_joint))
  }
  
  dat2 <- dat[valid_joint, , drop = FALSE]
  n_factors <- mod_joint@Model$nfact
  QMC_flag <- isTRUE(n_factors > 1)
  
  sc <- fscores(mod_joint, method = method, response.pattern = dat2, full.scores = TRUE, QMC = QMC_flag)
  
  theta_out <- rep(NA_real_, nrow(new_data))
  theta_out[valid_joint] <- sc[, 1]
  
  list(theta = theta_out, keep = valid_joint)
}

score_new_dataset_for_scenario <- function(new_data,
                                           scoring_engine,
                                           theta_source = c("joint", "ssq_engine"),
                                           scenario_obj = NULL,
                                           phq_items = paste0("PHQ9", sprintf("%02d", 1:9)),
                                           sex_var = "SEX",
                                           age_var = "AGE") {
  theta_source <- match.arg(theta_source)
  
  df <- as.data.frame(new_data)
  df <- .add_agegrp_if_needed(df, age_var = age_var)
  
  ssq_items <- scoring_engine$items
  df[, ssq_items] <- lapply(df[, ssq_items, drop = FALSE], .coerce_item_0based)
  df$Sum_SSQ <- rowSums(df[, ssq_items, drop = FALSE], na.rm = TRUE)
  
  if (theta_source == "joint") {
    if (is.null(scenario_obj) || is.null(scenario_obj$models) || is.null(scenario_obj$models$joint)) {
      stop("theta_source='joint' requested but scenario_obj$models$joint is missing")
    }
    if (!all(phq_items %in% names(df))) {
      stop("theta_source='joint' requested but PHQ items are not present in new_data")
    }
    df[, phq_items] <- lapply(df[, phq_items, drop = FALSE], .coerce_item_0based)
    
    theta_res <- .score_theta_joint(
      new_data = df,
      mod_joint = scenario_obj$models$joint,
      ssq_items = ssq_items,
      phq_items = phq_items,
      method = "EAP"
    )
  } else {
    theta_res <- .score_theta_ssq(df, scoring_engine)
  }
  
  df$Theta_Harmonized <- theta_res$theta
  df$Theta_Source_Used <- theta_source
  
  df$Cutoff_Theta_Applied <- .get_row_cutoff(df, scoring_engine, sex_var = sex_var)
  df$Depression_Binary <- as.integer(df$Theta_Harmonized >= df$Cutoff_Theta_Applied)
  
  df
}

# -----------------------------------------------------------------------------
# D. Reference outcome for supervised metrics/ROC
# -----------------------------------------------------------------------------

.compute_reference_outcome <- function(df,
                                       outcome_var = NULL,
                                       phq_items = paste0("PHQ9", sprintf("%02d", 1:9)),
                                       phq_sum_cut = 10) {
  if (!is.null(outcome_var) && outcome_var %in% names(df)) {
    return(as.integer(df[[outcome_var]]))
  }
  
  if (all(phq_items %in% names(df))) {
    phq_dat <- df[, phq_items, drop = FALSE]
    phq_dat[] <- lapply(phq_dat, .coerce_item_0based)
    df$Sum_PHQ <- rowSums(phq_dat, na.rm = TRUE)
    return(as.integer(df$Sum_PHQ >= phq_sum_cut))
  }
  
  rep(NA_integer_, nrow(df))
}

# -----------------------------------------------------------------------------
# E. Traditional ROC helpers
# -----------------------------------------------------------------------------

.roc_auc_ci <- function(y_true, predictor, conf.level = 0.95) {
  ok <- !is.na(y_true) & !is.na(predictor)
  y <- y_true[ok]
  x <- predictor[ok]
  if (length(unique(y)) < 2 || length(y) < 20) {
    return(list(auc = NA_real_, lcl = NA_real_, ucl = NA_real_, roc = NULL))
  }
  
  r <- pROC::roc(response = y, predictor = x, quiet = TRUE, direction = "auto")
  a <- as.numeric(pROC::auc(r))
  ci <- tryCatch(as.numeric(pROC::ci.auc(r, conf.level = conf.level)), error = function(e) c(NA, NA, NA))
  list(auc = a, lcl = ci[1], ucl = ci[3], roc = r)
}

.roc_threshold_youden <- function(roc_obj) {
  if (is.null(roc_obj)) return(NA_real_)
  cd <- pROC::coords(roc_obj, x = "best", best.method = "youden",
                     ret = c("threshold"), transpose = FALSE)
  suppressWarnings(as.numeric(cd["threshold"]))
}

# -----------------------------------------------------------------------------
# F. Main: score across run_ids + compute metrics
# -----------------------------------------------------------------------------

score_dataset_across_runids <- function(run_ids,
                                        new_data,
                                        out_dir = here::here("statistical_analysis/output/objects/irt"),
                                        outcome_var = NULL,
                                        phq_items = paste0("PHQ9", sprintf("%02d", 1:9)),
                                        phq_sum_cut = 10,
                                        conf.level = 0.95,
                                        n_boot_kappa = 500) {
  
  run_ids <- unique(as.character(run_ids))
  stopifnot(length(run_ids) > 0)
  
  all_rows <- list()
  
  add_row <- function(rid, file, approach, predictor, threshold, auc, auc_l, auc_u, met) {
    data.frame(
      Run_ID = rid,
      Scenario_File = file,
      Approach = approach,
      Predictor = predictor,
      Threshold = threshold,
      AUC = auc,
      AUC_CI_Low = auc_l,
      AUC_CI_High = auc_u,
      Accuracy = met$Accuracy,
      Accuracy_LCL = met$Accuracy_LCL,
      Accuracy_UCL = met$Accuracy_UCL,
      Sensitivity = met$Sensitivity,
      Sensitivity_LCL = met$Sensitivity_LCL,
      Sensitivity_UCL = met$Sensitivity_UCL,
      Specificity = met$Specificity,
      Specificity_LCL = met$Specificity_LCL,
      Specificity_UCL = met$Specificity_UCL,
      PPV = met$PPV,
      PPV_LCL = met$PPV_LCL,
      PPV_UCL = met$PPV_UCL,
      NPV = met$NPV,
      NPV_LCL = met$NPV_LCL,
      NPV_UCL = met$NPV_UCL,
      Cohen_Kappa = met$Cohen_Kappa,
      Cohen_Kappa_LCL = met$Cohen_Kappa_LCL,
      Cohen_Kappa_UCL = met$Cohen_Kappa_UCL,
      Observed_Prevalence = met$Observed_Prevalence,
      Predicted_Prevalence = met$Predicted_Prevalence,
      Predicted_Positive_N = met$Predicted_Positive_N,
      stringsAsFactors = FALSE
    )
  }
  
  for (rid in run_ids) {
    scn <- tryCatch(.load_scenario(out_dir, rid), error = function(e) NULL)
    if (is.null(scn)) {
      all_rows[[length(all_rows) + 1]] <- data.frame(
        Run_ID = rid,
        Scenario_File = NA_character_,
        Approach = "ERROR",
        Predictor = NA_character_,
        Threshold = NA_real_,
        AUC = NA_real_, AUC_CI_Low = NA_real_, AUC_CI_High = NA_real_,
        Accuracy = NA_real_, Accuracy_LCL = NA_real_, Accuracy_UCL = NA_real_,
        Sensitivity = NA_real_, Sensitivity_LCL = NA_real_, Sensitivity_UCL = NA_real_,
        Specificity = NA_real_, Specificity_LCL = NA_real_, Specificity_UCL = NA_real_,
        PPV = NA_real_, PPV_LCL = NA_real_, PPV_UCL = NA_real_,
        NPV = NA_real_, NPV_LCL = NA_real_, NPV_UCL = NA_real_,
        Cohen_Kappa = NA_real_, Cohen_Kappa_LCL = NA_real_, Cohen_Kappa_UCL = NA_real_,
        Observed_Prevalence = NA_real_, Predicted_Prevalence = NA_real_, Predicted_Positive_N = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    obj <- scn$obj
    engine <- obj$scoring_engine %||% obj$scoring_engine
    
    if (is.null(engine)) {
      all_rows[[length(all_rows) + 1]] <- data.frame(
        Run_ID = rid,
        Scenario_File = scn$file,
        Approach = "ERROR",
        Predictor = NA_character_,
        Threshold = NA_real_,
        AUC = NA_real_, AUC_CI_Low = NA_real_, AUC_CI_High = NA_real_,
        Accuracy = NA_real_, Accuracy_LCL = NA_real_, Accuracy_UCL = NA_real_,
        Sensitivity = NA_real_, Sensitivity_LCL = NA_real_, Sensitivity_UCL = NA_real_,
        Specificity = NA_real_, Specificity_LCL = NA_real_, Specificity_UCL = NA_real_,
        PPV = NA_real_, PPV_LCL = NA_real_, PPV_UCL = NA_real_,
        NPV = NA_real_, NPV_LCL = NA_real_, NPV_UCL = NA_real_,
        Cohen_Kappa = NA_real_, Cohen_Kappa_LCL = NA_real_, Cohen_Kappa_UCL = NA_real_,
        Observed_Prevalence = NA_real_, Predicted_Prevalence = NA_real_, Predicted_Positive_N = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    scored_ssq <- tryCatch(
      score_new_dataset_for_scenario(new_data, engine, theta_source = "ssq_engine", scenario_obj = obj, phq_items = phq_items),
      error = function(e) NULL
    )
    
    scored_joint <- tryCatch(
      score_new_dataset_for_scenario(new_data, engine, theta_source = "joint", scenario_obj = obj, phq_items = phq_items),
      error = function(e) NULL
    )
    
    df_truth <- scored_joint %||% scored_ssq
    if (is.null(df_truth)) next
    
    y_true <- .compute_reference_outcome(df_truth, outcome_var = outcome_var, phq_items = phq_items, phq_sum_cut = phq_sum_cut)
    
    # AUC is a property of the *continuous predictor* (Theta or SSQ sum) vs the reference outcome.
    # We compute it once here so it can be attached to BOTH the IRT-applied rows and the ROC rows.
    auc_theta_ssq <- list(auc = NA_real_, lcl = NA_real_, ucl = NA_real_)
    auc_theta_joint <- list(auc = NA_real_, lcl = NA_real_, ucl = NA_real_)
    auc_sum_ssq <- list(auc = NA_real_, lcl = NA_real_, ucl = NA_real_)
    
    if (!all(is.na(y_true))) {
      if (!is.null(scored_ssq)) {
        auc_theta_ssq <- .roc_auc_ci(y_true, scored_ssq$Theta_Harmonized, conf.level = conf.level)
      }
      if (!is.null(scored_joint)) {
        auc_theta_joint <- .roc_auc_ci(y_true, scored_joint$Theta_Harmonized, conf.level = conf.level)
      }
      auc_sum_ssq <- .roc_auc_ci(y_true, df_truth$Sum_SSQ, conf.level = conf.level)
    }
    
    if (!is.null(scored_ssq)) {
      met <- .compute_metrics_with_ci(y_true, scored_ssq$Depression_Binary, conf.level = conf.level, n_boot_kappa = n_boot_kappa)
      all_rows[[length(all_rows) + 1]] <- add_row(rid, scn$file, "IRT_applied_ssq_engine", "Theta_Harmonized", NA_real_, auc_theta_ssq$auc, auc_theta_ssq$lcl, auc_theta_ssq$ucl, met)
    }
    
    if (!is.null(scored_joint)) {
      met <- .compute_metrics_with_ci(y_true, scored_joint$Depression_Binary, conf.level = conf.level, n_boot_kappa = n_boot_kappa)
      all_rows[[length(all_rows) + 1]] <- add_row(rid, scn$file, "IRT_applied_joint", "Theta_Harmonized", NA_real_, auc_theta_joint$auc, auc_theta_joint$lcl, auc_theta_joint$ucl, met)
    }
    
    if (all(is.na(y_true))) next
    
    if (!is.null(scored_ssq)) {
      roc_th <- .roc_auc_ci(y_true, scored_ssq$Theta_Harmonized, conf.level = conf.level)
      thr <- .roc_threshold_youden(roc_th$roc)
      pred <- as.integer(scored_ssq$Theta_Harmonized >= thr)
      met <- .compute_metrics_with_ci(y_true, pred, conf.level = conf.level, n_boot_kappa = n_boot_kappa)
      all_rows[[length(all_rows) + 1]] <- add_row(rid, scn$file, "ROC_youden_ssq_engine", "Theta_Harmonized", thr, roc_th$auc, roc_th$lcl, roc_th$ucl, met)
    }
    
    if (!is.null(scored_joint)) {
      roc_th <- .roc_auc_ci(y_true, scored_joint$Theta_Harmonized, conf.level = conf.level)
      thr <- .roc_threshold_youden(roc_th$roc)
      pred <- as.integer(scored_joint$Theta_Harmonized >= thr)
      met <- .compute_metrics_with_ci(y_true, pred, conf.level = conf.level, n_boot_kappa = n_boot_kappa)
      all_rows[[length(all_rows) + 1]] <- add_row(rid, scn$file, "ROC_youden_joint", "Theta_Harmonized", thr, roc_th$auc, roc_th$lcl, roc_th$ucl, met)
    }
    
    roc_ssq <- auc_sum_ssq
    thr_ssq <- .roc_threshold_youden(roc_ssq$roc)
    pred_ssq <- as.integer(df_truth$Sum_SSQ >= thr_ssq)
    met_ssq <- .compute_metrics_with_ci(y_true, pred_ssq, conf.level = conf.level, n_boot_kappa = n_boot_kappa)
    all_rows[[length(all_rows) + 1]] <- add_row(rid, scn$file, "ROC_youden", "SSQ10_Sum", thr_ssq, roc_ssq$auc, roc_ssq$lcl, roc_ssq$ucl, met_ssq)
  }
  
  bind_rows(all_rows) %>% arrange(Run_ID, Predictor, Approach)
}

# -----------------------------------------------------------------------------
# G. Plotting helpers (Test Information + ROC curves)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("ggforce", quietly = TRUE)) install.packages("ggforce")
})

library(ggforce)

.pretty_approach <- function(x) {
  x <- as.character(x)
  x <- gsub("IRT_applied_joint", "IRT applied (Joint θ)", x, fixed = TRUE)
  x <- gsub("IRT_applied_ssq_engine", "IRT applied (SSQ-only θ)", x, fixed = TRUE)
  x <- gsub("ROC_youden_joint", "ROC Youden (Joint θ)", x, fixed = TRUE)
  x <- gsub("ROC_youden_ssq_engine", "ROC Youden (SSQ-only θ)", x, fixed = TRUE)
  x <- gsub("ROC_youden", "ROC Youden (SSQ Sum)", x, fixed = TRUE)
  x
}

.pretty_group_level <- function(level, facet_var) {
  lvl <- as.character(level)
  lvl <- gsub("Age_17_19", "Age 17–19", lvl, fixed = TRUE)
  lvl <- gsub("Age_20_24", "Age 20–24", lvl, fixed = TRUE)
  if (facet_var == "SEX_AGEGRP") {
    lvl <- gsub("__", " × ", lvl, fixed = TRUE)
    return(paste0(lvl))
  }
  if (facet_var == "AGEGRP") return(paste0(lvl))
  if (facet_var == "SEX") return(paste0(lvl))
  lvl
}

.pick_facet_var <- function(scored_df, engine) {
  mode <- engine$cutoff_apply_mode %||% "single"
  if (mode %in% c("none", "single")) return(NULL)
  if (mode %in% c("SEX", "AGEGRP", "SEX_AGEGRP")) return(mode)
  if (mode == "hierarchical") {
    hier <- engine$cutoff_apply_hierarchy %||% c("SEX_AGEGRP", "SEX", "AGEGRP", "none")
    hier <- hier[hier != "none"]
    for (g in hier) {
      if (g %in% names(scored_df)) return(g)
    }
    return(NULL)
  }
  NULL
}

# --- Test Information plot for a single scenario ---
plot_test_information_for_scenario <- function(run_id,
                                               out_dir = here::here("statistical_analysis/output/objects/irt"),
                                               theta_min = -4,
                                               theta_max = 4,
                                               step = 0.05,
                                               phq_items = paste0("PHQ9", sprintf("%02d", 1:9)),
                                               save_path = NULL) {
  scn <- .load_scenario(out_dir, run_id)
  obj <- scn$obj
  if (is.null(obj$models) || is.null(obj$models$joint)) stop("Scenario has no joint model; cannot compute test information.")
  
  mod_joint <- obj$models$joint
  engine <- obj$scoring_engine
  if (is.null(engine) || is.null(engine$items)) stop("Scenario has no scoring engine items; cannot identify SSQ items.")
  
  ssq_items <- engine$items
  n_dims <- mod_joint@Model$nfact
  
  t_seq <- seq(theta_min, theta_max, by = step)
  if (n_dims == 1) {
    Theta_grid <- matrix(t_seq)
  } else {
    Theta_grid <- matrix(0, nrow = length(t_seq), ncol = n_dims)
    Theta_grid[, 1] <- t_seq
  }
  
  phq_idx <- seq_along(phq_items)
  ssq_idx <- (length(phq_items) + 1):(length(phq_items) + length(ssq_items))
  
  info_phq <- tryCatch(mirt::testinfo(mod_joint, Theta_grid, which.items = phq_idx), error = function(e) rep(NA_real_, length(t_seq)))
  info_ssq <- tryCatch(mirt::testinfo(mod_joint, Theta_grid, which.items = ssq_idx), error = function(e) rep(NA_real_, length(t_seq)))
  
  dfp <- data.frame(
    Theta = rep(t_seq, 2),
    Information = c(info_phq, info_ssq),
    Scale = rep(c("PHQ-9", "SSQ-10"), each = length(t_seq))
  )
  
  # One vertical line at the overall (none/Overall) theta cut if available
  thr <- NA_real_
  if (!is.null(engine$cutoff_table)) {
    ct <- engine$cutoff_table
    if (all(c("Grouping", "Group_Level", "Cutoff_Theta") %in% names(ct))) {
      r <- ct[ct$Grouping == "none" & ct$Group_Level == "Overall", , drop = FALSE]
      if (nrow(r) == 1) thr <- suppressWarnings(as.numeric(r$Cutoff_Theta[1]))
    }
  }
  if (is.na(thr)) thr <- engine$threshold %||% engine$cutoff_theta %||% NA_real_
  
  p <- ggplot(dfp, aes(x = Theta, y = Information, color = Scale)) +
    geom_line(linewidth = 1) +
    theme_minimal(base_size = 13) +
    labs(
      #title = paste0("Test Information Functions (Run: ", run_id, ")"),
      x = "Latent Depression (Theta)",
      y = "Information",
      color = "Scale"
    ) +
    theme(legend.position = "right")
  
  if (!is.na(thr)) {
    ssq_df <- dfp[dfp$Scale == "SSQ-10", , drop = FALSE]
    y_end <- ssq_df$Information[which.min(abs(ssq_df$Theta - thr))]
    if (is.finite(y_end)) {
      p <- p + geom_segment(x = thr, xend = thr, y = 0, yend = y_end,
                            linetype = "dashed", linewidth = 0.8, color = "black")
    } else {
      p <- p + geom_vline(xintercept = thr, linetype = "dashed", linewidth = 0.8, color = "black")
    }
  }
  
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = 8.5, height = 5.5, dpi = 300)
  }
  p
}

# --- ROC curve plot for a single scenario (curves colored by approach; points overlay Youden + Applied cutoff)
# For group/hierarchical applied cutoffs, this function facets by the relevant group variable.
plot_roc_for_scenario <- function(run_id,
                                  new_data,
                                  out_dir = here::here("statistical_analysis/output/objects/irt"),
                                  outcome_var = NULL,
                                  phq_items = paste0("PHQ9", sprintf("%02d", 1:9)),
                                  phq_sum_cut = 10,
                                  conf.level = 0.95,
                                  save_path = NULL) {
  scn <- .load_scenario(out_dir, run_id)
  obj <- scn$obj
  engine <- obj$scoring_engine
  if (is.null(engine)) stop("Scenario has no scoring engine")
  
  # Score in BOTH ways (when possible)
  scored_ssq <- tryCatch(
    score_new_dataset_for_scenario(new_data, engine, theta_source = "ssq_engine", scenario_obj = obj, phq_items = phq_items),
    error = function(e) NULL
  )
  scored_joint <- tryCatch(
    score_new_dataset_for_scenario(new_data, engine, theta_source = "joint", scenario_obj = obj, phq_items = phq_items),
    error = function(e) NULL
  )
  
  # Build a single plotting frame that contains BOTH predictors (so curves don't accidentally reuse the same column)
  df_base <- as.data.frame(new_data)
  df_base <- .add_agegrp_if_needed(df_base, age_var = "AGE")
  
  ssq_items <- engine$items
  if (!all(ssq_items %in% names(df_base))) stop("New data missing SSQ items needed for ROC plotting.")
  df_base[, ssq_items] <- lapply(df_base[, ssq_items, drop = FALSE], .coerce_item_0based)
  df_base$Sum_SSQ <- rowSums(df_base[, ssq_items, drop = FALSE], na.rm = TRUE)
  
  # Reference outcome
  y_true <- .compute_reference_outcome(df_base, outcome_var = outcome_var, phq_items = phq_items, phq_sum_cut = phq_sum_cut)
  if (all(is.na(y_true)) || length(unique(y_true[!is.na(y_true)])) < 2) {
    stop("No valid reference outcome (y_true) available to compute ROC.")
  }
  
  # Attach both theta predictors (if available)
  df_base$Theta_SSQ <- if (!is.null(scored_ssq)) scored_ssq$Theta_Harmonized else NA_real_
  df_base$Theta_Joint <- if (!is.null(scored_joint)) scored_joint$Theta_Harmonized else NA_real_
  
  # Applied cutoff for THIS scenario (row-wise, depends on grouping strategy)
  df_base$Cutoff_Theta_Applied <- .get_row_cutoff(df_base, engine, sex_var = "SEX")
  
  # Create interaction group if needed for faceting
  if (!"SEX_AGEGRP" %in% names(df_base) && "SEX" %in% names(df_base) && "AGEGRP" %in% names(df_base)) {
    df_base$SEX_AGEGRP <- paste0(df_base$SEX, "__", df_base$AGEGRP)
  }
  
  # Decide whether to facet (group/hierarchical applied cutoffs)
  facet_var <- .pick_facet_var(df_base, engine)
  
  # helper to build ROC long data per facet
  build_curve_df <- function(df_sub, y_sub, predictor, approach_label, facet_label) {
    ok <- !is.na(y_sub) & !is.na(predictor)
    y <- y_sub[ok]
    x <- predictor[ok]
    if (length(unique(y)) < 2 || length(y) < 20) return(NULL)
    r <- pROC::roc(response = y, predictor = x, quiet = TRUE, direction = "auto")
    data.frame(
      FPR = 1 - r$specificities,
      TPR = r$sensitivities,
      Approach = approach_label,
      Facet = facet_label,
      stringsAsFactors = FALSE
    )
  }
  
  # helper to add points
  add_point <- function(y_sub, predictor, thr, approach_label, facet_label, point_label) {
    ok <- !is.na(y_sub) & !is.na(predictor)
    y <- y_sub[ok]
    x <- predictor[ok]
    if (length(unique(y)) < 2 || length(y) < 20) return(NULL)
    r <- pROC::roc(response = y, predictor = x, quiet = TRUE, direction = "auto")
    cd <- tryCatch(pROC::coords(r, x = thr, input = "threshold", ret = c("specificity", "sensitivity"), transpose = FALSE),
                   error = function(e) NULL)
    if (is.null(cd)) return(NULL)
    data.frame(
      FPR = 1 - as.numeric(cd["specificity"]),
      TPR = as.numeric(cd["sensitivity"]),
      Approach = approach_label,
      Facet = facet_label,
      Point = point_label,
      stringsAsFactors = FALSE
    )
  }
  
  # build per-facet datasets
  if (is.null(facet_var)) {
    splits <- list(Overall = df_base)
  } else {
    lv <- as.character(df_base[[facet_var]])
    lv[is.na(lv)] <- "Missing"
    splits <- split(df_base, lv, drop = TRUE)
    names(splits) <- vapply(names(splits), .pretty_group_level, character(1), facet_var = facet_var)
  }
  
  df_long <- data.frame()
  df_pts <- data.frame()
  
  for (facet_label in names(splits)) {
    df_sub <- splits[[facet_label]]
    y_sub <- .compute_reference_outcome(df_sub, outcome_var = outcome_var, phq_items = phq_items, phq_sum_cut = phq_sum_cut)
    if (all(is.na(y_sub)) || length(unique(y_sub[!is.na(y_sub)])) < 2) next
    
    # --- SSQ-only theta curve (portable predictor) ---
    if (!all(is.na(df_sub$Theta_SSQ))) {
      df_long <- dplyr::bind_rows(df_long,
                                  build_curve_df(df_sub, y_sub, df_sub$Theta_SSQ,
                                                 "IRT applied (SSQ-only θ)", facet_label))
      
      roc_ssq <- .roc_auc_ci(y_sub, df_sub$Theta_SSQ, conf.level = conf.level)$roc
      thr_y <- .roc_threshold_youden(roc_ssq)
      df_pts <- dplyr::bind_rows(df_pts,
                                 add_point(y_sub, df_sub$Theta_SSQ, thr_y,
                                           "IRT applied (SSQ-only θ)", facet_label, "Youden"))
      
      thr_appl <- suppressWarnings(as.numeric(stats::median(unique(df_sub$Cutoff_Theta_Applied), na.rm = TRUE)))
      if (is.finite(thr_appl)) {
        df_pts <- dplyr::bind_rows(df_pts,
                                   add_point(y_sub, df_sub$Theta_SSQ, thr_appl,
                                             "IRT applied (SSQ-only θ)", facet_label, "Applied"))
      }
    }
    
    # --- Joint theta curve (only if joint scoring succeeded) ---
    if (!all(is.na(df_sub$Theta_Joint))) {
      df_long <- dplyr::bind_rows(df_long,
                                  build_curve_df(df_sub, y_sub, df_sub$Theta_Joint,
                                                 "IRT applied (Joint θ)", facet_label))
      
      roc_j <- .roc_auc_ci(y_sub, df_sub$Theta_Joint, conf.level = conf.level)$roc
      thr_y <- .roc_threshold_youden(roc_j)
      df_pts <- dplyr::bind_rows(df_pts,
                                 add_point(y_sub, df_sub$Theta_Joint, thr_y,
                                           "IRT applied (Joint θ)", facet_label, "Youden"))
      
      thr_appl <- suppressWarnings(as.numeric(stats::median(unique(df_sub$Cutoff_Theta_Applied), na.rm = TRUE)))
      if (is.finite(thr_appl)) {
        df_pts <- dplyr::bind_rows(df_pts,
                                   add_point(y_sub, df_sub$Theta_Joint, thr_appl,
                                             "IRT applied (Joint θ)", facet_label, "Applied"))
      }
    }
    
    # --- SSQ sum curve ---
    if ("Sum_SSQ" %in% names(df_sub)) {
      df_long <- dplyr::bind_rows(df_long, build_curve_df(df_sub, y_sub, df_sub$Sum_SSQ, "ROC Youden (SSQ Sum)", facet_label))
      roc_s <- .roc_auc_ci(y_sub, df_sub$Sum_SSQ, conf.level = conf.level)$roc
      thr_y <- .roc_threshold_youden(roc_s)
      df_pts <- dplyr::bind_rows(df_pts, add_point(y_sub, df_sub$Sum_SSQ, thr_y, "ROC Youden (SSQ Sum)", facet_label, "Youden"))
    }
  }
  
  if (nrow(df_long) == 0) stop("No ROC curves could be created (check sample size and class balance in groups).")
  
  p <- ggplot(df_long, aes(x = FPR, y = TPR, color = Approach)) +
    geom_line(linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey40") +
    coord_equal() +
    theme_minimal(base_size = 13) +
    labs(
      #title = paste0("ROC Curves (Run: ", run_id, ")"),
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "Approach",
      shape = "Point"
    ) +
    theme(legend.position = "right")
  
  if (nrow(df_pts) > 0) {
    p <- p + geom_point(data = df_pts, aes(x = FPR, y = TPR, color = Approach, shape = Point),
                        size = 2.7, stroke = 0.8)
  }
  
  # Facet for group/hierarchical
  if (!is.null(facet_var)) {
    p <- p + ggforce::facet_col(facets = vars(Facet))
  }
  
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = 8.5, height = ifelse(is.null(facet_var), 5.5, 8.5), dpi = 300)
  }
  p
}

# --- Convenience: make plots for many run_ids and optionally save them ---
make_scenario_plots <- function(run_ids,
                                new_data,
                                out_dir = here::here("statistical_analysis/output/objects/irt"),
                                plot_dir = file.path(out_dir, "plots"),
                                outcome_var = NULL,
                                phq_sum_cut = 10,
                                conf.level = 0.95,
                                save = TRUE) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  run_ids <- unique(as.character(run_ids))
  
  plots <- list()
  for (rid in run_ids) {
    ti_path <- if (save) file.path(plot_dir, paste0("TestInfo__", rid, ".png")) else NULL
    roc_path <- if (save) file.path(plot_dir, paste0("ROC__", rid, ".png")) else NULL
    
    plots[[rid]] <- list(
      test_info = tryCatch(plot_test_information_for_scenario(rid, out_dir = out_dir, save_path = ti_path), error = function(e) e),
      roc = tryCatch(plot_roc_for_scenario(rid, new_data = new_data, out_dir = out_dir, outcome_var = outcome_var, phq_sum_cut = phq_sum_cut, conf.level = conf.level, save_path = roc_path), error = function(e) e)
    )
  }
  
  invisible(plots)
}


# -----------------------------------------------------------------------------
rank_batch_scenarios <- function(summary_df,
                                 evaluation = "none",
                                 min_specificity = NA_real_,
                                 min_sensitivity = NA_real_,
                                 maximize = c("Sensitivity", "Kappa", "Accuracy"),
                                 top_n = 20) {
  maximize <- match.arg(maximize)
  if (is.null(summary_df) || nrow(summary_df) == 0) return(summary_df)
  
  df <- summary_df
  if ("Evaluation" %in% names(df)) df <- df[df$Evaluation == evaluation, , drop = FALSE]
  
  # Coerce numeric columns (CSV round-trips can make them character)
  num_cols <- intersect(c("Accuracy", "Sensitivity", "Specificity", "Cohen_Kappa", "Cutoff_Theta_Overall"), names(df))
  for (cc in num_cols) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  
  if (!is.na(min_specificity) && "Specificity" %in% names(df)) df <- df[df$Specificity >= min_specificity, , drop = FALSE]
  if (!is.na(min_sensitivity) && "Sensitivity" %in% names(df)) df <- df[df$Sensitivity >= min_sensitivity, , drop = FALSE]
  
  if (nrow(df) == 0) return(df)
  
  # Map maximize target to available column names
  max_col <- maximize
  if (!max_col %in% names(df)) {
    if (identical(max_col, "Kappa") && "Cohen_Kappa" %in% names(df)) max_col <- "Cohen_Kappa"
  }
  
  ord <- order(df[[max_col]], decreasing = TRUE, na.last = TRUE)
  df <- df[ord, , drop = FALSE]
  head(df, top_n)
}


out_dir <- here::here("statistical_analysis/output/objects/irt")
batch_res <- readRDS(file.path(out_dir,"batch_res.rds"))
out_dir <- here::here("statistical_analysis/output/objects/irt/batch")
# -----------------------------------------------------------------------------

apply_hierarchy <- c("SEX_AGEGRP", "SEX", "AGEGRP", "none")
top <- rank_batch_scenarios(
  summary_df = batch_res$summary_overall,
  min_specificity = .95,
  maximize = "Sensitivity",
  top_n = 1000
) |> 
  #filter(grepl("1F",JOINT_Model)) |> 
  select(-c(PHQ_Model,SSQ_Model,Out_File)) |> 
  distinct(JOINT_Model,Apply_Strategy,Cutoff_Method,
           Cutoff_Theta_Overall,Accuracy,Sensitivity,Specificity,
           .keep_all = TRUE) |> 
  distinct(Cutoff_Theta_Overall,Accuracy,Sensitivity,Specificity,
           .keep_all = TRUE) |> 
  filter(Apply_Strategy!="hierarchical")|> 
  arrange(Cutoff_Method) |> 
  tibble::remove_rownames() |> 
  mutate(tag=stringr::str_to_upper(paste0(Cutoff_Method," ",JOINT_Model," (",Apply_Strategy,")")),
         tag=gsub("SEX","Gender",
                  gsub("AGEGRP","Age Group",
                       gsub("SEX_AGEGRP","Gender & Age Group",
                            gsub("F","-Factor Model",
                                 gsub("NONE","Overall",
                                      gsub("JOINT_|MODEL_|SPLIT_|_VALIDATED","",
                                           gsub("_AT_"," at ",
                                                tag)))))))
  )

####------------------
run_ids <- top$Run_ID
# New dataset must contain SSQ items for scoring; include PHQ items or outcome_var for metrics.
psychometric_data_path <- .default_psychometric_data_path()
if (!file.exists(psychometric_data_path) && !exists("dt_psychometric", inherits = TRUE)) {
  stop(
    "Data file not found at: ", psychometric_data_path,
    ". Set COLU_PSYCHOMETRIC_RDATA or place dt_psychometric.RData under data/inputs/."
  )
}

# --- Ensure grouping vars exist ---
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
load_and_clean_data <- function(ssq_idx, phq_idx,
                                data_path = .default_psychometric_data_path()) {
  
  if (!file.exists(data_path)) {
    if (exists("dt_psychometric", inherits = TRUE)) {
      dt_main <- get("dt_psychometric", inherits = TRUE)
    } else {
      stop(
        "Data file not found at: ", data_path,
        ". Set COLU_PSYCHOMETRIC_RDATA or place dt_psychometric.RData under data/inputs/."
      )
    }
  } else {
    loaded_objects <- load(file = data_path)
    if ("dt_psychometric" %in% loaded_objects) {
      dt_main <- dt_psychometric
    } else if ("dt_main" %in% loaded_objects) {
      dt_main <- dt_main
    } else {
      dt_main <- get(loaded_objects[1])
    }
  }
  
  ssq_items <- paste0("SSQ", sprintf("%02d", ssq_idx))
  phq_items <- paste0("PHQ9", sprintf("%02d", phq_idx))
  
  keep_cols <- intersect(c("USUBJID", "SEX", "AGE", ssq_items, phq_items), names(dt_main))
  df_irt <- dt_main[, keep_cols]
  
  df_irt <- .ensure_groups(df_irt)
  
  # Convert factor responses to 0-based integers
  for (v in ssq_items) {
    if (v %in% names(df_irt)) {
      df_irt[[v]] <- if (is.factor(df_irt[[v]])) as.numeric(df_irt[[v]]) - 1 else df_irt[[v]]
    }
  }
  for (v in phq_items) {
    if (v %in% names(df_irt)) {
      df_irt[[v]] <- if (is.factor(df_irt[[v]])) as.numeric(df_irt[[v]]) - 1 else df_irt[[v]]
    }
  }
  
  list(data = df_irt, ssq_items = ssq_items, phq_items = phq_items)
}

prep <- load_and_clean_data(
  ssq_idx = c(1:3, 8:14),
  phq_idx = 1:9,
  data_path = psychometric_data_path
)

dt_metrics_scenario <- score_dataset_across_runids(
  run_ids = run_ids,
  new_data = prep$data,
  out_dir = out_dir,
  outcome_var = NULL,          # or a column like "Depression_PHQ_Raw"
  phq_sum_cut = 10,
  conf.level = 0.95,
  n_boot_kappa = 500
)

dt..plot<-dt_metrics_scenario |> gather(
  measure,values, AUC:Cohen_Kappa_UCL
  ) |> 
  mutate(metric=gsub("_"," ",ifelse(grepl("CI|CL",measure),NA,measure)),
         measure=ifelse(!grepl("CI|CL",measure),"Statistic",measure),
         measure=ifelse(grepl("CI_Low|LCL",measure),"lower",measure),
         measure=ifelse(grepl("CI_High|UCL",measure),"upper",measure),
         metric=factor(metric,levels = unique(metric)),
         Predictor=factor(Predictor,levels = sort(unique(Predictor),decreasing = TRUE))
  ) |> 
  fill(metric,.direction = "down") |> 
  group_by(Approach,Predictor) |> 
  spread(measure,values) |> 
  mutate(Approach=.pretty_approach(Approach)) |> 
  ungroup() |> 
  right_join(
    top |> select(Run_ID,tag)
  )

s=1
p0<-plot_test_information_for_scenario(run_ids[s], out_dir = out_dir)+
  scale_color_uchicago()+
  ggplot2::guides(color = guide_legend(title = NULL, ncol = 1))+
  theme_forest(base_size = 16)+
  ggstatsplot::theme_ggstatsplot()+
  theme(text=element_text(size = 16))

p1<-dt..plot |> 
  filter(Run_ID==run_ids[s]) |> 
  ggplot()+
  
  aes(x=Statistic, 
      y=Approach,
      xmin = lower, 
      xmax = ifelse(round(lower,3)==round(upper,3),Statistic,upper),
      colour = paste0(interaction(Approach,
        sprintf("%.1f",Predicted_Prevalence*100),sep="; "),
        "%"),
      #linetype=Predictor,
      group = Predictor
  ) +
  
  geom_errorbarh(width=0.25)+
  geom_point()+
  # geom_text(aes(x=-.5,label = paste0(sprintf("%.1f",ifelse(!grepl("Cohen",metric),Statistic*100,Statistic))," (",
  #                                     sprintf("%.1f",ifelse(!grepl("Cohen",metric),lower*100,lower)),",",
  #                                     sprintf("%.1f",ifelse(!grepl("Cohen",metric),upper*100,upper)),ifelse(!grepl("Cohen",metric),")%",")")
  #                                     )
  #               ),color="black",hjust=0)+
  ggforce::facet_col(facets = ~metric)+
  scale_color_jama()+
  scale_x_continuous(breaks = seq(0,1,.25))+
  scale_linetype_manual(values = c("dashed","solid"))+
  labs(y="")+
  theme_forest(base_size = 16) +
  theme(axis.text.y = element_blank())+
  ggplot2::guides(color = guide_legend(title = NULL, ncol = 3))
p1

p2 <- plot_roc_for_scenario(run_ids[s], new_data = prep$data, out_dir = out_dir)+
  scale_color_jama()+
  theme_forest(base_size = 16)+
  guides(color = "none", fill = "none",shape=guide_legend(title = NULL, ncol = 1))+
  ggstatsplot::theme_ggstatsplot()+
  theme(text=element_text(size = 16))
