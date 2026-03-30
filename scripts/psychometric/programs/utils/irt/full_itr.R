# Script 03: Latent Harmonization via Item Response Theory (IRT)
# ==============================================================================
# PURPOSE:
# 1. Fit GRMs to SSQ-10 and PHQ-9 (supports exploratory and syntax-based models)
# 2. Check DIF / invariance (SEX, AGEGRP, and SEX_AGEGRP 4-level interaction)
# 3. Estimate latent trait (Theta) for all participants
# 4. Establish cut-offs on Theta (model-based TCC and empirical calibration)
# 5. Generate harmonized binary classification (overall or group-specific / hierarchical)
# 6. Export scoring engine to score NEW datasets (SSQ only)
# 7. Batch runner to compare multiple scenarios transparently
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
  if (!requireNamespace("mirt", quietly = TRUE)) install.packages("mirt")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("psych", quietly = TRUE)) install.packages("psych")
  if (!requireNamespace("gt", quietly = TRUE)) install.packages("gt")
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
})

library(here)
library(mirt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(psych)
library(gt)
library(stringr)

# ==============================================================================
# A. HELPERS
# ==============================================================================

# --- Robust output dir default (user requested) ---
.default_out_dir <- function() {
  here("statistical_analysis/output/objects/irt")
}

.resolve_repo_data_path <- function(env_var, relative_path) {
  override <- Sys.getenv(env_var, unset = "")
  if (nzchar(override)) return(override)
  here(relative_path)
}

.default_psychometric_data_path <- function() {
  .resolve_repo_data_path("COLU_PSYCHOMETRIC_RDATA", "data/inputs/dt_psychometric.RData")
}

.default_dreams_data_path <- function() {
  .resolve_repo_data_path("COLU_DREAMS_RDATA", "data/inputs/dt_dreams_multi_ssq.RData")
}

# --- Auto method chooser (helps QMCEM/MCEM for higher dims) ---
.choose_method <- function(nfact, method = c("auto", "EM", "QMCEM", "MCEM")) {
  method <- match.arg(method)
  if (method != "auto") return(method)
  if (is.na(nfact) || nfact <= 2) return("EM")
  return("QMCEM")
}

# --- Lavaan-like syntax to mirt.model ---
lavaan_to_mirt <- function(syntax) {
  # If numeric (e.g., 1 or 2 factors exploratory), return as is
  if (is.null(syntax) || is.numeric(syntax)) return(syntax)
  
  # Basic cleaning
  s <- syntax
  s <- gsub("=~", "=", s)
  s <- gsub("\\+", ",", s)
  
  # Parse factor names
  lines <- unlist(strsplit(s, "\n"))
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  
  factor_lines <- grep("=", lines, fixed = TRUE)
  factor_names <- sapply(strsplit(lines[factor_lines], "=", fixed = TRUE), function(x) trimws(x[1]))
  factor_names <- unique(factor_names)
  
  # If >1 factor and no explicit COV statement, add correlated factors
  if (length(factor_names) > 1 && !any(grepl("^COV\\s*=", lines))) {
    cov_stmt <- paste0("COV = ", paste(factor_names, collapse = "*"))
    s <- paste(s, cov_stmt, sep = "\n")
  }
  
  mirt.model(s)
}


# --- Ensure grouping vars exist ---

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

# --- Data loading ---
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

# --- Fit table helper ---
create_fit_table <- function(m2_obj, title_str) {
  if (is.null(m2_obj)) return(NULL)
  target_cols <- c("M2", "df", "p", "RMSEA", "RMSEA_5", "RMSEA_95", "SRMSR", "TLI", "CFI")
  valid_cols <- intersect(target_cols, names(m2_obj))
  if (length(valid_cols) == 0) return(NULL)
  
  fit_df <- as.data.frame(m2_obj)[, valid_cols, drop = FALSE]
  gt(fit_df) %>%
    tab_header(title = title_str) %>%
    fmt_number(columns = any_of(c("M2", "RMSEA", "RMSEA_5", "RMSEA_95", "SRMSR", "TLI", "CFI")), decimals = 3) %>%
    fmt_number(columns = any_of("p"), decimals = 4) %>%
    fmt_number(columns = any_of("df"), decimals = 0)
}

# --- Cache id helper (stable signature) ---
if (!exists("make_cache_id", mode = "function")) {
  make_cache_id <- function(x) {
    raw <- serialize(x, NULL)
    if (requireNamespace("digest", quietly = TRUE)) {
      return(digest::digest(raw, algo = "md5", serialize = FALSE))
    }
    r <- as.integer(raw)
    paste0("L", length(r), "_S", sum(r), "_W", sum((seq_along(r) %% 10007L) * r))
  }
}

# ==============================================================================
# B. MODEL FITTING + INVARIANCE
# ==============================================================================

fit_grm_models <- function(df_irt, ssq_items, phq_items,
                           model_phq = 1, model_ssq = 1, model_joint = 1,
                           method_phq = c("auto", "EM", "QMCEM", "MCEM"),
                           method_ssq = c("auto", "EM", "QMCEM", "MCEM"),
                           method_joint = c("auto", "EM", "QMCEM", "MCEM"),
                           technical_phq = list(),
                           technical_ssq = list(),
                           technical_joint = list()) {
  
  mirt_phq <- lavaan_to_mirt(model_phq)
  mirt_ssq <- lavaan_to_mirt(model_ssq)
  mirt_joint <- lavaan_to_mirt(model_joint)
  
  # PHQ
  phq_dat <- df_irt[, phq_items, drop = FALSE]
  phq_dat <- phq_dat[rowSums(!is.na(phq_dat)) > 0, , drop = FALSE]
  nfact_phq <- if (is.numeric(model_phq)) as.integer(model_phq) else NA_integer_
  method_phq <- .choose_method(nfact_phq, method_phq)
  mod_phq <- mirt(phq_dat, mirt_phq, itemtype = "graded", method = method_phq,
                  technical = technical_phq, verbose = FALSE)
  m2_phq <- tryCatch(M2(mod_phq, type = "C2", calcNull = FALSE), error = function(e) NULL)
  
  # SSQ
  ssq_dat <- df_irt[, ssq_items, drop = FALSE]
  ssq_dat <- ssq_dat[rowSums(!is.na(ssq_dat)) > 0, , drop = FALSE]
  nfact_ssq <- if (is.numeric(model_ssq)) as.integer(model_ssq) else NA_integer_
  method_ssq <- .choose_method(nfact_ssq, method_ssq)
  mod_ssq <- mirt(ssq_dat, mirt_ssq, itemtype = "graded", method = method_ssq,
                  technical = technical_ssq, verbose = FALSE)
  m2_ssq <- tryCatch(M2(mod_ssq, type = "C2", calcNull = FALSE), error = function(e) NULL)
  
  # JOINT
  joint_dat <- df_irt[, c(phq_items, ssq_items), drop = FALSE]
  joint_dat <- joint_dat[rowSums(!is.na(joint_dat)) > 0, , drop = FALSE]
  
  nfact_joint <- tryCatch({
    if (is.numeric(model_joint)) as.integer(model_joint) else {
      mm <- mirt_joint
      if (inherits(mm, "mirt.model")) mm@nfact else NA_integer_
    }
  }, error = function(e) NA_integer_)
  
  method_joint <- .choose_method(nfact_joint, method_joint)
  
  mod_joint <- mirt(joint_dat, mirt_joint, itemtype = "graded", method = method_joint,
                    technical = technical_joint, verbose = FALSE)
  
  list(
    models = list(phq = mod_phq, ssq = mod_ssq, joint = mod_joint),
    fit = list(phq_m2 = m2_phq, ssq_m2 = m2_ssq),
    tables = list(
      phq = create_fit_table(m2_phq, paste0("PHQ-9 Model Fit (", method_phq, ")")),
      ssq = create_fit_table(m2_ssq, paste0("SSQ-10 Model Fit (", method_ssq, ")"))
    ),
    methods_used = list(phq = method_phq, ssq = method_ssq, joint = method_joint)
  )
}

run_fitting_workflow <- function(df_irt, ssq_items, phq_items,
                                 phq_model, ssq_model, joint_model,
                                 method_phq, method_ssq, method_joint,
                                 technical_phq, technical_ssq, technical_joint,
                                 output_dir, cache_tag = NULL, overwrite = FALSE) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  sig <- list(
    phq_model = phq_model,
    ssq_model = ssq_model,
    joint_model = joint_model,
    method_phq = method_phq,
    method_ssq = method_ssq,
    method_joint = method_joint,
    technical_phq = technical_phq,
    technical_ssq = technical_ssq,
    technical_joint = technical_joint,
    tag = cache_tag
  )
  cid <- .hash_id(sig, 16)
  fit_file <- file.path(output_dir, paste0("03_IRT_fitted_models__", cid, ".rds"))
  
  if (file.exists(fit_file) && !isTRUE(overwrite)) {
    return(readRDS(fit_file))
  }
  
  res <- fit_grm_models(
    df_irt, ssq_items, phq_items,
    model_phq = phq_model,
    model_ssq = ssq_model,
    model_joint = joint_model,
    method_phq = method_phq,
    method_ssq = method_ssq,
    method_joint = method_joint,
    technical_phq = technical_phq,
    technical_ssq = technical_ssq,
    technical_joint = technical_joint
  )
  
  .saveRDS_atomic(res, fit_file)
  res
}

check_group_invariance <- function(item_data, group_vec, group_name, model_def,
                                   method_mg = c("auto", "EM", "QMCEM", "MCEM"),
                                   technical_mg = list()) {
  
  if (length(unique(na.omit(group_vec))) < 2) return(NULL)
  
  dat_grp <- data.frame(item_data, Group = as.factor(group_vec))
  dat_grp <- na.omit(dat_grp)
  
  mirt_model <- lavaan_to_mirt(model_def)
  nfact <- tryCatch({
    if (is.numeric(model_def)) as.integer(model_def) else {
      if (inherits(mirt_model, "mirt.model")) mirt_model@nfact else NA_integer_
    }
  }, error = function(e) NA_integer_)
  
  method_mg <- .choose_method(nfact, method_mg)
  
  # Try increasingly stable invariance configurations if we hit non-PD sigma issues
  inv_try <- list(
    c("slopes", "intercepts", "free_means", "free_var"),
    c("slopes", "intercepts", "free_means"),
    c("slopes", "intercepts")
  )
  
  last_err <- NULL
  mod_mg <- NULL
  used_inv <- NULL
  
  for (inv in inv_try) {
    tmp <- tryCatch({
      multipleGroup(
        dat_grp[, -ncol(dat_grp), drop = FALSE],
        model = mirt_model,
        group = dat_grp$Group,
        invariance = inv,
        method = method_mg,
        technical = technical_mg,
        verbose = FALSE
      )
    }, error = function(e) {
      last_err <<- e
      NULL
    })
    
    if (!is.null(tmp)) {
      mod_mg <- tmp
      used_inv <- inv
      break
    }
  }
  
  if (is.null(mod_mg)) {
    stop("multipleGroup failed for ", group_name, ": ", if (!is.null(last_err)) last_err$message else "unknown error")
  }
  
  dif_res <- tryCatch(DIF(mod_mg, which.par = c("a1", "d"), scheme = "drop", p.adjust = "fdr"),
                      error = function(e) NULL)
  
  coefs <- tryCatch(coef(mod_mg, simplify = TRUE), error = function(e) NULL)
  
  list(model = mod_mg, dif = dif_res, latent_pars = coefs, method = method_mg, invariance_used = used_inv)
}

summarize_invariance <- function(res_obj, group_label) {
  if (is.null(res_obj) || is.null(res_obj$latent_pars)) return(NULL)
  
  grp_names <- names(res_obj$latent_pars)
  if (is.null(grp_names) || length(grp_names) < 2) return(NULL)
  
  ref_grp <- grp_names[1]
  foc_grp <- grp_names[2]
  
  .get_mean_var1 <- function(group_obj) {
    mean1 <- NA_real_
    var1  <- NA_real_
    
    if (is.null(group_obj) || !is.list(group_obj)) return(list(mean1 = mean1, var1 = var1))
    
    # Typical structure from coef(mod_mg, simplify=TRUE): $means and $cov
    if (!is.null(group_obj$means)) {
      m <- suppressWarnings(as.numeric(group_obj$means))
      if (length(m) >= 1) mean1 <- m[1]
    }
    if (!is.null(group_obj$cov)) {
      cm <- tryCatch(as.matrix(group_obj$cov), error = function(e) NULL)
      if (!is.null(cm) && nrow(cm) >= 1 && ncol(cm) >= 1) var1 <- suppressWarnings(as.numeric(cm[1, 1]))
    }
    
    # Fallback: GroupPars (older representations)
    gp <- group_obj$GroupPars
    
    pick_by_name <- function(obj, key_regex) {
      if (is.null(obj)) return(NA_real_)
      if (is.atomic(obj) && !is.null(names(obj))) {
        hit <- grep(key_regex, names(obj))
        if (length(hit) >= 1) return(suppressWarnings(as.numeric(obj[hit[1]])))
      }
      if (is.list(obj) && !is.null(names(obj))) {
        hit <- grep(key_regex, names(obj))
        if (length(hit) >= 1) return(suppressWarnings(as.numeric(obj[[hit[1]]])))
      }
      if (is.matrix(obj) || is.data.frame(obj)) {
        rn <- rownames(obj); cn <- colnames(obj)
        if (!is.null(rn)) {
          hit <- grep(key_regex, rn)
          if (length(hit) >= 1) return(suppressWarnings(as.numeric(obj[hit[1], 1])))
        }
        if (!is.null(cn)) {
          hit <- grep(key_regex, cn)
          if (length(hit) >= 1) return(suppressWarnings(as.numeric(obj[1, hit[1]])))
        }
      }
      NA_real_
    }
    
    if (is.na(mean1)) {
      mean1 <- pick_by_name(gp, "^MEAN_1$")
      if (is.na(mean1)) mean1 <- pick_by_name(gp, "^MEAN(_|\\b)")
    }
    
    if (is.na(var1)) {
      var1 <- pick_by_name(gp, "^COV_11$")
      if (is.na(var1)) var1 <- pick_by_name(gp, "^COV_([0-9]+)\\1$")
    }
    
    list(mean1 = mean1, var1 = var1)
  }
  
  ref_obj <- res_obj$latent_pars[[ref_grp]]
  foc_obj <- res_obj$latent_pars[[foc_grp]]
  
  mv_ref <- .get_mean_var1(ref_obj)
  mv_foc <- .get_mean_var1(foc_obj)
  
  latent_mean_diff <- mv_foc$mean1 - mv_ref$mean1
  latent_var_ratio <- if (!is.na(mv_ref$var1) && mv_ref$var1 > 0) mv_foc$var1 / mv_ref$var1 else NA_real_
  
  dif_df <- res_obj$dif
  sig_items <- character(0)
  if (!is.null(dif_df)) {
    pcol <- if ("adj_p" %in% names(dif_df)) "adj_p" else if ("p" %in% names(dif_df)) "p" else NULL
    if (!is.null(pcol)) sig_items <- rownames(dif_df)[dif_df[[pcol]] < 0.05]
  }
  
  data.frame(
    Grouping = group_label,
    Method = if (!is.null(res_obj$method)) as.character(res_obj$method) else NA_character_,
    Reference_Group = ref_grp,
    Focal_Group = foc_grp,
    Latent_Mean_Diff = round(latent_mean_diff, 3),
    Latent_Var_Ratio = round(latent_var_ratio, 3),
    DIF_Items = if (length(sig_items) > 0) paste(sig_items, collapse = ", ") else "None",
    DIF_Count = length(sig_items),
    stringsAsFactors = FALSE
  )
}

create_dif_table <- function(res_obj, group_label) {
  if (is.null(res_obj) || is.null(res_obj$dif)) return(NULL)
  dif_df <- res_obj$dif
  dif_df$Item <- rownames(dif_df)
  dif_df$Grouping_Variable <- group_label
  cols_to_keep <- intersect(c("Grouping_Variable", "Item", "X2", "df", "p", "adj_p", "AIC", "BIC"), names(dif_df))
  
  gt(dif_df[, cols_to_keep, drop = FALSE]) %>%
    tab_header(title = paste("Detailed DIF Statistics:", group_label)) %>%
    fmt_number(columns = any_of(c("AIC", "BIC", "X2")), decimals = 2) %>%
    fmt_number(columns = any_of(c("p", "adj_p")), decimals = 4) %>%
    fmt_number(columns = any_of("df"), decimals = 0)
}

run_invariance_workflow <- function(df_irt, ssq_items, phq_items, joint_model_def,
                                    method_mg, technical_mg,
                                    output_dir, cache_tag = NULL, overwrite = FALSE) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  sig <- list(
    joint_model_def = joint_model_def,
    method_mg = method_mg,
    technical_mg = technical_mg,
    tag = cache_tag
  )
  cid <- .hash_id(sig, 16)
  inv_file <- file.path(output_dir, paste0("03_IRT_invariance_models__", cid, ".rds"))
  
  if (file.exists(inv_file) && !isTRUE(overwrite)) {
    return(readRDS(inv_file))
  }
  
  joint_data <- df_irt[, c(phq_items, ssq_items), drop = FALSE]
  valid_joint <- rowSums(!is.na(joint_data)) > 0
  clean_dat <- joint_data[valid_joint, , drop = FALSE]
  
  df_tmp <- df_irt[valid_joint, , drop = FALSE]
  df_tmp <- .ensure_groups(df_tmp)
  
  res_sex <- check_group_invariance(clean_dat, df_tmp$SEX, "SEX", joint_model_def, method_mg = method_mg, technical_mg = technical_mg)
  res_age <- check_group_invariance(clean_dat, df_tmp$AGEGRP, "AGEGRP", joint_model_def, method_mg = method_mg, technical_mg = technical_mg)
  res_int <- check_group_invariance(clean_dat, df_tmp$SEX_AGEGRP, "SEX_AGEGRP", joint_model_def, method_mg = method_mg, technical_mg = technical_mg)
  
  summ <- rbind(
    summarize_invariance(res_sex, "Biological Sex"),
    summarize_invariance(res_age, "Age Group"),
    summarize_invariance(res_int, "SEX x AGEGRP (4-level)")
  )
  
  tbl_summary <- if (!is.null(summ) && nrow(summ) > 0) {
    gt(summ) %>%
      tab_header(title = "Measurement Invariance (DIF) & Latent Impact") %>%
      fmt_number(columns = any_of(c("Latent_Mean_Diff", "Latent_Var_Ratio")), decimals = 3)
  } else {
    NULL
  }
  
  out <- list(
    sex = res_sex,
    age = res_age,
    sex_age = res_int,
    summary_df = summ,
    tables = list(summary = tbl_summary,
                  sex_dif = create_dif_table(res_sex, "Biological Sex"),
                  age_dif = create_dif_table(res_age, "Age Group"),
                  sex_age_dif = create_dif_table(res_int, "SEX x AGEGRP")),
    model_def = joint_model_def
  )
  
  .saveRDS_atomic(out, inv_file)
  out
}

# ==============================================================================
# C. LINKING + CUTOFF CALIBRATION
# ==============================================================================

.calc_metrics <- function(y_true, y_pred) {
  y_true <- as.integer(y_true)
  y_pred <- as.integer(y_pred)
  ok <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]
  if (length(y_true) == 0) return(NULL)
  
  tbl <- table(
    Pred = factor(y_pred, levels = c(0, 1)),
    True = factor(y_true, levels = c(0, 1))
  )
  
  TN <- tbl[1, 1]
  FN <- tbl[1, 2]
  FP <- tbl[2, 1]
  TP <- tbl[2, 2]
  total <- sum(tbl)
  
  acc <- (TP + TN) / total
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  ppv  <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  npv  <- if ((TN + FN) > 0) TN / (TN + FN) else 0
  kap <- tryCatch(psych::cohen.kappa(tbl)$kappa, error = function(e) NA_real_)
  
  list(
    Accuracy = acc,
    Sensitivity = sens,
    Specificity = spec,
    PPV = ppv,
    NPV = npv,
    Kappa = kap,
    TP = as.integer(TP),
    FP = as.integer(FP),
    TN = as.integer(TN),
    FN = as.integer(FN),
    Pred_Pos = as.integer(TP + FP),
    Pred_Neg = as.integer(TN + FN),
    True_Pos = as.integer(TP + FN),
    True_Neg = as.integer(TN + FP),
    N = as.integer(total)
  )
}

.select_cutoff_empirical <- function(theta, y_true,
                                     method = c("youden", "sens_at_spec", "spec_at_sens"),
                                     target_specificity = 0.95,
                                     target_sensitivity = 0.80) {
  
  method <- match.arg(method)
  
  ok <- !is.na(theta) & !is.na(y_true)
  theta <- theta[ok]
  y_true <- y_true[ok]
  
  if (length(theta) < 20) return(list(cutoff = NA_real_, summary_row = NULL))
  
  qs <- unique(as.numeric(stats::quantile(theta, probs = seq(0.01, 0.99, by = 0.01), na.rm = TRUE)))
  qs <- sort(unique(qs))
  
  res <- lapply(qs, function(cut) {
    y_pred <- ifelse(theta >= cut, 1L, 0L)
    m <- .calc_metrics(y_true, y_pred)
    if (is.null(m)) return(NULL)
    c(Cutoff_Theta = cut,
      Accuracy = m$Accuracy,
      Sensitivity = m$Sensitivity,
      Specificity = m$Specificity,
      Kappa = m$Kappa)
  })
  
  df <- as.data.frame(do.call(rbind, res))
  if (is.null(df) || nrow(df) == 0) return(list(cutoff = NA_real_, summary_row = NULL))
  
  if (method == "youden") {
    df$Youden <- df$Sensitivity + df$Specificity - 1
    best <- df[which.max(df$Youden), , drop = FALSE]
  } else if (method == "sens_at_spec") {
    keep <- df[df$Specificity >= target_specificity, , drop = FALSE]
    if (nrow(keep) == 0) {
      df$SpecDist <- abs(df$Specificity - target_specificity)
      best <- df[order(df$SpecDist, -df$Sensitivity, -df$Kappa), , drop = FALSE][1, , drop = FALSE]
    } else {
      best <- keep[order(-keep$Sensitivity, -keep$Kappa), , drop = FALSE][1, , drop = FALSE]
    }
  } else {
    keep <- df[df$Sensitivity >= target_sensitivity, , drop = FALSE]
    if (nrow(keep) == 0) {
      df$SensDist <- abs(df$Sensitivity - target_sensitivity)
      best <- df[order(df$SensDist, -df$Specificity, -df$Kappa), , drop = FALSE][1, , drop = FALSE]
    } else {
      best <- keep[order(-keep$Specificity, -keep$Kappa), , drop = FALSE][1, , drop = FALSE]
    }
  }
  
  list(
    cutoff = as.numeric(best$Cutoff_Theta),
    summary_row = best
  )
}

.compute_cutoff_model_tcc <- function(mod_joint, phq_items, phq_sum_cut = 10,
                                      theta_grid = seq(-4, 4, by = 0.01)) {
  
  phq_idx <- seq_along(phq_items)
  n_dims <- mod_joint@Model$nfact
  
  if (n_dims == 1) {
    Theta <- matrix(theta_grid)
  } else {
    Theta <- matrix(0, nrow = length(theta_grid), ncol = n_dims)
    Theta[, 1] <- theta_grid
  }
  
  tcc <- expected.test(mod_joint, Theta = Theta, which.items = phq_idx)
  idx <- which.min(abs(tcc - phq_sum_cut))
  theta_grid[idx]
}

.build_phq_theta_map <- function(mod_joint, phq_items, phq_sum_cut = 10,
                                 theta_grid = seq(-4, 4, by = 0.01)) {
  phq_idx <- seq_along(phq_items)
  n_dims <- mod_joint@Model$nfact
  
  if (n_dims == 1) {
    Theta <- matrix(theta_grid)
  } else {
    Theta <- matrix(0, nrow = length(theta_grid), ncol = n_dims)
    Theta[, 1] <- theta_grid
  }
  
  exp_sum <- expected.test(mod_joint, Theta = Theta, which.items = phq_idx)
  
  # Probability(PHQ sum >= phq_sum_cut) via scoreDist, with safe fallback
  prob_ge <- tryCatch({
    sd <- scoreDist(mod_joint, Theta = Theta, which.items = phq_idx)
    # sd$score: possible summed scores; sd$prob: matrix rows=Theta, cols=scores
    sc <- sd$score
    pr <- sd$prob
    idx_ge <- which(sc >= phq_sum_cut)
    if (length(idx_ge) == 0) rep(0, length(theta_grid)) else rowSums(pr[, idx_ge, drop = FALSE])
  }, error = function(e) {
    # Fallback: logistic approximation based on expected sum (crude, but stable)
    p <- plogis(exp_sum - phq_sum_cut)
    p
  })
  
  data.frame(
    Theta = theta_grid,
    PHQ_Expected_Sum = as.numeric(exp_sum),
    PHQ_Prob_GE10 = as.numeric(prob_ge),
    stringsAsFactors = FALSE
  )
}

.evaluate_cutoffs_by_group <- function(df, group_var, cutoff_method,
                                       target_specificity, target_sensitivity,
                                       cutoff_model_tcc,
                                       label_col = "Depression_PHQ_Raw") {
  
  if (identical(group_var, "none")) {
    levels <- "Overall"
    groups <- list(Overall = df)
  } else {
    if (!group_var %in% names(df)) return(NULL)
    groups <- split(df, df[[group_var]], drop = TRUE)
    levels <- names(groups)
  }
  
  out <- lapply(levels, function(lv) {
    d <- groups[[lv]]
    d <- d[!is.na(d$Theta_Harmonized) & !is.na(d[[label_col]]), , drop = FALSE]
    if (nrow(d) < 30) return(NULL)
    
    if (cutoff_method == "model_tcc") {
      cut <- cutoff_model_tcc
      y_pred <- ifelse(d$Theta_Harmonized >= cut, 1L, 0L)
      m <- .calc_metrics(d[[label_col]], y_pred)
      if (is.null(m)) return(NULL)
      
      data.frame(
        Grouping = group_var,
        Group_Level = ifelse(group_var == "none", "Overall", as.character(lv)),
        Cutoff_Method = cutoff_method,
        Cutoff_Theta = cut,
        Accuracy = m$Accuracy,
        Sensitivity = m$Sensitivity,
        Specificity = m$Specificity,
        PPV = m$PPV,
        NPV = m$NPV,
        Kappa = m$Kappa,
        N = m$N,
        stringsAsFactors = FALSE
      )
      
    } else {
      sel <- .select_cutoff_empirical(
        theta = d$Theta_Harmonized,
        y_true = d[[label_col]],
        method = cutoff_method,
        target_specificity = target_specificity,
        target_sensitivity = target_sensitivity
      )
      cut <- sel$cutoff
      if (is.na(cut)) return(NULL)
      
      y_pred <- ifelse(d$Theta_Harmonized >= cut, 1L, 0L)
      m <- .calc_metrics(d[[label_col]], y_pred)
      if (is.null(m)) return(NULL)
      
      data.frame(
        Grouping = group_var,
        Group_Level = ifelse(group_var == "none", "Overall", as.character(lv)),
        Cutoff_Method = cutoff_method,
        Cutoff_Theta = cut,
        Accuracy = m$Accuracy,
        Sensitivity = m$Sensitivity,
        Specificity = m$Specificity,
        PPV = m$PPV,
        NPV = m$NPV,
        Kappa = m$Kappa,
        N = m$N,
        stringsAsFactors = FALSE
      )
    }
  })
  
  do.call(rbind, out)
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

perform_linking_and_scoring <- function(mod_joint, df_irt, ssq_items, phq_items,
                                        cutoff_method = c("model_tcc", "youden", "sens_at_spec", "spec_at_sens"),
                                        target_specificity = 0.95,
                                        target_sensitivity = 0.80,
                                        cutoff_by_group_eval = c("none", "SEX", "AGEGRP", "SEX_AGEGRP"),
                                        cutoff_by_group_apply = c("none", "SEX", "AGEGRP", "SEX_AGEGRP", "hierarchical"),
                                        cutoff_apply_hierarchy = c("SEX_AGEGRP", "SEX", "AGEGRP", "none")) {
  
  cutoff_method <- match.arg(cutoff_method)
  cutoff_by_group_apply <- match.arg(cutoff_by_group_apply)
  
  df_irt <- .ensure_groups(df_irt)
  
  joint_data <- df_irt[, c(phq_items, ssq_items), drop = FALSE]
  valid_joint <- rowSums(!is.na(joint_data)) > 0
  
  # Use quasi-Monte Carlo integration for factor scores when model is high-dimensional
  nfact_joint <- tryCatch(mod_joint@Model$nfact, error = function(e) NA_integer_)
  use_qmc <- !is.na(nfact_joint) && nfact_joint >= 3
  theta_joint <- fscores(mod_joint, method = "EAP", full.scores = TRUE, QMC = use_qmc)
  df_irt$Theta_Harmonized <- NA_real_
  df_irt$Theta_Harmonized[valid_joint] <- theta_joint[, 1]
  
  df_irt$Sum_PHQ <- rowSums(df_irt[, phq_items, drop = FALSE], na.rm = TRUE)
  df_irt$Sum_SSQ <- rowSums(df_irt[, ssq_items, drop = FALSE], na.rm = TRUE)
  df_irt$Depression_PHQ_Raw <- ifelse(df_irt$Sum_PHQ >= 10, 1L, 0L)
  
  # Build a 1D mapping from Theta -> PHQ expected sum and P(PHQ>=10)
  phq_map <- .build_phq_theta_map(mod_joint, phq_items, phq_sum_cut = 10)
  
  # Attach mapped PHQ quantities to each participant (using Theta_Harmonized)
  df_irt$PHQ_Expected_Sum <- NA_real_
  df_irt$PHQ_Prob_GE10 <- NA_real_
  if (any(!is.na(df_irt$Theta_Harmonized))) {
    df_irt$PHQ_Expected_Sum <- approx(phq_map$Theta, phq_map$PHQ_Expected_Sum, xout = df_irt$Theta_Harmonized, rule = 2)$y
    df_irt$PHQ_Prob_GE10 <- approx(phq_map$Theta, phq_map$PHQ_Prob_GE10, xout = df_irt$Theta_Harmonized, rule = 2)$y
  }
  
  cutoff_model_tcc <- .compute_cutoff_model_tcc(mod_joint, phq_items, phq_sum_cut = 10)
  
  cutoff_by_group_eval <- unique(as.character(cutoff_by_group_eval))
  cutoff_by_group_eval <- cutoff_by_group_eval[cutoff_by_group_eval %in% c("none", "SEX", "AGEGRP", "SEX_AGEGRP")]
  if (length(cutoff_by_group_eval) == 0) cutoff_by_group_eval <- "none"
  
  eval_res <- lapply(cutoff_by_group_eval, function(g) {
    .evaluate_cutoffs_by_group(
      df = df_irt,
      group_var = g,
      cutoff_method = cutoff_method,
      target_specificity = target_specificity,
      target_sensitivity = target_sensitivity,
      cutoff_model_tcc = cutoff_model_tcc,
      label_col = "Depression_PHQ_Raw"
    )
  })
  
  cutoff_table <- do.call(rbind, eval_res)
  
  if (is.null(cutoff_table) || nrow(cutoff_table) == 0) {
    cutoff_table <- data.frame(
      Grouping = "none",
      Group_Level = "Overall",
      Cutoff_Method = cutoff_method,
      Cutoff_Theta = cutoff_model_tcc,
      Accuracy = NA_real_,
      Sensitivity = NA_real_,
      Specificity = NA_real_,
      PPV = NA_real_,
      NPV = NA_real_,
      Kappa = NA_real_,
      N = NA_integer_,
      stringsAsFactors = FALSE
    )
  }
  
  df_scored <- .apply_cutoffs(
    df_irt,
    cutoff_table = cutoff_table,
    apply_mode = cutoff_by_group_apply,
    hierarchy = cutoff_apply_hierarchy
  )
  
  m_app <- .calc_metrics(df_scored$Depression_PHQ_Raw, df_scored$Depression_Harmonized)
  
  metrics_df <- data.frame(
    Metric = c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "Cohen_Kappa"),
    Value = c(m_app$Accuracy, m_app$Sensitivity, m_app$Specificity, m_app$PPV, m_app$NPV, m_app$Kappa)
  )
  
  tbl_class <- gt(metrics_df) %>%
    tab_header(title = "Applied Cutoff: Harmonized (Theta) vs PHQ-9 Raw") %>%
    fmt_number(columns = "Value", decimals = 3) %>%
    fmt_percent(columns = "Value", rows = 1:5, decimals = 1)
  
  # Confusion matrix (Pred rows x True cols)
  conf_mat <- matrix(
    c(m_app$TN, m_app$FN,
      m_app$FP, m_app$TP),
    nrow = 2, byrow = TRUE,
    dimnames = list(Pred = c("0", "1"), True = c("0", "1"))
  )
  conf_df <- as.data.frame(conf_mat)
  conf_df$Pred <- rownames(conf_mat)
  conf_df <- conf_df[, c("Pred", "0", "1"), drop = FALSE]
  
  tbl_conf <- gt(conf_df) %>%
    tab_header(title = "Confusion Matrix (Applied Cutoff)") %>%
    cols_label(`0` = "True 0", `1` = "True 1")
  
  consolidated_df <- cutoff_table %>%
    mutate(
      Accuracy = as.numeric(Accuracy),
      Sensitivity = as.numeric(Sensitivity),
      Specificity = as.numeric(Specificity),
      PPV = as.numeric(PPV),
      NPV = as.numeric(NPV),
      Kappa = as.numeric(Kappa),
      Cutoff_Theta = as.numeric(Cutoff_Theta)
    )
  
  consolidated_gt <- gt(consolidated_df) %>%
    tab_header(title = "Cutoff Evaluation (All Requested Groupings)") %>%
    fmt_number(columns = any_of(c("Cutoff_Theta", "Kappa")), decimals = 3) %>%
    fmt_percent(columns = any_of(c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV")), decimals = 1)
  
  list(
    scores = df_scored,
    cutoff_model_tcc = cutoff_model_tcc,
    cutoff_table = cutoff_table,
    consolidated_df = consolidated_df,
    phq_theta_map = phq_map,
    tables = list(applied_metrics = tbl_class, confusion = tbl_conf, cutoff_eval = consolidated_gt),
    metrics = metrics_df
  )
}

# ==============================================================================
# D. SCORING ENGINE EXPORT (supports overall or group cutoffs)
# ==============================================================================

extract_scoring_engine <- function(mod_joint, ssq_items,
                                   cutoff_table,
                                   cutoff_by_group_apply,
                                   cutoff_apply_hierarchy,
                                   phq_items = paste0("PHQ9", sprintf("%02d", 1:9)),
                                   phq_sum_cut = 10) {
  
  pars <- mod2values(mod_joint)
  pars_ssq <- pars[pars$item %in% ssq_items, , drop = FALSE]
  nfact <- mod_joint@Model$nfact
  
  # Store mapping needed to produce PHQ expected/prob outputs for SSQ-only cohorts
  phq_map <- .build_phq_theta_map(mod_joint, phq_items = phq_items, phq_sum_cut = phq_sum_cut)
  
  list(
    parameters = pars_ssq,
    items = ssq_items,
    model_type = "graded",
    n_factors = nfact,
    cutoff_table = cutoff_table,
    cutoff_apply_mode = cutoff_by_group_apply,
    cutoff_apply_hierarchy = cutoff_apply_hierarchy,
    phq_theta_map = phq_map,
    phq_sum_cut = phq_sum_cut
  )
}

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

# ==============================================================================
# E. MAIN PIPELINE (ONE SCRIPT)
# ==============================================================================

run_irt_harmonization_pipeline <- function(
    phq_model = 1,
    ssq_model = 1,
    joint_model = 1,
    # Estimation
    method_phq = c("auto", "EM", "QMCEM", "MCEM"),
    method_ssq = c("auto", "EM", "QMCEM", "MCEM"),
    method_joint = c("auto", "EM", "QMCEM", "MCEM"),
    technical_phq = list(),
    technical_ssq = list(),
    technical_joint = list(),
    method_mg = c("auto", "EM", "QMCEM", "MCEM"),
    technical_mg = list(),
    # Cutoff calibration
    cutoff_method = c("model_tcc", "youden", "sens_at_spec", "spec_at_sens"),
    target_specificity = 0.95,
    target_sensitivity = 0.80,
    cutoff_by_group_eval = c("none", "SEX", "AGEGRP", "SEX_AGEGRP"),
    cutoff_by_group_apply = c("none", "SEX", "AGEGRP", "SEX_AGEGRP", "hierarchical"),
    cutoff_apply_hierarchy = c("SEX_AGEGRP", "SEX", "AGEGRP", "none"),
    # Sensitivity external
    run_sensitivity = FALSE,
    dreams_path = .default_dreams_data_path(),
    # Output
    output_dir = .default_out_dir(),
    save_results = TRUE,
    results_tag = NULL,
    export_engine = TRUE,
    engine_tag = NULL,
    export_dir = output_dir,
    save_sensitivity_scores = TRUE,
    overwrite_cache = FALSE
) {
  
  method_phq <- match.arg(method_phq)
  method_ssq <- match.arg(method_ssq)
  method_joint <- match.arg(method_joint)
  method_mg <- match.arg(method_mg)
  cutoff_method <- match.arg(cutoff_method)
  cutoff_by_group_apply <- match.arg(cutoff_by_group_apply)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("=== STARTING IRT HARMONIZATION PIPELINE ===")
  
  prep <- load_and_clean_data(ssq_idx = c(1:3, 8:14), phq_idx = 1:9)
  
  fit_res <- run_fitting_workflow(
    prep$data, prep$ssq_items, prep$phq_items,
    phq_model = phq_model,
    ssq_model = ssq_model,
    joint_model = joint_model,
    method_phq = method_phq,
    method_ssq = method_ssq,
    method_joint = method_joint,
    technical_phq = technical_phq,
    technical_ssq = technical_ssq,
    technical_joint = technical_joint,
    output_dir = output_dir,
    cache_tag = "fit",
    overwrite = overwrite_cache
  )
  
  inv_res <- run_invariance_workflow(
    prep$data, prep$ssq_items, prep$phq_items,
    joint_model_def = joint_model,
    method_mg = method_mg,
    technical_mg = technical_mg,
    output_dir = output_dir,
    cache_tag = "inv",
    overwrite = overwrite_cache
  )
  
  link_res <- perform_linking_and_scoring(
    mod_joint = fit_res$models$joint,
    df_irt = prep$data,
    ssq_items = prep$ssq_items,
    phq_items = prep$phq_items,
    cutoff_method = cutoff_method,
    target_specificity = target_specificity,
    target_sensitivity = target_sensitivity,
    cutoff_by_group_eval = cutoff_by_group_eval,
    cutoff_by_group_apply = cutoff_by_group_apply,
    cutoff_apply_hierarchy = cutoff_apply_hierarchy
  )
  
  scoring_engine <- extract_scoring_engine(
    mod_joint = fit_res$models$joint,
    ssq_items = prep$ssq_items,
    cutoff_table = link_res$cutoff_table,
    cutoff_by_group_apply = cutoff_by_group_apply,
    cutoff_apply_hierarchy = cutoff_apply_hierarchy,
    phq_items = prep$phq_items,
    phq_sum_cut = 10
  )
  
  if (isTRUE(export_engine)) {
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    
    eng_name <- if (!is.null(engine_tag) && nzchar(engine_tag)) {
      paste0("SSQ10_Harmonized_Scoring_Engine__", engine_tag, ".rds")
    } else {
      "SSQ10_Harmonized_Scoring_Engine.rds"
    }
    
    joint_name <- if (!is.null(engine_tag) && nzchar(engine_tag)) {
      paste0("Joint_IRT_Model_For_Scoring__", engine_tag, ".rds")
    } else {
      "Joint_IRT_Model_For_Scoring.rds"
    }
    
    saveRDS(scoring_engine, file.path(export_dir, eng_name))
    saveRDS(fit_res$models$joint, file.path(export_dir, joint_name))
  }
  
  sensitivity_scores <- NULL
  if (isTRUE(run_sensitivity) && file.exists(dreams_path)) {
    load(dreams_path)
    if (exists("dt_dreams_multi_ssq")) {
      sensitivity_scores <- tryCatch({
        s <- score_new_cohort(dt_dreams_multi_ssq, scoring_engine)
        if ("USUBJID" %in% names(dt_dreams_multi_ssq)) {
          cbind(dt_dreams_multi_ssq[, "USUBJID", drop = FALSE], s)
        } else {
          s
        }
      }, error = function(e) {
        message("Error scoring Dreams cohort: ", e$message)
        NULL
      })
      
      if (!is.null(sensitivity_scores) && isTRUE(save_sensitivity_scores)) {
        dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
        sens_name <- if (!is.null(results_tag) && nzchar(results_tag)) {
          paste0("Scored_Dreams_Cohort__", results_tag, ".rds")
        } else {
          "Scored_Dreams_Cohort.rds"
        }
        saveRDS(sensitivity_scores, file.path(export_dir, sens_name))
      }
    }
  }
  
  final_output <- list(
    input_structure = list(
      phq_model = phq_model,
      ssq_model = ssq_model,
      joint_model = joint_model,
      methods_used = fit_res$methods_used,
      cutoff_method = cutoff_method,
      target_specificity = target_specificity,
      target_sensitivity = target_sensitivity,
      cutoff_by_group_eval = cutoff_by_group_eval,
      cutoff_by_group_apply = cutoff_by_group_apply,
      cutoff_apply_hierarchy = cutoff_apply_hierarchy
    ),
    models = fit_res$models,
    fit_tables = fit_res$tables,
    invariance = inv_res,
    link_results = link_res,
    scores = link_res$scores,
    metrics_applied = link_res$metrics,
    scoring_engine = scoring_engine,
    sensitivity_scores = sensitivity_scores
  )
  
  if (isTRUE(save_results)) {
    save_name <- if (!is.null(results_tag) && nzchar(results_tag)) {
      paste0("03_IRT_Harmonization_Results__", results_tag, ".rds")
    } else {
      "03_IRT_Harmonization_Results.rds"
    }
    saveRDS(final_output, file.path(output_dir, save_name))
  }
  
  message("Pipeline complete.")
  final_output
}

# ------------------------------------------------------------------------------
# VALIDATED 2-FACTOR STRUCTURES (for this dataset)
# ------------------------------------------------------------------------------
# These MUST exist *before* batch_run_irt() is defined, because they are used as
# default arguments in that function.

if (!exists("phq_structure_validated", inherits = TRUE)) {
  phq_structure_validated <- "F1 =~ PHQ906 + PHQ907 + PHQ908 + PHQ909\nF2 =~ PHQ901 + PHQ902 + PHQ903 + PHQ904 + PHQ905"
}


if (!exists("ssq_structure_validated", inherits = TRUE)) {
  ssq_structure_validated <- "F1 =~ SSQ09 + SSQ10 + SSQ11 + SSQ12 + SSQ13 + SSQ14 + SSQ08\nF2 =~ SSQ01 + SSQ02 + SSQ03"
}


if (!exists("make_joint_validated_4f_structure", mode = "function")) {
  make_joint_validated_4f_structure <- function() {
    paste0(
      "PHQ_F1 =~ PHQ906 + PHQ907 + PHQ908 + PHQ909\n",
      "PHQ_F2 =~ PHQ901 + PHQ902 + PHQ903 + PHQ904 + PHQ905\n",
      "SSQ_F1 =~ SSQ09 + SSQ10 + SSQ11 + SSQ12 + SSQ13 + SSQ14 + SSQ08\n",
      "SSQ_F2 =~ SSQ01 + SSQ02 + SSQ03"
    )
  }
}


# ------------------------------------------------------------------------------
# CACHE ID HELPER (used by batch_run_irt)
# ------------------------------------------------------------------------------
# Generates a stable hash for any R object without requiring extra packages.
# (Uses md5 of a serialized representation)

if (!exists("make_cache_id", mode = "function")) {
  make_cache_id <- function(x) {
    # Stable hash for arbitrary R objects.
    # Prefer digest if available; otherwise use a deterministic fallback.
    raw <- serialize(x, NULL)
    
    if (requireNamespace("digest", quietly = TRUE)) {
      return(digest::digest(raw, algo = "md5", serialize = FALSE))
    }
    
    # Fallback (not cryptographic, but stable): length + weighted sums
    r <- as.integer(raw)
    s1 <- sum(r)
    s2 <- sum((seq_along(r) %% 10007L) * r)
    paste0("L", length(r), "_S", s1, "_W", s2)
  }
}

# ------------------------------------------------------------------------------
# HASH + ATOMIC SAVE (robust for parallel batch runs)
# ------------------------------------------------------------------------------
# We place these helpers here (near batch) to ensure they exist in worker sessions.

.saveRDS_atomic <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(file, ".", Sys.getpid(), ".tmp")
  saveRDS(object, tmp)
  file.rename(tmp, file)
}

make_cache_id <- function(x) {
  # Stable, non-empty hash for arbitrary R objects
  if (!requireNamespace("digest", quietly = TRUE)) install.packages("digest")
  digest::digest(x, algo = "xxhash64", serialize = TRUE)
}

.hash_id <- function(x, n = 16L) {
  substring(make_cache_id(x), 1L, n)
}

# ------------------------------------------------------------------------------
# batch_run_irt(): GRID SEARCH
# ------------------------------------------------------------------------------
# Runs a grid of scenarios and produces:
#   1) One .rds per scenario (full transparency; contains models, cutoffs, tables)
#   2) A single SUMMARY file (overall rows only) for quick comparison
#
# IMPORTANT:
# - During batch, we call run_irt_harmonization_pipeline(save_results=FALSE, export_engine=FALSE)
#   so we never overwrite canonical scoring engine artifacts.
# - Model fitting + invariance checks are still executed *per scenario*, but they are cached
#   via the cache_id logic, so identical scenarios re-use results.

batch_run_irt <- function(
    # Model candidates (IDs + defs)
  phq_models = list(
    `PHQ_1F` = 1,
    `PHQ_2F_validated` = phq_structure_validated
  ),
  ssq_models = list(
    `SSQ_1F` = 1,
    `SSQ_2F_validated` = ssq_structure_validated
  ),
  joint_models = list(
    `JOINT_1F` = 1,
    `JOINT_SPLIT_2F` = paste0(
      "F1 = PHQ901, PHQ902, PHQ903, PHQ904, PHQ905, PHQ906, PHQ907, PHQ908, PHQ909\n",
      "F2 = SSQ01, SSQ02, SSQ03, SSQ08, SSQ09, SSQ10, SSQ11, SSQ12, SSQ13, SSQ14\n",
      "COV = F1*F2"
    ),
    `JOINT_4F_validated` = make_joint_validated_4f_structure()
  ),
  cutoff_methods = c("model_tcc", "youden", "sens_at_spec", "spec_at_sens"),
  target_specificities = c(0.90, 0.95, 0.98),
  target_sensitivities = c(0.70, 0.80, 0.90),
  # Apply/eval
  apply_modes = c("single", "hierarchical"),
  apply_groups = c("none", "SEX", "AGEGRP", "SEX_AGEGRP"),
  eval_groupings = c("none", "SEX", "AGEGRP", "SEX_AGEGRP"),
  apply_hierarchy = c("SEX_AGEGRP", "SEX", "AGEGRP", "none"),
  # Estimation methods
  method_phq = "auto",
  method_ssq = "auto",
  method_joint = "auto",
  technical_phq = list(),
  technical_ssq = list(),
  technical_joint = list(),
  method_mg = "auto",
  technical_mg = list(),
  # Runtime
  run_sensitivity = FALSE,
  overwrite = FALSE,
  workers = 1L,
  output_dir = {
    if (requireNamespace("here", quietly = TRUE)) {
      here::here("statistical_analysis/output/objects/irt")
    } else {
      file.path(getwd(), "statistical_analysis/output/objects/irt")
    }
  },
  save_prefix = "IRT_batch"
) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- file.path(output_dir, "batch")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Build grid (all combinations)
  grid <- expand.grid(
    phq_model_id = names(phq_models),
    ssq_model_id = names(ssq_models),
    joint_model_id = names(joint_models),
    cutoff_method = cutoff_methods,
    target_specificity = target_specificities,
    target_sensitivity = target_sensitivities,
    cutoff_by_group_apply = apply_groups,
    apply_mode = apply_modes,
    stringsAsFactors = FALSE
  )
  
  # If apply_mode is 'single', hierarchy isn't used; if 'hierarchical', apply group is ignored.
  grid_targets <- grid
  
  results_overall <- list()
  failures <- list()
  
  run_one <- function(i) {
    p <- grid_targets[i, , drop = FALSE]
    
    # Choose apply strategy
    cutoff_by_group_apply <- if (identical(p$apply_mode, "hierarchical")) "hierarchical" else as.character(p$cutoff_by_group_apply)
    
    # Build stable run_id
    run_sig <- list(
      phq_model_id = p$phq_model_id,
      ssq_model_id = p$ssq_model_id,
      joint_model_id = p$joint_model_id,
      cutoff_method = p$cutoff_method,
      target_specificity = p$target_specificity,
      target_sensitivity = p$target_sensitivity,
      cutoff_by_group_apply = cutoff_by_group_apply,
      eval_groupings = eval_groupings,
      apply_hierarchy = apply_hierarchy,
      method_phq = method_phq,
      method_ssq = method_ssq,
      method_joint = method_joint,
      method_mg = method_mg
    )
    run_id <- .hash_id(run_sig, 16)
    
    out_file <- file.path(out_dir, paste0(save_prefix, "__", run_id, ".rds"))
    if (file.exists(out_file) && !isTRUE(overwrite)) {
      # Load and return summary only (fast)
      res <- readRDS(out_file)
    } else {
      res <- tryCatch({
        run_irt_harmonization_pipeline(
          phq_model = phq_models[[p$phq_model_id]],
          ssq_model = ssq_models[[p$ssq_model_id]],
          joint_model = joint_models[[p$joint_model_id]],
          method_phq = method_phq,
          method_ssq = method_ssq,
          method_joint = method_joint,
          technical_phq = technical_phq,
          technical_ssq = technical_ssq,
          technical_joint = technical_joint,
          method_mg = method_mg,
          technical_mg = technical_mg,
          cutoff_method = p$cutoff_method,
          target_sensitivity = as.numeric(p$target_sensitivity),
          target_specificity = as.numeric(p$target_specificity),
          cutoff_by_group_eval = eval_groupings,
          cutoff_by_group_apply = cutoff_by_group_apply,
          cutoff_apply_hierarchy = apply_hierarchy,
          run_sensitivity = run_sensitivity,
          # IMPORTANT for batch transparency: avoid overwriting canonical outputs
          save_results = FALSE,
          export_engine = FALSE,
          save_sensitivity_scores = FALSE,
          output_dir = output_dir
        )
      }, error = function(e) e)
      
      if (inherits(res, "error")) {
        return(list(status = "failed", run_id = run_id, error = res$message, params = p))
      }
      
      .saveRDS_atomic(res, out_file)
    }
    
    # Build ONE overall summary row per scenario (applied metrics + overall cutoff)
    ct <- res$link_results$cutoff_table
    cut_overall <- NA_real_
    if (!is.null(ct)) {
      row <- ct[ct$Grouping == "none" & ct$Group_Level == "Overall", , drop = FALSE]
      if (nrow(row) > 0) cut_overall <- as.numeric(row$Cutoff_Theta[1])
    }
    
    met <- res$link_results$metrics
    met_wide <- list(Accuracy = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_, PPV = NA_real_, NPV = NA_real_, Cohen_Kappa = NA_real_)
    if (!is.null(met) && all(c("Metric", "Value") %in% names(met))) {
      mv <- setNames(as.numeric(met$Value), as.character(met$Metric))
      if (!is.na(mv["Accuracy"])) met_wide$Accuracy <- mv["Accuracy"]
      if (!is.na(mv["Sensitivity"])) met_wide$Sensitivity <- mv["Sensitivity"]
      if (!is.na(mv["Specificity"])) met_wide$Specificity <- mv["Specificity"]
      if (!is.na(mv["PPV"])) met_wide$PPV <- mv["PPV"]
      if (!is.na(mv["NPV"])) met_wide$NPV <- mv["NPV"]
      if (!is.na(mv["Cohen_Kappa"])) met_wide$Cohen_Kappa <- mv["Cohen_Kappa"]
    }
    
    eval_summary <- data.frame(
      Run_ID = run_id,
      PHQ_Model = as.character(p$phq_model_id),
      SSQ_Model = as.character(p$ssq_model_id),
      JOINT_Model = as.character(p$joint_model_id),
      Cutoff_Method = as.character(p$cutoff_method),
      Apply_Strategy = as.character(cutoff_by_group_apply),
      Target_Specificity = as.numeric(p$target_specificity),
      Target_Sensitivity = as.numeric(p$target_sensitivity),
      Cutoff_Theta_Overall = cut_overall,
      Accuracy = met_wide$Accuracy,
      Sensitivity = met_wide$Sensitivity,
      Specificity = met_wide$Specificity,
      PPV = met_wide$PPV,
      NPV = met_wide$NPV,
      Cohen_Kappa = met_wide$Cohen_Kappa,
      Out_File = out_file,
      stringsAsFactors = FALSE
    )
    
    list(status = "ok", run_id = run_id, out_file = out_file, eval_summary = eval_summary)
  }
  
  # Run sequentially or in parallel
  idx <- seq_len(nrow(grid_targets))
  if (workers > 1L) {
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Parallel execution requested (workers > 1) but packages 'future' and/or 'future.apply' are not installed.")
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    out_list <- future.apply::future_lapply(idx, run_one, future.seed = TRUE)
  } else {
    out_list <- lapply(idx, run_one)
  }
  
  # Collect outputs
  for (x in out_list) {
    if (is.null(x)) next
    if (identical(x$status, "failed")) {
      failures[[x$run_id]] <- list(params = x$params, error = x$error)
      message("  [FAILED] ", x$run_id, " :: ", x$error)
    } else if (identical(x$status, "ok")) {
      if (!is.null(x$eval_summary) && nrow(x$eval_summary) > 0) {
        results_overall[[x$run_id]] <- x$eval_summary
      }
    }
  }
  
  overall_df <- do.call(rbind, results_overall)
  
  # Save consolidated summaries
  summary_rds <- file.path(out_dir, paste0(save_prefix, "__SUMMARY_overall.rds"))
  saveRDS(overall_df, summary_rds)
  
  summary_csv <- file.path(out_dir, paste0(save_prefix, "__SUMMARY_overall.csv"))
  try(write.csv(overall_df, summary_csv, row.names = FALSE), silent = TRUE)
  
  failures_rds <- file.path(out_dir, paste0(save_prefix, "__FAILURES.rds"))
  saveRDS(failures, failures_rds)
  
  # Optional: GT table
  tbl <- NULL
  if (!is.null(overall_df) && nrow(overall_df) > 0) {
    tbl <- gt(overall_df, groupname_col = "Run_ID") %>%
      tab_header(title = "Batch Cutoff Evaluation (Overall Rows Only)") %>%
      fmt_number(columns = any_of(c("Cutoff_Theta_Overall", "Cohen_Kappa")), decimals = 3) %>%
      fmt_percent(columns = any_of(c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV")), decimals = 1)
  }
  
  return(list(
    grid = grid_targets,
    summary_overall = overall_df,
    failures = failures,
    files = list(summary_rds = summary_rds, summary_csv = summary_csv, failures_rds = failures_rds),
    gt_overall = tbl
  ))
}

# ==============================================================================
# F. POST-BATCH: REBUILD SCENARIO REGISTRY (safe to run anytime)
# ==============================================================================
# This helper scans /batch/ for completed scenario .rds files, extracts the same
# overall-row summary used by batch_run_irt(), and writes scenario_registry.csv.
#
# It is designed to work even if a long batch was interrupted. It will only
# include scenarios that have already been saved.

rebuild_scenario_registry <- function(batch_dir, save_csv = TRUE, save_rds = TRUE, save_prefix = "IRT_batch") {
  if (is.null(batch_dir) || !dir.exists(batch_dir)) stop("batch_dir does not exist: ", batch_dir)
  
  # Only scenario files (exclude SUMMARY/FAILURES). Also ignore temp files.
  patt <- paste0("^", save_prefix, "__([a-f0-9]{10,64})\\.rds$")
  files <- list.files(batch_dir, pattern = patt, full.names = TRUE)
  
  files <- files[!grepl("__SUMMARY_", basename(files))]
  files <- files[!grepl("__FAILURES\\.rds$", basename(files))]
  files <- files[!grepl("\\.tmp$", basename(files))]
  
  if (length(files) == 0) {
    reg <- data.frame(
      Run_ID = character(),
      PHQ_Model = character(),
      SSQ_Model = character(),
      JOINT_Model = character(),
      Cutoff_Method = character(),
      Apply_Strategy = character(),
      Target_Specificity = numeric(),
      Target_Sensitivity = numeric(),
      Cutoff_Theta_Overall = numeric(),
      Accuracy = numeric(),
      Sensitivity = numeric(),
      Specificity = numeric(),
      PPV = numeric(),
      NPV = numeric(),
      Cohen_Kappa = numeric(),
      Out_File = character(),
      Status = character(),
      Error = character(),
      stringsAsFactors = FALSE
    )
    
    if (isTRUE(save_csv)) utils::write.csv(reg, file.path(batch_dir, "scenario_registry.csv"), row.names = FALSE)
    if (isTRUE(save_rds)) saveRDS(reg, file.path(batch_dir, "scenario_registry.rds"))
    return(reg)
  }
  
  extract_one <- function(f) {
    run_id <- sub(paste0("^", save_prefix, "__"), "", basename(f))
    run_id <- sub("\\.rds$", "", run_id)
    
    res <- tryCatch(readRDS(f), error = function(e) e)
    if (inherits(res, "error")) {
      return(data.frame(
        Run_ID = run_id,
        PHQ_Model = NA_character_,
        SSQ_Model = NA_character_,
        JOINT_Model = NA_character_,
        Cutoff_Method = NA_character_,
        Apply_Strategy = NA_character_,
        Target_Specificity = NA_real_,
        Target_Sensitivity = NA_real_,
        Cutoff_Theta_Overall = NA_real_,
        Accuracy = NA_real_,
        Sensitivity = NA_real_,
        Specificity = NA_real_,
        PPV = NA_real_,
        NPV = NA_real_,
        Cohen_Kappa = NA_real_,
        Out_File = f,
        Status = "failed_read",
        Error = res$message,
        stringsAsFactors = FALSE
      ))
    }
    
    # Recover parameters (best-effort)
    istruct <- res$input_structure
    phq_id <- NA_character_; ssq_id <- NA_character_; joint_id <- NA_character_
    cutoff_m <- NA_character_; apply_s <- NA_character_
    tgt_spec <- NA_real_; tgt_sens <- NA_real_
    
    if (!is.null(istruct) && is.list(istruct)) {
      phq_id <- if (!is.null(istruct$phq_model)) as.character(istruct$phq_model) else NA_character_
      ssq_id <- if (!is.null(istruct$ssq_model)) as.character(istruct$ssq_model) else NA_character_
      joint_id <- if (!is.null(istruct$joint_model)) as.character(istruct$joint_model) else NA_character_
      cutoff_m <- if (!is.null(istruct$cutoff_method)) as.character(istruct$cutoff_method) else NA_character_
      apply_s <- if (!is.null(istruct$cutoff_by_group_apply)) as.character(istruct$cutoff_by_group_apply) else NA_character_
      tgt_spec <- if (!is.null(istruct$target_specificity)) suppressWarnings(as.numeric(istruct$target_specificity)) else NA_real_
      tgt_sens <- if (!is.null(istruct$target_sensitivity)) suppressWarnings(as.numeric(istruct$target_sensitivity)) else NA_real_
    }
    
    # Overall cutoff
    cut_overall <- NA_real_
    ct <- res$link_results$cutoff_table
    if (!is.null(ct) && all(c("Grouping", "Group_Level", "Cutoff_Theta") %in% names(ct))) {
      row <- ct[ct$Grouping == "none" & ct$Group_Level == "Overall", , drop = FALSE]
      if (nrow(row) > 0) cut_overall <- suppressWarnings(as.numeric(row$Cutoff_Theta[1]))
    }
    
    # Applied metrics
    acc <- sens <- spec <- ppv <- npv <- kap <- NA_real_
    met <- res$link_results$metrics
    if (!is.null(met) && all(c("Metric", "Value") %in% names(met))) {
      mv <- setNames(suppressWarnings(as.numeric(met$Value)), as.character(met$Metric))
      if (!is.na(mv["Accuracy"])) acc <- mv["Accuracy"]
      if (!is.na(mv["Sensitivity"])) sens <- mv["Sensitivity"]
      if (!is.na(mv["Specificity"])) spec <- mv["Specificity"]
      if (!is.na(mv["PPV"])) ppv <- mv["PPV"]
      if (!is.na(mv["NPV"])) npv <- mv["NPV"]
      if (!is.na(mv["Cohen_Kappa"])) kap <- mv["Cohen_Kappa"]
    }
    
    data.frame(
      Run_ID = run_id,
      PHQ_Model = phq_id,
      SSQ_Model = ssq_id,
      JOINT_Model = joint_id,
      Cutoff_Method = cutoff_m,
      Apply_Strategy = apply_s,
      Target_Specificity = tgt_spec,
      Target_Sensitivity = tgt_sens,
      Cutoff_Theta_Overall = cut_overall,
      Accuracy = acc,
      Sensitivity = sens,
      Specificity = spec,
      PPV = ppv,
      NPV = npv,
      Cohen_Kappa = kap,
      Out_File = f,
      Status = "ok",
      Error = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  reg_list <- lapply(files, extract_one)
  reg <- do.call(rbind, reg_list)
  
  num_cols <- intersect(c("Target_Specificity", "Target_Sensitivity", "Cutoff_Theta_Overall", "Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "Cohen_Kappa"), names(reg))
  for (cc in num_cols) reg[[cc]] <- suppressWarnings(as.numeric(reg[[cc]]))
  
  if (isTRUE(save_csv)) utils::write.csv(reg, file.path(batch_dir, "scenario_registry.csv"), row.names = FALSE)
  if (isTRUE(save_rds)) saveRDS(reg, file.path(batch_dir, "scenario_registry.rds"))
  
  reg
}

# ==============================================================================
# D. POST-BATCH: RANK + EXPORT A SELECTED SCENARIO (for transparency → best pick)
# ==============================================================================

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

export_selected_scenario <- function(batch_dir,
                                     run_id,
                                     export_dir = {
                                       if (requireNamespace("here", quietly = TRUE)) {
                                         here::here("statistical_analysis/output/objects/irt")
                                       } else {
                                         file.path(getwd(), "statistical_analysis/output/objects/irt")
                                       }
                                     },
                                     tag = "SELECTED",
                                     overwrite = TRUE) {
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Find the scenario RDS produced by batch_run_irt()
  cand <- list.files(batch_dir, pattern = paste0("__", run_id, "\\.rds$"), full.names = TRUE)
  if (length(cand) == 0) stop("Could not find scenario file for run_id in batch_dir.")
  if (length(cand) > 1) cand <- cand[1]
  
  res <- readRDS(cand)
  if (is.null(res$scoring_engine) || is.null(res$models$joint)) {
    stop("Scenario object does not contain scoring_engine and/or joint model.")
  }
  
  eng_file <- file.path(export_dir, paste0("SSQ10_Harmonized_Scoring_Engine__", tag, ".rds"))
  joint_file <- file.path(export_dir, paste0("Joint_IRT_Model_For_Scoring__", tag, ".rds"))
  meta_file <- file.path(export_dir, paste0("SELECTED_SCENARIO_META__", tag, ".rds"))
  
  if (!overwrite) {
    if (file.exists(eng_file) || file.exists(joint_file) || file.exists(meta_file)) {
      stop("Export files already exist and overwrite = FALSE.")
    }
  }
  
  saveRDS(res$scoring_engine, eng_file)
  saveRDS(res$models$joint, joint_file)
  
  meta <- list(
    run_id = run_id,
    source_file = cand,
    input_structure = res$input_structure,
    threshold_applied = res$threshold,
    metrics_applied = res$metrics,
    batch_exported_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  saveRDS(meta, meta_file)
  
  invisible(list(engine = eng_file, joint_model = joint_file, meta = meta_file))
}

# ------------------------------------------------------------------------------
# Example batch run (uncomment to use)
# ------------------------------------------------------------------------------
# Notes:
# - This produces a transparent comparison table across ALL scenarios (overall rows only).
# - It DOES NOT overwrite a single scoring engine during the batch.
# - After reviewing the summary, pick the best run_id and export that scenario’s engine/model.
#
# # Choose parallel workers (leave at least 2 cores for the system)
# workers <- max(1L, min(8L, parallel::detectCores() - 2L))
#
# # Run the grid
# batch_res <- batch_run_irt(
#   cutoff_methods = c("model_tcc", "youden", "sens_at_spec", "spec_at_sens"),
#   target_specificities = c(0.90, 0.95, 0.98),
#   target_sensitivities = c(0.70, 0.80, 0.90),
#   apply_modes = c("single", "hierarchical"),
#   apply_groups = c("none", "SEX", "AGEGRP", "SEX_AGEGRP"),
#   # Evaluate ALL groupings for transparency in every scenario
#   eval_groupings = c("none", "SEX", "AGEGRP", "SEX_AGEGRP"),
#   # Use QMCEM automatically for 3+ dimensional models
#   method_joint = "auto",
#   method_mg = "auto",
#   workers = workers,
#   run_sensitivity = FALSE,
#   overwrite = TRUE
# )
#
# # View comparison table
# batch_res$gt_overall
#
# # Rank scenarios (example: maximize sensitivity subject to specificity ≥ 0.95)
# top <- rank_batch_scenarios(
#   summary_df = batch_res$summary_overall,
#   evaluation = "none",
#   min_specificity = 0.95,
#   maximize = "Sensitivity",
#   top_n = 10
# )
# top
#
# # Export the selected best scenario’s scoring engine + joint model
# # (batch_dir can be inferred from the summary file path)
# batch_dir <- dirname(batch_res$files$summary_rds)
# best_run_id <- top$Run_ID[1]
# export_selected_scenario(
#   batch_dir = batch_dir,
#   run_id = best_run_id,
#   export_dir = here("statistical_analysis/output/objects/irt"),
#   tag = "BEST_SPEC95_MAXSENS",
#   overwrite = TRUE
# )

# ==============================================================================
# Execute
# ==============================================================================
if (sys.nframe() == 0) {
  # Example: Run with default 1-factor models and Sensitivity Check enabled
  res <- run_irt_harmonization_pipeline(
    cutoff_method = "model_tcc",
    cutoff_by_group_eval = c("none", "SEX", "AGEGRP", "SEX_AGEGRP"),
    cutoff_by_group_apply = "none",
    run_sensitivity = TRUE
  )
}
