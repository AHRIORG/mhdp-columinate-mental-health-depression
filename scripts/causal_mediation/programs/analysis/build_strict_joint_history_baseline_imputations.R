#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(mice)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(paste0("^", flag, "="), "", hit[[1]])
}

m <- as.integer(arg_value("--m", "20"))
maxit <- as.integer(arg_value("--maxit", "20"))
seed <- as.integer(arg_value("--seed", "20260414"))

script_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][[1]])
script_arg <- gsub("~\\+~", " ", script_arg)
script_path <- normalizePath(script_arg, mustWork = TRUE)
workspace_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

input_dir <- file.path(workspace_root, "input")
listing_dir <- file.path(workspace_root, "output", "listings")
imputation_dir <- file.path(input_dir, "imputations", "strict_jmaes_final")
dir.create(imputation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(listing_dir, recursive = TRUE, showWarnings = FALSE)

subject_path <- file.path(input_dir, "jmaes_strict_subject_level_analysis.rds")
transition_path <- file.path(input_dir, "jmaes_strict_transition_bundle.rds")

subject <- readRDS(subject_path) %>%
  mutate(
    USUBJID = as.character(.data$USUBJID),
    SEX = factor(.data$SEX),
    RURAL = factor(.data$RURAL),
    MEDU = ordered(.data$MEDU, levels = c("None/Primary", "Secondary", "Completed Matric")),
    MHIV = factor(.data$MHIV, levels = c("Negative", "Positive")),
    MPTS = factor(.data$MPTS, levels = c("No partner", "Casual", "Regular", "Married", "Widowed/Separated")),
    MMIG = factor(.data$MMIG, levels = c("No", "Yes"))
  )

transition_bundle <- readRDS(transition_path)
overlap_ids <- as.character(transition_bundle$inputs$overlap_ids)

modal_lookup <- transition_bundle$saved_probs %>%
  transmute(
    USUBJID = as.character(.data$USUBJID),
    EARLY_MODAL = factor(.data$C1),
    MID_MODAL = factor(.data$C2),
    JOINT_MODAL = factor(.data$MLCJOINT)
  )

analysis_df <- subject %>%
  filter(.data$USUBJID %in% overlap_ids) %>%
  left_join(modal_lookup, by = "USUBJID") %>%
  arrange(.data$USUBJID)

impute_vars <- c("MAGE", "MEDU", "MHIV", "MPTS", "MMIG")
predictor_vars <- c(
  "SEX", "AGE", "RURAL", "CLYR",
  "DNKA", "ESXC", "FDSC", "GOVG", "VLNC", "SCSP", "SCHL", "DPBN",
  "JOINT_MODAL",
  impute_vars
)

missingness_summary <- tibble(
  variable = impute_vars,
  n_missing = vapply(impute_vars, function(var) sum(is.na(analysis_df[[var]])), numeric(1)),
  pct_missing = vapply(impute_vars, function(var) mean(is.na(analysis_df[[var]])), numeric(1))
) %>%
  mutate(pct_missing = round(100 * .data$pct_missing, 1))

write_csv(
  missingness_summary,
  file.path(listing_dir, "strict_jmaes_final_imputation_missingness.csv")
)

imp_data <- analysis_df %>%
  select("USUBJID", all_of(predictor_vars))

methods <- rep("", ncol(imp_data))
names(methods) <- names(imp_data)
methods[c("MAGE", "MEDU", "MHIV", "MPTS", "MMIG")] <- c("pmm", "cart", "cart", "cart", "cart")

predictor_matrix <- matrix(0, nrow = ncol(imp_data), ncol = ncol(imp_data))
rownames(predictor_matrix) <- colnames(predictor_matrix) <- names(imp_data)

for (target in impute_vars) {
  predictor_matrix[target, setdiff(predictor_vars, target)] <- 1
}

predictor_matrix[, "USUBJID"] <- 0
predictor_matrix["USUBJID", ] <- 0

ini <- mice(
  imp_data,
  m = 1,
  maxit = 0,
  method = methods,
  predictorMatrix = predictor_matrix,
  printFlag = FALSE
)

imp <- mice(
  imp_data,
  m = m,
  maxit = maxit,
  method = ini$method,
  predictorMatrix = ini$predictorMatrix,
  seed = seed,
  printFlag = TRUE
)

manifest <- vector("list", m)

fill_residual_missing <- function(data, donors, vars) {
  residual_fills <- 0L

  for (var in vars) {
    miss_idx <- which(is.na(data[[var]]))
    if (!length(miss_idx)) {
      next
    }

    for (idx in miss_idx) {
      donor_pool <- donors[[var]][!is.na(donors[[var]]) & donors$JOINT_MODAL == data$JOINT_MODAL[[idx]]]
      if (!length(donor_pool)) {
        donor_pool <- donors[[var]][!is.na(donors[[var]]) & donors$SEX == data$SEX[[idx]]]
      }
      if (!length(donor_pool)) {
        donor_pool <- donors[[var]][!is.na(donors[[var]])]
      }
      if (!length(donor_pool)) {
        next
      }
      data[[var]][idx] <- sample(donor_pool, size = 1)
      residual_fills <- residual_fills + 1L
    }
  }

  list(data = data, residual_fills = residual_fills)
}

for (i in seq_len(m)) {
  completed <- complete(imp, action = i) %>%
    as_tibble() %>%
    mutate(
      USUBJID = as.character(.data$USUBJID),
      SEX = factor(.data$SEX),
      RURAL = factor(.data$RURAL),
      MEDU = ordered(.data$MEDU, levels = c("None/Primary", "Secondary", "Completed Matric")),
      MHIV = factor(.data$MHIV, levels = c("Negative", "Positive")),
      MPTS = factor(.data$MPTS, levels = c("No partner", "Casual", "Regular", "Married", "Widowed/Separated")),
      MMIG = factor(.data$MMIG, levels = c("No", "Yes"))
    )

  residual_fill <- fill_residual_missing(completed, analysis_df, impute_vars)
  completed <- residual_fill$data

  output_path <- file.path(
    imputation_dir,
    sprintf("jmaes_strict_subject_level_final_imp_%02d.rds", i)
  )
  saveRDS(completed, output_path)

  manifest[[i]] <- tibble(
    imputation_id = i,
    output_path = output_path,
    n_subjects = nrow(completed),
    n_missing_after = sum(is.na(completed[, impute_vars])),
    n_residual_fills = residual_fill$residual_fills
  )
}

saveRDS(imp, file.path(imputation_dir, "jmaes_strict_baseline_mids.rds"))

write_csv(
  bind_rows(manifest),
  file.path(listing_dir, "strict_jmaes_final_imputation_manifest.csv")
)

write_csv(
  tibble(
    run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    script = script_path,
    workspace_root = workspace_root,
    source_subject_path = subject_path,
    source_transition_path = transition_path,
    n_overlap_subjects = length(overlap_ids),
    m_imputations = m,
    max_iterations = maxit,
    seed = seed
  ),
  file.path(listing_dir, "strict_jmaes_final_imputation_run_manifest.csv")
)

message("Wrote strict final imputation manifest to ", file.path(listing_dir, "strict_jmaes_final_imputation_manifest.csv"))
