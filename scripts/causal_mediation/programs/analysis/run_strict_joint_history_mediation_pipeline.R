#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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
seed <- as.integer(arg_value("--seed", "20260414"))
n_sim <- as.integer(arg_value("--sim", "300"))
block_sim <- as.integer(arg_value("--block-sim", "120"))
expanded_block_sim <- as.integer(arg_value("--expanded-block-sim", "0"))

script_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][[1]])
script_arg <- gsub("~\\+~", " ", script_arg)
script_path <- normalizePath(script_arg, mustWork = TRUE)
workspace_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

config_path <- file.path(workspace_root, "config", "adjustment_sets_final_inferential.csv")
imputation_dir <- file.path(workspace_root, "input", "imputations", "strict_jmaes_final")
listing_dir <- file.path(workspace_root, "output", "listings")
dir.create(listing_dir, recursive = TRUE, showWarnings = FALSE)

build_script <- file.path(workspace_root, "programs", "analysis", "build_strict_joint_history_baseline_imputations.R")
fit_script <- file.path(workspace_root, "programs", "analysis", "fit_joint_history_mediation.R")
pool_script <- file.path(workspace_root, "programs", "analysis", "pool_joint_history_mediation_results.R")

run_rscript <- function(script, args_vec) {
  status <- system2("Rscript", c(shQuote(script), shQuote(args_vec)))
  if (!identical(status, 0L)) {
    stop("Command failed for ", basename(script), " with status ", status)
  }
}

run_rscript(
  build_script,
  c(
    paste0("--m=", m),
    paste0("--seed=", seed)
  )
)

for (i in seq_len(m)) {
  subject_path <- file.path(
    imputation_dir,
    sprintf("jmaes_strict_subject_level_final_imp_%02d.rds", i)
  )
  output_prefix <- sprintf("strict_jmaes_final_imp_%02d", i)

  run_rscript(
    fit_script,
    c(
      "--branch=strict",
      paste0("--subject-input=", subject_path),
      paste0("--adjustment-config=", config_path),
      paste0("--output-prefix=", output_prefix),
      paste0("--sim=", n_sim),
      paste0("--block-sim=", block_sim),
      paste0("--expanded-block-sim=", expanded_block_sim),
      paste0("--seed=", seed + i)
    )
  )
}

run_rscript(
  pool_script,
  c(
    "--input-prefix=strict_jmaes_final_imp_",
    "--output-prefix=strict_jmaes_final",
    paste0("--m=", m)
  )
)

write_csv(
  tibble(
    run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    script = script_path,
    workspace_root = workspace_root,
    m_imputations = m,
    seed = seed,
    n_simulation_draws = n_sim,
    n_block_simulation_draws = block_sim,
    n_expanded_block_simulation_draws = expanded_block_sim,
    pooled_output_prefix = "strict_jmaes_final"
  ),
  file.path(listing_dir, "strict_jmaes_final_pipeline_manifest.csv")
)

message("Completed strict JMAES final inferential pipeline.")
