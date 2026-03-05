#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 0) {
  stop("Unable to resolve script path from commandArgs().")
}

this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(this_file), "..", ".."), winslash = "/", mustWork = TRUE)
objects_dir <- file.path(project_root, "_tools", "_objects")
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)

resolve_engine_bundle <- function(root) {
  env_path <- Sys.getenv("COLUMINATE_ENGINES_PATH", unset = "")
  candidates <- c(
    env_path,
    file.path(root, "_tools", "_objects", "irt_joint_models.rds"),
    file.path(
      root, "..", "OBJ00-Datasets Preperation", "data_management",
      "data_examples", "1_staging_snippets", "derived_models", "irt_joint_models.rds"
    )
  )
  candidates <- unique(candidates[nzchar(candidates)])
  idx <- which(file.exists(candidates))[1]
  if (is.na(idx)) {
    stop(
      "Unable to locate irt_joint_models.rds.\n",
      "Set COLUMINATE_ENGINES_PATH or ensure OBJ00 derived model exists."
    )
  }
  normalizePath(candidates[idx], winslash = "/", mustWork = TRUE)
}

source_path <- resolve_engine_bundle(project_root)
target_path <- normalizePath(
  file.path(objects_dir, "irt_joint_models.rds"),
  winslash = "/",
  mustWork = FALSE
)

if (!identical(source_path, target_path)) {
  file.copy(source_path, target_path, overwrite = TRUE)
}

bundle <- readRDS(target_path)
if (!is.list(bundle) || length(bundle) == 0) {
  stop("Copied bundle appears invalid at: ", target_path)
}

cat("Prepared API deploy bundle.\n")
cat("Source:", source_path, "\n")
cat("Target:", target_path, "\n")
cat("Engines:", length(bundle), "\n")
