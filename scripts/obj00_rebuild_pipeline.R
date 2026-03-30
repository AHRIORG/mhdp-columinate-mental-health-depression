#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# OBJ00 Public Rebuild Pipeline Runner
# Purpose:
#   1) Build OBJ00 core datasets from private inputs
#   2) Generate public-safe synthetic exemplar datasets
#   3) Validate reproducibility against ADaM reference objects when available
# -------------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (is.character(x) && length(x) == 1 && identical(x, "")) return(y)
  x
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0) return(default)
  val <- tolower(trimws(as.character(x)[1]))
  if (val %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

parse_args <- function(argv) {
  out <- list()
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    if (!grepl("=", kv, fixed = TRUE)) {
      out[[kv]] <- "true"
      next
    }
    parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    val <- paste(parts[-1], collapse = "=")
    out[[key]] <- val
  }
  out
}

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0) stop("Cannot resolve script path from commandArgs.")
  normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = TRUE)
}

find_public_repo_root <- function() {
  sdir <- dirname(script_path())
  root <- normalizePath(file.path(sdir, ".."), winslash = "/", mustWork = TRUE)
  probe <- file.path(
    root, "scripts", "obj00", "data_management", "scripts", "01_cleaning",
    "01b_build_core_datasets_candidate.R"
  )
  if (!file.exists(probe)) {
    stop("Could not find public OBJ00 script bundle from script location: ", sdir)
  }
  root
}

normalize_path_safe <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

rel_path <- function(from, to) {
  from <- normalizePath(from, winslash = "/", mustWork = TRUE)
  to <- normalizePath(to, winslash = "/", mustWork = FALSE)
  parts_from <- strsplit(from, "/", fixed = TRUE)[[1]]
  parts_to <- strsplit(to, "/", fixed = TRUE)[[1]]

  i <- 1
  limit <- min(length(parts_from), length(parts_to))
  while (i <= limit && identical(parts_from[i], parts_to[i])) i <- i + 1

  up <- rep("..", length(parts_from) - i + 1)
  down <- if (i <= length(parts_to)) parts_to[i:length(parts_to)] else character()
  rel <- file.path(c(up, down))
  if (identical(rel, "")) "." else rel
}

safe_unlink_dir <- function(path) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
}

probe_readable <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (file.access(path, 4) != 0) return(FALSE)
  ok <- tryCatch({
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    readBin(con, what = "raw", n = 16)
    TRUE
  }, error = function(e) FALSE)
  ok
}

resolve_candidate_file <- function(base_dir, rel_path_expected) {
  expected <- file.path(base_dir, rel_path_expected)
  if (probe_readable(expected)) return(expected)

  dir_part <- dirname(expected)
  fname <- basename(expected)
  stem <- sub("\\.[^.]+$", "", fname)
  ext <- sub("^.*(\\.[^.]+)$", "\\1", fname)

  local_candidates <- c(
    file.path(dir_part, fname),
    file.path(dir_part, paste0(stem, "-MACDBNW2N4T", ext))
  )
  local_candidates <- unique(local_candidates[file.exists(local_candidates)])
  for (cand in local_candidates) {
    if (probe_readable(cand)) return(cand)
  }

  patt <- paste0("^", stem, "(-[A-Za-z0-9]+)?", gsub("\\.", "\\\\.", ext), "$")
  recursive_hits <- list.files(
    base_dir, pattern = patt, recursive = TRUE, full.names = TRUE, ignore.case = FALSE
  )
  for (cand in unique(recursive_hits)) {
    if (probe_readable(cand)) return(cand)
  }

  NA_character_
}

stage_private_inputs <- function(obj00_root, private_data_dir, include_final_steps = TRUE) {
  required <- c(
    "raw_data/columinate_clean.dta",
    "raw_data/hdss_raw/SurveillanceEpisodesYrAgeDel_HSE_Labour_Partner_Parents_HIV.dta",
    "raw_data/hdss_raw/HouseholdMap.dta",
    "raw_data/hdss_raw/HSE_scores.csv"
  )
  if (isTRUE(include_final_steps)) {
    required <- c(required, "adam/dt_psychometric.RData", "adam/dt_dreams_multi_ssq.RData")
  }

  resolved <- setNames(
    vapply(required, function(rp) resolve_candidate_file(private_data_dir, rp), FUN.VALUE = character(1)),
    required
  )

  missing <- names(resolved)[is.na(resolved)]
  if (length(missing) > 0) {
    stop(
      "Could not resolve required input files from private data dir. Missing: ",
      paste0(file.path(private_data_dir, missing), collapse = "; ")
    )
  }

  stage_root <- file.path(obj00_root, ".tmp_obj00_pipeline_private_data")
  safe_unlink_dir(stage_root)
  dir.create(stage_root, recursive = TRUE, showWarnings = FALSE)

  for (rp in names(resolved)) {
    target <- file.path(stage_root, rp)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    ok <- file.symlink(from = resolved[[rp]], to = target)
    if (!isTRUE(ok)) stop("Failed to create symlink: ", target, " -> ", resolved[[rp]])
  }

  list(stage_root = stage_root, resolved = resolved)
}

load_rdata_object <- function(path, object_name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Object '", object_name, "' not found in ", path)
  }
  get(object_name, envir = env, inherits = FALSE)
}

get_key_md5 <- function(df, key_cols) {
  use_cols <- intersect(key_cols, names(df))
  if (length(use_cols) == 0) return(NA_character_)
  key_df <- df[, use_cols, drop = FALSE]
  ord <- do.call(order, c(unname(key_df), list(na.last = TRUE)))
  key_df <- key_df[ord, , drop = FALSE]
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(key_df, file = tmp, row.names = FALSE, na = "")
  as.character(tools::md5sum(tmp))
}

compare_to_reference <- function(results, ref_dir) {
  out <- list(
    reference_found = FALSE,
    childhood = list(),
    multistudies = list()
  )

  child_file <- file.path(ref_dir, "dt_childhood_exposure.RData")
  multi_file <- file.path(ref_dir, "dt_multistudies.RData")

  if (!file.exists(child_file) || !file.exists(multi_file)) return(out)
  out$reference_found <- TRUE

  ref_child <- load_rdata_object(child_file, "dt_childhood_exposure")
  ref_multi <- load_rdata_object(multi_file, "dt_multistudies")

  out$childhood$rows_match <- identical(nrow(ref_child), nrow(results$dt_childhood_exposure))
  out$childhood$cols_match <- identical(ncol(ref_child), ncol(results$dt_childhood_exposure))
  out$childhood$names_match <- identical(names(ref_child), names(results$dt_childhood_exposure))
  out$childhood$key_md5_ref <- get_key_md5(ref_child, c("USUBJID", "HHID", "HDSSYR", "EXPAGE"))
  out$childhood$key_md5_new <- get_key_md5(results$dt_childhood_exposure, c("USUBJID", "HHID", "HDSSYR", "EXPAGE"))
  out$childhood$key_md5_match <- identical(out$childhood$key_md5_ref, out$childhood$key_md5_new)

  out$multistudies$rows_match <- identical(nrow(ref_multi), nrow(results$dt_multistudies))
  out$multistudies$cols_match <- identical(ncol(ref_multi), ncol(results$dt_multistudies))
  out$multistudies$names_match <- identical(names(ref_multi), names(results$dt_multistudies))
  out$multistudies$key_md5_ref <- get_key_md5(ref_multi, c("USUBJID", "STUDY", "VISITDT", "DPBN"))
  out$multistudies$key_md5_new <- get_key_md5(results$dt_multistudies, c("USUBJID", "STUDY", "VISITDT", "DPBN"))
  out$multistudies$key_md5_match <- identical(out$multistudies$key_md5_ref, out$multistudies$key_md5_new)

  out
}

synth_values <- function(col_name, x, n) {
  upper <- toupper(col_name)

  if (inherits(x, "Date")) return(as.Date("2000-01-01") + seq_len(n) - 1L)
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) {
    base <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
    return(base + (seq_len(n) - 1L) * 86400)
  }
  if (is.factor(x)) {
    lev <- levels(x)
    if (length(lev) == 0) lev <- c("No", "Yes")
    return(factor(rep_len(lev, n), levels = lev))
  }
  if (is.logical(x)) return(rep_len(c(TRUE, FALSE), n))
  if (grepl("USUBJID", upper)) return(sprintf("SYNTH-%04d", seq_len(n)))
  if (grepl("HHID", upper)) return(sprintf("HH-%04d", ((seq_len(n) - 1L) %% 50L) + 1L))
  if (grepl("DATE|DT$", upper)) {
    return(format(as.Date("2000-01-01") + seq_len(n) - 1L, "%Y-%m-%d"))
  }
  if (grepl("SEX", upper)) return(rep_len(c("Female", "Male"), n))
  if (grepl("STUDY", upper)) return(rep_len(c("DREAMS", "Multilevel", "Isisekelo Sempilo"), n))
  if (grepl("BIN|OUTCOME|YES|NO|FLAG", upper)) return(rep_len(c("No", "Yes"), n))
  if (is.integer(x)) return(as.integer(seq_len(n)))
  if (is.numeric(x)) return(as.numeric(seq_len(n)))
  if (is.character(x)) return(sprintf("%s_%03d", col_name, seq_len(n)))
  rep_len(NA, n)
}

build_synthetic_example <- function(df, n = 20L) {
  n <- max(1L, as.integer(n))
  cols <- names(df)
  out <- as.data.frame(
    setNames(
      lapply(cols, function(col) synth_values(col, df[[col]], n)),
      cols
    ),
    stringsAsFactors = FALSE
  )
  out
}

compare_two_runs <- function(run_a, run_b) {
  list(
    childhood = list(
      rows_match = identical(nrow(run_a$dt_childhood_exposure), nrow(run_b$dt_childhood_exposure)),
      cols_match = identical(ncol(run_a$dt_childhood_exposure), ncol(run_b$dt_childhood_exposure)),
      names_match = identical(names(run_a$dt_childhood_exposure), names(run_b$dt_childhood_exposure)),
      key_md5_a = get_key_md5(run_a$dt_childhood_exposure, c("USUBJID", "HHID", "HDSSYR", "EXPAGE")),
      key_md5_b = get_key_md5(run_b$dt_childhood_exposure, c("USUBJID", "HHID", "HDSSYR", "EXPAGE"))
    ),
    multistudies = list(
      rows_match = identical(nrow(run_a$dt_multistudies), nrow(run_b$dt_multistudies)),
      cols_match = identical(ncol(run_a$dt_multistudies), ncol(run_b$dt_multistudies)),
      names_match = identical(names(run_a$dt_multistudies), names(run_b$dt_multistudies)),
      key_md5_a = get_key_md5(run_a$dt_multistudies, c("USUBJID", "STUDY", "VISITDT", "DPBN")),
      key_md5_b = get_key_md5(run_b$dt_multistudies, c("USUBJID", "STUDY", "VISITDT", "DPBN"))
    )
  )
}

write_examples <- function(results, out_dir, n_rows = 20L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  child_ex <- build_synthetic_example(results$dt_childhood_exposure, n = n_rows)
  multi_ex <- build_synthetic_example(results$dt_multistudies, n = n_rows)
  utils::write.csv(child_ex, file.path(out_dir, "dt_childhood_exposure_example_synth.csv"), row.names = FALSE)
  utils::write.csv(multi_ex, file.path(out_dir, "dt_multistudies_example_synth.csv"), row.names = FALSE)

  if (!is.null(results$dt_dreams_multi_ssq)) {
    dreams_ex <- build_synthetic_example(results$dt_dreams_multi_ssq, n = n_rows)
    utils::write.csv(dreams_ex, file.path(out_dir, "dt_dreams_multi_ssq_example_synth.csv"), row.names = FALSE)
  }

  manifest <- data.frame(
    dataset = c("dt_childhood_exposure", "dt_multistudies", "dt_dreams_multi_ssq"),
    source_rows = c(
      nrow(results$dt_childhood_exposure),
      nrow(results$dt_multistudies),
      if (is.null(results$dt_dreams_multi_ssq)) NA_integer_ else nrow(results$dt_dreams_multi_ssq)
    ),
    source_cols = c(
      ncol(results$dt_childhood_exposure),
      ncol(results$dt_multistudies),
      if (is.null(results$dt_dreams_multi_ssq)) NA_integer_ else ncol(results$dt_dreams_multi_ssq)
    ),
    example_rows = n_rows,
    synthetic = TRUE,
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(out_dir, "obj00_examples_manifest.csv"), row.names = FALSE)

  readme <- c(
    "# OBJ00 Synthetic Example Datasets",
    "",
    "These files are synthetic placeholders generated from dataset schema only.",
    "They do not contain participant-level records and are intended for public examples.",
    "",
    "Generated files:",
    "- dt_childhood_exposure_example_synth.csv",
    "- dt_multistudies_example_synth.csv",
    "- dt_dreams_multi_ssq_example_synth.csv (if available)",
    "- obj00_examples_manifest.csv"
  )
  writeLines(readme, con = file.path(out_dir, "README.md"))
}

write_validation_report <- function(report_file, cfg, results, compare_ref, compare_self) {
  lines <- c(
    "# OBJ00 Pipeline Validation Report",
    "",
    paste0("- Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- include_index_final_steps: ", cfg$include_index_final_steps),
    paste0("- save_to_private_adam: ", cfg$save_to_private_adam),
    paste0("- self_repro_check: ", cfg$self_repro_check),
    paste0("- Generated dt_childhood_exposure: ", nrow(results$dt_childhood_exposure), " x ", ncol(results$dt_childhood_exposure)),
    paste0("- Generated dt_multistudies: ", nrow(results$dt_multistudies), " x ", ncol(results$dt_multistudies)),
    paste0("- Best calibration: ", results$best_calibration %||% "N/A"),
    ""
  )

  if (!isTRUE(compare_ref$reference_found)) {
    lines <- c(lines, "## Reference comparison", "", "- ADaM reference objects were not found. Comparison skipped.")
  } else {
    lines <- c(
      lines,
      "## Reference comparison",
      "",
      paste0("- Childhood rows match: ", compare_ref$childhood$rows_match),
      paste0("- Childhood cols match: ", compare_ref$childhood$cols_match),
      paste0("- Childhood names match: ", compare_ref$childhood$names_match),
      paste0("- Childhood key hash match: ", compare_ref$childhood$key_md5_match),
      paste0("- Multistudies rows match: ", compare_ref$multistudies$rows_match),
      paste0("- Multistudies cols match: ", compare_ref$multistudies$cols_match),
      paste0("- Multistudies names match: ", compare_ref$multistudies$names_match),
      paste0("- Multistudies key hash match: ", compare_ref$multistudies$key_md5_match)
    )
  }

  if (is.null(compare_self)) {
    lines <- c(lines, "", "## Self reproducibility comparison", "", "- Self-check skipped.")
  } else {
    self_child_hash_match <- identical(compare_self$childhood$key_md5_a, compare_self$childhood$key_md5_b)
    self_multi_hash_match <- identical(compare_self$multistudies$key_md5_a, compare_self$multistudies$key_md5_b)
    lines <- c(
      lines,
      "",
      "## Self reproducibility comparison",
      "",
      paste0("- Childhood rows match (run1 vs run2): ", compare_self$childhood$rows_match),
      paste0("- Childhood cols match (run1 vs run2): ", compare_self$childhood$cols_match),
      paste0("- Childhood names match (run1 vs run2): ", compare_self$childhood$names_match),
      paste0("- Childhood key hash match (run1 vs run2): ", self_child_hash_match),
      paste0("- Multistudies rows match (run1 vs run2): ", compare_self$multistudies$rows_match),
      paste0("- Multistudies cols match (run1 vs run2): ", compare_self$multistudies$cols_match),
      paste0("- Multistudies names match (run1 vs run2): ", compare_self$multistudies$names_match),
      paste0("- Multistudies key hash match (run1 vs run2): ", self_multi_hash_match)
    )
  }

  writeLines(lines, con = report_file)
}

save_adam_outputs <- function(results, adam_dir) {
  dir.create(adam_dir, recursive = TRUE, showWarnings = FALSE)
  dt_childhood_exposure <- results$dt_childhood_exposure
  dt_multistudies <- results$dt_multistudies
  save(dt_childhood_exposure, file = file.path(adam_dir, "dt_childhood_exposure.RData"))
  save(dt_multistudies, file = file.path(adam_dir, "dt_multistudies.RData"))
  if (!is.null(results$dt_dreams_multi_ssq)) {
    dt_dreams_multi_ssq <- results$dt_dreams_multi_ssq
    save(dt_dreams_multi_ssq, file = file.path(adam_dir, "dt_dreams_multi_ssq.RData"))
  }
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  public_root <- find_public_repo_root()
  obj00_root <- file.path(public_root, "scripts", "obj00")
  release_dir <- file.path(public_root, "results")

  private_data_dir_raw <- args[["private-data-dir"]] %||% Sys.getenv("CO_LUMINATE_PRIVATE_DATA_DIR", unset = "")
  if (identical(private_data_dir_raw, "")) {
    stop(
      "Missing private data root. Provide --private-data-dir=/path/to/co-luminate ",
      "or set CO_LUMINATE_PRIVATE_DATA_DIR. ",
      "The script intentionally does not hardcode external/private machine paths."
    )
  }
  private_data_dir <- normalize_path_safe(private_data_dir_raw)
  examples_out <- normalize_path_safe(
    args[["examples-out"]] %||% file.path(public_root, "data", "data_examples")
  )
  validation_out <- normalize_path_safe(
    args[["validation-out"]] %||% file.path(release_dir, "validation")
  )
  adam_out <- normalize_path_safe(
    args[["adam-out"]] %||% file.path(release_dir, "intermediate", "adam")
  )

  score_engines_arg <- args[["score-engines-rds"]] %||%
    "data_management/data_examples/1_staging_snippets/derived_models/irt_joint_models.rds"
  if (!identical(score_engines_arg, "") && grepl("^/", score_engines_arg)) {
    score_engines_arg <- rel_path(obj00_root, score_engines_arg)
  }
  scoring_script_arg <- args[["scoring-script"]] %||% "data_management/scripts/02_harmonization/02b_irt_model.R"
  if (!identical(scoring_script_arg, "") && grepl("^/", scoring_script_arg)) {
    scoring_script_arg <- rel_path(obj00_root, scoring_script_arg)
  }

  cfg <- list(
    include_index_final_steps = as_bool(args[["include-final-steps"]], default = FALSE),
    save_to_private_adam = as_bool(args[["save-adam"]], default = FALSE),
    build_examples = as_bool(args[["build-examples"]], default = TRUE),
    self_repro_check = as_bool(args[["self-repro-check"]], default = TRUE),
    example_rows = as.integer(args[["example-rows"]] %||% 20L)
  )

  if (is.na(cfg$example_rows) || cfg$example_rows < 1L) cfg$example_rows <- 20L

  if (isTRUE(cfg$include_index_final_steps) && !file.exists(file.path(obj00_root, score_engines_arg))) {
    stop(
      "Final-step mode requires a valid scoring engine file. ",
      "Expected default at scripts/obj00/", score_engines_arg,
      " or pass --score-engines-rds=/path/to/irt_joint_models.rds"
    )
  }

  message("== OBJ00 Public Rebuild Pipeline ==")
  message("private_data_dir: ", private_data_dir)
  message("examples_out    : ", examples_out)
  message("validation_out  : ", validation_out)
  message("adam_out        : ", adam_out)

  staged <- stage_private_inputs(
    obj00_root = obj00_root,
    private_data_dir = private_data_dir,
    include_final_steps = cfg$include_index_final_steps
  )
  on.exit(safe_unlink_dir(staged$stage_root), add = TRUE)

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(obj00_root)
  suppressMessages(here::i_am("data_management/scripts/01_cleaning/01b_build_core_datasets_candidate.R"))

  builder_env <- new.env(parent = globalenv())
  source(
    file.path(obj00_root, "data_management", "scripts", "01_cleaning", "01b_build_core_datasets_candidate.R"),
    local = builder_env,
    chdir = TRUE
  )

  if (!exists("build_obj00_core_datasets", envir = builder_env, inherits = FALSE)) {
    stop("build_obj00_core_datasets() not found in candidate builder script.")
  }

  t0 <- Sys.time()
  results <- builder_env$build_obj00_core_datasets(list(
    private_dt_dir = rel_path(obj00_root, staged$stage_root),
    save_outputs = FALSE,
    verbose = TRUE,
    include_index_final_steps = cfg$include_index_final_steps,
    final_inputs_from_adam = cfg$include_index_final_steps,
    save_to_private_adam = FALSE,
    score_engines_rds = score_engines_arg,
    scoring_script = scoring_script_arg
  ))
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
  message("Build completed in ", elapsed, " seconds.")

  ref_compare <- compare_to_reference(results, ref_dir = file.path(staged$stage_root, "adam"))
  self_compare <- NULL

  if (isTRUE(cfg$self_repro_check)) {
    message("Running self reproducibility check (second build)...")
    run2 <- builder_env$build_obj00_core_datasets(list(
      private_dt_dir = rel_path(obj00_root, staged$stage_root),
      save_outputs = FALSE,
      verbose = FALSE,
      include_index_final_steps = cfg$include_index_final_steps,
      final_inputs_from_adam = cfg$include_index_final_steps,
      save_to_private_adam = FALSE,
      score_engines_rds = score_engines_arg,
      scoring_script = scoring_script_arg
    ))
    self_compare <- compare_two_runs(results, run2)
  }

  if (isTRUE(cfg$save_to_private_adam)) {
    save_adam_outputs(results, adam_out)
    message("Saved ADaM outputs to: ", adam_out)
  }

  if (isTRUE(cfg$build_examples)) {
    write_examples(results, out_dir = examples_out, n_rows = cfg$example_rows)
    message("Wrote synthetic examples to: ", examples_out)
  }

  dir.create(validation_out, recursive = TRUE, showWarnings = FALSE)
  report_file <- file.path(validation_out, "obj00_pipeline_validation_report.md")
  write_validation_report(report_file, cfg, results, ref_compare, self_compare)
  message("Wrote validation report: ", report_file)

  message("== Pipeline complete ==")
  if (isTRUE(ref_compare$reference_found)) {
    message(
      "Reference check summary | child(rows,cols,names,key): ",
      ref_compare$childhood$rows_match, ",",
      ref_compare$childhood$cols_match, ",",
      ref_compare$childhood$names_match, ",",
      ref_compare$childhood$key_md5_match
    )
    message(
      "Reference check summary | multi(rows,cols,names,key): ",
      ref_compare$multistudies$rows_match, ",",
      ref_compare$multistudies$cols_match, ",",
      ref_compare$multistudies$names_match, ",",
      ref_compare$multistudies$key_md5_match
    )
  } else {
    message("Reference ADaM files not found; comparison skipped.")
  }

  if (!is.null(self_compare)) {
    message(
      "Self-check summary | child(rows,cols,names,key): ",
      self_compare$childhood$rows_match, ",",
      self_compare$childhood$cols_match, ",",
      self_compare$childhood$names_match, ",",
      identical(self_compare$childhood$key_md5_a, self_compare$childhood$key_md5_b)
    )
    message(
      "Self-check summary | multi(rows,cols,names,key): ",
      self_compare$multistudies$rows_match, ",",
      self_compare$multistudies$cols_match, ",",
      self_compare$multistudies$names_match, ",",
      identical(self_compare$multistudies$key_md5_a, self_compare$multistudies$key_md5_b)
    )
  }
}

main()
