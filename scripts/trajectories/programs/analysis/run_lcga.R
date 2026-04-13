# ==============================================================================
# SCRIPT: Run Harmonized LCGA Models
# PURPOSE: Execute univariate and joint LCGA models from the public trajectory bundle.
# NOTES:
#   - Supports optional CLI filters: --spec=JPA,JPAC --time_label=time_0_5
#   - Supports --dry_run=true to write inputs without running Mplus
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(here)
  library(MplusAutomation)
  library(purrr)
})

find_trajectory_root <- function(start = getwd()) {
  is_trajectory_root <- function(path) {
    dir.exists(file.path(path, "programs", "analysis")) &&
      dir.exists(file.path(path, "programs", "utils"))
  }

  build_candidates <- function(path) {
    if (is.null(path) || !nzchar(path)) {
      return(character())
    }

    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    parents <- character()
    current <- path

    repeat {
      parents <- c(parents, current)
      next_path <- dirname(current)
      if (identical(next_path, current)) {
        break
      }
      current <- next_path
    }

    unique(c(
      parents,
      file.path(parents, "scripts", "trajectories"),
      file.path(parents, "trajectories")
    ))
  }

  candidate_roots <- unique(c(
    Sys.getenv("COLU_TRAJECTORY_ROOT", unset = ""),
    start,
    here::here(),
    file.path(start, "scripts", "trajectories"),
    file.path(here::here(), "scripts", "trajectories")
  ))

  for (candidate in unique(unlist(lapply(candidate_roots, build_candidates)))) {
    if (is_trajectory_root(candidate)) {
      resolved <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
      Sys.setenv(COLU_TRAJECTORY_ROOT = resolved)
      return(resolved)
    }
  }

  stop("Unable to locate the public trajectory bundle root from the current session.")
}

trajectory_root <- find_trajectory_root()
trajectory_path <- function(...) file.path(trajectory_root, ...)

source(trajectory_path("programs", "utils", "analysis_configuration.R"), local = TRUE)
source(trajectory_path("programs", "utils", "data_prep.R"), local = TRUE)
source(trajectory_path("programs", "utils", "lcga_selection.R"), local = TRUE)
source(trajectory_path("programs", "utils", "lcga_reporting_core.R"), local = TRUE)
source(trajectory_path("programs", "utils", "mplus_lcga.R"), local = TRUE)

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (grepl("^--", arg) && grepl("=", arg)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      val <- sub("^--[^=]+=(.*)$", "\\1", arg)
      out[[key]] <- val
    }
  }
  out
}

split_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    return(NULL)
  }
  trimws(unlist(strsplit(x, ",", fixed = TRUE)))
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  tolower(x) %in% c("true", "1", "yes", "y")
}

resolve_min_class_cutoff <- function(k, anchor_threshold = NULL, config = cfg) {
  mode <- tolower(config$min_class_pct_mode %||% "fixed")
  anchor_k <- suppressWarnings(as.integer(config$min_class_pct_anchor_k %||% NA))

  if (identical(mode, "anchor_k2")) {
    mode <- "anchor_k"
    anchor_k <- 2L
  }

  if (identical(mode, "anchor_k")) {
    if (!is.na(anchor_k) && k > anchor_k && !is.null(anchor_threshold)) {
      return(as.numeric(anchor_threshold))
    }
    return(as.numeric(config$min_class_pct))
  }

  as.numeric(config$min_class_pct)
}

list_canonical_lcga_outputs <- function(path = ".") {
  files <- list.files(
    path = path,
    pattern = "^lcga_model_class_[0-9]+\\.out$",
    full.names = TRUE,
    recursive = FALSE
  )
  files[order(as.integer(sub("^.*lcga_model_class_([0-9]+)\\.out$", "\\1", files)))]
}

read_canonical_lcga_models <- function(files) {
  if (length(files) == 0) {
    return(list())
  }

  stats::setNames(
    lapply(files, function(file) suppressWarnings(readModels(target = file, quiet = TRUE))),
    nm = basename(files)
  )
}

derive_structural_sensitivity_dir <- function(output_dir, profile_id) {
  dir_path <- file.path(output_dir, "structural_sensitivity", profile_id)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  dir_path
}

extract_summary_value <- function(model_obj, candidates) {
  if (is.null(model_obj$summaries)) {
    return(NA_real_)
  }

  summaries <- model_obj$summaries
  for (candidate in candidates) {
    if (candidate %in% names(summaries)) {
      value <- suppressWarnings(as.numeric(summaries[[candidate]][1]))
      if (!is.na(value)) {
        return(value)
      }
    }
  }

  NA_real_
}

model_estimated_successfully <- function(model_obj) {
  if (is.null(model_obj)) {
    return(FALSE)
  }

  metric <- extract_summary_value(model_obj, c("LL", "BIC", "AIC", "Observations"))
  !is.na(metric)
}

build_structural_sensitivity_summary_row <- function(
    structure_profile,
    model_obj,
    selected_k,
    directory,
    output_file,
    baseline = FALSE,
    error_message = NULL
) {
  structure_profile <- normalise_lcga_structure_profile(structure_profile)
  estimated <- model_estimated_successfully(model_obj)
  free_variances <- structure_profile$free_variances %||% character(0)
  free_covariances <- structure_profile$free_covariances %||% list()
  status <- if (!is.null(error_message)) {
    "Error"
  } else if (estimated) {
    "Estimated"
  } else if (file.exists(output_file)) {
    "Output available"
  } else {
    "Not run"
  }
  evidence <- if (!is.null(error_message)) {
    error_message
  } else if (estimated) {
    glue(
      "Free variances: {ifelse(length(free_variances) == 0, 'none', paste(free_variances, collapse = ', '))}; ",
      "free covariances: {ifelse(length(free_covariances) == 0, 'none', paste(vapply(free_covariances, paste, collapse = '-', character(1)), collapse = ', '))}."
    )
  } else if (file.exists(output_file)) {
    "Output file exists but summary metrics were not parsed."
  } else {
    "No output file detected."
  }

  tibble(
    ProfileID = structure_profile$id,
    ProfileLabel = structure_profile$label,
    Baseline = baseline,
    Description = structure_profile$description,
    SelectedK = as.integer(selected_k),
    AlternativeWithinClassHeterogeneity = !baseline && length(free_variances) > 0,
    AlternativeCovarianceStructure = !baseline && length(free_covariances) > 0,
    Status = status,
    BIC = extract_summary_value(model_obj, c("BIC")),
    Entropy = extract_summary_value(model_obj, c("Entropy")),
    MinClassPct = suppressWarnings(as.numeric(get_min_class_pct(model_obj))),
    OutputFile = output_file,
    SyntaxDirectory = directory,
    Evidence = evidence
  )
}

run_structural_sensitivity_audit <- function(
    best_selection,
    baseline_model,
    paths,
    lcga_spec,
    lp,
    df_mplus,
    time_scores,
    run_profile,
    dry_run = FALSE,
    config = cfg
) {
  profiles <- config$structural_sensitivity_profiles %||% list()
  profiles <- profiles[lengths(profiles) > 0]

  if (!isTRUE(config$structural_sensitivity_enabled) || length(profiles) == 0 || isTRUE(dry_run)) {
    return(NULL)
  }

  baseline_idx <- which(names(profiles) == "strict_lcga")[1]
  if (is.na(baseline_idx) || length(baseline_idx) == 0) {
    baseline_idx <- 1L
  }
  baseline_profile <- profiles[[baseline_idx]]
  baseline_profile <- normalise_lcga_structure_profile(baseline_profile)
  selected_k <- as.integer(best_selection$Classes)
  baseline_out <- file.path(paths$output_dir, glue("lcga_model_class_{selected_k}.out"))

  summary_rows <- list(
    build_structural_sensitivity_summary_row(
      structure_profile = baseline_profile,
      model_obj = baseline_model,
      selected_k = selected_k,
      directory = paths$output_dir,
      output_file = baseline_out,
      baseline = TRUE
    )
  )

  run_records <- list(
    strict_lcga = list(
      profile = baseline_profile,
      directory = paths$output_dir,
      output_file = baseline_out,
      model = baseline_model,
      baseline = TRUE
    )
  )

  alt_profiles <- profiles[setdiff(names(profiles), baseline_profile$id)]
  if (length(alt_profiles) == 0) {
    return(list(selected_k = selected_k, summary = bind_rows(summary_rows), runs = run_records))
  }

  sensitivity_run_profile <- list(
    tier = run_profile$tier,
    stage = "structural_sensitivity",
    processors = config$structural_sensitivity_processors %||% run_profile$processors,
    starts = if (identical(run_profile$stage %||% "screen", "final")) {
      config$structural_sensitivity_starts_final %||% run_profile$starts
    } else {
      config$structural_sensitivity_starts_screen %||% run_profile$starts
    }
  )

  current_wd <- getwd()
  on.exit(setwd(current_wd), add = TRUE)

  for (profile_id in names(alt_profiles)) {
    structure_profile <- normalise_lcga_structure_profile(alt_profiles[[profile_id]])
    sens_dir <- derive_structural_sensitivity_dir(paths$output_dir, structure_profile$id)
    out_file <- file.path(sens_dir, glue("lcga_model_class_{selected_k}.out"))
    model_obj <- NULL
    error_message <- NULL

    setwd(sens_dir)

    sens_model <- build_model_object(
      k = selected_k,
      lcga_spec = lcga_spec,
      lp = lp,
      df_mplus = df_mplus,
      time_scores = time_scores,
      run_profile = sensitivity_run_profile,
      dry_run = dry_run,
      structure_profile = structure_profile,
      include_selection_tests = FALSE,
      include_plot = FALSE
    )

    if (isTRUE(config$rerun_existing) || !file.exists(out_file)) {
      tryCatch({
        mplusModeler(
          sens_model,
          modelout = glue("lcga_model_class_{selected_k}.inp"),
          dataout = glue("analysis_data_class_{selected_k}.dat"),
          run = 1L,
          writeData = "always",
          hashfilename = FALSE
        )
      }, error = function(e) {
        error_message <<- e$message
      })
    }

    if (is.null(error_message) && file.exists(out_file)) {
      model_obj <- tryCatch(
        suppressWarnings(readModels(target = out_file, quiet = TRUE)),
        error = function(e) {
          error_message <<- e$message
          NULL
        }
      )
    } else if (is.null(error_message) && !file.exists(out_file)) {
      error_message <- "Sensitivity run did not produce an output file."
    }

    summary_rows[[length(summary_rows) + 1]] <- build_structural_sensitivity_summary_row(
      structure_profile = structure_profile,
      model_obj = model_obj,
      selected_k = selected_k,
      directory = sens_dir,
      output_file = out_file,
      baseline = FALSE,
      error_message = error_message
    )

    run_records[[structure_profile$id]] <- list(
      profile = structure_profile,
      directory = sens_dir,
      output_file = out_file,
      model = model_obj,
      baseline = FALSE,
      error_message = error_message
    )
  }

  list(
    selected_k = selected_k,
    summary = bind_rows(summary_rows),
    runs = run_records
  )
}

build_model_object <- function(
    k,
    lcga_spec,
    lp,
    df_mplus,
    time_scores,
    run_profile,
    dry_run = FALSE,
    structure_profile = NULL,
    include_selection_tests = TRUE,
    include_plot = TRUE
) {
  all_item_names <- lp$all_item_names
  series_tokens <- purrr::imap(lp$item_names_list, function(items, pfx) {
    purrr::map_chr(items, function(vn) {
      ag <- suppressWarnings(as.integer(sub(paste0("^", pfx), "", vn)))
      if (is.na(ag)) vn else paste0(vn, "(", ag, ")")
    })
  }) %>%
    unlist(use.names = FALSE)
  series_stmt <- mplus_series_statement(series_tokens)

  usevars_stmt <- mplus_wrap_tokens(c("USUBJID", all_item_names), prefix = "USEVARIABLES =", width = 88)
  cat_stmt <- mplus_wrap_tokens(all_item_names, prefix = "CATEGORICAL =", width = 88)

  class_blocks <- if (length(lcga_spec$outcome_prefixes) > 1) {
    get_joint_class_blocks(lp$item_names_list, time_scores, k, structure_profile = structure_profile)
  } else {
    get_univariate_class_blocks(lp$all_item_names, time_scores, k, structure_profile = structure_profile)
  }

  base_analysis <- glue(
    "TYPE = MIXTURE;\n",
    "ESTIMATOR = MLR;\n",
    "ALGORITHM = INTEGRATION;\n",
    "INTEGRATION = MONTECARLO;\n",
    "PROCESSORS = {run_profile$processors};"
  )

  mplusObject(
    TITLE = glue("LCGA - {k} Classes ({lcga_spec$name})"),
    VARIABLE = glue(
      "IDVARIABLE = USUBJID;\n",
      "CLASSES = c({k});\n",
      "{usevars_stmt}\n",
      "{cat_stmt}\n",
      "MISSING = ALL({cfg$missing_code});\n"
    ),
    ANALYSIS = glue("{base_analysis}\nSTARTS = {run_profile$starts};\n"),
    MODEL = class_blocks,
    OUTPUT = if (isTRUE(include_selection_tests) && k > 1) "TECH11 TECH14 SAMPSTAT;" else "SAMPSTAT;",
    PLOT = if (isTRUE(include_plot)) glue("TYPE = PLOT3;\n{series_stmt}") else "",
    SAVEDATA = if (k > 1 && !isTRUE(dry_run)) glue("FILE = savedata_class_{k}.dat; SAVE = CPROB BCHWEIGHTS;") else "",
    rdata = df_mplus,
    modelout = glue("lcga_model_class_{k}.inp")
  )
}

run_trajectory_lcga <- function(args = commandArgs(trailingOnly = TRUE)) {
  args <- parse_cli_args(args)
  spec_filter <- split_csv_arg(args$spec)
  time_filter <- split_csv_arg(args$time_label)
  dry_run <- parse_bool(args$dry_run, default = FALSE)

  DT <- load_long_data(cfg$rdata_path, cfg$data_object)

  selected_specs <- lcga_specs
  if (!is.null(spec_filter)) {
    selected_specs <- Filter(function(spec) spec$name %in% spec_filter, selected_specs)
  }

  selected_time_ranges <- cfg$time_ranges
  if (!is.null(time_filter)) {
    selected_time_ranges <- Filter(function(age_range) {
      time_label <- glue("time_{min(age_range)}_{max(age_range)}")
      time_label %in% time_filter
    }, selected_time_ranges)
  }

  if (length(selected_specs) == 0) {
    stop("No LCGA specifications matched the requested filter.")
  }

  if (length(selected_time_ranges) == 0) {
    stop("No time ranges matched the requested filter.")
  }

  for (age_range in selected_time_ranges) {
    ages <- as.integer(age_range)
    time_scores <- 0:(length(ages) - 1)
    time_label <- glue("time_{min(ages)}_{max(ages)}")

    for (lcga_spec in selected_specs) {
      is_joint <- length(lcga_spec$outcome_prefixes) > 1
      paths <- derive_lcga_paths(cfg$base_mplus_dir, lcga_spec$name, time_label, is_joint = is_joint)
      run_profile <- get_lcga_execution_profile(lcga_spec, ages, stage = cfg$lcga_fit_stage, config = cfg)
      min_class_anchor_threshold <- NULL
      min_class_anchor_k <- suppressWarnings(as.integer(cfg$min_class_pct_anchor_k %||% NA))
      min_class_mode <- tolower(cfg$min_class_pct_mode %||% "fixed")

      if (identical(min_class_mode, "anchor_k2")) {
        min_class_mode <- "anchor_k"
        min_class_anchor_k <- 2L
      }

      message(glue("\n>>> Running trajectory bundle: {lcga_spec$name} | {time_label}"))
      message(glue("    -> Execution tier: {run_profile$tier} | Stage: {run_profile$stage} | PROCESSORS = {run_profile$processors} | STARTS = {run_profile$starts}"))
      if (isTRUE(dry_run)) {
        message("    -> Dry run enabled: inputs will be written but Mplus will not be executed.")
      }
      if (identical(min_class_mode, "anchor_k") && !is.na(min_class_anchor_k)) {
        message(glue("    -> Min class cutoff mode: anchor at k = {min_class_anchor_k} (current fixed floor = {cfg$min_class_pct})"))
      }

      lp <- prep_wide_lcga(DT, cfg$id_var, cfg$time_var, ages, lcga_spec, vars_registry)
      df_lcga <- lp$df_lcga %>%
        mutate(USUBJID_ORIG = as.character(USUBJID))

      observation_filter <- apply_lcga_observation_filter(
        df_lcga = df_lcga,
        item_names_list = lp$item_names_list,
        min_observed_occasions = cfg$min_observed_occasions %||% NULL,
        mode = cfg$observed_occasions_mode %||% "any"
      )
      df_lcga <- observation_filter$df_lcga

      if (isTRUE(observation_filter$summary$Applied[[1]])) {
        message(glue(
          "    -> Observation filter applied: >= {observation_filter$summary$MinObservedOccasions[[1]]} occasions ",
          "({observation_filter$summary$ObservedOccasionsMode[[1]]}); retained ",
          "{observation_filter$summary$NRetained[[1]]} / {observation_filter$summary$NInput[[1]]} ",
          "({sprintf('%.1f', 100 * observation_filter$summary$PctRetained[[1]])}%)"
        ))
      }

      suppressWarnings({
        df_lcga$USUBJID <- as.numeric(df_lcga$USUBJID_ORIG)
      })
      if (any(is.na(df_lcga$USUBJID))) {
        df_lcga <- df_lcga %>% mutate(USUBJID = dplyr::dense_rank(USUBJID_ORIG))
      }

      all_item_names <- lp$all_item_names
      df_mplus <- df_lcga %>%
        mutate(across(all_of(all_item_names), ~ if_else(is.na(.x), cfg$missing_code, as.numeric(.x)))) %>%
        select(USUBJID, all_of(all_item_names))

      current_wd <- getwd()
      on.exit(setwd(current_wd), add = TRUE)
      setwd(paths$output_dir)

      for (k in cfg$k_classes) {
        model_obj <- build_model_object(
          k = k,
          lcga_spec = lcga_spec,
          lp = lp,
          df_mplus = df_mplus,
          time_scores = time_scores,
          run_profile = run_profile,
          dry_run = dry_run
        )

        out_file <- glue("lcga_model_class_{k}.out")

        if (isTRUE(cfg$rerun_existing) || !file.exists(out_file)) {
          tryCatch({
            mplusModeler(
              model_obj,
              modelout = glue("lcga_model_class_{k}.inp"),
              dataout = glue("analysis_data_class_{k}.dat"),
              run = if (isTRUE(dry_run)) 0L else 1L,
              writeData = "always",
              hashfilename = FALSE
            )
          }, error = function(e) {
            message(glue("Error K={k}: {e$message}"))
          })
        }

        if (!isTRUE(dry_run) && k > 1 && file.exists(out_file)) {
          res_k <- suppressWarnings(readModels(target = out_file, quiet = TRUE))
          if (!is.null(res_k$class_counts$modelEstimated)) {
            min_class_observed <- min(res_k$class_counts$modelEstimated$proportion, na.rm = TRUE)
            if (identical(min_class_mode, "anchor_k") && !is.na(min_class_anchor_k) && k == min_class_anchor_k) {
              min_class_anchor_threshold <- min_class_observed
              message(glue("    -> Anchored min class cutoff at k = {k}: {sprintf('%.3f', min_class_anchor_threshold)}"))
            }

            effective_min_class_cutoff <- resolve_min_class_cutoff(
              k = k,
              anchor_threshold = min_class_anchor_threshold,
              config = cfg
            )

            if (min_class_observed < effective_min_class_cutoff) {
              message(glue(
                "    -> Stopping after k = {k}: observed min class {sprintf('%.3f', min_class_observed)} ",
                "fell below cutoff {sprintf('%.3f', effective_min_class_cutoff)}"
              ))
              break
            }
          }
        }
      }

      canonical_out_files <- list_canonical_lcga_outputs(".")
      all_results <- read_canonical_lcga_models(canonical_out_files)
      fit_summary_aug <- tibble()
      structural_sensitivity <- NULL

      try({
        raw_summary <- purrr::imap_dfr(
          all_results,
          ~ if (!is.null(.x$summaries)) .x$summaries %>% mutate(model_id = .y) else NULL
        )
        if (nrow(raw_summary) > 0) {
          fit_summary_aug <- augment_fit_with_class_sizes(raw_summary, all_results)
        }
      }, silent = TRUE)

      fit_summary_norm <- normalise_fit_summary(fit_summary_aug)
      if (!isTRUE(dry_run) && nrow(fit_summary_norm) > 0) {
        best_selection <- tryCatch(
          select_best_solution(fit_summary_norm, spec_name = lcga_spec$name, time_label = time_label, config = cfg),
          error = function(e) {
            message(glue("    ! Structural sensitivity selection skipped: {e$message}"))
            NULL
          }
        )

        if (!is.null(best_selection) && !is.na(best_selection$model_id)) {
          baseline_model <- all_results[[best_selection$model_id]]
          structural_sensitivity <- tryCatch(
            run_structural_sensitivity_audit(
              best_selection = best_selection,
              baseline_model = baseline_model,
              paths = paths,
              lcga_spec = lcga_spec,
              lp = lp,
              df_mplus = df_mplus,
              time_scores = time_scores,
              run_profile = run_profile,
              dry_run = dry_run,
              config = cfg
            ),
            error = function(e) {
              message(glue("    ! Structural sensitivity audit error: {e$message}"))
              NULL
            }
          )

          if (!is.null(structural_sensitivity$summary) && nrow(structural_sensitivity$summary) > 0) {
            utils::write.csv(
              structural_sensitivity$summary,
              file = file.path(paths$listings_dir, glue("structural_sensitivity_{lcga_spec$name}_{time_label}.csv")),
              row.names = FALSE
            )
          }
        }
      }

      setwd(current_wd)

      saveRDS(
        list(
          meta = paths,
          run_profile = run_profile,
          all_results = all_results,
          fit_summary_aug = fit_summary_aug,
          structural_sensitivity = structural_sensitivity,
          observation_filter = observation_filter$summary,
          df_mplus = df_mplus,
          item_names_list = lp$item_names_list,
          meta_long = lp$outcome_meta,
          spec = lcga_spec,
          time_label = time_label,
          dry_run = dry_run
        ),
        paths$results_rds_path
      )
    }
  }

  message("\n>>> Trajectory LCGA execution complete.")
}

if (sys.nframe() == 0) {
  run_trajectory_lcga()
}
