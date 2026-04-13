# ==============================================================================
# UTILITY: Data Loading & Directory Management
# PURPOSE: Load longitudinal data, coerce outcome types for Mplus, and reshape
#          to wide LCGA-ready inputs in the public trajectory bundle.
# ==============================================================================

library(dplyr)
library(glue)
library(here)
library(purrr)
library(tidyr)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

if (!exists("find_trajectory_root", mode = "function", inherits = TRUE)) {
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
}

if (!exists("trajectory_path", mode = "function", inherits = TRUE)) {
  trajectory_path <- function(...) {
    file.path(find_trajectory_root(), ...)
  }
}

load_long_data <- function(rdata_path, object_name) {
  if (!file.exists(rdata_path)) {
    stop(glue("Error: Data file not found at {rdata_path}"))
  }

  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)

  if (!exists(object_name, envir = e)) {
    stop(glue("Error: Object '{object_name}' not found in RData."))
  }

  get(object_name, envir = e)
}

coerce_outcome <- function(x, type = c("binary", "ordinal", "continuous"), levels = NULL) {
  type <- match.arg(type)

  if (type == "continuous") {
    return(as.numeric(x))
  }

  if (!is.null(levels)) {
    x <- factor(x, levels = levels, ordered = identical(type, "ordinal"))
  } else {
    x <- if (identical(type, "ordinal")) {
      if (is.ordered(x)) x else ordered(x)
    } else {
      factor(x)
    }
  }

  as.integer(x) - 1L
}

derive_lcga_paths <- function(base_mplus_dir, analysis_name, time_label, is_joint = FALSE) {
  variant <- cfg$artifact_variant %||% NULL
  variant <- if (!is.null(variant) && nzchar(variant)) {
    gsub("[^A-Za-z0-9_-]+", "_", variant)
  } else {
    NULL
  }

  mplus_base <- if (!is.null(variant)) {
    file.path(base_mplus_dir, "variants", variant)
  } else {
    base_mplus_dir
  }

  analysis_root_mplus <- if (is_joint) {
    file.path(mplus_base, "lcga", "joint", analysis_name, time_label)
  } else {
    file.path(mplus_base, "lcga", "univariate", analysis_name, time_label)
  }

  output_base <- trajectory_path("output")
  if (!is.null(variant)) {
    output_base <- file.path(output_base, "variants", variant)
  }
  sub_path <- if (is_joint) {
    file.path("lcga", "joint", analysis_name)
  } else {
    file.path("lcga", "univariate", analysis_name)
  }

  paths <- list(
    mplus_exec = analysis_root_mplus,
    output_dir = file.path(analysis_root_mplus, "rautomation"),
    manual_dir = file.path(analysis_root_mplus, "manual"),
    objects_dir = file.path(output_base, "objects", sub_path),
    figures_dir = file.path(output_base, "figures", sub_path, time_label),
    tables_dir = file.path(output_base, "tables", sub_path, time_label),
    listings_dir = file.path(output_base, "listings", sub_path, time_label)
  )

  purrr::walk(paths, ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))

  paths$results_rds_path <- file.path(paths$objects_dir, glue("lcga_results_{time_label}.rds"))
  paths
}

count_lcga_observed_occasions <- function(df_lcga, item_names_list, mode = c("any", "complete")) {
  mode <- match.arg(mode)

  if (length(item_names_list) == 1) {
    items <- item_names_list[[1]]
    return(rowSums(!is.na(df_lcga[, items, drop = FALSE])))
  }

  ages <- sort(unique(as.integer(gsub(".*?(\\d+)$", "\\1", unlist(item_names_list, use.names = FALSE)))))
  by_age <- sapply(ages, function(age_value) {
    cols <- unlist(lapply(item_names_list, function(items) {
      grep(paste0(age_value, "$"), items, value = TRUE)
    }))

    non_missing_counts <- rowSums(!is.na(df_lcga[, cols, drop = FALSE]))
    if (identical(mode, "complete")) {
      non_missing_counts == length(cols)
    } else {
      non_missing_counts > 0
    }
  })

  rowSums(by_age)
}

apply_lcga_observation_filter <- function(
    df_lcga,
    item_names_list,
    min_observed_occasions = cfg$min_observed_occasions %||% NULL,
    mode = cfg$observed_occasions_mode %||% "any"
) {
  if (is.null(min_observed_occasions) || length(min_observed_occasions) == 0) {
    min_observed_occasions <- NA_integer_
  } else {
    min_observed_occasions <- suppressWarnings(as.integer(min_observed_occasions)[1])
  }

  if (is.na(min_observed_occasions) || min_observed_occasions <= 0) {
    return(list(
      df_lcga = df_lcga,
      summary = tibble(
        Applied = FALSE,
        MinObservedOccasions = NA_integer_,
        ObservedOccasionsMode = mode,
        NInput = nrow(df_lcga),
        NRetained = nrow(df_lcga),
        NRemoved = 0L,
        PctRetained = if (nrow(df_lcga) > 0) 1 else NA_real_
      ),
      observed_occasions = rep(NA_integer_, nrow(df_lcga))
    ))
  }

  observed_occasions <- count_lcga_observed_occasions(df_lcga, item_names_list, mode = mode)
  keep <- observed_occasions >= min_observed_occasions
  filtered <- df_lcga[keep, , drop = FALSE]

  list(
    df_lcga = filtered,
    summary = tibble(
      Applied = TRUE,
      MinObservedOccasions = min_observed_occasions,
      ObservedOccasionsMode = mode,
      NInput = nrow(df_lcga),
      NRetained = nrow(filtered),
      NRemoved = nrow(df_lcga) - nrow(filtered),
      PctRetained = if (nrow(df_lcga) > 0) nrow(filtered) / nrow(df_lcga) else NA_real_
    ),
    observed_occasions = observed_occasions
  )
}

prep_wide_one <- function(dt_long, id_var, time_var, ages, var_spec) {
  dt_sub <- dt_long %>% filter(.data[[time_var]] %in% ages)
  var_levels <- var_spec$levels %||% NULL

  long_std <- dt_sub %>%
    transmute(
      USUBJID = .data[[id_var]],
      AGE = .data[[time_var]],
      RAW = .data[[var_spec$variable]],
      Y = coerce_outcome(RAW, type = var_spec$type, levels = var_levels)
    )

  df_wide <- long_std %>%
    select(USUBJID, AGE, Y) %>%
    pivot_wider(names_from = AGE, values_from = Y, names_prefix = var_spec$wide_prefix) %>%
    select(any_of(c("USUBJID", paste0(var_spec$wide_prefix, ages))))

  item_names <- names(df_wide)[grepl(paste0("^", var_spec$wide_prefix), names(df_wide))]

  list(
    df_wide = df_wide,
    long_std = long_std,
    item_names = item_names,
    transform_levels = var_levels
  )
}

prep_wide_lcga <- function(dt_long, id_var, time_var, ages, lcga_spec, registry) {
  specs <- purrr::map(lcga_spec$outcome_prefixes, get_var_spec_by_prefix, registry = registry)
  preps <- purrr::map(specs, ~ prep_wide_one(dt_long, id_var, time_var, ages, .x))

  df_lcga <- purrr::reduce(purrr::map(preps, "df_wide"), full_join, by = "USUBJID")
  item_names_list <- purrr::map(preps, "item_names")
  names(item_names_list) <- purrr::map_chr(specs, ~ .x$wide_prefix)

  outcome_meta <- purrr::map2(preps, specs, function(p, s) {
    list(
      wide_prefix = s$wide_prefix,
      variable = s$variable,
      label = s$label,
      type = s$type,
      plot_target = s$plot_target %||% "lowest",
      transform_levels = p$transform_levels %||% s$levels %||% NULL
    )
  })

  list(
    df_lcga = df_lcga,
    preps = preps,
    specs = specs,
    item_names_list = item_names_list,
    all_item_names = unlist(item_names_list, use.names = FALSE),
    outcome_meta = outcome_meta
  )
}
