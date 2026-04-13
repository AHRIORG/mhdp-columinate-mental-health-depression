# ==============================================================================
# UTILITY: Harmonized LCGA Configuration
# PURPOSE: Define the canonical LCGA configuration, variable registry, and
#          model specifications for the public trajectory bundle.
# ==============================================================================

library(glue)
library(here)

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

if (!exists("repo_path", mode = "function", inherits = TRUE)) {
  repo_path <- function(...) {
    repo_root <- normalizePath(file.path(find_trajectory_root(), "..", ".."), winslash = "/", mustWork = TRUE)
    file.path(repo_root, ...)
  }
}

resolve_input_file <- function(env_var, default_rel_path) {
  candidates <- c(
    Sys.getenv(env_var, unset = ""),
    repo_path(default_rel_path)
  )

  candidates <- unique(candidates[nzchar(candidates)])

  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  candidates[[length(candidates)]]
}

existing_cfg <- if (exists("cfg", inherits = FALSE) && is.list(cfg)) cfg else NULL
existing_vars_registry <- if (exists("vars_registry", inherits = FALSE) && is.list(vars_registry)) vars_registry else NULL
existing_lcga_specs <- if (exists("lcga_specs", inherits = FALSE) && is.list(lcga_specs)) lcga_specs else NULL

childhood_rdata_path <- resolve_input_file("COLU_CHILDHOOD_RDATA", file.path("data", "inputs", "dt_childhood_exposure.RData"))
multistudies_rdata_path <- resolve_input_file("COLU_MULTISTUDIES_RDATA", file.path("data", "inputs", "dt_multistudies.RData"))

root_dir <- trajectory_path()
dirs_to_create <- c(
  trajectory_path(".cache", "R"),
  trajectory_path("data", "metadata"),
  trajectory_path("models"),
  trajectory_path("output", "tables"),
  trajectory_path("output", "figures"),
  trajectory_path("output", "objects"),
  trajectory_path("output", "listings"),
  trajectory_path("programs", "analysis"),
  trajectory_path("programs", "utils")
)
invisible(lapply(dirs_to_create, dir.create, recursive = TRUE, showWarnings = FALSE))
Sys.setenv(R_USER_CACHE_DIR = trajectory_path(".cache", "R"))

cfg <- list(
  rdata_path = childhood_rdata_path,
  data_object = "dt_childhood_exposure",
  multistudies_path = multistudies_rdata_path,
  multistudies_object = "dt_multistudies",
  id_var = "USUBJID",
  time_var = "EXPAGE",
  time_ranges = list(early = 0:5, mid = 6:12),
  k_classes = 1:15,
  base_mplus_dir = trajectory_path("models"),
  processors = 8,
  starts = "5000 500",
  rerun_existing = FALSE,
  missing_code = -9999,
  artifact_variant = NULL,
  min_observed_occasions = NULL,
  observed_occasions_mode = "any",
  selection_preset = "default",
  selection_test = "BLRT",
  selection_alpha = 0.05,
  min_class_pct = 0.05,
  min_class_pct_mode = "fixed",
  min_class_pct_anchor_k = NULL,
  bic_close_delta = 2,
  use_manual_overrides = FALSE,
  analysis_stage = "draft",
  lcga_fit_stage = NULL,
  lcga_processors_override = NULL,
  lcga_starts_override = NULL,
  lcga_execution_profiles = list(
    heavy_joint = list(
      processors = "10 2",
      starts_screen = "2000 200",
      starts_final = "5000 500"
    ),
    medium_joint = list(
      processors = "8",
      starts_screen = "1500 150",
      starts_final = "3000 300"
    ),
    light_univariate = list(
      processors = "8",
      starts_screen = "1000 100",
      starts_final = "2000 200"
    )
  ),
  manual_class_overrides = list(),
  composite_method = "mean",
  baseline_confounders = character(0),
  mediator_outcome_confounders = character(0),
  class_registry_path = trajectory_path("data", "metadata", "class_registry.csv"),
  structural_sensitivity_enabled = TRUE,
  structural_sensitivity_run_mode = "selected_k",
  structural_sensitivity_processors = "2",
  structural_sensitivity_starts_screen = "10 2",
  structural_sensitivity_starts_final = "20 5",
  structural_sensitivity_profiles = list(
    list(
      id = "strict_lcga",
      label = "Strict LCGA",
      description = "All within-class growth-factor variances and covariances fixed to zero.",
      free_variances = character(0),
      free_covariances = list()
    ),
    list(
      id = "free_is_variances",
      label = "Free I/S Variances",
      description = "Intercept and linear slope variances freed within class; covariance structure remains fixed.",
      free_variances = c("i", "s"),
      free_covariances = list(),
      variance_start_values = c(i = 0.25, s = 0.05)
    ),
    list(
      id = "free_is_covariance",
      label = "Free I/S Variances + Covariance",
      description = "Intercept and linear slope variances plus the within-outcome intercept-slope covariance are freed within class.",
      free_variances = c("i", "s"),
      free_covariances = list(c("i", "s")),
      variance_start_values = c(i = 0.25, s = 0.05),
      covariance_start_values = c(i_s = 0.00)
    )
  )
)

if (!is.null(existing_cfg)) {
  cfg <- utils::modifyList(cfg, existing_cfg, keep.null = TRUE)
}

if (is.null(cfg$lcga_fit_stage) || identical(cfg$lcga_fit_stage, "auto")) {
  cfg$lcga_fit_stage <- if (identical(cfg$analysis_stage, "final")) "final" else "screen"
}

normalise_structural_sensitivity_profiles <- function(profiles) {
  profiles <- profiles %||% list()

  if (length(profiles) == 0) {
    return(list())
  }

  normalised <- lapply(seq_along(profiles), function(idx) {
    profile <- profiles[[idx]]

    if (is.null(profile$id) || !nzchar(profile$id)) {
      profile$id <- paste0("profile_", idx)
    }
    if (is.null(profile$label) || !nzchar(profile$label)) {
      profile$label <- profile$id
    }
    if (is.null(profile$description) || !nzchar(profile$description)) {
      profile$description <- profile$label
    }

    profile$free_variances <- unique(tolower(as.character(profile$free_variances %||% character(0))))

    raw_covs <- profile$free_covariances %||% list()
    if (is.atomic(raw_covs) && !is.list(raw_covs)) {
      raw_covs <- list(raw_covs)
    }
    profile$free_covariances <- lapply(raw_covs, function(pair) {
      sort(tolower(as.character(pair)))
    })

    profile$variance_start_values <- profile$variance_start_values %||% c(i = 0.25, s = 0.05, q = 0.01)
    profile$covariance_start_values <- profile$covariance_start_values %||% c(i_s = 0.00, i_q = 0.00, q_s = 0.00)
    profile
  })

  names(normalised) <- vapply(normalised, `[[`, character(1), "id")
  normalised
}

cfg$structural_sensitivity_profiles <- normalise_structural_sensitivity_profiles(cfg$structural_sensitivity_profiles)

vars_registry <- existing_vars_registry %||% list(
  list(
    variable = "LIVEDWOM",
    wide_prefix = "LWOM",
    label = "Lived Without Mother",
    type = "binary",
    levels = c("No", "Yes"),
    plot_target = "highest"
  ),
  list(
    variable = "LIVEDWOF",
    wide_prefix = "LWOF",
    label = "Lived Without Father",
    type = "binary",
    levels = c("No", "Yes"),
    plot_target = "highest"
  ),
  list(
    variable = "HHCHLD14O",
    wide_prefix = "HC14",
    label = "Overcrowded (Children > 3)",
    type = "binary",
    levels = c("No", "Yes"),
    plot_target = "highest"
  ),
  list(
    variable = "HHSES",
    wide_prefix = "HSES",
    label = "Household Socioeconomic Status",
    type = "ordinal",
    levels = c("Low", "Middle", "High"),
    plot_target = "lowest"
  )
)

lcga_specs <- existing_lcga_specs %||% list(
  list(name = "ULWOM", label = "Univariate Lived Without Mother", outcome_prefixes = c("LWOM")),
  list(name = "ULWOF", label = "Univariate Lived Without Father", outcome_prefixes = c("LWOF")),
  list(name = "UHC14", label = "Univariate Overcrowded (Children > 3)", outcome_prefixes = c("HC14")),
  list(name = "UHSES", label = "Univariate Household Socioeconomic Status", outcome_prefixes = c("HSES")),
  list(name = "JPA", label = "Joint Parental Residence Patterns", outcome_prefixes = c("LWOM", "LWOF")),
  list(name = "JMAES", label = "Joint Maternal Absence + SES", outcome_prefixes = c("LWOM", "HSES")),
  list(name = "JPAC", label = "Joint Parental Absence + Overcrowding", outcome_prefixes = c("LWOM", "LWOF", "HC14")),
  list(name = "JPAES", label = "Joint Parental Absence + SES", outcome_prefixes = c("LWOM", "LWOF", "HSES")),
  list(name = "JPACES", label = "Joint Parental Absence + Overcrowding + SES", outcome_prefixes = c("LWOM", "LWOF", "HC14", "HSES"))
)

label_map <- setNames(
  vapply(vars_registry, function(x) x$label, character(1)),
  vapply(vars_registry, function(x) x$wide_prefix, character(1))
)

spec_map <- setNames(
  vars_registry,
  vapply(vars_registry, function(x) x$wide_prefix, character(1))
)

get_var_spec_by_prefix <- function(prefix, registry = vars_registry) {
  hits <- Filter(function(spec) identical(spec$wide_prefix, prefix), registry)
  if (length(hits) != 1) {
    stop(glue("Variable spec '{prefix}' was not found in analysis_configuration.R"))
  }
  hits[[1]]
}

get_lcga_spec <- function(name, specs = lcga_specs) {
  hits <- Filter(function(spec) identical(spec$name, name), specs)
  if (length(hits) != 1) {
    stop(glue("LCGA spec '{name}' was not found in analysis_configuration.R"))
  }
  hits[[1]]
}

is_joint_lcga_spec <- function(name, specs = lcga_specs) {
  length(get_lcga_spec(name, specs = specs)$outcome_prefixes) > 1
}

classify_lcga_execution_tier <- function(spec, ages) {
  if (length(spec$outcome_prefixes) <= 1) {
    return("light_univariate")
  }

  complexity_score <- length(spec$outcome_prefixes) * length(ages)
  if (complexity_score >= 24) {
    return("heavy_joint")
  }

  "medium_joint"
}

get_lcga_execution_profile <- function(spec, ages, stage = cfg$lcga_fit_stage, config = cfg) {
  stage <- stage %||% "screen"
  tier <- classify_lcga_execution_tier(spec, ages)
  profile <- config$lcga_execution_profiles[[tier]]

  if (is.null(profile)) {
    stop(glue("LCGA execution profile '{tier}' is not defined in analysis_configuration.R"))
  }

  starts_key <- switch(
    stage,
    screen = "starts_screen",
    final = "starts_final",
    stop(glue("Unknown lcga fit stage '{stage}'. Use 'screen' or 'final'."))
  )

  list(
    tier = tier,
    stage = stage,
    processors = config$lcga_processors_override %||% profile$processors %||% config$processors,
    starts = config$lcga_starts_override %||% profile[[starts_key]] %||% config$starts
  )
}
