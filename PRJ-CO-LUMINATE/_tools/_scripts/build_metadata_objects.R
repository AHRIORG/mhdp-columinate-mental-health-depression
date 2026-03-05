#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(utils)
})

project_root <- "."
objects_dir <- file.path(project_root, "_tools", "_objects")
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)

obj00_root <- file.path("..", "OBJ00-Datasets Preperation")
private_adam_roots <- c(
  file.path("..", "..", "_private_use", "data_management", "data", "co-luminate", "adam"),
  file.path("..", "_private_use", "data_management", "data", "co-luminate", "adam"),
  file.path("..", "..", "..", "_private_use", "data_management", "data", "co-luminate", "adam")
)
child_private_paths <- file.path(private_adam_roots, "dt_childhood_exposure.RData")
multi_private_paths <- file.path(private_adam_roots, "dt_multistudies.RData")

child_example_path <- file.path(
  obj00_root,
  "data_management",
  "data_examples",
  "0_raw_snippets",
  "dt_childhood_exposure_head.RData"
)
multi_example_rdata_path <- file.path(
  obj00_root,
  "data_management",
  "data_examples",
  "2_final_snippets",
  "dt_Multistudies_head.RData"
)
multi_example_csv_path <- file.path(
  obj00_root,
  "data_management",
  "data_examples",
  "2_final_snippets",
  "dt_Multistudies_head.csv"
)
master_dict_path <- file.path(
  obj00_root,
  "data_management",
  "docs",
  "dictionaries",
  "dt_master_dict.RData"
)
meta_helper_path <- file.path(
  obj00_root,
  "data_management",
  "scripts",
  "99_utils",
  "metadata_helper.R"
)

load_rdata_df <- function(path, preferred_name = NULL) {
  if (!file.exists(path)) {
    return(NULL)
  }

  src <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = src)

  if (!is.null(preferred_name) && preferred_name %in% loaded_names) {
    preferred_obj <- get(preferred_name, envir = src)
    if (is.data.frame(preferred_obj)) {
      return(preferred_obj)
    }
  }

  for (nm in loaded_names) {
    obj <- get(nm, envir = src)
    if (is.data.frame(obj)) {
      return(obj)
    }
  }

  NULL
}

load_first_df <- function(paths, preferred_name = NULL) {
  unique_paths <- unique(paths)
  for (path in unique_paths) {
    candidate <- tryCatch(
      load_rdata_df(path, preferred_name = preferred_name),
      error = function(e) NULL
    )
    if (!is.null(candidate)) {
      return(list(
        data = candidate,
        source = normalizePath(path, mustWork = FALSE)
      ))
    }
  }
  list(data = NULL, source = NA_character_)
}

child_loaded <- load_first_df(
  c(child_private_paths, child_example_path),
  preferred_name = "dt_childhood_exposure"
)
if (is.null(child_loaded$data)) {
  stop("Unable to load childhood exposure data from private or example paths.")
}
dt_chilhood_exposure <- child_loaded$data
data_source_child <- child_loaded$source

multi_loaded <- load_first_df(
  c(multi_private_paths, multi_example_rdata_path),
  preferred_name = "dt_multistudies"
)
if (is.null(multi_loaded$data)) {
  if (!file.exists(multi_example_csv_path)) {
    stop("Unable to load multistudies data from private paths or local examples.")
  }
  multi_loaded <- list(
    data = read.csv(multi_example_csv_path, check.names = FALSE),
    source = normalizePath(multi_example_csv_path, mustWork = FALSE)
  )
}
dt_multistudies <- multi_loaded$data
data_source_multi <- multi_loaded$source

find_col_name <- function(df, target_name) {
  idx <- match(toupper(target_name), toupper(names(df)))
  if (is.na(idx)) {
    return(NA_character_)
  }
  names(df)[idx]
}

count_unique_ids <- function(df, id_col) {
  if (is.na(id_col) || !(id_col %in% names(df))) {
    return(NA_integer_)
  }
  length(unique(df[[id_col]][!is.na(df[[id_col]])]))
}

derive_clyr_from_visitdt <- function(df) {
  visit_idx <- match("VISITDT", toupper(names(df)))

  if (is.na(visit_idx)) {
    if (!"CLYR" %in% names(df)) {
      df$CLYR <- NA_real_
    }
    return(df)
  }

  visit_values <- df[[visit_idx]]
  if (!inherits(visit_values, "Date")) {
    visit_values <- suppressWarnings(as.Date(visit_values))
  }

  visit_year <- suppressWarnings(as.numeric(format(visit_values, "%Y")))
  visit_year_min <- suppressWarnings(min(visit_year, na.rm = TRUE))

  if (is.finite(visit_year_min)) {
    df$CLYR <- visit_year - visit_year_min
  } else {
    df$CLYR <- NA_real_
  }

  df
}

dt_multistudies <- derive_clyr_from_visitdt(dt_multistudies)

id_col_child <- find_col_name(dt_chilhood_exposure, "USUBJID")
id_col_multi <- find_col_name(dt_multistudies, "USUBJID")
dpbn_col_multi <- find_col_name(dt_multistudies, "DPBN")

if (is.na(id_col_child) || is.na(id_col_multi)) {
  stop("Unable to align datasets: USUBJID was not found in both source datasets.")
}
if (is.na(dpbn_col_multi)) {
  stop("Unable to build metadata CONSORT: DPBN was not found in dt_multistudies.")
}

dt_multistudies_dpbn <- dt_multistudies[!is.na(dt_multistudies[[dpbn_col_multi]]), , drop = FALSE]
overlap_ids <- intersect(
  unique(dt_multistudies_dpbn[[id_col_multi]]),
  unique(dt_chilhood_exposure[[id_col_child]])
)
overlap_ids <- overlap_ids[!is.na(overlap_ids)]

dt_multistudies <- dt_multistudies_dpbn[
  dt_multistudies_dpbn[[id_col_multi]] %in% overlap_ids,
  ,
  drop = FALSE
]
dt_chilhood_exposure <- dt_chilhood_exposure[
  dt_chilhood_exposure[[id_col_child]] %in% overlap_ids,
  ,
  drop = FALSE
]

consort_counts <- data.frame(
  metric = c(
    "multi_records_loaded",
    "multi_participants_loaded",
    "multi_records_with_dpbn",
    "multi_participants_with_dpbn",
    "child_records_loaded",
    "child_participants_loaded",
    "overlap_participants",
    "multi_records_final",
    "multi_participants_final",
    "child_records_final",
    "child_participants_final",
    "multi_missing_dpbn_final"
  ),
  value = c(
    nrow(multi_loaded$data),
    count_unique_ids(multi_loaded$data, id_col_multi),
    nrow(dt_multistudies_dpbn),
    count_unique_ids(dt_multistudies_dpbn, id_col_multi),
    nrow(child_loaded$data),
    count_unique_ids(child_loaded$data, id_col_child),
    length(overlap_ids),
    nrow(dt_multistudies),
    count_unique_ids(dt_multistudies, id_col_multi),
    nrow(dt_chilhood_exposure),
    count_unique_ids(dt_chilhood_exposure, id_col_child),
    sum(is.na(dt_multistudies[[dpbn_col_multi]]))
  ),
  stringsAsFactors = FALSE
)

id_col_multi_loaded <- find_col_name(multi_loaded$data, "USUBJID")
study_col_multi_loaded <- find_col_name(multi_loaded$data, "STUDY")
dpbn_col_multi_loaded <- find_col_name(multi_loaded$data, "DPBN")

if (is.na(id_col_multi_loaded)) {
  id_col_multi_loaded <- id_col_multi
}
if (is.na(dpbn_col_multi_loaded)) {
  dpbn_col_multi_loaded <- dpbn_col_multi
}

ids_multi_loaded <- unique(as.character(multi_loaded$data[[id_col_multi_loaded]]))
ids_multi_loaded <- ids_multi_loaded[!is.na(ids_multi_loaded) & nzchar(ids_multi_loaded)]

ids_multi_with_dpbn <- unique(as.character(
  multi_loaded$data[[id_col_multi_loaded]][
    !is.na(multi_loaded$data[[id_col_multi_loaded]]) &
      !is.na(multi_loaded$data[[dpbn_col_multi_loaded]])
  ]
))
ids_multi_with_dpbn <- ids_multi_with_dpbn[!is.na(ids_multi_with_dpbn) & nzchar(ids_multi_with_dpbn)]

ids_child_loaded <- unique(as.character(child_loaded$data[[id_col_child]]))
ids_child_loaded <- ids_child_loaded[!is.na(ids_child_loaded) & nzchar(ids_child_loaded)]

ids_excl_missing_endpoint <- setdiff(ids_multi_loaded, ids_multi_with_dpbn)
ids_excl_no_sdoh_match <- setdiff(ids_multi_with_dpbn, ids_child_loaded)
ids_excl_no_eligible_mh <- setdiff(ids_child_loaded, ids_multi_with_dpbn)

get_study_list <- function(ids) {
  if (is.na(study_col_multi_loaded) || length(ids) == 0) {
    return(character(0))
  }

  study_vals <- as.character(multi_loaded$data[[study_col_multi_loaded]])
  id_vals <- as.character(multi_loaded$data[[id_col_multi_loaded]])
  keep <- !is.na(study_vals) & nzchar(study_vals) & !is.na(id_vals) & id_vals %in% ids
  sort(unique(study_vals[keep]))
}

get_study_breakdown <- function(ids) {
  if (is.na(study_col_multi_loaded) || length(ids) == 0) {
    return(data.frame(study = character(), n = integer(), stringsAsFactors = FALSE))
  }

  id_vals <- as.character(multi_loaded$data[[id_col_multi_loaded]])
  study_vals <- as.character(multi_loaded$data[[study_col_multi_loaded]])

  keep <- !is.na(id_vals) & nzchar(id_vals) & !is.na(study_vals) & nzchar(study_vals) & id_vals %in% ids
  if (!any(keep)) {
    return(data.frame(study = character(), n = integer(), stringsAsFactors = FALSE))
  }

  pairs <- unique(data.frame(
    id = id_vals[keep],
    study = study_vals[keep],
    stringsAsFactors = FALSE
  ))

  counts <- as.data.frame(table(pairs$study), stringsAsFactors = FALSE)
  names(counts) <- c("study", "n")
  counts$n <- as.integer(counts$n)
  counts <- counts[order(-counts$n, counts$study), , drop = FALSE]
  rownames(counts) <- NULL
  counts
}

collapse_studies <- function(studies) {
  if (length(studies) == 0) {
    return("No linked study identified")
  }
  paste(studies, collapse = ", ")
}

collapse_study_breakdown <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return("No linked study identified")
  }
  paste0(df$study, " (n=", df$n, ")", collapse = "; ")
}

consort_exclusion_studies <- data.frame(
  exclusion_key = c(
    "missing_depression_endpoint",
    "no_sdoh_match",
    "no_eligible_mental_health_match"
  ),
  participants = c(
    length(ids_excl_missing_endpoint),
    length(ids_excl_no_sdoh_match),
    length(ids_excl_no_eligible_mh)
  ),
  studies = c(
    collapse_studies(get_study_list(ids_excl_missing_endpoint)),
    collapse_studies(get_study_list(ids_excl_no_sdoh_match)),
    collapse_studies(get_study_list(ids_excl_no_eligible_mh))
  ),
  study_breakdown = c(
    collapse_study_breakdown(get_study_breakdown(ids_excl_missing_endpoint)),
    collapse_study_breakdown(get_study_breakdown(ids_excl_no_sdoh_match)),
    collapse_study_breakdown(get_study_breakdown(ids_excl_no_eligible_mh))
  ),
  stringsAsFactors = FALSE
)

load(master_dict_path)
source(meta_helper_path)
meta <- get_variable_metadata()
labels <- unlist(meta$labels)

format_min_max <- function(x) {
  numeric_x <- suppressWarnings(as.numeric(x))
  numeric_x <- numeric_x[is.finite(numeric_x)]
  if (length(numeric_x) == 0) {
    return("NA")
  }
  paste0(
    format(min(numeric_x), trim = TRUE, scientific = FALSE),
    " to ",
    format(max(numeric_x), trim = TRUE, scientific = FALSE)
  )
}

build_dictionary <- function(dict_df, set_column, labels_lookup, data_df) {
  cols <- c("variable", "label", "type", "value_label")
  out <- dict_df[dict_df[[set_column]] == 1, cols, drop = FALSE]
  out[] <- lapply(out, function(x) if (is.factor(x)) as.character(x) else x)
  out <- unique(out)
  out <- out[!(is.na(out$variable) | trimws(out$variable) == ""), , drop = FALSE]
  rownames(out) <- NULL

  missing_label <- is.na(out$label) | trimws(out$label) == ""
  out$label[missing_label] <- ifelse(
    out$variable[missing_label] %in% names(labels_lookup),
    labels_lookup[out$variable[missing_label]],
    out$variable[missing_label]
  )

  out$value_label[is.na(out$value_label)] <- ""
  numeric_rows <- which(tolower(out$type) == "numeric")
  for (i in numeric_rows) {
    var_name <- out$variable[i]
    range_text <- "Not available"
    if (var_name %in% names(data_df)) {
      range_text <- format_min_max(data_df[[var_name]])
    }
    if (trimws(out$value_label[i]) == "") {
      out$value_label[i] <- paste0("Min-Max: ", range_text)
    } else {
      out$value_label[i] <- paste0(out$value_label[i], " | Min-Max: ", range_text)
    }
  }
  out
}

extract_scalar <- function(text, pattern, default = "") {
  match_obj <- regexec(pattern, text, perl = TRUE)
  match_vec <- regmatches(text, match_obj)[[1]]
  if (length(match_vec) >= 2) {
    return(match_vec[2])
  }
  default
}

extract_c_vector <- function(text, name) {
  pattern <- paste0("(?s)", name, "\\s*<-\\s*c\\((.*?)\\)")
  match_obj <- regexec(pattern, text, perl = TRUE)
  match_vec <- regmatches(text, match_obj)[[1]]
  if (length(match_vec) < 2) {
    return(character())
  }
  block <- match_vec[2]
  vals <- regmatches(block, gregexpr("\"[^\"]+\"", block, perl = TRUE))[[1]]
  gsub("^\"|\"$", "", vals)
}

extract_registry_variables <- function(text) {
  hits <- regmatches(
    text,
    gregexpr("list\\(variable\\s*=\\s*\"[^\"]+\"", text, perl = TRUE)
  )[[1]]
  if (length(hits) == 0) {
    return(character())
  }
  gsub("^.*\"([^\"]+)\"$", "\\1", hits)
}

build_role_dictionary <- function(dict_df, role_spec, data_df, labels_lookup) {
  role_spec <- role_spec[role_spec$variable != "", , drop = FALSE]
  role_spec <- role_spec[order(role_spec$priority, role_spec$order_id), , drop = FALSE]
  role_spec <- role_spec[!duplicated(role_spec$variable), , drop = FALSE]

  present <- role_spec$variable %in% names(data_df)
  present_vars <- role_spec[present, , drop = FALSE]

  infer_type <- function(x) {
    if (is.factor(x)) return("factor")
    if (is.numeric(x)) return("numeric")
    if (inherits(x, "Date") || inherits(x, "POSIXct") || inherits(x, "POSIXt")) return("date")
    if (is.logical(x)) return("logical")
    "character"
  }

  infer_value_label <- function(x, type_name) {
    if (identical(type_name, "numeric")) {
      return(paste0("Min-Max: ", format_min_max(x)))
    }
    non_missing <- x[!is.na(x)]
    if (length(non_missing) == 0) return("")
    lvls <- unique(as.character(non_missing))
    lvls <- lvls[trimws(lvls) != ""]
    if (length(lvls) == 0) return("")
    if (length(lvls) > 12) {
      return(paste0(paste(lvls[1:12], collapse = "; "), "; ..."))
    }
    paste(lvls, collapse = "; ")
  }

  build_row <- function(var_name, role_name) {
    dict_idx <- match(var_name, dict_df$variable)
    if (!is.na(dict_idx)) {
      row <- dict_df[dict_idx, c("variable", "label", "type", "value_label"), drop = FALSE]
    } else {
      vec <- data_df[[var_name]]
      type_name <- infer_type(vec)
      row <- data.frame(
        variable = var_name,
        label = ifelse(
          var_name %in% names(labels_lookup),
          labels_lookup[[var_name]],
          var_name
        ),
        type = type_name,
        value_label = infer_value_label(vec, type_name),
        stringsAsFactors = FALSE
      )
    }
    row$analysis_role <- role_name
    row[, c("analysis_role", "variable", "label", "type", "value_label"), drop = FALSE]
  }

  if (nrow(present_vars) == 0) {
    out <- data.frame(
      analysis_role = character(),
      variable = character(),
      label = character(),
      type = character(),
      value_label = character(),
      stringsAsFactors = FALSE
    )
  } else {
    row_list <- lapply(seq_len(nrow(present_vars)), function(i) {
      build_row(present_vars$variable[i], present_vars$role[i])
    })
    out <- do.call(rbind, row_list)
  }
  rownames(out) <- NULL

  missing <- role_spec[!present, c("role", "variable"), drop = FALSE]
  names(missing) <- c("analysis_role", "variable")
  rownames(missing) <- NULL

  list(dictionary = out, missing = missing)
}

dt_chilhood_exposure_dictionary <- build_dictionary(
  dt_master_dict,
  "Childhood Exposure",
  labels,
  dt_chilhood_exposure
)
dt_multistudies_dictionary <- build_dictionary(
  dt_master_dict,
  "Multistudies",
  labels,
  dt_multistudies
)

if (!("CLYR" %in% dt_multistudies_dictionary$variable)) {
  dt_multistudies_dictionary <- rbind(
    dt_multistudies_dictionary,
    data.frame(
      variable = "CLYR",
      label = "Calendar Year",
      type = "numeric",
      value_label = paste0("Min-Max: ", format_min_max(dt_multistudies$CLYR)),
      stringsAsFactors = FALSE
    )
  )
}

analysis_cfg_path <- file.path(
  "..",
  "OBJ00-Datasets Preperation",
  "statistical_analysis",
  "programs",
  "utils",
  "analysis_configuration.R"
)
analysis_cfg_text <- paste(readLines(analysis_cfg_path, warn = FALSE), collapse = "\n")

cfg_id_var <- extract_scalar(analysis_cfg_text, "id_var\\s*=\\s*\"([^\"]+)\"")
cfg_time_var <- extract_scalar(analysis_cfg_text, "time_var\\s*=\\s*\"([^\"]+)\"")
cfg_registry_vars <- extract_registry_variables(analysis_cfg_text)

cfg_mediators <- extract_c_vector(analysis_cfg_text, "mediators_vars")
cfg_baseline <- extract_c_vector(analysis_cfg_text, "baseline_vars")
cfg_confounders <- extract_c_vector(analysis_cfg_text, "confounders_vars")
cfg_outcomes <- extract_c_vector(analysis_cfg_text, "outcome_var")

# Metadata presentation override: classify SCSP and SCHL as mediators.
role_override_to_mediator <- c("SCSP", "SCHL")
override_vars <- role_override_to_mediator[role_override_to_mediator %in% cfg_confounders]
cfg_mediators <- unique(c(cfg_mediators, override_vars))
cfg_confounders <- setdiff(cfg_confounders, role_override_to_mediator)

# Metadata presentation override: exclude selected parental-status fields
# from Baseline Factor role grouping.
baseline_exclude_vars <- c("MDTH", "FDTH", "MENM", "FENM", "MHMS", "FHMS")
cfg_baseline <- setdiff(cfg_baseline, baseline_exclude_vars)

dt_chilhood_role_spec <- rbind(
  data.frame(variable = cfg_id_var, role = "Identifier", priority = 1, order_id = 1),
  data.frame(variable = cfg_time_var, role = "Time Variable", priority = 2, order_id = 1),
  data.frame(
    variable = cfg_registry_vars,
    role = "Exposure (SDoH Trajectory Input)",
    priority = 3,
    order_id = seq_along(cfg_registry_vars)
  ),
  data.frame(
    variable = cfg_baseline,
    role = "Baseline Factor",
    priority = 4,
    order_id = seq_along(cfg_baseline)
  )
)

dt_multistudies_role_spec <- rbind(
  data.frame(variable = cfg_id_var, role = "Identifier", priority = 1, order_id = 1),
  data.frame(variable = cfg_outcomes, role = "Outcome", priority = 2, order_id = seq_along(cfg_outcomes)),
  data.frame(variable = cfg_mediators, role = "Mediator", priority = 3, order_id = seq_along(cfg_mediators)),
  data.frame(
    variable = cfg_confounders,
    role = "Mediator-Outcome Confounder",
    priority = 4,
    order_id = seq_along(cfg_confounders)
  ),
  data.frame(
    variable = cfg_baseline,
    role = "Baseline Factor",
    priority = 5,
    order_id = seq_along(cfg_baseline)
  )
)

dt_chilhood_role_outputs <- build_role_dictionary(
  dt_chilhood_exposure_dictionary,
  dt_chilhood_role_spec,
  dt_chilhood_exposure,
  labels
)
dt_multistudies_role_outputs <- build_role_dictionary(
  dt_multistudies_dictionary,
  dt_multistudies_role_spec,
  dt_multistudies,
  labels
)

dt_chilhood_exposure_role_dictionary <- dt_chilhood_role_outputs$dictionary
dt_multistudies_role_dictionary <- dt_multistudies_role_outputs$dictionary
dt_chilhood_exposure_role_missing <- dt_chilhood_role_outputs$missing
dt_multistudies_role_missing <- dt_multistudies_role_outputs$missing

derived_info <- data.frame(
  variable = c("DPBN", "ELCDV", "MDCDV"),
  source_variables = c(
    "PHQBIN, irt_joint_models.rds",
    "Childhood exposure trajectories (0-5 years)",
    "Childhood exposure trajectories (6-12 years)"
  ),
  derivation_note = c(
    "PHQBIN >= 10",
    "Initialized as NA_character_ during OBJ00 merge; later populated from early-years latent trajectory class assignment.",
    "Initialized as NA_character_ during OBJ00 merge; later populated from middle-late years latent trajectory class assignment."
  ),
  derivation_status = c(
    "Implemented in OBJ00 merge",
    "Placeholder in OBJ00; populated in class-derivation stage",
    "Placeholder in OBJ00; populated in class-derivation stage"
  ),
  stringsAsFactors = FALSE
)

dt_multistudies_dictionary$is_derived <- ifelse(
  dt_multistudies_dictionary$variable %in% derived_info$variable,
  "Yes",
  "No"
)
match_idx <- match(dt_multistudies_dictionary$variable, derived_info$variable)
dt_multistudies_dictionary$source_variables <- derived_info$source_variables[match_idx]
dt_multistudies_dictionary$derivation_note <- derived_info$derivation_note[match_idx]
dt_multistudies_dictionary$derivation_status <- derived_info$derivation_status[match_idx]
dt_multistudies_dictionary$source_variables[is.na(dt_multistudies_dictionary$source_variables)] <- ""
dt_multistudies_dictionary$derivation_note[is.na(dt_multistudies_dictionary$derivation_note)] <- ""
dt_multistudies_dictionary$derivation_status[is.na(dt_multistudies_dictionary$derivation_status)] <- ""

dt_multistudies_nonderived_dictionary <- dt_multistudies_dictionary[
  dt_multistudies_dictionary$is_derived == "No",
  c("variable", "label", "type", "value_label"),
  drop = FALSE
]

dt_multistudies_derived_dictionary <- dt_multistudies_dictionary[
  dt_multistudies_dictionary$is_derived == "Yes",
  c(
    "variable",
    "label",
    "type",
    "value_label",
    "source_variables",
    "derivation_note",
    "derivation_status"
  ),
  drop = FALSE
]

metadata_summary <- data.frame(
  Field = c(
    "Project",
    "Population",
    "Linked sample size",
    "Data source backbone",
    "Core mental health outcomes",
    "Priority social determinants",
    "Cohorts nested in AHRI HDSS",
    "Youth Co-Creators"
  ),
  Value = c(
    "CO-LUMINATE",
    "South African youth aged 13-24",
    format(count_unique_ids(dt_multistudies, id_col_multi), big.mark = ",", trim = TRUE),
    "AHRI HDSS",
    "Depression symptoms (SSQ-14, PHQ-9)",
    "SES and caregiver co-residency",
    "DREAMS, Multilevel, Isisekelo Sempilo, TasP, Thetha Nami",
    "11"
  ),
  stringsAsFactors = FALSE
)

object_snapshot <- data.frame(
  Dataset = c(
    "Childhood Social Exposure Panel",
    "Integrated Multi-study Mental Health Panel"
  ),
  Rows = c(nrow(dt_chilhood_exposure), nrow(dt_multistudies)),
  Columns = c(ncol(dt_chilhood_exposure), ncol(dt_multistudies)),
  stringsAsFactors = FALSE
)

saveRDS(
  dt_chilhood_exposure_dictionary,
  file.path(objects_dir, "dt_chilhood_exposure_dictionary.rds")
)
saveRDS(
  dt_multistudies_dictionary,
  file.path(objects_dir, "dt_multistudies_dictionary.rds")
)
saveRDS(
  dt_chilhood_exposure_role_dictionary,
  file.path(objects_dir, "dt_chilhood_exposure_role_dictionary.rds")
)
saveRDS(
  dt_multistudies_role_dictionary,
  file.path(objects_dir, "dt_multistudies_role_dictionary.rds")
)
saveRDS(
  dt_chilhood_exposure_role_missing,
  file.path(objects_dir, "dt_chilhood_exposure_role_missing.rds")
)
saveRDS(
  dt_multistudies_role_missing,
  file.path(objects_dir, "dt_multistudies_role_missing.rds")
)
saveRDS(
  dt_multistudies_nonderived_dictionary,
  file.path(objects_dir, "dt_multistudies_nonderived_dictionary.rds")
)
saveRDS(
  dt_multistudies_derived_dictionary,
  file.path(objects_dir, "dt_multistudies_derived_dictionary.rds")
)
saveRDS(metadata_summary, file.path(objects_dir, "metadata_summary.rds"))
saveRDS(object_snapshot, file.path(objects_dir, "object_snapshot.rds"))
saveRDS(consort_counts, file.path(objects_dir, "consort_counts.rds"))
saveRDS(consort_exclusion_studies, file.path(objects_dir, "consort_exclusion_studies.rds"))

save(
  dt_chilhood_exposure_dictionary,
  file = file.path(objects_dir, "dt_chilhood_exposure_dictionary.RData")
)
save(
  dt_multistudies_dictionary,
  file = file.path(objects_dir, "dt_multistudies_dictionary.RData")
)
save(
  dt_chilhood_exposure_role_dictionary,
  file = file.path(objects_dir, "dt_chilhood_exposure_role_dictionary.RData")
)
save(
  dt_multistudies_role_dictionary,
  file = file.path(objects_dir, "dt_multistudies_role_dictionary.RData")
)
save(
  dt_chilhood_exposure_role_missing,
  file = file.path(objects_dir, "dt_chilhood_exposure_role_missing.RData")
)
save(
  dt_multistudies_role_missing,
  file = file.path(objects_dir, "dt_multistudies_role_missing.RData")
)
save(
  dt_multistudies_nonderived_dictionary,
  file = file.path(objects_dir, "dt_multistudies_nonderived_dictionary.RData")
)
save(
  dt_multistudies_derived_dictionary,
  file = file.path(objects_dir, "dt_multistudies_derived_dictionary.RData")
)
save(metadata_summary, file = file.path(objects_dir, "metadata_summary.RData"))
save(object_snapshot, file = file.path(objects_dir, "object_snapshot.RData"))
save(consort_counts, file = file.path(objects_dir, "consort_counts.RData"))
save(consort_exclusion_studies, file = file.path(objects_dir, "consort_exclusion_studies.RData"))

cat("Created objects in:", objects_dir, "\n")
cat("Source - dt_chilhood_exposure:", data_source_child, "\n")
cat("Source - dt_multistudies:", data_source_multi, "\n")
cat(
  "Rows (final aligned) - dt_chilhood_exposure:",
  nrow(dt_chilhood_exposure),
  "| dt_multistudies:",
  nrow(dt_multistudies),
  "\n"
)
cat(
  "Participants (final aligned) - dt_chilhood_exposure:",
  count_unique_ids(dt_chilhood_exposure, id_col_child),
  "| dt_multistudies:",
  count_unique_ids(dt_multistudies, id_col_multi),
  "\n"
)
cat(
  "DPBN missing in final dt_multistudies:",
  sum(is.na(dt_multistudies[[dpbn_col_multi]])),
  "\n"
)
cat(
  "Dictionary rows - dt_chilhood_exposure:",
  nrow(dt_chilhood_exposure_dictionary),
  "| dt_multistudies:",
  nrow(dt_multistudies_dictionary),
  "\n"
)
cat(
  "dt_multistudies split - non-derived:",
  nrow(dt_multistudies_nonderived_dictionary),
  "| derived:",
  nrow(dt_multistudies_derived_dictionary),
  "\n"
)
cat(
  "Role-filtered rows - dt_chilhood_exposure:",
  nrow(dt_chilhood_exposure_role_dictionary),
  "| dt_multistudies:",
  nrow(dt_multistudies_role_dictionary),
  "\n"
)
