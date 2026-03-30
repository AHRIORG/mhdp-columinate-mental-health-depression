# ==============================================================================
# OBJ00 Candidate Build Script (Non-destructive)
# Purpose:
#   - Assess and re-implement the OBJ00 core data flow for:
#       1) dt_childhood_exposure
#       2) dt_multistudies
#   - Keep the existing production script unchanged.
#   - Provide a cleaner, modular script that can be validated in parallel.
# ------------------------------------------------------------------------------
# NOTE:
#   This script does NOT replace current pipelines.
#   To avoid accidental overwrite, saving outputs is OFF by default.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(haven)
  library(labelled)
  library(readr)
  library(here)
})

# ------------------------------------------------------------------------------
# Path helpers
# ------------------------------------------------------------------------------

detect_obj00_root <- function() {
  candidates <- c(
    here::here(),
    here::here("OBJ00-Datasets Preperation")
  )

  for (root in candidates) {
    probe <- file.path(root, "data_management", "scripts", "99_utils", "metadata_helper.R")
    if (file.exists(probe)) {
      return(root)
    }
  }

  stop("Unable to detect OBJ00 root with metadata_helper.R available.")
}

OBJ00_ROOT <- detect_obj00_root()
obj00_path <- function(...) file.path(OBJ00_ROOT, ...)

source(obj00_path("data_management", "scripts", "99_utils", "metadata_helper.R"))

# ------------------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------------------

calculate_mode <- function(x) {
  x_clean <- stats::na.omit(x)
  if (length(x_clean) == 0) return(NA)
  ux <- unique(x_clean)
  ux[which.max(tabulate(match(x_clean, ux)))]
}

private_path <- function(cfg, ...) {
  normalizePath(file.path(OBJ00_ROOT, cfg$private_dt_dir, ...), winslash = "/", mustWork = FALSE)
}

build_config <- function(config = list()) {
  defaults <- list(
    private_dt_dir = "__SET_PRIVATE_DATA_DIR__",
    ext_path = "raw_data/hdss_raw",
    adam_dir = "adam",
    save_outputs = FALSE,
    output_dir = "data_management/data_examples/2_final_snippets/candidate_outputs",
    save_format = c("rds", "rdata"),
    return_intermediates = FALSE,
    verbose = TRUE,
    include_index_final_steps = TRUE,
    final_inputs_from_adam = TRUE,
    dt_psychometric = NULL,
    dt_dreams_multi_ssq = NULL,
    score_engines_rds = "data_management/data_examples/1_staging_snippets/derived_models/irt_joint_models.rds",
    scoring_script = "data_management/scripts/02_harmonization/02b_irt_model.R",
    best_calibration = NULL,
    save_to_private_adam = FALSE
  )
  modifyList(defaults, config)
}

log_msg <- function(cfg, ...) {
  if (isTRUE(cfg$verbose)) message(...)
}

load_rdata_object <- function(path, object_name = NULL) {
  if (!file.exists(path)) {
    stop("RData file not found: ", path)
  }

  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)

  if (!is.null(object_name)) {
    if (!exists(object_name, envir = env, inherits = FALSE)) {
      stop("Object '", object_name, "' was not found in ", path)
    }
    return(get(object_name, envir = env, inherits = FALSE))
  }

  if (length(loaded) == 1) {
    return(get(loaded[[1]], envir = env, inherits = FALSE))
  }

  out <- mget(loaded, envir = env, inherits = FALSE)
  out
}

# ------------------------------------------------------------------------------
# Input loading
# ------------------------------------------------------------------------------

load_primary_inputs <- function(cfg) {
  raw_file <- private_path(cfg, "raw_data", "columinate_clean.dta")
  hse_file <- private_path(
    cfg,
    cfg$ext_path,
    "SurveillanceEpisodesYrAgeDel_HSE_Labour_Partner_Parents_HIV.dta"
  )
  hh_file <- private_path(cfg, cfg$ext_path, "HouseholdMap.dta")
  hses_file <- private_path(cfg, cfg$ext_path, "HSE_scores.csv")

  required_files <- c(raw_file, hse_file, hh_file, hses_file)
  missing <- required_files[!file.exists(required_files)]
  if (length(missing) > 0) {
    stop(
      "Missing required input file(s). Missing: ",
      paste0(missing, collapse = "; ")
    )
  }

  log_msg(cfg, "Loading private input datasets...")

  dt_raw <- haven::read_dta(raw_file) |> haven::as_factor()
  dt_hse <- haven::read_dta(hse_file) |> haven::as_factor()
  dt_hh <- haven::read_dta(hh_file) |> haven::as_factor()

  dt_hses <- readr::read_csv(hses_file, show_col_types = FALSE) |>
    mutate(
      hdssyr = lubridate::year(as.Date(VisitDate, format = "%d%b%Y")),
      ses3cat = factor(ses3cat, levels = c("Low", "Middle", "High"))
    )

  join_cols <- intersect(names(dt_hses), names(dt_hh))
  if (length(join_cols) > 0) {
    dt_hses <- left_join(dt_hses, dt_hh, by = join_cols)
  }

  dt_hses <- dt_hses |>
    select(any_of(c("HouseholdId", "hdssyr", "ses_score", "ses3cat"))) |>
    distinct()

  list(
    dt_raw = dt_raw,
    dt_hse = dt_hse,
    dt_hses = dt_hses
  )
}

# ------------------------------------------------------------------------------
# Childhood exposure flow
# ------------------------------------------------------------------------------

build_childhood_panel <- function(dt_raw, dt_hses) {
  dt_seed <- dt_raw |>
    mutate(year_of_birth = lubridate::year(dateofbirth)) |>
    select(
      IIntId, IndividualId, HouseholdId,
      sex, dateofbirth, year_of_birth,
      starts_with("p_"),
      starts_with("hh_")
    )

  vars <- names(dt_seed)
  suffix <- stringr::str_extract(vars, "[0-9.]+$")
  numeric_suffix <- suppressWarnings(as.numeric(suffix))
  common_vars <- vars[is.na(numeric_suffix)]
  episodes <- sort(unique(stats::na.omit(numeric_suffix)))

  DATA_episode <- NULL

  for (ep in episodes) {
    episode_vars <- vars[which(!is.na(numeric_suffix) & numeric_suffix == ep)]

    dt_ep <- dt_seed |>
      mutate(episode = ep) |>
      select(any_of(common_vars), episode, any_of(episode_vars)) |>
      rename_with(~ gsub(paste0(ep, "$"), "", .x)) |>
      select(-contains("3y")) |>
      distinct() |>
      rename(hdssyr = episode)

    if (is.null(DATA_episode)) {
      DATA_episode <- dt_ep
    } else {
      DATA_episode <- bind_rows(DATA_episode, dt_ep) |>
        arrange(IIntId, hdssyr) |>
        filter(
          !if_all(
            setdiff(names(DATA_episode), c(common_vars, "hdssyr")),
            is.na
          )
        )
    }
  }

  DATA_with_ses <- DATA_episode |>
    rename(ses_old = hh_ses) |>
    left_join(
      dt_hses |>
        rename(hh_ses = ses3cat, hh_sess = ses_score),
      by = c("HouseholdId", "hdssyr")
    ) |>
    select(-any_of(c("ses_old", "hh_sess"))) |>
    distinct()

  DATA <- DATA_with_ses |>
    mutate(
      birth_expmonth = round(
        as.numeric(difftime(
          as.Date(paste0(year_of_birth, "-12-31")),
          dateofbirth,
          units = "days"
        )) / 365.25 * 12,
        0
      ),
      expage = hdssyr - year_of_birth
    ) |>
    filter(!(hdssyr == year_of_birth & birth_expmonth < 7)) |>
    group_by(IIntId) |>
    mutate(
      minexpage = min(expage, na.rm = TRUE),
      maxexpage = max(expage, na.rm = TRUE)
    ) |>
    ungroup() |>
    relocate(minexpage, maxexpage, birth_expmonth, hdssyr, expage, .after = year_of_birth) |>
    filter(expage < 13)

  list(
    DATA = DATA,
    DATA_with_ses = DATA_with_ses
  )
}

build_parent_history <- function(dt_hh_dhs, dt_hse, dt_child_ids, parent_id) {
  stopifnot(parent_id %in% c("MotherId", "FatherId"))

  base_parent <- dt_hh_dhs |>
    select(all_of(parent_id), HouseholdId, ChildId, PDAT, CYBT, childage) |>
    filter(!is.na(.data[[parent_id]]))

  hse_parent <- dt_hse |>
    mutate(DoD = lubridate::year(as.Date(DoD))) |>
    select(
      YEPI = Episode,
      IndividualId, HouseholdId, hdssyr = CalendarYear,
      PAR_YOD = DoD,
      PAR_AGE = Age,
      PAR_DTH = Died,
      PAR_HMS = HHRelationshipTypeId,
      PAR_EMP = EmploymentStatus,
      PAR_PTS = PartnerStatus,
      PAR_IMG = InMigration,
      PAR_OMG = OutMigration,
      PAR_ENM = Enumeration,
      PAR_HIV = HIVStatus
    ) |>
    rename(!!parent_id := IndividualId) |>
    distinct()

  out <- base_parent |>
    full_join(hse_parent, by = c(parent_id, "HouseholdId")) |>
    filter(ChildId %in% dt_child_ids$ChildId) |>
    arrange(.data[[parent_id]], HouseholdId, hdssyr, YEPI) |>
    group_by(.data[[parent_id]]) |>
    fill(PAR_YOD, .direction = "downup") |>
    ungroup() |>
    filter(!if_all(starts_with("PAR_"), is.na)) |>
    distinct() |>
    group_by(.data[[parent_id]], HouseholdId, hdssyr) |>
    mutate(
      across(c(PAR_IMG, PAR_OMG, PAR_ENM), ~ ifelse(grepl("Yes", paste(as.character(.x), collapse = "")), "Yes", as.character(.x))),
      across(c(PAR_IMG, PAR_OMG, PAR_ENM), ~ ifelse(grepl("No", paste(.x, collapse = "")), "No", .x)),
      across(c(PAR_HIV), ~ ifelse(grepl("Positive|positive", paste(as.character(.x), collapse = "")), "Positive", as.character(.x))),
      across(c(PAR_HIV), ~ ifelse(grepl("Negative|negative", paste(.x, collapse = "")), "Negative", .x))
    ) |>
    fill(starts_with("PAR_"), .direction = "downup") |>
    arrange(PAR_HMS) |>
    mutate(
      PAR_HMS = as.character(PAR_HMS),
      PAR_HMS = ifelse(grepl("Head|Spouse", PAR_HMS), "Head/Spouse/Husband", PAR_HMS),
      PAR_HMS = ifelse(grepl("Son/Daughter-inlaw|Child|Grandchild", PAR_HMS), "Child", PAR_HMS),
      PAR_HMS = ifelse(grepl("Sibling|Grandparent|Parent", PAR_HMS), "Related", PAR_HMS),
      PAR_HMS = ifelse(grepl("Domestic|Unrelated", PAR_HMS), "Unrelated/Other", PAR_HMS)
    ) |>
    distinct() |>
    filter(!is.na(hdssyr)) |>
    arrange(.data[[parent_id]], HouseholdId, hdssyr, YEPI) |>
    filter(YEPI == max(YEPI, na.rm = TRUE)) |>
    ungroup() |>
    select(-YEPI)

  out
}

finalize_parent_baseline <- function(parent_history, parent_id, parent_prefix, parent_label) {
  parent_history |>
    mutate(
      PAR_HMS = factor(PAR_HMS, levels = c("Head/Spouse/Husband", "Child", "Related", "Unrelated/Other")),
      across(c(PAR_HIV, PAR_IMG, PAR_OMG, PAR_ENM), ~ factor(.x))
    ) |>
    filter(hdssyr <= CYBT + childage) |>
    filter(hdssyr >= CYBT + childage - 1) |>
    mutate(
      PAR_MIG = paste0(PAR_IMG, PAR_OMG),
      PAR_MIG = ifelse(grepl("Yes", PAR_MIG), "Yes", PAR_MIG),
      PAR_MIG = ifelse(!grepl("Yes", PAR_MIG) & grepl("No", PAR_MIG), "No", PAR_MIG),
      PAR_MIG = ifelse(!grepl("Yes|No", PAR_MIG), "Unknown/Missing", PAR_MIG),
      PAR_AGE = PAR_AGE - childage
    ) |>
    group_by(.data[[parent_id]], HouseholdId) |>
    filter(hdssyr == max(hdssyr, na.rm = TRUE)) |>
    ungroup() |>
    set_variable_labels(
      PAR_AGE = paste0(parent_label, " Age"),
      PAR_HMS = paste0(parent_label, " Head of HH relationship"),
      PAR_EMP = paste0(parent_label, " Employment Status"),
      PAR_PTS = paste0(parent_label, " Partner Status"),
      PAR_MIG = paste0(parent_label, " Ever Migrated (Previous 1 Year)"),
      PAR_ENM = paste0(parent_label, " Enumeration (Previous 1 Year)"),
      PAR_HIV = paste0(parent_label, " HIV Status")
    ) |>
    rename_with(~ gsub("PAR_", parent_prefix, .x), starts_with("PAR_"))
}

build_parent_attributes <- function(DATA, DATA_with_ses, dt_hse, meta) {
  dt_child_ids <- DATA_with_ses |>
    select(ChildId = IndividualId, HouseholdId) |>
    distinct()

  dt_hse_ids <- dt_hse |>
    select(IndividualId, HouseholdId, MotherId, FatherId) |>
    distinct() |>
    right_join(dt_child_ids |> rename(IndividualId = ChildId), by = c("IndividualId", "HouseholdId")) |>
    mutate(PDAT = ifelse(!if_all(c(MotherId, FatherId), is.na), "Yes", "No")) |>
    arrange(PDAT) |>
    rename(ChildId = IndividualId)

  dt_hh_dhs <- DATA |>
    mutate(childage = hdssyr - year_of_birth) |>
    select(ChildId = IndividualId, HouseholdId, hdssyr, CYBT = year_of_birth, childage) |>
    distinct() |>
    left_join(dt_hse_ids, by = c("ChildId", "HouseholdId")) |>
    group_by(ChildId) |>
    filter(childage == min(childage, na.rm = TRUE)) |>
    ungroup()

  mother_child_raw <- build_parent_history(dt_hh_dhs, dt_hse, dt_child_ids, parent_id = "MotherId")
  father_child_raw <- build_parent_history(dt_hh_dhs, dt_hse, dt_child_ids, parent_id = "FatherId")

  mother_child <- finalize_parent_baseline(
    parent_history = mother_child_raw,
    parent_id = "MotherId",
    parent_prefix = "M",
    parent_label = "Mother's"
  )

  father_child <- finalize_parent_baseline(
    parent_history = father_child_raw,
    parent_id = "FatherId",
    parent_prefix = "F",
    parent_label = "Father's"
  )

  dt_blvrs <- mother_child |>
    rename(IndividualId = ChildId, expage = childage) |>
    select(IndividualId, contains("AGE"), contains("HMS"), contains("EMP"), contains("PTS"), contains("MIG"), contains("ENM"), contains("HIV")) |>
    full_join(
      father_child |>
        rename(IndividualId = ChildId, expage = childage) |>
        select(IndividualId, contains("AGE"), contains("HMS"), contains("EMP"), contains("PTS"), contains("MIG"), contains("ENM"), contains("HIV")),
      by = c("IndividualId", "expage")
    ) |>
    distinct() |>
    select(-expage) |>
    set_variable_labels(.labels = meta$labels, .strict = FALSE)

  list(
    mother_child = mother_child,
    father_child = father_child,
    dt_blvrs = dt_blvrs
  )
}

build_dt_childhood_exposure <- function(DATA, dt_blvrs, meta) {
  DATA |>
    left_join(dt_blvrs, by = "IndividualId") |>
    select(-IndividualId) |>
    rename(!!!meta$rename_childhood) |>
    select(-any_of(c("MINEXPYR", "MAXEXPYR"))) |>
    mutate(
      across(c(HHCHLD14), ~ ifelse(EXPAGE <= 14 & (.x == 0 | is.na(.x)), 1, .x)),
      across(c(HHCHLD4), ~ ifelse(EXPAGE <= 4 & (.x == 0 | is.na(.x)), 1, .x))
    ) |>
    mutate(
      HHCHLD14O = HHCHLD14,
      HHCHLD4O = HHCHLD4,
      across(c(HHCHLD14O, HHCHLD4O), ~ ifelse(.x <= 3, 1, 2)),
      across(
        c(HHCHLD14O, HHCHLD4O),
        ~ factor(.x, levels = c(1, 2), labels = c("No", "Yes"))
      ),
      across(
        c(MEDU, FEDU),
        ~ factor(
          as.numeric(.x),
          levels = c(1, 2, 3, 4),
          labels = c("None/Primary", "None/Primary", "Secondary", "Completed Matric")
        )
      )
    ) |>
    set_variable_labels(.labels = meta$labels, .strict = FALSE) |>
    distinct()
}

# ------------------------------------------------------------------------------
# Multi-studies flow
# ------------------------------------------------------------------------------

build_dt_multistudies <- function(dt_raw, meta) {
  vars <- names(dt_raw)
  suffix <- str_extract(vars, "[0-9.]+$")
  numeric_suffix <- suppressWarnings(as.numeric(suffix))
  common_vars <- vars[is.na(numeric_suffix)]
  episodes <- sort(unique(stats::na.omit(numeric_suffix)))

  DATA1 <- NULL
  DATA2 <- NULL

  for (ep in episodes) {
    episode_vars <- vars[which(!is.na(numeric_suffix) & numeric_suffix == ep)]

    dt_ep <- dt_raw |>
      mutate(episode = ep) |>
      select(any_of(common_vars), episode, any_of(episode_vars)) |>
      rename_with(~ gsub(paste0(ep, "$"), "", .x))

    if (ep <= 10) {
      dt_ep <- dt_ep |> mutate(visityr = lubridate::year(visitdate))
      if (is.null(DATA1)) {
        DATA1 <- dt_ep
      } else {
        DATA1 <- bind_rows(DATA1, dt_ep) |>
          arrange(IndividualId, IIntId, visitdate) |>
          filter(!if_all(setdiff(names(DATA1), c(common_vars, "episode")), is.na))
      }
    } else {
      dt_ep <- dt_ep |> rename(visityr = episode)
      if (is.null(DATA2)) {
        DATA2 <- dt_ep
      } else {
        DATA2 <- bind_rows(DATA2, dt_ep) |>
          arrange(IndividualId, IIntId, visityr)
      }
    }
  }

  dt_analysis_multistudies <- left_join(DATA1, DATA2, by = intersect(names(DATA1), names(DATA2))) |>
    arrange(IndividualId, IIntId, visitdate) |>
    group_by(IndividualId, IIntId) |>
    mutate(
      neg_fill = if_else(hiv_status == "Negative", "Negative", NA_character_),
      pos_fill = if_else(hiv_status == "Positive", "Positive", NA_character_),
      pos_fill = if_else(hiv_status == "StatusUnknown", NA_character_, pos_fill)
    ) |>
    fill(pos_fill, .direction = "down") |>
    fill(neg_fill, .direction = "up") |>
    ungroup() |>
    mutate(hiv_status = factor(coalesce(neg_fill, pos_fill, as.character(hiv_status)))) |>
    select(-any_of(c("episode", "visityr", "neg_fill", "pos_fill")))

  cols_to_keep <- c(
    "IIntId", "IndividualId", "HouseholdId", "dateofbirth", "visitdate", "sex", "age",
    "urban_or_rural", "migration", "fatherstatus", "motherstatus", "orphanstatus",
    "clinic_data", "round", "study_name", "education", "employment", "governmentgrant",
    "everhadsex", "hiv_status", "violence", "social_support", "foodsec", "everpregnant", "sexual_behavior", "condomlesssex",
    "everdrankalcohol", "depression_score", "depression_outcome", "suicidal_ideation_score",
    "suicidal_ideation_outcome", "ssq14_score", "ssq14_outcome",
    "suicidal_ideation_score_source"
  )

  clean_special_factor <- function(x) {
    levels_clean <- unique(gsub("NoInformation|NotApplicable|Not Applicable", NA_character_, as.character(x)))
    levels_clean <- sort(levels_clean[!is.na(levels_clean)])
    factor(as.character(x), levels = levels_clean)
  }

  dt_multistudies <- dt_analysis_multistudies |>
    select(any_of(cols_to_keep)) |>
    rename(!!!meta$rename_multistudies) |>
    filter(grepl("baseline", VISIT, ignore.case = TRUE)) |>
    select(-VISIT) |>
    mutate(
      PIPSA = ifelse(grepl("PIPSA", MGRC), "Yes", "No"),
      EXTMG = ifelse(grepl("External migration", MGRC), "Yes", "No"),

      FDEC = ifelse(grepl("Deceased", FDEC), "Yes", "No"),
      MDEC = ifelse(grepl("Deceased", MDEC), "Yes", "No"),

      across(any_of(c("CDMLESSFL", "DNKA", "FDSC", "VLNC", "URBANCAT")), ~ factor(trimws(as.character(.x)), levels = sort(unique(trimws(as.character(.x)))))),

      URBAN = ifelse(URBANCAT == "Urban", "Yes", "No"),
      PURBN = ifelse(URBANCAT == "Peri-Urban", "Yes", "No"),
      RURAL = ifelse(URBANCAT == "Rural", "Yes", "No"),
      SCHL = ifelse(SCHL == "Inschool", "Yes", "No"),
      GOVG = factor(abs(as.numeric(GOVG) - 2), levels = c(0, 1), labels = c("No", "Yes")),

      across(any_of(c("PIPSA", "EXTMG", "SCHL", "ESXC", "EVRPREGCAT", "URBAN", "PURBN", "RURAL", "MDEC", "FDEC")), ~ factor(.x, levels = c("No", "Yes"))),

      PHQBIN = clean_special_factor(PHQBIN),
      SCIBIN = clean_special_factor(SCIBIN),
      SSQBIN = clean_special_factor(SSQBIN),
      SCSP = clean_special_factor(SCSP),

      SEXBHV = factor(as.character(SEXBHV), levels = c("EverHadSex", "EverPregnant", "NeverHadSex"), labels = c("Ever had Sex", "Ever Pregnant", "Never had Sex")),
      HIVS = factor(as.character(HIVS), levels = c("Negative", "Positive"), labels = c("Negative", "Positive")),
      EMPLOY = factor(as.character(EMPLOY), levels = sort(unique(gsub("NoInformation", NA_character_, as.character(EMPLOY)))))
    ) |>
    mutate(
      CLYR = as.numeric(lubridate::year(VISITDT)),
      CLYR = CLYR - min(CLYR, na.rm = TRUE)
    ) |>
    relocate(USUBJID, HHID) |>
    relocate(PIPSA, EXTMG, .after = MGRC) |>
    relocate(URBAN, PURBN, RURAL, .after = URBANCAT) |>
    set_variable_labels(.labels = meta$labels, .strict = FALSE) |>
    set_variable_labels(HIVS = "HIV Status") |>
    arrange(USUBJID, HHID)

  dt_multistudies
}

# ------------------------------------------------------------------------------
# Optional finalization from OBJ00 index.qmd
# ------------------------------------------------------------------------------

resolve_index_final_inputs <- function(cfg) {
  dt_psychometric <- cfg$dt_psychometric
  dt_dreams_multi_ssq <- cfg$dt_dreams_multi_ssq

  if (is.null(dt_psychometric) && isTRUE(cfg$final_inputs_from_adam)) {
    dt_psychometric <- load_rdata_object(
      private_path(cfg, cfg$adam_dir, "dt_psychometric.RData"),
      object_name = "dt_psychometric"
    )
  }

  if (is.null(dt_dreams_multi_ssq) && isTRUE(cfg$final_inputs_from_adam)) {
    dt_dreams_multi_ssq <- load_rdata_object(
      private_path(cfg, cfg$adam_dir, "dt_dreams_multi_ssq.RData"),
      object_name = "dt_dreams_multi_ssq"
    )
  }

  if (is.null(dt_psychometric) || is.null(dt_dreams_multi_ssq)) {
    stop(
      "Missing inputs for index finalization. Provide dt_psychometric and dt_dreams_multi_ssq ",
      "via config, or enable final_inputs_from_adam with available private adam objects."
    )
  }

  list(
    dt_psychometric = dt_psychometric,
    dt_dreams_multi_ssq = dt_dreams_multi_ssq
  )
}

ensure_scoring_function <- function(cfg) {
  if (exists("score_new_cohort", mode = "function")) {
    return(invisible(TRUE))
  }

  scoring_script <- obj00_path(cfg$scoring_script)
  if (!file.exists(scoring_script)) {
    stop("Scoring script was not found: ", scoring_script)
  }

  source(scoring_script)

  if (!exists("score_new_cohort", mode = "function")) {
    stop("score_new_cohort() is not available after sourcing: ", scoring_script)
  }

  invisible(TRUE)
}

apply_index_final_steps <- function(dt_multistudies_core, dt_childhood_exposure, meta, cfg) {
  inputs <- resolve_index_final_inputs(cfg)
  dt_psychometric <- inputs$dt_psychometric
  dt_dreams_multi_ssq <- inputs$dt_dreams_multi_ssq

  ensure_scoring_function(cfg)

  score_engine_file <- obj00_path(cfg$score_engines_rds)
  if (!file.exists(score_engine_file)) {
    stop("Score engines object not found: ", score_engine_file)
  }

  SCORE_ENGINES <- readRDS(score_engine_file)
  if (length(SCORE_ENGINES) == 0) {
    stop("Score engines object is empty: ", score_engine_file)
  }

  LS_CALB <- list()

  for (calib_name in names(SCORE_ENGINES)) {
    dt_dreams_multi_ssq_base <- dt_dreams_multi_ssq |>
      select(-any_of(c(
        "CALIB", "TDEPBIN",
        "Theta_Harmonized", "Depression_Binary", "PHQ_Expected_Sum", "PHQ_Prob_GE10"
      )))

    dt_dreams_multi_ssq_scored <- dt_dreams_multi_ssq_base |>
      bind_cols(
        dt_dreams_multi_ssq_base |>
          score_new_cohort(scoring_engine = SCORE_ENGINES[[calib_name]]) |>
          mutate(
            CALIB = calib_name,
            Depression_Binary = factor(Depression_Binary, levels = c(0, 1), labels = c("No", "Yes"))
          ) |>
          ungroup()
      ) |>
      set_variable_labels(
        Theta_Harmonized = "IRT Score",
        Depression_Binary = "Binary (theta cutoff based)",
        PHQ_Expected_Sum = "PHQ-equivalent Score",
        PHQ_Prob_GE10 = "Probabilistic PHQ>=10"
      ) |>
      select(
        CALIB, USUBJID, STUDY, AGE, SEX, SSQ01:SSQ03, SSQ08:SSQ14,
        any_of("Depression_Binary")
      ) |>
      rename(TDEPBIN = Depression_Binary) |>
      set_variable_labels(TDEPBIN = "theta cutoff based Depression (Binary)")

    dt_multistudies_depression <- dt_multistudies_core |>
      left_join(
        dt_dreams_multi_ssq_scored |>
          select(USUBJID, STUDY, TDEPBIN),
        by = c("USUBJID", "STUDY")
      ) |>
      left_join(
        dt_psychometric |>
          mutate(PHQBIN2 = ifelse(PHQSCR >= 10, "Yes", "No")) |>
          select(USUBJID, PHQBIN2),
        by = "USUBJID"
      ) |>
      mutate(
        PHQBIN = ifelse(!is.na(PHQBIN2), as.character(PHQBIN2), as.character(PHQBIN)),
        DPBN = ifelse(is.na(PHQBIN), as.character(TDEPBIN), as.character(PHQBIN)),
        DPBN = factor(DPBN, levels = c("No", "Yes"))
      ) |>
      select(USUBJID, STUDY, DPBN)

    LS_CALB[[calib_name]] <- list(
      dt_dreams_multi_ssq = dt_dreams_multi_ssq_scored,
      dt_multistudies_depression = dt_multistudies_depression
    )
  }

  best_calibration <- cfg$best_calibration
  if (is.null(best_calibration) || !best_calibration %in% names(LS_CALB)) {
    best_calibration <- names(LS_CALB)[[1]]
  }

  dt_multistudies_final <- dt_multistudies_core |>
    left_join(
      LS_CALB[[best_calibration]]$dt_multistudies_depression,
      by = c("USUBJID", "STUDY")
    ) |>
    mutate(
      ELCDV = NA_character_,
      MDCDV = NA_character_,
      GBCDV = NA_character_
    ) |>
    select(
      USUBJID,
      STUDY, VISITDT, SEX, AGE,
      URBANCAT, URBAN, PURBN, RURAL,
      MGRC, EXTMG, PIPSA,
      ORPH, MDEC, FDEC,
      SCHL, GOVG, SCSP,
      EVRPREGCAT, SEXBHV, CDMLESSFL, ESXC, HIVS,
      DNKA, VLNC, FDSC, CLYR,
      DPBN, ELCDV, MDCDV, GBCDV,
      PHQSCR, SSQSCR
    ) |>
    set_variable_labels(
      .labels = meta$labels,
      .strict = FALSE
    ) |>
    set_variable_labels(
      DPBN = "Depression Outcome",
      ELCDV = "Early Childhood Development (Ages 0-5)",
      MDCDV = "Middle-Late Childhood Development (Ages 6-12)",
      GBCDV = "Entire Childhood Development (Ages 0-12)"
    )

  # Consort-ready alignment:
  # 1) keep only participants with a defined DPBN endpoint
  # 2) ensure both datasets contain the same participants
  dt_multistudies_final <- dt_multistudies_final |>
    filter(!is.na(DPBN))

  id_overlap <- intersect(
    unique(dt_multistudies_final$USUBJID),
    unique(dt_childhood_exposure$USUBJID)
  )

  dt_multistudies_final <- dt_multistudies_final |>
    filter(USUBJID %in% id_overlap)

  dt_childhood_exposure <- dt_childhood_exposure |>
    filter(USUBJID %in% id_overlap)

  list(
    dt_childhood_exposure = dt_childhood_exposure,
    dt_multistudies = dt_multistudies_final,
    dt_multistudies_core = dt_multistudies_core,
    dt_dreams_multi_ssq = LS_CALB[[best_calibration]]$dt_dreams_multi_ssq,
    best_calibration = best_calibration,
    calibration_outputs = LS_CALB
  )
}

# ------------------------------------------------------------------------------
# Validation + optional save
# ------------------------------------------------------------------------------

validate_outputs <- function(dt_childhood_exposure, dt_multistudies) {
  req_child <- c("USUBJID", "HHID", "BRTHDT", "BRTHYR", "HDSSYR", "EXPAGE")
  req_multi_core <- c("USUBJID", "HHID", "VISITDT", "STUDY", "AGE", "PHQSCR", "PHQBIN")
  req_multi_final <- c("USUBJID", "STUDY", "VISITDT", "SEX", "AGE", "DPBN", "PHQSCR", "SSQSCR")

  miss_child <- setdiff(req_child, names(dt_childhood_exposure))
  req_multi <- if ("DPBN" %in% names(dt_multistudies)) req_multi_final else req_multi_core
  miss_multi <- setdiff(req_multi, names(dt_multistudies))

  if (length(miss_child) > 0) {
    stop("dt_childhood_exposure is missing required column(s): ", paste(miss_child, collapse = ", "))
  }

  if (length(miss_multi) > 0) {
    stop("dt_multistudies is missing required column(s): ", paste(miss_multi, collapse = ", "))
  }

  if ("DPBN" %in% names(dt_multistudies)) {
    if (any(is.na(dt_multistudies$DPBN))) {
      stop("dt_multistudies contains missing DPBN values after finalization.")
    }

    id_child <- sort(unique(dt_childhood_exposure$USUBJID))
    id_multi <- sort(unique(dt_multistudies$USUBJID))
    if (!identical(id_child, id_multi)) {
      stop(
        "Participant mismatch between dt_childhood_exposure and dt_multistudies. ",
        "n_child=", length(id_child),
        ", n_multi=", length(id_multi)
      )
    }
  }

  invisible(TRUE)
}

save_candidate_outputs <- function(results, cfg) {
  out_dir <- file.path(OBJ00_ROOT, cfg$output_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if ("rds" %in% cfg$save_format) {
    saveRDS(results$dt_childhood_exposure, file.path(out_dir, "dt_childhood_exposure_candidate.rds"))
    saveRDS(results$dt_multistudies, file.path(out_dir, "dt_multistudies_candidate.rds"))
    if (!is.null(results$dt_multistudies_core)) {
      saveRDS(results$dt_multistudies_core, file.path(out_dir, "dt_multistudies_core_candidate.rds"))
    }
    if (!is.null(results$dt_dreams_multi_ssq)) {
      saveRDS(results$dt_dreams_multi_ssq, file.path(out_dir, "dt_dreams_multi_ssq_candidate.rds"))
    }
  }

  if ("rdata" %in% cfg$save_format) {
    dt_childhood_exposure <- results$dt_childhood_exposure
    dt_multistudies <- results$dt_multistudies
    save(dt_childhood_exposure, file = file.path(out_dir, "dt_childhood_exposure_candidate.RData"))
    save(dt_multistudies, file = file.path(out_dir, "dt_multistudies_candidate.RData"))
    if (!is.null(results$dt_multistudies_core)) {
      dt_multistudies_core <- results$dt_multistudies_core
      save(dt_multistudies_core, file = file.path(out_dir, "dt_multistudies_core_candidate.RData"))
    }
    if (!is.null(results$dt_dreams_multi_ssq)) {
      dt_dreams_multi_ssq <- results$dt_dreams_multi_ssq
      save(dt_dreams_multi_ssq, file = file.path(out_dir, "dt_dreams_multi_ssq_candidate.RData"))
    }
  }
}

save_private_adam_outputs <- function(results, cfg) {
  adam_dir <- private_path(cfg, cfg$adam_dir)
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

# ------------------------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------------------------

build_obj00_core_datasets <- function(config = list()) {
  cfg <- build_config(config)
  meta <- get_variable_metadata()

  log_msg(cfg, "Starting candidate OBJ00 core build...")
  inputs <- load_primary_inputs(cfg)

  log_msg(cfg, "Building childhood exposure panel...")
  childhood <- build_childhood_panel(inputs$dt_raw, inputs$dt_hses)

  log_msg(cfg, "Building parental baseline attributes...")
  parents <- build_parent_attributes(
    DATA = childhood$DATA,
    DATA_with_ses = childhood$DATA_with_ses,
    dt_hse = inputs$dt_hse,
    meta = meta
  )

  log_msg(cfg, "Assembling dt_childhood_exposure...")
  dt_childhood_exposure <- build_dt_childhood_exposure(
    DATA = childhood$DATA,
    dt_blvrs = parents$dt_blvrs,
    meta = meta
  )

  log_msg(cfg, "Assembling dt_multistudies...")
  dt_multistudies_core <- build_dt_multistudies(
    dt_raw = inputs$dt_raw,
    meta = meta
  )

  results <- list()
  if (isTRUE(cfg$include_index_final_steps)) {
    log_msg(cfg, "Applying final index.qmd steps (depression calibration + final selection)...")
    results <- apply_index_final_steps(
      dt_multistudies_core = dt_multistudies_core,
      dt_childhood_exposure = dt_childhood_exposure,
      meta = meta,
      cfg = cfg
    )
  } else {
    results <- list(
      dt_childhood_exposure = dt_childhood_exposure,
      dt_multistudies = dt_multistudies_core
    )
  }

  validate_outputs(results$dt_childhood_exposure, results$dt_multistudies)

  if (isTRUE(cfg$return_intermediates)) {
    results$intermediates <- list(
      childhood = childhood,
      parents = parents,
      dt_multistudies_core = dt_multistudies_core
    )
  }

  if (isTRUE(cfg$save_outputs)) {
    log_msg(cfg, "Saving candidate outputs...")
    save_candidate_outputs(results, cfg)
  }

  if (isTRUE(cfg$save_to_private_adam)) {
    log_msg(cfg, "Saving outputs to private adam path...")
    save_private_adam_outputs(results, cfg)
  }

  log_msg(cfg, "Candidate OBJ00 core build completed.")
  results
}

if (sys.nframe() == 0) {
  message("Script loaded.")
  message("Run: results <- build_obj00_core_datasets()")
  message("No existing production script is modified by this file.")
}
