#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(jsonlite)
})

project_root <- "."
objects_dir <- file.path(project_root, "_tools", "_objects")
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)

private_adam_roots <- c(
  file.path("..", "..", "_private_use", "data_management", "data", "co-luminate", "adam"),
  file.path("..", "_private_use", "data_management", "data", "co-luminate", "adam"),
  file.path("..", "..", "..", "_private_use", "data_management", "data", "co-luminate", "adam")
)
child_private_paths <- file.path(private_adam_roots, "dt_childhood_exposure.RData")
multi_private_paths <- file.path(private_adam_roots, "dt_multistudies.RData")

min_cell_n <- 10
dimensions_cfg_path <- file.path(project_root, "_tools", "_scripts", "dashboard_dimensions.csv")

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

find_col_name <- function(df, target_name) {
  idx <- match(toupper(target_name), toupper(names(df)))
  if (is.na(idx)) {
    return(NA_character_)
  }
  names(df)[idx]
}

to_binary <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  yes_values <- c("yes", "1", "true", "positive", "depressed", "in school", "inschool", "female")
  no_values <- c("no", "0", "false", "negative", "not depressed", "not in school", "not inschool", "male")

  out <- rep(NA_real_, length(x_chr))
  out[x_chr %in% yes_values] <- 1
  out[x_chr %in% no_values] <- 0

  x_num <- suppressWarnings(as.numeric(x_chr))
  out[is.na(out) & x_num %in% c(0, 1)] <- x_num[is.na(out) & x_num %in% c(0, 1)]
  out
}

recode_govg_no_grant <- function(x) {
  out <- suppressWarnings(
    factor(abs(as.numeric(x) - 2), levels = c(0, 1), labels = c("No", "Yes"))
  )

  if (all(is.na(out))) {
    x_chr <- tolower(trimws(as.character(x)))
    out_chr <- rep(NA_character_, length(x_chr))
    out_chr[x_chr %in% c("yes", "1", "true")] <- "No"
    out_chr[x_chr %in% c("no", "0", "false")] <- "Yes"
    out <- factor(out_chr, levels = c("No", "Yes"))
  }

  out
}

extract_year <- function(x) {
  if (inherits(x, "Date")) {
    return(as.numeric(format(x, "%Y")))
  }
  if (inherits(x, "POSIXt")) {
    return(as.numeric(format(as.Date(x), "%Y")))
  }

  x_chr <- as.character(x)
  y <- suppressWarnings(as.numeric(format(as.Date(x_chr), "%Y")))
  if (!all(is.na(y))) {
    return(y)
  }

  y2 <- suppressWarnings(as.numeric(substr(x_chr, 1, 4)))
  y2[!(y2 >= 1900 & y2 <= 2100)] <- NA_real_
  y2
}

extract_date <- function(x) {
  if (inherits(x, "Date")) {
    return(as.Date(x))
  }
  if (inherits(x, "POSIXt")) {
    return(as.Date(x))
  }

  x_chr <- trimws(as.character(x))
  x_chr[x_chr == ""] <- NA_character_

  out <- suppressWarnings(as.Date(x_chr))
  missing_idx <- which(is.na(out) & !is.na(x_chr))
  if (length(missing_idx) > 0) {
    parsed <- suppressWarnings(lubridate::parse_date_time(
      x_chr[missing_idx],
      orders = c("Ymd", "ymd", "Y-m-d", "dmy", "mdy", "d-b-Y", "d B Y", "d/m/Y", "m/d/Y"),
      quiet = TRUE
    ))
    out[missing_idx] <- as.Date(parsed)
  }

  out
}

format_month_year <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return("Not available")
  }
  format(as.Date(x)[1], "%b %Y")
}

format_year <- function(x) {
  if (length(x) == 0 || all(is.na(x)) || !is.finite(x[1])) {
    return("Not available")
  }
  as.character(as.integer(round(x[1])))
}

clean_study <- function(x) {
  out <- trimws(as.character(x))
  out[out == "" | is.na(out)] <- "Unknown"
  out
}

clean_sex <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- ifelse(
    x_chr %in% c("female", "f", "2"),
    "Female",
    ifelse(
      x_chr %in% c("male", "m", "1"),
      "Male",
      "Unknown"
    )
  )
  out[is.na(x_chr) | x_chr == ""] <- "Unknown"
  out
}

band_age <- function(x) {
  out <- rep("Unknown", length(x))
  out[is.finite(x) & x >= 13 & x <= 15] <- "13-15"
  out[is.finite(x) & x >= 16 & x <= 18] <- "16-18"
  out[is.finite(x) & x >= 19 & x <= 21] <- "19-21"
  out[is.finite(x) & x >= 22 & x <= 24] <- "22-24"
  out[is.finite(x) & x > 24] <- "25+"
  factor(out, levels = c("13-15", "16-18", "19-21", "22-24", "25+", "Unknown"))
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  stats::median(x)
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE))
}

sanitize_colname <- function(x) {
  out <- tolower(gsub("[^a-zA-Z0-9]+", "_", x))
  out <- gsub("^_+|_+$", "", out)
  ifelse(out == "", "var", out)
}

clean_level <- function(x) {
  out <- trimws(as.character(x))
  out[out == "" | is.na(out)] <- "Unknown"
  out
}

if (!file.exists(dimensions_cfg_path)) {
  stop("Missing dashboard dimensions config: ", dimensions_cfg_path)
}
dimensions_cfg <- read.csv(dimensions_cfg_path, stringsAsFactors = FALSE, check.names = FALSE)
required_cfg_cols <- c("domain", "var", "label", "type")
if (!all(required_cfg_cols %in% names(dimensions_cfg))) {
  stop(
    "dashboard_dimensions.csv must contain columns: ",
    paste(required_cfg_cols, collapse = ", ")
  )
}
dimensions_cfg <- dimensions_cfg |>
  mutate(
    domain = trimws(domain),
    var = trimws(var),
    label = trimws(label),
    type = tolower(trimws(type))
  ) |>
  filter(var != "", label != "")

multi_loaded <- load_first_df(
  multi_private_paths,
  preferred_name = "dt_multistudies"
)
if (is.null(multi_loaded$data)) {
  stop("Unable to load private dt_multistudies.RData for dashboard aggregation.")
}
dt_multistudies <- multi_loaded$data
data_source_multi <- multi_loaded$source
n_input_rows_loaded <- nrow(dt_multistudies)

id_col_multi <- find_col_name(dt_multistudies, "USUBJID")
if (is.na(id_col_multi)) {
  dt_multistudies$USUBJID <- seq_len(nrow(dt_multistudies))
  id_col_multi <- "USUBJID"
} else if (id_col_multi != "USUBJID") {
  dt_multistudies$USUBJID <- dt_multistudies[[id_col_multi]]
}

study_col_multi <- find_col_name(dt_multistudies, "STUDY")
if (is.na(study_col_multi)) {
  stop("dt_multistudies is missing required column: STUDY")
}
if (study_col_multi != "STUDY") {
  dt_multistudies$STUDY <- dt_multistudies[[study_col_multi]]
}

dpbn_col_multi <- find_col_name(dt_multistudies, "DPBN")
if (is.na(dpbn_col_multi)) {
  stop("dt_multistudies is missing required column: DPBN")
}
if (dpbn_col_multi != "DPBN") {
  dt_multistudies$DPBN <- dt_multistudies[[dpbn_col_multi]]
}

dt_multistudies <- dt_multistudies |>
  filter(!is.na(DPBN))

n_rows_after_dpbn <- nrow(dt_multistudies)
n_overlap_ids <- NA_integer_
dt_childhood_overlap <- NULL

child_loaded <- load_first_df(
  child_private_paths,
  preferred_name = "dt_childhood_exposure"
)
data_source_child <- child_loaded$source

if (!is.null(child_loaded$data)) {
  dt_childhood <- child_loaded$data
  id_col_child <- find_col_name(dt_childhood, "USUBJID")

  if (!is.na(id_col_child)) {
    ids_multi <- unique(as.character(dt_multistudies$USUBJID))
    ids_child <- unique(as.character(dt_childhood[[id_col_child]]))
    overlap_ids <- intersect(ids_multi, ids_child)
    overlap_ids <- overlap_ids[!is.na(overlap_ids) & nzchar(overlap_ids)]
    n_overlap_ids <- length(overlap_ids)

    dt_childhood_overlap <- dt_childhood |>
      filter(as.character(.data[[id_col_child]]) %in% overlap_ids)

    dt_multistudies <- dt_multistudies |>
      filter(as.character(USUBJID) %in% overlap_ids)
  }
}

input_colnames <- names(dt_multistudies)

mediators <- c("DNKA", "ESXC", "FDSC", "GOVG", "VLNC", "SCHL", "SCSP")
mediator_labels <- c(
  DNKA = "Drank Alcohol",
  ESXC = "Ever Had Sex",
  FDSC = "Food Insecure",
  GOVG = "No Government Grant",
  VLNC = "Experienced Violence",
  SCHL = "In School",
  SCSP = "Has Social Support"
)

for (med in mediators) {
  if (!(med %in% names(dt_multistudies))) {
    dt_multistudies[[med]] <- NA
  }
}

cfg_vars <- unique(dimensions_cfg$var)
for (v in cfg_vars) {
  if (!(v %in% names(dt_multistudies))) {
    dt_multistudies[[v]] <- NA
  }
}

analysis_data <- dt_multistudies |>
  mutate(
    STUDY_clean = clean_study(STUDY),
    SEX_clean = clean_sex(SEX),
    AGE_num = suppressWarnings(as.numeric(AGE)),
    AGE_band = band_age(AGE_num),
    DPBN_bin = to_binary(DPBN),
    VISIT_date = extract_date(VISITDT),
    VISIT_year = extract_year(VISITDT),
    GOVG = recode_govg_no_grant(GOVG),
    DNKA_bin = to_binary(DNKA),
    ESXC_bin = to_binary(ESXC),
    FDSC_bin = to_binary(FDSC),
    GOVG_bin = to_binary(GOVG),
    VLNC_bin = to_binary(VLNC),
    SCHL_bin = to_binary(SCHL),
    SCSP_bin = to_binary(SCSP),
    RURAL_bin = to_binary(RURAL),
    EXTMG_bin = to_binary(EXTMG)
  ) |>
  filter(!is.na(DPBN_bin))

visit_dates_non_missing <- analysis_data$VISIT_date[!is.na(analysis_data$VISIT_date)]
if (length(visit_dates_non_missing) > 0) {
  visit_date_min <- min(visit_dates_non_missing)
  visit_date_max <- max(visit_dates_non_missing)
} else {
  visit_date_min <- as.Date(NA)
  visit_date_max <- as.Date(NA)
}
date_start_month_year <- format_month_year(visit_date_min)
date_end_month_year <- format_month_year(visit_date_max)
date_range_month_year <- if (identical(date_start_month_year, "Not available") ||
  identical(date_end_month_year, "Not available")) {
  "Not available"
} else {
  paste0(date_start_month_year, " to ", date_end_month_year)
}

sdoh_start_year <- "Not available"
sdoh_end_year <- "Not available"
sdoh_range_year <- "Not available"
sdoh_start_month_year <- "Not available"
sdoh_end_month_year <- "Not available"
sdoh_range_month_year <- "Not available"
if (!is.null(dt_childhood_overlap)) {
  hdssyr_col <- find_col_name(dt_childhood_overlap, "HDSSYR")
  if (!is.na(hdssyr_col)) {
    hdssyr_num <- suppressWarnings(as.numeric(as.character(dt_childhood_overlap[[hdssyr_col]])))
    hdssyr_num <- hdssyr_num[is.finite(hdssyr_num)]
    if (length(hdssyr_num) > 0) {
      hdssyr_min <- min(hdssyr_num, na.rm = TRUE)
      hdssyr_max <- max(hdssyr_num, na.rm = TRUE)
      sdoh_start_year <- format_year(hdssyr_min)
      sdoh_end_year <- format_year(hdssyr_max)
      sdoh_range_year <- paste0(sdoh_start_year, " to ", sdoh_end_year)
      sdoh_start_month_year <- paste("Jan", sdoh_start_year)
      sdoh_end_month_year <- paste("Dec", sdoh_end_year)
      sdoh_range_month_year <- paste0(sdoh_start_month_year, " to ", sdoh_end_month_year)
    }
  }
}

visit_year_min <- suppressWarnings(min(analysis_data$VISIT_year, na.rm = TRUE))
if (!is.finite(visit_year_min)) {
  analysis_data$CLYR <- NA_real_
} else {
  analysis_data$CLYR <- analysis_data$VISIT_year - visit_year_min
}

participant_snapshot <- analysis_data |>
  mutate(
    visit_ord = ifelse(is.na(VISIT_year), -Inf, VISIT_year)
  ) |>
  arrange(USUBJID, desc(visit_ord)) |>
  group_by(USUBJID) |>
  slice_head(n = 1) |>
  ungroup()

study_summary <- participant_snapshot |>
  group_by(study = STUDY_clean) |>
  summarise(
    participants = n(),
    depressed_n = sum(DPBN_bin == 1, na.rm = TRUE),
    prevalence_pct = round(100 * depressed_n / participants, 1),
    female_pct = round(100 * mean(SEX_clean == "Female", na.rm = TRUE), 1),
    median_age = round(safe_median(AGE_num), 1),
    age_q1 = round(safe_quantile(AGE_num, 0.25), 1),
    age_q3 = round(safe_quantile(AGE_num, 0.75), 1),
    .groups = "drop"
  ) |>
  filter(participants >= min_cell_n) |>
  arrange(desc(participants), study)

demog_summary <- participant_snapshot |>
  group_by(study = STUDY_clean, sex = SEX_clean, age_band = as.character(AGE_band)) |>
  summarise(
    n = n(),
    depressed_n = sum(DPBN_bin == 1, na.rm = TRUE),
    prevalence_pct = round(100 * depressed_n / n, 1),
    .groups = "drop"
  ) |>
  filter(n >= min_cell_n) |>
  arrange(study, sex, age_band)

year_summary <- analysis_data |>
  filter(!is.na(CLYR)) |>
  group_by(study = STUDY_clean, clyr = as.integer(CLYR)) |>
  summarise(
    n = n(),
    depressed_n = sum(DPBN_bin == 1, na.rm = TRUE),
    prevalence_pct = round(100 * depressed_n / n, 1),
    .groups = "drop"
  ) |>
  filter(n >= min_cell_n) |>
  arrange(study, clyr)

mediator_summary <- participant_snapshot |>
  transmute(
    study = STUDY_clean,
    DNKA = DNKA_bin,
    ESXC = ESXC_bin,
    FDSC = FDSC_bin,
    GOVG = GOVG_bin,
    VLNC = VLNC_bin,
    SCHL = SCHL_bin,
    SCSP = SCSP_bin
  ) |>
  tidyr::pivot_longer(
    cols = c("DNKA", "ESXC", "FDSC", "GOVG", "VLNC", "SCHL", "SCSP"),
    names_to = "mediator",
    values_to = "value"
  ) |>
  mutate(
    mediator = unname(mediator_labels[mediator])
  ) |>
  filter(!is.na(value)) |>
  group_by(study, mediator) |>
  summarise(
    n_complete = n(),
    yes_n = sum(value == 1, na.rm = TRUE),
    yes_pct = round(100 * yes_n / n_complete, 1),
    .groups = "drop"
  ) |>
  filter(n_complete >= min_cell_n) |>
  arrange(study, desc(yes_pct))

filter_schema <- dimensions_cfg |>
  mutate(
    source_present = var %in% input_colnames,
    available = FALSE,
    non_missing_n = 0L,
    distinct_levels = 0L,
    column_name = paste0("f_", sanitize_colname(var))
  )

filter_cube_base <- participant_snapshot |>
  transmute(
    study = STUDY_clean,
    sex = SEX_clean,
    age_band = as.character(AGE_band),
    DPBN_bin = DPBN_bin
  )

available_filter_cols <- character()

for (i in seq_len(nrow(filter_schema))) {
  var_name <- filter_schema$var[[i]]
  var_type <- filter_schema$type[[i]]
  out_col <- filter_schema$column_name[[i]]

  raw_vec <- participant_snapshot[[var_name]]

  if (var_type == "binary") {
    bin_col <- paste0(var_name, "_bin")
    if (bin_col %in% names(participant_snapshot)) {
      bin_vec <- participant_snapshot[[bin_col]]
    } else {
      bin_vec <- to_binary(raw_vec)
    }
    out_vec <- ifelse(
      is.na(bin_vec),
      "Unknown",
      ifelse(bin_vec == 1, "Yes", ifelse(bin_vec == 0, "No", "Unknown"))
    )
    non_missing_n <- sum(!is.na(bin_vec))
    distinct_levels <- dplyr::n_distinct(out_vec[out_vec != "Unknown"])
  } else {
    out_vec <- clean_level(raw_vec)
    non_missing_n <- sum(!(is.na(raw_vec) | trimws(as.character(raw_vec)) == ""))
    distinct_levels <- dplyr::n_distinct(out_vec[out_vec != "Unknown"])
  }

  filter_schema$non_missing_n[[i]] <- as.integer(non_missing_n)
  filter_schema$distinct_levels[[i]] <- as.integer(distinct_levels)

  usable <- isTRUE(filter_schema$source_present[[i]]) && non_missing_n > 0 && distinct_levels > 0
  filter_schema$available[[i]] <- usable

  if (usable) {
    filter_cube_base[[out_col]] <- out_vec
    available_filter_cols <- c(available_filter_cols, out_col)
  }
}

cube_group_cols <- c("study", "sex", "age_band", available_filter_cols)
filter_cube <- filter_cube_base |>
  group_by(across(all_of(cube_group_cols))) |>
  summarise(
    n = n(),
    depressed_n = sum(DPBN_bin == 1, na.rm = TRUE),
    prevalence_pct = round(100 * depressed_n / n, 1),
    .groups = "drop"
  ) |>
  filter(n >= min_cell_n)

filter_schema_dimension_status <- filter_schema

if (nrow(filter_cube) == 0) {
  filter_cube <- filter_cube_base |>
    group_by(study, sex, age_band) |>
    summarise(
      n = n(),
      depressed_n = sum(DPBN_bin == 1, na.rm = TRUE),
      prevalence_pct = round(100 * depressed_n / n, 1),
      .groups = "drop"
    ) |>
    filter(n >= min_cell_n)

  filter_schema$available <- FALSE
  filter_schema$column_name <- NA_character_
}

trajectory_schema <- filter_schema_dimension_status |>
  filter(tolower(domain) == "trajectory", available)

trajectory_trend_summary <- data.frame(
  trajectory_var = character(),
  trajectory_label = character(),
  study = character(),
  clyr = integer(),
  trajectory_class = character(),
  n = integer(),
  depressed_n = integer(),
  prevalence_pct = numeric(),
  stringsAsFactors = FALSE
)

trajectory_pairwise_summary <- data.frame(
  pair_id = character(),
  pair_label = character(),
  trajectory_var_x = character(),
  trajectory_label_x = character(),
  class_x = character(),
  trajectory_var_y = character(),
  trajectory_label_y = character(),
  class_y = character(),
  study = character(),
  n = integer(),
  depressed_n = integer(),
  prevalence_pct = numeric(),
  stringsAsFactors = FALSE
)

if (nrow(trajectory_schema) > 0) {
  trend_parts <- vector("list", nrow(trajectory_schema))
  for (i in seq_len(nrow(trajectory_schema))) {
    var_name <- trajectory_schema$var[[i]]
    var_label <- trajectory_schema$label[[i]]

    trend_parts[[i]] <- analysis_data |>
      transmute(
        study = STUDY_clean,
        clyr = as.integer(CLYR),
        DPBN_bin = DPBN_bin,
        trajectory_class = clean_level(.data[[var_name]])
      ) |>
      filter(!is.na(clyr)) |>
      group_by(study, clyr, trajectory_class) |>
      summarise(
        n = n(),
        depressed_n = sum(DPBN_bin == 1, na.rm = TRUE),
        prevalence_pct = round(100 * depressed_n / n, 1),
        .groups = "drop"
      ) |>
      filter(n >= min_cell_n) |>
      mutate(
        trajectory_var = var_name,
        trajectory_label = var_label
      ) |>
      select(
        trajectory_var, trajectory_label,
        study, clyr, trajectory_class,
        n, depressed_n, prevalence_pct
      )
  }

  trajectory_trend_summary <- bind_rows(trend_parts) |>
    arrange(trajectory_label, study, trajectory_class, clyr)

  if (nrow(trajectory_schema) > 1) {
    trajectory_pairs <- combn(trajectory_schema$var, 2, simplify = FALSE)
    pair_parts <- vector("list", length(trajectory_pairs))

    for (i in seq_along(trajectory_pairs)) {
      pair <- trajectory_pairs[[i]]
      var_x <- pair[[1]]
      var_y <- pair[[2]]
      label_x <- trajectory_schema$label[match(var_x, trajectory_schema$var)]
      label_y <- trajectory_schema$label[match(var_y, trajectory_schema$var)]

      pair_parts[[i]] <- participant_snapshot |>
        transmute(
          study = STUDY_clean,
          DPBN_bin = DPBN_bin,
          class_x = clean_level(.data[[var_x]]),
          class_y = clean_level(.data[[var_y]])
        ) |>
        group_by(study, class_x, class_y) |>
        summarise(
          n = n(),
          depressed_n = sum(DPBN_bin == 1, na.rm = TRUE),
          prevalence_pct = round(100 * depressed_n / n, 1),
          .groups = "drop"
        ) |>
        filter(n >= min_cell_n) |>
        mutate(
          pair_id = paste(var_x, var_y, sep = "__"),
          pair_label = paste(label_x, "x", label_y),
          trajectory_var_x = var_x,
          trajectory_label_x = label_x,
          trajectory_var_y = var_y,
          trajectory_label_y = label_y
        ) |>
        select(
          pair_id, pair_label,
          trajectory_var_x, trajectory_label_x, class_x,
          trajectory_var_y, trajectory_label_y, class_y,
          study, n, depressed_n, prevalence_pct
        )
    }

    trajectory_pairwise_summary <- bind_rows(pair_parts) |>
      arrange(pair_label, study, desc(n), class_x, class_y)
  }
}

kpi <- list(
  participants_total = nrow(participant_snapshot),
  records_total = nrow(analysis_data),
  studies_total = dplyr::n_distinct(participant_snapshot$STUDY_clean),
  depression_prev_pct = round(100 * mean(participant_snapshot$DPBN_bin == 1, na.rm = TRUE), 1),
  mh_start_month_year = date_start_month_year,
  mh_end_month_year = date_end_month_year,
  mh_range_month_year = date_range_month_year,
  sdoh_start_month_year = sdoh_start_month_year,
  sdoh_end_month_year = sdoh_end_month_year,
  sdoh_range_month_year = sdoh_range_month_year,
  sdoh_start_year = sdoh_start_year,
  sdoh_end_year = sdoh_end_year,
  sdoh_range_year = sdoh_range_year,
  min_cell_n = min_cell_n,
  filter_dimensions_available = sum(filter_schema$available),
  trajectory_dimensions_available = nrow(trajectory_schema)
)

dashboard_build_info <- list(
  source = data_source_multi,
  source_child = data_source_child,
  rendered_at = as.character(Sys.time()),
  n_input_rows = n_input_rows_loaded,
  n_rows_after_dpbn = n_rows_after_dpbn,
  n_overlap_ids = n_overlap_ids,
  n_rows_after_dpbn_and_overlap = nrow(dt_multistudies),
  mh_start_month_year = date_start_month_year,
  mh_end_month_year = date_end_month_year,
  mh_range_month_year = date_range_month_year,
  sdoh_start_month_year = sdoh_start_month_year,
  sdoh_end_month_year = sdoh_end_month_year,
  sdoh_range_month_year = sdoh_range_month_year,
  sdoh_start_year = sdoh_start_year,
  sdoh_end_year = sdoh_end_year,
  sdoh_range_year = sdoh_range_year,
  n_analysis_rows = nrow(analysis_data),
  n_participant_snapshot = nrow(participant_snapshot),
  min_cell_n = min_cell_n,
  n_filter_dimensions_configured = nrow(filter_schema),
  n_filter_dimensions_available = sum(filter_schema$available),
  n_filter_cube_rows = nrow(filter_cube),
  n_trajectory_dimensions_available = nrow(trajectory_schema),
  n_trajectory_trend_rows = nrow(trajectory_trend_summary),
  n_trajectory_pairwise_rows = nrow(trajectory_pairwise_summary)
)

saveRDS(kpi, file.path(objects_dir, "dashboard_kpi.rds"))
saveRDS(study_summary, file.path(objects_dir, "dashboard_study_summary.rds"))
saveRDS(demog_summary, file.path(objects_dir, "dashboard_demog_summary.rds"))
saveRDS(year_summary, file.path(objects_dir, "dashboard_year_summary.rds"))
saveRDS(mediator_summary, file.path(objects_dir, "dashboard_mediator_summary.rds"))
saveRDS(filter_schema, file.path(objects_dir, "dashboard_filter_schema.rds"))
saveRDS(filter_cube, file.path(objects_dir, "dashboard_filter_cube.rds"))
saveRDS(trajectory_trend_summary, file.path(objects_dir, "dashboard_trajectory_trend_summary.rds"))
saveRDS(trajectory_pairwise_summary, file.path(objects_dir, "dashboard_trajectory_pairwise_summary.rds"))
saveRDS(dashboard_build_info, file.path(objects_dir, "dashboard_build_info.rds"))

jsonlite::write_json(
  kpi,
  file.path(objects_dir, "dashboard_kpi.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  na = "null"
)
jsonlite::write_json(
  as.data.frame(study_summary),
  file.path(objects_dir, "dashboard_study_summary.json"),
  dataframe = "rows",
  pretty = TRUE,
  na = "null"
)
jsonlite::write_json(
  as.data.frame(demog_summary),
  file.path(objects_dir, "dashboard_demog_summary.json"),
  dataframe = "rows",
  pretty = TRUE,
  na = "null"
)
jsonlite::write_json(
  as.data.frame(year_summary),
  file.path(objects_dir, "dashboard_year_summary.json"),
  dataframe = "rows",
  pretty = TRUE,
  na = "null"
)
jsonlite::write_json(
  as.data.frame(mediator_summary),
  file.path(objects_dir, "dashboard_mediator_summary.json"),
  dataframe = "rows",
  pretty = TRUE,
  na = "null"
)
jsonlite::write_json(
  as.data.frame(filter_schema),
  file.path(objects_dir, "dashboard_filter_schema.json"),
  dataframe = "rows",
  pretty = TRUE,
  na = "null"
)
jsonlite::write_json(
  as.data.frame(filter_cube),
  file.path(objects_dir, "dashboard_filter_cube.json"),
  dataframe = "rows",
  pretty = TRUE,
  na = "null"
)
jsonlite::write_json(
  as.data.frame(trajectory_trend_summary),
  file.path(objects_dir, "dashboard_trajectory_trend_summary.json"),
  dataframe = "rows",
  pretty = TRUE,
  na = "null"
)
jsonlite::write_json(
  as.data.frame(trajectory_pairwise_summary),
  file.path(objects_dir, "dashboard_trajectory_pairwise_summary.json"),
  dataframe = "rows",
  pretty = TRUE,
  na = "null"
)
jsonlite::write_json(
  dashboard_build_info,
  file.path(objects_dir, "dashboard_build_info.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  na = "null"
)

save(kpi, file = file.path(objects_dir, "dashboard_kpi.RData"))
save(study_summary, file = file.path(objects_dir, "dashboard_study_summary.RData"))
save(demog_summary, file = file.path(objects_dir, "dashboard_demog_summary.RData"))
save(year_summary, file = file.path(objects_dir, "dashboard_year_summary.RData"))
save(mediator_summary, file = file.path(objects_dir, "dashboard_mediator_summary.RData"))
save(filter_schema, file = file.path(objects_dir, "dashboard_filter_schema.RData"))
save(filter_cube, file = file.path(objects_dir, "dashboard_filter_cube.RData"))
save(trajectory_trend_summary, file = file.path(objects_dir, "dashboard_trajectory_trend_summary.RData"))
save(trajectory_pairwise_summary, file = file.path(objects_dir, "dashboard_trajectory_pairwise_summary.RData"))
save(dashboard_build_info, file = file.path(objects_dir, "dashboard_build_info.RData"))

cat("Created dashboard aggregated objects in:", objects_dir, "\n")
cat("Source - dt_multistudies:", data_source_multi, "\n")
cat("Source - dt_childhood_exposure:", data_source_child, "\n")
cat(
  "Input rows:",
  n_input_rows_loaded,
  "| After DPBN filter:",
  n_rows_after_dpbn,
  "| After overlap filter:",
  nrow(dt_multistudies),
  "\n"
)
cat("Mental health date range (Month-Year):", date_range_month_year, "\n")
cat("SDoH trajectory range (Month-Year):", sdoh_range_month_year, "\n")
cat("Participant snapshot rows:", nrow(participant_snapshot), "\n")
cat("Study summary rows:", nrow(study_summary), "\n")
cat("Demographic summary rows:", nrow(demog_summary), "\n")
cat("Year summary rows:", nrow(year_summary), "\n")
cat("Mediator summary rows:", nrow(mediator_summary), "\n")
cat("Filter dimensions configured:", nrow(filter_schema), "\n")
cat("Filter dimensions available:", sum(filter_schema$available), "\n")
cat("Filter cube rows:", nrow(filter_cube), "\n")
cat("Trajectory dimensions available:", nrow(trajectory_schema), "\n")
cat("Trajectory trend rows:", nrow(trajectory_trend_summary), "\n")
cat("Trajectory cross-class rows:", nrow(trajectory_pairwise_summary), "\n")
cat("Suppression threshold (min cell n):", min_cell_n, "\n")
