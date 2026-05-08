`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

.package_root <- function() {
  root_opt <- getOption("portableIRTEngine.package_root", "")
  if (nzchar(root_opt) && dir.exists(root_opt)) {
    return(normalizePath(root_opt, mustWork = TRUE))
  }

  if (file.exists("DESCRIPTION")) {
    return(normalizePath(".", mustWork = TRUE))
  }

  stop(
    "Could not locate the portableIRTEngine package root. ",
    "Set options(portableIRTEngine.package_root = '/path/to/portable_irt_engine').",
    call. = FALSE
  )
}

.engines_dir <- function() {
  installed_dir <- system.file("engines", package = "portableIRTEngine")
  if (nzchar(installed_dir) && dir.exists(installed_dir)) {
    return(installed_dir)
  }

  source_dir <- file.path(.package_root(), "inst", "engines")
  if (dir.exists(source_dir)) {
    return(source_dir)
  }

  stop("Could not locate the packaged engines directory.", call. = FALSE)
}

.extdata_dir <- function() {
  installed_dir <- system.file("extdata", package = "portableIRTEngine")
  if (nzchar(installed_dir) && dir.exists(installed_dir)) {
    return(installed_dir)
  }

  source_dir <- file.path(.package_root(), "inst", "extdata")
  if (dir.exists(source_dir)) {
    return(source_dir)
  }

  stop("Could not locate the packaged extdata directory.", call. = FALSE)
}

.standardize_sex <- function(x) {
  if (is.null(x)) {
    return(x)
  }

  x_chr <- trimws(as.character(x))
  x_low <- tolower(x_chr)

  mapped <- ifelse(
    x_low %in% c("male", "m", "boy"),
    "Male",
    ifelse(
      x_low %in% c("female", "f", "girl"),
      "Female",
      ifelse(x_chr %in% c("", "NA", "Na", "na"), NA_character_, x_chr)
    )
  )

  mapped
}

.standardize_age_group <- function(x) {
  if (is.null(x)) {
    return(x)
  }

  x_chr <- trimws(as.character(x))
  x_low <- tolower(x_chr)

  mapped <- ifelse(
    x_low %in% c("age_17_19", "17-19", "17_19", "17 to 19"),
    "Age_17_19",
    ifelse(
      x_low %in% c("age_20_24", "20-24", "20_24", "20 to 24"),
      "Age_20_24",
      ifelse(x_chr %in% c("", "NA", "Na", "na"), NA_character_, x_chr)
    )
  )

  mapped
}

.ensure_groups <- function(df) {
  if ("SEX" %in% names(df)) {
    df$SEX <- .standardize_sex(df$SEX)
  }

  if ("AGEGRP" %in% names(df)) {
    df$AGEGRP <- .standardize_age_group(df$AGEGRP)
  }

  if (!"AGEGRP" %in% names(df) && "AGE" %in% names(df)) {
    age_num <- suppressWarnings(as.numeric(df$AGE))
    df$AGEGRP <- ifelse(is.na(age_num), NA_character_, ifelse(age_num < 20, "Age_17_19", "Age_20_24"))
  }

  if (!"SEX_AGEGRP" %in% names(df)) {
    sx <- if ("SEX" %in% names(df)) as.character(df$SEX) else NA_character_
    ag <- if ("AGEGRP" %in% names(df)) as.character(df$AGEGRP) else NA_character_
    df$SEX_AGEGRP <- ifelse(is.na(sx) | is.na(ag), NA_character_, paste0(sx, "__", ag))
  }

  df
}
