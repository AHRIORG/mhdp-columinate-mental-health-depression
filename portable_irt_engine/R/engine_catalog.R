list_engines <- function() {
  catalog_path <- file.path(.engines_dir(), "engine_catalog.csv")
  out <- utils::read.csv(catalog_path, stringsAsFactors = FALSE)
  out$default <- as.logical(out$default)
  if ("available_in_app" %in% names(out)) {
    out$available_in_app <- as.logical(out$available_in_app)
  }
  if ("scoring_ready" %in% names(out)) {
    out$scoring_ready <- as.logical(out$scoring_ready)
  }
  out
}

default_engine_id <- function() {
  catalog <- list_engines()
  default_rows <- catalog[catalog$default %in% TRUE, , drop = FALSE]

  if (nrow(default_rows) == 0) {
    stop("No default engine is marked in the engine catalog.", call. = FALSE)
  }

  default_rows$engine_id[1]
}

engine_metadata <- function(engine_id = NULL) {
  catalog <- list_engines()
  engine_id <- engine_id %||% default_engine_id()
  meta <- catalog[catalog$engine_id == engine_id, , drop = FALSE]

  if (nrow(meta) == 0) {
    stop("Unknown engine_id: ", engine_id, call. = FALSE)
  }

  meta
}

engine_provenance <- function(engine_id = NULL) {
  meta <- engine_metadata(engine_id)
  meta_file <- meta$meta_file_name[1]

  if (is.na(meta_file) || !nzchar(meta_file)) {
    return(NULL)
  }

  readRDS(file.path(.engines_dir(), meta_file))
}

engine_input_schema <- function(engine_id = NULL) {
  engine <- load_engine(engine_id)

  data.frame(
    field = c(engine$items, "SEX", "AGE", "AGEGRP"),
    required = c(rep(TRUE, length(engine$items)), FALSE, FALSE, FALSE),
    expected_values = c(
      rep("0/1 or Yes/No", length(engine$items)),
      "Male/Female (optional)",
      "Numeric age in years (optional)",
      "Age_17_19 / Age_20_24 (optional)"
    ),
    stringsAsFactors = FALSE
  )
}

load_engine <- function(engine_id = NULL) {
  meta <- engine_metadata(engine_id)
  engine <- readRDS(file.path(.engines_dir(), meta$file_name[1]))
  attr(engine, "engine_id") <- meta$engine_id[1]
  attr(engine, "engine_label") <- meta$label[1]
  attr(engine, "engine_metadata") <- meta
  engine
}
