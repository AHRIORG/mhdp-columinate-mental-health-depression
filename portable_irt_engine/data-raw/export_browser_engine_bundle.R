suppressPackageStartupMessages({
  library(jsonlite)
})

source_root <- normalizePath(file.path(getwd(), "."), mustWork = TRUE)

read_catalog <- function(root) {
  out <- utils::read.csv(
    file.path(root, "inst", "engines", "engine_catalog.csv"),
    stringsAsFactors = FALSE
  )
  out$default <- as.logical(out$default)
  out$available_in_app <- as.logical(out$available_in_app)
  out$scoring_ready <- as.logical(out$scoring_ready)
  out$n_factors <- as.integer(out$n_factors)
  out
}

load_engine_object <- function(root, file_name) {
  readRDS(file.path(root, "inst", "engines", file_name))
}

engine_label_browser <- function(row) {
  scope_map <- c(
    none = "Overall",
    SEX = "Sex",
    AGEGRP = "Age group",
    SEX_AGEGRP = "Sex + age group",
    hierarchical = "Hierarchical"
  )
  method_map <- c(
    sens_at_spec = "Sensitivity at target specificity",
    spec_at_sens = "Specificity at target sensitivity",
    youden = "Youden index threshold",
    model_tcc = "Model-based TCC threshold"
  )

  scope <- unname(scope_map[row$cutoff_apply_mode])
  method <- unname(method_map[row$cutoff_method])
  scope <- ifelse(is.na(scope), row$cutoff_apply_mode, scope)
  method <- ifelse(is.na(method), row$cutoff_method, method)
  paste("SSQ-10", paste0(row$n_factors, "-factor"), scope, method, sep = ", ")
}

compact_parameters <- function(pars_df) {
  items <- unique(as.character(pars_df$item))
  out <- lapply(items, function(item_name) {
    item_rows <- pars_df[pars_df$item == item_name, , drop = FALSE]
    a <- suppressWarnings(as.numeric(item_rows$value[item_rows$name == "a1"][1]))
    d <- suppressWarnings(as.numeric(item_rows$value[item_rows$name == "d1"][1]))
    list(a = a, d = d)
  })
  names(out) <- items
  out
}

export_browser_engine_bundle <- function(
  root = ".",
  out_dir = file.path("..", "website", "_tools", "irt-scoring-static", "engines")
) {
  root <- normalizePath(root, mustWork = TRUE)
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  catalog <- read_catalog(root)
  catalog <- catalog[
    catalog$available_in_app %in% TRUE &
      catalog$scoring_ready %in% TRUE &
      catalog$n_factors == 1L,
    ,
    drop = FALSE
  ]
  if (nrow(catalog) == 0) {
    stop("No production-safe 1-factor engines found for browser export.", call. = FALSE)
  }

  engines <- lapply(seq_len(nrow(catalog)), function(i) {
    row <- catalog[i, , drop = FALSE]
    engine <- load_engine_object(root, row$file_name[[1]])

    list(
      engine_id = row$engine_id[[1]],
      label = engine_label_browser(row),
      source_label = row$label[[1]],
      description = row$description[[1]],
      notes = row$notes[[1]],
      default = isTRUE(row$default[[1]]),
      n_factors = as.integer(row$n_factors[[1]]),
      model_type = as.character(engine$model_type %||% "graded"),
      items = unname(as.character(engine$items)),
      parameters = compact_parameters(engine$parameters),
      cutoff_apply_mode = as.character(engine$cutoff_apply_mode %||% "none"),
      cutoff_method = as.character(row$cutoff_method[[1]]),
      cutoff_table = engine$cutoff_table,
      theta_map = engine$phq_theta_map,
      metrics_applied = engine$cutoff_table
    )
  })

  default_engine <- catalog$engine_id[which(catalog$default)[1]]
  if (length(default_engine) == 0 || is.na(default_engine)) {
    default_engine <- catalog$engine_id[[1]]
  }

  bundle <- list(
    bundle_name = "portable-irt-browser-engines",
    generated_at = as.character(Sys.time()),
    default_engine_id = default_engine,
    n_engines = length(engines),
    engines = engines
  )

  jsonlite::write_json(
    bundle,
    path = file.path(out_dir, "production_engines.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    digits = NA,
    na = "null"
  )

  example_src <- file.path(root, "inst", "extdata", "example_batch.csv")
  example_dst <- file.path(out_dir, "example_batch.csv")
  if (file.exists(example_src)) {
    file.copy(example_src, example_dst, overwrite = TRUE)
  }

  message("Wrote browser engine bundle to: ", out_dir)
  invisible(file.path(out_dir, "production_engines.json"))
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

if (sys.nframe() == 0) {
  export_browser_engine_bundle()
}
