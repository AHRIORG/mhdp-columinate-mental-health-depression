sanitize_engine_metadata <- function(engine_dir = file.path("inst", "engines")) {
  meta_files <- list.files(
    engine_dir,
    pattern = "_meta[.]rds$",
    full.names = TRUE
  )

  if (length(meta_files) == 0) {
    stop("No metadata RDS files found under: ", engine_dir, call. = FALSE)
  }

  for (meta_file in meta_files) {
    meta <- readRDS(meta_file)

    if (!is.null(meta$source_file) && length(meta$source_file) == 1) {
      meta$source_file <- basename(meta$source_file)
    }

    meta$source_scope <- "public_release_sanitized"
    meta$source_note <- paste(
      "Absolute local source paths were removed during public-release packaging.",
      "The source_file field now stores only the original artifact basename."
    )

    saveRDS(meta, meta_file)
    message("Sanitized metadata: ", meta_file)
  }

  invisible(meta_files)
}

if (sys.nframe() == 0) {
  sanitize_engine_metadata()
}
