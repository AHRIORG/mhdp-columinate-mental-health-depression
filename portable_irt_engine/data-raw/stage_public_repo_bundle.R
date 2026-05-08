stage_public_repo_bundle <- function(
  source_root = ".",
  public_repo_root = "../../../public/mhdp-columinate-mental-health-depression",
  target_subdir = "portable_irt_engine"
) {
  source_root <- normalizePath(source_root, mustWork = TRUE)
  public_repo_root <- normalizePath(public_repo_root, mustWork = TRUE)
  target_root <- file.path(public_repo_root, target_subdir)

  dir.create(target_root, recursive = TRUE, showWarnings = FALSE)

  include_paths <- c(
    ".Rbuildignore",
    ".dockerignore",
    ".gitignore",
    "DESCRIPTION",
    "DEPLOYMENT.md",
    "Dockerfile",
    "LICENSE",
    "NAMESPACE",
    "README.md",
    "PUBLIC_RELEASE_MANIFEST.md",
    "R",
    "deploy",
    "inst",
    "shiny-app",
    "tests",
    file.path("data-raw", "sanitize_engine_metadata.R"),
    file.path("data-raw", "stage_public_repo_bundle.R"),
    file.path("data-raw", "rebuild_engine_bundle.R")
  )

  copy_path <- function(rel_path) {
    src <- file.path(source_root, rel_path)
    dst <- file.path(target_root, rel_path)

    if (!file.exists(src)) {
      stop("Missing release path: ", src, call. = FALSE)
    }

    if (dir.exists(src)) {
      dir.create(dst, recursive = TRUE, showWarnings = FALSE)
      files <- list.files(src, all.files = FALSE, no.. = TRUE, full.names = TRUE, recursive = TRUE)
      if (length(files) == 0) {
        return(invisible(NULL))
      }
      rel_files <- substring(files, nchar(src) + 2L)
      for (i in seq_along(files)) {
        out_file <- file.path(dst, rel_files[i])
        dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
        ok <- file.copy(files[i], out_file, overwrite = TRUE, copy.date = TRUE)
        if (!ok) stop("Failed to copy file: ", files[i], call. = FALSE)
      }
    } else {
      dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
      ok <- file.copy(src, dst, overwrite = TRUE, copy.date = TRUE)
      if (!ok) stop("Failed to copy file: ", src, call. = FALSE)
    }
  }

  invisible(lapply(include_paths, copy_path))

  message("Staged public bundle at: ", target_root)
  invisible(target_root)
}

if (sys.nframe() == 0) {
  stage_public_repo_bundle()
}
