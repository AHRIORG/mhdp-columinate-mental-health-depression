# Public-safe helpers for the pathway screening release bundle.

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
})

find_pathway_screening_root <- function(start = getwd()) {
  is_bundle_root <- function(path) {
    dir.exists(file.path(path, "programs", "analysis")) &&
      dir.exists(file.path(path, "programs", "utils"))
  }

  build_candidates <- function(path) {
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
    unique(c(parents, file.path(parents, "scripts", "pathway_screening")))
  }

  candidates <- unique(c(
    Sys.getenv("COLU_PATHWAY_SCREENING_ROOT", unset = ""),
    start,
    getwd()
  ))

  for (candidate in unique(unlist(lapply(candidates[nzchar(candidates)], build_candidates)))) {
    if (is_bundle_root(candidate)) {
      resolved <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
      Sys.setenv(COLU_PATHWAY_SCREENING_ROOT = resolved)
      return(resolved)
    }
  }

  stop("Unable to locate the pathway screening bundle root.")
}

pathway_screening_path <- function(...) {
  file.path(find_pathway_screening_root(), ...)
}

pathway_repo_root <- function() {
  normalizePath(file.path(find_pathway_screening_root(), "..", ".."), winslash = "/", mustWork = TRUE)
}

default_pathway_release_dir <- function() {
  repo_root <- pathway_repo_root()
  private_release <- file.path(repo_root, "results", "release", "pathway_screening")
  public_release <- file.path(repo_root, "results", "pathway_screening")
  if (dir.exists(file.path(repo_root, "results", "release"))) {
    private_release
  } else {
    public_release
  }
}

clean_public_text <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "NA" | x == "NaN"] <- "-"
  user_path_marker <- paste0("/", "Users", "/")
  cloud_path_marker <- paste0("One", "Drive")
  x <- str_replace_all(x, paste0(user_path_marker, "[^,;\\\"]+"), "restricted input")
  x <- str_replace_all(x, paste0(cloud_path_marker, "[^,;\\\"]*"), "restricted input")
  private_repo_marker <- paste0("private", "/", "CO-LUMINATE")
  private_parent_marker <- paste0("org_repos", "/", "private")
  objective_markers <- c(
    paste0("OBJ", "00-Datasets"),
    paste0("OBJ", "01-Clustering"),
    paste0("OBJ", "02-Interactive Pathways"),
    paste0("OBJ", "03-Psychometric"),
    paste0("OBJ", "04-Causal Mediation")
  )
  objective_replacements <- c(
    "data-management workflow",
    "class-registry workflow",
    "pathway-inventory workflow",
    "depression-harmonization workflow",
    "causal-mediation workflow"
  )
  project_markers <- c(
    paste0("PRJ", "LCGA"),
    paste0("PRJ", "PSM"),
    paste0("PRJ", "CMA")
  )
  project_replacements <- c(
    "trajectory modelling",
    "pathway screening",
    "causal mediation"
  )
  x <- str_replace_all(x, private_repo_marker, "private source workspace")
  x <- str_replace_all(x, private_parent_marker, "restricted workspace")
  for (i in seq_along(objective_markers)) {
    x <- str_replace_all(x, paste0(objective_markers[[i]], "[^,;\\\"]*"), objective_replacements[[i]])
  }
  for (i in seq_along(project_markers)) {
    x <- str_replace_all(x, project_markers[[i]], project_replacements[[i]])
  }
  x
}

drop_private_path_columns <- function(data) {
  keep <- vapply(names(data), function(name) {
    values <- as.character(data[[name]])
    private_repo_marker <- paste0("private", "/", "CO-LUMINATE")
    private_parent_marker <- paste0("org_repos", "/", "private")
    user_path_marker <- paste0("/", "Users", "/")
    cloud_path_marker <- paste0("One", "Drive")
    path_like_values <- any(str_detect(values, paste0(user_path_marker, "|", cloud_path_marker, "|", private_repo_marker, "|", private_parent_marker, "|\\.rds")), na.rm = TRUE)
    path_like_name <- str_detect(name, "(^|_)(result_path|source_path|output_path|workspace_root|script|inventory_source)$")
    !(path_like_values || path_like_name)
  }, logical(1))
  data[, keep, drop = FALSE]
}

write_public_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data[] <- lapply(data, clean_public_text)
  readr::write_csv(data, path, na = "-")
  invisible(path)
}
