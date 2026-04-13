# ==============================================================================
# UTILITY: Canonical Class Registry
# PURPOSE: Provide a bundle-level source of truth for LCGA class labels,
#          archetypes, colors, and narrative metadata in the public trajectory bundle.
# NOTES:
#   - Source_Project and Source_File are public-safe provenance labels only.
#   - Any unpublished curation logic should be materialized into the registry
#     before handoff rather than referenced here by private project path/name.
# ==============================================================================

library(dplyr)
library(glue)
library(here)
library(readr)
library(tibble)

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

class_registry_schema <- c(
  "Model_Family",
  "Variable",
  "Period",
  "Classes",
  "LatentClass",
  "Assigned_Class_Label",
  "Archetype_Group",
  "Display_Label",
  "Narrative_Label",
  "Color_Name",
  "Color_Hex",
  "Rationale",
  "Legacy_Display_Label",
  "Legacy_Archetype_Group",
  "Legacy_Rationale",
  "Source_Project",
  "Source_File",
  "Status",
  "Decision_Source"
)

cb_dict <- tribble(
  ~name, ~hex,
  "gray", "#BBBBBB",
  "indigo", "#332288",
  "light_gray", "#DDDDBD",
  "red", "#D73027",
  "bluish_green", "#009E73",
  "teal", "#44AA99",
  "green", "#228833",
  "orange", "#E69F00",
  "vermillion", "#D55E00",
  "pink", "#EE6677",
  "purple", "#AA3377",
  "lavender", "#CC79A7",
  "wine", "#882255",
  "sky", "#88CCEE",
  "blue", "#4477AA",
  "yellow", "#CCBB44",
  "missing", "#000000"
) %>%
  mutate(hex = toupper(hex))

get_default_period_label <- function(model_family = NULL) {
  if (is.null(model_family)) {
    return(NA_character_)
  }

  switch(
    toupper(model_family),
    "LCA" = "Cross-sectional",
    NA_character_
  )
}

get_colorblind_palette <- function(include_missing = FALSE) {
  palette_tbl <- cb_dict
  if (!isTRUE(include_missing)) {
    palette_tbl <- palette_tbl %>% filter(name != "missing")
  }

  palette_tbl %>%
    distinct(hex, .keep_all = TRUE) %>%
    mutate(name = make.unique(name, sep = "_"))
}

default_class_palette <- function(n_classes) {
  if (is.null(n_classes) || n_classes < 1) {
    return(tibble(
      LatentClass = integer(),
      Color_Name = character(),
      Color_Hex = character()
    ))
  }

  palette_tbl <- get_colorblind_palette(include_missing = FALSE)
  palette_idx <- rep(seq_len(nrow(palette_tbl)), length.out = n_classes)

  tibble(
    LatentClass = seq_len(n_classes),
    Color_Name = palette_tbl$name[palette_idx],
    Color_Hex = palette_tbl$hex[palette_idx]
  )
}

default_class_labels <- function(n_classes, default_prefix = "Class") {
  if (is.null(n_classes) || n_classes < 1) {
    return(tibble(
      LatentClass = integer(),
      Assigned_Class_Label = character(),
      Archetype_Group = character(),
      Display_Label = character(),
      Narrative_Label = character(),
      Rationale = character(),
      Legacy_Display_Label = character(),
      Legacy_Archetype_Group = character(),
      Legacy_Rationale = character(),
      Source_Project = character(),
      Source_File = character(),
      Status = character(),
      Decision_Source = character()
    ))
  }

  labels <- paste(default_prefix, seq_len(n_classes))

  tibble(
    LatentClass = seq_len(n_classes),
    Assigned_Class_Label = labels,
    Archetype_Group = labels,
    Display_Label = labels,
    Narrative_Label = labels,
    Rationale = NA_character_,
    Legacy_Display_Label = NA_character_,
    Legacy_Archetype_Group = NA_character_,
    Legacy_Rationale = NA_character_,
    Source_Project = NA_character_,
    Source_File = NA_character_,
    Status = "default",
    Decision_Source = "default"
  )
}

empty_class_registry <- function() {
  out <- tibble()
  for (col in class_registry_schema) {
    out[[col]] <- vector("character", length = 0)
  }
  out$Classes <- integer()
  out$LatentClass <- integer()
  out
}

get_class_registry_path <- function(path = NULL) {
  path %||% trajectory_path("data", "metadata", "class_registry.csv")
}

write_empty_class_registry <- function(path = get_class_registry_path(), overwrite = FALSE) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop(glue("Registry already exists at {path}. Set overwrite = TRUE to replace it."))
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(empty_class_registry(), path, na = "")
  invisible(path)
}

normalise_class_registry <- function(registry) {
  registry <- as_tibble(registry)

  for (col in class_registry_schema) {
    if (!col %in% names(registry)) {
      registry[[col]] <- NA
    }
  }

  registry %>%
    mutate(
      Model_Family = ifelse(is.na(Model_Family), NA_character_, toupper(as.character(Model_Family))),
      Variable = as.character(Variable),
      Period = as.character(Period),
      Classes = suppressWarnings(as.integer(Classes)),
      LatentClass = suppressWarnings(as.integer(LatentClass)),
      Assigned_Class_Label = as.character(Assigned_Class_Label),
      Archetype_Group = as.character(Archetype_Group),
      Display_Label = as.character(Display_Label),
      Narrative_Label = as.character(Narrative_Label),
      Color_Name = as.character(Color_Name),
      Color_Hex = toupper(as.character(Color_Hex)),
      Rationale = as.character(Rationale),
      Legacy_Display_Label = as.character(Legacy_Display_Label),
      Legacy_Archetype_Group = as.character(Legacy_Archetype_Group),
      Legacy_Rationale = as.character(Legacy_Rationale),
      Source_Project = as.character(Source_Project),
      Source_File = as.character(Source_File),
      Status = as.character(Status),
      Decision_Source = as.character(Decision_Source)
    ) %>%
    distinct(Model_Family, Variable, Period, Classes, LatentClass, .keep_all = TRUE) %>%
    arrange(Model_Family, Variable, Period, Classes, LatentClass)
}

load_class_registry <- function(path = get_class_registry_path(), strict = FALSE) {
  if (!file.exists(path)) {
    if (isTRUE(strict)) {
      stop(glue("Class registry not found at {path}"))
    }
    return(empty_class_registry())
  }

  readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_character())) %>%
    normalise_class_registry()
}

validate_class_registry <- function(registry = NULL, path = get_class_registry_path()) {
  registry <- registry %||% load_class_registry(path = path, strict = FALSE)
  registry <- normalise_class_registry(registry)

  issues <- list()

  dupes <- registry %>%
    count(Model_Family, Variable, Period, Classes, LatentClass, name = "n") %>%
    filter(n > 1) %>%
    mutate(issue = "duplicate_key")

  if (nrow(dupes) > 0) {
    issues[["duplicate_key"]] <- dupes
  }

  missing_display <- registry %>%
    filter(
      !is.na(Model_Family),
      !is.na(Variable),
      !is.na(Classes),
      !is.na(LatentClass),
      is.na(Display_Label) | !nzchar(Display_Label)
    ) %>%
    mutate(issue = "missing_display_label")

  if (nrow(missing_display) > 0) {
    issues[["missing_display_label"]] <- missing_display
  }

  missing_colors <- registry %>%
    filter(
      !is.na(Model_Family),
      !is.na(Variable),
      !is.na(Classes),
      !is.na(LatentClass),
      is.na(Color_Hex) | !nzchar(Color_Hex)
    ) %>%
    mutate(issue = "missing_color")

  if (nrow(missing_colors) > 0) {
    issues[["missing_color"]] <- missing_colors
  }

  if (length(issues) == 0) {
    return(tibble())
  }

  bind_rows(issues)
}

resolve_class_metadata <- function(
  variable,
  classes,
  period = NULL,
  model_family = "LCGA",
  registry = NULL,
  registry_path = get_class_registry_path(),
  default_prefix = "Class"
) {
  if (is.null(classes) || is.na(classes) || classes < 1) {
    return(tibble(
      Model_Family = character(),
      Variable = character(),
      Period = character(),
      Classes = integer(),
      LatentClass = integer(),
      Assigned_Class_Label = character(),
      Archetype_Group = character(),
      Display_Label = character(),
      Narrative_Label = character(),
      Color_Name = character(),
      Color_Hex = character(),
      Rationale = character(),
      Legacy_Display_Label = character(),
      Legacy_Archetype_Group = character(),
      Legacy_Rationale = character(),
      Source_Project = character(),
      Source_File = character(),
      Status = character(),
      Decision_Source = character()
    ))
  }

  registry_norm <- normalise_class_registry(registry %||% load_class_registry(path = registry_path, strict = FALSE))
  model_family_norm <- toupper(as.character(model_family %||% NA_character_))
  period_norm <- period %||% get_default_period_label(model_family_norm)

  matched <- registry_norm %>%
    filter(
      Variable == variable,
      Classes == as.integer(classes),
      if (!is.na(model_family_norm)) Model_Family == model_family_norm else TRUE,
      if (!is.na(period_norm)) Period == period_norm else TRUE
    ) %>%
    select(-any_of(c("Variable", "Classes", "Period", "Model_Family")))

  defaults <- default_class_labels(classes, default_prefix = default_prefix) %>%
    left_join(default_class_palette(classes), by = "LatentClass")

  defaults %>%
    left_join(matched, by = "LatentClass", suffix = c("_default", "")) %>%
    transmute(
      Model_Family = model_family_norm,
      Variable = variable,
      Period = period_norm,
      Classes = as.integer(classes),
      LatentClass = LatentClass,
      Assigned_Class_Label = coalesce(Assigned_Class_Label, Assigned_Class_Label_default),
      Archetype_Group = coalesce(Archetype_Group, Archetype_Group_default),
      Display_Label = coalesce(Display_Label, Display_Label_default, Assigned_Class_Label, Archetype_Group),
      Narrative_Label = coalesce(Narrative_Label, Narrative_Label_default, Display_Label, Assigned_Class_Label),
      Color_Name = coalesce(Color_Name, Color_Name_default),
      Color_Hex = coalesce(Color_Hex, Color_Hex_default),
      Rationale = coalesce(Rationale, Rationale_default),
      Legacy_Display_Label = coalesce(Legacy_Display_Label, Legacy_Display_Label_default),
      Legacy_Archetype_Group = coalesce(Legacy_Archetype_Group, Legacy_Archetype_Group_default),
      Legacy_Rationale = coalesce(Legacy_Rationale, Legacy_Rationale_default),
      Source_Project = coalesce(Source_Project, Source_Project_default),
      Source_File = coalesce(Source_File, Source_File_default),
      Status = coalesce(Status, Status_default),
      Decision_Source = coalesce(Decision_Source, Decision_Source_default)
    ) %>%
    arrange(LatentClass)
}

get_class_label_map <- function(
  variable,
  classes,
  period = NULL,
  model_family = "LCGA",
  registry = NULL,
  registry_path = get_class_registry_path(),
  default_prefix = "Class",
  label_col = c("Display_Label", "Assigned_Class_Label", "Archetype_Group")
) {
  label_col <- match.arg(label_col)

  class_meta <- resolve_class_metadata(
    variable = variable,
    classes = classes,
    period = period,
    model_family = model_family,
    registry = registry,
    registry_path = registry_path,
    default_prefix = default_prefix
  )

  setNames(class_meta[[label_col]], class_meta$LatentClass)
}

get_class_color_map <- function(
  variable,
  classes,
  period = NULL,
  model_family = "LCGA",
  registry = NULL,
  registry_path = get_class_registry_path(),
  default_prefix = "Class",
  key = c("label", "class")
) {
  key <- match.arg(key)

  class_meta <- resolve_class_metadata(
    variable = variable,
    classes = classes,
    period = period,
    model_family = model_family,
    registry = registry,
    registry_path = registry_path,
    default_prefix = default_prefix
  )

  if (key == "label") {
    return(setNames(class_meta$Color_Hex, class_meta$Display_Label))
  }

  setNames(class_meta$Color_Hex, class_meta$LatentClass)
}
