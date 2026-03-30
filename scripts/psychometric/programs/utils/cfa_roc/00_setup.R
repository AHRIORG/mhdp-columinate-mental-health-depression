############################################################
# SSQ / PHQ Monte Carlo CFA Validation
# 00_setup.R - Packages, data loading, configuration
############################################################

load_required_packages <- function() {
  pkgs <- c(
    "here", "psych", "lavaan", "semTools", "dplyr",
    "tidyr", "pROC", "ggplot2", "ggsci", "gt"
  )
  for (p in pkgs) {
    if (!require(p, character.only = TRUE)) {
      install.packages(p)
      library(p, character.only = TRUE)
    }
  }
}

load_main_data <- function(path = NULL, object_name = NULL) {
  # 1) If a data.frame object already exists in the environment, use it
  if (!is.null(object_name) && exists(object_name, inherits = TRUE)) {
    dt_main <- get(object_name)
    if (!is.data.frame(dt_main)) {
      stop(paste("Object", object_name, "exists but is not a dataframe."))
    }
    message(paste("Using dataframe already in environment:", object_name))
    return(dt_main)
  }
  
  # 2) Otherwise, a path must be provided
  if (is.null(path)) {
    stop("Either 'path' must be provided or 'object_name' must refer to a dataframe in the environment.")
  }
  
  if (!file.exists(path)) {
    stop(paste("Data file not found at:", path))
  }
  
  loaded_objs <- load(file = path)
  
  # If an explicit object_name was requested and is in the RData, use that
  if (!is.null(object_name) && object_name %in% loaded_objs) {
    dt_main <- get(object_name)
    message(paste("Using dataframe from RData file:", object_name))
    
  } else if (exists(loaded_objs, inherits = FALSE)) {
    # As per your requested version: if RData only contains one object
    dt_main <- get(loaded_objs)
    
  } else {
    # Search for a data.frame among loaded objects
    found_df <- FALSE
    for (obj_name in loaded_objs) {
      obj <- get(obj_name)
      if (is.data.frame(obj)) {
        dt_main <- obj
        found_df <- TRUE
        message(paste("Using dataframe:", obj_name))
        break
      }
    }
    if (!found_df) {
      if (length(loaded_objs) > 0) {
        dt_main <- get(loaded_objs[1])
        warning(paste("No dataframe found. Using first object:", loaded_objs[1]))
      } else {
        stop("No objects found in RData file.")
      }
    }
  }
  
  if (is.null(dt_main) || !is.data.frame(dt_main)) {
    stop("Loaded object is not a valid dataframe.")
  }
  
  message(paste("Data loaded successfully. Rows:", nrow(dt_main), "Cols:", ncol(dt_main)))
  dt_main
}

create_item_sets <- function() {
  list(
    "SSQ-14 (Original)"                        = paste0("SSQ", sprintf("%02d", 1:14)),
    "SSQ-10 (Theoretical Depressive Symptoms)" = paste0("SSQ", sprintf("%02d", c(1:3, 8:14))),
    "SSQ-08 (Dixon Chibanda et.al)"           = paste0("SSQ", sprintf("%02d", c(1, 8:14))),
    "PHQ-09 (Original)"                       = paste0("PHQ9", sprintf("%02d", 1:9))
  )
}

analysis_config <- function(item_sets_list,
                            target_set_name       = "SSQ-10 (Theoretical Depressive Symptoms)",
                            dataset_label         = "Isisekelo Sempilo",
                            manual_method         = "Parallel",  # or NULL for auto
                            manual_n_factors      = 3,            # or NULL for auto
                            n_mc_samples          = 100,
                            min_items_per_factor  = 3,
                            manual_merge_pairs    = list(c(3, 1))) {   # F3 -> F1 by default
  if (!target_set_name %in% names(item_sets_list)) {
    stop(paste("target_set_name not found in item_sets_list:", target_set_name))
  }
  
  target_items <- item_sets_list[[target_set_name]]
  
  message(paste("Target Set:", target_set_name))
  message(paste("Items:", paste(target_items, collapse = ", ")))
  
  list(
    target_set_name       = target_set_name,
    target_items          = target_items,
    dataset_label         = dataset_label,
    manual_method         = manual_method,
    manual_n_factors      = manual_n_factors,
    n_mc_samples          = n_mc_samples,
    min_items_per_factor  = min_items_per_factor,
    manual_merge_pairs    = manual_merge_pairs
  )
}
