# ==============================================================================
# UTILITY: LCGA Reporting Core
# PURPOSE: Shared helpers for trajectory reporting scripts.
# ==============================================================================

library(dplyr)
library(glue)
library(gt)
library(stringr)
library(tibble)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

period_label_from_time_label <- function(time_label) {
  dplyr::case_when(
    identical(time_label, "time_0_5") ~ "Early Childhood (0–5 Years)",
    identical(time_label, "time_6_12") ~ "Middle–Late Childhood (6–12 Years)",
    grepl("^time_\\d+_\\d+$", time_label) ~ {
      ages <- as.integer(str_match(time_label, "^time_(\\d+)_(\\d+)$")[, 2:3])
      glue("Childhood Window ({ages[1]}–{ages[2]} Years)")
    },
    TRUE ~ time_label
  )
}

ages_from_time_label <- function(time_label) {
  age_match <- str_match(time_label, "^time_(\\d+)_(\\d+)$")
  if (any(is.na(age_match[1, 2:3]))) {
    stop(glue("Unable to parse ages from time label '{time_label}'."))
  }

  seq.int(as.integer(age_match[1, 2]), as.integer(age_match[1, 3]))
}

resolve_report_time_labels <- function(time_labels = NULL, config = cfg) {
  if (!is.null(time_labels) && length(time_labels) > 0) {
    return(unique(as.character(time_labels)))
  }

  unique(vapply(
    config$time_ranges,
    function(age_range) {
      ages <- as.integer(age_range)
      glue("time_{min(ages)}_{max(ages)}")
    },
    FUN.VALUE = character(1)
  ))
}

filter_lcga_specs <- function(specs, spec_names = NULL) {
  if (is.null(spec_names) || length(spec_names) == 0) {
    return(specs)
  }

  wanted <- unique(as.character(spec_names))
  Filter(function(spec) spec$name %in% wanted, specs)
}

normalise_fit_summary <- function(fit_df) {
  if (is.null(fit_df) || nrow(fit_df) == 0) {
    return(tibble())
  }

  fit_df <- as_tibble(fit_df)

  if (!"Classes" %in% names(fit_df) && "NLatentClasses" %in% names(fit_df)) {
    fit_df$Classes <- fit_df$NLatentClasses
  }
  if (!"BLRT_p" %in% names(fit_df) && "BLRT_PValue" %in% names(fit_df)) {
    fit_df$BLRT_p <- fit_df$BLRT_PValue
  }
  if (!"BLRT_p" %in% names(fit_df) && "T14_BLRT_PValue" %in% names(fit_df)) {
    fit_df$BLRT_p <- fit_df$T14_BLRT_PValue
  }
  if (!"VLMR_p" %in% names(fit_df) && "T11_VLMR_PValue" %in% names(fit_df)) {
    fit_df$VLMR_p <- fit_df$T11_VLMR_PValue
  }
  if (!"model_id" %in% names(fit_df) && "filename" %in% names(fit_df)) {
    fit_df$model_id <- fit_df$filename
  }

  for (col in c("Classes", "BIC", "Entropy", "BLRT_p", "VLMR_p", "MinClassPct", "model_id")) {
    if (!col %in% names(fit_df)) {
      fit_df[[col]] <- NA
    }
  }

  fit_df %>%
    mutate(
      Classes = suppressWarnings(as.integer(Classes)),
      BIC = suppressWarnings(as.numeric(BIC)),
      Entropy = suppressWarnings(as.numeric(Entropy)),
      BLRT_p = suppressWarnings(as.numeric(BLRT_p)),
      VLMR_p = suppressWarnings(as.numeric(VLMR_p)),
      MinClassPct = suppressWarnings(as.numeric(MinClassPct))
    ) %>%
    arrange(Classes)
}

resolve_manual_k <- function(spec_name, time_label = NULL, config = cfg) {
  if (!isTRUE(config$use_manual_overrides)) {
    return(NULL)
  }

  key_window <- if (!is.null(time_label)) paste0(spec_name, "_", time_label) else NULL

  config$manual_class_overrides[[key_window]] %||%
    config$manual_class_overrides[[spec_name]] %||%
    NULL
}

select_best_solution <- function(fit_df, spec_name, time_label = NULL, config = cfg) {
  manual_k <- resolve_manual_k(spec_name, time_label = time_label, config = config)

  if (!is.null(manual_k)) {
    matched <- fit_df %>% filter(Classes == as.integer(manual_k))
    if (nrow(matched) == 0) {
      stop(glue("Manual class override for {spec_name} ({manual_k}) was not found in fit summary."))
    }

    row <- matched %>% slice(1)
    return(list(
      model_row = row,
      model_name = if ("model_name" %in% names(row)) row$model_name[[1]] else NA,
      model_id = if ("model_id" %in% names(row)) row$model_id[[1]] else NA,
      Classes = row$Classes[[1]],
      BIC = row$BIC[[1]],
      Entropy = row$Entropy[[1]],
      test = config$selection_test,
      p_col = if (identical(config$selection_test, "BLRT")) "BLRT_p" else "VLMR_p",
      min_class_pct = row$MinClassPct[[1]],
      bic_close_delta = config$bic_close_delta,
      alpha = config$selection_alpha,
      filtered_table = fit_df
    ))
  }

  select_lcga_model(
    fit_df = fit_df,
    preset = config$selection_preset,
    test = config$selection_test,
    alpha = config$selection_alpha,
    min_class_pct = config$min_class_pct,
    bic_close_delta = config$bic_close_delta
  )
}

extract_lcga_class_assignments <- function(model_obj, class_labels = NULL, id_var = "USUBJID") {
  if (is.null(model_obj$savedata)) {
    return(NULL)
  }

  savedata <- tibble::as_tibble(model_obj$savedata)
  class_col <- grep("^C$|^MLCC$", names(savedata), value = TRUE)[1]
  prob_cols <- sort(grep("^CPROB\\d+$", names(savedata), value = TRUE))

  if (is.na(class_col) || length(prob_cols) == 0) {
    return(NULL)
  }

  n_classes <- length(prob_cols)
  label_map <- class_labels %||% setNames(paste("Class", seq_len(n_classes)), as.character(seq_len(n_classes)))
  if (is.null(names(label_map))) {
    names(label_map) <- as.character(seq_along(label_map))
  }

  savedata %>%
    mutate(Assigned_Class = as.integer(.data[[class_col]])) %>%
    select(any_of(id_var), Assigned_Class, all_of(prob_cols)) %>%
    mutate(
      Assigned_Class_Label = factor(
        Assigned_Class,
        levels = seq_len(n_classes),
        labels = unname(label_map[as.character(seq_len(n_classes))])
      )
    )
}

summarise_lcga_avepp <- function(class_assign_data) {
  if (is.null(class_assign_data)) {
    return(list(by_class = tibble(), mean_diag = NA_real_, min_diag = NA_real_))
  }

  prob_cols <- sort(grep("^CPROB\\d+$", names(class_assign_data), value = TRUE))
  if (length(prob_cols) == 0 || !"Assigned_Class" %in% names(class_assign_data)) {
    return(list(by_class = tibble(), mean_diag = NA_real_, min_diag = NA_real_))
  }

  by_class <- class_assign_data %>%
    group_by(Assigned_Class, Assigned_Class_Label) %>%
    group_modify(function(.x, .y) {
      prob_col <- paste0("CPROB", .y$Assigned_Class[[1]])
      tibble(
        Count = nrow(.x),
        DiagonalAvePP = if (prob_col %in% names(.x)) mean(.x[[prob_col]], na.rm = TRUE) else NA_real_
      )
    }) %>%
    ungroup()

  list(
    by_class = by_class,
    mean_diag = mean(by_class$DiagonalAvePP, na.rm = TRUE),
    min_diag = min(by_class$DiagonalAvePP, na.rm = TRUE)
  )
}

trajectory_y_label <- function(meta) {
  if ((meta$type %||% "") == "continuous") {
    return("Value")
  }

  target_label <- meta$plot_target %||% "highest"
  lvls <- meta$transform_levels %||% meta$levels

  if (!is.null(lvls) && length(lvls) > 0) {
    if (identical(target_label, "lowest")) {
      target_label <- lvls[1]
    }
    if (identical(target_label, "highest")) {
      target_label <- tail(lvls, 1)
    }
  }

  glue("Probability ({target_label})")
}

save_gt_html <- function(gt_tbl, path) {
  if (is.null(gt_tbl)) {
    return(invisible(NULL))
  }
  gt::gtsave(gt_tbl, path)
  invisible(path)
}

apply_class_color_scale <- function(plot_obj, class_color_map = NULL, aes_name = c("color", "fill")) {
  if (is.null(plot_obj) || is.null(class_color_map) || length(class_color_map) == 0) {
    return(plot_obj)
  }

  aes_name <- match.arg(aes_name)
  if (identical(aes_name, "fill")) {
    return(plot_obj + ggplot2::scale_fill_manual(values = class_color_map, drop = FALSE))
  }

  plot_obj + ggplot2::scale_color_manual(values = class_color_map, drop = FALSE)
}
