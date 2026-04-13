# ==============================================================================
# UTILITY: Joint LCGA Plotting Core
# PURPOSE: Joint-process growth-factor extraction and composite trajectory
#          helpers for the public trajectory bundle.
# ==============================================================================

library(dplyr)
library(glue)
library(ggplot2)
library(purrr)
library(scales)
library(tibble)

extract_joint_gf_means <- function(unstd, j, classes = NULL) {
  target_pattern <- glue("^[ISQ][._]?{j}$")

  gf <- unstd %>%
    filter(grepl("Means", paramHeader, ignore.case = TRUE), !is.na(LatentClass)) %>%
    mutate(paramU = toupper(trimws(param))) %>%
    filter(grepl(target_pattern, paramU)) %>%
    mutate(letter = substr(paramU, 1, 1)) %>%
    select(LatentClass, letter, est) %>%
    tidyr::pivot_wider(names_from = letter, values_from = est) %>%
    mutate(across(any_of(c("I", "S", "Q")), as.numeric))

  if (nrow(gf) == 0) {
    return(gf)
  }

  if (!("Q" %in% names(gf))) {
    gf$Q <- 0
  }

  if (!all(c("I", "S") %in% names(gf))) {
    warning(glue("Found some parameters for process {j} but missing I or S growth factors."))
    return(tibble())
  }

  if (!is.null(classes)) {
    gf <- gf %>% mutate(LatentClass = factor(LatentClass, levels = classes))
  } else {
    gf <- gf %>% mutate(LatentClass = factor(LatentClass))
  }

  gf
}

plot_joint_composite <- function(traj_out_a, method = "mean") {
  if (length(traj_out_a$data) == 0) {
    return(NULL)
  }

  smooth_all <- purrr::imap_dfr(traj_out_a$data, ~ .x$smooth %>% mutate(prefix = .y))
  points_all <- purrr::imap_dfr(traj_out_a$data, ~ .x$points %>% mutate(prefix = .y))

  if (nrow(smooth_all) == 0) {
    return(NULL)
  }

  calc_composite <- function(vals) {
    if (identical(method, "mean")) {
      return(mean(vals, na.rm = TRUE))
    }
    1 - prod(1 - vals, na.rm = TRUE)
  }

  comp_smooth <- smooth_all %>%
    group_by(LatentClass, Age, t) %>%
    summarise(val = calc_composite(val), .groups = "drop")

  comp_points <- points_all %>%
    group_by(LatentClass, Age, t) %>%
    summarise(val = calc_composite(val), .groups = "drop")

  title_str <- if (identical(method, "mean")) {
    "Composite: Mean Probability"
  } else {
    "Composite: Prob(Any Outcome)"
  }

  x_breaks <- sort(unique(comp_points$Age))
  if (length(x_breaks) == 0) {
    x_breaks <- waiver()
  }

  ggplot(comp_smooth, aes(x = Age, y = val, color = LatentClass, group = LatentClass)) +
    geom_line(linewidth = 1.5) +
    geom_point(data = comp_points, aes(x = Age, y = val, color = LatentClass), size = 3) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_x_continuous(breaks = x_breaks) +
    theme_minimal(base_size = 14) +
    labs(
      title = title_str,
      y = "Composite Probability",
      x = "Child Age (Years)",
      color = "Latent Class"
    ) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
}

plot_joint_trajectories <- function(
  best_model,
  outcome_meta,
  item_names_list,
  ages,
  time_scores,
  grid_step = 0.1,
  class_label_map = NULL
) {
  unstd <- best_model$parameters$unstandardized
  if (is.null(unstd)) {
    return(list(plots = list(), data = list()))
  }

  plots <- list()
  data_out <- list()

  t_grid <- seq(min(time_scores), max(time_scores), by = grid_step)
  age_grid <- t_grid + min(ages)

  for (j in seq_along(outcome_meta)) {
    meta <- outcome_meta[[j]]
    prefix <- meta$wide_prefix
    items <- item_names_list[[prefix]]
    if (is.null(items)) {
      next
    }

    gf_means <- extract_joint_gf_means(unstd, j = j)
    if (nrow(gf_means) == 0) {
      next
    }

    thresholds <- NULL
    if (meta$type %in% c("binary", "ordinal")) {
      thresholds <- extract_thresholds(unstd, items, ages, target = meta$plot_target)
    }

    points_data <- compute_trajectory_data(
      gf_means,
      ages,
      time_scores,
      thresholds,
      meta$type,
      target = meta$plot_target
    )

    thr_grid <- NULL
    if (!is.null(thresholds)) {
      tau_interp <- approx(x = thresholds$Age, y = thresholds$tau, xout = age_grid, rule = 2)$y
      thr_grid <- tibble(Age = age_grid, tau = tau_interp)
    }

    smooth_data <- compute_trajectory_data(
      gf_means,
      age_grid,
      t_grid,
      thr_grid,
      meta$type,
      target = meta$plot_target
    )

    if (!is.null(class_label_map) && length(class_label_map) > 0) {
      points_data <- points_data %>%
        mutate(LatentClass = factor(as.character(LatentClass), levels = names(class_label_map), labels = unname(class_label_map)))
      smooth_data <- smooth_data %>%
        mutate(LatentClass = factor(as.character(LatentClass), levels = names(class_label_map), labels = unname(class_label_map)))
    }

    p <- plot_lcga_trajectory(
      points_data = points_data,
      smooth_data = smooth_data,
      var_label = meta$label,
      y_label = trajectory_y_label(meta),
      is_percent = meta$type != "continuous"
    )

    plots[[prefix]] <- p
    data_out[[prefix]] <- list(points = points_data, smooth = smooth_data, meta = meta)
  }

  list(plots = plots, data = data_out)
}

collect_joint_model_trajectory_data <- function(
    all_results,
    outcome_index,
    item_names,
    ages,
    time_scores,
    var_type,
    plot_target = "highest",
    category_levels = NULL
) {
  model_rows <- purrr::map_dfr(names(all_results), function(model_id) {
    model_obj <- all_results[[model_id]]
    if (is.null(model_obj$parameters$unstandardized)) {
      return(tibble())
    }

    model_k <- suppressWarnings(as.integer(stringr::str_extract(model_id, "\\d+")))
    gf_means <- extract_joint_gf_means(model_obj$parameters$unstandardized, j = outcome_index)
    if (nrow(gf_means) == 0) {
      return(tibble())
    }

    thresholds <- if (tolower(var_type) %in% c("binary", "ordinal")) {
      extract_thresholds(model_obj$parameters$unstandardized, item_names, ages, target = plot_target)
    } else {
      NULL
    }

    points_data <- compute_trajectory_data(gf_means, ages, time_scores, thresholds, var_type, target = plot_target) %>%
      mutate(
        LatentClass = factor(LatentClass, levels = sort(unique(LatentClass)), labels = paste("Class", sort(unique(LatentClass)))),
        ModelK = factor(glue("K={model_k}"))
      )

    t_grid <- seq(min(time_scores), max(time_scores), by = 0.1)
    age_grid <- t_grid + min(ages)
    thr_grid <- if (!is.null(thresholds)) {
      tibble(Age = age_grid, tau = approx(x = thresholds$Age, y = thresholds$tau, xout = age_grid, rule = 2)$y)
    } else {
      NULL
    }

    smooth_data <- compute_trajectory_data(gf_means, age_grid, t_grid, thr_grid, var_type, target = plot_target) %>%
      mutate(
        LatentClass = factor(LatentClass, levels = sort(unique(LatentClass)), labels = paste("Class", sort(unique(LatentClass)))),
        ModelK = factor(glue("K={model_k}"))
      )

    bind_rows(
      points_data %>% mutate(DataType = "points"),
      smooth_data %>% mutate(DataType = "smooth")
    )
  })

  list(
    points = model_rows %>% filter(DataType == "points") %>% select(-DataType),
    smooth = model_rows %>% filter(DataType == "smooth") %>% select(-DataType)
  )
}
