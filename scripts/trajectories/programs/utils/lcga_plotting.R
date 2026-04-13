# ==============================================================================
# UTILITY: LCGA Plotting Core
# PURPOSE: Core parameter extraction, trajectory calculation, and light plotting
#          helpers for univariate LCGA work in the public trajectory bundle.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(glue)
library(purrr)
library(scales)
library(stringr)
library(tibble)
library(tidyr)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

extract_univariate_gf_means <- function(unstd) {
  gf <- unstd %>%
    filter(grepl("Means", paramHeader, ignore.case = TRUE), !is.na(LatentClass)) %>%
    mutate(paramU = toupper(trimws(param))) %>%
    filter(grepl("^[ISQ]", paramU)) %>%
    mutate(letter = substr(paramU, 1, 1)) %>%
    select(LatentClass, letter, est) %>%
    group_by(LatentClass, letter) %>%
    slice(1) %>%
    ungroup() %>%
    pivot_wider(names_from = letter, values_from = est) %>%
    mutate(across(any_of(c("I", "S", "Q")), as.numeric))

  if (nrow(gf) == 0) {
    return(gf)
  }

  if (!("Q" %in% names(gf))) {
    gf$Q <- 0
  }

  if (!all(c("I", "S") %in% names(gf))) {
    warning("Found some univariate parameters but missing I or S growth factors.")
    return(tibble())
  }

  gf
}

extract_thresholds <- function(param_df, item_names, ages, target = "lowest") {
  item_map <- tibble(item = item_names, Age = ages)
  pattern <- paste0("^(", paste(item_names, collapse = "|"), ")")

  thr <- param_df %>%
    filter(paramHeader == "Thresholds") %>%
    filter(grepl(pattern, param, ignore.case = TRUE)) %>%
    mutate(
      var = param,
      item = sub("(\\$|\\.|#).*", "", var),
      thr_num = suppressWarnings(as.integer(stringr::str_match(var, "\\$(\\d+)")[, 2])),
      thr_num = ifelse(is.na(thr_num), 1L, thr_num)
    ) %>%
    inner_join(item_map, by = "item")

  if (nrow(thr) == 0) {
    return(NULL)
  }

  thr_sel <- if (identical(target, "highest")) {
    thr %>% group_by(Age) %>% filter(thr_num == max(thr_num)) %>% ungroup()
  } else {
    thr %>% filter(thr_num == 1)
  }

  thr_final <- thr_sel %>%
    select(Age, tau = est) %>%
    mutate(tau = as.numeric(tau)) %>%
    distinct()

  tibble(Age = ages) %>%
    left_join(thr_final, by = "Age") %>%
    mutate(tau = tidyr::replace_na(tau, 0))
}

extract_all_thresholds <- function(param_df, item_names, ages) {
  item_map <- tibble(item = item_names, Age = ages)
  pattern <- paste0("^(", paste(item_names, collapse = "|"), ")")

  thr <- param_df %>%
    filter(paramHeader == "Thresholds") %>%
    filter(grepl(pattern, param, ignore.case = TRUE)) %>%
    mutate(
      var = param,
      item = sub("(\\$|\\.|#).*", "", var),
      thr_num = suppressWarnings(as.integer(stringr::str_match(var, "\\$(\\d+)")[, 2])),
      thr_num = ifelse(is.na(thr_num), 1L, thr_num)
    ) %>%
    inner_join(item_map, by = "item") %>%
    select(Age, k = thr_num, tau = est) %>%
    mutate(tau = as.numeric(tau)) %>%
    arrange(Age, k)

  if (nrow(thr) == 0) {
    return(NULL)
  }

  thr
}

compute_trajectory_data <- function(gf_means, ages, time_scores, thresholds = NULL, var_type, target = "lowest") {
  df <- tidyr::expand_grid(
    LatentClass = unique(gf_means$LatentClass),
    idx = seq_along(ages)
  ) %>%
    mutate(
      Age = ages[idx],
      t = time_scores[idx]
    ) %>%
    left_join(gf_means, by = "LatentClass")

  if (var_type %in% c("binary", "ordinal")) {
    if (!is.null(thresholds)) {
      df <- df %>% left_join(thresholds, by = "Age")
    } else {
      df <- df %>% mutate(tau = 0)
    }
  }

  df %>%
    mutate(across(any_of(c("I", "S", "Q", "tau")), as.numeric)) %>%
    mutate(
      eta = I + S * t + Q * (t^2),
      val = if (identical(var_type, "continuous")) {
        eta
      } else if (var_type %in% c("binary", "ordinal") && identical(target, "highest")) {
        1 - pnorm(tau - eta)
      } else if (var_type %in% c("binary", "ordinal")) {
        pnorm(tau - eta)
      } else {
        eta
      }
    )
}

compute_probs_all_cats <- function(gf_means, ages, time_scores, thresholds, cat_levels) {
  base <- tidyr::expand_grid(
    LatentClass = unique(gf_means$LatentClass),
    idx = seq_along(ages)
  ) %>%
    mutate(Age = ages[idx], t = time_scores[idx]) %>%
    left_join(gf_means, by = "LatentClass") %>%
    mutate(eta = I + S * t + Q * (t^2))

  calc_row <- function(eta, taus) {
    if (length(taus) == 0) {
      taus <- 0
    }
    taus <- sort(unique(taus))
    k <- length(taus) + 1
    p <- numeric(k)

    p[1] <- pnorm(taus[1] - eta)
    if (k > 2) {
      for (i in 2:(k - 1)) {
        p[i] <- pnorm(taus[i] - eta) - pnorm(taus[i - 1] - eta)
      }
    }
    p[k] <- 1 - pnorm(taus[length(taus)] - eta)
    p
  }

  results <- list()

  for (ag in unique(base$Age)) {
    taus_ag <- thresholds %>% filter(near(Age, ag)) %>% pull(tau)
    sub_base <- base %>% filter(near(Age, ag))
    probs_mat <- t(sapply(sub_base$eta, calc_row, taus = taus_ag))

    n_cols <- ncol(probs_mat)
    curr_levels <- if (!is.null(cat_levels) && length(cat_levels) >= n_cols) {
      cat_levels[1:n_cols]
    } else {
      paste0("Cat_", 0:(n_cols - 1))
    }
    colnames(probs_mat) <- curr_levels

    sub_res <- bind_cols(sub_base %>% select(LatentClass, Age, t), as_tibble(probs_mat)) %>%
      pivot_longer(cols = all_of(curr_levels), names_to = "Category", values_to = "Prob")

    results[[paste0(ag)]] <- sub_res
  }

  bind_rows(results) %>%
    mutate(Category = factor(Category, levels = cat_levels))
}

infer_positive_category_label <- function(var_label, category_levels = NULL) {
  if (!is.null(category_levels) && length(category_levels) > 0) {
    last_level <- as.character(tail(category_levels, 1))
    if (!tolower(last_level) %in% c("yes", "1", "true")) {
      return(last_level)
    }
  }

  var_label_lower <- tolower(var_label %||% "")

  case_when(
    str_detect(var_label_lower, "overcrowd") ~ "Overcrowding",
    str_detect(var_label_lower, "without|absence") ~ "Absence",
    TRUE ~ "Yes"
  )
}

prepare_distribution_plot_data <- function(
  smooth_data,
  points_data = NULL,
  var_type = NULL,
  var_label = NULL,
  category_levels = NULL
) {
  facet_categories <- TRUE
  subtitle <- NULL

  if (!is.null(var_type) && tolower(var_type) == "binary") {
    keep_category <- if (!is.null(category_levels) && length(category_levels) > 0) {
      as.character(tail(category_levels, 1))
    } else {
      as.character(tail(levels(smooth_data$Category), 1))
    }

    display_label <- infer_positive_category_label(var_label, category_levels)

    smooth_data <- smooth_data %>%
      filter(as.character(Category) == keep_category) %>%
      mutate(Category = factor(display_label, levels = display_label))

    if (!is.null(points_data)) {
      points_data <- points_data %>%
        filter(as.character(Category) == keep_category) %>%
        mutate(Category = factor(display_label, levels = display_label))
    }

    facet_categories <- FALSE
    subtitle <- glue("Probability of {display_label}")
  } else if (!is.null(smooth_data$Category)) {
    category_levels <- levels(factor(smooth_data$Category))
    display_levels <- glue("Probability of {category_levels}")

    smooth_data <- smooth_data %>%
      mutate(Category = factor(glue("Probability of {as.character(Category)}"), levels = display_levels))

    if (!is.null(points_data)) {
      points_data <- points_data %>%
        mutate(Category = factor(glue("Probability of {as.character(Category)}"), levels = display_levels))
    }
  }

  list(
    smooth_data = smooth_data,
    points_data = points_data,
    facet_categories = facet_categories,
    subtitle = subtitle
  )
}

plot_lcga_trajectory <- function(points_data, smooth_data, var_label, y_label, is_percent = TRUE) {
  x_breaks <- sort(unique(points_data$Age))
  if (length(x_breaks) == 0) {
    x_breaks <- waiver()
  }

  p <- ggplot() +
    geom_line(
      data = smooth_data,
      aes(x = Age, y = val, color = LatentClass, group = LatentClass),
      linewidth = 1.2
    ) +
    geom_point(
      data = points_data,
      aes(x = Age, y = val, color = LatentClass),
      size = 3
    ) +
    scale_x_continuous(breaks = x_breaks) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
    labs(
      title = glue("Latent Trajectories: {var_label}"),
      y = y_label,
      x = "Child Age (Years)",
      color = "Latent Class"
    )

  if (is_percent) {
    p <- p + scale_y_continuous(labels = scales::percent, limits = c(0, 1))
  }

  p
}

plot_category_distribution <- function(
  smooth_data,
  var_label,
  points_data = NULL,
  facet_categories = TRUE,
  subtitle = NULL
) {
  x_breaks <- if (!is.null(points_data) && nrow(points_data) > 0) {
    sort(unique(points_data$Age))
  } else {
    seq(ceiling(min(smooth_data$Age, na.rm = TRUE)), floor(max(smooth_data$Age, na.rm = TRUE)), 1)
  }

  p <- ggplot(smooth_data, aes(x = Age, y = Prob, color = LatentClass, group = LatentClass)) +
    geom_line(linewidth = 1.2)

  if (!is.null(points_data) && nrow(points_data) > 0) {
    p <- p +
      geom_point(
        data = points_data,
        aes(x = Age, y = Prob, color = LatentClass),
        size = 3
      )
  }

  if (isTRUE(facet_categories)) {
    p <- p + facet_wrap(~Category, scales = "fixed")
  }

  p +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_x_continuous(breaks = x_breaks) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      strip.text = element_text(size = 13, hjust = 0.5),
      panel.border = element_rect(color = "grey80", fill = NA)
    ) +
    labs(
      title = glue("{var_label}"),
      subtitle = subtitle,
      y = "Probability",
      x = "Child Age (Years)",
      color = "Class"
    )
}

compute_target_indicator <- function(value, var_type, plot_target = "highest", category_levels = NULL) {
  if (is.na(value)) {
    return(NA_real_)
  }

  if (tolower(var_type %||% "") == "continuous") {
    return(as.numeric(value))
  }

  numeric_value <- suppressWarnings(as.numeric(value))
  if (is.na(numeric_value)) {
    return(NA_real_)
  }

  max_code <- if (!is.null(category_levels) && length(category_levels) > 0) {
    length(category_levels) - 1
  } else {
    numeric_value
  }

  if (tolower(plot_target %||% "highest") == "lowest") {
    as.numeric(numeric_value == 0)
  } else {
    as.numeric(numeric_value == max_code)
  }
}

prepare_observed_target_data <- function(
    df_mplus,
    item_names,
    ages,
    class_assign_data = NULL,
    id_var = "USUBJID",
    var_type,
    plot_target = "highest",
    category_levels = NULL,
    class_levels = NULL,
    missing_code = -9999
) {
  if (is.null(df_mplus) || length(item_names) == 0) {
    return(tibble())
  }

  observed_df <- df_mplus %>%
    select(any_of(c(id_var, item_names))) %>%
    pivot_longer(cols = all_of(item_names), names_to = "item", values_to = "raw_value") %>%
    mutate(Age = as.integer(stringr::str_extract(item, "\\d+$"))) %>%
    filter(!is.na(raw_value), raw_value != missing_code) %>%
    mutate(
      val = purrr::map_dbl(
        raw_value,
        compute_target_indicator,
        var_type = var_type,
        plot_target = plot_target,
        category_levels = category_levels
      )
    ) %>%
    filter(!is.na(val)) %>%
    select(any_of(id_var), Age, val)

  if (!is.null(class_assign_data) && nrow(class_assign_data) > 0) {
    observed_df <- observed_df %>%
      left_join(
        class_assign_data %>%
          select(any_of(id_var), Assigned_Class_Label) %>%
          rename(LatentClass = Assigned_Class_Label),
        by = id_var
      ) %>%
      filter(!is.na(LatentClass))
  } else {
    observed_df <- observed_df %>% mutate(LatentClass = factor("Observed"))
  }

  observed_df %>%
    mutate(
      LatentClass = if (!is.null(class_levels) && length(class_levels) > 0) {
        factor(LatentClass, levels = class_levels)
      } else {
        factor(LatentClass)
      }
    )
}

plot_estimated_plus_observed_by_class <- function(
    observed_data,
    points_data,
    smooth_data,
    var_label,
    y_label,
    is_percent = TRUE,
    spaghetti_n = 100
) {
  if (nrow(observed_data) == 0 || nrow(points_data) == 0 || nrow(smooth_data) == 0) {
    return(NULL)
  }

  id_col <- names(observed_data)[1]
  sampled_ids <- observed_data %>%
    distinct(LatentClass, SubjectID = .data[[id_col]]) %>%
    group_by(LatentClass) %>%
    slice_head(n = spaghetti_n) %>%
    ungroup()

  join_by <- c("LatentClass", stats::setNames("SubjectID", id_col))
  observed_sample <- observed_data %>%
    semi_join(sampled_ids, by = join_by)

  x_breaks <- sort(unique(points_data$Age))
  if (length(x_breaks) == 0) {
    x_breaks <- waiver()
  }

  p <- ggplot() +
    geom_line(
      data = observed_sample,
      aes(x = Age, y = val, group = .data[[id_col]]),
      color = "grey75",
      alpha = 0.20,
      linewidth = 0.4
    ) +
    geom_point(
      data = observed_sample,
      aes(x = Age, y = val),
      color = "grey70",
      alpha = 0.25,
      size = 0.7
    ) +
    geom_line(
      data = smooth_data,
      aes(x = Age, y = val, color = LatentClass, group = LatentClass),
      linewidth = 1.4
    ) +
    geom_point(
      data = points_data,
      aes(x = Age, y = val, color = LatentClass),
      size = 2.4
    ) +
    facet_wrap(~LatentClass) +
    scale_x_continuous(breaks = x_breaks) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = glue("Estimated Means With Observed Trajectories: {var_label}"),
      subtitle = "Observed lines are sampled within assigned class; colored line is the fitted class-specific mean",
      y = y_label,
      x = "Child Age (Years)",
      color = "Latent Class"
    )

  if (is_percent) {
    p <- p + scale_y_continuous(labels = scales::percent, limits = c(0, 1))
  }

  p
}

plot_model_trajectory_grid <- function(
    points_all,
    smooth_all,
    var_label,
    y_label,
    is_percent = TRUE
) {
  if (nrow(points_all) == 0 || nrow(smooth_all) == 0) {
    return(NULL)
  }

  x_breaks <- sort(unique(points_all$Age))
  if (length(x_breaks) == 0) {
    x_breaks <- waiver()
  }

  p <- ggplot() +
    geom_line(
      data = smooth_all,
      aes(x = Age, y = val, color = LatentClass, group = interaction(ModelK, LatentClass)),
      linewidth = 1.1
    ) +
    geom_point(
      data = points_all,
      aes(x = Age, y = val, color = LatentClass),
      size = 2
    ) +
    facet_wrap(~ModelK) +
    scale_x_continuous(breaks = x_breaks) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = glue("Estimated Mean Trajectories Across Fitted Models: {var_label}"),
      subtitle = "Each panel shows one fitted class solution",
      y = y_label,
      x = "Child Age (Years)",
      color = "Latent Class"
    )

  if (is_percent) {
    p <- p + scale_y_continuous(labels = scales::percent, limits = c(0, 1))
  }

  p
}

collect_univariate_model_trajectory_data <- function(
    all_results,
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
    gf_means <- extract_univariate_gf_means(model_obj$parameters$unstandardized)
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

plot_class_sizes <- function(class_assign_data, class_color_map = NULL) {
  if (is.null(class_assign_data) || nrow(class_assign_data) == 0) {
    return(NULL)
  }

  size_df <- class_assign_data %>%
    count(Assigned_Class_Label, name = "n") %>%
    mutate(prop = n / sum(n))

  p <- ggplot(size_df, aes(x = Assigned_Class_Label, y = prop, fill = Assigned_Class_Label)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), vjust = -0.3, size = 4) +
    scale_y_continuous(labels = scales::percent, limits = c(0, min(1, max(size_df$prop) + 0.1))) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 20, hjust = 1),
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = "Class Size Distribution",
      y = "Share of Assigned Sample"
    )

  if (!is.null(class_color_map) && length(class_color_map) > 0) {
    p <- p + scale_fill_manual(values = class_color_map, drop = FALSE)
  }

  p
}
