# ==============================================================================
# UTILITY: Class Overlap Reporting
# PURPOSE: Summarize posterior overlap and estimated-profile separation for
#          selected LCGA solutions in the public trajectory bundle.
# ==============================================================================

library(dplyr)
library(glue)
library(gt)
library(scales)
library(tibble)
library(tidyr)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

infer_selected_model_id <- function(res, k) {
  model_id <- paste0("lcga_model_class_", as.integer(k), ".out")
  if (!model_id %in% names(res$all_results)) {
    stop(glue("Selected model id {model_id} was not found in all_results."))
  }
  model_id
}

build_joint_class_feature_matrix <- function(res, model_obj, ages, time_scores) {
  age_grid <- seq(min(ages), max(ages), by = 0.25)
  t_grid <- age_grid - min(ages)

  feature_rows <- purrr::map_dfr(seq_along(res$meta_long), function(j) {
    meta <- res$meta_long[[j]]
    gf_means <- extract_joint_gf_means(model_obj$parameters$unstandardized, j = j)
    if (nrow(gf_means) == 0) {
      return(tibble())
    }

    thresholds <- if (tolower(meta$type) %in% c("binary", "ordinal")) {
      extract_thresholds(
        model_obj$parameters$unstandardized,
        res$item_names_list[[meta$wide_prefix]],
        ages,
        target = meta$plot_target
      )
    } else {
      NULL
    }

    thr_grid <- if (!is.null(thresholds)) {
      tibble(
        Age = age_grid,
        tau = approx(x = thresholds$Age, y = thresholds$tau, xout = age_grid, rule = 2)$y
      )
    } else {
      NULL
    }

    compute_trajectory_data(
      gf_means = gf_means,
      ages = age_grid,
      time_scores = t_grid,
      thresholds = thr_grid,
      var_type = meta$type,
      target = meta$plot_target
    ) %>%
      transmute(
        LatentClass = as.integer(as.character(LatentClass)),
        Feature = glue("{meta$wide_prefix}_A{formatC(Age, format = 'f', digits = 2)}"),
        Value = as.numeric(val)
      )
  })

  feature_rows %>%
    pivot_wider(names_from = Feature, values_from = Value) %>%
    arrange(LatentClass)
}

build_univariate_class_feature_matrix <- function(res, model_obj, ages, time_scores) {
  age_grid <- seq(min(ages), max(ages), by = 0.25)
  t_grid <- age_grid - min(ages)
  meta <- res$meta_long[[1]]
  gf_means <- extract_univariate_gf_means(model_obj$parameters$unstandardized)

  thresholds <- if (tolower(meta$type) %in% c("binary", "ordinal")) {
    extract_thresholds(
      model_obj$parameters$unstandardized,
      names(res$df_mplus)[names(res$df_mplus) != "USUBJID"],
      ages,
      target = meta$plot_target
    )
  } else {
    NULL
  }

  thr_grid <- if (!is.null(thresholds)) {
    tibble(
      Age = age_grid,
      tau = approx(x = thresholds$Age, y = thresholds$tau, xout = age_grid, rule = 2)$y
    )
  } else {
    NULL
  }

  compute_trajectory_data(
    gf_means = gf_means,
    ages = age_grid,
    time_scores = t_grid,
    thresholds = thr_grid,
    var_type = meta$type,
    target = meta$plot_target
  ) %>%
    transmute(
      LatentClass = as.integer(as.character(LatentClass)),
      Feature = glue("{meta$wide_prefix}_A{formatC(Age, format = 'f', digits = 2)}"),
      Value = as.numeric(val)
    ) %>%
    pivot_wider(names_from = Feature, values_from = Value) %>%
    arrange(LatentClass)
}

compute_lcga_estimate_distance <- function(res, k, ages, time_scores) {
  model_id <- infer_selected_model_id(res, k)
  model_obj <- res$all_results[[model_id]]

  feature_mat <- if (length(res$meta_long) > 1) {
    build_joint_class_feature_matrix(res, model_obj, ages, time_scores)
  } else {
    build_univariate_class_feature_matrix(res, model_obj, ages, time_scores)
  }

  if (nrow(feature_mat) < 2) {
    return(tibble(
      Class_i = integer(),
      Class_j = integer(),
      EstimateDistance = numeric()
    ))
  }

  mat <- as.matrix(feature_mat %>% select(-LatentClass))
  rownames(mat) <- feature_mat$LatentClass
  dist_long <- as.matrix(dist(mat, method = "euclidean")) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Class_i") %>%
    pivot_longer(-Class_i, names_to = "Class_j", values_to = "EstimateDistance") %>%
    mutate(
      Class_i = as.integer(Class_i),
      Class_j = as.integer(Class_j)
    )

  as_tibble(dist_long)
}

compute_class_overlap_metrics <- function(savedata_post, class_labels = NULL) {
  if ("Assigned_Class" %in% names(savedata_post) && !"AssignedClass" %in% names(savedata_post)) {
    savedata_post <- savedata_post %>% rename(AssignedClass = Assigned_Class)
  }
  prob_cols <- grep("^CPROB\\d+$", names(savedata_post), value = TRUE)
  n_classes <- length(prob_cols)
  stopifnot(n_classes > 1)

  if (is.null(class_labels) && "Assigned_Class_Label" %in% names(savedata_post)) {
    inferred <- savedata_post %>%
      distinct(AssignedClass, Assigned_Class_Label) %>%
      arrange(AssignedClass)
    class_labels <- setNames(as.character(inferred$Assigned_Class_Label), inferred$AssignedClass)
  }

  label_map <- class_labels %||% setNames(paste("Class", seq_len(n_classes)), seq_len(n_classes))

  overlap_df <- expand_grid(
    AssignedClass = seq_len(n_classes),
    CompetingClass = seq_len(n_classes)
  ) %>%
    rowwise() %>%
    mutate(
      MeanPosterior = mean(
        savedata_post %>%
          filter(AssignedClass == .env$AssignedClass) %>%
          pull(paste0("CPROB", .env$CompetingClass)),
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    mutate(
      AssignedLabel = unname(label_map[as.character(AssignedClass)]),
      CompetingLabel = unname(label_map[as.character(CompetingClass)])
    )

  prob_mat <- as.matrix(savedata_post[, prob_cols, drop = FALSE])
  top1 <- max.col(prob_mat, ties.method = "first")
  top1_prob <- prob_mat[cbind(seq_len(nrow(prob_mat)), top1)]
  prob_mat_masked <- prob_mat
  prob_mat_masked[cbind(seq_len(nrow(prob_mat_masked)), top1)] <- -Inf
  top2 <- max.col(prob_mat_masked, ties.method = "first")
  top2_prob <- prob_mat[cbind(seq_len(nrow(prob_mat)), top2)]

  margin_df <- savedata_post %>%
    transmute(
      AssignedClass = AssignedClass,
      AssignedLabel = unname(label_map[as.character(AssignedClass)]),
      TopClass = top1,
      TopLabel = unname(label_map[as.character(top1)]),
      TopProb = top1_prob,
      SecondClass = top2,
      SecondLabel = unname(label_map[as.character(top2)]),
      SecondProb = top2_prob,
      Margin = top1_prob - top2_prob
    )

  summary_df <- margin_df %>%
    group_by(AssignedClass, AssignedLabel) %>%
    summarise(
      NAssigned = n(),
      MeanTopProb = mean(TopProb, na.rm = TRUE),
      MeanSecondProb = mean(SecondProb, na.rm = TRUE),
      MeanMargin = mean(Margin, na.rm = TRUE),
      P10Margin = quantile(Margin, probs = 0.10, na.rm = TRUE),
      TopCompetitorClass = as.integer(names(sort(table(SecondClass), decreasing = TRUE)[1])),
      .groups = "drop"
    ) %>%
    mutate(
      TopCompetitorLabel = unname(label_map[as.character(TopCompetitorClass)])
    )

  list(
    overlap = overlap_df,
    margins = margin_df,
    summary = summary_df
  )
}

plot_posterior_overlap_heatmap <- function(overlap_df, fill_limits = c(0, 1)) {
  ggplot(overlap_df, aes(x = CompetingLabel, y = AssignedLabel, fill = MeanPosterior)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = percent(MeanPosterior, accuracy = 0.1)), size = 3) +
    scale_fill_gradient(low = "#F3F6FA", high = "#0B7285", limits = fill_limits, labels = percent) +
    labs(
      title = "Posterior Overlap Heatmap",
      subtitle = "Diagonal cells are AvePP; off-diagonals show competing-class posterior mass.",
      x = "Competing Class",
      y = "Assigned Class",
      fill = "Mean posterior"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank()
    )
}

plot_posterior_margin_distribution <- function(margin_df, color_map = NULL) {
  p <- ggplot(margin_df, aes(x = AssignedLabel, y = Margin, fill = AssignedLabel)) +
    geom_violin(alpha = 0.7, trim = FALSE, color = NA) +
    geom_boxplot(width = 0.18, outlier.alpha = 0.15, fill = "white", color = "#333333") +
    coord_flip() +
    scale_y_continuous(labels = number_format(accuracy = 0.01), limits = c(0, 1)) +
    labs(
      title = "Posterior Margin Distribution",
      subtitle = "Margin = highest posterior minus second-highest posterior for each assigned case.",
      x = NULL,
      y = "Posterior margin"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )

  if (!is.null(color_map)) {
    p <- p + scale_fill_manual(values = color_map, drop = FALSE)
  }

  p
}

build_class_overlap_summary_table <- function(
    summary_df,
    distance_df = NULL,
    analysis_label,
    selected_k
) {
  summary_out <- summary_df

  if (!is.null(distance_df) && nrow(distance_df) > 0) {
    summary_out <- summary_out %>%
      left_join(
        distance_df %>%
          rename(
            AssignedClass = Class_i,
            TopCompetitorClass = Class_j
          ),
        by = c("AssignedClass", "TopCompetitorClass")
      )
  }

  gt(summary_out) %>%
    tab_header(
      title = md(glue("**Class Overlap Summary - {analysis_label}**")),
      subtitle = glue("Selected solution: K = {selected_k}")
    ) %>%
    fmt_percent(columns = c(MeanTopProb, MeanSecondProb, MeanMargin, P10Margin), decimals = 1) %>%
    fmt_number(columns = any_of("EstimateDistance"), decimals = 3) %>%
    cols_label(
      AssignedLabel = "Assigned Class",
      NAssigned = "n assigned",
      MeanTopProb = "Mean top posterior",
      MeanSecondProb = "Mean second posterior",
      MeanMargin = "Mean margin",
      P10Margin = "10th pct margin",
      TopCompetitorLabel = "Main competing class",
      EstimateDistance = "Estimate distance"
    ) %>%
    tab_source_note("Estimate distance is the Euclidean distance between class-specific estimated trajectory profiles.")
}
