# ==============================================================================
# SCRIPT: Report Harmonized Joint LCGA Results
# PURPOSE: Generate first-pass tables and figures for the public joint trajectory bundle.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(gt)
  library(here)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

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

trajectory_root <- find_trajectory_root()
trajectory_path <- function(...) file.path(trajectory_root, ...)

source(trajectory_path("programs", "utils", "analysis_configuration.R"), local = TRUE)
source(trajectory_path("programs", "utils", "class_registry.R"), local = TRUE)
source(trajectory_path("programs", "utils", "data_prep.R"), local = TRUE)
source(trajectory_path("programs", "utils", "lcga_plotting.R"), local = TRUE)
source(trajectory_path("programs", "utils", "lcga_plotting_joint.R"), local = TRUE)
source(trajectory_path("programs", "utils", "lcga_reporting_core.R"), local = TRUE)
source(trajectory_path("programs", "utils", "class_overlap_reporting.R"), local = TRUE)
source(trajectory_path("programs", "utils", "reporting_standards.R"), local = TRUE)
source(trajectory_path("programs", "utils", "lcga_selection.R"), local = TRUE)

run_trajectory_joint_report <- function(spec_names = NULL, time_labels = NULL) {
  dt_multi <- tryCatch(
    load_long_data(cfg$multistudies_path, cfg$multistudies_object),
    error = function(e) NULL
  )
  all_joint_reports <- list()
  joint_specs <- filter_lcga_specs(lcga_specs, spec_names = spec_names)
  joint_specs <- Filter(function(spec) length(spec$outcome_prefixes) > 1, joint_specs)
  report_time_labels <- resolve_report_time_labels(time_labels, config = cfg)

  for (time_label in report_time_labels) {
    ages <- ages_from_time_label(time_label)
    time_scores <- 0:(length(ages) - 1)
    period_label <- period_label_from_time_label(time_label)

    for (lcga_spec in joint_specs) {
      paths <- derive_lcga_paths(cfg$base_mplus_dir, lcga_spec$name, time_label, is_joint = TRUE)
      rds_file <- paths$results_rds_path
      if (!file.exists(rds_file)) {
        next
      }

      message(glue("\n>>> Reporting trajectory bundle (joint): {lcga_spec$label} ({time_label})"))
      res <- readRDS(rds_file)
      fit_df <- normalise_fit_summary(res$fit_summary_aug)
      if (nrow(fit_df) == 0) {
        next
      }

      best_selection <- tryCatch(
        select_best_solution(fit_df, spec_name = lcga_spec$name, time_label = time_label, config = cfg),
        error = function(e) {
          message("    ! Selection error: ", e$message)
          NULL
        }
      )

      if (is.null(best_selection) || is.na(best_selection$model_id)) {
        next
      }

      best_k <- best_selection$Classes
      best_model <- res$all_results[[best_selection$model_id]]
      if (is.null(best_model)) {
        next
      }

      class_meta <- resolve_class_metadata(
        variable = lcga_spec$name,
        classes = best_k,
        period = period_label,
        model_family = "LCGA"
      )
      label_map_vec <- setNames(class_meta$Display_Label, class_meta$LatentClass)
      color_map_vec <- setNames(class_meta$Color_Hex, class_meta$Display_Label)
      manual_override_key <- NULL
      if (isTRUE(cfg$use_manual_overrides)) {
        time_key <- paste0(lcga_spec$name, "_", time_label)
        manual_override_key <- if (!is.null(cfg$manual_class_overrides[[time_key]])) {
          time_key
        } else if (!is.null(cfg$manual_class_overrides[[lcga_spec$name]])) {
          lcga_spec$name
        } else {
          NULL
        }
      }

      class_assign_data <- extract_lcga_class_assignments(best_model, class_labels = label_map_vec, id_var = cfg$id_var)
      avail_tbl <- NULL
      desc_tbl <- NULL
      avepp_tbl <- NULL
      est_tbl <- NULL
      p_comp <- NULL
      model_grid_plots <- list()
      observed_plots <- list()
      final_solution_tbl <- NULL
      missingness_tbl <- NULL
      structural_sensitivity_tbl <- NULL
      sensitivity_tbl <- NULL
      overlap_tbl <- NULL
      checklist_tbl <- NULL
      smart_tbl <- NULL
      overlap_heatmap <- NULL
      overlap_margin_plot <- NULL
      traj_out <- list(plots = list(), data = list())
      sensitivity_summary <- tibble()
      syntax_dir <- res$meta$output_dir %||% paths$output_dir

      if (!is.null(res$df_mplus)) {
        avail_tbl <- map_dfr(ages, function(ag) {
          cols <- unlist(lapply(res$item_names_list, function(it) grep(paste0(ag, "$"), it, value = TRUE)))
          n_obs <- res$df_mplus %>%
            select(all_of(cols)) %>%
            filter(if_any(everything(), ~ . != cfg$missing_code)) %>%
            nrow()
          tibble(age = ag, n_total = nrow(res$df_mplus), n_obs = n_obs, pct_obs = n_obs / n_total)
        }) %>%
          gt() %>%
          tab_header(title = md(glue("**Data Availability - {lcga_spec$label}**"))) %>%
          fmt_percent(columns = pct_obs, decimals = 1)

        save_gt_html(avail_tbl, file.path(paths$tables_dir, glue("00_Data_Availability_{lcga_spec$name}.html")))

        desc_tbl <- map_dfr(names(res$item_names_list), function(prefix) {
          meta <- Filter(function(x) identical(x$wide_prefix, prefix), res$meta_long)[[1]]
          cat_levels <- meta$transform_levels %||% paste0("Cat_", 0:4)
          items <- res$item_names_list[[prefix]]
          res$df_mplus %>%
            select(all_of(items)) %>%
            pivot_longer(everything(), names_to = "item", values_to = "value") %>%
            filter(!is.na(value), value != cfg$missing_code) %>%
            mutate(
              variable = meta$label,
              age = as.integer(str_extract(item, "\\d+$")),
              category = factor(value, levels = 0:(length(cat_levels) - 1), labels = cat_levels[seq_len(length(cat_levels))])
            ) %>%
            count(variable, age, category, name = "n") %>%
            group_by(variable, age) %>%
            mutate(pct = n / sum(n), stat = glue("{n} ({sprintf('%.1f', pct * 100)}%)")) %>%
            ungroup() %>%
            select(variable, age, category, stat)
        }) %>%
          pivot_wider(names_from = category, values_from = stat, values_fill = "0 (0.0%)") %>%
          gt(groupname_col = "variable") %>%
          tab_header(title = md(glue("**Descriptive Statistics - {lcga_spec$label}**")))

        save_gt_html(desc_tbl, file.path(paths$tables_dir, glue("00_Descriptive_Statistics_{lcga_spec$name}.html")))
      }

      fit_tbl <- fit_df %>%
        select(Classes, BIC, Entropy, BLRT_p, VLMR_p, MinClassPct) %>%
        gt() %>%
        tab_header(title = md(glue("**Fit Summary - {lcga_spec$label}**"))) %>%
        fmt_number(columns = c(BIC), decimals = 1) %>%
        fmt_number(columns = c(Entropy, BLRT_p, VLMR_p), decimals = 3) %>%
        fmt_percent(columns = MinClassPct, decimals = 1) %>%
        tab_style(style = cell_fill(color = "lightgreen"), locations = cells_body(rows = Classes == best_k))

      save_gt_html(fit_tbl, file.path(paths$tables_dir, glue("01_fit_summary_{lcga_spec$name}.html")))

      if (!is.null(class_assign_data)) {
        avepp_summary <- summarise_lcga_avepp(class_assign_data)
        avepp_tbl <- avepp_summary$by_class %>%
          gt() %>%
          tab_header(title = md(glue("**Average Posterior Probabilities - {lcga_spec$label}**"))) %>%
          fmt_number(columns = DiagonalAvePP, decimals = 3)

        save_gt_html(avepp_tbl, file.path(paths$tables_dir, glue("02_AvePP_{lcga_spec$name}.html")))

        p_sizes <- plot_class_sizes(class_assign_data, class_color_map = color_map_vec)
        if (!is.null(p_sizes)) {
          ggsave(file.path(paths$figures_dir, glue("04_class_sizes_{lcga_spec$name}.png")), p_sizes, width = 10, height = 6)
        }
      }

      if (!is.null(best_model$parameters$unstandardized)) {
        est_rows <- purrr::map_dfr(seq_along(res$meta_long), function(j) {
          meta <- res$meta_long[[j]]
          gf <- extract_joint_gf_means(best_model$parameters$unstandardized, j = j)
          if (nrow(gf) == 0) {
            return(tibble())
          }

          gf %>%
            mutate(
              Outcome = meta$label,
              LatentClass = factor(as.character(LatentClass), levels = names(label_map_vec), labels = unname(label_map_vec))
            )
        })

        if (nrow(est_rows) > 0) {
          est_tbl <- est_rows %>%
            select(Outcome, LatentClass, I, S, Q) %>%
            gt(groupname_col = "Outcome") %>%
            tab_header(title = md(glue("**Growth Factor Means - {lcga_spec$label}**"))) %>%
            fmt_number(columns = any_of(c("I", "S", "Q")), decimals = 3)

          save_gt_html(est_tbl, file.path(paths$tables_dir, glue("04_Estimates_{lcga_spec$name}.html")))
        }

        traj_out <- plot_joint_trajectories(
          best_model = best_model,
          outcome_meta = res$meta_long,
          item_names_list = res$item_names_list,
          ages = ages,
          time_scores = time_scores,
          class_label_map = label_map_vec
        )

        if (length(traj_out$plots) > 0) {
          for (nm in names(traj_out$plots)) {
            p_curr <- apply_class_color_scale(traj_out$plots[[nm]], color_map_vec, aes_name = "color")
            ggsave(file.path(paths$figures_dir, glue("02_trajectory_{nm}.png")), p_curr, width = 12, height = 6)
          }

          p_comp <- plot_joint_composite(traj_out, method = cfg$composite_method)
          p_comp <- apply_class_color_scale(p_comp, color_map_vec, aes_name = "color")
          ggsave(file.path(paths$figures_dir, glue("03_composite_{lcga_spec$name}.png")), p_comp, width = 12, height = 6)

          for (nm in names(traj_out$plots)) {
            j_idx <- which(vapply(res$meta_long, function(x) identical(x$wide_prefix, nm), logical(1)))
            if (length(j_idx) != 1) {
              next
            }

            meta <- res$meta_long[[j_idx]]
            y_label <- trajectory_y_label(meta)

            model_traj_data <- collect_joint_model_trajectory_data(
              all_results = res$all_results,
              outcome_index = j_idx,
              item_names = res$item_names_list[[nm]],
              ages = ages,
              time_scores = time_scores,
              var_type = meta$type,
              plot_target = meta$plot_target,
              category_levels = meta$transform_levels %||% meta$levels
            )
            p_model_grid <- plot_model_trajectory_grid(
              points_all = model_traj_data$points,
              smooth_all = model_traj_data$smooth,
              var_label = meta$label,
              y_label = y_label,
              is_percent = meta$type != "continuous"
            )
            if (!is.null(p_model_grid)) {
              model_grid_plots[[nm]] <- p_model_grid
              n_models <- dplyr::n_distinct(model_traj_data$points$ModelK)
              ggsave(
                file.path(paths$figures_dir, glue("03b_all_model_trajectories_{nm}.png")),
                p_model_grid,
                width = 16,
                height = max(6, 3 * ceiling(n_models / 3))
              )
            }

            if (!is.null(class_assign_data) && !is.null(res$df_mplus) && !is.null(traj_out$data[[nm]])) {
              observed_data <- prepare_observed_target_data(
                df_mplus = res$df_mplus,
                item_names = res$item_names_list[[nm]],
                ages = ages,
                class_assign_data = class_assign_data,
                id_var = cfg$id_var,
                var_type = meta$type,
                plot_target = meta$plot_target,
                category_levels = meta$transform_levels %||% meta$levels,
                class_levels = unname(label_map_vec),
                missing_code = cfg$missing_code
              )
              p_observed <- plot_estimated_plus_observed_by_class(
                observed_data = observed_data,
                points_data = traj_out$data[[nm]]$points,
                smooth_data = traj_out$data[[nm]]$smooth,
                var_label = meta$label,
                y_label = y_label,
                is_percent = meta$type != "continuous",
                spaghetti_n = 100
              )
              if (!is.null(p_observed)) {
                p_observed <- apply_class_color_scale(p_observed, color_map_vec, aes_name = "color")
                observed_plots[[nm]] <- p_observed
                ggsave(
                  file.path(paths$figures_dir, glue("03c_estimated_plus_observed_{nm}.png")),
                  p_observed,
                  width = 16,
                  height = max(7, 3 * ceiling(best_k / 3))
                )
              }
            }
          }
        }
      }

      tryCatch({
        final_solution_df <- build_lcga_final_solution_summary(
          best_model = best_model,
          class_assign_data = class_assign_data,
          outcome_meta = res$meta_long,
          analysis_label = lcga_spec$label,
          is_joint = TRUE
        )
        final_solution_tbl <- generate_lcga_final_solution_summary_table(
          final_solution_df,
          analysis_label = lcga_spec$label
        )
        save_gt_html(final_solution_tbl, file.path(paths$tables_dir, glue("06a_Final_Solution_Summary_{lcga_spec$name}.html")))
      }, error = function(e) {
        message("    ! Final Solution Summary Error: ", e$message)
      })

      tryCatch({
        candidate_vars <- unique(c(
          cfg$baseline_confounders %||% character(),
          cfg$mediator_outcome_confounders %||% character(),
          "SEX",
          "DEPBIN"
        ))
        missingness_df <- build_missingness_correlates_summary(
          df_mplus = res$df_mplus,
          item_names_list = res$item_names_list,
          outcome_meta = res$meta_long,
          id_var = cfg$id_var,
          missing_code = cfg$missing_code,
          dt_multi = dt_multi,
          candidate_vars = candidate_vars
        )
        missingness_tbl <- generate_missingness_correlates_table(
          missingness_df,
          analysis_label = lcga_spec$label
        )
        save_gt_html(missingness_tbl, file.path(paths$tables_dir, glue("06b_Missingness_Correlates_{lcga_spec$name}.html")))
      }, error = function(e) {
        message("    ! Missingness Correlates Error: ", e$message)
      })

      tryCatch({
        structural_sensitivity_tbl <- generate_lcga_structural_sensitivity_table(
          structural_sensitivity_audit = res$structural_sensitivity %||% NULL,
          analysis_label = lcga_spec$label
        )
        if (!is.null(structural_sensitivity_tbl)) {
          save_gt_html(structural_sensitivity_tbl, file.path(paths$tables_dir, glue("07a_Structural_Sensitivity_{lcga_spec$name}.html")))
        }
      }, error = function(e) {
        message("    ! Structural Sensitivity Error: ", e$message)
      })

      tryCatch({
        empirical_probs <- collect_empirical_category_probs(
          res$df_mplus,
          res$item_names_list,
          res$meta_long,
          ages,
          cfg$missing_code
        )
        model_probs <- compute_model_implied_category_probs(
          best_model,
          res$item_names_list,
          res$meta_long,
          ages,
          time_scores,
          is_joint = TRUE
        )
        sensitivity_out <- summarise_lcga_sensitivity(empirical_probs, model_probs, res$meta_long)
        sensitivity_summary <- sensitivity_out$summary
        sensitivity_tbl <- generate_lcga_sensitivity_table(
          fit_df = fit_df,
          best_k = best_k,
          class_assign_data = class_assign_data,
          sensitivity_summary = sensitivity_summary,
          analysis_label = lcga_spec$label
        )
        save_gt_html(sensitivity_tbl, file.path(paths$tables_dir, glue("07_Misspecification_Sensitivity_{lcga_spec$name}.html")))
      }, error = function(e) {
        message("    ! Sensitivity Table Error: ", e$message)
      })

      if (!is.null(class_assign_data) && nrow(class_assign_data) > 0) {
        tryCatch({
          overlap_metrics <- compute_class_overlap_metrics(class_assign_data)
          distance_df <- compute_lcga_estimate_distance(res, best_k, ages, time_scores)
          overlap_tbl <- build_class_overlap_summary_table(
            summary_df = overlap_metrics$summary,
            distance_df = distance_df,
            analysis_label = lcga_spec$label,
            selected_k = best_k
          )
          save_gt_html(overlap_tbl, file.path(paths$tables_dir, glue("06c_Class_Overlap_Summary_{lcga_spec$name}.html")))

          overlap_heatmap <- plot_posterior_overlap_heatmap(overlap_metrics$overlap)
          ggsave(
            file.path(paths$figures_dir, glue("03d_posterior_overlap_heatmap_{lcga_spec$name}.png")),
            overlap_heatmap,
            width = 14,
            height = max(7, 0.55 * best_k + 4)
          )

          overlap_margin_plot <- plot_posterior_margin_distribution(
            overlap_metrics$margins,
            color_map = color_map_vec
          )
          ggsave(
            file.path(paths$figures_dir, glue("03e_posterior_margin_distribution_{lcga_spec$name}.png")),
            overlap_margin_plot,
            width = 14,
            height = max(7, 0.45 * best_k + 4)
          )
        }, error = function(e) {
          message("    ! Class Overlap Error: ", e$message)
        })
      }

      tryCatch({
        checklist_tbl <- generate_lcga_grolts_checklist_table(
          analysis_label = lcga_spec$label,
          time_label = period_label,
          ages = ages,
          time_var = cfg$time_var,
          outcome_meta = res$meta_long,
          fit_df = fit_df,
          best_k = best_k,
          class_assign_data = class_assign_data,
          run_profile = res$run_profile,
          syntax_dir = syntax_dir,
          availability_available = !is.null(avail_tbl),
          descriptive_available = !is.null(desc_tbl),
          sensitivity_available = !is.null(sensitivity_tbl),
          estimates_available = !is.null(est_tbl),
          weighted_available = FALSE,
          final_plot_available = length(traj_out$plots) > 0,
          per_model_plots_available = length(model_grid_plots) > 0,
          combined_observed_plot_available = length(observed_plots) > 0,
          missingness_correlates_available = !is.null(missingness_tbl),
          numerical_summary_available = !is.null(final_solution_tbl),
          structural_sensitivity_audit = res$structural_sensitivity %||% NULL,
          custom_labels_applied = any(!grepl("^Class\\s+[0-9]+$", class_meta$Display_Label)),
          manual_override_key = manual_override_key
        )
        save_gt_html(checklist_tbl, file.path(paths$tables_dir, glue("08_GRoLTS_Checklist_{lcga_spec$name}.html")))
      }, error = function(e) {
        message("    ! GRoLTS Checklist Error: ", e$message)
      })

      tryCatch({
        smart_tbl <- generate_lca_smart_checklist_table(
          analysis_label = lcga_spec$label,
          indicator_vars = lcga_spec$outcome_prefixes,
          fit_summary_aug = fit_df,
          best_k = best_k,
          report_tables = list(
            availability = avail_tbl,
            descriptive = desc_tbl,
            fit = fit_tbl,
            avepp = avepp_tbl %||% NULL,
            estimates = est_tbl,
            weighted = NULL
          ),
          starts = res$run_profile$starts %||% cfg$starts,
          syntax_dir = syntax_dir,
          profile_plot_available = length(traj_out$plots) > 0,
          class_labels = class_meta$Display_Label,
          manual_override_key = manual_override_key
        )
        save_gt_html(smart_tbl, file.path(paths$tables_dir, glue("09_SMART_LCA_Checklist_{lcga_spec$name}.html")))
      }, error = function(e) {
        message("    ! SMART-LCA Checklist Error: ", e$message)
      })

      all_joint_reports[[paste(lcga_spec$name, time_label, sep = "_")]] <- list(
        spec = lcga_spec,
        time_label = time_label,
        best_k = best_k,
        results_path = rds_file,
        tables = list(
          availability = avail_tbl,
          descriptive = desc_tbl,
          fit = fit_tbl,
          avepp = avepp_tbl %||% NULL,
          estimates = est_tbl,
          final_solution = final_solution_tbl,
          missingness = missingness_tbl,
          structural_sensitivity = structural_sensitivity_tbl,
          overlap = overlap_tbl,
          sensitivity = sensitivity_tbl,
          grolts = checklist_tbl,
          smart_lca = smart_tbl
        ),
        diagnostics = list(
          sensitivity_summary = sensitivity_summary
        ),
        figures = list(
          trajectories = traj_out$plots,
          composite = p_comp,
          all_models = model_grid_plots,
          observed = observed_plots,
          overlap_heatmap = overlap_heatmap,
          overlap_margin = overlap_margin_plot
        )
      )
    }
  }

  saveRDS(
    all_joint_reports,
    trajectory_path("output", "objects", glue("all_joint_report_{cfg$analysis_stage}.rds"))
  )
  message("\n>>> Trajectory joint reporting complete.")
}

if (sys.nframe() == 0) {
  run_trajectory_joint_report()
}
