# ==============================================================================
# UTILITY: Mplus LCGA Syntax Generation
# PURPOSE: Functions to generate Mplus input blocks for univariate and joint
#          latent class growth analysis models.
# ==============================================================================

library(glue)
library(purrr)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

mplus_wrap_tokens <- function(tokens, prefix, width = 88, indent = 2, end = ";") {
  tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
  if (length(tokens) == 0) {
    return(glue("{prefix} {end}"))
  }

  cont_prefix <- paste0(strrep(" ", indent))
  lines <- character(0)
  current <- paste(prefix, tokens[1])

  if (nchar(current) > width) {
    lines <- c(lines, prefix)
    current <- paste(cont_prefix, tokens[1])
  }

  if (length(tokens) > 1) {
    for (tok in tokens[-1]) {
      candidate <- paste(current, tok)
      if (nchar(candidate) > width) {
        lines <- c(lines, current)
        current <- paste(cont_prefix, tok)
      } else {
        current <- candidate
      }
    }
  }

  lines <- c(lines, current)
  lines[length(lines)] <- paste0(lines[length(lines)], end)
  paste(lines, collapse = "\n")
}

mplus_growth_statement_factors <- function(item_names, time_scores, factors = c("i", "s", "q"), width = 88) {
  prefix <- glue("  {factors[1]} {factors[2]} {factors[3]} |")
  toks <- paste0(item_names, "@", time_scores)
  mplus_wrap_tokens(toks, prefix = prefix, width = width, indent = 4, end = ";")
}

mplus_series_statement <- function(tokens, width = 88) {
  mplus_wrap_tokens(tokens, prefix = "  SERIES =", width = width, indent = 4, end = ";")
}

normalise_lcga_structure_profile <- function(profile = NULL) {
  profile <- profile %||% list(
    id = "strict_lcga",
    label = "Strict LCGA",
    description = "All within-class growth-factor variances and covariances fixed to zero.",
    free_variances = character(0),
    free_covariances = list(),
    variance_start_values = c(i = 0.25, s = 0.05, q = 0.01),
    covariance_start_values = c(i_s = 0.00, i_q = 0.00, q_s = 0.00)
  )

  profile$id <- profile$id %||% "strict_lcga"
  profile$label <- profile$label %||% profile$id
  profile$description <- profile$description %||% profile$label
  profile$free_variances <- unique(tolower(as.character(profile$free_variances %||% character(0))))

  raw_covs <- profile$free_covariances %||% list()
  if (is.atomic(raw_covs) && !is.list(raw_covs)) {
    raw_covs <- list(raw_covs)
  }
  profile$free_covariances <- lapply(raw_covs, function(pair) {
    sort(tolower(as.character(pair)))
  })

  profile$variance_start_values <- profile$variance_start_values %||% c(i = 0.25, s = 0.05, q = 0.01)
  profile$covariance_start_values <- profile$covariance_start_values %||% c(i_s = 0.00, i_q = 0.00, q_s = 0.00)
  profile
}

get_factor_family <- function(factor_name) {
  substr(tolower(trimws(factor_name)), 1, 1)
}

render_free_or_fixed_variance_line <- function(factor_name, structure_profile) {
  structure_profile <- normalise_lcga_structure_profile(structure_profile)
  family <- get_factor_family(factor_name)

  if (family %in% structure_profile$free_variances) {
    start_value <- unname(structure_profile$variance_start_values[[family]] %||% 0.10)
    return(glue("  {factor_name}*{format(start_value, trim = TRUE, scientific = FALSE)};"))
  }

  glue("  {factor_name}@0;")
}

is_free_covariance_pair <- function(factor_a, factor_b, structure_profile) {
  pair_family <- sort(c(get_factor_family(factor_a), get_factor_family(factor_b)))
  free_pairs <- structure_profile$free_covariances %||% list()
  any(vapply(free_pairs, function(x) identical(sort(x), pair_family), logical(1)))
}

render_free_or_fixed_covariance_line <- function(factor_a, factor_b, structure_profile) {
  structure_profile <- normalise_lcga_structure_profile(structure_profile)
  pair_key <- paste(sort(c(get_factor_family(factor_a), get_factor_family(factor_b))), collapse = "_")

  if (is_free_covariance_pair(factor_a, factor_b, structure_profile)) {
    start_value <- unname(structure_profile$covariance_start_values[[pair_key]] %||% 0.00)
    return(glue("  {factor_a} WITH {factor_b}*{format(start_value, trim = TRUE, scientific = FALSE)};"))
  }

  glue("  {factor_a} WITH {factor_b}@0;")
}

build_growth_structure_lines <- function(factor_groups, structure_profile = NULL) {
  structure_profile <- normalise_lcga_structure_profile(structure_profile)
  factor_groups <- factor_groups[!vapply(factor_groups, is.null, logical(1))]
  all_factors <- unlist(factor_groups, use.names = FALSE)

  variance_lines <- vapply(all_factors, render_free_or_fixed_variance_line, character(1), structure_profile = structure_profile)

  covariance_lines <- character(0)
  for (group in factor_groups) {
    if (length(group) < 2) {
      next
    }

    for (a in seq_len(length(group) - 1)) {
      for (b in (a + 1):length(group)) {
        covariance_lines <- c(
          covariance_lines,
          render_free_or_fixed_covariance_line(group[a], group[b], structure_profile = structure_profile)
        )
      }
    }
  }

  if (length(factor_groups) > 1) {
    for (g1 in seq_len(length(factor_groups) - 1)) {
      for (g2 in (g1 + 1):length(factor_groups)) {
        for (factor_a in factor_groups[[g1]]) {
          for (factor_b in factor_groups[[g2]]) {
            covariance_lines <- c(covariance_lines, glue("  {factor_a} WITH {factor_b}@0;"))
          }
        }
      }
    }
  }

  paste(c(variance_lines, covariance_lines), collapse = "\n")
}

get_univariate_class_blocks <- function(item_names, time_scores, k, structure_profile = NULL) {
  structure_profile <- normalise_lcga_structure_profile(structure_profile)
  growth_spec <- paste0(item_names, "@", time_scores)
  growth_wrapped <- paste(
    strwrap(paste(growth_spec, collapse = " "), width = 75, prefix = "    ", initial = "    "),
    collapse = "\n"
  )

  base <- glue("%OVERALL%\n  i s q | \n{growth_wrapped};\n\n")
  mean_line <- "  [i s q];"
  structure_lines <- build_growth_structure_lines(list(c("i", "s", "q")), structure_profile = structure_profile)

  if (k == 1) {
    return(glue("{base}{mean_line}\n{structure_lines}\n"))
  }

  class_parts <- purrr::map_chr(seq_len(k), function(cc) {
    glue("%c#{cc}%\n{mean_line}\n{structure_lines}\n")
  }) %>%
    paste(collapse = "\n")

  glue("{base}{class_parts}")
}

get_joint_class_blocks <- function(item_names_list, time_scores, k, structure_profile = NULL) {
  structure_profile <- normalise_lcga_structure_profile(structure_profile)
  n_out <- length(item_names_list)
  factor_triplets <- purrr::map(seq_len(n_out), ~ paste0(c("i", "s", "q"), .x))
  all_factors <- unlist(factor_triplets, use.names = FALSE)

  overall_growth <- purrr::map_chr(seq_along(item_names_list), function(j) {
    mplus_growth_statement_factors(item_names_list[[j]], time_scores, factors = factor_triplets[[j]])
  }) %>%
    paste(collapse = "\n")

  mean_line <- sub("^  \\[", "  [", mplus_wrap_tokens(all_factors, prefix = "  [", width = 88, indent = 2, end = "];"))
  structure_lines <- build_growth_structure_lines(factor_triplets, structure_profile = structure_profile)
  base <- glue("%OVERALL%\n{overall_growth}\n\n")

  if (k == 1) {
    return(glue("{base}{mean_line}\n{structure_lines}\n"))
  }

  class_parts <- purrr::map_chr(seq_len(k), function(cc) {
    glue("%c#{cc}%\n{mean_line}\n{structure_lines}\n")
  }) %>%
    paste(collapse = "\n")

  glue("{base}{class_parts}")
}
