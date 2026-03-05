#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(stats)
  library(utils)
  library(dplyr)
  library(tidyr)
  library(ggsci)
  library(ggplot2)
  library(lubridate)
  library(patchwork)
  library(broom)
  library(scales)
})

project_root <- "."
objects_dir <- file.path(project_root, "_tools", "_objects")
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)

private_adam_roots <- c(
  file.path("..", "..", "_private_use", "data_management", "data", "co-luminate", "adam"),
  file.path("..", "_private_use", "data_management", "data", "co-luminate", "adam"),
  file.path("..", "..", "..", "_private_use", "data_management", "data", "co-luminate", "adam")
)
multi_private_paths <- file.path(private_adam_roots, "dt_multistudies.RData")

load_rdata_df <- function(path, preferred_name = NULL) {
  if (!file.exists(path)) {
    return(NULL)
  }

  src <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = src)

  if (!is.null(preferred_name) && preferred_name %in% loaded_names) {
    preferred_obj <- get(preferred_name, envir = src)
    if (is.data.frame(preferred_obj)) {
      return(preferred_obj)
    }
  }

  for (nm in loaded_names) {
    obj <- get(nm, envir = src)
    if (is.data.frame(obj)) {
      return(obj)
    }
  }

  NULL
}

load_first_df <- function(paths, preferred_name = NULL) {
  unique_paths <- unique(paths)
  for (path in unique_paths) {
    candidate <- tryCatch(
      load_rdata_df(path, preferred_name = preferred_name),
      error = function(e) NULL
    )
    if (!is.null(candidate)) {
      return(list(
        data = candidate,
        source = normalizePath(path, mustWork = FALSE)
      ))
    }
  }
  list(data = NULL, source = NA_character_)
}

to_binary <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  yes_values <- c("yes", "1", "true", "positive", "depressed", "in school", "inschool")
  no_values <- c("no", "0", "false", "negative", "not depressed", "not in school", "not inschool")

  out <- rep(NA_real_, length(x_chr))
  out[x_chr %in% yes_values] <- 1
  out[x_chr %in% no_values] <- 0
  out
}

recode_govg_no_grant <- function(x) {
  # Primary mapping used in OBJ00 script.
  out <- suppressWarnings(
    factor(abs(as.numeric(x) - 2), levels = c(0, 1), labels = c("No", "Yes"))
  )

  if (all(is.na(out))) {
    x_chr <- tolower(trimws(as.character(x)))
    out_chr <- rep(NA_character_, length(x_chr))
    out_chr[x_chr %in% c("yes", "1", "true")] <- "No"
    out_chr[x_chr %in% c("no", "0", "false")] <- "Yes"
    out <- factor(out_chr, levels = c("No", "Yes"))
  }

  out
}

extract_year <- function(x) {
  if (inherits(x, "Date")) {
    return(as.numeric(format(x, "%Y")))
  }
  if (inherits(x, "POSIXt")) {
    return(as.numeric(format(as.Date(x), "%Y")))
  }

  x_chr <- as.character(x)
  y <- suppressWarnings(as.numeric(format(as.Date(x_chr), "%Y")))
  if (!all(is.na(y))) {
    return(y)
  }

  y2 <- suppressWarnings(as.numeric(substr(x_chr, 1, 4)))
  y2[!(y2 >= 1900 & y2 <= 2100)] <- NA_real_
  y2
}

format_or_ci <- function(est, lcl, ucl) {
  if (is.na(est) || is.na(lcl) || is.na(ucl)) {
    return("NA")
  }
  sprintf("%.2f (%.2f, %.2f)", est, lcl, ucl)
}

first_non_missing <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) {
    return(NA_real_)
  }
  y[[1]]
}

mediators <- c("DNKA", "ESXC", "FDSC", "GOVG", "VLNC", "SCHL", "SCSP")
mediator_labels <- c(
  DNKA = "Drank Alcohol",
  ESXC = "Ever Had Sex",
  FDSC = "Food Insecure",
  GOVG = "No Government Grant",
  VLNC = "Experienced Violence",
  SCHL = "In School",
  SCSP = "Has Social Support"
)

covariates <- c("SEX", "AGE", "HIVS", "VSTY", "EXTMG", "RURAL", "ORPH")
required_vars <- c("DPBN", "VISITDT", covariates[covariates != "VSTY"], mediators)

multi_loaded <- load_first_df(
  multi_private_paths,
  preferred_name = "dt_multistudies"
)
if (is.null(multi_loaded$data)) {
  stop("Unable to load private dt_multistudies.RData for mediation preliminary analysis.")
}
dt_multistudies <- multi_loaded$data
data_source_multi <- multi_loaded$source

missing_required <- setdiff(required_vars, names(dt_multistudies))
if (length(missing_required) > 0) {
  stop("Missing required columns in dt_multistudies: ", paste(missing_required, collapse = ", "))
}

analysis_data <- dt_multistudies |>
  mutate(
    GOVG = recode_govg_no_grant(GOVG),
    DPBN_bin = to_binary(DPBN),
    VSTY = extract_year(VISITDT)
  ) |>
  mutate(
    VSTY = if (all(is.na(VSTY))) NA_real_ else VSTY - min(VSTY, na.rm = TRUE)
  ) |>
  select(
    USUBJID, STUDY,
    DPBN_bin,
    SEX, AGE, VSTY, HIVS, RURAL, EXTMG, ORPH,
    all_of(mediators)
  ) |>
  filter(!is.na(DPBN_bin))

extract_term_row <- function(fit_tidy, variable) {
  out <- fit_tidy |>
    filter(grepl(paste0("^", variable), term)) |>
    arrange(p.value)

  if (nrow(out) == 0) {
    return(data.frame(
      variable = variable,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p.value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  out <- out[1, c("estimate", "conf.low", "conf.high", "p.value"), drop = FALSE]
  out$variable <- variable
  out[, c("variable", "estimate", "conf.low", "conf.high", "p.value"), drop = FALSE]
}

unadj_results <- lapply(mediators, function(med) {
  f <- as.formula(paste("DPBN_bin ~", med))
  fit <- glm(formula = f, family = binomial, data = analysis_data)
  fit_tidy <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  extract_term_row(fit_tidy, med)
})
unadj_results <- bind_rows(unadj_results) |>
  rename(
    unadj_estimate = estimate,
    unadj_conf.low = conf.low,
    unadj_conf.high = conf.high,
    unadj_p.value = p.value
  )

base_results <- lapply(mediators, function(med) {
  f <- as.formula(paste("DPBN_bin ~", med, "+ SEX + AGE + HIVS + VSTY + EXTMG + RURAL + ORPH"))
  fit <- glm(formula = f, family = binomial, data = analysis_data)
  fit_tidy <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  extract_term_row(fit_tidy, med)
})
base_results <- bind_rows(base_results) |>
  rename(
    base_estimate = estimate,
    base_conf.low = conf.low,
    base_conf.high = conf.high,
    base_p.value = p.value
  )

full_formula <- as.formula(
  paste(
    "DPBN_bin ~",
    paste(c(mediators, covariates), collapse = " + ")
  )
)
full_fit <- glm(formula = full_formula, family = binomial, data = analysis_data)
full_tidy <- tidy(full_fit, conf.int = TRUE, exponentiate = TRUE)
adj_results <- lapply(mediators, function(med) extract_term_row(full_tidy, med))
adj_results <- bind_rows(adj_results) |>
  rename(
    adj_estimate = estimate,
    adj_conf.low = conf.low,
    adj_conf.high = conf.high,
    adj_p.value = p.value
  )

summary_df <- data.frame(
  variable = mediators,
  label = unname(mediator_labels[mediators]),
  stringsAsFactors = FALSE
) |>
  left_join(unadj_results, by = "variable") |>
  left_join(base_results, by = "variable") |>
  left_join(adj_results, by = "variable") |>
  mutate(
    unadj_or_ci = mapply(format_or_ci, unadj_estimate, unadj_conf.low, unadj_conf.high),
    base_or_ci = mapply(format_or_ci, base_estimate, base_conf.low, base_conf.high),
    adj_or_ci = mapply(format_or_ci, adj_estimate, adj_conf.low, adj_conf.high)
  )

mediator_bin <- analysis_data |>
  transmute(across(all_of(mediators), to_binary))

pair_grid <- t(combn(mediators, 2))

pair_assoc_df <- lapply(seq_len(nrow(pair_grid)), function(i) {
  var_a <- pair_grid[i, 1]
  var_b <- pair_grid[i, 2]
  a <- mediator_bin[[var_a]]
  b <- mediator_bin[[var_b]]
  keep <- complete.cases(a, b)
  a <- a[keep]
  b <- b[keep]

  if (length(a) < 3 || length(unique(a)) < 2 || length(unique(b)) < 2) {
    return(data.frame(
      mediator_a = var_a,
      mediator_b = var_b,
      label_a = unname(mediator_labels[var_a]),
      label_b = unname(mediator_labels[var_b]),
      n_complete = length(a),
      or = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p_value = NA_real_,
      cramer_v = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  tab <- table(factor(a, levels = c(0, 1)), factor(b, levels = c(0, 1)))
  ft <- fisher.test(tab)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  n <- sum(tab)

  data.frame(
    mediator_a = var_a,
    mediator_b = var_b,
    label_a = unname(mediator_labels[var_a]),
    label_b = unname(mediator_labels[var_b]),
    n_complete = n,
    or = unname(ft$estimate),
    conf.low = ft$conf.int[[1]],
    conf.high = ft$conf.int[[2]],
    p_value = ft$p.value,
    cramer_v = sqrt(as.numeric(chi$statistic) / n),
    stringsAsFactors = FALSE
  )
})
pair_assoc_df <- bind_rows(pair_assoc_df) |>
  arrange(desc(cramer_v), p_value)

pair_cor_df <- lapply(seq_len(nrow(pair_grid)), function(i) {
  var_a <- pair_grid[i, 1]
  var_b <- pair_grid[i, 2]
  a <- mediator_bin[[var_a]]
  b <- mediator_bin[[var_b]]
  keep <- complete.cases(a, b)
  a <- a[keep]
  b <- b[keep]

  if (length(a) < 3 || length(unique(a)) < 2 || length(unique(b)) < 2) {
    return(data.frame(
      mediator_a = var_a,
      mediator_b = var_b,
      label_a = unname(mediator_labels[var_a]),
      label_b = unname(mediator_labels[var_b]),
      n_complete = length(a),
      correlation = NA_real_,
      p_value = NA_real_,
      abs_correlation = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  ct <- suppressWarnings(cor.test(a, b, method = "pearson"))
  data.frame(
    mediator_a = var_a,
    mediator_b = var_b,
    label_a = unname(mediator_labels[var_a]),
    label_b = unname(mediator_labels[var_b]),
    n_complete = length(a),
    correlation = unname(ct$estimate),
    p_value = ct$p.value,
    abs_correlation = abs(unname(ct$estimate)),
    stringsAsFactors = FALSE
  )
})
pair_cor_df <- bind_rows(pair_cor_df) |>
  arrange(desc(abs_correlation), p_value)

pair_combined_df <- pair_assoc_df |>
  rename(
    n_assoc = n_complete,
    or_p_value = p_value
  ) |>
  full_join(
    pair_cor_df |>
      rename(
        n_corr = n_complete,
        rho = correlation,
        rho_p_value = p_value
      ) |>
      select(mediator_a, mediator_b, n_corr, rho, rho_p_value, abs_correlation),
    by = c("mediator_a", "mediator_b")
  ) |>
  mutate(
    label_a = unname(mediator_labels[mediator_a]),
    label_b = unname(mediator_labels[mediator_b]),
    n_complete = dplyr::coalesce(n_assoc, n_corr)
  ) |>
  select(
    mediator_a, mediator_b, label_a, label_b, n_complete,
    or, conf.low, conf.high, or_p_value, cramer_v,
    rho, rho_p_value, abs_correlation
  ) |>
  arrange(desc(abs_correlation), desc(cramer_v), or_p_value, rho_p_value)

corr_matrix <- suppressWarnings(
  cor(mediator_bin, use = "pairwise.complete.obs", method = "pearson")
)
dimnames(corr_matrix) <- list(
  unname(mediator_labels[colnames(corr_matrix)]),
  unname(mediator_labels[colnames(corr_matrix)])
)

or_matrix <- matrix(
  NA_real_,
  nrow = length(mediators),
  ncol = length(mediators),
  dimnames = list(
    unname(mediator_labels[mediators]),
    unname(mediator_labels[mediators])
  )
)
for (i in seq_len(nrow(pair_assoc_df))) {
  a_idx <- match(pair_assoc_df$mediator_a[[i]], mediators)
  b_idx <- match(pair_assoc_df$mediator_b[[i]], mediators)
  if (is.na(a_idx) || is.na(b_idx)) {
    next
  }
  or_val <- pair_assoc_df$or[[i]]
  or_matrix[a_idx, b_idx] <- or_val
  or_matrix[b_idx, a_idx] <- or_val
}

forest_df <- bind_rows(
  summary_df |>
    transmute(
      variable, label,
      model = "Unadjusted",
      estimate = unadj_estimate,
      conf.low = unadj_conf.low,
      conf.high = unadj_conf.high
    ),
  summary_df |>
    transmute(
      variable, label,
      model = "Confounders Adjusted",
      estimate = base_estimate,
      conf.low = base_conf.low,
      conf.high = base_conf.high
    ),
  summary_df |>
    transmute(
      variable, label,
      model = "Fully Adjusted",
      estimate = adj_estimate,
      conf.low = adj_conf.low,
      conf.high = adj_conf.high
    )
) |>
  filter(!is.na(estimate))

vars_order <- summary_df |>
  arrange(adj_estimate) |>
  pull(variable)
label_order <- unname(mediator_labels[vars_order])
label_order <- unique(label_order[!is.na(label_order)])

forest_df <- forest_df |>
  mutate(
    label = factor(label, levels = label_order),
    model = factor(model, levels = c(
      "Fully Adjusted",
      "Confounders Adjusted",
      "Unadjusted"
    ))
  ) |>
  arrange(label)

forest_plot <- ggplot(
  forest_df,
  aes(x = estimate, y = label, color = model, group = model)
) +
  geom_vline(xintercept = 1, linetype = "solid", color = "gray40", alpha = 0.5) +
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    position = position_dodge(width = 0.7),
    width = 0.2,
    linewidth = 0.8,
    orientation = "y",
    alpha = 0.8
  ) +
  geom_point(
    position = position_dodge(width = 0.7),
    size = 2,
    shape = 21,
    fill = "white",
    stroke = 1.2
  ) +
  geom_text(
    aes(x = -2.95, label = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)),
    position = position_dodge(width = 0.7),
    hjust = 0,
    size = 3.3,
    show.legend = FALSE,
    color = "black"
  ) +
  annotate(
    "text",
    x = -2.95,
    y = length(levels(forest_df$label)) + 1,
    label = "OR (95% CI)",
    fontface = "bold",
    hjust = 0,
    size = 3.8
  ) +
  scale_x_continuous(breaks = c(0, 1, 2, 3), limits = c(-3, 3)) +
  scale_color_jama() +
  labs(
    title = "Forest Plot of Potential Mediators",
    subtitle = "Odds Ratios for Depression (95% CI)",
    x = "Odds Ratio",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray30", size = 11),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  ) +
  coord_cartesian(clip = "off")

# Match OBJ00 panel structure for "B" panel (mediator configurations by pattern).
panel_mediators <- c("FDSC", "GOVG", "ESXC", "DNKA", "VLNC")
panel_order <- vars_order[vars_order %in% panel_mediators]
if (length(panel_order) == 0) {
  panel_order <- panel_mediators
}
panel_label_levels <- unname(mediator_labels[panel_order])

pattern_source <- analysis_data |>
  mutate(across(all_of(panel_mediators), to_binary)) |>
  mutate(
    pattern = apply(
      as.data.frame(pick(all_of(panel_mediators))),
      1,
      function(x) paste0(ifelse(is.na(x), "M", as.integer(x)), collapse = "")
    ),
    nfactor = rowSums(!is.na(pick(all_of(panel_mediators))))
  )

pattern_df <- pattern_source |>
  group_by(pattern) |>
  summarise(
    n = n(),
    npat = n_distinct(USUBJID),
    depressed_n = sum(DPBN_bin == 1, na.rm = TRUE),
    prevalence = depressed_n / n,
    nfactor = max(nfactor, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(mpat = grepl("M", pattern)) |>
  filter(depressed_n > 0) |>
  arrange(mpat, desc(npat), desc(prevalence), desc(nfactor), pattern)

# if (nrow(pattern_df) > 20) {
#   pattern_df <- pattern_df[1:20, , drop = FALSE]
# }
pattern_df <- pattern_df |>
  mutate(pattern_id = factor(seq_len(nrow(pattern_df))))

pattern_lookup <- pattern_df |>
  select(pattern, pattern_id)

dt_upset <- pattern_source |>
  left_join(pattern_lookup, by = "pattern") |>
  filter(!is.na(pattern_id)) |>
  select(
    USUBJID, STUDY, DPBN_bin, pattern, pattern_id, nfactor,
    all_of(panel_mediators),
    RURAL, HIVS, EXTMG, SCSP, SCHL, SEX
  )

pattern_long <- dt_upset |>
  select(pattern_id, all_of(panel_mediators)) |>
  pivot_longer(cols = all_of(panel_mediators), names_to = "variable", values_to = "status_raw") |>
  mutate(
    status = factor(
      ifelse(is.na(status_raw), "Missing", ifelse(status_raw == 1, "Yes", "No")),
      levels = c("Yes", "No", "Missing")
    ),
    mediators = factor(unname(mediator_labels[variable]), levels = panel_label_levels),
    mediator_num = as.numeric(mediators)
  )

dt_segs <- pattern_long |>
  filter(status == "Yes") |>
  group_by(pattern_id) |>
  summarise(
    ymin = min(mediator_num, na.rm = TRUE),
    ymax = max(mediator_num, na.rm = TRUE),
    .groups = "drop"
  )

char_vars <- c("RURAL", "HIVS", "EXTMG", "SCSP", "SCHL", "SEX")
char_labels <- c(
  RURAL = "Reside in Rural Area",
  HIVS = "HIV Positive",
  EXTMG = "External Migration",
  SCSP = "Has Social Support",
  SCHL = "In School",
  SEX = "Females"
)

to_char_yes <- function(variable, value) {
  if (variable == "SEX") {
    x_chr <- tolower(trimws(as.character(value)))
    ifelse(
      x_chr %in% c("female", "f", "2", "woman", "women", "girl"),
      1,
      ifelse(
        x_chr %in% c("male", "m", "1", "man", "men", "boy"),
        0,
        NA_real_
      )
    )
  } else {
    to_binary(value)
  }
}

chars_long <- dt_upset |>
  select(pattern_id, all_of(char_vars)) |>
  pivot_longer(cols = all_of(char_vars), names_to = "variable", values_to = "value") |>
  mutate(yes = mapply(to_char_yes, variable, value)) |>
  group_by(pattern_id, variable) |>
  summarise(
    vals = if (all(is.na(yes))) NA_real_ else mean(yes, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(charct = unname(char_labels[variable])) |>
  group_by(charct) |>
  mutate(char_order = sum(vals, na.rm = TRUE)) |>
  ungroup() |>
  arrange(desc(char_order))
chars_long$charct <- factor(chars_long$charct, levels = unique(chars_long$charct))

p_chars <- ggplot(chars_long, aes(x = pattern_id, y = charct, fill = vals)) +
  geom_tile(color = "white") +
  theme_void() +
  scale_fill_gradientn(
    colours = c("blue", "dodgerblue", "grey89", "orange", "red"),
    rescaler = ~ scales::rescale_mid(.x, mid = 0.5),
    limits = c(0, 1),
    labels = scales::percent
  ) +
  theme(
    axis.text.y = element_text(hjust = 1),
    legend.position = "bottom",
    legend.key.width = grid::unit(1, "cm")
  )
mid_pat<-pattern_df$pattern_id[median(as.numeric(pattern_df$pattern_id))]
p_top <- ggplot(pattern_df, aes(x = pattern_id, y = npat)) +
  geom_vline(xintercept = pattern_df$pattern_id,linewidth=0.25,linetype='dashed',color="gray85") +
  geom_label(
    data = pattern_df[which.max(pattern_df$npat), , drop = FALSE],
    aes(x = mid_pat, 
        y = (npat), label = "Number of Participants"),
    inherit.aes = FALSE,
    size = 5, vjust = 1, hjust = 0.5
  ) +
  geom_col(fill = "black") +
  geom_text(aes(label = npat), angle = 90, vjust = 0.5, hjust = 1, size = 3,color="white") +
  theme_void() +
  theme(axis.text.y = element_text(hjust = 1))

p_middle <- ggplot(pattern_df, aes(x = pattern_id, y = -prevalence)) +
  geom_vline(xintercept = pattern_df$pattern_id,linewidth=0.25,linetype='dashed',color="gray85") +
  geom_label(
    data = pattern_df[which.max(pattern_df$prevalence), , drop = FALSE],
    aes(x = mid_pat, y = -(prevalence - 0.03), label = "Prevalence of Depression"),
    inherit.aes = FALSE,
    size = 5
  ) +
  geom_col(aes(fill = "Depressed")) +
  theme_void() +
  scale_y_continuous(labels = function(x) scales::percent(abs(x))) +
  scale_fill_manual(values = c("Depressed" = "purple")) +
  theme(axis.text.y = element_text(hjust = 1)) +
  guides(fill = "none")

p_bottom <- ggplot(pattern_long, aes(x = pattern_id, y = mediator_num)) +
  geom_segment(
    data = dt_segs,
    aes(x = pattern_id, xend = pattern_id, y = ymin, yend = ymax),
    linewidth = 2,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_point(color = "white", size = 5)+
  geom_point(data=pattern_long |> filter(status!="Missing"),
             aes(color = status), size = 3) +
  geom_point(data=pattern_long |> filter(status=="Missing"),
             aes(shape = status), size = 3) +
  scale_y_continuous(
    breaks = seq_along(panel_label_levels),
    labels = panel_label_levels
  ) +
  theme_void() +
  scale_color_manual(values = c("Yes" = "black", "No" = "grey50", "Missing" = "gray95")) +
  scale_shape_manual(values = c("Missing" = 13)) +
  theme(
    axis.text.y = element_text(hjust = 1, face = "bold"),
    legend.justification = "left"
  )
p_bottom
combined_caption <- paste(
  "Unadjusted model is a bivariate association.",
  "Confounder adjusted model adjust for Sex, Age, HIV Positive status, Year of Visit, External migration, Rural residence, and orphanhood.",
  "Fully adjusted model includes all mediators simultaneously in addition to the confounder covariates.",
  "Panel B excludes potential mediators not associated with depression and shows combinations where there is at least one depressed individual.",
  sep = "\n"
)

combined_plot <- (
  (forest_plot + labs(title = NULL, subtitle = NULL, caption = combined_caption)) |
    (p_chars / p_top / p_middle / p_bottom)
) +
  plot_layout(guides = "collect", widths = c(1, 3)) +
  plot_annotation(tag_levels = list(c("A", "B", "", ""))) &
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.title = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot"
  )

sig_adj <- summary_df |>
  filter(!is.na(adj_p.value), adj_p.value < 0.05) |>
  arrange(desc(adj_estimate))
sig_text <- if (nrow(sig_adj) == 0) {
  "No mediator retained statistical significance in the fully adjusted model."
} else {
  paste(
    apply(sig_adj, 1, function(r) {
      paste0(
        r[["label"]], " (aOR ", format_or_ci(
          as.numeric(r[["adj_estimate"]]),
          as.numeric(r[["adj_conf.low"]]),
          as.numeric(r[["adj_conf.high"]])
        ), ")"
      )
    }),
    collapse = "; "
  )
}

nonsig_adj <- summary_df |>
  filter(is.na(adj_p.value) | adj_p.value >= 0.05) |>
  pull(label)
nonsig_text <- if (length(nonsig_adj) == 0) {
  "All candidate mediators were significant in the fully adjusted model."
} else {
  paste("Non-significant in the fully adjusted model:", paste(nonsig_adj, collapse = ", "), ".")
}

interpretation <- c(
  paste(
    "Fully adjusted associations indicate the strongest independent mediator candidates are:",
    sig_text
  ),
  paste(
    nonsig_text,
    "These results support prioritizing robust signals in mediation pathway testing while retaining weaker variables as supporting covariates where scientifically justified."
  )
)

plot_file <- file.path(objects_dir, "mediation_prelim_combined_plot.png")
ggsave(
  filename = plot_file,
  plot = combined_plot,
  width = 17,
  height = 9,
  dpi = 300
)

saveRDS(summary_df, file.path(objects_dir, "mediation_prelim_summary.rds"))
saveRDS(forest_df, file.path(objects_dir, "mediation_prelim_forest_data.rds"))
saveRDS(pattern_df, file.path(objects_dir, "mediation_prelim_pattern_summary.rds"))
saveRDS(pair_assoc_df, file.path(objects_dir, "mediation_prelim_pair_assoc.rds"))
saveRDS(pair_cor_df, file.path(objects_dir, "mediation_prelim_pair_cor.rds"))
saveRDS(pair_combined_df, file.path(objects_dir, "mediation_prelim_pair_combined.rds"))
saveRDS(corr_matrix, file.path(objects_dir, "mediation_prelim_corr_matrix.rds"))
saveRDS(or_matrix, file.path(objects_dir, "mediation_prelim_or_matrix.rds"))
saveRDS(interpretation, file.path(objects_dir, "mediation_prelim_interpretation.rds"))
saveRDS(
  list(
    source = data_source_multi,
    rendered_at = as.character(Sys.time()),
    n_analysis = nrow(analysis_data),
    n_patterns = nrow(pattern_df),
    n_pairs = nrow(pair_combined_df)
  ),
  file.path(objects_dir, "mediation_prelim_build_info.rds")
)

save(summary_df, file = file.path(objects_dir, "mediation_prelim_summary.RData"))
save(forest_df, file = file.path(objects_dir, "mediation_prelim_forest_data.RData"))
save(pattern_df, file = file.path(objects_dir, "mediation_prelim_pattern_summary.RData"))
save(pair_assoc_df, file = file.path(objects_dir, "mediation_prelim_pair_assoc.RData"))
save(pair_cor_df, file = file.path(objects_dir, "mediation_prelim_pair_cor.RData"))
save(pair_combined_df, file = file.path(objects_dir, "mediation_prelim_pair_combined.RData"))
save(corr_matrix, file = file.path(objects_dir, "mediation_prelim_corr_matrix.RData"))
save(or_matrix, file = file.path(objects_dir, "mediation_prelim_or_matrix.RData"))
save(interpretation, file = file.path(objects_dir, "mediation_prelim_interpretation.RData"))

cat("Created mediation preliminary objects in:", objects_dir, "\n")
cat("Source - dt_multistudies:", data_source_multi, "\n")
cat("Analysis rows used (no row-level data saved):", nrow(analysis_data), "\n")
cat("Mediators summarized:", nrow(summary_df), "\n")
cat("Pattern rows summarized:", nrow(pattern_df), "\n")
cat("Pairwise mediator rows (combined):", nrow(pair_combined_df), "\n")
