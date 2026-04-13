# ==============================================================================
# UTILITY: Model Selection & Fit Summaries
# PURPOSE: Functions to compare LCGA models and select the preferred solution.
# ==============================================================================

library(dplyr)
library(glue)
library(purrr)

get_min_class_pct <- function(model_obj) {
  if (!is.null(model_obj$savedata) && is.data.frame(model_obj$savedata)) {
    cols <- names(model_obj$savedata)
    class_col <- cols[toupper(cols) %in% c("C", "MLCC")]
    if (length(class_col) > 0) {
      cc <- model_obj$savedata[[class_col[1]]]
      cc <- cc[!is.na(cc)]
      if (length(cc) > 0) {
        p <- prop.table(table(cc))
        return(as.numeric(min(p)))
      }
    }
  }

  NA_real_
}

augment_fit_with_class_sizes <- function(fit_df, all_results) {
  key_col <- if ("model_id" %in% names(fit_df)) "model_id" else "model_name"

  class_metrics <- tibble::tibble(
    !!key_col := names(all_results),
    MinClassPct = purrr::map_dbl(all_results, get_min_class_pct)
  )

  left_join(fit_df, class_metrics, by = key_col)
}

select_lcga_model <- function(
  fit_df,
  preset = c("default", "bic_only", "bic_entropy", "bic_vlmr"),
  test = c("BLRT", "VLMR"),
  alpha = 0.05,
  min_class_pct = 0.05,
  bic_close_delta = 2
) {
  preset <- match.arg(preset)
  test <- match.arg(test)

  if (!"Classes" %in% names(fit_df) && "NLatentClasses" %in% names(fit_df)) {
    fit_df <- fit_df %>% mutate(Classes = NLatentClasses)
  }
  if (!"BLRT_p" %in% names(fit_df) && "BLRT_PValue" %in% names(fit_df)) {
    fit_df <- fit_df %>% mutate(BLRT_p = BLRT_PValue)
  }
  if (!"VLMR_p" %in% names(fit_df) && "T11_VLMR_PValue" %in% names(fit_df)) {
    fit_df <- fit_df %>% mutate(VLMR_p = T11_VLMR_PValue)
  }

  require_test <- preset %in% c("default", "bic_vlmr")
  prefer_entropy <- preset %in% c("default", "bic_entropy", "bic_vlmr")
  if (preset == "bic_vlmr") {
    test <- "VLMR"
  }
  p_col <- if (test == "BLRT") "BLRT_p" else "VLMR_p"

  needed <- c("Classes", "BIC", "Entropy")
  missing_needed <- setdiff(needed, names(fit_df))
  if (length(missing_needed) > 0) {
    stop(glue("Selection failed: fit_df missing columns: {paste(missing_needed, collapse = ', ')}"))
  }

  df <- fit_df %>%
    mutate(
      BIC = suppressWarnings(as.numeric(BIC)),
      Entropy = suppressWarnings(as.numeric(Entropy)),
      BLRT_p = suppressWarnings(as.numeric(BLRT_p)),
      VLMR_p = suppressWarnings(as.numeric(VLMR_p))
    )

  if ("MinClassPct" %in% names(df)) {
    df <- df %>% filter(is.na(MinClassPct) | MinClassPct >= min_class_pct)
  }

  if (require_test && p_col %in% names(df)) {
    df2 <- df %>% filter(!is.na(.data[[p_col]]) & .data[[p_col]] < alpha)
    if (nrow(df2) > 0) {
      df <- df2
    }
  }

  if (nrow(df) == 0) {
    stop("Selection failed: no candidate models available after filtering.")
  }

  df <- df %>% arrange(BIC)
  best0 <- df %>% slice(1)
  best <- best0

  if (prefer_entropy && is.finite(bic_close_delta) && bic_close_delta > 0) {
    cand <- df %>% filter(BIC <= best0$BIC + bic_close_delta)
    if (nrow(cand) > 1) {
      best <- cand %>% arrange(desc(Entropy), BIC) %>% slice(1)
    }
  }

  list(
    model_row = best,
    model_name = if ("model_name" %in% names(best)) best$model_name[[1]] else NA,
    model_id = if ("model_id" %in% names(best)) best$model_id[[1]] else NA,
    Classes = best$Classes[[1]],
    BIC = best$BIC[[1]],
    Entropy = best$Entropy[[1]],
    test = test,
    p_col = p_col,
    min_class_pct = if ("MinClassPct" %in% names(best)) best$MinClassPct[[1]] else NA_real_,
    bic_close_delta = bic_close_delta,
    alpha = alpha,
    filtered_table = df
  )
}
