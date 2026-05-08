.normalize_binary_response <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    raw <- trimws(x)
    low <- tolower(raw)

    mapped <- ifelse(
      low %in% c("", "na", "n/a", "missing"),
      NA_character_,
      ifelse(
        low %in% c("yes", "y", "true", "1"),
        "1",
        ifelse(low %in% c("no", "n", "false", "0"), "0", raw)
      )
    )

    x <- suppressWarnings(as.numeric(mapped))
    return(x)
  }

  suppressWarnings(as.numeric(x))
}

validate_ssq10_data <- function(data, engine_id = NULL, engine = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  engine <- engine %||% load_engine(engine_id)
  required_items <- engine$items
  cleaned <- data

  missing_columns <- setdiff(required_items, names(cleaned))
  issues <- character()

  if (length(missing_columns) > 0) {
    issues <- c(
      issues,
      paste("Missing required SSQ-10 item columns:", paste(missing_columns, collapse = ", "))
    )
  }

  present_items <- intersect(required_items, names(cleaned))
  if (length(present_items) > 0) {
    for (item in present_items) {
      cleaned[[item]] <- .normalize_binary_response(cleaned[[item]])
      bad_idx <- which(!is.na(cleaned[[item]]) & !cleaned[[item]] %in% c(0, 1))
      if (length(bad_idx) > 0) {
        issues <- c(
          issues,
          paste0(
            item, " contains values outside the supported 0/1 range in row(s): ",
            paste(utils::head(bad_idx, 10), collapse = ", ")
          )
        )
      }
    }
  }

  cleaned <- .ensure_groups(cleaned)

  items_answered <- integer(nrow(cleaned))
  all_items_missing <- logical(nrow(cleaned))
  if (length(present_items) > 0) {
    items_answered <- rowSums(!is.na(cleaned[, present_items, drop = FALSE]))
    all_items_missing <- items_answered == 0
  }

  if (any(all_items_missing)) {
    issues <- c(
      issues,
      paste(
        "Rows with all SSQ-10 items missing:",
        paste(utils::head(which(all_items_missing), 10), collapse = ", ")
      )
    )
  }

  list(
    valid = length(issues) == 0,
    issues = unique(issues),
    data = cleaned,
    required_items = required_items,
    items_answered = items_answered,
    all_items_missing = all_items_missing
  )
}
