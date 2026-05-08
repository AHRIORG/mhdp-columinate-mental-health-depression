.apply_cutoffs <- function(df, cutoff_table,
                           apply_mode = c("none", "SEX", "AGEGRP", "SEX_AGEGRP", "hierarchical"),
                           hierarchy = c("SEX_AGEGRP", "SEX", "AGEGRP", "none")) {
  apply_mode <- match.arg(apply_mode)
  df$Depression_Binary <- NA_integer_
  df$Cutoff_Theta_Applied <- NA_real_
  df$Cutoff_Grouping_Used <- NA_character_
  df$Cutoff_Group_Level_Used <- NA_character_

  overall_idx <- cutoff_table$Grouping == "none" & cutoff_table$Group_Level == "Overall"
  overall_cut <- cutoff_table$Cutoff_Theta[overall_idx][1] %||% NA_real_

  get_cut <- function(grouping, level) {
    idx <- cutoff_table$Grouping == grouping & cutoff_table$Group_Level == level
    if (!any(idx)) {
      return(list(
        cutoff = NA_real_,
        grouping = NA_character_,
        level = NA_character_
      ))
    }
    list(
      cutoff = as.numeric(cutoff_table$Cutoff_Theta[idx][1]),
      grouping = grouping,
      level = as.character(level)
    )
  }

  overall_cut_info <- list(
    cutoff = as.numeric(overall_cut),
    grouping = "none",
    level = "Overall"
  )

  assign_cutoff <- function(idx, cut_info) {
    if (length(idx) == 0 || is.na(cut_info$cutoff)) {
      return()
    }

    df$Cutoff_Theta_Applied[idx] <<- cut_info$cutoff
    df$Cutoff_Grouping_Used[idx] <<- cut_info$grouping
    df$Cutoff_Group_Level_Used[idx] <<- cut_info$level
    df$Depression_Binary[idx] <<- ifelse(df$Theta_Harmonized[idx] >= cut_info$cutoff, 1L, 0L)
  }

  if (apply_mode == "none") {
    if (!is.na(overall_cut_info$cutoff)) {
      idx <- which(!is.na(df$Theta_Harmonized))
      assign_cutoff(idx, overall_cut_info)
    }
    return(df)
  }

  if (apply_mode == "hierarchical") {
    for (i in seq_len(nrow(df))) {
      th <- df$Theta_Harmonized[i]
      if (is.na(th)) {
        next
      }

      cut_info <- list(cutoff = NA_real_, grouping = NA_character_, level = NA_character_)
      for (grouping in hierarchy) {
        if (grouping == "none") {
          cut_info <- overall_cut_info
        } else if (grouping %in% names(df)) {
          level <- as.character(df[[grouping]][i])
          if (!is.na(level)) {
            cut_info <- get_cut(grouping, level)
          }
        }

        if (!is.na(cut_info$cutoff)) {
          break
        }
      }

      assign_cutoff(i, cut_info)
    }

    return(df)
  }

  grouping <- apply_mode
  if (!grouping %in% names(df)) {
    if (!is.na(overall_cut_info$cutoff)) {
      idx <- which(!is.na(df$Theta_Harmonized))
      assign_cutoff(idx, overall_cut_info)
    }
    return(df)
  }

  for (level in unique(stats::na.omit(df[[grouping]]))) {
    cut_info <- get_cut(grouping, as.character(level))
    if (is.na(cut_info$cutoff)) {
      cut_info <- overall_cut_info
    }

    idx <- which(df[[grouping]] == level & !is.na(df$Theta_Harmonized))
    assign_cutoff(idx, cut_info)
  }

  if (!is.na(overall_cut_info$cutoff)) {
    idx <- which(is.na(df$Depression_Binary) & !is.na(df$Theta_Harmonized))
    assign_cutoff(idx, overall_cut_info)
  }

  df
}

apply_engine_cutoff <- function(theta_data, engine_id = NULL, engine = NULL) {
  if (!"Theta_Harmonized" %in% names(theta_data)) {
    stop("`theta_data` must include a Theta_Harmonized column.", call. = FALSE)
  }

  engine <- engine %||% load_engine(engine_id)
  df <- .ensure_groups(theta_data)
  meta <- attr(engine, "engine_metadata")
  out <- .apply_cutoffs(
    df = df,
    cutoff_table = engine$cutoff_table,
    apply_mode = engine$cutoff_apply_mode %||% "none",
    hierarchy = engine$cutoff_apply_hierarchy %||% c("none")
  )

  out$Cutoff_Apply_Mode <- engine$cutoff_apply_mode %||% "none"
  out$Engine_ID <- attr(engine, "engine_id") %||% engine_id %||% NA_character_
  out$Engine_Label <- attr(engine, "engine_label") %||% out$Engine_ID
  out$Cutoff_Method <- if (!is.null(meta) && "cutoff_method" %in% names(meta)) {
    meta$cutoff_method[[1]]
  } else {
    NA_character_
  }

  out
}
