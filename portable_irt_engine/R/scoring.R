.template_scoring_data <- function(engine) {
  n_items <- length(engine$items)
  template_data <- as.data.frame(
    rbind(
      rep(0, n_items),
      rep(1, n_items),
      rep(c(0, 1), length.out = n_items),
      rep(c(1, 0), length.out = n_items)
    ),
    stringsAsFactors = FALSE
  )
  names(template_data) <- engine$items
  template_data
}

.rebuild_fixed_parameter_model <- function(engine) {
  template_data <- .template_scoring_data(engine)
  n_items <- ncol(template_data)
  item_types_vec <- rep(engine$model_type, n_items)
  n_factors <- engine$n_factors %||% 1L

  template_pars <- mirt::mirt(
    template_data,
    n_factors,
    itemtype = item_types_vec,
    pars = "values",
    verbose = FALSE
  )

  saved_pars <- engine$parameters

  for (i in seq_len(nrow(template_pars))) {
    item_name <- template_pars$item[i]
    par_name <- template_pars$name[i]
    if (identical(item_name, "GROUP")) {
      next
    }

    idx <- saved_pars$item == item_name & saved_pars$name == par_name
    if (sum(idx) == 1) {
      template_pars$value[i] <- saved_pars$value[idx]
      template_pars$est[i] <- FALSE
    }
  }

  mirt::mirt(
    template_data,
    n_factors,
    itemtype = item_types_vec,
    pars = template_pars,
    verbose = FALSE,
    calcNull = FALSE
  )
}

.map_theta_outputs <- function(theta, engine) {
  out <- data.frame(
    Theta_Harmonized = theta,
    PHQ_Expected_Sum = NA_real_,
    PHQ_Prob_GE10 = NA_real_
  )

  if (!is.null(engine$phq_theta_map) && length(theta) > 0) {
    tm <- engine$phq_theta_map
    out$PHQ_Expected_Sum <- stats::approx(tm$Theta, tm$PHQ_Expected_Sum, xout = theta, rule = 2)$y
    out$PHQ_Prob_GE10 <- stats::approx(tm$Theta, tm$PHQ_Prob_GE10, xout = theta, rule = 2)$y
  }

  out
}

score_dataset <- function(data, engine_id = NULL, engine = NULL) {
  engine <- engine %||% load_engine(engine_id)
  validation <- validate_ssq10_data(data, engine = engine)

  if (!validation$valid) {
    stop(paste(validation$issues, collapse = "\n"), call. = FALSE)
  }

  cleaned <- validation$data
  dat_score <- cleaned[, engine$items, drop = FALSE]
  keep_rows <- rowSums(!is.na(dat_score)) > 0
  total_ssq_score <- rowSums(dat_score, na.rm = TRUE)
  total_ssq_score[validation$items_answered == 0] <- NA_real_

  scored <- data.frame(
    SSQ10_Total_Score = total_ssq_score,
    Theta_Harmonized = rep(NA_real_, nrow(cleaned)),
    Depression_Binary = rep(NA_integer_, nrow(cleaned)),
    PHQ_Expected_Sum = rep(NA_real_, nrow(cleaned)),
    PHQ_Prob_GE10 = rep(NA_real_, nrow(cleaned)),
    Cutoff_Theta_Applied = rep(NA_real_, nrow(cleaned)),
    Cutoff_Grouping_Used = rep(NA_character_, nrow(cleaned)),
    Cutoff_Group_Level_Used = rep(NA_character_, nrow(cleaned)),
    Cutoff_Apply_Mode = rep(NA_character_, nrow(cleaned)),
    Cutoff_Method = rep(NA_character_, nrow(cleaned)),
    Engine_ID = rep(NA_character_, nrow(cleaned)),
    Engine_Label = rep(NA_character_, nrow(cleaned)),
    Items_Answered = validation$items_answered,
    stringsAsFactors = FALSE
  )

  if (!any(keep_rows)) {
    return(cbind(cleaned, scored))
  }

  dat_score_clean <- dat_score[keep_rows, , drop = FALSE]
  mod_fixed <- .rebuild_fixed_parameter_model(engine)

  nfact_fixed <- tryCatch(mod_fixed@Model$nfact, error = function(e) NA_integer_)
  use_qmc <- !is.na(nfact_fixed) && nfact_fixed >= 3
  theta <- mirt::fscores(
    mod_fixed,
    method = "EAP",
    full.scores = TRUE,
    QMC = use_qmc,
    response.pattern = dat_score_clean
  )[, 1]

  mapped <- .map_theta_outputs(theta, engine)
  scored$Theta_Harmonized[keep_rows] <- mapped$Theta_Harmonized
  scored$PHQ_Expected_Sum[keep_rows] <- mapped$PHQ_Expected_Sum
  scored$PHQ_Prob_GE10[keep_rows] <- mapped$PHQ_Prob_GE10

  tmp <- cleaned
  tmp$Theta_Harmonized <- scored$Theta_Harmonized
  tmp <- apply_engine_cutoff(tmp, engine = engine)
  scored$Depression_Binary <- tmp$Depression_Binary
  scored$Cutoff_Theta_Applied <- tmp$Cutoff_Theta_Applied
  scored$Cutoff_Grouping_Used <- tmp$Cutoff_Grouping_Used
  scored$Cutoff_Group_Level_Used <- tmp$Cutoff_Group_Level_Used
  scored$Cutoff_Apply_Mode <- tmp$Cutoff_Apply_Mode
  scored$Cutoff_Method <- tmp$Cutoff_Method
  scored$Engine_ID <- tmp$Engine_ID
  scored$Engine_Label <- tmp$Engine_Label

  cbind(cleaned, scored)
}

score_person <- function(data, engine_id = NULL, engine = NULL) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  if (nrow(data) != 1) {
    stop("score_person() expects exactly one row of input.", call. = FALSE)
  }

  score_dataset(data, engine_id = engine_id, engine = engine)
}
