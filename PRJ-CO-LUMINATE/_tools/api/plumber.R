#* @apiTitle CO-LUMINATE IRT Calibration API
#* @apiDescription Local API scaffold for scoring SSQ responses with calibrated IRT engines.

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

detect_project_root <- function() {
  candidates <- c(
    getwd(),
    file.path(getwd(), "PRJ-CO-LUMINATE"),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "PRJ-CO-LUMINATE"),
    file.path(getwd(), "..", ".."),
    file.path(getwd(), "..", "..", "PRJ-CO-LUMINATE")
  )
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))

  for (root in candidates) {
    if (file.exists(file.path(root, "_quarto.yml")) &&
        dir.exists(file.path(root, "_tools"))) {
      return(root)
    }
  }

  stop("Unable to detect PRJ-CO-LUMINATE project root from current working directory.")
}

project_root <- detect_project_root()

resolve_engine_bundle <- function(root) {
  env_path <- Sys.getenv("COLUMINATE_ENGINES_PATH", unset = "")
  candidates <- c(
    env_path,
    file.path(root, "_tools", "_objects", "irt_joint_models.rds"),
    file.path(
      root, "..", "OBJ00-Datasets Preperation", "data_management",
      "data_examples", "1_staging_snippets", "derived_models", "irt_joint_models.rds"
    )
  )
  candidates <- unique(candidates[nzchar(candidates)])
  idx <- which(file.exists(candidates))[1]
  if (is.na(idx)) {
    stop(
      "Unable to locate irt_joint_models.rds. ",
      "Set COLUMINATE_ENGINES_PATH or place bundle in _tools/_objects."
    )
  }
  normalizePath(candidates[idx], winslash = "/", mustWork = TRUE)
}

resolve_scoring_script <- function(root) {
  candidates <- c(
    Sys.getenv("COLUMINATE_SCORING_SCRIPT", unset = ""),
    file.path(root, "_tools", "api", "scoring_runtime.R"),
    file.path(root, "..", "OBJ00-Datasets Preperation", "data_management", "scripts", "02_harmonization", "02b_irt_model.R")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  idx <- which(file.exists(candidates))[1]
  if (is.na(idx)) {
    stop("Unable to locate 02b_irt_model.R. Set COLUMINATE_SCORING_SCRIPT.")
  }
  normalizePath(candidates[idx], winslash = "/", mustWork = TRUE)
}

engine_bundle_path <- resolve_engine_bundle(project_root)
scoring_script_path <- resolve_scoring_script(project_root)
ENGINES <- readRDS(engine_bundle_path)

if (!is.list(ENGINES) || length(ENGINES) == 0) {
  stop("Engine bundle is empty or invalid: ", engine_bundle_path)
}

if (!exists("score_new_cohort", mode = "function")) {
  source(scoring_script_path)
}
if (!exists("score_new_cohort", mode = "function")) {
  stop("score_new_cohort() was not loaded from scoring script.")
}

default_engine_id <- names(ENGINES)[[1]]

#* CORS filter for browser-based clients
#* @filter cors
function(req, res) {
  origin <- Sys.getenv("COLUMINATE_API_CORS_ORIGIN", unset = "*")
  res$setHeader("Access-Control-Allow-Origin", origin)
  res$setHeader("Access-Control-Allow-Methods", "GET,POST,OPTIONS")
  res$setHeader("Access-Control-Allow-Headers", "Content-Type, Authorization")
  res$setHeader("Access-Control-Max-Age", "86400")

  if (identical(toupper(req$REQUEST_METHOD %||% ""), "OPTIONS")) {
    res$status <- 200L
    return(list())
  }
  plumber::forward()
}

extract_auth_header <- function(req) {
  headers <- req$HEADERS %||% req$headers %||% list()
  get_header <- function(x, key) {
    nms <- names(x)
    if (is.null(nms)) return(NULL)
    idx <- match(tolower(key), tolower(nms))
    if (is.na(idx)) return(NULL)
    x[[idx]]
  }

  header <- get_header(headers, "authorization") %||%
    req$HTTP_AUTHORIZATION %||%
    req$Authorization
  as.character(header %||% "")
}

extract_bearer_token <- function(header) {
  if (!nzchar(header)) return("")
  if (!grepl("^Bearer\\s+", header, ignore.case = TRUE)) return("")
  sub("^Bearer\\s+", "", header, ignore.case = TRUE)
}

is_public_endpoint <- function(path) {
  public_paths <- c("/health", "/engines", "/contract", "/openapi.json")
  if (path %in% public_paths) return(TRUE)
  startsWith(path, "/__docs__")
}

#* Optional Bearer auth filter for secured deployments
#* @filter auth
function(req, res) {
  required_key <- Sys.getenv("COLUMINATE_API_KEY", unset = "")
  if (!nzchar(required_key)) {
    return(plumber::forward())
  }

  method <- toupper(req$REQUEST_METHOD %||% "")
  path <- req$PATH_INFO %||% req$PATH_INFO_RAW %||% req$PATH %||% ""
  if (method == "OPTIONS" || is_public_endpoint(path)) {
    return(plumber::forward())
  }

  token <- extract_bearer_token(extract_auth_header(req))
  if (!nzchar(token) || !identical(token, required_key)) {
    res$status <- 401L
    return(list(error = "Unauthorized. Provide a valid Bearer token."))
  }

  plumber::forward()
}

coerce_item_value <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    return(NA_real_)
  }
  if (is.factor(x)) {
    return(as.numeric(x) - 1)
  }
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  x_chr <- trimws(tolower(as.character(x)))
  if (x_chr %in% c("yes", "true")) return(1)
  if (x_chr %in% c("no", "false")) return(0)
  out <- suppressWarnings(as.numeric(x_chr))
  if (is.na(out)) return(NA_real_)
  out
}

extract_body_json <- function(req) {
  if (!is.null(req$body) && length(req$body) > 0) {
    return(req$body)
  }
  txt <- req$postBody %||% ""
  if (!nzchar(txt)) return(NULL)
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

error_response <- function(res, status = 400L, message = "Bad request", details = NULL) {
  res$status <- status
  out <- list(error = message)
  if (!is.null(details)) out$details <- details
  out
}

engine_metadata <- function(engine_id, engine_obj) {
  ct <- engine_obj$cutoff_table %||% data.frame()
  overall <- if (nrow(ct) > 0) {
    ct[ct$Grouping == "none" & ct$Group_Level == "Overall", , drop = FALSE]
  } else {
    data.frame()
  }
  list(
    engine_id = engine_id,
    n_items = length(engine_obj$items %||% character()),
    items = engine_obj$items %||% character(),
    n_factors = as.integer(engine_obj$n_factors %||% NA_integer_),
    model_type = as.character(engine_obj$model_type %||% NA_character_),
    cutoff_apply_mode = as.character(engine_obj$cutoff_apply_mode %||% NA_character_),
    cutoff_groupings = if (nrow(ct) > 0) sort(unique(as.character(ct$Grouping))) else character(),
    overall_cutoff_theta = if (nrow(overall) > 0) as.numeric(overall$Cutoff_Theta[1]) else NA_real_
  )
}

engine_catalog <- lapply(names(ENGINES), function(id) engine_metadata(id, ENGINES[[id]]))
names(engine_catalog) <- names(ENGINES)

get_engine <- function(engine_id = NULL) {
  eid <- engine_id %||% default_engine_id
  if (!(eid %in% names(ENGINES))) {
    return(NULL)
  }
  ENGINES[[eid]]
}

parse_single_record <- function(payload, engine) {
  item_names <- engine$items
  top_level_items <- payload[item_names]
  top_level_items <- top_level_items[!vapply(top_level_items, is.null, logical(1))]

  responses <- payload$responses
  if (is.null(responses)) responses <- list()

  merged <- as.list(top_level_items)
  merged[names(responses)] <- responses

  missing_items <- setdiff(item_names, names(merged))
  if (length(missing_items) > 0) {
    stop("Missing required SSQ items: ", paste(missing_items, collapse = ", "))
  }

  row <- as.list(setNames(vector("list", length(item_names)), item_names))
  for (nm in item_names) {
    row[[nm]] <- coerce_item_value(merged[[nm]])
  }

  if (!is.null(payload$sex)) row$SEX <- as.character(payload$sex)
  if (!is.null(payload$age)) row$AGE <- suppressWarnings(as.numeric(payload$age))
  if (!is.null(payload$agegrp)) row$AGEGRP <- as.character(payload$agegrp)
  if (!is.null(payload$sex_agegrp)) row$SEX_AGEGRP <- as.character(payload$sex_agegrp)

  as.data.frame(row, stringsAsFactors = FALSE)
}

parse_batch_records <- function(payload, engine) {
  item_names <- engine$items

  records <- payload$records
  if (is.null(records)) {
    if (is.data.frame(payload)) {
      records <- split(payload, seq_len(nrow(payload)))
    } else {
      stop("Batch payload must include `records`.")
    }
  }

  if (is.data.frame(records)) {
    df <- records
  } else if (is.list(records)) {
    rows <- lapply(records, function(rec) {
      as.data.frame(rec, stringsAsFactors = FALSE)
    })
    df <- dplyr::bind_rows(rows)
  } else {
    stop("`records` must be an array of objects.")
  }

  missing_items <- setdiff(item_names, names(df))
  if (length(missing_items) > 0) {
    stop("Missing required SSQ items in batch: ", paste(missing_items, collapse = ", "))
  }

  for (nm in item_names) {
    df[[nm]] <- vapply(df[[nm]], coerce_item_value, numeric(1))
  }
  if ("AGE" %in% names(df)) {
    df$AGE <- suppressWarnings(as.numeric(df$AGE))
  }
  df
}

build_anchor_rows <- function(engine, reference_df = NULL) {
  item_names <- engine$items
  pars <- engine$parameters

  max_by_item <- vapply(item_names, function(itm) {
    if (is.null(pars) || !all(c("item", "name") %in% names(pars))) return(1)
    n_thresh <- sum(pars$item == itm & grepl("^d[0-9]+$", pars$name))
    max(1, n_thresh)
  }, numeric(1))

  low <- as.list(setNames(rep(0, length(item_names)), item_names))
  high <- as.list(setNames(as.numeric(max_by_item), item_names))

  if (!is.null(reference_df)) {
    if ("SEX" %in% names(reference_df)) {
      low$SEX <- as.character(reference_df$SEX[[1]])
      high$SEX <- as.character(reference_df$SEX[[1]])
    }
    if ("AGE" %in% names(reference_df)) {
      low$AGE <- suppressWarnings(as.numeric(reference_df$AGE[[1]]))
      high$AGE <- suppressWarnings(as.numeric(reference_df$AGE[[1]]))
    }
  }

  dplyr::bind_rows(
    as.data.frame(low, stringsAsFactors = FALSE),
    as.data.frame(high, stringsAsFactors = FALSE)
  )
}

score_with_engine <- function(input_df, engine) {
  if (nrow(input_df) == 0) {
    stop("No records to score.")
  }

  anchor <- build_anchor_rows(engine, input_df)
  aug <- dplyr::bind_rows(input_df, anchor)
  scored <- score_new_cohort(aug, engine)

  scored <- as.data.frame(scored, stringsAsFactors = FALSE)
  scored <- scored[seq_len(nrow(input_df)), , drop = FALSE]
  rownames(scored) <- NULL
  scored
}

api_contract <- list(
  version = "0.1.0",
  auth = list(
    type = "Bearer",
    header = "Authorization",
    required_when = "COLUMINATE_API_KEY is set"
  ),
  endpoints = list(
    list(
      path = "/health",
      method = "GET",
      description = "API health check."
    ),
    list(
      path = "/engines",
      method = "GET",
      description = "List available scoring engines and metadata."
    ),
    list(
      path = "/score/single",
      method = "POST",
      description = "Score one respondent.",
      body_example = list(
        engine_id = default_engine_id,
        responses = list(
          SSQ01 = 1, SSQ02 = 0, SSQ03 = 1, SSQ08 = 0, SSQ09 = 1,
          SSQ10 = 0, SSQ11 = 1, SSQ12 = 0, SSQ13 = 1, SSQ14 = 0
        ),
        sex = "Female",
        age = 19
      )
    ),
    list(
      path = "/score/batch",
      method = "POST",
      description = "Score multiple respondents.",
      body_example = list(
        engine_id = default_engine_id,
        records = list(
          list(
            SSQ01 = 1, SSQ02 = 0, SSQ03 = 1, SSQ08 = 0, SSQ09 = 1,
            SSQ10 = 0, SSQ11 = 1, SSQ12 = 0, SSQ13 = 1, SSQ14 = 0,
            SEX = "Female", AGE = 19
          ),
          list(
            SSQ01 = 0, SSQ02 = 0, SSQ03 = 0, SSQ08 = 0, SSQ09 = 0,
            SSQ10 = 0, SSQ11 = 0, SSQ12 = 0, SSQ13 = 0, SSQ14 = 0,
            SEX = "Male", AGE = 22
          )
        )
      )
    )
  )
)

#* Health check
#* @get /health
function() {
  list(
    status = "ok",
    api = "CO-LUMINATE IRT Calibration API",
    engines_loaded = length(ENGINES),
    default_engine_id = default_engine_id,
    engine_bundle_path = engine_bundle_path,
    scoring_script_path = scoring_script_path,
    timestamp = as.character(Sys.time())
  )
}

#* List available scoring engines
#* @get /engines
function() {
  list(
    default_engine_id = default_engine_id,
    engines = unname(engine_catalog)
  )
}

#* Get API contract
#* @get /contract
function() {
  api_contract
}

#* Score one respondent
#* @post /score/single
function(req, res) {
  payload <- tryCatch(extract_body_json(req), error = function(e) e)
  if (inherits(payload, "error") || is.null(payload)) {
    return(error_response(res, 400L, "Invalid or missing JSON body."))
  }

  engine_id <- payload$engine_id %||% payload$engine %||% default_engine_id
  engine <- get_engine(engine_id)
  if (is.null(engine)) {
    return(error_response(res, 400L, "Unknown engine_id.", list(engine_id = engine_id)))
  }

  row_df <- tryCatch(parse_single_record(payload, engine), error = function(e) e)
  if (inherits(row_df, "error")) {
    return(error_response(res, 400L, "Single-record payload validation failed.", row_df$message))
  }

  scored <- tryCatch(score_with_engine(row_df, engine), error = function(e) e)
  if (inherits(scored, "error")) {
    return(error_response(res, 422L, "Scoring failed.", scored$message))
  }

  list(
    engine_id = engine_id,
    output = as.list(scored[1, , drop = FALSE])
  )
}

#* Score batch respondents
#* @post /score/batch
function(req, res) {
  payload <- tryCatch(extract_body_json(req), error = function(e) e)
  if (inherits(payload, "error") || is.null(payload)) {
    return(error_response(res, 400L, "Invalid or missing JSON body."))
  }

  engine_id <- payload$engine_id %||% payload$engine %||% default_engine_id
  engine <- get_engine(engine_id)
  if (is.null(engine)) {
    return(error_response(res, 400L, "Unknown engine_id.", list(engine_id = engine_id)))
  }

  records_df <- tryCatch(parse_batch_records(payload, engine), error = function(e) e)
  if (inherits(records_df, "error")) {
    return(error_response(res, 400L, "Batch payload validation failed.", records_df$message))
  }

  max_batch_n <- suppressWarnings(as.integer(Sys.getenv("COLUMINATE_API_MAX_BATCH", "5000")))
  if (!is.na(max_batch_n) && nrow(records_df) > max_batch_n) {
    return(error_response(
      res, 413L, "Batch too large.",
      list(max_batch_n = max_batch_n, received_n = nrow(records_df))
    ))
  }

  scored <- tryCatch(score_with_engine(records_df, engine), error = function(e) e)
  if (inherits(scored, "error")) {
    return(error_response(res, 422L, "Scoring failed.", scored$message))
  }

  list(
    engine_id = engine_id,
    n_records = nrow(scored),
    output = scored
  )
}
