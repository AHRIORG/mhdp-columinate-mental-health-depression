suppressPackageStartupMessages({
  library(jsonlite)
})

root <- normalizePath(file.path(getwd(), "."), mustWork = TRUE)
options(portableIRTEngine.package_root = root)

source(file.path(root, "R", "helpers.R"))
source(file.path(root, "R", "engine_catalog.R"))
source(file.path(root, "R", "input_validation.R"))
source(file.path(root, "R", "cutoffs.R"))
source(file.path(root, "R", "scoring.R"))
source(file.path(root, "data-raw", "export_browser_engine_bundle.R"))

browser_dir <- file.path(root, "..", "website", "_tools", "irt-scoring-static", "engines")
bundle_path <- export_browser_engine_bundle(root = root, out_dir = browser_dir)
example_path <- file.path(browser_dir, "example_batch.csv")
node_script <- file.path(root, "..", "website", "_tools", "irt-scoring-static", "smoke_test.mjs")

catalog <- list_engines()
catalog <- catalog[
  catalog$available_in_app %in% TRUE &
    catalog$scoring_ready %in% TRUE &
    catalog$n_factors == 1L,
  ,
  drop = FALSE
]

if (nrow(catalog) == 0) {
  stop("No production-safe 1-factor engines found for smoke testing.", call. = FALSE)
}

example_batch <- utils::read.csv(example_path, stringsAsFactors = FALSE)

compare_frame <- function(engine_id) {
  engine <- load_engine(engine_id)
  r_scored <- score_dataset(example_batch, engine = engine)
  js_stderr <- tempfile(fileext = ".log")
  on.exit(unlink(js_stderr), add = TRUE)
  js_raw <- system2(
    "node",
    c(shQuote(node_script), shQuote(bundle_path), shQuote(example_path), shQuote(engine_id)),
    stdout = TRUE,
    stderr = js_stderr
  )
  js_status <- attr(js_raw, "status") %||% 0L
  if (!identical(js_status, 0L)) {
    err_text <- if (file.exists(js_stderr)) paste(readLines(js_stderr, warn = FALSE), collapse = "\n") else ""
    stop("Browser smoke runner failed for ", engine_id, if (nzchar(err_text)) paste0(": ", err_text), call. = FALSE)
  }
  js_parsed <- jsonlite::fromJSON(paste(js_raw, collapse = "\n"), simplifyDataFrame = TRUE)
  js_scored <- as.data.frame(js_parsed$scored, stringsAsFactors = FALSE)

  keep_cols <- c(
    "participant_id",
    "SEX",
    "AGE",
    "AGEGRP",
    "SSQ10_Total_Score",
    "Theta_Harmonized",
    "Depression_Binary",
    "PHQ_Expected_Sum",
    "PHQ_Prob_GE10",
    "Cutoff_Theta_Applied",
    "Cutoff_Grouping_Used",
    "Cutoff_Group_Level_Used",
    "Cutoff_Apply_Mode",
    "Cutoff_Method",
    "Engine_ID",
    "Items_Answered"
  )

  list(
    r = r_scored[, keep_cols, drop = FALSE],
    js = js_scored[, keep_cols, drop = FALSE]
  )
}

numeric_cols <- c(
  "SSQ10_Total_Score",
  "Theta_Harmonized",
  "PHQ_Expected_Sum",
  "PHQ_Prob_GE10",
  "Cutoff_Theta_Applied",
  "Items_Answered"
)
exact_cols <- c(
  "participant_id",
  "SEX",
  "AGE",
  "AGEGRP",
  "Depression_Binary",
  "Cutoff_Grouping_Used",
  "Cutoff_Group_Level_Used",
  "Cutoff_Apply_Mode",
  "Cutoff_Method",
  "Engine_ID"
)

tolerance <- 0.02
summary_rows <- list()

for (engine_id in catalog$engine_id) {
  compared <- compare_frame(engine_id)

  numeric_diff <- sapply(numeric_cols, function(col) {
    max(abs(as.numeric(compared$r[[col]]) - as.numeric(compared$js[[col]])), na.rm = TRUE)
  })
  numeric_diff[!is.finite(numeric_diff)] <- 0

  exact_match <- sapply(exact_cols, function(col) {
    identical(as.character(compared$r[[col]]), as.character(compared$js[[col]]))
  })

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    engine_id = engine_id,
    max_theta_abs_diff = numeric_diff[["Theta_Harmonized"]],
    max_phq_sum_abs_diff = numeric_diff[["PHQ_Expected_Sum"]],
    max_phq_prob_abs_diff = numeric_diff[["PHQ_Prob_GE10"]],
    max_cutoff_abs_diff = numeric_diff[["Cutoff_Theta_Applied"]],
    exact_fields_ok = all(exact_match),
    stringsAsFactors = FALSE
  )

  if (any(numeric_diff > tolerance) || !all(exact_match)) {
    print(compared$r)
    print(compared$js)
    stop(
      paste0(
        "Browser smoke test failed for ", engine_id,
        ". Max numeric diffs: ",
        paste(names(numeric_diff), signif(numeric_diff, 6), collapse = "; "),
        ". Exact fields ok: ", all(exact_match)
      ),
      call. = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
summary_path <- file.path(root, "..", "website", "_tools", "irt-scoring-static", "engines", "browser_smoke_summary.csv")
utils::write.csv(summary_df, summary_path, row.names = FALSE, na = "")

print(summary_df)
message("Browser smoke test passed for ", nrow(summary_df), " engine(s).")
message("Summary written to: ", summary_path)
