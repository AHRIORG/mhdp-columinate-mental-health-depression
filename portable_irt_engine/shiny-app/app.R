app_dir <- normalizePath(getwd(), mustWork = TRUE)
pkg_root <- normalizePath(file.path(app_dir, ".."), mustWork = TRUE)
options(portableIRTEngine.package_root = pkg_root)

r_files <- list.files(file.path(pkg_root, "R"), pattern = "[.][Rr]$", full.names = TRUE)
for (r_file in sort(r_files)) {
  source(r_file, local = globalenv())
}

library(shiny)

engine_catalog <- list_engines()
selectable_catalog <- subset(
  engine_catalog,
  is.na(available_in_app) | available_in_app %in% TRUE
)

app_tabs <- data.frame(
  value = c("overview", "manual", "batch", "results", "privacy", "documentation"),
  label = c("Overview", "Manual Scoring", "Batch Scoring", "Results", "Privacy", "Documentation"),
  stringsAsFactors = FALSE
)

glossary_abbreviations <- data.frame(
  Abbreviation = c("SSQ-10", "PHQ-9", "IRT", "PPV", "NPV", "ROC", "TCC", "AUC", "CSV"),
  Meaning = c(
    "Shona Symptom Questionnaire-10",
    "Patient Health Questionnaire-9",
    "Item Response Theory",
    "Positive Predictive Value",
    "Negative Predictive Value",
    "Receiver Operating Characteristic",
    "Test Characteristic Curve",
    "Area Under the Curve",
    "Comma-Separated Values"
  ),
  stringsAsFactors = FALSE
)

glossary_terms <- data.frame(
  Term = c(
    "Harmonisation",
    "Harmonized theta",
    "Sensitivity",
    "Specificity",
    "Accuracy",
    "Depression prevalence",
    "Cutoff threshold",
    "Threshold approach",
    "Threshold scope",
    "Youden index threshold",
    "Model-based TCC threshold",
    "Positive predictive value",
    "Negative predictive value",
    "Cohen's kappa"
  ),
  Definition = c(
    "The process of making scores from different questionnaires comparable on a shared scale.",
    "The underlying depression score estimated by the harmonised IRT model.",
    "The proportion of truly depressed participants correctly identified by the engine.",
    "The proportion of truly non-depressed participants correctly identified by the engine.",
    "The overall proportion of classifications that are correct.",
    "The proportion of scored participants classified as depressed in the current scored dataset.",
    "The theta value above which the engine classifies a participant as depressed.",
    "The rule used to choose the decision threshold, for example prioritising sensitivity or specificity.",
    "The grouping level at which thresholds are applied, such as overall, by sex, by age group, or by sex plus age group.",
    "A threshold chosen to balance sensitivity and specificity using the Youden index.",
    "A threshold chosen from the model-implied test characteristic curve rather than a directly optimised empirical rule.",
    "The proportion of participants classified as depressed who are truly depressed in the reference sample.",
    "The proportion of participants classified as non-depressed who are truly non-depressed in the reference sample.",
    "A measure of agreement between the engine classification and the reference classification after accounting for agreement expected by chance."
  ),
  stringsAsFactors = FALSE
)

engine_label <- function(engine_id) {
  match_row <- engine_catalog[engine_catalog$engine_id == engine_id, , drop = FALSE]
  if (nrow(match_row) == 0) {
    return(engine_id)
  }

  model_label <- if (!is.na(match_row$n_factors[[1]])) paste0(match_row$n_factors[[1]], "-factor") else NULL
  scope_label <- pretty_cutoff_mode(match_row$cutoff_apply_mode[[1]])
  method_label <- pretty_cutoff_method(match_row$cutoff_method[[1]])

  pieces <- c("SSQ-10", model_label, scope_label, method_label)
  pieces <- pieces[!is.na(pieces) & nzchar(pieces)]
  paste(pieces, collapse = ", ")
}

format_metric <- function(value, digits = 2, fallback = "Not available") {
  if (length(value) == 0) {
    return(fallback)
  }

  out <- rep(fallback, length(value))
  ok <- !is.na(value)
  if (any(ok)) {
    out[ok] <- format(round(as.numeric(value[ok]), digits), nsmall = digits, trim = TRUE)
  }

  if (length(out) == 1) out[[1]] else out
}

format_percent <- function(value, digits = 1, fallback = "Not available") {
  if (length(value) == 0) {
    return(fallback)
  }

  out <- rep(fallback, length(value))
  ok <- !is.na(value)
  if (any(ok)) {
    out[ok] <- sprintf(paste0("%.", digits, "f%%"), 100 * as.numeric(value[ok]))
  }

  if (length(out) == 1) out[[1]] else out
}

engine_change_note <- function(engine_id) {
  meta <- engine_metadata(engine_id)
  if (nrow(meta) == 0) {
    return(NULL)
  }

  if (!is.na(meta$n_factors[[1]]) && meta$n_factors[[1]] == 1) {
    return(
      paste(
        "This packaged engine uses the shared 1-factor SSQ-10 scorer.",
        "Across the current production engines, switching the engine usually changes",
        "the stored cutoff strategy or subgroup threshold rather than the latent theta itself."
      )
    )
  }

  paste(
    "This engine uses a distinct factor structure, so both the scoring model and the cutoff",
    "strategy may differ from the default engine."
  )
}

pretty_age_group <- function(value) {
  if (length(value) == 0) {
    return(character())
  }

  out <- as.character(value)
  out[is.na(out) | !nzchar(out)] <- NA_character_
  out <- gsub("^Age_", "", out)
  out <- gsub("_", "-", out)
  out
}

pretty_cutoff_method <- function(value) {
  mapping <- c(
    sens_at_spec = "Sensitivity at target specificity",
    spec_at_sens = "Specificity at target sensitivity",
    youden = "Youden index threshold",
    model_tcc = "Model-based TCC threshold"
  )

  out <- as.character(value)
  idx <- match(out, names(mapping))
  out[!is.na(idx)] <- unname(mapping[idx[!is.na(idx)]])
  out
}

engine_metrics_lookup <- function(engine_id) {
  meta <- tryCatch(engine_provenance(engine_id), error = function(e) NULL)
  if (is.null(meta) || !"metrics_applied" %in% names(meta) || !is.data.frame(meta$metrics_applied)) {
    return(setNames(numeric(), character()))
  }

  metrics_df <- meta$metrics_applied
  vals <- suppressWarnings(as.numeric(metrics_df$Value))
  names(vals) <- metrics_df$Metric
  vals
}

format_metric_percent <- function(value, digits = 1, fallback = "Not available") {
  if (length(value) == 0 || all(is.na(value))) {
    return(fallback)
  }
  format_percent(value, digits = digits, fallback = fallback)
}

engine_performance_note <- function(engine_id) {
  metrics <- engine_metrics_lookup(engine_id)
  if (length(metrics) == 0) {
    return("Packaged reference performance metrics are not available for this engine.")
  }

  paste(
    "Reference performance in the packaging sample:",
    "sensitivity", format_metric_percent(metrics[["Sensitivity"]]),
    ", specificity", format_metric_percent(metrics[["Specificity"]]),
    ", accuracy", format_metric_percent(metrics[["Accuracy"]]),
    "."
  )
}

pretty_cutoff_mode <- function(value) {
  mapping <- c(
    none = "Overall",
    SEX = "Sex",
    AGEGRP = "Age group",
    SEX_AGEGRP = "Sex + age group",
    hierarchical = "Hierarchical"
  )

  out <- as.character(value)
  idx <- match(out, names(mapping))
  out[!is.na(idx)] <- unname(mapping[idx[!is.na(idx)]])
  out
}

engine_choices <- setNames(
  selectable_catalog$engine_id,
  vapply(selectable_catalog$engine_id, engine_label, character(1))
)

cutoff_target_label <- function(grouping, level) {
  if (length(grouping) == 0 || is.na(grouping) || !nzchar(grouping)) {
    return("Not available")
  }

  if (identical(grouping, "none")) {
    return("Overall")
  }

  if (length(level) == 0 || is.na(level) || !nzchar(level)) {
    return(pretty_cutoff_mode(grouping))
  }

  if (identical(grouping, "AGEGRP")) {
    return(pretty_age_group(level))
  }

  if (identical(grouping, "SEX")) {
    return(as.character(level))
  }

  if (identical(grouping, "SEX_AGEGRP")) {
    parts <- strsplit(as.character(level), "__", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      return(paste(parts[[1]], pretty_age_group(parts[[2]])))
    }
  }

  level_text <- gsub("__", " ", as.character(level), fixed = TRUE)
  level_text <- gsub("Age_", "", level_text, fixed = TRUE)
  level_text <- gsub("_", "-", level_text, fixed = TRUE)
  level_text
}

preview_id_columns <- function(df) {
  candidates <- names(df)[grepl("(^id$|_id$|participant|record|study|case)", names(df), ignore.case = TRUE)]
  candidates <- candidates[!grepl("^engine(_id)?$", candidates, ignore.case = TRUE)]
  unique(candidates)[seq_len(min(length(unique(candidates)), 2))]
}

compact_result_preview <- function(df, max_rows = 10) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    return(df)
  }

  id_cols <- preview_id_columns(df)
  core_cols <- c(
    "SEX",
    "AGE",
    "AGEGRP",
    "Items_Answered",
    "SSQ10_Total_Score",
    "Theta_Harmonized",
    "Depression_Binary",
    "PHQ_Expected_Sum",
    "PHQ_Prob_GE10",
    "Cutoff_Theta_Applied"
  )

  keep_cols <- unique(c(id_cols, core_cols))
  keep_cols <- keep_cols[keep_cols %in% names(df)]

  out <- utils::head(df[, keep_cols, drop = FALSE], max_rows)

  if (all(c("Cutoff_Grouping_Used", "Cutoff_Group_Level_Used") %in% names(df))) {
    out$Cutoff_Target <- vapply(
      seq_len(nrow(out)),
      function(i) {
        cutoff_target_label(df$Cutoff_Grouping_Used[i], df$Cutoff_Group_Level_Used[i])
      },
      character(1)
    )
  }

  if ("Depression_Binary" %in% names(out)) {
    out$Depression_Binary <- ifelse(
      is.na(out$Depression_Binary),
      NA_character_,
      ifelse(out$Depression_Binary == 1, "Yes", "No")
    )
  }

  if ("AGEGRP" %in% names(out)) {
    out$AGEGRP <- pretty_age_group(out$AGEGRP)
  }

  if ("PHQ_Prob_GE10" %in% names(out)) {
    out$PHQ_Prob_GE10 <- ifelse(
      is.na(out$PHQ_Prob_GE10),
      NA_character_,
      format_percent(out$PHQ_Prob_GE10)
    )
  }

  for (col in c("SSQ10_Total_Score", "Theta_Harmonized", "PHQ_Expected_Sum", "Cutoff_Theta_Applied")) {
    if (col %in% names(out)) {
      out[[col]] <- round(as.numeric(out[[col]]), 3)
    }
  }

  friendly_names <- c(
    SEX = "Sex",
    AGE = "Age",
    AGEGRP = "Age Group",
    Items_Answered = "Items Answered",
    SSQ10_Total_Score = "Total SSQ-10 Score",
    Theta_Harmonized = "Harmonized Theta",
    Depression_Binary = "Depression",
    PHQ_Expected_Sum = "Expected PHQ-9 Sum",
    PHQ_Prob_GE10 = "Probability PHQ-9 >= 10",
    Cutoff_Theta_Applied = "Applied Cutoff Theta",
    Cutoff_Target = "Cutoff Target",
    Engine_Label = "Engine"
  )

  rename_idx <- match(names(out), names(friendly_names))
  names(out) <- ifelse(
    is.na(rename_idx),
    names(out),
    unname(friendly_names[rename_idx])
  )

  out
}

metric_card <- function(label, value, note = NULL) {
  div(
    class = "stat-card",
    span(label),
    strong(value),
    if (!is.null(note)) tags$small(note)
  )
}

collapsible_panel <- function(title, copy = NULL, body) {
  tags$details(
    class = "panel collapsible-panel",
    tags$summary(
      class = "collapsible-summary",
      div(
        class = "collapsible-summary-copy",
        h3(title),
        if (!is.null(copy)) p(class = "panel-copy", copy)
      ),
      span(class = "collapsible-indicator", "Open")
    ),
    div(class = "collapsible-body", body)
  )
}

make_tab_button <- function(value, label, active = FALSE) {
  tags$button(
    type = "button",
    class = paste("tab-btn", if (active) "active"),
    `aria-selected` = if (active) "true" else "false",
    onclick = sprintf(
      "Shiny.setInputValue('goto_tab', '%s', {priority: 'event'})",
      value
    ),
    label
  )
}

state_badge <- function(text, class_name = "info") {
  tags$strong(class = paste("state-badge", class_name), text)
}

example_batch_path <- function() {
  normalizePath(
    file.path(pkg_root, "inst", "extdata", "example_batch.csv"),
    mustWork = TRUE
  )
}

manual_input_df <- function(input, engine) {
  age_value <- suppressWarnings(as.numeric(input$manual_age))
  if (length(age_value) == 0 || is.na(age_value)) {
    age_value <- NA_real_
  }

  row <- list(
    SEX = if (nzchar(input$manual_sex)) input$manual_sex else NA_character_,
    AGE = age_value
  )

  for (item in engine$items) {
    value <- input[[paste0("manual_", item)]]
    row[[item]] <- if (is.null(value) || value == "") NA_integer_ else as.integer(value)
  }

  as.data.frame(row, stringsAsFactors = FALSE)
}

ui <- fluidPage(
  tags$head(
    tags$title("Portable SSQ-10 IRT Engine"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$link(rel = "stylesheet", href = "styles.css")
  ),
  div(
    class = "app-shell portable-app",
    tags$header(
      class = "topbar",
      div(
        p(class = "eyebrow", "CO-LUMINATE Harmonisation Toolkit"),
        h1("Portable SSQ-10 IRT Engine"),
        p(
          class = "lede",
          paste(
            "Score SSQ-10 responses against the packaged harmonised depression engine,",
            "review outputs in a browser-friendly workspace, and download results",
            "without retaining raw uploaded data."
          )
        )
      ),
      div(
        class = "topbar-actions",
        span(class = "status-pill", "In-memory scoring"),
        span(class = "status-pill status-pill--soft", "No raw data retention"),
        tags$button(
          type = "button",
          class = "primary-btn",
          onclick = "Shiny.setInputValue('goto_tab', 'batch', {priority: 'event'})",
          "Open Batch Scoring"
        )
      )
    ),
    tags$section(
      class = "summary-banner panel",
      div(
        span(class = "mini-label", "Selected engine"),
        strong(textOutput("banner_engine", inline = TRUE))
      ),
      div(
        span(class = "mini-label", "Current workspace"),
        uiOutput("banner_tab")
      ),
      div(
        span(class = "mini-label", "Latest run"),
        p(class = "sentence-text", textOutput("banner_latest_summary", inline = TRUE))
      ),
      div(
        span(class = "mini-label", "Privacy mode"),
        p(
          class = "sentence-text",
          paste(
            "Uploaded data are scored only for the active session and are not",
            "written back to GitHub or project storage."
          )
        )
      )
    ),
    uiOutput("tab_nav"),
    tabsetPanel(
      id = "app_tabs",
      type = "hidden",
      selected = "overview",
      tabPanel(
        title = "Overview",
        value = "overview",
        div(
          class = "panel panel-intro",
          h2("Overview"),
          p(
            "This workspace packages the current Objective 3 SSQ-10 harmonisation engine into a reusable browser tool. Start here for orientation, then move into manual scoring, batch scoring, and result review."
          )
        ),
        div(
          class = "workflow-strip",
          div(
            class = "workflow-step",
            strong("1. Review the engine"),
            "Confirm which packaged engine is active, what input fields are expected, and how the stored cutoff is applied."
          ),
          div(
            class = "workflow-step",
            strong("2. Score one participant"),
            "Use the manual scoring view to test a single SSQ-10 profile and inspect all harmonised outputs before wider use."
          ),
          div(
            class = "workflow-step",
            strong("3. Score a cohort"),
            "Upload a CSV, validate the schema, score the file in memory, and download the resulting harmonised dataset."
          ),
          div(
            class = "workflow-step",
            strong("4. Review the outputs"),
            "Inspect the latest session summary, the most recent result preview, and the interpretation of the returned fields."
          )
        ),
        div(
          class = "overview-grid",
          tags$section(
            class = "panel stats-panel",
            h3("Current engine snapshot"),
            div(class = "stats-grid", uiOutput("overview_stats")),
            div(
              class = "overview-note",
              span(class = "mini-label", "Current selection"),
              strong(textOutput("overview_current", inline = TRUE)),
              p(class = "sentence-text", textOutput("overview_sentence", inline = TRUE))
            )
          ),
          tags$section(
            class = "panel guide-card",
            h3("Manual Scoring"),
            p(
              "Use the Manual Scoring view to enter one participant at a time, confirm how the packaged engine behaves, and inspect the harmonised metrics immediately."
            ),
            tags$ul(
              class = "guide-list",
              tags$li("Select the engine version and optional subgroup fields."),
              tags$li("Enter each SSQ-10 item as a binary response."),
              tags$li("Review theta, PHQ-oriented outputs, and binary classification.")
            ),
            tags$button(
              type = "button",
              class = "secondary-btn",
              onclick = "Shiny.setInputValue('goto_tab', 'manual', {priority: 'event'})",
              "Go to Manual Scoring"
            )
          ),
          tags$section(
            class = "panel guide-card",
            h3("Batch Scoring"),
            p(
              "Use the Batch Scoring view for cohort files. The app checks the schema, scores the data in memory, and returns a downloadable CSV for the same active session."
            ),
            tags$ul(
              class = "guide-list",
              tags$li("Choose the engine and upload a CSV file."),
              tags$li("Review validation notes before downloading."),
              tags$li("Preview the first scored rows before sharing outputs.")
            ),
            tags$button(
              type = "button",
              class = "secondary-btn",
              onclick = "Shiny.setInputValue('goto_tab', 'batch', {priority: 'event'})",
              "Go to Batch Scoring"
            )
          ),
          tags$section(
            class = "panel guide-card",
            h3("Documentation"),
            p(
              "Use the documentation and privacy views to confirm the expected schema, returned fields, and the privacy-by-design boundaries of the deployment."
            ),
            tags$ul(
              class = "guide-list",
              tags$li("Inspect the packaged engine catalog."),
              tags$li("Review the expected item schema and outputs."),
              tags$li("Confirm the no-retention handling for raw uploads.")
            ),
            tags$button(
              type = "button",
              class = "secondary-btn",
              onclick = "Shiny.setInputValue('goto_tab', 'documentation', {priority: 'event'})",
              "Go to Documentation"
            )
          )
        ),
        div(
          class = "overview-stack",
          tags$section(
            class = "panel",
            h3("Why use this engine"),
            p(
              class = "panel-copy",
              "Used appropriately, the engine helps researchers produce more comparable depression measures from SSQ-10 data while keeping the scoring workflow reproducible and transparent."
            ),
            tags$ul(
              class = "guide-list",
              tags$li("It places SSQ-10 responses onto a harmonised depression scale rather than relying only on a raw total score."),
              tags$li("It returns both continuous and binary outputs, which supports different reporting needs across studies."),
              tags$li("It makes the applied threshold explicit, including subgroup-specific targets where those are available."),
              tags$li("It supports batch scoring with the same packaged logic, which improves consistency across participants and cohorts."),
              tags$li("It helps researchers reuse existing SSQ-10 data without rerunning the full psychometric workflow each time.")
            )
          ),
          tags$section(
            class = "panel",
            h3("Important limitations"),
            p(
              class = "panel-copy",
              "This tool is intended to support harmonised research scoring. It is not a clinical diagnosis tool, and its outputs should be interpreted alongside study context, local governance, and appropriate subject-matter oversight."
            ),
            tags$ul(
              class = "guide-list",
              tags$li("The app scores binary SSQ-10 responses only. It does not accept SSQ-14 or PHQ-9 responses as direct input."),
              tags$li("Current production-safe engines share the same 1-factor SSQ-10 scorer, so changing engines will often change the applied threshold rather than the harmonised theta itself."),
              tags$li("Packaged performance metrics such as sensitivity and specificity come from the packaging sample and may not generalise unchanged to other cohorts, settings, languages, age ranges, or treatment-exposed groups."),
              tags$li("If subgroup variables needed for a subgroup-specific threshold are missing or unavailable, the app falls back to the stored overall threshold and records the applied target in the output."),
              tags$li("The binary depression classification is a harmonised research indicator and should not be used on its own for diagnosis, triage, or treatment decisions.")
            )
          ),
          collapsible_panel(
            title = "Packaged engines",
            copy = "The production default is set in the engine catalog and surfaced throughout the interface.",
            body = tagList(
              div(
                class = "data-status",
                "Only production-safe engines are exposed in the scoring views. Audit-only exports remain packaged for provenance but are hidden from live scoring until repaired."
              ),
              div(class = "table-shell", tableOutput("engine_table"))
            )
          ),
          collapsible_panel(
            title = "Expected input schema",
            copy = "SEX and AGE remain optional, but the SSQ-10 items are required for scoring.",
            body = div(class = "table-shell", tableOutput("schema_table"))
          )
        )
      ),
      tabPanel(
        title = "Manual Scoring",
        value = "manual",
        div(
          class = "panel panel-intro",
          h2("Manual Scoring"),
          p(
            "Use this view to score one participant at a time. This is useful for testing item combinations, reviewing the behaviour of the stored cutoff, and demonstrating the harmonised outputs without uploading a file."
          ),
          p(
            class = "panel-copy",
            "These outputs support harmonised research interpretation and should not be treated as a standalone clinical diagnosis."
          )
        ),
        div(
          class = "tool-grid tool-grid--manual",
          tags$section(
            class = "panel",
            h3("Participant input"),
            p(
              class = "panel-copy",
              "Choose the engine, optionally add subgroup fields, then capture the binary SSQ-10 responses."
            ),
            div(
              class = "form-grid",
              selectInput(
                "manual_engine",
                "Engine version",
                choices = engine_choices,
                selected = default_engine_id(),
                selectize = FALSE
              ),
              selectInput(
                "manual_sex",
                "SEX (optional)",
                choices = stats::setNames(c("", "Female", "Male"), c("", "Female", "Male")),
                selectize = FALSE
              ),
              textInput("manual_age", "AGE (optional)", value = "")
            ),
            h3("SSQ-10 items"),
            p(
              class = "panel-copy",
              "Record each item as 0 for No or 1 for Yes. Blank values are treated as missing."
            ),
            uiOutput("manual_item_inputs"),
            div(
              class = "canvas-actions",
              actionButton("score_manual", "Score participant", class = "primary-btn"),
              tags$button(
                type = "button",
                class = "secondary-btn",
                onclick = "Shiny.setInputValue('goto_tab', 'results', {priority: 'event'})",
                "View Latest Results"
              )
            )
          ),
          tags$section(
            class = "panel builder-card builder-card--result",
            h3("Scoring result"),
            p(
              class = "panel-copy",
              "The result card stays session-local. The in-app table shows a compact reader-friendly view, while downloaded files retain the full scored dataset."
            ),
            uiOutput("manual_highlights"),
            uiOutput("manual_status"),
            div(class = "table-shell", tableOutput("manual_result"))
          )
        ),
        tags$section(
          class = "panel",
          h3("Current profile across engines"),
          p(
            class = "panel-copy",
            "This comparison shows how the current manual profile behaves across all production-safe engines. It is especially useful when theta stays constant and the change is happening in the applied cutoff."
          ),
          uiOutput("manual_compare_note"),
          div(class = "table-shell", tableOutput("manual_engine_compare"))
        )
      ),
      tabPanel(
        title = "Batch Scoring",
        value = "batch",
        div(
          class = "panel panel-intro",
          h2("Batch Scoring"),
          p(
            "Upload a CSV cohort file, validate the schema, score the dataset in memory, and download the harmonised output for the current session."
          ),
          p(
            class = "panel-copy",
            "The returned binary classification is a research-use indicator derived from the packaged engine and should not replace diagnostic assessment."
          )
        ),
        div(
          class = "tool-grid tool-grid--batch",
          tags$section(
            class = "panel",
            h3("Upload data"),
            p(
              class = "panel-copy",
              "Use the example file if you want to test the workflow before uploading a real dataset."
            ),
            selectInput(
              "batch_engine",
              "Engine version",
              choices = engine_choices,
              selected = default_engine_id(),
              selectize = FALSE
            ),
            fileInput("batch_file", "Upload CSV", accept = c(".csv")),
            div(
              class = "canvas-actions",
              actionButton("score_batch", "Score uploaded file", class = "primary-btn"),
              downloadButton("download_batch", "Download scored CSV", class = "secondary-btn"),
              downloadButton("download_example", "Download example CSV", class = "ghost-btn")
            )
          ),
          tags$section(
            class = "panel builder-card builder-card--result",
            h3("Validation and preview"),
            p(
              class = "panel-copy",
              "Validation messages are shown here before you download outputs. The preview shows the first ten scored rows with a compact set of reader-friendly scoring columns."
            ),
            uiOutput("batch_highlights"),
            uiOutput("batch_status"),
            div(class = "table-shell", tableOutput("batch_preview"))
          )
        )
      ),
      tabPanel(
        title = "Results",
        value = "results",
        div(
          class = "panel panel-intro",
          h2("Results"),
          p(
            "This session summary shows the most recent successful scoring run and a preview of the returned data structure."
          )
        ),
        div(
          class = "results-grid",
          tags$section(
            class = "panel",
            h3("Latest summary"),
            div(class = "data-status", textOutput("latest_summary_text")),
            p(
              class = "panel-copy",
              "Only the latest result from the current app session is shown here. Closing the session clears the in-memory state. The preview is intentionally compact so the scoring outputs remain readable in the browser."
            ),
            div(class = "table-shell", tableOutput("latest_result_preview"))
          )
        )
      ),
      tabPanel(
        title = "Privacy",
        value = "privacy",
        div(
          class = "panel panel-intro",
          h2("Privacy"),
          p(
            "The intended deployment pattern is privacy-first: raw uploads stay within the active session and only user-requested exports are returned."
          )
        ),
        div(
          class = "panel-grid panel-grid--two",
          tags$section(
            class = "panel",
            h3("What the app does"),
            tags$ul(
              class = "guide-list",
              tags$li("Loads the selected engine and validates the uploaded schema."),
              tags$li("Scores the data in memory during the active session."),
              tags$li("Returns harmonised outputs directly to the user."),
              tags$li("Generates downloadable files only when the user requests them.")
            )
          ),
          tags$section(
            class = "panel",
            h3("What the app does not do"),
            tags$ul(
              class = "guide-list",
              tags$li("It does not write raw uploaded datasets back into the project repository."),
              tags$li("It does not persist session inputs to GitHub or project storage."),
              tags$li("It should not log raw item-level responses in production telemetry."),
              tags$li("It is not designed to replace governance for regulated or restricted data environments.")
            )
          )
        )
      ),
      tabPanel(
        title = "Documentation",
        value = "documentation",
        div(
          class = "panel panel-intro",
          h2("Documentation"),
          p(
            "Use this view to inspect the expected input contract and the interpretation of the returned harmonised fields."
          )
        ),
        div(
          class = "panel-grid panel-grid--two",
          tags$section(
            class = "panel",
            h3("Scoring outputs"),
            tags$ul(
              class = "guide-list",
              tags$li("SSQ10_Total_Score is the observed raw sum of endorsed SSQ-10 items."),
              tags$li("Theta_Harmonized is the underlying harmonised depression score."),
              tags$li("Depression_Binary applies the stored engine-specific cutoff strategy."),
              tags$li("PHQ_Expected_Sum maps the harmonised theta back into a PHQ-oriented scale."),
              tags$li("PHQ_Prob_GE10 estimates the probability of being at or above the PHQ-9 >= 10 threshold."),
              tags$li("The returned output also includes engine and cutoff metadata so result changes remain auditable.")
            )
          ),
          tags$section(
            class = "panel",
            h3("Expected input schema"),
            div(class = "table-shell", tableOutput("schema_table_doc"))
          )
        ),
        tags$section(
          class = "panel",
          h3("Benefits and intended strengths"),
          tags$ul(
            class = "guide-list",
            tags$li("The app provides a reproducible way to convert SSQ-10 responses into a harmonised depression metric."),
            tags$li("It reports multiple outputs together, including raw SSQ-10 total score, harmonised theta, PHQ-oriented summaries, and the binary classification."),
            tags$li("The packaged threshold metadata makes the scoring decision auditable rather than opaque."),
            tags$li("The same engine logic can be applied to single records or whole cohorts, which supports consistent scoring workflows."),
            tags$li("The in-memory scoring design supports privacy-conscious use by avoiding routine storage of raw uploaded data.")
          )
        ),
        tags$section(
          class = "panel",
          h3("Limitations and restrictions"),
          tags$ul(
            class = "guide-list",
            tags$li("This app is designed for harmonised SSQ-10 scoring. It is not a general-purpose mental health screening platform."),
            tags$li("The direct input contract assumes binary SSQ-10 item responses. Other questionnaires should be converted upstream only if a validated mapping exists."),
            tags$li("The current production-safe engine family is based on the validated 1-factor SSQ-10 scorer. Multidimensional audit exports remain packaged for provenance but are not exposed for live scoring until repaired."),
            tags$li("Threshold performance summaries are reference-sample properties, not guarantees of identical performance in new datasets."),
            tags$li("Subgroup-specific thresholds depend on correctly supplied subgroup fields. Where those fields are missing, unmatched, or outside the stored threshold table, the app falls back to the overall threshold."),
            tags$li("Treatment exposure, setting, language, or cohort differences may change how items function. Additional validation may be needed before using the engine in a new population."),
            tags$li("The returned depression classification should not be interpreted as a diagnosis and should not be used as the sole basis for clinical or administrative decisions.")
          )
        ),
        tags$section(
          class = "panel",
          h3("Glossary: Abbreviations"),
          p(
            class = "panel-copy",
            "These abbreviations appear across the engine table, performance summaries, and scored outputs."
          ),
          div(class = "table-shell", tableOutput("glossary_abbrev_table"))
        ),
        tags$section(
          class = "panel",
          h3("Glossary: Statistical Terms"),
          p(
            class = "panel-copy",
            "These short definitions are intended to make the app easier to read for non-specialist users."
          ),
          div(class = "table-shell", tableOutput("glossary_terms_table"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  default_engine <- load_engine(default_engine_id())
  latest_result <- reactiveVal(NULL)
  latest_summary <- reactiveVal("No scoring run yet.")
  latest_engine_id <- reactiveVal(default_engine_id())
  current_tab <- reactiveVal("overview")
  manual_scoring_started <- reactiveVal(FALSE)
  batch_scoring_started <- reactiveVal(FALSE)

  observeEvent(input$goto_tab, {
    if (!is.null(input$goto_tab) && input$goto_tab %in% app_tabs$value) {
      current_tab(input$goto_tab)
      updateTabsetPanel(session, "app_tabs", selected = input$goto_tab)
    }
  }, ignoreInit = TRUE)

  output$tab_nav <- renderUI({
    div(
      class = "tab-nav",
      role = "tablist",
      `aria-label` = "Portable IRT sections",
      lapply(seq_len(nrow(app_tabs)), function(i) {
        make_tab_button(
          value = app_tabs$value[[i]],
          label = app_tabs$label[[i]],
          active = identical(current_tab(), app_tabs$value[[i]])
        )
      })
    )
  })

  output$banner_engine <- renderText({
    tab_value <- current_tab()
    if (identical(tab_value, "manual") && !is.null(input$manual_engine)) {
      return(engine_label(input$manual_engine))
    }
    if (identical(tab_value, "batch") && !is.null(input$batch_engine)) {
      return(engine_label(input$batch_engine))
    }
    engine_label(latest_engine_id())
  })

  output$banner_tab <- renderUI({
    tab_value <- current_tab()
    tab_label <- app_tabs$label[match(tab_value, app_tabs$value)]
    badge_class <- switch(
      tab_value,
      manual = "reviewed",
      batch = "approved",
      results = "info",
      privacy = "draft",
      documentation = "info",
      "info"
    )
    state_badge(tab_label, badge_class)
  })

  output$banner_latest_summary <- renderText({
    latest_summary()
  })

  output$overview_stats <- renderUI({
    default_row <- engine_catalog[engine_catalog$engine_id == default_engine_id(), , drop = FALSE]
    tagList(
      metric_card("Engines", nrow(engine_catalog)),
      metric_card("SSQ-10 items", length(default_engine$items)),
      metric_card("Default model", paste0(default_row$n_factors[[1]], "-factor")),
      metric_card("Returned fields", "Core + audit", "Scores plus cutoff and engine metadata")
    )
  })

  output$overview_current <- renderText({
    engine_label(latest_engine_id())
  })

  output$overview_sentence <- renderText({
    latest_summary()
  })

  output$engine_table <- renderTable({
    display_catalog <- engine_catalog
    metrics_list <- lapply(display_catalog$engine_id, engine_metrics_lookup)
    get_metric <- function(metric_name) {
      vapply(
        metrics_list,
        function(x) if (metric_name %in% names(x)) format_metric_percent(x[[metric_name]]) else "Not available",
        character(1)
      )
    }

    data.frame(
      Engine = vapply(display_catalog$engine_id, engine_label, character(1)),
      Availability = ifelse(
        is.na(display_catalog$available_in_app) | display_catalog$available_in_app %in% TRUE,
        "Available",
        "Audit only"
      ),
      Model = paste0(display_catalog$n_factors, "-factor"),
      `Threshold approach` = pretty_cutoff_method(display_catalog$cutoff_method),
      `Threshold scope` = pretty_cutoff_mode(display_catalog$cutoff_apply_mode),
      Accuracy = get_metric("Accuracy"),
      Sensitivity = get_metric("Sensitivity"),
      Specificity = get_metric("Specificity"),
      PPV = get_metric("PPV"),
      NPV = get_metric("NPV"),
      `Cohen's kappa` = vapply(
        metrics_list,
        function(x) if ("Cohen_Kappa" %in% names(x)) format_metric(x[["Cohen_Kappa"]], digits = 3) else "Not available",
        character(1)
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)

  output$schema_table <- renderTable({
    engine_input_schema(default_engine_id())
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$schema_table_doc <- renderTable({
    engine_input_schema(default_engine_id())
  }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)

  output$glossary_abbrev_table <- renderTable({
    glossary_abbreviations
  }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)

  output$glossary_terms_table <- renderTable({
    glossary_terms
  }, striped = TRUE, bordered = TRUE, spacing = "s", rownames = FALSE)

  output$manual_item_inputs <- renderUI({
    engine_id <- if (is.null(input$manual_engine)) default_engine_id() else input$manual_engine
    engine <- load_engine(engine_id)
    div(
      class = "item-grid",
      lapply(engine$items, function(item) {
        selectInput(
          inputId = paste0("manual_", item),
          label = item,
          choices = stats::setNames(c("", "0", "1"), c("", "0 = No", "1 = Yes")),
          selected = "",
          selectize = FALSE
        )
      })
    )
  })

  observeEvent(input$score_manual, {
    manual_scoring_started(TRUE)
  }, ignoreInit = TRUE)

  manual_scored <- reactive({
    req(manual_scoring_started())
    engine <- load_engine(input$manual_engine)
    person <- manual_input_df(input, engine)
    score_person(person, engine = engine)
  })

  observeEvent(manual_scored(), {
    latest_result(manual_scored())
    latest_summary(paste("Manual scoring complete using", engine_label(input$manual_engine), "."))
    latest_engine_id(input$manual_engine)
  }, ignoreInit = TRUE)

  output$manual_highlights <- renderUI({
    if (!isTRUE(manual_scoring_started())) {
      return(div(class = "empty-note", "Score a participant to populate the harmonised outputs."))
    }
    result <- manual_scored()[1, , drop = FALSE]
    div(
      class = "stats-grid score-highlight-grid",
      metric_card("SSQ-10 total", format_metric(result$SSQ10_Total_Score, digits = 0)),
      metric_card("Theta", format_metric(result$Theta_Harmonized, digits = 3)),
      metric_card("Depression", ifelse(result$Depression_Binary == 1, "Yes", "No")),
      metric_card("PHQ expected sum", format_metric(result$PHQ_Expected_Sum, digits = 2)),
      metric_card("Prob. PHQ >= 10", format_percent(result$PHQ_Prob_GE10)),
      metric_card("Cutoff theta", format_metric(result$Cutoff_Theta_Applied, digits = 3)),
      metric_card("Cutoff target", cutoff_target_label(result$Cutoff_Grouping_Used, result$Cutoff_Group_Level_Used)),
      metric_card("Threshold scope", pretty_cutoff_mode(result$Cutoff_Apply_Mode)),
      metric_card("Threshold approach", pretty_cutoff_method(result$Cutoff_Method))
    )
  })

  output$manual_status <- renderUI({
    if (!isTRUE(manual_scoring_started())) {
      return(NULL)
    }
    div(
      class = "data-status",
      paste(
        "Manual scoring is active in memory for this session.",
        "After the first run, changing the engine or the item inputs recalculates the result immediately.",
        engine_change_note(input$manual_engine),
        engine_performance_note(input$manual_engine)
      )
    )
  })

  output$manual_result <- renderTable({
    req(manual_scoring_started())
    compact_result_preview(manual_scored(), max_rows = 1)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$manual_compare_note <- renderUI({
    req(manual_scoring_started())
    compare <- manual_engine_compare_data()
    theta_values <- compare$Theta_Harmonized
    theta_same <- length(theta_values) > 1 && all(abs(theta_values - theta_values[1]) < 1e-8)
    dep_levels <- unique(stats::na.omit(compare$Depression_Binary))

    if (theta_same) {
      return(
        div(
          class = "data-status",
          paste(
            "Theta is currently identical across the selectable engines for this profile.",
            "That is expected for the packaged 1-factor family.",
            "What differs is the cutoff strategy, and in some profiles that can still change the binary classification."
          )
        )
      )
    }

    if (length(dep_levels) > 1) {
      return(
        div(
          class = "data-status",
          "The engines are producing different classifications for this profile, so the comparison table below is showing a meaningful engine-specific decision difference."
        )
      )
    }

    div(
      class = "data-status",
      "The current profile is stable across the selectable engines right now. The table below still shows which cutoff and grouping rule each engine applied."
    )
  })

  manual_engine_compare_data <- reactive({
    req(manual_scoring_started())
    base_engine <- load_engine(input$manual_engine)
    person <- manual_input_df(input, base_engine)

    out <- lapply(selectable_catalog$engine_id, function(engine_id) {
      engine <- load_engine(engine_id)
      result <- score_person(person, engine = engine)[1, , drop = FALSE]
      meta <- engine_metadata(engine_id)

      data.frame(
        Selected = ifelse(identical(engine_id, input$manual_engine), "Yes", ""),
        Engine = engine_label(engine_id),
        Model = paste0(meta$n_factors[[1]], "-factor"),
        `Threshold Approach` = pretty_cutoff_method(meta$cutoff_method[[1]]),
        `Threshold Scope` = pretty_cutoff_mode(result$Cutoff_Apply_Mode[[1]]),
        Theta_Harmonized = result$Theta_Harmonized[[1]],
        Cutoff_Theta_Applied = result$Cutoff_Theta_Applied[[1]],
        Cutoff_Target = cutoff_target_label(
          result$Cutoff_Grouping_Used[[1]],
          result$Cutoff_Group_Level_Used[[1]]
        ),
        Depression_Binary = ifelse(
          is.na(result$Depression_Binary[[1]]),
          NA_character_,
          ifelse(result$Depression_Binary[[1]] == 1, "Yes", "No")
        ),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, out)
  })

  output$manual_engine_compare <- renderTable({
    req(manual_scoring_started())
    compare <- manual_engine_compare_data()
    compare$Theta_Harmonized <- round(compare$Theta_Harmonized, 3)
    compare$Cutoff_Theta_Applied <- round(compare$Cutoff_Theta_Applied, 3)
    compare
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  observeEvent(input$score_batch, {
    batch_scoring_started(TRUE)
  }, ignoreInit = TRUE)

  batch_raw <- reactive({
    req(batch_scoring_started())
    req(input$batch_file)
    utils::read.csv(input$batch_file$datapath, stringsAsFactors = FALSE)
  })

  batch_validation <- reactive({
    req(batch_scoring_started())
    req(input$batch_file)
    validate_ssq10_data(batch_raw(), engine = load_engine(input$batch_engine))
  })

  batch_scored <- reactive({
    req(batch_scoring_started())
    validation <- batch_validation()
    if (!validation$valid) {
      return(NULL)
    }

    engine <- load_engine(input$batch_engine)
    score_dataset(batch_raw(), engine = engine)
  })

  observeEvent(batch_scored(), {
    result <- batch_scored()
    if (is.null(result)) {
      return()
    }

    latest_result(result)
    latest_summary(
      paste(
        "Batch scoring complete for", nrow(result), "row(s) using",
        engine_label(input$batch_engine), "."
      )
    )
    latest_engine_id(input$batch_engine)
  }, ignoreInit = TRUE)

  output$batch_highlights <- renderUI({
    if (!isTRUE(batch_scoring_started())) {
      return(div(class = "empty-note", "Upload and score a file to populate the batch summary."))
    }

    validation <- batch_validation()
    if (!validation$valid || is.null(batch_scored())) {
      return(NULL)
    }

    result <- batch_scored()
    prevalence <- mean(result$Depression_Binary, na.rm = TRUE)
    div(
      class = "stats-grid score-highlight-grid",
      metric_card("Rows scored", nrow(result)),
      metric_card("Mean theta", format_metric(mean(result$Theta_Harmonized, na.rm = TRUE), digits = 3)),
      metric_card("Depression prevalence", sprintf("%.1f%%", 100 * prevalence)),
      metric_card("Engine", paste0(load_engine(input$batch_engine)$n_factors, "-factor")),
      metric_card("Threshold scope", pretty_cutoff_mode(unique(result$Cutoff_Apply_Mode)[1])),
      metric_card("Threshold approach", pretty_cutoff_method(unique(result$Cutoff_Method)[1]))
    )
  })

  output$batch_status <- renderUI({
    if (!isTRUE(batch_scoring_started())) {
      return(NULL)
    }

    validation <- batch_validation()
    if (validation$valid) {
      return(
        div(
          class = "data-status",
          paste(
            "Validation passed for", nrow(batch_raw()),
            "row(s). Data were scored in memory only.",
            "After the first run, changing the engine recalculates the loaded file for the active session.",
            engine_performance_note(input$batch_engine)
          )
        )
      )
    }

    div(
      class = "data-status data-status--warning",
      tags$strong("Validation needs attention."),
      tags$ul(class = "guide-list", lapply(validation$issues, tags$li))
    )
  })

  output$batch_preview <- renderTable({
    req(batch_scoring_started())
    req(batch_validation()$valid)
    req(batch_scored())
    compact_result_preview(batch_scored(), max_rows = 10)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$download_batch <- downloadHandler(
    filename = function() {
      paste0("portable_irt_scores_", Sys.Date(), ".csv")
    },
    content = function(file) {
      scored <- batch_scored()
      req(!is.null(scored))
      utils::write.csv(scored, file, row.names = FALSE)
    }
  )

  output$download_example <- downloadHandler(
    filename = function() {
      "portable_irt_example_batch.csv"
    },
    content = function(file) {
      file.copy(example_batch_path(), file)
    }
  )

  output$latest_summary_text <- renderText({
    latest_summary()
  })

  output$latest_result_preview <- renderTable({
    req(latest_result())
    compact_result_preview(latest_result(), max_rows = 10)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
}

shinyApp(ui = ui, server = server)
