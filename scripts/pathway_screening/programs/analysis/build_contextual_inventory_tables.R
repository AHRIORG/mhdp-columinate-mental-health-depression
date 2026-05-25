#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

script_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][[1]])
script_path <- normalizePath(script_arg, mustWork = TRUE)
source(file.path(dirname(dirname(script_path)), "utils", "pathway_screening_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- grep(glue("^--{name}="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(glue("^--{name}="), "", hit[[1]])
}

inventory_path <- arg_value(
  "inventory-path",
  pathway_screening_path("data", "inputs", "pathway_inventory.md")
)
output_dir <- arg_value(
  "output-dir",
  file.path(default_pathway_release_dir(), "tables")
)

if (!file.exists(inventory_path)) {
  stop("Inventory source not found: ", inventory_path)
}

lines <- readLines(inventory_path, warn = FALSE)
start_idx <- which(trimws(lines) == "## Pathway Inventory")
if (length(start_idx) == 0) {
  stop("Could not locate '## Pathway Inventory' in the supplied source file.")
}

inventory_lines <- lines[(start_idx[[1]] + 1):length(lines)]
inventory_lines <- inventory_lines[nzchar(trimws(inventory_lines))]

current_domain <- NA_character_
domain_order <- 0L
inventory_rows <- list()

for (line in inventory_lines) {
  line_trim <- trimws(line)

  if (str_detect(line_trim, "^###\\s+")) {
    domain_order <- domain_order + 1L
    current_domain <- str_remove(line_trim, "^###\\s+")
    next
  }

  if (!str_detect(line_trim, "^-\\s+")) {
    next
  }

  pathway_chain <- str_remove(line_trim, "^-\\s+")
  pathway_chain <- str_replace(pathway_chain, "\\.$", "")
  pathway_lower <- tolower(pathway_chain)

  pathway_valence <- case_when(
    str_detect(pathway_lower, "positive mental health$") ~ "Positive mental health",
    str_detect(pathway_lower, "poor mental health$") ~ "Poor mental health",
    TRUE ~ "Mixed or unspecified"
  )

  captured_domain <- case_when(
    str_detect(pathway_lower, "food|crop|income|subsistence|poverty|financial") ~ "Food insecurity and material deprivation",
    str_detect(pathway_lower, "violence|crime|unsafe|drug|alcohol|tavern|abandoned buildings|domestic violence") ~ "Violence exposure and unsafe environments",
    str_detect(pathway_lower, "alcohol|drug|substance") ~ "Alcohol use and substance-related risk",
    str_detect(pathway_lower, "sexual|truancy|peer|risky behaviours") ~ "Sexual debut and broader adolescent risk ecology",
    str_detect(pathway_lower, "grant") ~ "Social protection and government grant conditions",
    str_detect(pathway_lower, "school|learning center|punishment|fencing|monitoring") ~ "School climate and educational pathway conditions",
    str_detect(pathway_lower, "migration|family separation|loneliness") ~ "Family separation and migration-related stress",
    str_detect(pathway_lower, "water|house|housing|hygiene") ~ "Housing quality and water insecurity",
    str_detect(pathway_lower, "church|spiritual|prayer|witchcraft|ritual") ~ "Religion, spirituality, and cultural support or strain",
    str_detect(pathway_lower, "clinic|hospital|healthcare|treatment|care") ~ "Service access, care barriers, and stigma",
    str_detect(pathway_lower, "pollution|mine|weather|erosion|environment|dumping") ~ "Environmental degradation and climate stress",
    str_detect(pathway_lower, "social capital|community cohesion|excluded|support|sports|family decision-making") ~ "Social connectedness and community inclusion or exclusion",
    TRUE ~ "Cross-domain pathway"
  )

  operationalization_status <- case_when(
    captured_domain == "Food insecurity and material deprivation" ~ "Directly represented by food insecurity",
    captured_domain == "Violence exposure and unsafe environments" ~ "Directly represented by violence exposure",
    captured_domain == "Alcohol use and substance-related risk" ~ "Directly represented by alcohol use",
    captured_domain == "Sexual debut and broader adolescent risk ecology" ~ "Partially represented by ever having had sexual intercourse",
    captured_domain == "Social protection and government grant conditions" ~ "Partially represented by lack of government grant receipt",
    captured_domain == "School climate and educational pathway conditions" ~ "Partially represented by school disengagement",
    captured_domain == "Social connectedness and community inclusion or exclusion" ~ "Partially represented by low social support",
    TRUE ~ "Not directly harmonized in the present analysis"
  )

  inventory_rows[[length(inventory_rows) + 1L]] <- tibble(
    domain_order = domain_order,
    source_domain = current_domain,
    pathway_chain = pathway_chain,
    pathway_valence = pathway_valence,
    captured_domain = captured_domain,
    operationalization_status = operationalization_status
  )
}

inventory_long <- bind_rows(inventory_rows) |>
  mutate(row_id = row_number()) |>
  select(row_id, domain_order, source_domain, captured_domain, pathway_chain, pathway_valence, operationalization_status)

inventory_summary <- inventory_long |>
  count(captured_domain, operationalization_status, name = "n_pathways") |>
  arrange(desc(n_pathways), captured_domain)

manifest <- tibble(
  build_timestamp = format(Sys.time(), tz = "Africa/Johannesburg"),
  inventory_source = basename(inventory_path),
  long_rows = nrow(inventory_long),
  summary_rows = nrow(inventory_summary)
)

write_public_csv(inventory_long, file.path(output_dir, "contextual_pathway_inventory_long.csv"))
write_public_csv(inventory_summary, file.path(output_dir, "contextual_pathway_inventory_summary.csv"))
write_public_csv(manifest, file.path(output_dir, "contextual_pathway_inventory_manifest.csv"))

message("Contextual pathway inventory tables written to ", output_dir)
