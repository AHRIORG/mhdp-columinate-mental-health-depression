#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glue)
  library(readr)
})

script_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][[1]])
source(file.path(dirname(dirname(normalizePath(script_arg, mustWork = TRUE))), "utils", "pathway_screening_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- grep(glue("^--{name}="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(glue("^--{name}="), "", hit[[1]])
}

input_dir <- arg_value(
  "input-dir",
  pathway_screening_path("data", "inputs", "pathway_screening_raw")
)
output_dir <- arg_value(
  "output-dir",
  file.path(default_pathway_release_dir(), "tables")
)

if (!dir.exists(input_dir)) {
  stop("Input directory not found: ", input_dir)
}

csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  stop("No CSV files found in: ", input_dir)
}

for (input_file in csv_files) {
  data <- read_csv(input_file, col_types = cols(.default = col_character()), show_col_types = FALSE)
  data <- drop_private_path_columns(data)
  write_public_csv(data, file.path(output_dir, basename(input_file)))
}

message("Sanitised ", length(csv_files), " pathway-screening table(s) into ", output_dir)
