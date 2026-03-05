#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(plumber)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 0) {
  stop("Unable to resolve script path from commandArgs().")
}

this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(this_file), "..", ".."), winslash = "/", mustWork = TRUE)
setwd(project_root)

api_file <- file.path("_tools", "api", "plumber.R")
if (!file.exists(api_file)) {
  stop("Missing API file: ", api_file)
}

port <- suppressWarnings(as.integer(Sys.getenv("COLUMINATE_API_PORT", unset = "8088")))
if (is.na(port)) {
  port <- 8088L
}
host <- Sys.getenv("COLUMINATE_API_HOST", unset = "127.0.0.1")

cat("Starting CO-LUMINATE IRT API\n")
cat("Project root:", project_root, "\n")
cat("URL: http://", host, ":", port, "\n", sep = "")
cat("Contract: http://", host, ":", port, "/contract\n", sep = "")

pr <- plumber::plumb(api_file)
pr$run(host = host, port = port)
