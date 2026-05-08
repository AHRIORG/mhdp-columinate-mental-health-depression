options(
  portableIRTEngine.package_root = normalizePath(
    file.path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
)

for (r_file in sort(list.files(file.path("..", "..", "R"), pattern = "[.][Rr]$", full.names = TRUE))) {
  source(r_file, local = environment())
}

test_that("validation accepts the packaged example batch", {
  example_path <- file.path("..", "..", "inst", "extdata", "example_batch.csv")
  example_df <- utils::read.csv(example_path, stringsAsFactors = FALSE)
  validation <- validate_ssq10_data(example_df)
  expect_true(validation$valid)
})

test_that("validation flags missing required item columns", {
  example_path <- file.path("..", "..", "inst", "extdata", "example_batch.csv")
  example_df <- utils::read.csv(example_path, stringsAsFactors = FALSE)
  example_df$SSQ14 <- NULL
  validation <- validate_ssq10_data(example_df)
  expect_false(validation$valid)
  expect_true(any(grepl("Missing required SSQ-10 item columns", validation$issues, fixed = TRUE)))
})
