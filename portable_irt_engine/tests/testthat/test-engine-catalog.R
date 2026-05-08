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

test_that("engine catalog exposes a default engine", {
  catalog <- list_engines()
  expect_true(nrow(catalog) >= 1)
  expect_true(any(catalog$default))
  expect_identical(default_engine_id(), catalog$engine_id[catalog$default][1])
})

test_that("default engine loads with required components", {
  engine <- load_engine()
  expect_true(all(c("parameters", "items", "cutoff_table", "phq_theta_map") %in% names(engine)))
  expect_true(length(engine$items) > 0)
})
