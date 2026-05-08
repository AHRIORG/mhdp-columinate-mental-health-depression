app_dir <- normalizePath(file.path(getwd(), "shiny-app"), mustWork = TRUE)
host <- Sys.getenv("PORTABLE_IRT_HOST", unset = "0.0.0.0")
port <- suppressWarnings(as.integer(Sys.getenv("PORT", unset = Sys.getenv("PORTABLE_IRT_PORT", unset = "8080"))))

if (is.na(port)) {
  port <- 8080L
}

shiny::runApp(
  appDir = app_dir,
  host = host,
  port = port,
  launch.browser = FALSE
)
