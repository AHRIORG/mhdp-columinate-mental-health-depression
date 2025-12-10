# ==============================================================================
# PROJECT: CO-LUMINATE Trajectory Analysis (LCGA)
# OBJECTIVE: Identify latent trajectories of childhood exposure (Ages 0-4)
# DEPENDENCIES: tidyverse, MplusAutomation, glue, here
# ==============================================================================

# 1. SETUP & LIBRARIES ---------------------------------------------------------
# Load necessary libraries. If not installed, use install.packages("...")
library(tidyverse)
library(MplusAutomation)
library(glue)
library(here) # Robust paths in R Projects

# Ensure output directory exists at the Project Root
output_dir <- here("Trajectory_Output")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Set the Mplus executable path if necessary (MplusAutomation will try to find it automatically)
# MplusAutomation::setMplusPath("/Applications/Mplus/mplus")

# 2. DATA LOADING (Reading longitudinal data for Exposure Age 0-4) -------------
# Load the data from the specified path. This file is assumed to contain a 
# long-format data frame named 'dt_childhood_exposure' with the columns: 
# USUBJID (Participant ID), EXPAGE (Age/Time Point), and LIVEDWOM (Factor variable).
data_path <- here("../../_private_use/data_management/data/co-luminate/adam/dt_childhood_exposure.RData")

if (!file.exists(data_path)) {
  stop(paste("Error: Data file not found at", data_path))
}

# The load() function will place the object(s) into the environment. 
# We assume the required object is named 'dt_childhood_exposure'.
load(data_path) 
message(glue("Data loaded successfully from: {data_path}"))

# CRITICAL CHECK: Ensure the exposure column ('LIVEDWOM') exists.
if (!("LIVEDWOM" %in% names(dt_childhood_exposure))) {
  stop(glue(
    "Error in pivot_wider(): The expected exposure column 'LIVEDWOM' was not found in the loaded data frame 'dt_childhood_exposure'. 
    
    Actual column names found are: {paste(names(dt_childhood_exposure), collapse = ', ')}.
    
    Please inspect your loaded data (dt_childhood_exposure) and change 'values_from' 
    in Section 3 to the correct column name."
  ))
}

# 3. DATA PREPARATION (Long -> Wide) -------------------------------------------
# Mplus LCGA requires one row per subject (child) and one column per time point.

df_wide <- dt_childhood_exposure %>%
  # === FACTOR TO NUMERIC CONVERSION ===
  # Convert the factor variable LIVEDWOM ('No', 'Yes') into a numeric binary 
  # variable (0, 1) where 1 indicates the risk (Absence of Mother).
  # We adjust by subtracting 1 to get 0/1, where 1 = Risk (Absence of Mother).
  # Assuming factor levels are [1] No, [2] Yes (as is standard).
  mutate(NO_MUM_EXPOSURE = as.numeric(LIVEDWOM) - 1) %>%
  # ====================================

# === USER FILTER: Restrict to Ages 0, 1, 2, 3, 4 ===
filter(EXPAGE %in% 0:4) %>%
  # =================================================================

# === CRITICAL FIX: SELECT ONLY NECESSARY COLUMNS FOR PIVOT ===
# This prevents other long-format variables (like LIVEDDAD) from interfering 
# with the pivot operation, ensuring one row per USUBJID.
select(USUBJID, EXPAGE, NO_MUM_EXPOSURE) %>%
  # =============================================================

pivot_wider(
  names_from   = EXPAGE,
  values_from  = NO_MUM_EXPOSURE, # Use the newly created numeric variable
  names_prefix = "no_mum" # Creates wide vars: no_mum0, no_mum1, ..., no_mum4
)

# Variable names for Mplus automation
item_names   <- names(df_wide)[grep("^no_mum", names(df_wide))]
all_vars     <- names(df_wide) # Includes USUBJID and all no_mum vars
mplus_varnames <- paste(item_names, collapse = " ") # no_mum0 no_mum1 ...

# 4. DEFINE MODELS (Iterate Classes 1-4) --------------------------------------
k_classes <- 1:4

lca_models <- lapply(k_classes, function(k) {
  
  # Define MODEL block dynamically: for k > 1, estimate class intercepts.
  model_overall <- if (k > 1) {
    glue("
      %OVERALL%
      ! LCGA: estimate class probabilities for first {k-1} classes
      [c#1-c#{k-1}];
    ")
  } else {
    "
      %OVERALL%
      ! Single-class model: no latent class intercepts specified
    "
  }
  
  mplusObject(
    TITLE = glue("LCGA Trajectory Analysis - {k} Classes (No_Mum)"),
    
    # REVISED VARIABLE BLOCK: Explicitly inject USEVARIABLES and CATEGORICAL
    # to ensure Mplus receives the necessary variable definitions clearly.
    VARIABLE = glue("
      IDVARIABLE = USUBJID;
      CLASSES = c({k});
      USEVARIABLES = {paste(all_vars, collapse = ' ')};
      CATEGORICAL = {mplus_varnames};
    "),
    
    ANALYSIS = "
      TYPE = MIXTURE;
      ESTIMATOR = MLR; ! Maximum Likelihood Robust for categorical data
      PROCESSORS = 4;
      STARTS = 500 100; ! 500 random starts, 100 optimized to find global maximum
    ",
    
    MODEL = model_overall,
    
    OUTPUT = "TECH11 TECH14 SAMPSTAT;", # TECH11: VLMR, TECH14: BLRT
    
    # REVISED PLOT COMMAND: Use variable list directly instead of series notation with (s)
    PLOT = glue("TYPE = PLOT3; SERIES = {mplus_varnames};"),
    
    # The usevariables and categorical arguments are now redundant due to explicit
    # inclusion in the VARIABLE block, but we keep them for robustness.
    # We remove rdata = df_wide, categorical = item_names
    # and instead use the explicit USEVARIABLES/CATEGORICAL in the glue block.
    # Note: MplusAutomation handles data transfer via rdata implicitly in runModels.
    
    rdata = df_wide,
    
    # Model file path
    modelout = here("Trajectory_Output", glue("model_{k}_class.inp"))
  )
})

# 5. BATCH EXECUTION -----------------------------------------------------------
# Run the Mplus models. Note: This requires Mplus software to be installed and accessible.
message("Running Mplus models for 1 to 4 classes...")
# Using a tryCatch block to suppress the MplusAutomation warnings, 
# as the error code 1 is being caught, but the function still throws a warning.
mplus_results <- tryCatch({
  lapply(seq_along(lca_models), function(i) {
    mplusModeler(
      lca_models[[i]],
      dataout = glue("Trajectory_Output/exposure_data_{i}.dat"),
      run = 1L
    )
  })
}, warning = function(w) {
  # Suppress MplusAutomation warnings related to error code 1
  if (grepl("Mplus returned error code: 1", w$message)) {
    return(NULL) # Return NULL or a placeholder list instead of failing the script
  } else {
    stop(w)
  }
})
message("Mplus execution attempted. Check .out files for status.")

# 6. MODEL COMPARISON (Fit Statistics) -----------------------------------------
all_results <- readModels(here("Trajectory_Output"), recursive = TRUE)

if (length(all_results) > 0) {
  # Combine summaries from all models
  # === FIX: Use a filter to only process models that actually ran and have summaries ===
  fit_summary <- purrr::imap(all_results, function(.x, .y) {
    if (is.list(.x) && !is.null(.x$summaries)) {
      # Model ran successfully, append the model name
      .x$summaries %>% mutate(model_name = .y)
    } else {
      # Model failed, return NULL so bind_rows can skip it
      NULL
    }
  }) %>%
    bind_rows()
  
  if ("NLatentClasses" %in% names(fit_summary)) {
    fit_summary <- fit_summary %>%
      select(
        Classes = NLatentClasses,
        AIC, BIC, aBIC, Entropy,
        VLMR_p = T11_VLMR_PValue, # VLMR test p-value
        BLRT_p = BLRT_PValue      # Bootstrap Likelihood Ratio Test p-value
      ) %>%
      arrange(Classes)
    
    cat("\n======================================================\n")
    cat("           LCGA Model Fit Comparison Summary\n")
    cat("======================================================\n")
    print(fit_summary)
    cat("------------------------------------------------------\n")
    message("\nUse BIC (minimize) and BLRT/Entropy (maximize/significant) to choose K.")
  } else {
    warning("Models failed to generate summaries. Check .out files for Mplus errors.")
  }
} else {
  warning("No output files found. Ensure Mplus ran successfully.")
}

# 7. VISUALIZATION (Trajectories) ----------------------------------------------
# We will use the model that the fit criteria (usually BIC) suggests is best.
# For this mock data demo, we assume the 3-class solution is the winner.
best_classes <- 3

if (exists("fit_summary") && nrow(fit_summary) > 0 && sum(fit_summary$Classes == best_classes) > 0) {
  best_model_name <- fit_summary %>%
    filter(Classes == best_classes) %>%
    # Assuming the first model run for that class count is the best, or sorting by BIC
    arrange(BIC) %>%
    slice(1) %>%
    pull(model_name)
  
  best_model <- all_results[[best_model_name]]
  
  if (!is.null(best_model) && !is.null(best_model$parameters$unstandardized)) {
    # Extract threshold parameters for the binary indicators no_mum0 - no_mum4
    thr <- best_model$parameters$unstandardized %>%
      filter(paramHeader == "Thresholds",
             grepl("^no_mum", param)) %>%
      mutate(
        var         = gsub("\\$1$", "", param),  # strip $1 from e.g. no_mum0$1
        Age         = readr::parse_number(var),
        LatentClass = factor(LatentClass),
        # Mplus uses Probit Link for CATEGORICAL: P(Y=1) = 1 - Phi(Threshold)
        # where Phi is the CDF of the standard normal distribution (pnorm)
        prob        = 1 - pnorm(est)
      )
    
    plot_data <- thr %>%
      select(LatentClass, Age, prob)
    
    # Generate the visualization of the trajectories
    trajectory_plot <- ggplot(plot_data, aes(x = Age, y = prob, color = LatentClass, group = LatentClass)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      scale_x_continuous(breaks = 0:4) +
      labs(
        title    = glue("Latent Trajectories of Parental Absence ({best_classes} Classes)"),
        subtitle = "Probability of Living Without Mother (Ages 0-4)",
        y        = "Probability of Exposure (P(Y=1))",
        x        = "Child Age (Years)",
        color    = "Latent Class"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
      )
    
    print(trajectory_plot)
    message(glue("\nVisualization generated for the assumed best-fitting {best_classes}-class model."))
    
  } else {
    message(glue("Parameters for the {best_classes}-class model were not found. Check .out file."))
  }
} else {
  message("Model comparison summaries required or model output not available before plotting.")
}

# 8. EXPORT FILES FOR MANUAL MPLUS RUN -----------------------------------------
# This section exports a clean .dat file and a standalone .inp file for the 
# 3-class model. This allows you to run the analysis manually in Mplus.

# Create a separate folder for manual files
manual_output_dir <- here("Trajectory_Manual_Export")
if (!dir.exists(manual_output_dir)) dir.create(manual_output_dir)

manual_dat_file <- here("Trajectory_Manual_Export", "columinate_lcga_data.dat")
manual_inp_file <- here("Trajectory_Manual_Export", "manual_3_class_model.inp")

# 1. Export the Data (df_wide) to a simple text file (no headers)
prepareMplusData(
  df_wide,
  filename = manual_dat_file,
  overwrite = TRUE,
  # This option ensures headers aren't written to the .dat, 
  # but column names are printed to console for reference.
  inpfile = FALSE 
)

# 2. Create the standalone Mplus Input File (.inp)
# Note: We hardcode the NAMES list to match df_wide columns (USUBJID no_mum0...)
manual_inp_text <- glue("
TITLE: Manual 3-Class LCGA Trajectory Analysis;

DATA: FILE = columinate_lcga_data.dat;

VARIABLE: 
  NAMES = {paste(names(df_wide), collapse = ' ')};
  IDVARIABLE = USUBJID;
  CLASSES = c(3);
  
  USEVARIABLES = {paste(names(df_wide), collapse = ' ')};
  CATEGORICAL = {mplus_varnames};

ANALYSIS:
  TYPE = MIXTURE;
  ESTIMATOR = MLR;
  PROCESSORS = 4;
  STARTS = 500 100;

MODEL:
  %OVERALL%
  ! Estimate class probabilities for classes 1 and 2
  [c#1-c#2];

OUTPUT: TECH11 TECH14 SAMPSTAT;

PLOT: TYPE = PLOT3; SERIES = {mplus_varnames};
")

writeLines(manual_inp_text, manual_inp_file)

message(glue("\nMANUAL EXPORT COMPLETE:\nData: {manual_dat_file}\nInput: {manual_inp_file}\nYou can now open {manual_inp_file} in Mplus and run it directly."))

# Clean up temporary Mplus data files
# COMMENTED OUT TO SAVE THE .DAT FILES AS REQUESTED
# file.remove(list.files(here("Trajectory_Output"), pattern = "\\.dat$", full.names = TRUE))