# ==============================================================================
# PROJECT: CO-LUMINATE Trajectory Analysis (LCGA) - Using Mplus Variable Names
# OBJECTIVE: Identify latent trajectories of "No Mother" exposure (Ages 0-4)
# DEPENDENCIES: tidyverse, MplusAutomation, glue, here
# ==============================================================================

# 1. SETUP & LIBRARIES ---------------------------------------------------------
library(tidyverse)
library(MplusAutomation)
library(glue)
library(here) 

# Ensure output directory exists at the Project Root
output_dir <- here("Trajectory_Output_Nomum")
if (!dir.exists(output_dir)) dir.create(output_dir)

# 2. DATA LOADING (Reading longitudinal data for Exposure Age 1-5) ----
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
# Use the Mplus variable names: USUBJID for ID, LWOM0-LWOM4 for DVs (5 time points).

df_wide <- dt_childhood_exposure %>%
  # === FACTOR TO NUMERIC CONVERSION ===
  # Convert the factor variable LIVEDWOM ('No', 'Yes') into a numeric binary 
  # variable (0, 1) where 1 indicates the risk (Absence of Mother).
  # We adjust by subtracting 1 to get 0/1, where 1 = Risk (Absence of Mother).
  # Assuming factor levels are [1] No, [2] Yes (as is standard).
  mutate(NO_MUM_EXPOSURE = as.numeric(LIVEDWOM) - 1) %>%
  # ====================================

# === USER FILTER: Restrict to Ages 0, 1, 2, 3, 4 ===
filter(EXPAGE %in% 1:5) %>%
  # =================================================================

# === CRITICAL FIX: SELECT ONLY NECESSARY COLUMNS FOR PIVOT ===
# This prevents other long-format variables (like LIVEDDAD) from interfering 
# with the pivot operation, ensuring one row per USUBJID.
select(USUBJID, EXPAGE, NO_MUM_EXPOSURE) %>%
  # =============================================================

pivot_wider(
  names_from   = EXPAGE,
  values_from  = NO_MUM_EXPOSURE, # Use the newly created numeric variable
  # NOTE: We use "LWOM" prefix (short for Lived Without Mother) instead of "LIVEDWOM"
  # because Mplus often truncates variable names > 8 chars, causing duplicates 
  # (e.g., LIVEDWOM0 and LIVEDWOM1 both become LIVEDWOM).
  names_prefix = "LWOM" # Creates wide vars: LWOM0, LWOM1, ..., LWOM4
) %>%
  
  # === REORDER COLUMNS EXPLICITLY ===
  # Ensure columns are in the correct chronological order (0, 1, 2, 3, 4)
  # This matches the expected input for the SERIES command.
  select(USUBJID, LWOM1, LWOM2, LWOM3, LWOM4, LWOM5)


# Variable names for Mplus automation
item_names   <- names(df_wide)[grep("^LWOM", names(df_wide))] # LWOM0 to LWOM4
all_vars     <- names(df_wide) # USUBJID and all LWOM vars
# The variables are now expected to be: LWOM0, LWOM1, LWOM2, LWOM3, LWOM4
mplus_varnames <- paste(item_names, collapse = " ") 

# Construct explicit series string for PLOT command: LWOM0(0) LWOM1(1)...
# This aligns the variable name with the specific time point.
mplus_series_string <- paste0(item_names, "(", 0:4, ")", collapse = " ")

# 4. DEFINE MODELS (Iterate Classes 1-4) --------------------------------------
k_classes <- 1:4

lca_models <- lapply(k_classes, function(k) {
  
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
    
    # VARIABLE BLOCK: Updated based on successful manual run
    # Removed MISSING = .; because MplusAutomation automatically adds it
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
      STARTS = 500 100; ! 500 random starts, 100 optimized
    ",
    
    MODEL = model_overall,
    
    OUTPUT = "TECH11 TECH14 SAMPSTAT;", # TECH11: VLMR, TECH14: BLRT
    
    # PLOT COMMAND: Updated to use explicit time scores (e.g., LWOM0(0))
    PLOT = glue("TYPE = PLOT3; SERIES = {mplus_series_string};"),
    
    # SAVEDATA COMMAND: Essential for Adjudication (Assigning individuals to classes)
    # SAVE = CPROB saves the posterior probabilities and the most likely class membership.
    SAVEDATA = glue("FILE = savedata_{k}.dat; SAVE = CPROB;"),
    
    rdata = df_wide,
    
    # Model file path
    # FIX: Use ONLY the filename here, not the full path.
    modelout = glue("model_{k}_class.inp") 
  )
})

# 5. BATCH EXECUTION -----------------------------------------------------------
message("Running Mplus models for 1 to 4 classes (LWOM trajectory, 5 time points: Ages 0-4)...")

# Store the current working directory to restore it later
current_wd <- getwd()

# CHANGE WORKING DIRECTORY: 
# We switch to the output directory so Mplus only sees filenames, not long paths.
setwd(output_dir)

# Using simple lapply with error handling to run models. 
mplus_results <- lapply(seq_along(lca_models), function(i) {
  tryCatch({
    mplusModeler(
      lca_models[[i]],
      # FIX: Use relative filename for data output
      dataout = glue("exposure_data_{i}.dat"),
      run = 1L
    )
  }, error = function(e) {
    message(glue("Error running Model {i}: {e$message}"))
    return(NULL)
  })
})

# ALWAYS restore the original working directory
setwd(current_wd)

message("Mplus execution attempted. Checking for output files.")

# 6. MODEL COMPARISON (Fit Statistics) -----------------------------------------
all_results <- readModels(here("Trajectory_Output_Nomum"), recursive = TRUE)

if (length(all_results) > 0) {
  # === FIX: Use a filter to only process models that actually ran and have summaries ===
  fit_summary <- purrr::imap(all_results, function(.x, .y) {
    # Check if .x is a list (valid object) AND has a 'summaries' element
    if (is.list(.x) && !is.null(.x$summaries)) {
      # Model ran successfully, append the model name
      .x$summaries %>% mutate(model_name = .y)
    } else {
      # Model failed or invalid object, return NULL so bind_rows can skip it
      NULL
    }
  }) %>%
    bind_rows()
  # ====================================================================================
  
  if ("NLatentClasses" %in% names(fit_summary) && nrow(fit_summary) > 0) {
    fit_summary <- fit_summary %>%
      select(
        model_name, # === FIX: Include model_name so it can be pulled later ===
        Classes = NLatentClasses,
        AIC, BIC, aBIC, Entropy,
        VLMR_p = T11_VLMR_PValue,
        BLRT_p = BLRT_PValue
      ) %>%
      arrange(Classes)
    
    cat("\n======================================================\n")
    cat("  LCGA Model Fit Comparison Summary (LWOM Trajectories)\n")
    cat("======================================================\n")
    print(fit_summary)
    cat("------------------------------------------------------\n")
    message("\nUse BIC (minimize) and BLRT/Entropy (maximize/significant) to choose K.")
  } else {
    warning("No successful model summaries found. Check .out files for Mplus errors.")
  }
} else {
  warning("No output files found. Ensure Mplus ran successfully.")
}

# 7. VISUALIZATION (Trajectories) ----------------------------------------------
best_classes <- 4

if (exists("fit_summary") && nrow(fit_summary) > 0 && sum(fit_summary$Classes == best_classes) > 0) {
  best_model_name <- fit_summary %>%
    filter(Classes == best_classes) %>%
    arrange(BIC) %>%
    slice(1) %>%
    pull(model_name)
  
  best_model <- all_results[[best_model_name]]
  
  if (!is.null(best_model) && !is.null(best_model$parameters$unstandardized)) {
    
    # === DEBUGGING: INSPECT PARAMETER NAMES ===
    message("\n--- DEBUG: INSPECTING PARAMETER NAMES ---")
    all_params <- best_model$parameters$unstandardized$param
    print(head(unique(all_params), 20))
    message("-----------------------------------------\n")
    # ==========================================
    
    # Extract threshold parameters for the binary indicators LWOM0 - LWOM4
    thr <- best_model$parameters$unstandardized %>%
      filter(paramHeader == "Thresholds",
             # FIX: Added ignore.case = TRUE because Mplus output variables are UPPERCASE
             # Regex now matches the shorter prefix 'LWOM'
             grepl("LWOM", param, ignore.case = TRUE)) %>% 
      mutate(
        # Remove '$1' OR '.Cat.1' suffix if present to get variable name
        var         = gsub("(\\$1|\\.Cat\\.1)$", "", param), 
        Age         = readr::parse_number(var),
        LatentClass = factor(LatentClass),
        # P(Y=1) = 1 - Phi(Threshold)
        prob        = 1 - pnorm(est)
      )
    
    # Check if we successfully extracted data
    if (nrow(thr) == 0) {
      warning("Plot data is empty! Check the debug output above for parameter name mismatches.")
    } else {
      print(head(thr)) # Show preview of plot data
    }
    
    plot_data <- thr %>%
      select(LatentClass, Age, prob)
    
    # Generate the visualization of the trajectories
    trajectory_plot <- ggplot(plot_data, aes(x = Age, y = prob, color = LatentClass, group = LatentClass)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      # X-axis breaks now reflect Age 0 to 4
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

# 8. CLASS ADJUDICATION (Assigning Individuals to Classes) ---------------------
# This section extracts the class assignment (Most Likely Class) for each
# individual from the best-fitting model.

if (exists("best_model") && !is.null(best_model$savedata)) {
  
  message(glue("\nExtracting class assignments for the {best_classes}-class model..."))
  
  # === ROBUST COLUMN DETECTION ===
  # Inspect the actual column names in savedata to find the class membership column.
  savedata_cols <- names(best_model$savedata)
  message(glue("Columns in savedata: {paste(savedata_cols, collapse = ', ')}"))
  
  # Look for potential class column names: "C" or "MLCC"
  potential_cols <- c("C", "MLCC")
  class_col_name <- savedata_cols[toupper(savedata_cols) %in% potential_cols]
  
  if (length(class_col_name) == 0) {
    stop(glue("The Class Membership column (neither 'C' nor 'MLCC') was found in the savedata. Please check the 'Columns in savedata' message above."))
  }
  
  # Use the first match found (prioritizing the order in savedata_cols)
  class_col_name <- class_col_name[1]
  message(glue("Using column '{class_col_name}' for class membership."))
  
  class_assignments <- best_model$savedata %>%
    # Use all_of() with the detected name to handle case sensitivity safely
    select(USUBJID, Assigned_Class = all_of(class_col_name), starts_with("CPROB")) %>%
    mutate(USUBJID = as.character(USUBJID)) # Ensure ID matches original format if needed
  
  # Display the first few rows
  print(head(class_assignments))
  
  # Calculate Average Posterior Probabilities (AvePP)
  # This is a measure of classification quality. Ideally > 0.70 for all classes.
  avepp <- class_assignments %>%
    group_by(Assigned_Class) %>%
    summarise(
      Count = n(),
      AvePP = mean(get(paste0("CPROB", unique(Assigned_Class))))
    ) %>%
    mutate(Proportion = Count / sum(Count))
  
  cat("\nClassification Quality (Average Posterior Probabilities):\n")
  print(avepp)
  
  # Merge back with original data if needed for further analysis
  # df_adjudicated <- left_join(df_wide, class_assignments, by = "USUBJID")
  
} else {
  message("No savedata found in the best model object. Ensure SAVEDATA = ...; was run.")
}


# 9. EXPORT FILES FOR MANUAL MPLUS RUN -----------------------------------------
# This section exports a clean .dat file and a standalone .inp file for the 
# 3-class model. This allows you to run the analysis manually in Mplus.

# Create a separate folder for manual files
manual_output_dir <- here("Trajectory_Manual_Export_Nomum")
if (!dir.exists(manual_output_dir)) dir.create(manual_output_dir)

manual_dat_file <- here("Trajectory_Manual_Export_Nomum", "columinate_lcga_nomum_data.dat")
manual_inp_file <- here("Trajectory_Manual_Export_Nomum", "manual_3_class_nomum_model.inp")

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
# Updated to match the successful manual run syntax (MISSING and SERIES)
manual_inp_text <- glue("
TITLE: Manual 3-Class LCGA Trajectory Analysis (No Mother);

DATA: FILE = columinate_lcga_nomum_data.dat;

VARIABLE: 
  NAMES = {paste(names(df_wide), collapse = ' ')};
  IDVARIABLE = USUBJID;
  CLASSES = c(3);
  
  USEVARIABLES = {paste(names(df_wide), collapse = ' ')};
  CATEGORICAL = {mplus_varnames};
  MISSING = .;

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

PLOT: TYPE = PLOT3; SERIES = {mplus_series_string};

SAVEDATA: FILE = savedata_manual.dat; SAVE = CPROB;
")

writeLines(manual_inp_text, manual_inp_file)

message(glue("\nMANUAL EXPORT COMPLETE:\nData: {manual_dat_file}\nInput: {manual_inp_file}\nYou can now open {manual_inp_file} in Mplus and run it directly."))

# Clean up temporary Mplus data files
# COMMENTED OUT TO SAVE THE .DAT FILES AS REQUESTED
# file.remove(list.files(here("Trajectory_Output_Nomum"), pattern = "\\.dat$", full.names = TRUE))