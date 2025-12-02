# R Script for Data Processing and Preparation
# -------------------------------------------
# This script loads the raw screening data, applies inclusion/exclusion
# criteria, cleans and recodes questionnaire items, calculates total scores,
# and saves a final analysis-ready dataset.
#
# Author: Marothi Peter LETSOALO (Generative AI improved - Gemini Pro)
# GenAI Purpose: Used to document the script and optimize the process
# Date: September 11, 2025
# -------------------------------------------

# -----------------------------------------------------------------------------
# Step 1: Setup - Load necessary packages
# -----------------------------------------------------------------------------
# Ensure you have these packages installed: install.packages(c("tidyverse", "here", "labelled"))
library(tidyverse)
library(here) # For robust file path management
library(labelled) # For handling variable labels

cat("Step 1: Packages loaded successfully.\n")

# Helper function to calculate the mode
calculate_mode <- function(x) {
  # Exclude NA values from mode calculation
  x_clean <- na.omit(x)
  if (length(x_clean) == 0) return(NA)
  ux <- unique(x_clean)
  ux[which.max(tabulate(match(x, ux)))]
}

# -----------------------------------------------------------------------------
# Step 2: Data Loading
# -----------------------------------------------------------------------------
# Load the RData file
load(here("data_management/transformed_data", "dt_analysis_screening.RData"))

# Check if the data exists in the environment before proceeding
if (!exists("dt_analysis_screening")) {
  stop("The raw dataset 'dt_analysis_screening' was not found in the environment. Please load it first.")
}

cat("Step 2: Raw data from .RData file loaded successfully.\n")

# -----------------------------------------------------------------------------
# Step 3: Disposition - Identify participants for exclusion and Export
# -----------------------------------------------------------------------------
# This section creates a disposition table that flags participants who do not
# meet the study's inclusion criteria.

dt_disposition_ <-
  dt_analysis_screening %>%
  mutate(
    # --- Calculate missing counts first ---
    ssq_missing_count = rowSums(is.na(select(., starts_with("SSQ"))) | select(., starts_with("SSQ")) == "Prefer not to answer"),
    phq_missing_count = rowSums(is.na(select(., starts_with("PHQ90"))) | select(., starts_with("PHQ90")) %in% c("Prefer not to answer", "Not applicable")),
    
    # --- Define exclusion flags ---
    EXCAGE = AGE < 17 | AGE > 24 | is.na(AGE),
    EXCVST = is.na(VISITDT),
    EXCNMV_SSQ = if_all(starts_with("SSQ"), ~ .x %in% c(NA, "Prefer not to answer")),
    EXCNMV_PHQ9 = if_all(starts_with("PHQ90"), ~ .x %in% c(NA, "Prefer not to answer", "Not applicable")),
    # New flags for excessive missingness (>20%)
    EXCMIS_SSQ = ssq_missing_count >= 3,
    EXCMIS_PHQ9 = phq_missing_count >= 2
  ) |>
  rowwise() |>
  mutate(
    # Create a consolidated reason for exclusion
    EXCREAS = paste(
      c(
        if (EXCAGE) "Age not 17-24",
        if (EXCVST) "Missing Visit Date or Missing Age",
        if (EXCNMV_SSQ) "All SSQ items are missing",
        if (EXCNMV_PHQ9) "All PHQ-9 items are missing",
        if (EXCMIS_SSQ && !EXCNMV_SSQ) "Excessive missing data on SSQ-14 (>= 3 items)",
        if (EXCMIS_PHQ9 && !EXCNMV_PHQ9) "Excessive missing data on PHQ-9 (>= 2 items)"
      ),
      collapse = "; "
    )
  ) |>
  ungroup() |>
  mutate(
    # Create a final EXCLUDE flag and clean up the reason text
    EXCLUDE = EXCREAS != "",
    EXCREAS = if_else(EXCREAS == "", "Not Excluded", EXCREAS)
  ) |> 
  set_variable_labels(
      EXCAGE = "Age not 17-24",
      EXCVST = "Missing Visit Date or Missing Age",
      EXCNMV_SSQ = "All SSQ items are missing",
      EXCNMV_PHQ9 = "All PHQ-9 items are missing",
      EXCMIS_SSQ = "Excessive missing data on SSQ-14 (>= 3 items)",
      EXCMIS_PHQ9 = "Excessive missing data on PHQ-9 (>= 2 items)"
    )

cat("Step 3: Disposition logic updated and applied successfully.\n")
cat("Summary of Exclusions:\n")
print(table(dt_disposition_$EXCREAS))

# --- Export the disposition table ---
# Create the output directory if it doesn't exist
if (!dir.exists(here("statistical_analysis/data/adam"))) {
  dir.create(here("statistical_analysis/data/adam"), recursive = TRUE)
}
dt_disposition<-dt_disposition_ |>
  select(USUBJID,VISITDT,STUDYDY, contains("EXC"))
save(dt_disposition,
     file = here("statistical_analysis/data/adam", "dt_disposition.RData"))
cat("\nStep 3b: Disposition table 'dt_disposition.RData' has been exported to statistical_analysis/data/adam.\n")


# -----------------------------------------------------------------------------
# Step 4: Data Cleaning, Recoding, Imputation, and Score Calculation
# -----------------------------------------------------------------------------
# This section converts factors to numeric, imputes minimal missingness,
# and calculates final total scores.

dt_scored <- dt_disposition_ |>
  # --- Recode SSQ-14 items to numeric (0/1) ---
  mutate(across(
    starts_with("SSQ"),
    ~ case_when(
      . == "Yes" ~ 1,
      . == "No" ~ 0,
      TRUE ~ NA_real_ # Handles NA, "Prefer not to answer", etc.
    ),
    .names = "{.col}N"
  )) |>
  # --- Recode PHQ-9 items to numeric (0-3) ---
  mutate(across(
    starts_with("PHQ90"),
    ~ case_when(
      . == "Absolutely not" | . == "Not at all" ~ 0,
      . == "Several days" ~ 1,
      . == "More than half days" | . == "More than half the days" ~ 2,
      . == "Almost daily" | . == "Nearly every day" ~ 3,
      TRUE ~ NA_real_
    ),
    .names = "{.col}N"
  )) |>
  # --- Perform row-wise modal imputation for minimal missingness ---
  rowwise() |>
  mutate(
    # Impute SSQ items
    ssq_mode = calculate_mode(c_across(ends_with("N") & starts_with("SSQ"))),
    across(ends_with("N") & starts_with("SSQ"), ~if_else(is.na(.), ssq_mode, .)),
    # Impute PHQ items
    phq_mode = calculate_mode(c_across(ends_with("N") & starts_with("PHQ90"))),
    across(ends_with("N") & starts_with("PHQ90"), ~if_else(is.na(.), phq_mode, .))
  ) |>
  # --- Calculate Total Scores after imputation ---
  mutate(
    SSQSCR = sum(c_across(ends_with("N") & starts_with("SSQ"))),
    PHQSCR = sum(c_across(ends_with("N") & starts_with("PHQ90")))
  ) |>
  ungroup() # End rowwise operations

cat("\nStep 4: Item recoding, modal imputation, and score calculation complete.\n")


# -----------------------------------------------------------------------------
# Step 5: Create Final Analysis-Ready Dataset and Add Labels
# -----------------------------------------------------------------------------

# Combine all processed data into a final dataframe for analysis
analysis_ready_data <- dt_scored |>
  # Apply the exclusion criteria defined in Step 3
  filter(!EXCLUDE) |>
  # Select the final set of variables for modeling and analysis
  select(
    USUBJID,
    AGE,
    SEX,
    FI01,
    FI02,
    ends_with("N"), # All the numeric SSQ and PHQ items
    SSQSCR,
    PHQSCR,
    # Also include original categorical items as requested
    starts_with("SSQ"),
    starts_with("PHQ")
  ) |>
  select(-contains("_"))

# --- Preserve and Modify Variable Labels ---
# This loop runs on the final data object to prevent labels from being stripped
# by intermediate dplyr operations.
ssq_phq_vars <- names(dt_analysis_screening)[grepl("^SSQ|^PHQ90", names(dt_analysis_screening))]
numeric_vars <- paste0(ssq_phq_vars, "N")

# Loop through, copy the original label, and append "(Numeric)"
for (i in seq_along(ssq_phq_vars)) {
  original_var <- ssq_phq_vars[i]
  numeric_var <- numeric_vars[i]
  # Check if the numeric variable exists in the final dataset
  if(numeric_var %in% names(analysis_ready_data)) {
    original_label <- var_label(dt_analysis_screening[[original_var]])
    # Handle cases where the original label might be NULL
    if (is.null(original_label)) {
      original_label <- original_var
    }
    # Assign the new label to the column in the final dataframe
    var_label(analysis_ready_data[[numeric_var]]) <- paste(original_label, "(Numeric)")
  }
}

# Label the total score variables
var_label(analysis_ready_data$SSQSCR) <- "SSQ-14 Total Score"
var_label(analysis_ready_data$PHQSCR) <- "PHQ-9 Total Score"


cat("\nVariable labels added successfully.\n")

# Display the structure of the final dataset
cat("\n\nStructure of the final analysis-ready dataset:\n")
str(analysis_ready_data)

# Display the first few rows of the final dataset
cat("\n\nFirst 6 rows of the final analysis-ready dataset:\n")
print(head(analysis_ready_data))

# Display a label for one of the new numeric variables to confirm success
cat("\n\nChecking label for SSQ01N:\n")
print(var_label(analysis_ready_data$SSQ01N))

# Display a label for one of the new score variables to confirm success
cat("\n\nChecking label for PHQSCR:\n")
print(var_label(analysis_ready_data$PHQSCR))


# Save the final dataset in both .RData and .csv formats
save(analysis_ready_data, file = here("statistical_analysis/data/adam", "analysis_ready_data.RData"))
write.csv(analysis_ready_data, here("statistical_analysis/data/adam", "analysis_ready_data.csv"), row.names = FALSE)

cat("\nStep 5: 'analysis_ready_data.RData' and 'analysis_ready_data.csv' have been successfully saved to data/adam.\n")
cat("\nScript finished.\n")

