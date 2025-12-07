# Script 04: Demographics and Item-Level Analysis
# ==============================================================================
# PURPOSE: 
# 1. Generate Table 1: Participant Characteristics (Age, Sex, PHQ-9 Status + SES/Health).
# 2. Generate Table 2: SSQ-14 Item Endorsement by Depression Status & 
#    Item-Total Correlations (Point-Biserial) with PHQ-9 and SSQ-14 Total.
# ==============================================================================

# Load necessary libraries
if (!require("here")) install.packages("here")
if (!require("dplyr")) install.packages("dplyr")
if (!require("gtsummary")) install.packages("gtsummary")
if (!require("gt")) install.packages("gt")
if (!require("psych")) install.packages("psych") # For point-biserial correlation

library(here)
library(dplyr)
library(gtsummary)
library(gt)
library(psych)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------

# Load Main Data
data_path <- here("../../_private_use/data_management/data/co-luminate/adam/dt_psychometric.RData")

if (file.exists(data_path)) {
  load(file = data_path)
  if (!exists("dt_main")) dt_main <- get(ls()[1]) 
  message("Data loaded successfully.")
} else {
  stop(paste("Data file not found at:", data_path))
}

# Define Variables
ssq14_items <- paste0("SSQ", sprintf("%02d", 1:14))
phq_total_var <- "PHQSCR" 
ssq_total_var <- "SSQSCR" # New variable for Item-Total Correlation
age_var <- "AGE"
sex_var <- "SEX"

# Additional Demographic/Clinical Variables
extra_demog_vars <- c("URBANCAT", "SCHOOL", "EMPLOY", "FOODSEC", "GOVTGRNT", 
                      "ORPHSTAT", "HIVSTRESC", "VIOLENCE", "DRNKALC")

# Items to highlight (Excluded in SSQ-10)
excluded_items <- c("SSQ04", "SSQ05", "SSQ06", "SSQ07")

# -----------------------------------------------------------------------------
# 2. Data Preparation
# -----------------------------------------------------------------------------

# Check which extra vars exist in dataset to avoid errors
existing_extras <- intersect(extra_demog_vars, names(dt_main))
if(length(existing_extras) < length(extra_demog_vars)) {
  warning("Some requested demographic variables are missing from the dataset.")
}

# Create Analysis Dataset
df_demog <- dt_main %>%
  select(all_of(c(age_var, sex_var, phq_total_var, ssq_total_var, ssq14_items, existing_extras))) %>%
  mutate(
    # Create Age Groups
    Age_Group = case_when(
      .data[[age_var]] < 20 ~ "17-19",
      .data[[age_var]] >= 20 ~ "20-24",
      TRUE ~ NA_character_
    ),
    # Create Depression Status (PHQ-9 >= 10)
    Depression_Status = factor(ifelse(.data[[phq_total_var]] >= 10, "Depressed (PHQ-9 >= 10)", "Non-Depressed (< 10)"),
                               levels = c("Non-Depressed (< 10)", "Depressed (PHQ-9 >= 10)")),
    
    # Adjust Employment: Students are 'Not Applicable' for strict unemployment definition
    EMPLOY = if ("SCHOOL" %in% names(.) & "EMPLOY" %in% names(.)) {
      case_when(
        SCHOOL == "Inschool" & EMPLOY == "Not employed" ~ "Student (Not Applicable)",
        is.na(EMPLOY) ~ "Status Unknown",
        TRUE ~ as.character(EMPLOY)
      )
    } else {
      EMPLOY
    }
  )

# Ensure SSQ items are numeric 0/1 for correlation, but Factors for Summary Table
# We keep a numeric copy for correlations
df_numeric <- df_demog
df_numeric[ssq14_items] <- lapply(df_numeric[ssq14_items], function(x) {
  if (is.factor(x) || is.character(x)) as.numeric(x == "Yes") else x
})

# -----------------------------------------------------------------------------
# 3. Table 1: Participant Characteristics (Stacked Groups)
# -----------------------------------------------------------------------------

# Define groups of variables (intersect with available data)
vars_demog <- intersect(c("Age_Group", sex_var, "URBANCAT"), names(df_demog))
vars_ses   <- intersect(c("SCHOOL", "EMPLOY", "FOODSEC", "GOVTGRNT"), names(df_demog))
vars_fam   <- intersect(c("ORPHSTAT"), names(df_demog))
vars_hlth  <- intersect(c("HIVSTRESC", "VIOLENCE", "DRNKALC"), names(df_demog))

# Define labels as a NAMED LIST (using = not ~) to allow robust subsetting
var_labels <- list(
  Age_Group = "Age Group (Years)",
  SEX = "Biological Sex",
  URBANCAT = "Household Location",
  SCHOOL = "Currently in School",
  EMPLOY = "Employment Status",
  FOODSEC = "Food Insecurity",
  GOVTGRNT = "Receives Government Grant",
  ORPHSTAT = "Orphanhood Status",
  HIVSTRESC = "HIV Status",
  VIOLENCE = "Experienced Violence",
  DRNKALC = "Ever Drank Alcohol"
)

# Function to generate sub-tables
make_tbl <- function(vars) {
  if(length(vars) == 0) return(NULL)
  
  # Ensure we only pick labels for variables that actually exist in 'vars'
  relevant_labels <- var_labels[names(var_labels) %in% vars]
  
  df_demog %>%
    select(all_of(vars), Depression_Status) %>%
    tbl_summary(
      by = Depression_Status,
      label = relevant_labels, 
      statistic = list(all_categorical() ~ "{n} ({p}%)"),
      missing = "ifany"
    ) %>%
    add_overall() %>%
    italicize_levels() # Italicize levels, Labels remain Normal (default)
}

# Create sub-tables
t1_demog <- make_tbl(vars_demog)
t1_ses   <- make_tbl(vars_ses)
t1_fam   <- make_tbl(vars_fam)
t1_hlth  <- make_tbl(vars_hlth)

# Stack them
# Filter out NULL tables if variables were missing
tbl_list <- list(t1_demog, t1_ses, t1_fam, t1_hlth)
group_headers <- c("Demographics", "Socio-economics", "Family/Vulnerability", "Health & Risk")

# Keep only valid tables and headers
valid_idx <- !sapply(tbl_list, is.null)
tbl_list <- tbl_list[valid_idx]
group_headers <- group_headers[valid_idx]

table1 <- tbl_stack(tbl_list, group_header = group_headers) %>%
  modify_caption("**Participant Characteristics**") %>%
  # Apply Bold to Group Headers (using gtsummary native function before as_gt)
  modify_table_styling(
    columns = label,
    rows = label %in% group_headers,
    text_format = "bold"
  ) %>%
  as_gt()%>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups(groups = everything())
  ) %>%
  tab_source_note(
    source_note = "Note: Unemployment is defined as being willing/available to work but jobless. Participants currently in school are classified as 'Student (Not Applicable)' to align with labour force participation definitions."
  ) %>%
  tab_options(table.font.names = "Arial")

# print(table1) # Suppressed

# -----------------------------------------------------------------------------
# 4. Table 2: Item Endorsement and Correlations
# -----------------------------------------------------------------------------

# A. Calculate Frequencies by Group using gtsummary
table2_freq <- df_demog %>%
  select(all_of(ssq14_items), Depression_Status) %>%
  tbl_summary(
    by = Depression_Status,
    statistic = list(all_categorical() ~ "{n} ({p}%)"), # Show N and %
    label = list(
      SSQ01 ~ "SSQ01: Thinking too much",
      SSQ02 ~ "SSQ02: Difficulty concentrating",
      SSQ03 ~ "SSQ03: Irritability",
      SSQ04 ~ "SSQ04: Nightmares",
      SSQ05 ~ "SSQ05: Hallucinations",
      SSQ06 ~ "SSQ06: Stomach ache",
      SSQ07 ~ "SSQ07: Frightened easily",
      SSQ08 ~ "SSQ08: Poor sleep",
      SSQ09 ~ "SSQ09: Crying",
      SSQ10 ~ "SSQ10: Feeling tired",
      SSQ11 ~ "SSQ11: Suicidal thoughts",
      SSQ12 ~ "SSQ12: Generally unhappy",
      SSQ13 ~ "SSQ13: Work or School Lagging Behind",
      SSQ14 ~ "SSQ14: Difficulty deciding"
    ),
    missing = "no" # Exclude missing rows from denominator usually
  ) %>%
  add_overall() %>%
  bold_labels() 

# B. Calculate Point-Biserial Correlations
# 1. Item vs PHQ-9 (Concurrent Validity)
# 2. Item vs SSQ-14 (Internal Consistency)
cor_results <- data.frame(Item = ssq14_items, 
                          r_pb_phq = NA, p_val_phq = NA,
                          r_pb_ssq = NA, p_val_ssq = NA)

for (item in ssq14_items) {
  # PHQ-9 Correlation
  res_phq <- cor.test(df_numeric[[item]], df_numeric[[phq_total_var]])
  cor_results[cor_results$Item == item, "r_pb_phq"] <- res_phq$estimate # Store raw, format later
  
  p_phq <- res_phq$p.value
  cor_results[cor_results$Item == item, "p_val_phq"] <- if(p_phq < 0.001) "<0.001" else sprintf("%.3f", p_phq)
  
  # SSQ-14 Correlation
  res_ssq <- cor.test(df_numeric[[item]], df_numeric[[ssq_total_var]])
  cor_results[cor_results$Item == item, "r_pb_ssq"] <- res_ssq$estimate # Store raw, format later
  
  p_ssq <- res_ssq$p.value
  cor_results[cor_results$Item == item, "p_val_ssq"] <- if(p_ssq < 0.001) "<0.001" else sprintf("%.3f", p_ssq)
}

# Create a tibble for merging with gtsummary
# Note: gtsummary merges by "variable" name
cor_tibble <- cor_results %>%
  rename(variable = Item) %>%
  mutate(
    # Force 3 digits for correlations using sprintf
    Corr_PHQ_Label = paste0(sprintf("%.3f", r_pb_phq), " (p", ifelse(p_val_phq=="<0.001", "", "="), p_val_phq, ")"),
    Corr_SSQ_Label = paste0(sprintf("%.3f", r_pb_ssq), " (p", ifelse(p_val_ssq=="<0.001", "", "="), p_val_ssq, ")")
  ) %>%
  select(variable, Corr_PHQ_Label, Corr_SSQ_Label)

# C. Merge and Finalize Table 2
table2_final <- table2_freq %>%
  modify_table_body(
    ~ .x %>%
      left_join(cor_tibble, by = "variable")
  ) %>%
  modify_header(
    label = "**SSQ-14 Item**",
    stat_0 = "**Total Sample**\n(N = {N})",
    stat_1 = "**Non-Depressed**\n(PHQ < 10)\n(N = {n})",
    stat_2 = "**Depressed**\n(PHQ >= 10)\n(N = {n})",
    Corr_PHQ_Label = "**Correlation**\n(with PHQ-9)",
    Corr_SSQ_Label = "**Correlation**\n(with SSQ-14)"
  ) %>%
  modify_fmt_fun(
    c(Corr_PHQ_Label, Corr_SSQ_Label) ~ function(x) x # Prevent auto-formatting errors
  ) %>%
  as_gt() %>%
  # Update Header: Removed "Table 2." prefix
  tab_header(
    title = md("**Distribution of SSQ-14 Items by Depression Status and Correlations**")
  ) %>%
  # Highlight the excluded items
  tab_style(
    style = list(
      cell_fill(color = "#F9E3E3"), # Light Red/Pink background
      cell_text(style = "italic")
    ),
    locations = cells_body(
      rows = variable %in% excluded_items
    )
  ) %>%
  tab_source_note(
    source_note = "Note: Data presented as n (%). Correlations are point-biserial coefficients between binary items and continuous total scores. Items highlighted in red (SSQ04-07) are those excluded from the theoretical SSQ-10 subset."
  ) %>%
  tab_options(table.font.names = "Arial")

# print(table2_final) # Suppressed

# -----------------------------------------------------------------------------
# 5. Save Tables
# -----------------------------------------------------------------------------
output_dir <- here("statistical_analysis/output/tables")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

gtsave(table1, filename = file.path(output_dir, "Table1_Demographics.html"))
gtsave(table2_final, filename = file.path(output_dir, "Table2_ItemAnalysis.html"))

# Also save as objects for RMarkdown/Quarto inclusion
saveRDS(list(table1 = table1, table2 = table2_final), 
        file = file.path(output_dir, "04_Demographic_Tables.rds"))

message(paste("Tables saved to:", output_dir))