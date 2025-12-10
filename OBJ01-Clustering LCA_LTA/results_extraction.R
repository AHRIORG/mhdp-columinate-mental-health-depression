# ==============================================================================
# PROJECT: CO-LUMINATE - LCGA Result Extraction
# OBJECTIVE: Extract plots, fit stats, and class assignments for Quarto/RMarkdown.
# INPUT: Mplus .out files
# OUTPUT: .rds file containing list of results for MULTIPLE models.
# ==============================================================================

# 1. SETUP ---------------------------------------------------------------------
library(tidyverse)
library(MplusAutomation)
library(glue)
library(here)
library(gt) # For beautiful tables

# === USER CONFIGURATION: SELECT YOUR MODELS HERE ===
# You can now specify multiple classes to extract (e.g., c(3, 4))
target_classes <- c(3, 4) 
output_dir     <- here("Trajectory_Output_Nomum") 
# ==================================================

# Initialize master results list
final_output <- list(
  fit_comparison_table = NULL,
  models = list() # Will hold specific results for each k
)

# 2. GLOBAL FIT STATISTICS TABLE (Comparison of All Models) --------------------
message("Reading all models in directory to generate comparison table...")

# Read all .out files in the directory to construct the comparison
all_models_list <- readModels(output_dir, quiet = TRUE, recursive = FALSE)

if (length(all_models_list) > 0) {
  
  # Extract summaries into a single dataframe
  fit_data <- map(all_models_list, "summaries") %>%
    bind_rows() %>%
    select(
      Classes = NLatentClasses,
      LogLikelihood = LL,
      AIC, BIC, aBIC, 
      Entropy,
      VLMR_p = T11_VLMR_PValue,
      BLRT_p = BLRT_PValue
    ) %>%
    arrange(Classes)
  
  # Create gt table comparing all models
  tbl_fit <- fit_data %>%
    gt() %>%
    tab_header(
      title = "LCGA Model Fit Comparison",
      subtitle = "Information Criteria and Classification Accuracy across all fitted models"
    ) %>%
    cols_label(
      Classes = "Latent Classes",
      LogLikelihood = "Log-Likelihood",
      aBIC = "Sample-Size Adj. BIC",
      VLMR_p = "VLMR p-value",
      BLRT_p = "BLRT p-value"
    ) %>%
    fmt_number(
      columns = c(LogLikelihood, AIC, BIC, aBIC),
      decimals = 2
    ) %>%
    fmt_number(
      columns = c(Entropy),
      decimals = 3
    ) %>%
    fmt_number(
      columns = c(VLMR_p, BLRT_p),
      decimals = 3
    ) %>%
    sub_missing(
      columns = everything(),
      missing_text = "-"
    ) %>%
    tab_footnote(
      footnote = "Lower values indicate better model fit.",
      locations = cells_column_labels(columns = c(AIC, BIC, aBIC))
    ) %>%
    tab_footnote(
      footnote = "Values > 0.80 indicate good separation between classes.",
      locations = cells_column_labels(columns = c(Entropy))
    ) %>%
    # Highlight rows corresponding to the user-selected models
    tab_style(
      style = list(
        cell_fill(color = "#e6f3ff"),
        cell_text(weight = "bold")
      ),
      locations = cells_body(rows = Classes %in% target_classes)
    ) %>%
    tab_options(
      table.width = pct(100)
    )
  
  final_output$fit_comparison_table <- tbl_fit
}

# 3. ITERATE THROUGH TARGET CLASSES --------------------------------------------
# Function to extract results for a single k
extract_single_model <- function(k, dir) {
  
  message(glue("\n--- Processing {k}-Class Model ---"))
  
  model_filename <- glue("exposure_data_{target_class}.out")
  model_path     <- file.path(dir, model_filename)
  
  if (!file.exists(model_path)) {
    warning(glue("File not found: {model_filename}. Skipping."))
    return(NULL)
  }
  
  # Read model with "what='all'" to ensure everything is parsed
  model_res <- readModels(model_path, quiet = TRUE, what = "all")
  
  # Initialize list for this specific model
  model_artifacts <- list(
    k = k,
    plot = NULL,
    avepp_table = NULL,
    data = NULL
  )
  
  # --- A. Trajectory Plot Extraction ---
  if (!is.null(model_res$parameters$unstandardized)) {
    message(glue("[{k}-Class] Parameters found. Generating plot..."))
    
    plot_data <- model_res$parameters$unstandardized %>%
      filter(paramHeader == "Thresholds") %>%
      filter(grepl("LWOM|LIVEDWOM", param, ignore.case = TRUE))
    
    if (nrow(plot_data) > 0) {
      plot_data <- plot_data %>%
        mutate(
          var_clean   = gsub("(\\$1|\\.Cat\\.1)$", "", param),
          Age         = readr::parse_number(var_clean),
          LatentClass = factor(LatentClass),
          prob        = 1 - pnorm(est) 
        ) %>%
        select(LatentClass, Age, prob)
      
      p <- ggplot(plot_data, aes(x = Age, y = prob, color = LatentClass, group = LatentClass)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3) +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        scale_x_continuous(breaks = unique(plot_data$Age)) +
        labs(
          title    = glue("Latent Trajectories ({k} Classes)"),
          subtitle = "Probability of Living Without Mother (Ages 0-4)",
          y        = "Probability",
          x        = "Child Age (Years)",
          color    = "Class"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "bottom",
          plot.title = element_text(face = "bold")
        )
      
      model_artifacts$plot <- p
      
      # Save PNG
      ggsave(here(dir, glue("plot_{k}_class_trajectory.png")), 
             plot = p, width = 8, height = 6, bg = "white")
    } else {
      warning(glue("[{k}-Class] Threshold parameters not found matching 'LWOM' or 'LIVEDWOM'. Plot skipped."))
    }
  } else {
    warning(glue("[{k}-Class] 'unstandardized' parameters are NULL. Check Mplus output file for errors."))
  }
  
  # --- B. Diagnostics & Adjudication ---
  if (!is.null(model_res$savedata)) {
    message(glue("[{k}-Class] SAVEDATA found. Extracting diagnostics..."))
    
    cols <- names(model_res$savedata)
    # Print columns to help debugging if extraction fails
    # print(glue("Columns: {paste(cols, collapse=', ')}"))
    
    class_col <- cols[toupper(cols) %in% c("C", "MLCC")][1]
    
    if (!is.na(class_col)) {
      message(glue("[{k}-Class] Found class column: {class_col}"))
      
      # AvePP Table
      avepp_data <- model_res$savedata %>%
        rename(Assigned_Class = all_of(class_col)) %>%
        group_by(Assigned_Class) %>%
        summarise(
          Count = n(),
          AvePP = mean(get(paste0("CPROB", unique(Assigned_Class))))
        ) %>%
        mutate(Proportion = Count / sum(Count))
      
      tbl_avepp <- avepp_data %>%
        gt() %>%
        tab_header(
          title = glue("Classification Diagnostics ({k}-Class)"),
          subtitle = "Average Posterior Probabilities (AvePP)"
        ) %>%
        fmt_number(columns = c(AvePP), decimals = 3) %>%
        fmt_percent(columns = c(Proportion), decimals = 1) %>%
        tab_style(
          style = list(cell_text(weight = "bold")),
          locations = cells_body(columns = c(AvePP), rows = AvePP > 0.7)
        )
      
      model_artifacts$avepp_table <- tbl_avepp
      
      # Adjudicated Data
      adjudicated_data <- model_res$savedata %>%
        select(USUBJID, 
               LCGA_Class = all_of(class_col), 
               starts_with("CPROB")) %>%
        mutate(
          USUBJID = as.character(USUBJID),
          LCGA_Model = paste0(k, "_Class")
        )
      
      model_artifacts$data <- adjudicated_data
      
      # Save CSV
      write_csv(adjudicated_data, here(dir, glue("adjudicated_{k}class.csv")))
    } else {
      warning(glue("[{k}-Class] Could not find Class column (C or MLCC) in savedata."))
    }
  } else {
    warning(glue("[{k}-Class] 'savedata' is NULL. Mplus may not have generated the save file or readModels couldn't parse it."))
  }
  
  return(model_artifacts)
}

# Run extraction for all target classes
model_results_list <- map(target_classes, ~ extract_single_model(.x, output_dir))
names(model_results_list) <- paste0("k", target_classes)

# Combine into final object
final_output$models <- model_results_list

# 4. SAVE LIST OBJECT (.rds) ---------------------------------------------------
rds_filename <- here(output_dir, "lcga_results_multimodel.rds")
saveRDS(final_output, file = rds_filename)

message(glue("\nSUCCESS: Multi-model results bundled and saved to:\n{rds_filename}"))
message("Structure: final_output$fit_comparison_table, final_output$models$k3$plot, etc.")