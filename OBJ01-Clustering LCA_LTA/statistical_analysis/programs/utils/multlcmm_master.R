# ==============================================================================
# PROJECT: CO-LUMINATE Trajectory Analysis (multlcmm)
# OBJECTIVE: Identify JOINT latent trajectories from multiple socioeconomic variables
# METHOD: multlcmm (Multivariate Latent Class Linear Mixed Models)
# OUTPUT: .rds list object containing tables, latent process plots, and data
# ==============================================================================

# 1. SETUP & LIBRARIES ---------------------------------------------------------
library(tidyverse)
library(lcmm)
library(here)
library(gt)
library(glue)

# Define output directory
output_dir <- here("Trajectory_Output_multlcmm")
if (!dir.exists(output_dir)) dir.create(output_dir)

models_dir <- here(output_dir, "model_objects")
if (!dir.exists(models_dir)) dir.create(models_dir)

# 2. CONFIGURATION -------------------------------------------------------------
# Define the set of variables to model JOINTLY
# Note: Order matters for the 'link' argument in multlcmm

multivariate_config <- list(
  # General Settings
  id_var      = "USUBJID",
  time_var    = "EXPAGE",
  time_range  = 1:2,
  num_classes = 2,
  rerun_models = FALSE,
  
  # List of Variables to include in the joint model
  variables = list(
    list(name = "LIVEDWOM", type = "binary",  link = "thresholds", label = "No Mother")#,
    #list(name = "LIVEDWOF", type = "binary",  link = "thresholds", label = "No Father")#,
    # list(name = "HHSES",    type = "ordinal", link = "thresholds", label = "Household SES"),
    # list(name = "EDUCMOTH", type = "ordinal", link = "thresholds", label = "Mother Edu"),
    # list(name = "HHCHLD4",  type = "continuous", link = "linear",  label = "Children <4y")
  )
)

# 3. GLOBAL DATA LOADING -------------------------------------------------------
data_path <- here("../../_private_use/data_management/data/co-luminate/adam/dt_childhood_exposure.RData")
if (!file.exists(data_path)) stop("Data file not found.")
load(data_path)
message(glue("Data loaded successfully from: {data_path}"))

# 4. ANALYSIS FUNCTION DEFINITION ----------------------------------------------
run_multivariate_lcga <- function(config, input_data, out_dir, mod_dir) {
  
  # --- Setup ---
  id_var     <- config$id_var
  time_var   <- config$time_var
  time_range <- config$time_range
  vars_info  <- config$variables
  
  # Extract names and links for the formula
  var_names <- map_chr(vars_info, "name")
  var_links <- map_chr(vars_info, "link")
  var_types <- map_chr(vars_info, "type")
  
  message(glue("\n========================================================"))
  message(glue("STARTING MULTIVARIATE LCGA"))
  message(glue("Variables: {paste(var_names, collapse = ', ')}"))
  message(glue("Links: {paste(var_links, collapse = ', ')}"))
  message(glue("========================================================"))
  
  # --- Data Preparation ---
  # multlcmm needs a dataframe where each outcome is a column, 
  # and rows are repeated for time points.
  
  # 1. Filter Time
  dt_prep <- input_data %>%
    filter(!!sym(time_var) %in% time_range) %>%
    mutate(
      ID_NUM = as.numeric(as.factor(!!sym(id_var))),
      TIME_NUM = as.numeric(!!sym(time_var))
    )
  
  # 2. Process Specific Variable Types (Factors -> Numeric)
  # multlcmm requires numeric inputs even for categorical/ordinal (1, 2, 3...)
  for (i in seq_along(vars_info)) {
    v <- vars_info[[i]]
    v_col <- v$name
    
    if (v$type %in% c("binary", "ordinal")) {
      # Convert factors to integer levels (1, 2, 3...)
      if (is.factor(dt_prep[[v_col]])) {
        dt_prep[[v_col]] <- as.numeric(dt_prep[[v_col]])
      } else {
        # Ensure 0/1 binary becomes 1/2 for thresholds link if needed, 
        # though lcmm handles 0/1. Standardizing to 1-based is safer for ordinal.
        # Check if min is 0
        if(min(dt_prep[[v_col]], na.rm=TRUE) == 0) {
          dt_prep[[v_col]] <- dt_prep[[v_col]] + 1
        }
      }
    } else {
      dt_prep[[v_col]] <- as.numeric(dt_prep[[v_col]])
    }
  }
  
  dt_prep <- as.data.frame(dt_prep)
  
  # --- Helper Functions ---
  calc_entropy <- function(m) {
    if (m$ng == 1) return(NA)
    if (!is.null(m$entropy)) return(m$entropy)
    # Fallback calc if needed
    if (!is.null(m$pprob)) {
      probs <- m$pprob[, grep("^prob", colnames(m$pprob))]
      probs[probs < 1e-9] <- 1e-9
      n <- nrow(probs)
      k <- ncol(probs)
      num <- sum(-apply(probs * log(probs), 1, sum))
      den <- n * log(k)
      return(1 - (num / den))
    }
    return(NA)
  }
  
  # --- Model Formula Construction ---
  # Fixed: Y1 + Y2 ... ~ 1 + TIME + I(TIME^2)
  outcome_str <- paste(var_names, collapse = " + ")
  fixed_fml   <- as.formula(paste(outcome_str, "~ 1 + TIME_NUM + I(TIME_NUM^2)"))
  random_fml  <- ~ 1 + TIME_NUM # Random intercept & slope on the LATENT PROCESS
  mixture_fml <- ~ 1 + TIME_NUM + I(TIME_NUM^2) # Classes differ in shape
  
  # Initialize Lists
  results_list <- list()
  temp_models  <- list()
  
  # --- Model Loop ---
  for (k in 1:config$num_classes) {
    message(glue("Processing k={k}..."))
    
    model_path <- here(mod_dir, glue("multlcga_k{k}.rds"))
    
    tryCatch({
      if (file.exists(model_path) && !config$rerun_models) {
        message("  -> Loading existing model...")
        m <- readRDS(model_path)
      } else {
        message("  -> Estimating new model...")
        
        if (k == 1) {
          m <- multlcmm(
            fixed = fixed_fml,
            random = random_fml,
            subject = "ID_NUM",
            data = dt_prep,
            link = var_links,
            ng = 1,
            verbose = FALSE
          )
        } else {
          m <- multlcmm(
            fixed = fixed_fml,
            mixture = mixture_fml,
            random = random_fml,
            subject = "ID_NUM",
            data = dt_prep,
            link = var_links,
            ng = k,
            B = temp_models[[1]], # Initialize with 1-class estimates
            verbose = FALSE
          )
        }
        saveRDS(m, model_path)
      }
      
      if (k == 1) temp_models[[1]] <- m
      
      # --- Store Results Components ---
      k_res <- list(model_obj = m)
      
      # 1. Plot Latent Process (The underlying factor trajectory)
      # We predict the Latent Process (Ydep) over time
      dat_pred <- data.frame(TIME_NUM = seq(min(time_range), max(time_range), length.out = 50))
      
      # predictL computes the predicted mean latent process for each class
      pred_latent <- predictL(m, dat_pred, var.time = "TIME_NUM")
      
      plot_data <- pred_latent$pred %>%
        as.data.frame() %>%
        mutate(TIME = dat_pred$TIME_NUM)
      
      # Reshape
      if (k == 1) {
        # Single column usually named 'y_pred' or similar
        col_name <- names(plot_data)[1] 
        plot_data <- plot_data %>% 
          rename(Value = all_of(col_name)) %>%
          mutate(Class = "1")
      } else {
        plot_data <- plot_data %>%
          pivot_longer(cols = starts_with("class"), names_to = "Class", values_to = "Value") %>%
          mutate(Class = str_extract(Class, "[0-9]+"))
      }
      
      p <- ggplot(plot_data, aes(x = TIME, y = Value, color = Class, group = Class)) +
        geom_line(linewidth = 1.2) +
        labs(
          title = NULL,
          subtitle = glue("{k} Classes (Joint Latent Process)"),
          y = "Latent Process Level (Standardized)",
          x = "Time"
        ) +
        theme_minimal() + theme(legend.position = "bottom")
      
      ggsave(here(out_dir, glue("plot_mult_latent_k{k}.png")), plot = p, width=8, height=6)
      k_res$plot_latent <- p
      
      # 2. Diagnostics
      if (k > 1) {
        class_assignments <- m$pprob %>% 
          rename(ID_NUM = ID_NUM) %>% 
          select(ID_NUM, class, starts_with("prob"))
        
        avepp <- class_assignments %>% 
          group_by(class) %>% 
          summarise(Count = n(), AvePP = mean(get(paste0("prob", unique(class))))) %>% 
          mutate(Proportion = Count / sum(Count))
        
        k_res$diagnostics_table <- avepp %>%
          gt() %>%
          tab_header(title = glue("Diagnostics k={k}")) %>%
          fmt_number(columns = c(AvePP), decimals = 3, use_seps=FALSE) %>%
          fmt_percent(columns = c(Proportion), decimals = 1)
        
        # Merge back original IDs
        id_map <- dt_prep %>% select(ID_NUM, RAW_ID) %>% distinct()
        k_res$adjudicated_data <- class_assignments %>%
          left_join(id_map, by = "ID_NUM") %>%
          rename(!!sym(id_var) := RAW_ID) %>%
          select(!!sym(id_var), everything(), -ID_NUM)
      }
      
      results_list[[paste0("k", k)]] <- k_res
      
    }, error = function(e) message(glue("Error k={k}: {e$message}")))
  }
  
  # --- Summary Table ---
  fit_data <- map_df(names(results_list), function(n) {
    m <- results_list[[n]]$model_obj
    if (is.null(m)) return(NULL)
    tibble(Classes = m$ng, AIC = m$AIC, BIC = m$BIC, Entropy = calc_entropy(m), Converged = (m$conv == 1))
  })
  
  if (nrow(fit_data) > 0) {
    final_table <- fit_data %>% 
      arrange(Classes) %>%
      gt() %>%
      tab_header(title = "Multivariate LCGA Fit Comparison") %>%
      fmt_number(columns = c(AIC, BIC), decimals = 2, use_seps=FALSE) %>%
      fmt_number(columns = c(Entropy), decimals = 3, use_seps=FALSE) %>%
      sub_missing(missing_text = "-") %>%
      tab_style(style = list(cell_fill(color = "#e6f3ff"), cell_text(weight = "bold")), locations = cells_body(rows = BIC == min(BIC)))
    
    results_list$fit_summary_table <- final_table
  }
  
  # Save Master Object
  saveRDS(results_list, here(out_dir, "multivariate_lcga_results.rds"))
  message("Multivariate analysis complete.")
  return(results_list)
}

# 5. EXECUTION -----------------------------------------------------------------
results <- run_multivariate_lcga(multivariate_config, dt_childhood_exposure, output_dir, models_dir)