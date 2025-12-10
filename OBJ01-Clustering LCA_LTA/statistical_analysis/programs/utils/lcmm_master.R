# ==============================================================================
# PROJECT: CO-LUMINATE Trajectory Analysis (lcmm)
# OBJECTIVE: Flexible LCGA for Binary, Ordinal, Continuous, or Count data
# OUTPUT: .rds list object containing gt tables, ggplots, and model data
# ==============================================================================

# 1. SETUP & LIBRARIES ---------------------------------------------------------
library(tidyverse)
library(lcmm)
library(here)
library(gt)
library(glue)

# 2. ANALYSIS FUNCTION DEFINITION ----------------------------------------------
run_single_lcga <- function(config, input_data, out_dir, mod_dir) {
  
  # --- A. Unpack Configuration & Setup ---
  target_var  <- config$target_var
  id_var      <- config$id_var      
  time_var    <- config$time_var    
  var_type    <- config$var_type
  use_linear  <- config$use_linear_for_binary
  time_range  <- config$time_range  
  num_classes <- config$num_classes
  rerun       <- config$rerun_models
  
  # Robustness: Handle missing time_range
  if (is.null(time_range)) {
    if (!time_var %in% names(input_data)) stop(glue("Time variable '{time_var}' not found in data."))
    obs_times <- na.omit(input_data[[time_var]])
    time_range <- sort(unique(obs_times))
    message(glue("Notice: 'time_range' derived from data: {min(time_range)} to {max(time_range)}."))
  }
  
  # Determine Display Label
  if (!is.null(config$var_label)) {
    display_label <- config$var_label
  } else {
    data_label <- attr(input_data[[target_var]], "label")
    display_label <- if (!is.null(data_label)) data_label else target_var
  }
  
  # Determine Y-Axis Label
  if (!is.null(config$y_axis_label)) {
    y_label <- config$y_axis_label
  } else {
    y_label <- if (var_type == "binary") "Probability" else "Outcome Level"
  }
  
  # Determine X-Axis Label (Time)
  if (!is.null(config$time_var_label)) {
    time_label <- config$time_var_label
  } else {
    time_label <- time_var # Default to variable name if not provided
  }
  
  message(glue("\n========================================================"))
  message(glue("STARTING ANALYSIS FOR: {display_label} ({target_var})"))
  message(glue("Type: {var_type} | Time: {time_var} ({min(time_range)}-{max(time_range)}) | Max K: {num_classes}"))
  message(glue("========================================================"))
  
  # --- B. Data Preparation ---
  dt_long <- input_data %>%
    rename(
      RAW_ID = all_of(id_var),
      TIME   = all_of(time_var),
      Y_RAW  = all_of(target_var)
    ) %>%
    filter(TIME %in% time_range) %>%
    mutate(
      ID = as.numeric(as.factor(RAW_ID)),
      Y_OUTCOME = case_when(
        var_type == "binary" & use_linear & is.factor(Y_RAW) ~ as.numeric(Y_RAW) - 1,
        var_type == "binary" & use_linear & !is.factor(Y_RAW) ~ as.numeric(Y_RAW),
        var_type == "ordinal" & is.factor(Y_RAW) ~ as.numeric(Y_RAW),
        var_type == "ordinal" & !is.factor(Y_RAW) ~ as.numeric(Y_RAW), 
        var_type == "binary" & !use_linear ~ as.numeric(Y_RAW), 
        TRUE ~ as.numeric(Y_RAW)
      )
    ) %>%
    select(RAW_ID, ID, TIME, Y_OUTCOME) %>%
    arrange(ID, TIME) %>%
    as.data.frame()
  
  # Determine Model Method
  use_gaussian <- (var_type %in% c("continuous", "count")) || (var_type == "binary" & use_linear)
  method_label <- if(use_gaussian) "hlme (Gaussian)" else "lcmm (Thresholds)"
  
  # Define Interpretation Note based on Method/Type
  longitudinal_note <- if (use_gaussian) {
    if (var_type == "binary") "Linear changes in Probability (0-1)."
    else "Linear changes in the Outcome Unit."
  } else {
    "Effects on the Latent Process Scale (SDs)."
  }
  
  # --- C. Helper Functions ---
  calc_entropy <- function(m) {
    if (m$ng == 1) return(NA)
    if (!is.null(m$entropy)) return(m$entropy)
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
  
  # --- UPDATED: WIDE FORMAT ESTIMATES WITH GROUPING & DUPLICATE FIX ---
  extract_estimates_table <- function(m, title_label, note) {
    if (is.null(m)) return(NULL)
    
    coefs <- m$best
    n_coef <- length(coefs)
    
    # 1. Variance Matrix Handling
    if (is.null(m$V)) {
      message(glue("  [Debug] {title_label}: m$V is NULL. SEs cannot be computed."))
      return(NULL)
    }
    
    if (is.matrix(m$V)) {
      if (nrow(m$V) != n_coef) return(NULL)
      variances <- diag(m$V)
    } else {
      expected_tri_len <- n_coef * (n_coef + 1) / 2
      if (length(m$V) == expected_tri_len) {
        V_mat <- matrix(0, nrow = n_coef, ncol = n_coef)
        V_mat[upper.tri(V_mat, diag = TRUE)] <- m$V
        variances <- diag(V_mat)
      } else if (length(m$V) == n_coef) {
        variances <- m$V
      } else {
        return(NULL)
      }
    }
    
    se <- sqrt(variances)
    z <- coefs / se
    p_val <- 2 * (1 - pnorm(abs(z)))
    
    # 2. Section Parsing (Assign Types)
    n_prob <- m$N[1] 
    n_fix  <- m$N[2] 
    n_vc   <- m$N[3] 
    n_cor  <- m$N[4] 
    
    types <- c(
      rep("Class-Membership Model", n_prob),
      rep("Longitudinal Model", n_fix),
      rep("Random Effects & Residuals", n_vc + n_cor)
    )
    if(length(types) < n_coef) types <- c(types, rep("Other", n_coef - length(types)))
    
    param_names <- names(coefs)
    if (is.null(param_names)) param_names <- paste0("Param_", 1:n_coef)
    
    df_clean <- data.frame(
      Raw_Name = param_names,
      Estimate = coefs,
      SE = se,
      P = p_val,
      Type = types,
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        Class_Extract = str_extract(Raw_Name, " class[0-9]+$"),
        Class_Num = str_extract(Class_Extract, "[0-9]+"),
        
        Term = str_remove(Raw_Name, " class[0-9]+$"),
        
        # FIX: Ensure unique naming for variance components to prevent pivot duplicates
        Term_Clean = case_when(
          Term == "intercept" ~ "Intercept",
          Term == "TIME" ~ "Linear Slope",
          Term == "I(TIME^2)" ~ "Quadratic Slope",
          # Preserve suffix for varcov/stderr to distinguish multiple random effects
          str_detect(Term, "stderr") ~ str_replace(Term, "stderr", "Residual SE"),
          str_detect(Term, "varcov") ~ str_replace(Term, "varcov", "Random Effect"),
          TRUE ~ Term
        ),
        
        Class_Num = if_else(is.na(Class_Num), "1", Class_Num),
        
        Stars = case_when(
          P < 0.001 ~ "***",
          P < 0.01 ~ "**",
          P < 0.05 ~ "*",
          TRUE ~ ""
        ),
        Cell_Value = glue::glue("{Stars}{formatC(Estimate, format='f', digits=3)} ({formatC(SE, format='f', digits=3)})")
      )
    
    # Check for duplicates before pivoting (Safety check)
    if (any(duplicated(df_clean[, c("Type", "Term_Clean", "Class_Num")]))) {
      # Fallback: append original index to name if still duplicate
      df_clean <- df_clean %>% 
        group_by(Type, Term_Clean, Class_Num) %>%
        mutate(Term_Clean = if(n() > 1) paste0(Term_Clean, " (", row_number(), ")") else Term_Clean) %>%
        ungroup()
    }
    
    df_wide <- df_clean %>%
      select(Type, Term_Clean, Class_Num, Cell_Value) %>%
      pivot_wider(
        names_from = Class_Num,
        values_from = Cell_Value,
        names_prefix = "Class "
      ) %>%
      mutate(
        Type_Order = case_when(
          Type == "Class-Membership Model" ~ 1,
          Type == "Longitudinal Model" ~ 2,
          Type == "Random Effects & Residuals" ~ 3,
          TRUE ~ 4
        ),
        Term_Order = case_when(
          Term_Clean == "Intercept" ~ 1,
          Term_Clean == "Linear Slope" ~ 2,
          Term_Clean == "Quadratic Slope" ~ 3,
          TRUE ~ 4
        )
      ) %>%
      arrange(Type_Order, Term_Order) %>%
      select(-Type_Order, -Term_Order)
    
    df_wide %>%
      gt(groupname_col = "Type") %>%
      tab_header(title = glue("Model Estimates: {title_label}")) %>%
      cols_label(Term_Clean = "Parameter") %>%
      sub_missing(missing_text = "") %>%
      tab_source_note(
        source_note = md("*** p<0.001, ** p<0.01, * p<0.05. Values are Coefficient (SE).")
      ) %>%
      tab_source_note(
        source_note = glue("Interpretation: Class-Membership = Log-Odds; Longitudinal = {note}")
      ) %>%
      cols_align(
        align = "right",
        columns = starts_with("Class ")
      ) %>%
      cols_align(
        align = "left",
        columns = "Term_Clean"
      ) %>%
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_row_groups()
      ) %>%
      tab_options(
        table.width = pct(100),
        row.striping.include_table_body = FALSE
      )
  }
  
  # ============================================================================
  # PHASE 1: MODEL ACQUISITION (Load or Fit)
  # ============================================================================
  fitted_models <- list()
  temp_models   <- list() 
  
  # 1. Baseline LGCM
  message("\n[Phase 1] Checking Baseline LGCM...")
  base_path <- here(mod_dir, glue("baseline_lgcm_{target_var}_std.rds"))
  
  tryCatch({
    if (file.exists(base_path) && !rerun) {
      message("  -> Loading existing Baseline model...")
      m_baseline <- readRDS(base_path)
    } else {
      message("  -> Estimating new Baseline model...")
      if (use_gaussian) {
        m_baseline <- hlme(fixed = Y_OUTCOME ~ 1 + TIME, random = ~ 1 + TIME, subject = "ID", data = dt_long, ng = 1, verbose = FALSE)
      } else {
        m_baseline <- lcmm(fixed = Y_OUTCOME ~ 1 + TIME, random = ~ 1 + TIME, subject = "ID", data = dt_long, link = "thresholds", ng = 1, verbose = FALSE)
      }
      saveRDS(m_baseline, base_path)
    }
    fitted_models[["baseline"]] <- m_baseline
  }, error = function(e) message(glue("  -> Failed Baseline: {e$message}")))
  
  # 2. LCGA Models
  for (k in 1:num_classes) {
    message(glue("[Phase 1] Checking k={k} Model..."))
    model_path <- here(mod_dir, glue("lcga_k{k}_{target_var}_std.rds"))
    
    tryCatch({
      if (file.exists(model_path) && !rerun) {
        message("  -> Loading existing k={k} model...")
        m <- readRDS(model_path)
      } else {
        message("  -> Estimating new k={k} model...")
        if (use_gaussian) {
          if (k == 1) m <- hlme(fixed = Y_OUTCOME ~ 1 + TIME + I(TIME^2), random = ~ 1, subject = "ID", data = dt_long, ng = 1, verbose = FALSE)
          else m <- hlme(fixed = Y_OUTCOME ~ 1 + TIME + I(TIME^2), mixture = ~ 1 + TIME + I(TIME^2), random = ~ 1, subject = "ID", data = dt_long, ng = k, B = temp_models[[1]], verbose = FALSE)
        } else {
          if (k == 1) m <- lcmm(fixed = Y_OUTCOME ~ 1 + TIME + I(TIME^2), random = ~ 1, subject = "ID", data = dt_long, link = "thresholds", ng = 1, verbose = FALSE)
          else m <- lcmm(fixed = Y_OUTCOME ~ 1 + TIME + I(TIME^2), mixture = ~ 1 + TIME + I(TIME^2), random = ~ 1, subject = "ID", data = dt_long, link = "thresholds", ng = k, B = temp_models[[1]], verbose = FALSE)
        }
        saveRDS(m, model_path)
      }
      
      fitted_models[[paste0("k", k)]] <- m
      if (k == 1) temp_models[[1]] <- m
      
    }, error = function(e) message(glue("  -> Failed k={k}: {e$message}")))
  }
  
  # ============================================================================
  # PHASE 2: RESULTS GENERATION (Reporting)
  # ============================================================================
  message("\n[Phase 2] Generating Reports (Tables & Plots)...")
  
  lcga_results <- list(
    config = config,
    display_label = display_label,
    fit_summary_table = NULL,
    models = list()
  )
  
  # --- Report: Baseline LGCM ---
  if (!is.null(fitted_models[["baseline"]])) {
    m <- fitted_models[["baseline"]]
    
    set.seed(123)
    dt_subset <- dt_long %>% filter(ID %in% sample(unique(ID), min(50, length(unique(ID)))))
    
    dat_pred_base <- data.frame(TIME = seq(min(time_range), max(time_range), length.out = 50))
    pred_base <- predictY(m, dat_pred_base, var.time = "TIME")
    plot_data_mean <- pred_base$pred %>% as.data.frame() %>% mutate(TIME = dat_pred_base$TIME)
    pred_col <- setdiff(names(plot_data_mean), "TIME")[1]
    plot_data_mean <- plot_data_mean %>% rename(Value = all_of(pred_col))
    if (use_gaussian & var_type == "binary") plot_data_mean$Value <- pmin(pmax(plot_data_mean$Value, 0), 1)
    
    p <- ggplot() +
      geom_line(data = dt_subset, aes(x = TIME, y = Y_OUTCOME, group = ID), color = "grey70", alpha = 0.6, linewidth = 0.5, na.rm = TRUE) + # FIX: na.rm=TRUE
      geom_line(data = plot_data_mean, aes(x = TIME, y = Value), color = "blue", linewidth = 1.5) +
      labs(title = NULL, subtitle = "Baseline LGCM (1-Class, Random Slope)", y = y_label, x = time_label) +
      theme_minimal()
    
    if (use_gaussian & var_type == "binary") p <- p + scale_y_continuous(labels = scales::percent, limits = c(0, 1))
    
    ggsave(here(out_dir, glue("plot_baseline_{target_var}.png")), plot = p, width=8, height=6)
    
    lcga_results$models[["baseline_lgcm"]] <- list(
      model_obj = m,
      plot = p,
      estimates_table = extract_estimates_table(m, "Baseline LGCM (1-Class)", longitudinal_note)
    )
  }
  
  # --- Report: LCGA Models ---
  for (k in 1:num_classes) {
    m <- fitted_models[[paste0("k", k)]]
    if (is.null(m)) next
    
    k_res <- list(model_obj = m)
    
    dat_pred <- data.frame(TIME = seq(min(time_range), max(time_range), length.out = 50))
    pred_vals <- predictY(m, dat_pred, var.time = "TIME")
    plot_data <- pred_vals$pred %>% as.data.frame() %>% mutate(TIME = dat_pred$TIME)
    
    if (k == 1) {
      pred_col <- setdiff(names(plot_data), "TIME")[1]
      plot_data <- plot_data %>% rename(Value = all_of(pred_col)) %>% mutate(Class = "1")
    } else {
      plot_data <- plot_data %>% pivot_longer(cols = starts_with("Ypred_class"), names_to = "Class", values_to = "Value") %>% mutate(Class = str_replace(Class, "Ypred_class", ""))
    }
    
    if (use_gaussian & var_type == "binary") plot_data$Value <- pmin(pmax(plot_data$Value, 0), 1)
    
    p <- ggplot(plot_data, aes(x = TIME, y = Value, color = Class, group = Class)) +
      geom_line(linewidth = 1.2, na.rm = TRUE) + # FIX: na.rm=TRUE
      labs(title = NULL, subtitle = glue("{k} Classes"), y = y_label, x = time_label) +
      theme_minimal() + theme(legend.position = "bottom")
    
    if (use_gaussian & var_type == "binary") p <- p + scale_y_continuous(labels = scales::percent, limits = c(0, 1))
    k_res$plot <- p
    ggsave(here(out_dir, glue("plot_{target_var}_k{k}.png")), plot = p, width=8, height=6)
    
    # 2. Estimates Table (Stacked Wide) with Interpretation Note
    k_res$estimates_table <- extract_estimates_table(m, glue("{k}-Class Solution"), longitudinal_note)
    
    if (k > 1) {
      class_assignments <- m$pprob %>% rename(ID = ID) %>% select(ID, class, starts_with("prob"))
      avepp <- class_assignments %>% group_by(class) %>% summarise(Count = n(), AvePP = mean(get(paste0("prob", unique(class))))) %>% mutate(Proportion = Count / sum(Count))
      
      k_res$diagnostics_table <- avepp %>% 
        select(class, Count, Proportion, AvePP) %>%
        gt() %>% 
        tab_header(title = glue("{k}-Class Classification Diagnostics and Average Posterior Probabilities")) %>%
        cols_label(class = "Class", Count = "Count", Proportion = "Proportion", AvePP = "Avg. Posterior Prob.") %>%
        fmt_number(columns = c(AvePP), decimals = 3, use_seps=FALSE) %>% 
        fmt_percent(columns = c(Proportion), decimals = 1) %>% 
        tab_footnote(footnote = "AvePP > 0.70 indicates reliable classification.", locations = cells_column_labels(columns = c(AvePP))) %>%
        tab_style(style = list(cell_text(weight = "bold")), locations = cells_body(columns = c(AvePP), rows = AvePP > 0.7))
      
      k_res$adjudicated_data <- class_assignments %>%
        left_join(dt_long %>% select(ID, RAW_ID) %>% distinct(), by = "ID") %>%
        rename(!!sym(id_var) := RAW_ID) %>%
        select(!!sym(id_var), everything(), -ID)
    }
    
    lcga_results$models[[paste0("k", k)]] <- k_res
  }
  
  # --- D. Summary Table ---
  fit_data <- map_df(names(lcga_results$models), function(n) {
    m <- lcga_results$models[[n]]$model_obj
    if (is.null(m)) return(NULL)
    mod_label <- if (n == "baseline_lgcm") "Baseline LGCM (1-Class, Random Slope)" else paste0(m$ng, "-Class LCGA (Quad)")
    tibble(Model = mod_label, Classes = m$ng, AIC = m$AIC, BIC = m$BIC, Entropy = calc_entropy(m), Converged = (m$conv == 1))
  })
  
  if (nrow(fit_data) > 0) {
    fit_data <- fit_data %>% arrange(Classes, desc(Model))
    lcga_results$fit_summary_table <- fit_data %>% gt() %>% 
      tab_header(title = "LCGA Model Fit Comparison", subtitle = display_label) %>%
      tab_source_note(source_note = glue("Method: {method_label}")) %>%
      fmt_number(columns = c(AIC, BIC), decimals = 2, use_seps = FALSE) %>% 
      fmt_number(columns = c(Entropy), decimals = 3, use_seps = FALSE) %>% 
      sub_missing(missing_text = "-") %>% 
      tab_footnote(footnote = "Lower values indicate better fit.", locations = cells_column_labels(columns = c(AIC, BIC))) %>%
      tab_footnote(footnote = "Values > 0.80 indicate good separation (N/A for 1 class).", locations = cells_column_labels(columns = c(Entropy))) %>%
      tab_style(style = list(cell_fill(color = "#e6f3ff"), cell_text(weight = "bold")), locations = cells_body(rows = BIC == min(BIC)))
  }
  
  saveRDS(lcga_results, here(out_dir, glue("lcga_results_{target_var}.rds")))
  message(glue("Analysis complete for {target_var}."))
  return(lcga_results)
}

# 4. CONFIGURATION & EXECUTION -------------------------------------------------

# --- A. Global Settings (Directories & Data) ---
output_dir <- here("Trajectory_Output_lcmm")
if (!dir.exists(output_dir)) dir.create(output_dir)

models_dir <- here(output_dir, "model_objects")
if (!dir.exists(models_dir)) dir.create(models_dir)

data_path <- here("../../_private_use/data_management/data/co-luminate/adam/dt_childhood_exposure.RData")
if (!file.exists(data_path)) stop("Data file not found.")
load(data_path)
message(glue("Data loaded successfully from: {data_path}"))

# --- B. Analysis Configurations ---
analysis_configs <- list(
  # Config 1: Lived Without Mother (Binary)
  list(
    target_var  = "LIVEDWOM",
    var_label   = "Lived Without Mother",
    id_var      = "USUBJID",
    time_var    = "EXPAGE",
    time_var_label = "Child Age (Years)",
    var_type    = "binary",
    use_linear_for_binary = TRUE,
    time_range  = 1:5,
    num_classes = 4,
    rerun_models = FALSE
  ),
  # Config 2: Lived Without Father (Binary)
  list(
    target_var  = "LIVEDWOF",
    var_label   = "Lived Without Father",
    id_var      = "USUBJID",
    time_var    = "EXPAGE",
    time_var_label = "Child Age (Years)",
    var_type    = "binary",
    use_linear_for_binary = TRUE,
    time_range  = 1:5,
    num_classes = 4,
    rerun_models = FALSE
  ),
  # Config 3: Household Socioeconomic Status (Ordinal)
  list(
    target_var  = "HHSES",
    var_label   = "Household Socioeconomic Status",
    id_var      = "USUBJID",
    time_var    = "EXPAGE",
    time_var_label = "Child Age (Years)",
    var_type    = "ordinal",
    use_linear_for_binary = FALSE, 
    time_range  = 1:5,
    num_classes = 4,
    rerun_models = FALSE
  ),
  # Config 4: Mother’s Level of Education (Ordinal)
  list(
    target_var  = "EDUCMOTH",
    var_label   = "Mother’s Level of Education",
    id_var      = "USUBJID",
    time_var    = "EXPAGE",
    time_var_label = "Child Age (Years)",
    var_type    = "ordinal",
    use_linear_for_binary = FALSE, 
    time_range  = 1:5,
    num_classes = 4,
    rerun_models = FALSE
  ),
  # Config 5: Father’s Level of Education (Ordinal)
  list(
    target_var  = "EDUCFATH",
    var_label   = "Father’s Level of Education",
    id_var      = "USUBJID",
    time_var    = "EXPAGE",
    time_var_label = "Child Age (Years)",
    var_type    = "ordinal",
    use_linear_for_binary = FALSE, 
    time_range  = 1:5,
    num_classes = 4,
    rerun_models = FALSE
  ),
  # Config 6: Number of Children Aged 0–4y in Household
  list(
    target_var  = "HHCHLD4",
    var_label   = "Number of Children Aged 0–4y in Household",
    id_var      = "USUBJID",
    time_var    = "EXPAGE",
    time_var_label = "Child Age (Years)",
    var_type    = "continuous",
    use_linear_for_binary = FALSE, 
    time_range  = 1:5,
    num_classes = 4,
    rerun_models = FALSE
  ),
  # Config 7: Number of Children Aged 0–14y in Household
  list(
    target_var  = "HHCHLD14",
    var_label   = "Number of Children Aged 0–14y in Household",
    id_var      = "USUBJID",
    time_var    = "EXPAGE",
    time_var_label = "Child Age (Years)",
    var_type    = "continuous",
    use_linear_for_binary = FALSE, 
    time_range  = 1:5,
    num_classes = 4,
    rerun_models = FALSE
  )
)

# --- C. Run Batch Analysis ---
batch_results <- map(analysis_configs, ~ run_single_lcga(.x, dt_childhood_exposure, output_dir, models_dir))
names(batch_results) <- map_chr(analysis_configs, "target_var")