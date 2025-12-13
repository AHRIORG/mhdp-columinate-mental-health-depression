# ==============================================================================
# PROJECT: CO-LUMINATE - Multi-Trajectory Profile Analysis
# OBJECTIVE: Cluster individuals based on their univariate trajectory class memberships
# METHOD: Latent Class Analysis (LCA) using poLCA
# INPUT: Adjudicated CSVs from 'Trajectory_Output_lcmm'
# OUTPUT: Profile plots and summary tables
# ==============================================================================

# 1. SETUP ---------------------------------------------------------------------
# Check and install poLCA if missing (specific to this analysis)
if (!requireNamespace("poLCA", quietly = TRUE)) install.packages("poLCA")

library(tidyverse)
library(poLCA)
library(here)
library(gt)
library(glue)

# Directory containing univariate results
input_dir <- here("Trajectory_Output_lcmm")
output_dir <- here("Trajectory_Output_Profiles")
if (!dir.exists(output_dir)) dir.create(output_dir)

# 2. LOAD & MERGE ADJUDICATED DATA ---------------------------------------------
message("Loading adjudicated class assignments from univariate models...")

# Find all adjudicated csv files
file_list <- list.files(input_dir, pattern = "adjudicated_.*\\.csv", full.names = TRUE)

if (length(file_list) == 0) stop("No adjudicated data files found. Run the univariate script first.")

# Read and Merge
# We assume standard naming: adjudicated_VARNAME.csv
merged_data <- file_list %>%
  map(function(f) {
    # Extract variable name from filename (e.g., adjudicated_LIVEDWOM.csv -> LIVEDWOM)
    var_name <- str_extract(basename(f), "(?<=adjudicated_).*(?=\\.csv)")
    
    read_csv(f, show_col_types = FALSE) %>%
      select(USUBJID, LCGA_Class) %>% # Keep only ID and Class
      rename(!!var_name := LCGA_Class) # Rename Class col to Var Name
  }) %>%
  reduce(full_join, by = "USUBJID")

# Clean Data for LCA
# poLCA requires positive integers (1, 2, 3...), no NAs allowed for included variables
lca_data <- merged_data %>%
  drop_na() %>% # Listwise deletion for screening
  mutate(across(-USUBJID, ~ as.numeric(as.factor(.)))) # Ensure 1-based integers

message(glue("Data merged. N = {nrow(lca_data)} individuals with complete trajectory data."))
print(head(lca_data))

# 3. DEFINE LCA FORMULA --------------------------------------------------------
# Use all columns except USUBJID as indicators
indicators <- setdiff(names(lca_data), "USUBJID")
f_lca <- as.formula(paste("cbind(", paste(indicators, collapse = ","), ") ~ 1"))

message(glue("LCA Formula: {paste(deparse(f_lca), collapse='')}"))

# 4. RUN LCA MODELS (1 to 5 Profiles) ------------------------------------------
set.seed(2025)
lca_results <- list()
fit_stats <- data.frame()

for (k in 1:5) {
  message(glue("Estimating {k}-Profile Solution..."))
  
  tryCatch({
    # Run poLCA (nrep=5 to avoid local maxima)
    mod <- poLCA(f_lca, lca_data, nclass = k, nrep = 5, verbose = FALSE, graphs = FALSE)
    
    # Store Results
    lca_results[[paste0("k", k)]] <- mod
    
    # Calculate Entropy (poLCA doesn't output it directly in simple summary)
    # Entropy formula: 1 - sum(-p*ln(p)) / (N*ln(K))
    post_probs <- mod$posterior
    numerator <- sum(-rowSums(post_probs * log(post_probs + 1e-9)))
    entropy <- 1 - (numerator / (nrow(post_probs) * log(k)))
    if(k==1) entropy <- NA
    
    fit_stats <- rbind(fit_stats, data.frame(
      Profiles = k,
      AIC = mod$aic,
      BIC = mod$bic,
      Likelihood = mod$llik,
      Entropy = entropy
    ))
    
  }, error = function(e) message(glue("Error k={k}: {e$message}")))
}

# 5. FIT SUMMARY ---------------------------------------------------------------
tbl_fit <- fit_stats %>%
  gt() %>%
  tab_header(title = "Multi-Trajectory Profile Fit Comparison") %>%
  fmt_number(columns = c(AIC, BIC, Likelihood), decimals = 2) %>%
  fmt_number(columns = c(Entropy), decimals = 3) %>%
  sub_missing(missing_text = "-") %>%
  tab_style(
    style = list(cell_fill(color = "#e6f3ff"), cell_text(weight = "bold")),
    locations = cells_body(rows = BIC == min(BIC))
  )

# Save
saveRDS(list(fit = fit_stats, models = lca_results), here(output_dir, "lca_profile_results.rds"))

# 6. VISUALIZE BEST SOLUTION ---------------------------------------------------
best_k <- fit_stats$Profiles[which.min(fit_stats$BIC)]
best_mod <- lca_results[[paste0("k", best_k)]]

message(glue("\nBest fit based on BIC: {best_k} Profiles. Generating visualization..."))

# Extract parameters for plotting
# Reshape probabilities: Profile -> Variable -> Class
probs_long <- best_mod$probs %>%
  map2_dfr(names(best_mod$probs), ~ {
    as.data.frame(.x) %>%
      mutate(Profile = row_number()) %>%
      pivot_longer(cols = -Profile, names_to = "Class_Label", values_to = "Probability") %>%
      mutate(Variable = .y)
  }) %>%
  mutate(
    # Clean up Class Labels (usually "Pr(1)", "Pr(2)")
    Class = parse_number(Class_Label)
  )

# Plot: Heatmap of Class Probabilities per Profile
p_profiles <- ggplot(probs_long, aes(x = as.factor(Class), y = Probability, fill = as.factor(Class))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Profile ~ Variable, labeller = label_both) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = glue("Multi-Trajectory Profiles ({best_k} Groups)"),
    subtitle = "Probability of belonging to univariate trajectory classes within each Profile",
    x = "Univariate Trajectory Class (1=Low, 2=High, etc.)",
    y = "Probability",
    fill = "Original Class"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(here(output_dir, glue("plot_profiles_k{best_k}.png")), plot = p_profiles, width = 12, height = 8)

message(glue("DONE. Results saved to {output_dir}"))