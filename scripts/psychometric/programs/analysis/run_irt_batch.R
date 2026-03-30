# ============================================================
# FULL SCENARIO BATCH RUNNER
# ============================================================

# 0) Source Latent Harmonization via Item Response Theory (IRT) Scrtpt
source(here("statistical_analysis/programs/utils/irt/full_itr.R"))

# 1) Output folder (as requested)
output_dir <- here::here("statistical_analysis/output/objects/irt")

# 1a) Data preflight
psychometric_data_path <- Sys.getenv(
  "COLU_PSYCHOMETRIC_RDATA",
  unset = here::here("data/inputs/dt_psychometric.RData")
)
if (!file.exists(psychometric_data_path) && !exists("dt_psychometric", inherits = TRUE)) {
  stop(
    "Data file not found at: ", psychometric_data_path,
    ". Set COLU_PSYCHOMETRIC_RDATA or place dt_psychometric.RData under data/inputs/."
  )
}

# 2) Parallel workers: leave at least 2 cores free for the system
cores <- parallel::detectCores()
if (is.na(cores) || cores < 2) cores <- 2L
workers <- max(1L, min(8L, cores - 2L))   # change 8 -> 10 if you want more parallelism
workers

# Optional: reduce "future" max globals issues in some environments
options(future.globals.maxSize = 8 * 1024^3)  # 8GB

# 3) Define the FULL scenario grid inputs
phq_models <- list(
  PHQ_1F = 1,
  PHQ_2F_validated = phq_structure_validated
)

ssq_models <- list(
  SSQ_1F = 1,
  SSQ_2F_validated = ssq_structure_validated
)

joint_models <- list(
  JOINT_1F = 1,
  JOINT_SPLIT_2F = paste0(
    "F1 = PHQ901, PHQ902, PHQ903, PHQ904, PHQ905, PHQ906, PHQ907, PHQ908, PHQ909\n",
    "F2 = SSQ01, SSQ02, SSQ03, SSQ08, SSQ09, SSQ10, SSQ11, SSQ12, SSQ13, SSQ14\n",
    "COV = F1*F2"
  ),
  JOINT_4F_validated = make_joint_validated_4f_structure()
)

cutoff_methods <- c("model_tcc", "youden", "sens_at_spec", "spec_at_sens")
target_specificities <- c(0.90, 0.95)
target_sensitivities <- c(0.95, 0.98)

apply_modes <- c("single", "hierarchical")
apply_groups <- c("none", "SEX", "AGEGRP", "SEX_AGEGRP")

# Evaluate all groupings for transparency
eval_groupings <- c("none", "SEX", "AGEGRP", "SEX_AGEGRP")

apply_hierarchy <- c("SEX_AGEGRP", "SEX", "AGEGRP", "none")

# 4) How many scenarios is this?
n_scenarios <-
  length(phq_models) * length(ssq_models) * length(joint_models) *
  length(cutoff_methods) * length(target_specificities) * length(target_sensitivities) *
  length(apply_modes) * length(apply_groups)

message("Total scenarios to run: ", n_scenarios)
# Expected with these defaults: 2*2*3*4*3*3*2*4 = 3456 scenarios

# 5) Run the full batch (RESUME-SAFE)
# - overwrite=FALSE means: keep existing scenarios; skip those already saved
batch_res <- batch_run_irt(
  phq_models = phq_models,
  ssq_models = ssq_models,
  joint_models = joint_models,
  
  cutoff_methods = cutoff_methods,
  target_specificities = target_specificities,
  target_sensitivities = target_sensitivities,
  
  apply_modes = apply_modes,
  apply_groups = apply_groups,
  eval_groupings = eval_groupings,
  apply_hierarchy = apply_hierarchy,
  
  method_phq = "auto",
  method_ssq = "auto",
  method_joint = "auto",
  method_mg = "auto",
  
  run_sensitivity = FALSE,       # keep FALSE during batch
  overwrite = FALSE,             # IMPORTANT: resume-safe
  workers = workers,
  
  output_dir = output_dir,
  save_prefix = "IRT_batch"
)
saveRDS(batch_res,file = file.path(output_dir,"batch_res.rds"))


