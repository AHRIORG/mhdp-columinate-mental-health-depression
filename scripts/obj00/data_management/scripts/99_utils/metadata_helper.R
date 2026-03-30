# ==============================================================================
# METADATA HELPER FUNCTION
# Purpose: centralized storage for variable renaming and labeling metadata
# ==============================================================================

get_variable_metadata <- function() {
  list(
    # --- Renaming Mappings (New Name = Old Name) ---
    rename_childhood = c(
      "USUBJID" = "IIntId",
      "HHID" = "HouseholdId",
      "SEX" = "sex",
      "BRTHDT" = "dateofbirth",
      "BRTHYR" = "year_of_birth",
      "HDSSYR" = "hdssyr",
      "MINEXPYR" = "minexpage",
      "MAXEXPYR" = "maxexpage",
      "EXPAGE" = "expage",
      "BRTHEXPM" = "birth_expmonth",
      "LIVEDWOM" = "p_without_her",
      "LIVEDWOF" = "p_without_him",
      "MEDU" = "p_edu_mother",
      "FEDU" = "p_edu_father",
      "HHCHLD14" = "hh_child_14num",
      "HHCHLD4" = "hh_child_4num",
      "HHSES" = "hh_ses"
    ),
    
    rename_multistudies = c(
      "USUBJID" = "IIntId",
      "HHID" = "HouseholdId",
      "BRTHDT" = "dateofbirth",
      "VISITDT" = "visitdate",
      "SEX" = "sex",
      "AGE" = "age",
      "URBANCAT" = "urban_or_rural",
      "MGRC" = "migration",
      "ORPH" = "orphanstatus",
      "FDEC" = "fatherstatus",
      "MDEC" = "motherstatus",
      "CLINDATA" = "clinic_data",
      "VISIT" = "round",
      "STUDY" = "study_name",
      "SCHL" = "education",
      "EMPLOY" = "employment",
      "GOVG" = "governmentgrant",
      "SCSP" = "social_support",
      "ESXC" = "everhadsex",
      "HIVS" = "hiv_status",
      "VLNC" = "violence",
      "FDSC" = "foodsec",
      "EVRPREGCAT" = "everpregnant",
      "SEXBHV" = "sexual_behavior",
      "CDMLESSFL" = "condomlesssex",
      "DNKA" = "everdrankalcohol",
      "PHQSCR" = "depression_score",
      "PHQBIN" = "depression_outcome",
      "SCISCR" = "suicidal_ideation_score",
      "SCIBIN" = "suicidal_ideation_outcome",
      "SSQSCR" = "ssq14_score",
      "SSQBIN" = "ssq14_outcome",
      "SCISRC" = "suicidal_ideation_score_source"
    ),
    
    rename_isisekelo = c(
      "USUBJID" = "IIntId",
      "HHID" = "HHIntId",
      "ARMCD" = "arm",
      "BRTHDT" = "pd_dob",
      "AGE" = "pd_age",
      "VISITDT" = "vstdt",
      "STUDYDY" = "stddy",
      "SEX" = "pd_sex",
      "SSQ01" = "mv_ment_health01",
      "SSQ02" = "mv_ment_health02",
      "SSQ03" = "mv_ment_health03",
      "SSQ04" = "mv_ment_health04",
      "SSQ05" = "mv_ment_health05",
      "SSQ06" = "mv_ment_health06",
      "SSQ07" = "mv_ment_health07",
      "SSQ08" = "mv_ment_health08",
      "SSQ09" = "mv_ment_health09",
      "SSQ10" = "mv_ment_health10",
      "SSQ11" = "mv_ment_health11",
      "SSQ12" = "mv_ment_health12",
      "SSQ13" = "mv_ment_health13",
      "SSQ14" = "mv_ment_health14",
      "PHQ901" = "mv_phq9_ment01",
      "PHQ902" = "mv_phq9_ment02",
      "PHQ903" = "mv_phq9_ment03",
      "PHQ904" = "mv_phq9_ment04",
      "PHQ905" = "mv_phq9_ment05",
      "PHQ906" = "mv_phq9_ment06",
      "PHQ907" = "mv_phq9_ment07",
      "PHQ908" = "mv_phq9_ment08",
      "PHQ909" = "mv_phq9_ment09"
    ),
    
    # --- Variable Labels ---
    # Using a named list ensures we can easily check for duplicates by name
    labels = list(
      "USUBJID" = "Unique Subject Identifier",
      "HHID" = "Unique Household Identifier",
      "CLYR" = "Calendar year",
      "BRTHYR" = "Year of Birth",
      "HDSSYR" = "HDSS Data Collection Year",
      "MINEXPYR" = "Minimum Exposure Age",
      "MAXEXPYR" = "Maximum Exposure Age",
      "EXPAGE" = "Exposure Age",
      "BRTHEXPM" = "Number of Exposure Months During Year of Birth",
      "LIVEDWOM" = "Lived Without Mother",
      "LIVEDWOF" = "Lived Without Father",
      "MEDU" = "Mother's Education",
      "FEDU" = "Father's Education",
      "HHCHLD14" = "Number of children (0–14y)",
      "HHCHLD14O" = "Overcrowded (0–14y children > 3)",
      "HHCHLD4" = "Number of children (0–4y)",
      "HHCHLD4O" = "Overcrowded (0–4y children > 3)",
      "HHSES" = "Household Socioeconomic Status",
      "BRTHDT" = "Date of Birth",
      "VISITDT" = "Visit Date",
      "SEX" = "Sex",
      "AGE" = "Age (Years)",
      "URBANCAT" = "Household Urbanicity Category",
      "URBAN"="Urban",
      "PURBN"="Peri-Urban",
      "RURAL"="Rural",
      "MGRC" = "Migration Status",
      "PIPSA"="PIPSA",
      "EXTMG"="External migration",
      "FDEC" = "Father Deceased",
      "MDEC" = "Mother Deceased",
      "ORPH" = "Orphanhood Status",
      "CLINDATA" = "Has Clinic Data",
      "VISIT" = "Visit Name",
      "STUDY" = "Study Name",
      "SCHL" = "Currently in School",
      "EMPLOY" = "Current Employment",
      "GOVG" = "Receives Government Grant",
      "ESXC" = "Ever Had Sex",
      "HIVS" = "HIV Test Result",
      "VLNC" = "Ever Experienced Violence",
      "FDSC" = "Food Insecurity",
      "EVRPREGCAT" = "Ever Been Pregnant",
      "SEXBHV" = "Sexual Behavior Category",
      "DNKA" = "Ever Drank Alcohol",
      "SCSP" = "Social Support",
      "CDMLESSFL" = "Condomless Sex",
      "HSV2STRESC" = "HSV-2 Test Result",
      "CHLSTRESC" = "Chlamydia Result",
      "GONSTRESC" = "Gonorrhea Result",
      "PHQSCR" = "PHQ-9 Total Score",
      "PHQBIN" = "PHQ-9 ≥ 10",
      "SCISCR" = "Suicidal Ideation Score",
      "SCIBIN" = "Suicidal Ideation",
      "SSQSCR" = "SSQ-14 Score",
      "SSQBIN" = "SSQ-14 ≥ 9",
      "SCISRC" = "Suicidal Score Source",
      "TRISTRESC" = "Trichomonas Result",
      "ARMCD" = "Study Arm",
      "STUDYDY" = "Study Day",
      "SSQ01" = "SSQ01-Thinking Deeply or About Many Things",
      "SSQ02" = "SSQ02-Trouble Concentrating",
      "SSQ03" = "SSQ03-Lose Temper or Get Annoyed",
      "SSQ04" = "SSQ04-Nightmares or Bad Dreams",
      "SSQ05" = "SSQ05-See or Hear Things Others Could Not",
      "SSQ06" = "SSQ06-Stomach Aching",
      "SSQ07" = "SSQ07-Frightened by Trivial Things",
      "SSQ08" = "SSQ08-Trouble Sleeping or Lose Sleep",
      "SSQ09" = "SSQ09-Cried or Wanted to Cry",
      "SSQ10" = "SSQ10-Feeling Run Down or Tired",
      "SSQ11" = "SSQ11-At Times Feel Like Committing Suicide",
      "SSQ12" = "SSQ12-Generally Unhappy",
      "SSQ13" = "SSQ13-Work or School Lagging Behind",
      "SSQ14" = "SSQ14-Problems in Deciding What to Do",
      "PHQ901" = "PHQ901-Little Interest or Pleasure",
      "PHQ902" = "PHQ902-Feeling Down, Depressed, Hopeless",
      "PHQ903" = "PHQ903-Trouble Sleeping",
      "PHQ904" = "PHQ904-Feeling Tired or Having Little Energy",
      "PHQ905" = "PHQ905-Poor Appetite or Overeating",
      "PHQ906" = "PHQ906-Feeling Bad About Yourself",
      "PHQ907" = "PHQ907-Trouble Concentrating",
      "PHQ908" = "PHQ908-Moving or Speaking Slowly/Fidgety",
      "PHQ909" = "PHQ909-Thoughts of Self-Harm"
    )
  )
}

#' Compare and Extract Missing Metadata
#'
#' This function inspects a labelled dataframe, identifies variable labels that
#' are present in the data but missing from the central metadata 'labels' list,
#' and returns them. This helps in keeping the metadata helper updated.
#'
#' @param df A dataframe containing labelled variables.
#' @param current_metadata The list object returned by get_variable_metadata().
#'
#' @return A named list of new labels found in 'df' that are not in 'current_metadata'.
compare_and_extract_metadata <- function(df, current_metadata) {
  # 1. Extract labels from the dataframe
  df_labels <- var_label(df)
  
  # Remove NULLs (unlabelled variables)
  df_labels <- df_labels[!sapply(df_labels, is.null)]
  
  if (length(df_labels) == 0) {
    message("No labels found in the provided dataframe.")
    return(invisible(list()))
  }
  
  # 2. Get existing label keys from metadata
  existing_keys <- names(current_metadata$labels)
  
  # 3. Identify new keys
  new_keys <- setdiff(names(df_labels), existing_keys)
  
  if (length(new_keys) == 0) {
    message("All labels in dataframe are already present in metadata.")
    return(invisible(list()))
  }
  
  # 4. Extract the new labels
  new_labels <- df_labels[new_keys]
  
  # 5. Output guidance
  message(sprintf("Found %d new variable label(s) not in metadata.", length(new_labels)))
  message("You can copy the list below to update get_variable_metadata()$labels:")
  
  # Return as a named list (can be printed easily)
  return(new_labels)
}