# Project path package
require(here)
require(fs)

# Functions and Packages
## Functions
{
  parent_dirs<-dir_ls(path_home())
  func_cloud_dirs<-file.path(
    parent_dirs[grepl("Library",parent_dirs)],
    "Mobile Documents/com~apple~CloudDocs/WCLOUD/pipelines_and_functions"
  )
  source(file.path(func_cloud_dirs,"data_processing_and_description_funcs.R"))
}
# Packages

pkgs<-c("tidyverse","labelled","haven")
func_install_and_load(pkgs)

# Read data
## Raw
dt_.raw <- read_dta("data_management/raw_data/columinate_clean.dta") |> 
  as_factor()

# Dictionary: Data variable and description
dictionary<-func_dictionary(dt_.raw)

# Episodes

vars<-names(dt_.raw)
common_vars<-vars[which(is.na(as.numeric(str_extract(vars, "[\\d\\.]+$"))))]
episodes<-sort(as.numeric(na.omit(unique(as.numeric(str_extract(vars, "[\\d\\.]+$"))))))
for(episode. in episodes){
  episode_vars<-vars[which(as.numeric(str_extract(vars, "[\\d\\.]+$"))==episode.)]
  dt.<-dt_.raw |> 
    mutate(episode=episode.) |> 
    select(any_of(common_vars),episode,any_of(episode_vars)) |> 
    rename_with(~ gsub(paste0(episode.,"$"),"",.x))
  
  if(episode.==min(episodes) | episode.==min(episodes[which(episodes>10)])){
    if(episode.<=10){
      DATA1<-dt. |> mutate(visityr=year(visitdate))
    } else{
      DATA2<-dt. |> rename(visityr=episode)
    }
  } else{
    if(episode.<=10){
      DATA1<-bind_rows(DATA1,dt.) |> 
        arrange(IndividualId,IIntId,visitdate) |> 
        filter(!if_all(names(DATA1)[which(!names(DATA1)%in%c(common_vars, "episode"))], is.na))
    } else{
      DATA2<-bind_rows(DATA2,dt. |> rename(visityr=episode)) |> 
        arrange(IndividualId,IIntId,visityr)# |> 
      #filter(!if_all(names(DATA2)[which(!names(DATA2)%in%c(common_vars, "visityr"))], is.na))
    }
    
  }
}

varlabs<-dictionary |> 
  group_by(row_number()) |> 
  mutate(
    rmv=paste0(as.numeric(str_extract(variable, "[\\d\\.]+$")),"$"),
    variable=gsub(rmv,"",variable),
    rmv=paste0("^",as.numeric(str_extract(label, "^[\\d\\.]+"))),
    label=trimws(gsub(rmv,"",label)),
    label=ifelse(variable=="round","Visit name",label)
  ) |> 
  ungroup() |> 
  select(variable,label) |> 
  distinct()


dt_analysis_multistudies<- left_join(DATA1, DATA2) |>
  # Arrange by individual and time to ensure correct filling order
  arrange(IndividualId, IIntId, visitdate) |>
  # Operations are performed within each individual's timeline
  group_by(IndividualId, IIntId) |>
  mutate(
    neg_fill = if_else(hiv_status == "Negative", "Negative", NA_character_),
    pos_fill = if_else(hiv_status == "Positive", "Positive", NA_character_)
  ) |>
  fill(neg_fill, .direction = "up") |>
  fill(pos_fill, .direction = "down") |>
  ungroup() |>
  mutate(
    hiv_status = factor(coalesce(neg_fill, pos_fill, as.character(hiv_status)))
  ) |>
  select(-c(episode, visityr, neg_fill, pos_fill)) |>
  set_variable_labels(
    .labels = setNames(as.list(varlabs$label), varlabs$variable), .strict = FALSE
  )
dt_dictionary_multistudies<-func_dictionary(dt_analysis_multistudies)

save(dt_analysis_multistudies,file=here("data_management/transformed_data/dt_analysis_multistudies.RData"))
save(dt_dictionary_multistudies,file=here("data_management/transformed_data/dt_dictionary_multistudies.RData"))

# Isisekelo Screening Data
isisekelo_semp <- read_dta("data_management/raw_data/isisekelo semp.dta") |> 
  as_factor()

dictionary<-func_dictionary(isisekelo_semp)

selected_vars<-dictionary |> 
  filter(
    grepl("^arm$|mv_ment|mv_phq9|IIntId|HHIntId|pd_|visitdate|foodfreq",variable)
  ) |> 
  pull(variable)

dt_analysis_screening.<-isisekelo_semp |> 
  select(
    any_of(selected_vars)
  ) |> 
  mutate(pd_age=round(as.numeric(difftime(ca_visitdate,pd_dob,units="days"))/365.5,1),
         stddy=as.numeric(difftime(ca_visitdate,min(ca_visitdate,na.rm=TRUE),units="days")),
         brtdt=pd_dob) |> 
  rename(vstdt=ca_visitdate) |> 
  relocate(brtdt,vstdt,stddy,.after=pd_age)



# --- Rename Variables to CDISC Standard ---
var_labels <- list(
  USUBJID = "Unique Subject Identifier",
  HHID = "Household Identifier",
  ARMCD = "Study Arm",
  BRTHDT = "Date of Birth",
  AGE = "Age",
  VISITDT = "Visit Date",
  STUDYDY = "Study Day",
  SEX = "Sex",
  
  SSQ01 = "Thinking deeply or about many things (SSQ)",
  SSQ02 = "Trouble concentrating (SSQ)",
  SSQ03 = "Lose temper or get annoyed (SSQ)",
  SSQ04 = "Nightmares or bad dreams (SSQ)",
  SSQ05 = "See or hear things others could not (SSQ)",
  SSQ06 = "Stomach aching (SSQ)",
  SSQ07 = "Frightened by trivial things (SSQ)",
  SSQ08 = "Trouble sleeping or lose sleep (SSQ)",
  SSQ09 = "Cried or wanted to cry (SSQ)",
  SSQ10 = "Feeling run down or tired (SSQ)",
  SSQ11 = "At times feel like committing suicide (SSQ)",
  SSQ12 = "Generally unhappy (SSQ)",
  SSQ13 = "Work or school lagging behind (SSQ)",
  SSQ14 = "Problems in deciding what to do (SSQ)",
  
  PHQ901 = "Little interest or pleasure (PHQ)",
  PHQ902 = "Feeling down, depressed, hopeless (PHQ)",
  PHQ903 = "Trouble sleeping (PHQ)",
  PHQ904 = "Feeling tired or having little energy (PHQ)",
  PHQ905 = "Poor appetite or overeating (PHQ)",
  PHQ906 = "Feeling bad about yourself (PHQ)",
  PHQ907 = "Trouble concentrating (PHQ)",
  PHQ908 = "Moving or speaking slowly/fidgety (PHQ)",
  PHQ909 = "Thoughts of self-harm (PHQ)",
  
  FI01 = "FI: Household cut food size in last 12 mos",
  FI02 = "FI: Frequency of cutting food size"
)

dt_analysis_screening <- dt_analysis_screening. |>
  rename(
    # Demographics and Identifiers
    USUBJID = IIntId,
    HHID = HHIntId,
    ARMCD = arm,
    # Assuming pd_dob and brtdt are the same; choosing one.
    BRTHDT = pd_dob, 
    AGE = pd_age,
    VISITDT = vstdt,
    STUDYDY = stddy,
    SEX = pd_sex,
    
    # General Mental Health Screen (SSQ)
    SSQ01 = mv_ment_health01,
    SSQ02 = mv_ment_health02,
    SSQ03 = mv_ment_health03,
    SSQ04 = mv_ment_health04,
    SSQ05 = mv_ment_health05,
    SSQ06 = mv_ment_health06,
    SSQ07 = mv_ment_health07,
    SSQ08 = mv_ment_health08,
    SSQ09 = mv_ment_health09,
    SSQ10 = mv_ment_health10,
    SSQ11 = mv_ment_health11,
    SSQ12 = mv_ment_health12,
    SSQ13 = mv_ment_health13,
    SSQ14 = mv_ment_health14,
    
    # Patient Health Questionnaire-9 (PHQ-9)
    PHQ901 = mv_phq9_ment01,
    PHQ902 = mv_phq9_ment02,
    PHQ903 = mv_phq9_ment03,
    PHQ904 = mv_phq9_ment04,
    PHQ905 = mv_phq9_ment05,
    PHQ906 = mv_phq9_ment06,
    PHQ907 = mv_phq9_ment07,
    PHQ908 = mv_phq9_ment08,
    PHQ909 = mv_phq9_ment09,
    
    # Food Insecurity (FI)
    FI01 = ie_foodfreq,
    FI02 = ie_foodfreqskip
  ) |> 
  # Select variables to keep, removing redundant ones like 'brtdt'
  select(-brtdt) |>
  mutate(
    ARMCD=paste("Arm",ARMCD), 
    across(starts_with("PHQ9"), 
           ~ factor(.x,
                    levels=c("Absolutely not",
                             "Several days",
                             "More than half days",
                             "Almost daily")
           )
    ),
    FI01=factor(as.character(FI01),levels=unique(as.character(FI01))),
    FI02=ifelse(grepl("No|Prefer not to answer",FI01),"Not Applicable",as.character(FI02)),
    FI02=factor(as.character(FI02),levels=unique(as.character(FI02))),
    across(starts_with("SSQ"), 
           ~ factor(as.character(.x),
                    levels=unique(as.character(.x))
                    )
           )
    ) |> 
  set_variable_labels(
    .labels = var_labels,.strict=FALSE
    )



dt_dictionary_screening<-func_dictionary(dt_analysis_screening)

save(dt_analysis_screening,file=here("data_management/transformed_data/dt_analysis_screening.RData"))
save(dt_dictionary_screening,file=here("data_management/transformed_data/dt_dictionary_screening.RData"))




