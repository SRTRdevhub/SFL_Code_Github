
# Cleans the SAF for use in the survival from listing model


LungDataCleaning <- function(root, 
                             organ, 
                             organ_ab, 
                             current_saf, 
                             cohort_start, 
                             cohort_end, 
                             max_followup){
  
  
  ### Reading the candidate data
  cand <- read.csv(file.path(root,
                             "Data Sets",
                             paste0(organ_ab, 
                                    "_candidate_file_", 
                                    current_saf,
                                    ".csv")),
                   stringsAsFactors = FALSE,
                   header = TRUE)

  
  ### We only want to evaluate candidates with the original listing. This matters in lung because OPTN intentionally double lists certain pediatric 
  ###   candidates to ensure they have adult priority. For this reason, I am removing rows that point to an original PX_ID
  cand <- cand %>%
    subset(is.na(CAND_PL_ORIG_PX_ID))

  
  ### Based on the cohort definition, these candidates cannot be in the cohort
  cand <- cand %>%
    subset(as.Date(CAN_LISTING_DT, format = "%m/%d/%y") >= 
             (as.Date(cohort_start, format = "%Y%m%d") - (365 * max_followup))) %>%
    subset(WL_ORG == "LU: Lung") %>%
    mutate(can_male = ifelse(CAN_GENDER == "", NA,
                             ifelse(CAN_GENDER == "M", 1, 0)),
           can_male_missing = ifelse(CAN_GENDER == "", 1, 0),
           can_age = CAN_AGE_IN_MONTHS_AT_LISTING / 12,
           can_age_missing = ifelse(is.na(CAN_AGE_IN_MONTHS_AT_LISTING), 1, 0),
           can_race = ifelse(CAN_RACE_SRTR == "", NA,
                             ifelse(CAN_RACE_SRTR == "ASIAN: Asian", "Asian",
                                    ifelse(CAN_RACE_SRTR == "BLACK: Black", "Black",
                                           ifelse(CAN_RACE_SRTR == "WHITE: White", "White", "Other")))),
           can_race_missing = ifelse(CAN_RACE_SRTR == "", 1, 0),
           can_citizen_status = ifelse(CAN_CITIZENSHIP %in% c("", "998"), NA,
                                       ifelse(CAN_CITIZENSHIP == "1: US Citizen", "Citizen",
                                              ifelse(CAN_CITIZENSHIP %in% c("2: RESIDENT ALIEN", 
                                                                            "4: Non-US Citizen/US Resident"), 
                                                     "NonCitizen_Resident",
                                                     ifelse(CAN_CITIZENSHIP %in% c("3: NON-RESIDENT ALIEN, Year Entered US",
                                                                                   "5: Non-US Citizen/Non-US Resident, Traveled to US for Reason Other Than Transplant",
                                                                                   "6: Non-US Citizen/Non-US Resident, Traveled to US for Transplant"),
                                                            "NonCitizen_NonResident",
                                                            "Other")))),
           can_citizen_status_missing = ifelse(CAN_CITIZENSHIP %in% c("", "998"), 1, 0),
           can_education = ifelse(CAN_EDUCATION %in% c("", "998: UNKNOWN"), NA,
                                  ifelse(CAN_EDUCATION %in% c("996: N/A (< 5 YRS OLD)", "1: NONE", "2: GRADE SCHOOL (0-8)"), 
                                         "LessThanHS",
                                         ifelse(CAN_EDUCATION == "3: HIGH SCHOOL (9-12) or GED",
                                                "HighSchool",
                                                ifelse(CAN_EDUCATION == "4: ATTENDED COLLEGE/TECHNICAL SCHOOL",
                                                       "SomeCollege",
                                                       ifelse(CAN_EDUCATION == "5: ASSOCIATE/BACHELOR DEGREE",
                                                              "AssociateBachelor",
                                                              ifelse(CAN_EDUCATION == "6: POST-COLLEGE GRADUATE DEGREE",
                                                                     "Graduate",
                                                                     "Other")))))),
           can_education_missing = ifelse(CAN_EDUCATION %in% c("", "998: UNKNOWN"), 1, 0),
           can_work_for_income = ifelse(CAN_WORK_INCOME %in% c("", "U"), NA,
                                        ifelse(CAN_WORK_INCOME == "Y", 
                                               1,
                                               0)),
           can_work_for_income_missing = ifelse(CAN_WORK_INCOME %in% c("", "U"), 1, 0),
           can_prior_hl = ifelse(CAN_PREV_HL %in% 1, 1, 0),
           can_prior_hr = ifelse(CAN_PREV_HR %in% 1, 1, 0),
           can_prior_in = ifelse(CAN_PREV_IN %in% 1, 1, 0),
           can_prior_ki = ifelse(CAN_PREV_KI %in% 1, 1, 0),
           can_prior_kp = ifelse(CAN_PREV_KP %in% 1, 1, 0),
           can_prior_li = ifelse(CAN_PREV_LI %in% 1, 1, 0),
           can_prior_lu = ifelse(CAN_PREV_LU %in% 1, 1, 0),
           can_prior_pa = ifelse(CAN_PREV_PA %in% 1, 1, 0),
           can_payer = ifelse(CAN_PRIMARY_PAY == "", NA,
                              ifelse(CAN_PRIMARY_PAY == "1: Private insurance", 
                                     "Private",
                                     ifelse(CAN_PRIMARY_PAY %in% c("13: Public insurance - Medicare Unspecified",
                                                                   "14: US/State Govt Agency",
                                                                   "3: Public insurance - Medicare FFS (Fee for Service)",
                                                                   "4: Public insurance - Medicare & Choice",
                                                                   "5: Public insurance - CHIP (Children's Health Insurance Program)",
                                                                   "6: Public insurance - Department of VA",
                                                                   "7: Public insurance - Other government"),
                                            "PublicNonMedicaid",
                                            ifelse(CAN_PRIMARY_PAY == "2: Public insurance - Medicaid",
                                                   "PublicMedicaid",
                                                   "Other")))),
           can_payer_missing = ifelse(CAN_PRIMARY_PAY == "", 1, 0),
           can_weight_kg = ifelse(is.na(CAN_WGT_KG), 
                                  NA,
                                  CAN_WGT_KG),
           can_weight_kg_missing = ifelse(is.na(CAN_WGT_KG), 1, 0),
           can_height_cm = ifelse(is.na(CAN_HGT_CM),
                                  NA,
                                  CAN_HGT_CM),
           can_height_cm_missing = ifelse(is.na(CAN_HGT_CM), 1, 0),
           can_bmi = ifelse(is.na(CAN_BMI),
                            NA,
                            CAN_BMI),
           can_bmi_missing = ifelse(is.na(CAN_BMI), 1, 0),
           
           can_blood_type = ifelse(CAN_ABO %in% c("", "UNK"), NA,
                                   ifelse(CAN_ABO %in% c("A: A", "A1: A1", "A2: A2"), "A",
                                          ifelse(CAN_ABO %in% c("AB: AB", "A1B: A1B", "A2B: A2B"), "AB",
                                                 ifelse(CAN_ABO %in% c("B: B"), "B",
                                                        ifelse(CAN_ABO %in% c("O: O"), "O", "Unexpected"))))),
           can_blood_type_missing = ifelse(CAN_ABO %in% c("", "UNK"), 1, 0),
           can_diabetes_type = ifelse(CAN_DIAB_TY %in% c("", "998: Diabetes Status Unknown"),
                                      NA,
                                      ifelse(CAN_DIAB_TY == "1: No", 
                                             "None",
                                             ifelse(CAN_DIAB_TY == "2: Type I", 
                                                    "Type1",
                                                    ifelse(CAN_DIAB_TY == "3: Type II", 
                                                           "Type2",
                                                           "Other")))),
           can_diabetes_type_missing = ifelse(CAN_DIAB_TY %in% c("", "998: Diabetes Status Unknown"), 1, 0),
           can_prior_cancer = ifelse(CAN_MALIG %in% c("", "U"), NA,
                                     ifelse(CAN_MALIG == "Y", 1, 0)),
           can_prior_cancer_missing = ifelse(CAN_MALIG %in% c("", "U"), 1, 0),
           
           
           can_on_life_support = ifelse(CAN_LIFE_SUPPORT %in% c("", NA), NA,
                                        ifelse(CAN_LIFE_SUPPORT == "Y", 1, 0)),
           
           
           can_had_resist_infect = ifelse(CAN_RESIST_INFECT %in% c("", "U", NA), NA,
                                          ifelse(CAN_RESIST_INFECT %in% "Y", 1, 0)),
           
           can_pa_bp_diast = CAN_PULM_ART_DIAST,
           can_pa_bp_syst = CAN_PULM_ART_SYST,
           can_pa_mean = CAN_PULM_ART_MEAN,
           can_hemo_pcw_mean = CAN_PCW_MEAN,
           can_hemo_pco2 = CAN_PCO2,
           can_smoking_status = ifelse(CAN_HIST_CIGARETTE %in% c("U", "", NA), NA,
                                       ifelse(CAN_HIST_CIGARETTE %in% "Y", 1, 0)),
           can_prior_cardiac_surg = ifelse(CAN_CARDIAC_SURG %in% c("", "U", NA), NA,
                                           ifelse(CAN_CARDIAC_SURG %in% "Y", 1, 0)),
           can_prior_congen_surg = ifelse(CAN_PRIOR_CONGEN_SURG %in% c("", "U", NA), NA,
                                          ifelse(CAN_PRIOR_CONGEN_SURG %in% "Y", 1, 0)),
           
           
           can_peds_cog_develop = ifelse(CAN_COGNITIVE_DEVELOP == "", 
                                         NA,
                                         ifelse(CAN_COGNITIVE_DEVELOP == "998: Not Assessed",
                                                "NotAssessed",
                                                ifelse(CAN_COGNITIVE_DEVELOP == "1: Definite Cognitive delay/impairment",
                                                       "Definite",
                                                       ifelse(CAN_COGNITIVE_DEVELOP == "2: Probable Cognitive delay/impairment",
                                                              "Probable",
                                                              ifelse(CAN_COGNITIVE_DEVELOP == "3: Questionable Cognitive delay/impairment",
                                                                     "Questionable",
                                                                     "None"))))),
           can_peds_motor_delay = ifelse(CAN_MOTOR_DEVELOP == "",
                                         NA,
                                         ifelse(CAN_MOTOR_DEVELOP == "998: Not Assessed", 
                                                "NotAssessed",
                                                ifelse(CAN_MOTOR_DEVELOP == "1: Definite Motor delay/impairment",
                                                       "Definite",
                                                       ifelse(CAN_MOTOR_DEVELOP == "2: Probable Motor delay/impairment",
                                                              "Probable",
                                                              ifelse(CAN_MOTOR_DEVELOP == "3: Questionable Motor delay/impairment",
                                                                     "Questionable",
                                                                     "None"))))),
           can_peds_progress = ifelse(CAN_ACADEMIC_PROGRESS %in% c("", "998: Status Unknown"),
                                      NA,
                                      ifelse(CAN_ACADEMIC_PROGRESS == "996: Not Applicable < 5 years old/ High School graduate or GED",
                                             "NotApplicable",
                                             ifelse(CAN_ACADEMIC_PROGRESS == "1: Within One Grade Level of Peers",
                                                    "Within1Level",
                                                    ifelse(CAN_ACADEMIC_PROGRESS == "2: Delayed Grade Level",
                                                           "Delayed",
                                                           "SpecialEducation")))),
           can_peds_level = ifelse(CAN_ACADEMIC_LEVEL %in% c("", "998: Status Unknown"),
                                   NA,
                                   ifelse(CAN_ACADEMIC_LEVEL == "996: Not Applicable < 5 years old/ High School graduate or GED",
                                          "NotApplicable",
                                          ifelse(CAN_ACADEMIC_LEVEL == "1: Full academic load",
                                                 "FullLoad",
                                                 ifelse(CAN_ACADEMIC_LEVEL == "2: Reduced academic load",
                                                        "ReducedLoad",
                                                        "CannotParticipate")))),
           ctr_id = paste0(CAN_LISTING_CTR_CD, 
                           CAN_LISTING_CTR_TY))
  
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  ### Importing the stathist file for LAS values
  
  lu_las <- read.csv(file.path(root,
                               "Data Sets",
                               paste0("lu_alloc_scores_",
                                      current_saf,
                                      ".csv")),
                     stringsAsFactors = FALSE,
                     header = TRUE)
  
  
  ### Merging the lu_las and cand files
  cand_las <- merge(cand,
                    lu_las,
                    by = "PX_ID",
                    suffixes = c("", "_las"))

  
  ### Selecting the appropriate rows
  cand_las <- cand_las %>%
    mutate(time_from_listing = as.numeric(as.Date(CANHX_LAS_SCORE_BEGIN_DT, format = "%m/%d/%Y") - 
                                            as.Date(CAN_LISTING_DT, format = "%m/%d/%y"))) %>%
    subset(time_from_listing <= 1) %>% 
    group_by(PX_ID) %>%
    mutate(max_time_from_listing = max(time_from_listing)) %>%
    ungroup() %>%
    subset(time_from_listing == max_time_from_listing) %>%
    mutate(las_at_listing = CANHX_LAS_SCORE_CALC)
  
  
  ### Making sure the merge and subsetting worked appropriately
  if(nrow(cand_las) != nrow(cand)){ stop("Merge with LAS stathist file failed.") }
  
  
  ### Cleaning up the workspace
  cand <- cand_las
  
  
  rm(list = c("lu_las", "cand_las"))
  gc()

  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  ### Importing the stathist file for the components of LAS
  
  lu_las_hist_alloc <- read.csv(file.path(root,
                                          "Data Sets",
                                          paste0("lu_las_hist_alloc_",
                                                 current_saf,
                                                 ".csv")),
                                stringsAsFactors = FALSE,
                                header = TRUE)
  
  
  ### Merging the lu_las and cand files
  cand_las_hist_alloc <- merge(cand,
                               lu_las_hist_alloc,
                               by = "PX_ID",
                               suffixes = c("", "_las_component"))
  
  
  ### Selecting the appropriate rows
  cand_las_hist_alloc <- cand_las_hist_alloc %>%
    mutate(days_from_listing = as.numeric(as.Date(CANHX_BEGIN_DT, format = "%m/%d/%Y") - 
                                            as.Date(CAN_LISTING_DT, format = "%m/%d/%y"))) %>%
    group_by(PX_ID) %>%
    mutate(min_days_from_listing = min(days_from_listing)) %>%
    ungroup() %>%
    subset(days_from_listing == min_days_from_listing) %>%
    mutate(fvc_at_listing = CANHX_FVC_ABSOLUTE,
           six_min_walk_at_listing = CANHX_SIX_MIN_WALK,
           cardiac_index_at_listing = CANHX_CI,
           cvp_int_at_listing = CANHX_CVP_INT,
           pco2_at_listing = CANHX_PCO2,
           disease_grouping = CANHX_GROUPING)
  
  
  ### Checking the merge
  if(any(duplicated(cand_las_hist_alloc$PX_ID))){ stop("Unexpected duplicated PX_IDs") }
  if(nrow(cand_las_hist_alloc) != nrow(cand)){ stop("Merge with LAS component stathist file failed.") }
  
  
  ### Cleaning the workspace
  cand <- cand_las_hist_alloc
  
  rm(list = c("lu_las_hist_alloc", "cand_las_hist_alloc"))
  gc()
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  ### Importing the pediatric priority file to obtain more information on candidates aged less than 12
  
  lu_ped_priority <- read.csv(file.path(root,
                                        "Data Sets",
                                        paste0("lu_ped_priority_",
                                               current_saf,
                                               ".csv")),
                              stringsAsFactors = FALSE,
                              header = TRUE)
  
  
  ### Merging with the candidate file
  ped_priority <- merge(cand,
                        lu_ped_priority,
                        by = "PX_ID")
  
  
  ### Selecting candidates listed at less than 12 years old and selecting the first row within 30 days of listing
  ped_priority <- ped_priority %>%
    subset(can_age < 12) %>%
    mutate(days_from_listing = as.numeric(as.Date(CAN_PL_EFFECTIVE_DT, format = "%m/%d/%Y") -
                                            as.Date(CAN_LISTING_DT, format = "%m/%d/%y"))) %>%
    subset(days_from_listing <= 30) %>%
    group_by(PX_ID) %>%
    mutate(min_days_from_listing = min(days_from_listing)) %>%
    ungroup() %>%
    subset(days_from_listing == min_days_from_listing) %>%
    group_by(PX_ID) %>%
    summarise(max_ped_priority_at_listing = min(CAN_PL_PRIORITY))
  
  
  ### Merging back with the candidate file
  cand_ped_priority <- merge(cand,
                             ped_priority,
                             by = "PX_ID",
                             all.x = TRUE)
  
  
  ### Checking the merge
  if(nrow(cand_ped_priority) != nrow(cand)){ stop("Merge with the pediatric priority file failed.") }
  
  
  ### Setting missing values for pediatric allocation priority to 2 due to the implicit assumptions of allocation
  cand_ped_priority <- cand_ped_priority %>%
    mutate(max_ped_priority_at_listing = ifelse(is.na(max_ped_priority_at_listing) & {can_age < 12},
                                                1,
                                                max_ped_priority_at_listing))
  
  
  ### Cleaning the workspace
  cand <- cand_ped_priority
  
  rm(list = c("lu_ped_priority", "ped_priority", "cand_ped_priority"))
  gc()
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  ### Identifying candidates listed for kipa organs
  
  kipa_candidates <- read.csv(file.path(root,
                                        "Data Sets",
                                        paste0("kipa_multi_candidates_",
                                               current_saf,
                                               ".csv")),
                              stringsAsFactors = FALSE,
                              header = TRUE)
  
  
  ### merging the kipa candidate file and cand
  cand_kipa_cand <- merge(cand,
                          kipa_candidates,
                          by = "PERS_ID",
                          all.x = TRUE,
                          suffixes = c("", "_kipa_cands"))
  
  
  cand_kipa_cand <- cand_kipa_cand %>%
    mutate(kipa_ctr_id = ifelse(is.na(CAN_LISTING_CTR_CD_kipa_cands),
                                NA,
                                paste0(CAN_LISTING_CTR_CD_kipa_cands, 
                                       CAN_LISTING_CTR_TY_kipa_cands)),
           listed_at_ctr = ifelse(is.na(kipa_ctr_id),
                                  0,
                                  ifelse(kipa_ctr_id == ctr_id,
                                         1, 
                                         0)),
           listed_at_listing = ifelse(is.na(CAN_LISTING_DT_kipa_cands),
                                      0,
                                      ifelse({CAN_LISTING_DT == ""} | {CAN_LISTING_DT_kipa_cands == ""}, 
                                             0,
                                             ifelse((as.Date(CAN_LISTING_DT_kipa_cands, format = "%m/%d/%y") - 30) <= as.Date(CAN_LISTING_DT, format = "%m/%d/%y") & 
                                                      {{CAN_REM_DT_kipa_cands == ""} | 
                                                          {as.Date(CAN_LISTING_DT, format = "%m/%d/%y") <= as.Date(CAN_REM_DT_kipa_cands, format = "%m/%d/%y")}},
                                                    1, 
                                                    0)))) %>%
    group_by(PX_ID) %>%
    summarise(listed_for_ki = 1 * {sum({listed_at_ctr == 1} & {listed_at_listing == 1} & {WL_ORG_kipa_cands %in% "KI: Kidney"}) > 0},
              listed_for_pa = 1 * {sum({listed_at_ctr == 1} & {listed_at_listing == 1} & {WL_ORG_kipa_cands %in% "PA: Pancreas"}) > 0},
              listed_for_kp = 1 * {sum({listed_at_ctr == 1} & {listed_at_listing == 1} & {WL_ORG_kipa_cands %in% "KP: Kidney-Pancreas"}) > 0})
    
  
  
  ### merging back with the cand file and cleaning the workspace
  cand_kipa <- merge(cand,
                     cand_kipa_cand,
                     by = "PX_ID")
  
  
  if(nrow(cand_kipa) != nrow(cand)){ stop("Process for determining multilisted status in KIPA failed.") }
  
  
  cand <- cand_kipa
  
  
  rm(list = c("cand_kipa", "cand_kipa_cand", "kipa_candidates"))
  
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  ### Identifying candidates listed for non-liver liin organs
  
  liin_candidates <- read.csv(file.path(root,
                                        "Data Sets",
                                        paste0("liin_multi_candidates_",
                                               current_saf,
                                               ".csv")),
                              stringsAsFactors = FALSE,
                              header = TRUE)
  
  
  ### merging the kipa candidate file and cand
  cand_liin_cand <- merge(cand,
                          liin_candidates,
                          by = "PERS_ID",
                          all.x = TRUE,
                          suffixes = c("", "_liin_cands"))
  
  
  cand_liin_cand <- cand_liin_cand %>%
    mutate(liin_ctr_id = ifelse(is.na(CAN_LISTING_CTR_CD_liin_cands),
                                NA,
                                paste0(CAN_LISTING_CTR_CD_liin_cands, 
                                       CAN_LISTING_CTR_TY_liin_cands)),
           listed_at_ctr = ifelse(is.na(liin_ctr_id),
                                  0,
                                  ifelse(liin_ctr_id == ctr_id,
                                         1, 
                                         0)),
           listed_at_listing = ifelse(is.na(CAN_LISTING_DT_liin_cands),
                                      0,
                                      ifelse({CAN_LISTING_DT == ""} | {CAN_LISTING_DT_liin_cands == ""}, 
                                             0,
                                             ifelse((as.Date(CAN_LISTING_DT_liin_cands, format = "%m/%d/%y") - 30) <= as.Date(CAN_LISTING_DT, format = "%m/%d/%y") & 
                                                      {{CAN_REM_DT_liin_cands == ""} | 
                                                          {as.Date(CAN_LISTING_DT, format = "%m/%d/%y") <= as.Date(CAN_REM_DT_liin_cands, format = "%m/%d/%y")}},
                                                    1, 
                                                    0)))) %>%
    group_by(PX_ID) %>%
    summarise(listed_for_in = 1 * {sum({listed_at_ctr == 1} & {listed_at_listing == 1} & {WL_ORG_liin_cands %in% "IN: Intestine"}) > 0},
              listed_for_li = 1 * {sum({listed_at_ctr == 1} & {listed_at_listing == 1} & {WL_ORG_liin_cands %in% "LI: Liver"}) > 0})
  
  
  
  ### merging back with the cand file and cleaning the workspace
  cand_liin <- merge(cand,
                     cand_liin_cand,
                     by = "PX_ID")
  
  
  if(nrow(cand_liin) != nrow(cand)){ stop("Process for determining multilisted status in LIIN failed.") }
  
  
  cand <- cand_liin
  
  
  rm(list = c("cand_liin", "cand_liin_cand", "liin_candidates"))
  gc()
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  ### Identifying candidates listed for thor organs
  
  thor_candidates <- read.csv(file.path(root,
                                        "Data Sets",
                                        paste0("thor_multi_candidates_",
                                               current_saf,
                                               ".csv")),
                              stringsAsFactors = FALSE,
                              header = TRUE)
  
  
  ### Removing the candidates listed for lung transplant
  thor_candidates <- thor_candidates %>%
    subset(WL_ORG != "LU: Lung")
  
  
  ### merging the thor candidate file and cand
  cand_thor_cand <- merge(cand,
                          thor_candidates,
                          by = "PERS_ID",
                          all.x = TRUE,
                          suffixes = c("", "_thor_cands"))
  
  
  cand_thor_cand <- cand_thor_cand %>%
    mutate(thor_ctr_id = ifelse(is.na(CAN_LISTING_CTR_CD_thor_cands),
                                NA,
                                paste0(CAN_LISTING_CTR_CD_thor_cands, 
                                       CAN_LISTING_CTR_TY_thor_cands)),
           listed_at_ctr = ifelse(is.na(thor_ctr_id),
                                  0,
                                  ifelse(thor_ctr_id == ctr_id,
                                         1, 
                                         0)),
           listed_at_listing = ifelse(is.na(CAN_LISTING_DT_thor_cands),
                                      0,
                                      ifelse({CAN_LISTING_DT == ""} | {CAN_LISTING_DT_thor_cands == ""}, 
                                             0,
                                             ifelse((as.Date(CAN_LISTING_DT_thor_cands, format = "%m/%d/%y") - 30) <= as.Date(CAN_LISTING_DT, format = "%m/%d/%y") & 
                                                      {{CAN_REM_DT_thor_cands == ""} | 
                                                          {as.Date(CAN_LISTING_DT, format = "%m/%d/%y") <= as.Date(CAN_REM_DT_thor_cands, format = "%m/%d/%y")}},
                                                    1, 
                                                    0)))) %>%
    group_by(PX_ID) %>%
    summarise(listed_for_hl = 1 * {sum({listed_at_ctr == 1} & {listed_at_listing == 1} & {WL_ORG_thor_cands %in% "HL: Heart/Lung"}) > 0},
              listed_for_hr = 1 * {sum({listed_at_ctr == 1} & {listed_at_listing == 1} & {WL_ORG_thor_cands %in% "HR: Heart"}) > 0})
  
  
  
  ### merging back with the cand file and cleaning the workspace
  cand_thor <- merge(cand,
                     cand_thor_cand,
                     by = "PX_ID")
  
  
  if(nrow(cand_thor) != nrow(cand)){ stop("Process for determining multilisted status in KIPA failed.") }
  
  
  cand <- cand_thor
  
  
  rm(list = c("cand_thor", "cand_thor_cand", "thor_candidates"))
  gc()
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  ### Determining distance between candidate ZIP code and center ZIP code
  institutes <- read.csv(file.path(root,
                                   "Data Sets",
                                   paste0("institutes_", current_saf, ".csv")),
                         stringsAsFactors = FALSE,
                         header = TRUE)
  
  
  ### Merging the institutes file with the cand file to get the ZIP code for the center
  cand_ctr_zip <- institutes %>%
    mutate(ctr_id = paste0(CTR_CD, substr(CTR_TY, 1, 3)),
           ctr_zip = zip) %>%
    select(ctr_id, 
           ctr_zip) %>%
    merge(cand, 
          by = "ctr_id",
          all.y = TRUE)
  
  
  # Checking the merge
  if(nrow(cand) != nrow(cand_ctr_zip)){ stop("Merge with institutions file failed.") }
  
  
  # Updating the cand file and removing cand_ctr_zip
  cand <- cand_ctr_zip
  
  rm(list = "cand_ctr_zip"); gc()
  
  
  ### Reading the zipcode file
  zipcodes <- read.csv(file.path(root,
                                 "Data Sets",
                                 paste0("zipcodes_", current_saf, ".csv")),
                       stringsAsFactors = FALSE,
                       header = TRUE)
  
  
  ### Merging to get candidate latitude and longitude
  cand_zip <- zipcodes %>%
    rename(can_zip = ZIP) %>%
    merge(cand,
          by = "can_zip",
          all.y = TRUE)
  
  
  # Checking the merge
  if(nrow(cand) != nrow(cand_zip)){ stop("Merge with zipcodes file failed for candidate zip.") }
  
  
  # Updating the cand file and removing cand_zip
  cand <- cand_zip
  
  rm(list = "cand_zip"); gc()
  
  
  ### Merging to get center latitude and longitude
  cand_ctr_zip <- zipcodes %>%
    rename(ctr_zip = ZIP,
           ctr_latitude = can_latitude,
           ctr_longitude = can_longitude) %>%
    merge(cand,
          by = "ctr_zip",
          all.y = TRUE)
  
  
  # Checking the merge
  if(nrow(cand) != nrow(cand_ctr_zip)){ stop("Merge with zipcodes file failed for center zip.") }
  
  
  # Updating the cand file and removing cand_ctr_zip
  cand <- cand_ctr_zip
  
  rm(list = "cand_ctr_zip"); gc()
  
  
  ### Finally calculating the distance between the ZIP codes
  xs <- cbind(cand$can_longitude, cand$can_latitude)
  ys <- cbind(cand$ctr_longitude, cand$ctr_latitude)
  
  
  cand$distance_from_ctr <- sapply(1:nrow(xs), 
                                   function(x){ 
                                     if(is.na(xs[x,1]) | is.na(ys[x,1])){ tmp <- NA }
                                     if({! is.na(xs[x,1])} & {! is.na(ys[x,1])}){ 
                                       tmp <- spDists(x = xs[x,, drop = FALSE], 
                                                      y = ys[x,, drop = FALSE], 
                                                      longlat = TRUE, 
                                                      diagonal = TRUE) 
                                     }
                                     
                                     return(tmp)
                                   })
  
  cand$log2_distance_from_ctr <- log2(cand$distance_from_ctr + 1)

  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  ### Calculating the important followup variables for survival from listing
  cand <- create_outcome_vars(cand,
                              max_followup,
                              cohort_start,
                              cohort_end)
  
  
  return(cand)
}


