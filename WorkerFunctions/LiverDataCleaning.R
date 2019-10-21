
# Cleans the SAF for use in the survival from listing model


LiverDataCleaning <- function(root, 
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
  
  
  ### Based on the cohort definition, these candidates cannot be in the cohort
  cand <- cand %>%
    subset(as.Date(CAN_LISTING_DT, format = "%m/%d/%y") >= 
             (as.Date(cohort_start, format = "%Y%m%d") - (365 * max_followup))) %>%
    subset(WL_ORG == "LI: Liver") %>%
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
           can_prior_abdom_surg = ifelse(CAN_PREV_ABDOM_SURG %in% c("", "U", NA), NA,
                                         ifelse(CAN_PREV_ABDOM_SURG == "Y", 1, 0)),
           can_bacterial_perit = ifelse(CAN_BACTERIA_PERIT %in% c("", "U", NA), NA,
                                        ifelse(CAN_BACTERIA_PERIT == "Y", 1, 0)),
           can_portal_vein_thrombosis = ifelse(CAN_PORTAL_VEIN %in% c("", "U", ""), NA,
                                               ifelse(CAN_PORTAL_VEIN == "Y", 1, 0)),
           can_tipss = ifelse(CAN_TIPSS %in% c("", "U", NA), NA,
                              ifelse(CAN_TIPSS == "Y", 1, 0)),
           
           
           
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
  
  
  
  ### Reading the diagnosis mapping
  dgn_map <- read.csv(file.path(root,
                                "Diagnosis.csv"),
                      stringsAsFactors = FALSE,
                      header = TRUE)
  
  
  dgn_map$organ_name <- sapply(dgn_map$TABLE.LABEL, 
                               function(x){ strsplit(x, 
                                                     split = ": ")[[1]][1] })
  
  dgn_map <- dgn_map %>%
    subset(organ_name %in% c("Liver ", "Liver", "Other Specify"))
  
  
  # Cleaning the SAS formats to better match with the diagnosis dictionary
  cand <- cand %>% 
    mutate(can_dgn_clean_char = ifelse(substr(CAN_DGN, 1, 3) == "999",
                                       substr(CAN_DGN, 1, 3),
                                       substr(CAN_DGN, 1, 4)),
           can_dgn_clean = as.numeric(can_dgn_clean_char))
  
  
  cand$can_dgn_grp <- dgn_map$TABLE.LABEL[match(cand$can_dgn_clean, dgn_map$OPTN.CODE)]
  
  
  # There are three diagnsoses not currently included in the dictionary: We will move them to the other category (after checking for unexpected diagnoses)
  dictionary_miss_dgn <- names(table(cand$CAN_DGN[is.na(cand$can_dgn_grp)]))
  expected_miss       <- dictionary_miss_dgn %in% c("", "3072: KI:HEPATORENAL SYNDROME", "3073: KI:LITHIUM TOXICITY", "3074: KI:HIV NEPHROPATHY")
  
  if(any(! expected_miss)){ stop("Diagnosis dictionary did not include every diagnosis") }
  
  cand <- cand %>%
    mutate(can_dgn_grp = ifelse(is.na(can_dgn_grp) & {CAN_DGN == ""},
                                NA,
                                ifelse(is.na(can_dgn_grp) & {CAN_DGN != ""},
                                       "Liver : Other",
                                       can_dgn_grp)),
           can_dgn_grp_col = ifelse(is.na(can_dgn_grp),
                                    NA,
                                    ifelse(can_dgn_grp == "Liver : Acute Hepatic Necrosis",
                                           "AHN",
                                           ifelse(can_dgn_grp %in% "Liver : Biliary Atresia",
                                                  "BiliaryAtresia",
                                                  ifelse(can_dgn_grp == "Liver : Cholestatic Liver Disease/Cirrhosis",
                                                         "CholestaticCirrhosis",
                                                         ifelse(can_dgn_grp == "Liver : Malignant Neoplasms",
                                                                "MalignantNeoplasms",
                                                                ifelse(can_dgn_grp == "Liver : Metabolic Diseases",
                                                                       "MetabolicDiseases",
                                                                       ifelse(can_dgn_grp == "Liver : Non-Cholestatic Cirrhosis",
                                                                              "NonCholestaticCirrhosis",
                                                                              ifelse(can_dgn_grp %in% c("Liver : Other", 
                                                                                                        "Other Specify"),
                                                                                     "Other",
                                                                                     "unexpected level")))))))))
  
  
  if(any(cand$can_dgn_grp_col %in% "unexpected level")){ stop("Unexpected diagonsis level") }
  
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  ################################################################################################################################################
  
  ### Importing the stathist file
  
  li_stathist <- read.csv(file.path(root,
                                    "Data Sets",
                                    paste0("li_alloc_meld_",
                                           current_saf,
                                           ".csv")),
                          stringsAsFactors = FALSE,
                          header = TRUE)
  
  
  
  ### Selecting the appropriate rows from the stathist file
  li_stathist <- li_stathist %>%
    subset(PX_ID %in% cand$PX_ID) %>%
    mutate(time_from_listing = as.numeric(as.Date(CANHX_BEGIN_DT, format = "%m/%d/%Y") - as.Date(CAN_LISTING_DT, format = "%m/%d/%y"))) %>%
    subset(time_from_listing <= 30) %>%
    group_by(PX_ID) %>%
    mutate(min_time_from_listing = min(time_from_listing)) %>%
    ungroup() %>%
    subset(time_from_listing == min_time_from_listing)
    
  
  
  
  ### Merging the cand and stathist files
  cand_stathist <- merge(cand,
                         li_stathist,
                         by = 'PX_ID',
                         all.x = TRUE,
                         suffixes = c("", "_stathist"))
  
    
  if(nrow(cand) != nrow(cand_stathist)){ stop("Merge with the stathist file failed.") }
  
  
  
  ### Cleaning some of the stathist variables
  cand_stathist <- cand_stathist %>%
    mutate(lab_meld = as.numeric(substr(CANHX_SRTR_LAB_MELD, 1, 4)) - 6200,
           log2_cand_inr = ifelse(is.na(CANHX_INR),
                                  NA,
                                  log2(CANHX_INR)),
           can_ascites = ifelse(CANHX_ASCITES %in% c("", NA),
                                NA,
                                ifelse(CANHX_ASCITES == "1: Absent",
                                       1,
                                       ifelse(CANHX_ASCITES == "2: Slight",
                                              2,
                                              ifelse(CANHX_ASCITES == "3: Moderate",
                                                     3,
                                                     NA)))),
           log2_cand_bili = ifelse(is.na(CANHX_BILI),
                                   NA,
                                   log2(CANHX_BILI + 1)),
           log2_cand_s_cr = ifelse(is.na(CANHX_SERUM_CREAT),
                                   NA,
                                   log2(CANHX_SERUM_CREAT + 1)),
           cand_hcc_exc = ifelse({CANHX_EXC_DIAG_HCC1 %in% 1} | 
                                   {CANHX_EXC_DIAG_HCC2 %in% 1} | 
                                   {CANHX_EXC_DIAG_HCC_NOPOLICY %in% 1}, 1, 0),
           cand_other_exc = ifelse(CANHX_EXC_DIAG_OTHER %in% 1, 1, 0))
  
  
  ### Cleaning up the workspace
  cand <- cand_stathist
  
  
  rm(list = c("cand_stathist", "li_stathist"))
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
  
  
  ### Removing the candidates listed for liver transplant
  liin_candidates <- liin_candidates %>%
    subset(WL_ORG != "LI: Liver")
  
  
  
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
    summarise(listed_for_in = 1 * {sum({listed_at_ctr == 1} & {listed_at_listing == 1} & {WL_ORG_liin_cands %in% "IN: Intestine"}) > 0})
  
  
  
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
              listed_for_hr = 1 * {sum({listed_at_ctr == 1} & {listed_at_listing == 1} & {WL_ORG_thor_cands %in% "HR: Heart"}) > 0},
              listed_for_lu = 1 * {sum({listed_at_ctr == 1} & {listed_at_listing == 1} & {WL_ORG_thor_cands %in% "LU: Lung"}) > 0})
  
  
  
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
  
  ### Adding ecological risk-factors
  
  ruca <- read.csv(file.path(root,
                             "Data Sets",
                             "Ecological_Data",
                             "cleaned_ruca.csv"),
                   stringsAsFactors = FALSE,
                   header = TRUE)
  
  
  ### Merging with the cand file
  cand_ruca <- ruca %>%
    rename(can_zip = ZIPCODEN) %>%
    merge(cand,
          by = "can_zip",
          all.y = TRUE)
  
  
  ### Puerto Rico is typically not in the RUCA database
  
  
  
  hh_income <- read.csv(file.path(root,
                                  "Data Sets",
                                  "Ecological_Data",
                                  "cleaned_hh_data.csv"),
                        stringsAsFactors = FALSE,
                        header = TRUE)
  
  
  ### Merging with the cand file
  cand_hh <- hh_income %>%
    rename(can_zip = zip) %>%
    merge(cand,
          by = "can_zip",
          all.y = TRUE)
  
  
  ### How should foreign nationals be handled for median?
  
  rm(list = c("cand_ruca", "cand_hh", "ruca", "hh_income"))
  gc()
  
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


