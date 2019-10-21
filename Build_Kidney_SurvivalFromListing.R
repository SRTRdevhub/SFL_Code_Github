###
### This function will build the survival from listing model for kidney candidates
###


### Loading the necessary packages
library(dplyr)
library(sp)
library(AUC)
library(lattice)
library(mice)
library(glmnet)
library(survival)



### Directory information (requires updating)
root <- file.path("SFL_Code_Github")


### Organ information
organ        <- "Kidney"
organ_ab     <- "ki"
organ_saf_ab <- "kipa"
current_saf  <- "SAF1808"



# Type of model
model_type <- "Build"


# Cohort definitions:
#   - 'cohort_start' and 'cohort_end': The period-prevalent cohort period for the survival from listing
#   - 'max_followup': The maximum followup (in years) allowed
cohort_start <- "20150701"
cohort_end   <- "20170630"
max_followup <- 7



### The cutpoints for the piecewise hazard function
follow_up_cuts <- 365 * 0:{max_followup - 1}



# Reading the control table
control <- read.csv(file.path(root,
                              organ,
                              "Control",
                              paste0("Control_", model_type, ".csv")),
                    stringsAsFactors = FALSE,
                    header = TRUE)



### Creating a vector from the control table that list the candidate variables that need extracting from the SAF

cand_vars <- control %>% 
  subset(SAF_EXPORT) %>% 
  select(var_name) %>% 
  unlist()



### These are objects necessary to correctly run the multiple imputation.  They will ultimately depend on the variables derived in '...DataCleaning.R'

mi_num <- 10

adult_mi_covariates <- control$var_name[control$adult_build_include]
adult_mi_cov_method <- control$MI_TYPE[control$adult_build_include]

peds_mi_covariates <- control$var_name[control$ped_build_include]
peds_mi_cov_method <- control$MI_TYPE[control$ped_build_include]


keep_vars <- control$var_name[control$MI_KEEP]

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

###  Sourcing the necessary worker functions

source(file.path(root, "WorkerFunctions", paste0(organ, "DataCleaning.R")))
source(file.path(root, "WorkerFunctions", "sas_saf_grab.R"))
source(file.path(root, "WorkerFunctions", "create_var_for_levels.R"))
source(file.path(root, "WorkerFunctions", "create_dummy_variables.R"))
source(file.path(root, "WorkerFunctions", "create_linear_splines.R"))
source(file.path(root, "WorkerFunctions", "create_timevarying_vars.R"))
source(file.path(root, "WorkerFunctions", "create_followup_vars.R"))
source(file.path(root, "WorkerFunctions", "mi.R"))
source(file.path(root, "WorkerFunctions", "estimate_model_i.R"))
source(file.path(root, "WorkerFunctions", "piecewise_cv_fold.R"))
source(file.path(root, "WorkerFunctions", "create_coxzph_figs.R"))
source(file.path(root, "WorkerFunctions", "create_outcome_vars.R"))

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

### Exporting the SAF data required for the survival from listing model

if(! file.exists(file.path(root, "Data Sets", paste0(organ_ab, "_candidate_file_", current_saf, ".csv")))){
  saf_extract <- sas_saf_grab(file.path(root, "WorkerFunctions"), 
                              file.path(root, "Data Sets"), 
                              current_saf, 
                              organ_ab, 
                              organ_saf_ab, 
                              cand_vars)
}


###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

### Cleaning the data
if(! file.exists(file.path(root, organ, "MI_Output", paste0("mi_", organ_ab, "_adult_data_", mi_num, ".Rdata")))){
  
  sfl_data <- KidneyDataCleaning(root, 
                                 organ, 
                                 organ_ab, 
                                 current_saf, 
                                 cohort_start, 
                                 cohort_end, 
                                 max_followup)
  
  
  ### Trimming the continuous variables included in the model
  trim_vars <- control$var_name[{control$adult_build_include | control$ped_build_include} & {control$MI_TYPE == "continuous"}]
  
  
  # Removing age from trim_vars
  trim_vars <- trim_vars[trim_vars != "can_age"]
  
  
  for(k in trim_vars){
    
    adult_boo <- {sfl_data$can_age >= 18}
    peds_boo  <- {sfl_data$can_age < 18}
    
    
    adult_q1  <- quantile(sfl_data[[k]][adult_boo], 
                          probs = 1/100,
                          na.rm = TRUE)
    adult_q99 <- quantile(sfl_data[[k]][adult_boo], 
                          probs = 99/100,
                          na.rm = TRUE)
    
    
    peds_q1  <- quantile(sfl_data[[k]][peds_boo], 
                         probs = 1/100,
                         na.rm = TRUE)
    peds_q99 <- quantile(sfl_data[[k]][peds_boo], 
                         probs = 99/100,
                         na.rm = TRUE)
    
    
    sfl_data[[k]][adult_boo] <- ifelse(sfl_data[[k]][adult_boo] < adult_q1,
                                       adult_q1,
                                       ifelse(sfl_data[[k]][adult_boo] > adult_q99,
                                              adult_q99,
                                              sfl_data[[k]][adult_boo]))
    
    
    sfl_data[[k]][peds_boo] <- ifelse(sfl_data[[k]][peds_boo] < peds_q1,
                                       peds_q1,
                                       ifelse(sfl_data[[k]][peds_boo] > peds_q99,
                                              peds_q99,
                                              sfl_data[[k]][peds_boo]))
  }
  
  
  sfl_data_peds <- sfl_data %>%
    subset(can_age < 18)
  
  
  sfl_data_adult <- sfl_data %>%
    subset(can_age >= 18)
  
  
  # Removing the primary slf_data file
  rm(list = c("sfl_data", "adult_boo", "peds_boo"))
  gc()
  
  
  ### Saving the raw datasets for (1) identifying missingness for the LBV adjustment, and (2) the appropriate reference level
  save(sfl_data_peds,
       file = file.path(root,
                        organ,
                        "sfl_data_peds.Rdata"))
  
  
  save(sfl_data_adult,
       file = file.path(root,
                        organ,
                        "sfl_data_adult.Rdata"))
  
  
  mi(sfl_data_peds,
     mi_num,
     peds_mi_covariates,
     peds_mi_cov_method,
     keep_vars,
     "peds",
     file.path(root, organ, "MI_Output"))
  
  
  
  mi(sfl_data_adult,
     mi_num,
     adult_mi_covariates,
     adult_mi_cov_method,
     keep_vars,
     "adult",
     file.path(root, organ, "MI_Output"))
  
  
  # Removing the data sets to clear space
  rm(list = c("sfl_data_adult", "sfl_data_peds"))
  gc()
}

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

### Importing the Results of the multiple imputation

if(! file.exists(file.path(root, organ, "piecewise_separate", paste0("adult_cv_prob_", mi_num, ".csv")))){
  
  model_time <- NULL
  
  for(m in 1:mi_num){
    
    cox_time <- estimate_model_i(root,
                                 organ,
                                 organ_ab,
                                 "adult",
                                 m,
                                 follow_up_cuts,
                                 adult_mi_covariates,
                                 adult_mi_cov_method,
                                 control,
                                 "coxph_constant")
    
    
    pemc_time <- estimate_model_i(root,
                                  organ,
                                  organ_ab,
                                  "adult",
                                  m,
                                  follow_up_cuts,
                                  adult_mi_covariates,
                                  adult_mi_cov_method,
                                  control,
                                  "piecewise_constant")
    
    
    pemti_time <- estimate_model_i(root,
                                   organ,
                                   organ_ab,
                                   "adult",
                                   m,
                                   follow_up_cuts,
                                   adult_mi_covariates,
                                   adult_mi_cov_method,
                                   control,
                                   "piecewise_varying")
    
    
    pemtd_time <- estimate_model_i(root,
                                   organ,
                                   organ_ab,
                                   "adult",
                                   m,
                                   follow_up_cuts,
                                   adult_mi_covariates,
                                   adult_mi_cov_method,
                                   control,
                                   "piecewise_separate")
    
    
    
    m_time <- c(cox_time,
                pemc_time,
                pemti_time,
                pemtd_time)
    
    
    
    write.table(matrix(m_time,
                       nrow = 1),
                file.path(root, 
                          organ,
                          paste0(organ_ab, "_model_estimation_times.csv")),
                row.names = FALSE,
                col.names = FALSE,
                sep = ",",
                append = TRUE)
    
    
    model_time <- rbind(model_time,
                        m_time)
  }
}


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

### Running the unpenalized PEM with time-varying effects

model_time <- NULL

for(m in 1:10){
  
  glm_time <- estimate_model_i(root,
                               organ,
                               organ_ab,
                               "adult",
                               m,
                               follow_up_cuts,
                               adult_mi_covariates,
                               adult_mi_cov_method,
                               control,
                               "gam_pem")
  
  
  model_time <- c(model_time,
                  glm_time)
}


write.table(matrix(model_time, ncol = 1),
            file.path(root, 
                      organ,
                      paste0(organ_ab, "_glm_model_estimation_times.csv")),
            row.names = FALSE,
            col.names = FALSE,
            sep = ",")









































