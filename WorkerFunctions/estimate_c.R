

estimate_c <- function(root,
                       organ,
                       organ_ab,
                       age_type,
                       m,
                       model_type){
  
  ### root: The base directory for the model building
  ### organ: The organ name (e.g., Kidney) for model building
  ### organ_ab: The organ abbreviation (e.g., ki) for model building
  ### age_type: Identifies whether the model is for pediatric or adult candidates
  ### m: Identifies the iteration of the multiple imputation
  ### model_type: Specifying the type of model to fit
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Loading the appropriate iteration of MI data
  load(file.path(root, 
                 organ, 
                 "MI_Output", 
                 paste0("mi_", organ_ab, "_", age_type, "_data_", m, ".Rdata")))
  
  
  curr_data <- mi_data
  
  
  rm(list = "mi_data")
  gc()
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Importing the cross-validated probabilities. The structure of the data set depends on the model_type.
  
  if(model_type == "coxph_constant"){
    
    cv_data <- read.csv(file.path(root, 
                                  organ, 
                                  model_type,
                                  paste0(age_type, "_cv_prob_", m, ".csv")),
                        stringsAsFactors = FALSE,
                        header = FALSE, 
                        col.names = c("PX_ID", "cv_exp", "surv_1yr", "surv_2yr"))
  }
  if(model_type != "coxph_constant"){
    
    cv_data <- read.csv(file.path(root, 
                                  organ, 
                                  model_type,
                                  paste0(age_type, "_cv_prob_", m, ".csv")),
                        stringsAsFactors = FALSE,
                        header = FALSE, 
                        col.names = c("PX_ID", "cv_exp", "cv_deviance", "surv_1yr", "surv_2yr"))
  }
  
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Data manipulation before estimating the Brier Score
  
  pred_data <- merge(cv_data,
                     curr_data,
                     by = "PX_ID",
                     suffixes = c("_cv", ""))
  
  
  ### Determining the 'right_time' equivalent for the maximum follow-up date and creating an indicator for whether the probability requires estimation
  pred_data <- pred_data %>%
    mutate(time_in_cohort = right_time - left_time + 1,
           Delta_1yr = ifelse(time_in_cohort >= (1 * 365 + 1), 
                              1,
                              as.numeric(as.character(event))),
           Delta_2yr = ifelse(time_in_cohort >= (2 * 365 + 1), 
                              1,
                              as.numeric(as.character(event))),
           Event_1yr = ifelse(time_in_cohort > (1 * 365 + 1), 
                              0,
                              as.numeric(as.character(event))),
           Event_2yr = ifelse(time_in_cohort > (2 * 365 + 1), 
                              0,
                              as.numeric(as.character(event))))
  

  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Estimating the Brier Score at the end of follow-up, which depends on the model_type
  
  c_1yr <- timeROC::timeROC(pred_data$time_in_cohort,
                            as.numeric(as.character(pred_data$event)),
                            marker = 1 - pred_data$surv_1yr,
                            cause = 1,
                            times  = c(365))$AUC[2]
  
  
  c_2yr <- timeROC::timeROC(pred_data$time_in_cohort,
                            as.numeric(as.character(pred_data$event)),
                            marker = 1 - pred_data$surv_1yr,
                            cause = 1,
                            times  = c(2 * 365))$AUC[2]

  
  return(c(c_1yr, c_2yr))
}




