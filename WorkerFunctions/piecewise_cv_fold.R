###
### Purpoose: Perform a fold of the cross-validation for a piecewise exponential model
###
### Author: Andrew Wey
###

piecewise_cv_fold <- function(curr_imputed,
                              curr_vars,
                              curr_model_statement,
                              curr_model_matrix,
                              curr_outcome,
                              curr_offset,
                              curr_ids,
                              curr_fold_vec,
                              curr_fold,
                              piecewise_type,
                              root,
                              organ,
                              age_type,
                              mi_i){
  
  ### curr_imputed: The data set containing the currently imputed values
  ### curr_vars: A vector of covariate names (before creating time-varying components)
  ### curr_model_statement: The current version of the model statement
  ### curr_model_matrix: The model matrix for the full data set.
  ### curr_outcome: The outcome vector for the full data set. The vector has the same order as curr_model_matrix.
  ### curr_offset: The offset vector for the full data set. The vector has the same order as curr_model_matrix.
  ### curr_fold_vec: The vector that identifies the fold for each unique PX_ID.
  ### curr_fold: The current fold being estimated.
  ### piecewise_type: Identifies the type of piecewise exponential model being estimated. Specifically, whether the covariates change over time.
  ### root: The base directory for the model building. Determines directory information for the fold predictions.
  ### organ: The organ name (e.g., Kidney) for model building. Determines directory information for the fold predictions.
  ### age_type: Identifies whether the model is for pediatric or adult candidates. Determines directory information for the fold predictions.
  ### mi_i: Identifies the iteration of the multiple imputation. Determines directory information for the fold predictions.
  
  
  if(! piecewise_type %in% c("piecewise_constant", "piecewise_varying")){ stop("Incorrect piecewise_type argument for piecewise_cv_fold function") }
  
  
  ###################################################################################################################################################
  ###################################################################################################################################################
  ###################################################################################################################################################
  ###################################################################################################################################################
  
  ### Fitting the model for cross-validation
  
  # Model matrix for current CV iteration
  curr_cv_matrix <- curr_model_matrix[curr_fold_vec != curr_fold,]
  
  # Model outcome for current CV iteration
  curr_cv_outcome <- curr_outcome[curr_fold_vec != curr_fold]
  
  # Model offset for current CV iteration
  curr_cv_offset <- curr_offset[curr_fold_vec != curr_fold]
  
  
  ### Time to fit the actual model(!)
  curr_cv_k <- cv.glmnet(x = curr_cv_matrix,
                         y = curr_cv_outcome,
                         family = "poisson",
                         offset = curr_cv_offset)
  
  
  # Cleaning up the workspace
  rm(list = c("curr_cv_matrix", "curr_cv_outcome", "curr_cv_offset"))
  gc()
  
  ###################################################################################################################################################
  ###################################################################################################################################################
  ###################################################################################################################################################
  ###################################################################################################################################################
  
  ### Estimating the cross-validated number of expected events and deviance
  
  # Model matrix for current CV iteration
  cv_matrix <- curr_model_matrix[curr_fold_vec == curr_fold,]
  
  # Model outcome for current CV iteration
  cv_outcome <- curr_outcome[curr_fold_vec == curr_fold]
  
  # Model offset for current CV iteration
  cv_offset <- curr_offset[curr_fold_vec == curr_fold]
  
  
  # Linear predictor component of the predicted values
  cv_pred_lp <- {cv_matrix %*% curr_cv_k$glmnet.fit$beta[,curr_cv_k$glmnet.fit$lambda == curr_cv_k$lambda.min, drop = FALSE]}[,1]
  
  
  # Intercept component of the predicted values
  cv_pred_int <- curr_cv_k$glmnet.fit$a0[curr_cv_k$glmnet.fit$lambda == curr_cv_k$lambda.min]
  
  # The final predicted cumulative hazards for each segment of the piecewise hazard
  cv_chz <- exp(cv_pred_int + cv_pred_lp + cv_offset)
  
  
  # creating a dataset with multiple rows per PX_ID
  cv_df <- data.frame(PX_ID = curr_ids[curr_fold_vec == curr_fold],
                      cv_exp = cv_chz,
                      cv_deviance = 2 * (ifelse(cv_outcome == 0, 0, log(1 / cv_chz)) - (cv_outcome - cv_chz)))
  
  # condensing the data set into a single row for each px_id
  cv_df <- cv_df %>%
    group_by(PX_ID) %>%
    summarise(cv_exp = sum(cv_exp),
              cv_deviance = sum(cv_deviance))
  
  
  rm(list = c("cv_matrix", "cv_outcome", "cv_offset", "cv_pred_lp", "cv_pred_int", "cv_chz"))
  gc()

  ###################################################################################################################################################
  ###################################################################################################################################################
  ###################################################################################################################################################
  ###################################################################################################################################################
  
  ### Estimating the cross-validated probability of survival 1-year after start of follow-up
  
  ### Constructing a data set to predict the probabilities for PX_ID not included in the CV
  cv_1yr_data <- curr_imputed %>%
    subset(PX_ID %in% curr_ids[curr_fold_vec == curr_fold]) %>%
    mutate(right_time = left_time + 1 * 365)
  
  
  ### Creating new followup variables for the modified right_time variable
  cv_update <- create_followup_vars(cv_1yr_data,
                                    follow_up_cuts)
  
  cv_1yr_data <- cv_update$new_data
  
  
  ### Reshaping the cv_1yr_data
  cv_1yr_long <- reshape(cv_1yr_data,
                         varying = list(paste0("days_at_risk_",
                                               follow_up_cuts,
                                               "_",
                                               c(follow_up_cuts[-1], "")),
                                        paste0("event_",
                                               follow_up_cuts,
                                               "_",
                                               c(follow_up_cuts[-1], ""))),
                         v.names = c("days_at_risk",
                                     "event"),
                         times = paste0(follow_up_cuts,
                                        "_",
                                        c(follow_up_cuts[-1], "")),
                         direction = "long")
  
  
  # Cleaning up the workspace
  rm(list = c("cv_1yr_data", "cv_update"))
  gc()
  
  
  cv_1yr_long <- cv_1yr_long %>%
    subset(days_at_risk > 0) %>%
    mutate(intercept = 1)
  
  
  ### If the model is piecewise_varying, create separate variables for each piece of the baseline hazard
  if(piecewise_type == "piecewise_varying"){
    
    ### Creating separate variables for each piece of the baseline hazard
    tmp_Data <- create_timevarying_vars(curr_data = cv_1yr_long,
                                        curr_vars = curr_vars,
                                        bl_cuts = follow_up_cuts)
    
    
    # Updating the data set
    cv_1yr_long <- tmp_Data$new_data
    
    
    # Updating the variables for the model
    curr_vars <- tmp_Data$new_vars
    
    
    # Cleaning up the workspacke
    rm(list = "tmp_Data")
    gc()
  }
  
  
  # Model matrix for LASSO fitting
  cv_iter_model_matrix <- model.matrix(as.formula(curr_model_statement),
                                       data = cv_1yr_long)
  
  # Model offset vector for LASSO fitting
  cv_iter_1yr_offset <- log(cv_1yr_long$days_at_risk)
  
  # Linear predictor component of the predicted values
  cv_iter_1yr_lp <- {cv_iter_model_matrix %*% curr_cv_k$glmnet.fit$beta[,curr_cv_k$glmnet.fit$lambda == curr_cv_k$lambda.min]}[,1]
  
  # Intercept component of the predicted values
  cv_iter_1yr_int <- curr_cv_k$glmnet.fit$a0[curr_cv_k$glmnet.fit$lambda == curr_cv_k$lambda.min]
  
  # The final predicted cumulative hazards for each segment of the piecewise hazard
  cv_iter_1yr_chz <- exp(cv_iter_1yr_int + cv_iter_1yr_lp + cv_iter_1yr_offset)
  
  # The final cross-validated linear predictor for the c-statistic
  cv_iter_1yr <- data.frame(PX_ID = cv_1yr_long$PX_ID,
                            chz_1yr = cv_iter_1yr_chz) %>%
    group_by(PX_ID) %>%
    summarise(surv_1yr = exp(-sum(chz_1yr)))
  
  ###################################################################################################################################################
  ###################################################################################################################################################
  ###################################################################################################################################################
  ###################################################################################################################################################
  
  ### Estimating the cross-validated probability of survival 2-years after start of follow-up
  
  ### Constructing a data set to predict the probabilities for PX_ID not included in the CV
  cv_2yr_data <- curr_imputed %>%
    subset(PX_ID %in% curr_ids[curr_fold_vec == curr_fold]) %>%
    mutate(right_time = left_time + 2 * 365)
  
  
  ### Creating new followup variables for the modified right_time variable
  cv_update <- create_followup_vars(cv_2yr_data,
                                    follow_up_cuts)
  
  cv_2yr_data <- cv_update$new_data
  
  
  ### Reshaping the cv_2yr_data
  cv_2yr_long <- reshape(cv_2yr_data,
                       varying = list(paste0("days_at_risk_",
                                             follow_up_cuts,
                                             "_",
                                             c(follow_up_cuts[-1], "")),
                                      paste0("event_",
                                             follow_up_cuts,
                                             "_",
                                             c(follow_up_cuts[-1], ""))),
                       v.names = c("days_at_risk",
                                   "event"),
                       times = paste0(follow_up_cuts,
                                      "_",
                                      c(follow_up_cuts[-1], "")),
                       direction = "long")
  
  
  # Cleaning up the workspace
  rm(list = c("cv_2yr_data", "cv_update"))
  gc()
  
  
  cv_2yr_long <- cv_2yr_long %>%
    subset(days_at_risk > 0) %>%
    mutate(intercept = 1)
  
  
  ### If the model is piecewise_varying, create separate variables for each piece of the baseline hazard
  if(piecewise_type == "piecewise_varying"){
    
    ### Creating separate variables for each piece of the baseline hazard
    tmp_Data <- create_timevarying_vars(curr_data = cv_2yr_long,
                                        curr_vars = curr_vars,
                                        bl_cuts = follow_up_cuts)
    
    
    # Updating the data set
    cv_2yr_long <- tmp_Data$new_data
    
    
    # Updating the variables for the model
    curr_vars <- tmp_Data$new_vars
    
    
    # Cleaning up the workspacke
    rm(list = "tmp_Data")
    gc()
  }
  
  
  # Model matrix for LASSO fitting
  cv_iter_model_matrix <- model.matrix(as.formula(curr_model_statement),
                                       data = cv_2yr_long)
  
  # Model offset vector for LASSO fitting
  cv_iter_2yr_offset <- log(cv_2yr_long$days_at_risk)
  
  # Linear predictor component of the predicted values
  cv_iter_2yr_lp <- {cv_iter_model_matrix %*% curr_cv_k$glmnet.fit$beta[,curr_cv_k$glmnet.fit$lambda == curr_cv_k$lambda.min]}[,1]
  
  # Intercept component of the predicted values
  cv_iter_2yr_int <- curr_cv_k$glmnet.fit$a0[curr_cv_k$glmnet.fit$lambda == curr_cv_k$lambda.min]
  
  # The final predicted cumulative hazards for each segment of the piecewise hazard
  cv_iter_2yr_chz <- exp(cv_iter_2yr_int + cv_iter_2yr_lp + cv_iter_2yr_offset)
  
  # The final cross-validated linear predictor for the c-statistic
  cv_iter_2yr <- data.frame(PX_ID = cv_2yr_long$PX_ID,
                            chz_2yr = cv_iter_2yr_chz) %>%
    group_by(PX_ID) %>%
    summarise(surv_2yr = exp(-sum(chz_2yr)))
  
  
  ### creating the final data set for the cross-validated data
  cv_fnl <- merge(cv_df,
                  cv_iter_1yr,
                  by = "PX_ID") %>%
    merge(cv_iter_2yr,
          by = "PX_ID")
  
  
  if(any(c(nrow(cv_df), nrow(cv_iter_1yr), nrow(cv_iter_2yr)) != nrow(cv_fnl))){ stop("merge of CV datasets failed.") }
  
  
  ### Writing data for keeping track of the predicted values
  write.table(cv_fnl,
              file = file.path(root,
                               organ,
                               piecewise_type,
                               paste0(age_type, "_cv_prob_", mi_i, ".csv")),
              col.names = FALSE,
              row.names = FALSE,
              append = TRUE,
              sep = ",")
  
  
  ### Cleaning the workspace
  rm(list = c("curr_cv_k", "cv_iter_model_matrix", "cv_iter_offset", "cv_iter_pred_lp", "cv_df", 
              "cv_iter_pred_int", "cv_iter_chz", "cv_iter_pred", "cv_pred_long"))
  gc()
  
  
  return(NULL)
}







