###
### Purpose: Estimate a survival from listing model for the ith iteration of multiple imputation
###
### Author: Andrew Wey
###


estimate_model_i <- function(root,
                             organ,
                             organ_ab,
                             age_type,
                             mi_i,
                             follow_up_cuts,
                             curr_mi_covariates,
                             curr_mi_cov_method,
                             curr_control,
                             model_type){
  
  ### root: The base directory for the model building
  ### organ: The organ name (e.g., Kidney) for model building
  ### organ_ab: The organ abbreviation (e.g., ki) for model building
  ### age_type: Identifies whether the model is for pediatric or adult candidates
  ### mi_i: Identifies the iteration of the multiple imputation
  ### follow_up_cuts: The cutpoints for splitting the baseline hazard and creating time-varying effects
  ### curr_mi_covariates: The covariates included in the current model
  ### curr_mi_cov_method: The covariate type (categorical/continuous) for covariates included in the current model
  ### curr_control: The control table for the current model (organ-specific)
  ### model_type: Specifying the type of model to fit
  
  
  ### Checks to make sure the function arguments are reasonable
  if({! age_type %in% c("peds", "adult")} & {organ_ab != "pa"}){
    stop("Incorrect argument to age_type")
  }
  if({! age_type %in% c("pa", "kp")} & {organ_ab == "pa"}){
    stop("Incorrect argument to age_type")
  }
  if(! model_type %in% c("piecewise_constant", "piecewise_varying", "piecewise_separate", "coxph_constant", "gam_pem")){
    stop("Incorrect argument to model_type")
  }
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Loading the raw data; used to select the appropriate reference levels for categorical variables
  load(file.path(root,
                 organ,
                 paste0("sfl_data_",
                        age_type,
                        ".Rdata")))
  
  
  ### Assign the data set to a generic name
  assign("raw_data",
         get(paste0("sfl_data_",
                    age_type)))
  
  
  ### Removing the original raw data set from the workspace
  rm(list = paste0("sfl_data_", age_type))
  gc()
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Loading the appropriate iteration of MI data
  load(file.path(root, 
                 organ, 
                 "MI_Output", 
                 paste0("mi_", organ_ab, "_", age_type, "_data_", mi_i, ".Rdata")))
  
  
  curr_data <- mi_data
  
  
  ### The MI changes the event status to a factor. In general, the following code expects event status as numeric.  So we switch it back.
  curr_data$event <- as.numeric(as.character(curr_data$event))

  
  rm(list = "mi_data")
  gc()
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Additional cleaning of binary categorical variables for easier use in LASSO fitting 
  
  dummy <- create_dummy_variables(curr_data,
                                  raw_data,
                                  curr_mi_covariates,
                                  curr_mi_cov_method,
                                  curr_control)
  
  # Updated data
  curr_imputed <- dummy$new_data
  
  # Updated variables
  curr_vars <- dummy$new_vars
  
  
  # If model_type is coxph_constant, then we need to track the type of the variables
  if(model_type %in% c("coxph_constant", "gam_pem")){
    
    curr_vars_type <- rep("categorical", 
                          length(curr_vars))
  }
  
  
  rm(list = c("dummy", "curr_data"))
  gc()
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Creating linear splines for the continuous variables; only necessary for non-coxph versions
  
  if(model_type %in% c("piecewise_constant", "piecewise_varying", "piecewise_separate")){
    
    # for loop to create the splines for each continuous variable
    for(k in 1:length(curr_mi_covariates)){
      
      covar <- curr_mi_covariates[k]
      
      if(curr_mi_cov_method[k] == "continuous"){
        
        curr_spline <- create_linear_splines(curr_imputed,
                                             covar)
        
        
        curr_imputed <- curr_spline$new_data
        
        
        curr_vars <- c(curr_vars,
                       curr_spline$new_vars)
        
        rm(list = "curr_spline")
        gc()
      }
    }
  }
  
  # If the model_type is coxph_constant, then the continuous variables need to enter the curr_vars vector
  if(model_type %in% c("coxph_constant", "gam_pem")){
    
    # for loop to create the splines for each continuous variable
    for(k in 1:length(curr_mi_covariates)){
      
      covar <- curr_mi_covariates[k]
      
      if(curr_mi_cov_method[k] == "continuous"){
        
        curr_vars <- c(curr_vars,
                       covar)
        
        curr_vars_type <- c(curr_vars_type,
                            "continuous")
      }
    }
  }

  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Creating the multiple followup variables and then transforming to a long format; only necessary for non-coxph versions
  
  if(model_type %in% c("piecewise_constant", "piecewise_varying", "piecewise_separate", "gam_pem")){
    
        
    curr_update <- create_followup_vars(curr_imputed,
                                        follow_up_cuts)
    
    
    curr_imputed <- curr_update$new_data
    
    
    rm(list = c("curr_update"))
    gc()
  }

  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Reshaping to a long format. Only necessary for the piecewise constant hazard versions
  
  if(model_type %in% c("piecewise_constant", "piecewise_varying", "piecewise_separate", "gam_pem")){
    
    model_data <- reshape(curr_imputed,
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
    
    
    ### The intercept is used to create the different components of the piecewise baseline hazard
    model_data <- model_data %>%
      subset(days_at_risk > 0) %>%
      mutate(intercept = 1)
    
    
    
    ### Creating a vector that has the base covariates (before adding the time-varying components); necessary for the cross-validation later
    static_vars <- curr_vars
    
    
    ### If the model is piecewise_varying, create separate variables for each piece of the baseline hazard
    if(model_type %in% c("piecewise_varying", "gam_pem")){
      
      ### Creating separate variables for each piece of the baseline hazard
      tmp_Data <- create_timevarying_vars(curr_data = model_data,
                                          curr_vars = curr_vars,
                                          curr_vars_type = curr_vars_type,
                                          bl_cuts = follow_up_cuts)
      
      
      # Updating the data set
      model_data <- tmp_Data$new_data
      
      
      # Updating the variables for the model
      curr_vars <- tmp_Data$new_vars
      
      
      # Updating the variable types for the time-varying effect
      curr_vars_type <- tmp_Data$new_vars_type
      
      
      # Cleaning up the workspacke
      rm(list = "tmp_Data")
      gc()
    }
  }
  
  
  ### Standardizing the data set name for the cox proportional hazards model
  if(model_type %in% c("coxph_constant")){
    
    model_data <- curr_imputed
    
    ### Cleaning the work space for potential memory issues 
    rm(list = c("curr_imputed"))
    gc()
  }
  
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Fitting the gam_pem model
  
  if(model_type %in% c("gam_pem")){
    
    ### Fitting the cross-validation models: We need to complete cross-validation and
    ###   the rows for a given candidate need to be kept together
    # Create indicators for appropriate folds
    uniq_ids <- curr_imputed$PX_ID
    
    
    # Checking for duplicated PX_IDs
    if(any(duplicated(uniq_ids))){ stop("Duplicated PX_IDs") }
    
    
    if((length(uniq_ids) %/% 10) > 0){
      
      fold_ind <- rep(1:10,
                      length(uniq_ids) %/% 10)
      
      if((length(uniq_ids) %% 10) > 0){
        fold_ind <- c(fold_ind,
                      1:{length(uniq_ids) %% 10})
      }
    }
    if((length(uniq_ids) %/% 10) == 0){
      
      fold_ind <- 1:{length(uniq_ids) %% 10}
    }
    
    # sample the fold indicators with replacement (seed is set to current year)
    set.seed(as.numeric(substr(Sys.time(), 1, 4)))
    
    fold_vec <- sample(fold_ind,
                       length(fold_ind),
                       replace = FALSE)
    
    
    # extending fold_vec to match model_data
    model_fold_vec <- fold_vec[match(model_data$PX_ID, uniq_ids)]
    
    
    # Clean workspace
    rm(list = c("fold_ind", "fold_vec"))
    gc()
    
    
    ### gam_pem is not penalized. So we do not want to include effects for the first interval. Otherwise, the model is over-parameterized
    first_interval <- paste0(follow_up_cuts,
                             "_",
                             c(follow_up_cuts[-1], ""))[1]
    
    
    model_statement <- paste0("event ~ ",
                              paste(paste0(ifelse(curr_vars_type[! grepl(paste0("fup_", first_interval), curr_vars)] == "categorical",
                                                  "(",
                                                  "("),
                                           curr_vars[! grepl(paste0("fup_", first_interval), curr_vars)],
                                           ")"),
                                    collapse = " + "),
                              " + offset(log(days_at_risk))")
    
    ### Fit the overall model
    start_time <- Sys.time()
    all_mi_i <- glm(as.formula(model_statement),
                    data = model_data,
                    model = FALSE,
                    family = poisson)
    end_time <- Sys.time()
    
    
    # If it does not exist, create the directory for the model fit
    if(! dir.exists(file.path(root, organ, model_type))){
      dir.create(file.path(root, organ, model_type))
    }
    
    
    # Creating an object with the model coefficients
    all_mi_i_coef <- coef(all_mi_i)
    
    
    # Save the model output
    save(all_mi_i_coef,
         file = file.path(root, 
                          organ,
                          model_type,
                          paste0(age_type, "_all_fit_", mi_i, ".Rdata")))    
    
    rm(list = 'all_mi_i')
    gc()
    
    
    ### Running the cross-validation
    for(k in 1:10){
      
      ###################################################################################################################################################
      ###################################################################################################################################################

      ### Fitting the model for cross-validation

      ### Time to fit the actual model(!)
      curr_cv_k <- glm(as.formula(model_statement),
                       data = model_data,
                       family = poisson,
                       subset = model_fold_vec != k)

      ###################################################################################################################################################
      ###################################################################################################################################################
      
      ### Estimating the cross-validated number of expected events and deviance
      
      # The predicted cumulative hazards for each segment of the piecewise hazard
      cv_chz <- predict(curr_cv_k,
                        newdata = model_data[model_fold_vec == k,],
                        type = "response")
      
      
      # creating a dataset with multiple rows per PX_ID
      cv_df <- data.frame(PX_ID = model_data$PX_ID[model_fold_vec == k],
                          cv_exp = cv_chz,
                          cv_deviance = 2 * (ifelse(model_data$event[model_fold_vec == k] == 0, 
                                                    0, 
                                                    log(1 / cv_chz)) - 
                                               (model_data$event[model_fold_vec == k] - cv_chz)))
      
      
      # condensing the data set into a single row for each px_id
      cv_df <- cv_df %>%
        group_by(PX_ID) %>%
        summarise(cv_exp = sum(cv_exp),
                  cv_deviance = sum(cv_deviance))
      
      
      rm(list = c("cv_chz"))
      gc()
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      ### Estimating the cross-validated probability of survival 1-year after start of follow-up
      
      ### Constructing a data set to predict the probabilities for PX_ID not included in the CV
      cv_1yr_data <- curr_imputed %>%
        subset(PX_ID %in% model_data$PX_ID[model_fold_vec == k]) %>%
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
      
      
      ### Creating separate variables for each piece of the baseline hazard
      tmp_Data <- create_timevarying_vars(curr_data = cv_1yr_long,
                                          curr_vars = static_vars,
                                          bl_cuts = follow_up_cuts)
      
      
      # Updating the data set
      cv_1yr_long <- tmp_Data$new_data

      
      # Cleaning up the workspacke
      rm(list = "tmp_Data")
      gc()
      
      
      # The final predicted cumulative hazards for each segment of the piecewise hazard
      cv_iter_1yr_chz <- predict(curr_cv_k,
                                 newdata = cv_1yr_long,
                                 type = "response")
      
      
      # The final cross-validated linear predictor for the C-statistic
      cv_iter_1yr <- data.frame(PX_ID = cv_1yr_long$PX_ID,
                                chz_1yr = cv_iter_1yr_chz) %>%
        group_by(PX_ID) %>%
        summarise(surv_1yr = exp(-sum(chz_1yr)))
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      ### Estimating the cross-validated probability of survival 2-year after start of follow-up
      
      ### Constructing a data set to predict the probabilities for PX_ID not included in the CV
      cv_2yr_data <- curr_imputed %>%
        subset(PX_ID %in% model_data$PX_ID[model_fold_vec == k]) %>%
        mutate(right_time = left_time + 2 * 365)
      
      
      ### Creating new followup variables for the modified right_time variable
      cv_update <- create_followup_vars(cv_2yr_data,
                                        follow_up_cuts)
      
      cv_2yr_data <- cv_update$new_data
      
      
      ### Reshaping the cv_1yr_data
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
      
      
      ### Creating separate variables for each piece of the baseline hazard
      tmp_Data <- create_timevarying_vars(curr_data = cv_2yr_long,
                                          curr_vars = static_vars,
                                          bl_cuts = follow_up_cuts)
      
      
      # Updating the data set
      cv_2yr_long <- tmp_Data$new_data
      
      
      # Cleaning up the workspacke
      rm(list = "tmp_Data")
      gc()
      
      
      # The final predicted cumulative hazards for each segment of the piecewise hazard
      cv_iter_2yr_chz <- predict(curr_cv_k,
                                 newdata = cv_2yr_long,
                                 type = "response")
      
      
      # The final cross-validated linear predictor for the C-statistic
      cv_iter_2yr <- data.frame(PX_ID = cv_2yr_long$PX_ID,
                                chz_2yr = cv_iter_2yr_chz) %>%
        group_by(PX_ID) %>%
        summarise(surv_2yr = exp(-sum(chz_2yr)))
      
      
      ### Creating the final data set for the cross-validated data
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
                                   model_type,
                                   paste0(age_type, "_cv_prob_", mi_i, ".csv")),
                  col.names = FALSE,
                  row.names = FALSE,
                  append = TRUE,
                  sep = ",")
    }
  }
  
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Fitting the piecewise_constant or piecewise_varying model
  
  if(model_type %in% c("piecewise_constant", "piecewise_varying")){
    
    
    ### Fitting the cross-validation models: We need to complete cross-validation and
    ###   the rows for a given candidate need to be kept together
    # Create indicators for appropriate folds
    uniq_ids <- curr_imputed$PX_ID
    
    
    # Checking for duplicated PX_IDs
    if(any(duplicated(uniq_ids))){ stop("Duplicated PX_IDs") }
    
    
    if((length(uniq_ids) %/% 10) > 0){
      
      fold_ind <- rep(1:10,
                      length(uniq_ids) %/% 10)
      
      if((length(uniq_ids) %% 10) > 0){
        fold_ind <- c(fold_ind,
                      1:{length(uniq_ids) %% 10})
      }
    }
    if((length(uniq_ids) %/% 10) == 0){
      
      fold_ind <- 1:{length(uniq_ids) %% 10}
    }
    
    # sample the fold indicators with replacement (seed is set to current year)
    set.seed(as.numeric(substr(Sys.time(), 1, 4)))
    
    fold_vec <- sample(fold_ind,
                       length(fold_ind),
                       replace = FALSE)
    
    
    # extending fold_vec to match model_data
    model_fold_vec <- fold_vec[match(model_data$PX_ID, uniq_ids)]
    
    
    # Clean workspace
    rm(list = c("fold_ind", "fold_vec"))
    gc()
    
    
    
    ### The model fitting component of the code
    model_statement <- paste("~ -1 + ",
                             paste(curr_vars, 
                                   collapse = " + "))
    
    
    # Model matrix for LASSO fitting
    lasso_model_matrix <- model.matrix(as.formula(model_statement),
                                       data = model_data)
    
    
    # Model outcome vector for LASSO fitting
    lasso_outcome <- model_data$event
    
    # Model offset vector for LASSO fitting
    lasso_offset <- log(model_data$days_at_risk)
    
    # The PX_IDs
    lasso_pxs <- model_data$PX_ID
    
    # Cleaning workspace
    rm(list = "model_data")
    gc()
    
    
    
    ### Time to fit the actual model(!)
    start_time <- Sys.time()
    all_mi_i <- cv.glmnet(x = lasso_model_matrix,
                          y = lasso_outcome,
                          family = "poisson",
                          offset = lasso_offset,
                          nfolds = 10,
                          foldid = model_fold_vec)
    end_time <- Sys.time()
    
    
    # If it does not exist, create the directory for the model fit
    if(! dir.exists(file.path(root, organ, model_type))){
      dir.create(file.path(root, organ, model_type))
    }
    
    
    # Save the model output
    save(all_mi_i,
         file = file.path(root, 
                          organ,
                          model_type,
                          paste0(age_type, "_all_fit_", mi_i, ".Rdata")))
    
    
    rm(list = "all_mi_i")
    gc()
    
    
    # Looping through the CV folds
    for(k in 1:10){
      piecewise_cv_fold(curr_imputed,
                        static_vars,
                        model_statement, 
                        lasso_model_matrix,
                        lasso_outcome,
                        lasso_offset,
                        lasso_pxs,
                        model_fold_vec,
                        k,
                        model_type,
                        root,
                        organ,
                        age_type,
                        mi_i)
    }
  }

  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Fitting the piecewise_separate model
  
  if(model_type == "piecewise_separate"){

    ### Fitting the cross-validation models: We need to complete cross-validation and
    ###   the rows for a given candidate need to be kept together
    # Create indicators for appropriate folds
    uniq_ids <- curr_imputed$PX_ID
    
    
    # Checking for duplicated PX_IDs
    if(any(duplicated(uniq_ids))){ stop("Duplicated PX_IDs") }
    
    
    if((length(uniq_ids) %/% 10) > 0){
      
      fold_ind <- rep(1:10,
                      length(uniq_ids) %/% 10)
      
      if((length(uniq_ids) %% 10) > 0){
        fold_ind <- c(fold_ind,
                      1:{length(uniq_ids) %% 10})
      }
    }
    if((length(uniq_ids) %/% 10) == 0){
      
      fold_ind <- 1:{length(uniq_ids) %% 10}
    }
    
    # sample the fold indicators with replacement (seed is set to current year)
    set.seed(as.numeric(substr(Sys.time(), 1, 4)))
    
    fold_vec <- sample(fold_ind,
                       length(fold_ind),
                       replace = FALSE)
    
    
    # extending fold_vec to match model_data
    model_fold_vec <- fold_vec[match(model_data$PX_ID, uniq_ids)]
    
    
    # Clean workspace
    rm(list = c("fold_ind"))
    gc()
    
    
    ### The piecewise_separate model fits a separate LASSO model for each of the follow-up periods
    intervals <- paste0(follow_up_cuts, "_", c(follow_up_cuts[-1], ""))
    
    
    # The model statement does not depend on the given interval
    model_statement <- paste("~ -1 + ",
                             paste(curr_vars, 
                                   collapse = " + "))
    
    
    # loop through each interval
    start_time <- Sys.time()
    for(r in intervals){
      
      # Subsetting model data to the current interval
      model_data_curr_interval <- model_data %>%
        subset(time == r)
      
      
      # Model matrix for LASSO fitting
      lasso_model_matrix <- model.matrix(as.formula(model_statement),
                                         data = model_data_curr_interval)
      
      # Model outcome vector for LASSO fitting
      lasso_outcome <- model_data_curr_interval$event
      
      # Model offset vector for LASSO fitting
      lasso_offset <- log(model_data_curr_interval$days_at_risk)
      
      # The PX_IDs
      lasso_pxs <- model_data_curr_interval$PX_ID
      
      # The CV vector for the current interval
      lasso_folds <- model_fold_vec[model_data$time == r]
      
      
      # Cleaning workspace
      rm(list = "model_data_curr_interval")
      gc()
      
      
      ### Time to fit the actual model(!)
      all_mi_i <- cv.glmnet(x = lasso_model_matrix,
                            y = lasso_outcome,
                            family = "poisson",
                            offset = lasso_offset,
                            nfolds = 10,
                            foldid = lasso_folds)
      
      
      # If it does not exist, create the directory for the model fit
      if(! dir.exists(file.path(root, organ, model_type))){
        dir.create(file.path(root, organ, model_type))
      }
      
      
      # Save the model output
      save(all_mi_i,
           file = file.path(root, 
                            organ,
                            model_type,
                            paste0(age_type, "_fit_", r, "_", mi_i, ".Rdata")))
      
      
      rm(list = "all_mi_i")
      gc()
    }
    end_time <- Sys.time()
    
    #########################################################################################################
    #########################################################################################################
    
    ### looping through 10 cross-validation goals
    
    for(k in 1:10){
      
      ### Constructing a data set to predict the survival 1 year after the start of follow-up
      cv_1yr_data <- curr_imputed %>%
        subset(PX_ID %in% uniq_ids[fold_vec == k]) %>%
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
      
      ################################################################################################################################################
      ################################################################################################################################################
      
      ### Constructing a data set to predict the survival 2 years after the start of follow-up
      cv_2yr_data <- curr_imputed %>%
        subset(PX_ID %in% uniq_ids[fold_vec == k]) %>%
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
      
      
      
      ### loop through each interval: need to track the predicted for px_id and interval combination
      # The data frame for tracking the output for each interval
      cv_df_long <- NULL
      cv_iter_1yr_long <- NULL
      cv_iter_2yr_long <- NULL
      
      for(r in intervals){
        
        # Subsetting model data to the current interval
        model_data_curr_interval <- model_data %>%
          subset({PX_ID %in% uniq_ids[fold_vec != k]} & 
                   {time == r})
        
        
        # Model matrix for LASSO fitting
        lasso_model_matrix <- model.matrix(as.formula(model_statement),
                                           data = model_data_curr_interval)
        
        # Model outcome vector for LASSO fitting
        lasso_outcome <- model_data_curr_interval$event
        
        # Model offset vector for LASSO fitting
        lasso_offset <- log(model_data_curr_interval$days_at_risk)
        
        # The PX_IDs
        lasso_pxs <- model_data_curr_interval$PX_ID
        
        # Cleaning workspace
        rm(list = "model_data_curr_interval")
        gc()
        
        
        ### Time to fit the actual model(!)
        curr_interval_r_cv_k <- cv.glmnet(x = lasso_model_matrix,
                                          y = lasso_outcome,
                                          family = "poisson",
                                          offset = lasso_offset)
        
        #################################################################################################################################
        #################################################################################################################################
        
        ### Estimating the cross-validated number of expected events and deviance
        
        # Subsetting model data to the current interval
        model_data_curr_interval <- model_data %>%
          subset({PX_ID %in% uniq_ids[fold_vec == k]} & 
                   {time == r})
        
        
        # Model matrix for LASSO fitting
        cv_matrix <- model.matrix(as.formula(model_statement),
                                  data = model_data_curr_interval)
        
        # Model outcome vector for LASSO fitting
        cv_outcome <- model_data_curr_interval$event
        
        # Model offset vector for LASSO fitting
        cv_offset <- log(model_data_curr_interval$days_at_risk)
        
        
        # Linear predictor component of the predicted values
        cv_pred_lp <- {cv_matrix %*% curr_interval_r_cv_k$glmnet.fit$beta[,curr_interval_r_cv_k$glmnet.fit$lambda == curr_interval_r_cv_k$lambda.min]}[,1]
        
        # Intercept component of the predicted values
        cv_pred_int <- curr_interval_r_cv_k$glmnet.fit$a0[curr_interval_r_cv_k$glmnet.fit$lambda == curr_interval_r_cv_k$lambda.min]
        
        # The final predicted cumulative hazards for each segment of the piecewise hazard
        cv_chz <- exp(cv_pred_int + cv_pred_lp + cv_offset)
        
        
        # creating a dataset with multiple rows per PX_ID
        cv_df <- data.frame(PX_ID = model_data_curr_interval$PX_ID,
                            cv_exp = cv_chz,
                            cv_deviance = 2 * (ifelse(cv_outcome == 0, 0, log(1 / cv_chz)) - (cv_outcome - cv_chz)))
        
        # condensing the data set into a single row for each px_id
        cv_df <- cv_df %>%
          group_by(PX_ID) %>%
          summarise(cv_exp = sum(cv_exp),
                    cv_deviance = sum(cv_deviance))
        
        
        cv_df_long <- rbind(cv_df_long,
                            cv_df)
        
        
        rm(list = c("cv_matrix", "cv_outcome", "cv_offset", "cv_pred_lp", "cv_pred_int", "cv_chz", "cv_df"))
        gc()
        
        #################################################################################################################################
        #################################################################################################################################
        
        ### Selecting the rows of the cv_clong that correspond to the current follow-up interval
        curr_cv_1yr_long <- cv_1yr_long %>%
          subset(time == r)
        
        
        # Model matrix for LASSO fitting
        cv_iter_model_matrix <- model.matrix(as.formula(model_statement),
                                             data = curr_cv_1yr_long)
        
        # Model offset vector for LASSO fitting
        cv_iter_1yr_offset <- log(curr_cv_1yr_long$days_at_risk)
        
        # Linear predictor component of the predicted values
        cv_iter_1yr_lp <- {cv_iter_model_matrix %*% curr_interval_r_cv_k$glmnet.fit$beta[,curr_interval_r_cv_k$glmnet.fit$lambda == curr_interval_r_cv_k$lambda.min]}[,1]
        
        # Intercept component of the predicted values
        cv_iter_1yr_int <- curr_interval_r_cv_k$glmnet.fit$a0[curr_interval_r_cv_k$glmnet.fit$lambda == curr_interval_r_cv_k$lambda.min]
        
        # The final predicted cumulative hazards for each segment of the piecewise hazard
        cv_iter_1yr_chz <- exp(cv_iter_1yr_int + cv_iter_1yr_lp + cv_iter_1yr_offset)
        
        # The final predicted linear predictor for a given PX_ID and follow-up
        cv_iter_1yr <- data.frame(PX_ID = curr_cv_1yr_long$PX_ID,
                                  chz_1yr = cv_iter_1yr_chz)
        
        # binding cv_iter_c to cv_iter_2yr_long; the data frame that will have multiple rows per px_id
        cv_iter_1yr_long <- rbind(cv_iter_1yr_long,
                                  cv_iter_1yr)
        
        ### Cleaning the workspace
        rm(list = c("cv_iter_1yr"))
        gc()
        
        #################################################################################################################################
        #################################################################################################################################
        
        ### Selecting the rows of the cv_clong that correspond to the current follow-up interval
        curr_cv_2yr_long <- cv_2yr_long %>%
          subset(time == r)
        
        
        # Model matrix for LASSO fitting
        cv_iter_model_matrix <- model.matrix(as.formula(model_statement),
                                             data = curr_cv_2yr_long)
        
        # Model offset vector for LASSO fitting
        cv_iter_2yr_offset <- log(curr_cv_2yr_long$days_at_risk)
        
        # Linear predictor component of the predicted values
        cv_iter_2yr_lp <- {cv_iter_model_matrix %*% curr_interval_r_cv_k$glmnet.fit$beta[,curr_interval_r_cv_k$glmnet.fit$lambda == curr_interval_r_cv_k$lambda.min]}[,1]
        
        # Intercept component of the predicted values
        cv_iter_2yr_int <- curr_interval_r_cv_k$glmnet.fit$a0[curr_interval_r_cv_k$glmnet.fit$lambda == curr_interval_r_cv_k$lambda.min]
        
        # The final predicted cumulative hazards for each segment of the piecewise hazard
        cv_iter_2yr_chz <- exp(cv_iter_2yr_int + cv_iter_2yr_lp + cv_iter_2yr_offset)
        
        # The final predicted linear predictor for a given PX_ID and follow-up
        cv_iter_2yr <- data.frame(PX_ID = curr_cv_2yr_long$PX_ID,
                                  chz_2yr = cv_iter_2yr_chz)
        
        # binding cv_iter_c to cv_iter_2yr_long; the data frame that will have multiple rows per px_id
        cv_iter_2yr_long <- rbind(cv_iter_2yr_long,
                                  cv_iter_2yr)
        
        ### Cleaning the workspace
        rm(list = c("curr_interval_r_cv_k", "cv_iter_2yr"))
        gc()
      }
      
      
      cv_iter_1yr <- cv_iter_1yr_long  %>%
        group_by(PX_ID) %>%
        summarise(surv_1yr = exp(-sum(chz_1yr)))
      
      
      cv_iter_2yr <- cv_iter_2yr_long  %>%
        group_by(PX_ID) %>%
        summarise(surv_2yr = exp(-sum(chz_2yr)))
      
      
      cv_df <- cv_df_long %>%
        group_by(PX_ID) %>%
        summarise(cv_exp = sum(cv_exp),
                  cv_deviance = sum(cv_deviance))
      
      
      
      ### creating the final data set for the cross-validated data
      cv_fnl <- merge(cv_df,
                      cv_iter_1yr,
                      by = "PX_ID")%>%
        merge(cv_iter_2yr,
              by = "PX_ID")
      
      
      if(any(c(nrow(cv_df), nrow(cv_iter_1yr), nrow(cv_iter_2yr)) != nrow(cv_fnl))){ stop("merge of CV datasets failed.") }
      
      
      ### Writing data for keeping track of the predicted values
      write.table(cv_fnl,
                  file = file.path(root,
                                   organ,
                                   model_type,
                                   paste0(age_type, "_cv_prob_", mi_i, ".csv")),
                  col.names = FALSE,
                  row.names = FALSE,
                  append = TRUE,
                  sep = ",")
    }
  }
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  
  ### Fitting the coxph_constant model
  
  if(model_type %in% c("coxph_constant")){
    
    
    ### The model statement for the coxph depends on the variable type.  model_prep prepares each variable for entry into coxph
    model_prep <- ifelse(curr_vars_type == "continuous",
                         paste0("pspline(",
                                curr_vars,
                                ")"),
                         curr_vars)
    
    
    model_statement <- paste("Surv(left_time, I(right_time + 1), event) ~",
                             paste(model_prep,
                                   collapse = " + "))
    
    
    ### We need to complete cross-validation for coxph
    # Create indicators for appropriate folds
    if((nrow(model_data) %/% 10) > 0){
      
      fold_ind <- rep(1:10,
                      nrow(model_data) %/% 10)
      
      if((nrow(model_data) %% 10) > 0){
        fold_ind <- c(fold_ind,
                      1:{nrow(model_data) %% 10})
      }
    }
    if((nrow(model_data) %/% 10) == 0){
      
      fold_ind <- 1:{nrow(model_data) %% 10}
    }
    
    # sample the fold indicators with replacement (seed is set to current year)
    set.seed(as.numeric(substr(Sys.time(), 1, 4)))
    
    fold_vec <- sample(fold_ind,
                       length(fold_ind),
                       replace = FALSE)
    
    # Clean workspace
    rm(list = "fold_ind")
    gc()
    
    
    ### Fit the overall model
    start_time <- Sys.time()
    all_mi_i <- coxph(as.formula(model_statement),
                      data = model_data)
    end_time <- Sys.time()
    
    
    ### Creating figures to investigat potential non-proportional hazards
    create_coxzph_figs(root,
                       organ,
                       age_type,
                       all_mi_i,
                       model_data,
                       mi_i,
                       follow_up_cuts)
    
    
    # If it does not exist, create the directory for the model fit
    if(! dir.exists(file.path(root, organ, model_type))){
      dir.create(file.path(root, organ, model_type))
    }
    
    
    # Save the model output
    save(all_mi_i,
         file = file.path(root, 
                          organ,
                          model_type,
                          paste0(age_type, "_all_fit_", mi_i, ".Rdata")))    
    
    rm(list = 'all_mi_i')
    gc()
    
    
    ### Fitting the cross-validation models
    
    # Creating a vector to track the cross-validated estimated probabilities of survival 1 or 2 years after listing
    cv_exps  <- rep(NA, nrow(model_data))
    cv_1yr_probs <- rep(NA, nrow(model_data))
    cv_2yr_probs <- rep(NA, nrow(model_data))
    
    for(k in 1:10){
      
      cv_k <- coxph(as.formula(model_statement),
                    data = model_data,
                    subset = fold_vec != k)
      
      
      ### The expected for calculating deviance
      cv_k_exp <- predict(cv_k,
                          newdata = model_data[fold_vec == k, ],
                          type = "expected")
      
      
      ### Saving the expected for the deviance
      cv_exps[fold_vec == k]  <- cv_k_exp

      #####################################################################################################################
      #####################################################################################################################
      
      ### I need a modified right-time variable to calculate probability of survival 1-year after follow-up
      
      pred_data <- model_data[fold_vec == k, ]
      pred_data <- pred_data %>%
        mutate(right_time = left_time + 1 * 365)
      
      
      # The baseline cumulative hazard function
      bl_hazard <- basehaz(cv_k)
      
      
      # The linear predictors for individal
      cv_k_lps <- predict(cv_k,
                          newdata = pred_data,
                          type = "lp")
      
      
      # Calculating the expected events at the end of the cohort period
      cv_k_chf_left  <- sapply(pred_data$left_time, function(x){ 
        if(sum(bl_hazard$time <= x) == 0){ fnl <- 0 }
        if(sum(bl_hazard$time <= x) >  0){ fnl <- bl_hazard$hazard[sum(bl_hazard$time <= x)] }
        
        return(fnl)
      })
      
      cv_k_chf_right <- sapply(pred_data$right_time, function(x){ 
        if(sum(bl_hazard$time <= x) == 0){ fnl <- 0 }
        if(sum(bl_hazard$time <= x) >  0){ fnl <- bl_hazard$hazard[sum(bl_hazard$time <= x)] }
        
        return(fnl)
      })
      
      cv_k_e <- (cv_k_chf_right - cv_k_chf_left) * exp(cv_k_lps)
      
      
      ### Saving the cross-validated probability of survival 1-year after start of follow-up
      cv_1yr_probs[fold_vec == k] <- exp(-cv_k_e)
      
      #####################################################################################################################
      #####################################################################################################################
      
      ### I need a modified right-time variable to calculate probability of survival 2-years after follow-up
      
      pred_data <- model_data[fold_vec == k, ]
      pred_data <- pred_data %>%
        mutate(right_time = left_time + 2 * 365)
      
      
      # The baseline cumulative hazard function
      bl_hazard <- basehaz(cv_k)
      
      
      # The linear predictors for individal
      cv_k_lps <- predict(cv_k,
                          newdata = pred_data,
                          type = "lp")
      
      
      # Calculating the expected events at the end of the cohort period
      cv_k_chf_left  <- sapply(pred_data$left_time, function(x){ 
        if(sum(bl_hazard$time <= x) == 0){ fnl <- 0 }
        if(sum(bl_hazard$time <= x) >  0){ fnl <- bl_hazard$hazard[sum(bl_hazard$time <= x)] }
        
        return(fnl)
      })
      
      cv_k_chf_right <- sapply(pred_data$right_time, function(x){ 
        if(sum(bl_hazard$time <= x) == 0){ fnl <- 0 }
        if(sum(bl_hazard$time <= x) >  0){ fnl <- bl_hazard$hazard[sum(bl_hazard$time <= x)] }
        
        return(fnl)
      })
      
      cv_k_e <- (cv_k_chf_right - cv_k_chf_left) * exp(cv_k_lps)
      
      
      ### Saving the cross-validated probability of survival 2-years after start of follow-up
      cv_2yr_probs[fold_vec == k] <- exp(-cv_k_e)

      
      rm(list = "cv_k")
      gc()
    }
    
    

    write.table(data.frame(px_id = model_data$PX_ID,
                           cv_exp = cv_exps,
                           cv_1yr_prob = cv_1yr_probs,
                           cv_2yr_prob = cv_2yr_probs),
                file = file.path(root,
                                 organ,
                                 model_type,
                                 paste0(age_type, "_cv_prob_", mi_i, ".csv")),
                row.names = FALSE,
                col.names = FALSE,
                sep = ",")
  }
  
  return(model_time = as.numeric(difftime(end_time, start_time, units = "mins")))
}





