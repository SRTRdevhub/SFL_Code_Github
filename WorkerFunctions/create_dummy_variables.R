
# Use create_var_for_levels and control table information to create the additional dummy variables for model fitting

create_dummy_variables <- function(dataset,
                                   raw_dataset, 
                                   covariates,
                                   cov_method,
                                   control_table){
  
  ### dataset: The multiply imputed dataset with the categorical variables
  ### raw_dataset: The dataset before multiple imputation (includes missing values)
  ### covariates: The variables included in the model
  ### cov_method: The type of variable
  ### control_table: The control table for the current model
  
  
  # A vector of variable names that will be used in model fitting
  new_variables <- NULL
  
  
  for(k in 1:length(covariates)){
    
    covar <- covariates[k]
    
    if(cov_method[k] %in% c("binary", "categorical")){
      
      cat_levels <- unlist(control_table[control_table$var_name == covar, 
                                         grepl("level_", names(control_table)) & 
                                           {! grepl("nice_level_", names(control_table))}])
      
      
      cat_levels <- cat_levels[! cat_levels %in% c("", NA)]
      
      
      if(sum(! is.na(dataset[[covar]])) == 0){ stop("all missing") }
      
      
      tmp_update <- create_var_for_levels(dataset,
                                          raw_dataset,
                                          covar,
                                          cat_levels)
      
      
      # Updated dataset
      dataset <- tmp_update$new_data
      
      # Variables to include in the model
      new_variables <- c(new_variables,
                         tmp_update$new_vars)
      
      rm(list = c("tmp_update"))
      gc()
    }
  }
  
  return(list(new_data = dataset,
              new_vars = new_variables))
}



















