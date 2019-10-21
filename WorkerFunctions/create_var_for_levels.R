
### Function to create variables for each level of a binary or categorical variable


create_var_for_levels <- function(curr_data, 
                                  raw_data,
                                  covar, 
                                  exp_levels){
  
  
  ### curr_data: The multiply imputed dataset
  ### raw_data: The dataset before multiple imputation; used to consistently determine the reference level
  ### covar: The covariate name in curr_data and raw_data to create different levels for
  ### exp_levels: The expected levels for the covariate based on information in the control table
  
  
  curr_table  <- table(curr_data[[covar]])
  
  
  if(sum(! dimnames(curr_table)[[1]] %in% exp_levels) > 0){ stop(paste("Unexpected level in", covar)) }
  

  level_ns <- sapply(exp_levels, 
                     function(xx){ sum(raw_data[[covar]] == xx, na.rm = TRUE) })
  
  
  ref_level <- exp_levels[level_ns == max(level_ns)][1]

  
  new_levels <- exp_levels[! exp_levels %in% ref_level]
  new_names  <- paste(covar, new_levels, sep = ".")
  
  
  for(l in exp_levels){
    curr_data[[paste(covar, l, sep = ".")]] <- ifelse(curr_data[[covar]] == l, 1, 0)
  }
  
  
  return(list(new_data = curr_data,
              new_vars = new_names))
}

