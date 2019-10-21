###
### Purpose: Create the variables that allow covariate effects to change over the follow-up period.  Expects data in a long format.
###
### Author: Andrew Wey
###


create_timevarying_vars <- function(curr_data,
                                    curr_vars,
                                    curr_vars_type = NULL,
                                    bl_cuts){
  
  ### curr_data: A long data set for the current iteration of MI
  ### curr_vars: The current variables for the survival from listing model.  These are the variables that will be transformed into time-varying components.
  ### curr_vars_type: The variable types for the current variables. For backwards compatibility, this can be NULL
  ### bl_cuts: The current cutpoints for the baseline hazard
  
  
  ### Creating new vectors with variable names
  var_name_vec <- c("intercept", curr_vars)
  
  
  ### If curr_vars_type is not NULL, then start tracking variable types
  if(! is.null(curr_vars_type)){
    var_type_vec <- c("categorical", curr_vars_type)
  }
  
  
  if(any(duplicated(var_name_vec))){
    stop("Duplicated names in variable name vector")
  }
  
  
  ### Creating separate variables for each piece of the baseline hazard
  new_vars <- NULL
  new_vars_type <- NULL
  
  
  for(k in var_name_vec){
    
    # Adding the primary variable to the covariate list
    if(k != "intercept"){
      
      new_vars <- c(new_vars, 
                    k)
      
      
      # If curr_vars_type is not NULL, then track the variable type
      if(! is.null(curr_vars_type)){
        new_vars_type <- c(new_vars_type,
                           var_type_vec[which(var_name_vec == k)])
      }
    }
    
    
    for(l in 1:length(bl_cuts)){
      
      if(l < length(bl_cuts)){
        
        lh_cut <- bl_cuts[l]
        rh_cut <- bl_cuts[l + 1]
        
        
        curr_time <- paste0(lh_cut,
                            "_",
                            rh_cut)
      }
      if(l == length(bl_cuts)){
        
        lh_cut <- bl_cuts[l]
        
        
        curr_time <- paste0(lh_cut,
                            "_")
      }
      
      
      ### We should know if the time period is not in the data
      if(! curr_time %in% unique(curr_data$time)){ warning(paste(curr_time, "was not an observed time period.")) }
      
      
      curr_data[[paste0(k, ".fup_", curr_time)]] <- ifelse(curr_data$time == curr_time,
                                                           curr_data[[k]],
                                                           0)
      
      
      # Updating the new variable list
      new_vars <- c(new_vars,
                    paste0(k, ".fup_", curr_time))
      
      
      # If curr_vars_type is not NULL, then track the variable type
      if(! is.null(curr_vars_type)){
        new_vars_type <- c(new_vars_type,
                           var_type_vec[which(var_name_vec == k)])
      }
    }
  }
  
  
  ### Constructing the final list, depends on curr_vars_type
  if(! is.null(curr_vars_type)){
    
    fnl_list <- list(new_data = curr_data,
                     new_vars = new_vars,
                     new_vars_type = new_vars_type)
  } else {
    
    fnl_list <- list(new_data = curr_data,
                     new_vars = new_vars)
  }
  
  
  return(fnl_list)
}






