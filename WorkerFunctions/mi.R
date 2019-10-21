

mi <- function(current_data, n_mi, mi_covariates, mi_cov_method, keep_vars, data_subset, current_dir){
  
  
  if(mean(keep_vars %in% names(current_data)) < 1){ stop("Data set does not contain all 'keep_vars'") }
  
  
  output.path <- current_dir
  
  
  ### Creating outcome variables for use in MI
  current_data <- current_data %>%
    mutate(log_days_at_risk_time = log(right_time - left_time + 1))
  
  
  
  mi_method_dictionary <- rbind(c("continuous", "binary", "categorical"), c("pmm", "logreg", "polyreg"))
  mi_method            <- mi_method_dictionary[2, match(mi_cov_method, mi_method_dictionary[1,])]
  
  
  mi_data    <- as.data.frame(current_data[, c("log_days_at_risk_time", "event", mi_covariates)])
  mi_factors <- which(c("continuous", "binary", mi_cov_method) != "continuous")
  
  
  current_data <- current_data[, c(keep_vars)]
  
  
  for(j in mi_factors){ mi_data[,j] <- as.factor(mi_data[,j]) }
  
  
  set.seed(1)
  mi_output <- mice(mi_data, 
                    m = n_mi, 
                    method = c("pmm", "logreg", mi_method))
  
  
  for(j in 1:n_mi){
    mi_data <- cbind(current_data, complete(mi_output, j), MI = j)
    save(mi_data, file = file.path(output.path, paste0("mi_", organ_ab, "_", data_subset, "_data_", j, ".Rdata")))
    
    rm(list = c("mi_data"))
  }
  
  return(NULL)
}
