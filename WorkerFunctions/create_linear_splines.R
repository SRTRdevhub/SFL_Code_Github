
# Create linear splines for continuous variables


create_linear_splines <- function(dataset, continuous_var){
  
  ### dataset: The dataset for the linear splines
  ### continuous_var: The continuous variable to create linear splines for
  
  
  # A vector to track the names of the linear splines
  new_variables <- NULL
  
  
  pknots <- pretty(range(dataset[[continuous_var]]),
                   n = 11)
  
  
  lh_knots <- pknots[-1]
  rh_knots <- pknots[-length(pknots)]
  
  for(l in 1:length(pknots)){
    
    
    if(pknots[l] >= 0){ knot_suffix <- pknots[l] }
    if(pknots[l] <  0){ knot_suffix <- paste0("n", abs(pknots[l])) }
    

    if(l == 1){
      
      # Right-hand linear spline
      dataset[[paste0(continuous_var, "_rh_", knot_suffix)]] <- {dataset[[continuous_var]] > pknots[l]} * (dataset[[continuous_var]] - pknots[l])
      
      # Adding the spline names to vector of new variables
      new_variables <- c(new_variables,
                         paste0(continuous_var, "_rh_", knot_suffix))
    }
    if(l == length(pknots)){
      
      # Left-hand linear spline
      dataset[[paste0(continuous_var, "_lh_", knot_suffix)]] <- {dataset[[continuous_var]] < pknots[l]} * (pknots[l] - dataset[[continuous_var]])
      
      # Adding the spline names to vector of new variables
      new_variables <- c(new_variables,
                         paste0(continuous_var, "_lh_", knot_suffix))
    }
    if(! l %in% c(1, length(pknots))){
      
      # Left-hand linear spline
      dataset[[paste0(continuous_var, "_lh_", knot_suffix)]] <- {dataset[[continuous_var]] < pknots[l]} * (pknots[l] - dataset[[continuous_var]])
      
      # Right-hand linear spline
      dataset[[paste0(continuous_var, "_rh_", knot_suffix)]] <- {dataset[[continuous_var]] > pknots[l]} * (dataset[[continuous_var]] - pknots[l])
      
      
      # Adding the spline names to vector of new variables
      new_variables <- c(new_variables,
                         paste0(continuous_var, "_lh_", knot_suffix),
                         paste0(continuous_var, "_rh_", knot_suffix))
    }
  }
  
  
  return(list(new_data = dataset,
              new_vars = new_variables))
}









