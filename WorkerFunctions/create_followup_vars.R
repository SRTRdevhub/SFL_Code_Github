###
### Purpose: Create the followup variables (event and followup time) for each portion of the baseline hazard. Expects data in a wide format.
###
### Author: Andrew Wey
###



create_followup_vars <- function(curr_data,
                                 bl_cuts){
  
  ### curr_data: Data set for creating the followup variables
  ### bl_cuts: The cutpoints for creating the baseline hazard

  for(l in 1:length(bl_cuts)){
    
    if(l < length(bl_cuts)){
      
      lh_cut <- bl_cuts[l]
      rh_cut <- bl_cuts[l + 1]
      
      
      curr_data[[paste0("days_at_risk_", lh_cut, "_", rh_cut)]] <- ifelse({curr_data$left_time < rh_cut} & {curr_data$right_time >= lh_cut},
                                                                          ifelse(curr_data$left_time >= lh_cut,
                                                                                 ifelse(curr_data$right_time < rh_cut,
                                                                                        curr_data$right_time - curr_data$left_time + 1,
                                                                                        rh_cut - curr_data$left_time),
                                                                                 ifelse(curr_data$right_time < rh_cut,
                                                                                        curr_data$right_time - lh_cut + 1,
                                                                                        rh_cut - lh_cut)),
                                                                          0)
      
      
      curr_data[[paste0("event_", lh_cut, "_", rh_cut)]] <- ifelse({curr_data$right_time >= lh_cut} & {curr_data$right_time < rh_cut},
                                                                   curr_data$event, 
                                                                   0)
    }
    if(l == length(bl_cuts)){
      
      lh_cut <- bl_cuts[l]
      

      curr_data[[paste0("days_at_risk_", lh_cut, "_")]] <- ifelse(curr_data$right_time >= lh_cut,
                                                                  ifelse(curr_data$left_time >= lh_cut,
                                                                         curr_data$right_time - curr_data$left_time + 1,
                                                                         curr_data$right_time - lh_cut + 1),
                                                                  0)
      
      
      curr_data[[paste0("event_", lh_cut, "_")]] <- ifelse(curr_data$right_time >= lh_cut,
                                                           curr_data$event,
                                                           0)
    }
  }
  
  
  return(list(new_data = curr_data))
}







