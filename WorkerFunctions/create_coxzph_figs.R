
### A function that creates plots to investigate the non-proportional hazards of survival from listing models



create_coxzph_figs <- function(root,
                               organ,
                               age_type,
                               model_obj,
                               model_data,
                               mi_i,
                               follow_up_cuts){
  
  ### root: The base location for the survival from listing models
  ### organ: The current organ 
  ### model_obj: The fitted coxph_constant model object for an MI iteration
  ### model_data: The data for the fitted coxph_constant model
  ### mi_i: The current MI iteration; used in naming the output files
  ### follow_up_cuts: The cutpoints for the piecewise exponential models
  
  
  ### The zph object
  coxzph_object <- cox.zph(model_obj,
                           transform = 'identity')
  
  
  ### The plotting of zph uses a scaled version of the timezp variable. Need to determine the appropriate factor
  scaling_fac <- max(coxzph_object$x) / max(model_data$right_time)
  
  
  ### Creating the vector for the names of each variable
  var_names <- dimnames(coxzph_object$y)[[2]]
  
  
  ### If it does not exist, creating the directory
  if(! dir.exists(file.path(root, organ, "zph_check"))){
    dir.create(file.path(root, organ, "zph_check"))
  }
  
  
  ### Creating the PDF and looping through each name
  pdf(file.path(root, organ, "zph_check", paste0("coxph_constant_", age_type, "_mi", mi_i, ".pdf")))
  
  for(i in 1:length(var_names)){
    
    if(! is.na(coxzph_object$table[i,1])){
      
      plot(coxzph_object[i],
           resid = FALSE)
      title(main = var_names[i])
      for(ct in follow_up_cuts){ abline(v = ct, lwd = 2, lty = 3) }
    }
  }
  
  dev.off()
  
  
  return(NULL)
}








