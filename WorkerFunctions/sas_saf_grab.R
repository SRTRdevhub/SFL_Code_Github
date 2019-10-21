

sas_saf_grab <- function(sasdir, data_dir, current_saf, organ_ab, organ_saf_ab, recip_vars){
  
  ## sasdir is the directory where the SAS match run extraction file exists 
  
  
  ## Read the SAS code
  if(file.exists(file.path(sasdir, "saf_extract.sas"))){ sascode <- readLines(file.path(sasdir, "saf_extract.sas")) } 
  if(! file.exists(file.path(sasdir, "saf_extract.sas"))){ stop("SAS file for extracting SAF data does not exist") } 
  
  
  sascoder <- c("OPTIONS NOCENTER NODATE NONOTES NONUMBER",
                "  LINESIZE=80 PAGESIZE=5000 MISSING=.", 
                "  FORMCHAR = '|----|+|---+=|-/<>*';",
                "",
                sascode,
                "",
                "QUIT;",
                "")
  
  
  # This section simply replaces the generic holders WITHIN the SAS code with the appropriate information for the current cohort
  sascoder <- gsub(pattern = "ORGAN_AB", replacement = organ_ab, x = sascoder)
  sascoder <- gsub(pattern = "ORGAN_SAF_AB", replacement = organ_saf_ab, x = sascoder)
  sascoder <- gsub(pattern = "DATA_DIR", replacement = data_dir, x = sascoder)
  sascoder <- gsub(pattern = "CURRENT_SAF", replacement = current_saf, x = sascoder)
  sascoder <- gsub(pattern = "RECIP_VARS", replacement = paste(recip_vars, collapse = " "), x = sascoder)
  
  
  # Adding suffixes to the csv files to ensure the data sets are kept
  sascoder <- gsub(pattern = ".csv", replacement = paste0("_", current_saf, ".csv"), x = sascoder)
  
  
  writeLines(sascoder, file.path(sasdir, paste0("saf_extract_", organ_ab, "_R.sas")))
  
  ## Run the code
  shell(paste0("cd /d ", gsub("/","\\\\",sasdir), " & ", "sas ", paste0("saf_extract_", organ_ab, "_R.sas")))
  
  ## Read the LST file 
  if(file.exists(paste0(sasdir, "/mr_extract", organ_ab,"_R.lst"))) {
    saslst <- readLines(paste0(sasdir, "/mr_extract", organ_ab, "_R.lst")) } else {
      saslst <- NULL }
  
  ## Read the LOG file 
  if(file.exists(paste0(sasdir, "/mr_extract", organ_ab, "_R.log"))) {
    saslog <- readLines(paste0(sasdir, "/mr_extract", organ_ab, "_R.log")) } else {
      saslog <- NULL }
  
  
  
  ######################################################################################################################################################
  ######################################################################################################################################################
  ######################################################################################################################################################
  #
  # Some organs require specific information that is not contained in the candidate SAFs.  I decided to code those differences for each organ
  #   with different SAS functions rather than making accommendations within the control table...
  #   

  if(organ_ab %in% c("ki", "li", "lu", "hr", "pa")){
    
    sascoder_specific <- readLines(file.path(sasdir, paste0("saf_extract_", organ_ab, "_specific.sas")))
    
    sascoder_specific <- gsub(pattern = "DATA_DIR", replacement = data_dir, x = sascoder_specific)
    sascoder_specific <- gsub(pattern = "CURRENT_SAF", replacement = current_saf, x = sascoder_specific)
    
    writeLines(sascoder_specific, file.path(sasdir, paste0("saf_extract_", organ_ab, "_specific_R.sas")))
    
    shell(paste0("cd /d ", gsub("/","\\\\",sasdir), " & ", "sas ", paste0("saf_extract_", organ_ab, "_specific_R.sas")))
  }
  
  
  
  return(list(code=sascoder,lst=saslst,log=saslog))
}
