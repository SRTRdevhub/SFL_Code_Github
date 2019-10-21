#
# A function to create the outcome variables for each of the survival from listing models
#


create_outcome_vars <- function(curr_data,
                                max_followup,
                                cohort_start,
                                cohort_end){
  
  
  ### Calculating the important followup variables for survival from listing
  curr_data <- curr_data %>% 
    mutate(max_followup_censor_dt = as.Date(CAN_LISTING_DT, format = "%m/%d/%y") + (365 * max_followup),
           numeric_listing_dt = as.numeric(as.Date(CAN_LISTING_DT, format = "%m/%d/%y")),
           numeric_removal_dt = ifelse(CAN_REM_DT == "",
                                       Inf,
                                       as.numeric(as.Date(CAN_REM_DT, format = "%m/%d/%y"))),
           death_dt = ifelse({CAN_DEATH_DT != ""} & {PERS_ALL_DEATH_DT != ""},
                             ifelse(as.Date(CAN_DEATH_DT, format = "%m/%d/%y") <= as.Date(PERS_ALL_DEATH_DT, format = "%m/%d/%y"),
                                    CAN_DEATH_DT,
                                    PERS_ALL_DEATH_DT),
                             ifelse({CAN_DEATH_DT != ""} & {PERS_ALL_DEATH_DT == ""},
                                    CAN_DEATH_DT,
                                    ifelse({CAN_DEATH_DT == ""} & {PERS_ALL_DEATH_DT != ""},
                                           PERS_ALL_DEATH_DT,
                                           ifelse(PERS_RESTRICT_DEATH_DT != "",
                                                  PERS_RESTRICT_DEATH_DT,
                                                  "")))),
           death_dt_restricted = ifelse({CAN_DEATH_DT != ""} & {PERS_ALL_DEATH_DT != ""},
                                        0,
                                        ifelse({CAN_DEATH_DT != ""} & {PERS_ALL_DEATH_DT == ""},
                                               0,
                                               ifelse({CAN_DEATH_DT == ""} & {PERS_ALL_DEATH_DT != ""},
                                                      0,
                                                      ifelse(PERS_RESTRICT_DEATH_DT != "",
                                                             1,
                                                             0))))) %>%
    subset({as.Date(CAN_LISTING_DT, format = "%m/%d/%y") <= as.Date(cohort_end, format = "%Y%m%d")} & 
             {max_followup_censor_dt >= as.Date(cohort_start, format = "%Y%m%d")} &
             {{death_dt == ""} | 
                 {as.Date(CAN_LISTING_DT, format = "%m/%d/%y") <= as.Date(death_dt, format = "%m/%d/%y")}} & 
             {{death_dt == ""} | 
                 {as.Date(death_dt, format = "%m/%d/%y") >= as.Date(cohort_start, format = "%Y%m%d")}}) %>%
    group_by(PERS_ID, 
             ctr_id) %>%
    mutate(min_listing_dt = min(numeric_listing_dt)) %>%
    ungroup() %>%
    subset(numeric_listing_dt == min_listing_dt) %>%
    mutate(event_dt = ifelse(death_dt == "", 
                             ifelse(as.Date(cohort_end, format = "%Y%m%d") <= max_followup_censor_dt,
                                    format(as.Date(cohort_end, format = "%Y%m%d"), format = "%m/%d/%y"),
                                    format(max_followup_censor_dt, format = "%m/%d/%y")),
                             ifelse(as.Date(death_dt, format = "%m/%d/%y") > as.Date(cohort_end, format = "%Y%m%d"),
                                    ifelse(as.Date(cohort_end, format = "%Y%m%d") <= max_followup_censor_dt,
                                           format(as.Date(cohort_end, format = "%Y%m%d"), format = "%m/%d/%y"),
                                           format(max_followup_censor_dt, format = "%m/%d/%y")),
                                    ifelse(as.Date(death_dt, format = "%m/%d/%y") <= max_followup_censor_dt,
                                           death_dt,
                                           format(max_followup_censor_dt, format = "%m/%d/%y")))),
           start_dt = ifelse(as.Date(CAN_LISTING_DT, format = "%m/%d/%y") >= as.Date(cohort_start, format = "%Y%m%d"),
                             CAN_LISTING_DT,
                             format(as.Date(cohort_start, format = "%Y%m%d"), format = "%m/%d/%y")),
           
           left_time = as.numeric(as.Date(start_dt, format = "%m/%d/%y") - as.Date(CAN_LISTING_DT, format = "%m/%d/%y")),
           right_time = as.numeric(as.Date(event_dt, format = "%m/%d/%y") - as.Date(CAN_LISTING_DT, format = "%m/%d/%y")),
           event = ifelse(event_dt == death_dt, 1, 0)) %>%
    subset({as.Date(event_dt, format = "%m/%d/%y") >= as.Date(cohort_start, format = "%Y%m%d")} &  # removes the candidates relisted prior to start
             {as.Date(event_dt, format = "%m/%d/%y") >= as.Date(start_dt, format = "%m/%d/%y")})

  
  return(curr_data)
}




