library(rstan) # rstan_2.32.7 
library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
# library(aghq) # aghq_0.4.1
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(forcats)
library(bayesplot)
library(ggplot2)


load('./data/work_case_adm_test_age.RData')


### Get wastewater reproduction 
#############################################################
load(file="./wastewater model/results_model_CAN_phachist_ONfilt.RData")
source("./functions_general/prep_data_covid_with_fittedwastewater.R")

weights_datadic = data.frame(site_id = c("TAB","THC","THU","TNT"), 
                             weight = c(0.499/(0.177+0.232+0.064+0.499),
                                        0.177/(0.177+0.232+0.064+0.499),
                                        0.232/(0.177+0.232+0.064+0.499),
                                        0.064/(0.177+0.232+0.064+0.499)))

start_dates = results$v_u_fixed %>% 
  group_by(site_id) %>% 
  summarise(date_inc = min(sample_date))

data_foranalysis_full <- prep_data(case_data = case_age %>% group_by(earliest_week_end_date) %>% summarise(totcase = sum(case)),
                                   y_var = "totcase",
                                   results = results,
                                   AR=TRUE,
                                   weight_ratio = TRUE,
                                   weights_datadic = weights_datadic,
                                   start_dates = start_dates, 
                                   add_lag = 0)


case_age<- case_age %>% 
  arrange(-desc(age))

case_age_merged <- 
  case_age %>% 
  mutate(age = fct_collapse(age, "0 to 39"= c("0 to 4","5 to 11", "12 to 39"))) %>% 
  group_by(age,earliest_week_end_date, earliest_week_start_date, week_start_date) %>% 
  summarise(case = sum(case),
            hosp = sum(hosp))

###############

data_foranalysis_full$analysis_d2 %>%
  filter(earliest_week_end_date>= "2020-10-10") # 2020-10-17

case_age %>% filter(!is.na(hosp))

# date_range = data_foranalysis_full$analysis_d2 %>% filter(!is.na(ratio)) %$% earliest_week_end_date %>% range()

date_range = data_foranalysis_full$analysis_d2 %>% filter(!is.na(ratio))  %$% earliest_week_end_date %>% range()


init = case_age_merged$case[case_age_merged$earliest_week_end_date == (date_range[1]-7)] 


case_merged <-case_age_merged %>% 
  filter(earliest_week_end_date >= date_range[1] & earliest_week_end_date <= (date_range[2]))

case_age_matrix <- case_merged  %>% dcast(earliest_week_end_date~age, value.var = 'case') 


save(file = "./age epidemic curve Toronto/epidemic_model_data_noadm.RData",list = c("case_age_matrix",
                                                                                    "date_range",
                                                                                    "data_foranalysis_full",
                                                                                    "case_merged",
                                                                                    "case_age_merged"))
                                                                                        