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
                                   start_dates = start_dates)


data_foranalysis_full$analysis_d[[1]] %>% 
  arrange(-desc(earliest_week_end_date)) %>% 
  filter(!is.na(ratio_cases)) %$% range(earliest_week_end_date)

ratio = data_foranalysis_full$analysis_d[[1]] %>% 
  arrange(-desc(earliest_week_end_date)) %>% 
  filter(!is.na(ratio_cases)) %$% ratio_v_u_fixed

####


case_age<- case_age %>% 
  arrange(-desc(age))

case_age %>% 
  filter(!is.na(hosp)) %>% 
  mutate(rate = hosp/case) %>% 
  ggplot(aes(earliest_week_end_date, log(rate)))+ 
  geom_point()+ 
  facet_wrap(~age)


case_age_merged <- 
  case_age %>% 
  mutate(age = fct_collapse(age, "0 to 39"= c("0 to 4","5 to 11", "12 to 39"))) %>% 
  group_by(age,earliest_week_end_date, earliest_week_start_date, week_start_date) %>% 
  summarise(case = sum(case),
            hosp = sum(hosp))



case_age_merged %>% 
  filter(!is.na(hosp)) %>% 
  mutate(prop = hosp/case) %>% 
  ggplot(aes(earliest_week_end_date, log(prop)))+ 
  geom_point()+ 
  facet_wrap(~age)

# What we see is increasing rate of admissions among reported
# Not surprising

###############

date_range = data_foranalysis_full$analysis_d2 %>% filter(!is.na(ratio)) %>% 
  filter(earliest_week_end_date >= "2021-10-16")%$% earliest_week_end_date %>% range()


init = case_age_merged$case[case_age_merged$earliest_week_end_date == (date_range[1]-7)] 


case_merged <-case_age_merged %>% 
  filter(earliest_week_end_date >= date_range[1] & earliest_week_end_date <= (date_range[2]))

case_age_matrix <- case_merged  %>% dcast(earliest_week_end_date~age, value.var = 'case') 

hosp_age_matrix <- case_merged  %>% dcast(earliest_week_end_date~age, value.var = 'hosp') 


(case_age_matrix - hosp_age_matrix) %>% summary()

save(file = "./age epidemic curve Toronto/forcomputer/epidemic_model_nonlag_data.RData", list = c("case_age_matrix",
                                                                                         "hosp_age_matrix",
                                                                                         "date_range",
                                                                                         "data_foranalysis_full",
                                                                                         "case_merged",
                                                                                         "case_age_merged"))
