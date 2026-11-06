library(readxl)
library(magrittr)
library(tidyverse)
library(janitor)
library(lubridate)
library(forcats)
# Check whether total = sum (cases by age)
load("./data/work_d_toronto_2024_08_08.RData")
case_age <- read_excel("./data/Weekly case counts or rates of COVID-19 by age group in Toronto Public Health.xlsx",
                       col_types = c("date", "numeric", "date", "date", "text"),
                       skip = 1)

test_age <- read_excel("./data/SARS-CoV-2 weekly total tests, positive tests, or percent positivity by age group in Toronto Public Health.xlsx", 
                       skip = 1)

hosp_age <- read_excel("./data/COVID-19 hospital admissions by age group in Toronto Public Health-3.xlsx", 
                       skip = 1)

case_age <- case_age %>% 
  clean_names() %>% 
  mutate(age_group_report = factor(age_group_report)) %>% 
  rename(case = number_of_cases,
         age=age_group_report) %>% 
  mutate_at(vars(grep("_date", colnames(.), value = TRUE)), ymd)

test_age <- test_age %>% 
  clean_names() %>% 
  mutate(selected_age_groups = factor(selected_age_groups)) %>% 
  rename(tests = total_number_of_tests,
         age=selected_age_groups) %>% 
  mutate_at(vars(grep("_date", colnames(.), value = TRUE)), ymd)

hosp_age <- hosp_age %>% 
  clean_names() %>% 
  mutate(age_group_ba_format = factor(age_group_ba_format)) %>% 
  rename(hosp = number,
         age=age_group_ba_format) %>% 
  mutate_at(vars(grep("_date", colnames(.), value = TRUE)), ymd)

which(!((case_age %>% 
  group_by(earliest_week_end_date) %>% 
  summarise(case = sum(case, na.rm = TRUE)) %>% 
filter(earliest_week_end_date <= "2024-06-01" & earliest_week_end_date>="2020-03-07")%$% case) == work_d_toronto$number_of_cases))

## good. 

### Look into the unknown age groups
case_age %>% filter(age=="Unknown") %>% 
  filter(!is.na(case)) %$% case %>% sum() #49 cases with unknown age groups.


## Impute the missing ages

case_age <- case_age %>% 
  filter(earliest_week_end_date>= "2020-03-22") %>% 
  mutate(case = ifelse(is.na(case), 0, case)) %>% 
  group_by(week_start_date) %>%
  mutate(unknown_cases = case[age == "Unknown"]) %>%
  filter(age != "Unknown") %>% 
  mutate(case_prop = case/sum(case)) %>% 
  mutate(case_toadd = ceiling(unknown_cases*case_prop)) %>% 
  arrange(earliest_week_end_date,desc(case_prop)) %>% 
  mutate(case_toadd_rollcumsum = cumsum(case_toadd),
         diff = unknown_cases-case_toadd_rollcumsum) %>% 
  group_by(earliest_week_end_date) %>% 
  # filter(!any(diff==0))
  mutate(case_imp = ifelse(diff>=0, case_toadd+case, case)) %>% 
  dplyr::select(1:6,case_imp) 

case_age %>% 
  group_by(earliest_week_end_date) %>% 
  summarise(case = sum(case),
            case_imp = sum(case_imp),
            unknown_cases = unknown_cases[1]) %>% 
  filter((case + unknown_cases)!=case_imp)

hosp_age$hosp %>% summary()
test_age$tests %>% summary()


### Filter the dates to match wastewaeter

case_age <- case_age %>% 
  filter(earliest_week_end_date <="2024-04-30") 

c("<1", "1 to 4","5 to 11","12 to 19", "20 to 39", "40 to 59", "60 to 79", "80+")  

hosp_age <- hosp_age %>% 
  rename(earliest_week_start_date =earliest_overalldatetable_week_start_date,
         earliest_week_end_date = week_end_date) %>% 
  filter(earliest_week_end_date <="2024-04-30") 

c("0 to 4",   "5 to 11",  "12 to 17", "18 to 39", "40 to 59", "60 to 79", "80+")

test_age$age %>% unique()
c("<1", "1 to 4",   "5 to 11",  "12 to 19", "20 to 64", "65+" )
####### Align ages %>% unique()

## Fix hosp age first

case_age <- case_age %>% 
  mutate(age = fct_collapse(age, `0 to 4` = c("<1","1 to 4")),
         age = fct_collapse(age, "12 to 39"= c("12 to 19","20 to 39"))) %>% 
  group_by(age, earliest_week_end_date,earliest_week_start_date, week_start_date) %>% 
  summarise(case = sum(case))

hosp_age <- hosp_age %>% 
  mutate(age = fct_collapse(age, "12 to 39"= c("12 to 17", "18 to 39"))) %>% 
  group_by(age,earliest_week_start_date, earliest_week_end_date) %>% 
  summarise(hosp = sum(hosp))

hosp_age$age <- factor(hosp_age$age, levels = c("0 to 4","5 to 11",  "12 to 39", "40 to 59", "60 to 79", "80+" ))
case_age$age <- factor(case_age$age, levels =  c("0 to 4","5 to 11",  "12 to 39", "40 to 59", "60 to 79", "80+" ))

test_age <- test_age %>% 
  mutate(age = factor(age, levels = c("<1", "1 to 4",   "5 to 11",  "12 to 19", "20 to 64", "65+" )))

#########
case_age$earliest_week_end_date %>% range()
hosp_age$earliest_week_end_date %>% range()

case_age <- case_age %>% 
  left_join(hosp_age %>% ungroup() %>%  dplyr::select("age","earliest_week_end_date","hosp"), by = c("age","earliest_week_end_date"))

save(file = './data/work_case_adm_test_age.RData',  list= c("case_age","test_age"))
