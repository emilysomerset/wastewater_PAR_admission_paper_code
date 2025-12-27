library(rstan) # rstan_2.32.7 
library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
# library(aghq) # aghq_0.4.1
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(forcats)
library(bayesplot)
library(ggplot2)
library(lubridate)

load("./age epidemic curve Toronto/epidemic_model_data_noadm.RData")
load("./age epidemic curve Toronto/results_nolag_toronto_thinned_noadm.RData")
load("./data/seroprevalence_by_age_canada.RData")

sero_keep_can <- x %>% 
  filter(interpretation=="infection-induced") %>% 
  rename("earliest_week_end_date"="samplingdate") 

sero_keep_can <- sero_keep_can %>% 
  rename(age = age_grp) %>% 
  mutate(age = factor(age, unique(sero_keep_can$age_grp),
                      labels = c(
                        "1 to 19",
                        "17 to 24",
                        "18 to 34",
                        "22 to 39",
                        "25 to 39",
                        "20 to 59",
                        "35 to 49",
                        "40 to 49",
                        "40 to 59",
                        "50 to 59",
                        "50 to 64",
                        "60+",
                        "60 to 69",
                        "50+",
                        "65+",
                        "70 to 93"
                      )))

load("./data/seroprevalence_by_age_provincial.RData")

sero_keep_ont <- x %>% 
  filter(interpretation=="infection-induced") %>% 
  rename("earliest_week_end_date"="samplingdate") %>% 
  filter(geo=="ON") %>% 
  rename("age"= age_grp) %>% 
  mutate(age = factor(age, levels = c(
    "'6 months to 2 years'",
    "'2 to 5 years'",
    "'5 to 10 years'",
    "'10 to 17 years'",
    "'17-24'",
    "'18-34'",
    "'15-50'",
    "'25-39'",
    "'35-49'",
    "'40-59'",
    "'50-64'",
    "'60+'",
    "'65+'"
  ),
  labels = c(
    "6 months to 2 years",
    "2 to 5 years",
    "5 to 10 years",
    "10 to 17 years",
    "17 to 24",
    "18 to 34",
    "15 to 50",
    "25 to 39",
    "35 to 49",
    "40 to 59",
    "50 to 64",
    "60+",
    "65+"
  )))



final_summary$age <- factor(final_summary$age, levels = 1:4, labels = levels(case_age_merged$age))
final_summary$earliest_week_end_date <- (final_summary$time-1)*7 + min(case_merged$earliest_week_end_date)

age_count_df <- data.frame(age = c("0 to 39","40 to 59", "60 to 79", "80+","60+",
                                   "60 to 69", "70+","70 to 79"),
                           cc = c(384295+134810+185650+243955+246785+213810,
                                  185820+175875+ 184060+191380,
                                  170935+141550+118910+81880,
                                  62790+71855,
                                  170935+141550+118910+81880 + 62790+71855,
                                  170935 + 141550, 
                                  118910+81880+62790+71855,
                                  118910+81880))

count_60plus = 6.65/100*age_count_df$cc[age_count_df$age=="60 to 69"] + 5.74/100*age_count_df$cc[age_count_df$age=="70+"]
count_60plus_upr = 8.99/100*age_count_df$cc[age_count_df$age=="60 to 69"] + 8.71/100*age_count_df$cc[age_count_df$age=="70+"]
count_60plus_lwr = 4.375/100*age_count_df$cc[age_count_df$age=="60 to 69"] + 3.67/100*age_count_df$cc[age_count_df$age=="70+"]
inc_60plus = count_60plus/age_count_df$cc[age_count_df$age == "60+"]
inc_60plus_upr = count_60plus_upr/age_count_df$cc[age_count_df$age == "60+"]
inc_60plus_lwr = count_60plus_lwr/age_count_df$cc[age_count_df$age == "60+"]

# df2 = data.frame(earliest_week_end_date=rep(lubridate::ymd("2021-03-01"),6),
#                  seroprev_est = c(4.72, 5.54, inc_60plus*100,6.65,5.74,5.74),
#                  seroprev_lo95 = c(2.40,3.73, inc_60plus_lwr*100,4.375,3.67,3.67),
#                  seroprev_hi95 = c(7.96, 7.58,inc_60plus_upr*100,8.99,8.71,8.71),
#                  age = c("0 to 39","40 to 59","60+","60 to 69", "70+","70+")) %>%
#   mutate(project = "Ab-C")


df2 = data.frame(earliest_week_end_date=rep(lubridate::ymd("2021-03-01"),4),
                 seroprev_est = c(4.72, 5.54,6.65,5.74)/100,
                 seroprev_lo95 = c(2.40,3.73, 4.375,3.67)/100,
                 seroprev_hi95 = c(7.96, 7.58,8.99,8.71)/100,
                 age = c("0 to 39","40 to 59","60 to 69", "70+")) %>%
  mutate(project = "Ab-C")

sero_keep_can <- sero_keep_can %>% 
  filter(earliest_week_end_date>= min(final_summary$earliest_week_end_date)) %>% 
  dplyr::select(project, earliest_week_end_date, seroprev_est, seroprev_lo95, seroprev_hi95,age) %>% 
  full_join(df2)

sero_keep_ont <- sero_keep_ont %>% 
  filter(earliest_week_end_date>= min(final_summary$earliest_week_end_date)) %>% 
  dplyr::select(project, earliest_week_end_date, seroprev_est, seroprev_lo95, seroprev_hi95,age) %>% 
  full_join(df2)

sero_keep_can <- sero_keep_can %>% 
  mutate(age = factor(age, levels = c(
    "1 to 19",
    "17 to 24",
    "18 to 34",
    "20 to 59",
    "22 to 39",
    "25 to 39",
    "35 to 49",
    "40 to 49",
    "40 to 59",
    "50 to 59",
    "50 to 64",
    "50+",
    "60 to 69",
    "60+",
    "65+",
    "70 to 93",
    "70+"
  )))

sero_keep_ont <- sero_keep_ont %>% 
  mutate(age = factor(age, levels = c(
    "6 months to 2 years",
    "0 to 39",
    "2 to 5 years",
    "5 to 10 years",
    "10 to 17 years",
    "15 to 50",
    "17 to 24",
    "18 to 34",
    "25 to 39",
    "35 to 49",
    "40 to 59",
    "50 to 64",
    "60 to 69",
    "60+",
    "65+",
    "70+"
  )))
library(grid)

mytheme <- theme_bw() +
  theme(
    legend.position = "top",        
    legend.justification = "left",  
    legend.box = "vertical",        
    legend.box.just = "left",       
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.margin = margin(t = -9, r = 0, b = -3, l = 0),
    legend.box.margin = margin(t = 0, r = 0, b = -7, l = 0),
    axis.title.y = element_text(size = 11),
    axis.title.x = element_blank(),
    axis.ticks.x.top = element_blank(),
    legend.key.width = unit(0.4, "cm")
  ) 

myguides <- guides(
  shape = guide_legend(title = "Serosurvey study", nrow = 1, direction = "horizontal"),
  col   = guide_legend(title = "Serosurvey age", nrow = 1, direction = "horizontal",
                       override.aes = list(linetype = 1, shape = NA, size = 1.2))
)


age_ser = c("60+", "60 to 69", "70+","65+")
gg_60_79=final_summary %>% 
  filter(age == "60 to 79") %>% 
  filter(earliest_week_end_date<"2021-12-01") %>%
  left_join(age_count_df, by = c("age")) %>%
  mutate(cum_med = X_cumsum_med/cc*100,
         cum_lwr = X_cumsum_q025/cc*100,
         cum_upr = X_cumsum_q975/cc*100) %>%
  ggplot(aes(earliest_week_end_date, cum_med))+
  geom_ribbon(aes(ymin =  cum_lwr, ymax = cum_upr), alpha = 0.2)+
  geom_line()+ 
  scale_y_continuous(name = expression(C[60*"-"*79 * "," ~ t]~"(%)"), breaks = scales::pretty_breaks(n=8))+
  geom_line(data=case_merged %>% group_by(age) %>%
              mutate(case = cumsum(case)) %>%
              left_join(age_count_df, by = "age") %>%
              mutate(case = case/cc*100) %>%
              filter(earliest_week_end_date<"2021-12-01")%>% 
              filter(age == "60 to 79")  ,
            aes(earliest_week_end_date,case), linetype = "dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=8), name = "Date")+
  geom_point(data = sero_keep_ont %>% 
               filter(earliest_week_end_date<"2021-12-01") %>% 
               filter(age %in% age_ser)   , 
             aes(earliest_week_end_date, seroprev_est*100, shape = factor(project),col = age),inherit.aes = FALSE)+
  geom_errorbar(data = sero_keep_ont %>%  
                  filter(earliest_week_end_date<"2021-12-01") %>% 
                  filter(age %in% age_ser)  , aes(earliest_week_end_date, ymin = seroprev_lo95*100, ymax = seroprev_hi95*100,col = age),inherit.aes = FALSE)+
  geom_segment(data = sero_keep_ont %>% filter(project == "Ab-C")%>%
                 filter(age %in% c("60 to 69","70+")) , aes(x = ymd("2021-02-01"), xend = ymd("2021-03-31"), y = seroprev_lo95*100, col = age)) +
  geom_segment(data = sero_keep_ont %>% filter(project == "Ab-C")%>%
                 filter(age %in% c("60 to 69","70+")) , aes(x = ymd("2021-02-01"), xend = ymd("2021-03-31"), y = seroprev_hi95*100, col = age)) +
  mytheme + myguides+
  scale_shape_manual(values = c(`Ab-C` = 19, CBS = 17, ANTE = 15))  # filled shapes

age_ser = c("35 to 49",
            "40 to 59",
            "50 to 64")

gg_40_59=final_summary %>% 
  filter(age == "40 to 59") %>% 
  filter(earliest_week_end_date<"2021-12-01") %>%
  left_join(age_count_df, by = c("age")) %>%
  mutate(cum_med = X_cumsum_med/cc*100,
         cum_lwr = X_cumsum_q025/cc*100,
         cum_upr = X_cumsum_q975/cc*100) %>%
  ggplot(aes(earliest_week_end_date, cum_med))+
  geom_ribbon(aes(ymin =  cum_lwr, ymax = cum_upr), alpha = 0.2)+
  geom_line()+ 
  scale_y_continuous(name = expression(C[40*"-"*59 * "," ~ t]~"(%)"), breaks = scales::pretty_breaks(n=8))+
  geom_line(data=case_merged %>% group_by(age) %>%
              mutate(case = cumsum(case)) %>%
              left_join(age_count_df, by = "age") %>%
              mutate(case = case/cc*100) %>%
              filter(earliest_week_end_date<"2021-12-01")%>% 
              filter(age == "40 to 59")  ,
            aes(earliest_week_end_date,case), linetype = "dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=8), name = "Date")+
  geom_point(data = sero_keep_ont %>% 
               filter(earliest_week_end_date<"2021-12-01") %>% 
               filter(age %in% age_ser)   , 
             aes(earliest_week_end_date, seroprev_est*100, shape = factor(project),col = age),inherit.aes = FALSE)+
  geom_errorbar(data = sero_keep_ont %>%  
                  filter(earliest_week_end_date<"2021-12-01") %>% 
                  filter(age %in% age_ser)  , aes(earliest_week_end_date, ymin = seroprev_lo95*100, ymax = seroprev_hi95*100,col = age),inherit.aes = FALSE)+
  geom_segment(data = sero_keep_ont %>% filter(project == "Ab-C")%>%
                 filter(age %in% c("40 to 59")) , aes(x = ymd("2021-02-01"), xend = ymd("2021-03-31"), y = seroprev_lo95*100, col = age)) +
  geom_segment(data = sero_keep_ont %>% filter(project == "Ab-C")%>%
                 filter(age %in% c("40 to 59")) , aes(x = ymd("2021-02-01"), xend = ymd("2021-03-31"), y = seroprev_hi95*100, col = age))  +
  mytheme +myguides+
  scale_shape_manual(values = c(`Ab-C` = 19, CBS = 17, ANTE = 15))  # filled shapes


age_ser = c("15 to 50",
            "17 to 24",
            "18 to 34",
            "25 to 39",
            "0 to 39")

gg_0_39=final_summary %>% 
  filter(age == "0 to 39") %>% 
  filter(earliest_week_end_date<"2021-12-01") %>%
  left_join(age_count_df, by = c("age")) %>%
  mutate(cum_med = X_cumsum_med/cc*100,
         cum_lwr = X_cumsum_q025/cc*100,
         cum_upr = X_cumsum_q975/cc*100) %>%
  ggplot(aes(earliest_week_end_date, cum_med))+
  geom_ribbon(aes(ymin =  cum_lwr, ymax = cum_upr), alpha = 0.2)+
  geom_line()+ 
  scale_y_continuous(name = expression(C[0*"-"*39 * "," ~ t]~"(%)"), breaks = scales::pretty_breaks(n=8))+
  geom_line(data=case_merged %>% group_by(age) %>%
              mutate(case = cumsum(case)) %>%
              left_join(age_count_df, by = "age") %>%
              mutate(case = case/cc*100) %>%
              filter(earliest_week_end_date<"2021-12-01")%>% 
              filter(age == "0 to 39")  ,
            aes(earliest_week_end_date,case), linetype = "dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=8), name = "Date")+
  geom_point(data = sero_keep_ont %>% 
               filter(earliest_week_end_date<"2021-12-01") %>% 
               filter(age %in% age_ser)   , 
             aes(earliest_week_end_date, seroprev_est*100, shape = factor(project),col = age),inherit.aes = FALSE)+
  geom_errorbar(data = sero_keep_ont %>%  
                  filter(earliest_week_end_date<"2021-12-01") %>% 
                  filter(age %in% age_ser)  , aes(earliest_week_end_date, ymin = seroprev_lo95*100, ymax = seroprev_hi95*100,col = age),inherit.aes = FALSE)+
  geom_segment(data = sero_keep_ont %>% filter(project == "Ab-C")%>%
                 filter(age %in% c("0 to 39")) , aes(x = ymd("2021-02-01"), xend = ymd("2021-03-31"), y = seroprev_lo95*100, col = age)) +
  geom_segment(data = sero_keep_ont %>% filter(project == "Ab-C")%>%
                 filter(age %in% c("0 to 39")) , aes(x = ymd("2021-02-01"), xend = ymd("2021-03-31"), y = seroprev_hi95*100, col = age))  +
  mytheme+ myguides+
  scale_shape_manual(values = c(`Ab-C` = 19, CBS = 17, ANTE = 15))  # filled shapes




age_ser = c(    "60+",
                "65+",
                "70+")

gg_80plus=final_summary %>% 
  filter(age == "80+") %>% 
  filter(earliest_week_end_date<"2021-12-01") %>%
  left_join(age_count_df, by = c("age")) %>%
  mutate(cum_med = X_cumsum_med/cc*100,
         cum_lwr = X_cumsum_q025/cc*100,
         cum_upr = X_cumsum_q975/cc*100) %>%
  ggplot(aes(earliest_week_end_date, cum_med))+
  geom_ribbon(aes(ymin =  cum_lwr, ymax = cum_upr), alpha = 0.2)+
  geom_line()+ 
  scale_y_continuous(name = expression(C[80*"+"* "," ~ t]~"(%)"), breaks = scales::pretty_breaks(n=8))+
  geom_line(data=case_merged %>% group_by(age) %>%
              mutate(case = cumsum(case)) %>%
              left_join(age_count_df, by = "age") %>%
              mutate(case = case/cc*100) %>%
              filter(earliest_week_end_date<"2021-12-01")%>% 
              filter(age == "80+")  ,
            aes(earliest_week_end_date,case), linetype = "dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=8), name = "Date")+
  geom_point(data = sero_keep_ont %>% 
               filter(earliest_week_end_date<"2021-12-01") %>% 
               filter(age %in% age_ser)   , 
             aes(earliest_week_end_date, seroprev_est*100, shape = factor(project),col = age),inherit.aes = FALSE)+
  geom_errorbar(data = sero_keep_ont %>%  
                  filter(earliest_week_end_date<"2021-12-01") %>% 
                  filter(age %in% age_ser)  , aes(earliest_week_end_date, ymin = seroprev_lo95*100, ymax = seroprev_hi95*100,col = age),inherit.aes = FALSE)+
  geom_segment(data = sero_keep_ont %>% filter(project == "Ab-C")%>%
                 filter(age %in% c("70+")) , aes(x = ymd("2021-02-01"), xend = ymd("2021-03-31"), y = seroprev_lo95*100, col = age)) +
  geom_segment(data = sero_keep_ont %>% filter(project == "Ab-C")%>%
                 filter(age %in% c("70+")) , aes(x = ymd("2021-02-01"), xend = ymd("2021-03-31"), y = seroprev_hi95*100, col = age))  +
  mytheme+ myguides+
  scale_shape_manual(values = c(`Ab-C` = 19, CBS = 17, ANTE = 15))  # filled shapes



ggsave(filename = paste0("./age epidemic curve Toronto/plots/cumu_older.pdf"),
       plot = gg_60_79, 
       device = "pdf",
       width = 4, 
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto/plots/cumu_adult.pdf"),
       plot = gg_40_59, 
       device = "pdf",
       width = 4, 
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto/plots/cumu_young.pdf"),
       plot = gg_0_39, 
       device = "pdf",
       width = 4, 
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto/plots/cumu_oldest.pdf"),
       plot = gg_80plus, 
       device = "pdf",
       width = 4, 
       height = 8/3,
       dpi = 300)


#####################
final_summary %>%
  mutate_if(is.numeric, function(x){x*100}) %>% 
  mutate(
    age1_CI = sprintf("%.2f [%.2f, %.2f]", age1_med, age1_q025, age1_q975),
    age2_CI = sprintf("%.2f [%.2f, %.2f]", age2_med, age2_q025, age2_q975),
    age3_CI = sprintf("%.2f [%.2f, %.2f]", age3_med, age3_q025, age3_q975),
    age4_CI = sprintf("%.2f [%.2f, %.2f]", age4_med, age4_q025, age4_q975)
  ) %>%
  select(age, ends_with("_CI")) %>% 
  slice(1:4)
###################

ggyoung <- final_summary %>% 
  left_join(case_merged, by = c("earliest_week_end_date","age")) %>% 
  filter(age == "0 to 39") %>% 
  ggplot(aes(earliest_week_end_date, value_int_med))+ 
  geom_ribbon(aes(ymin = value_int_q025,ymax = value_int_q975),alpha  = 0.3)+
  # geom_ribbon(aes(ymin = value_int_q10,ymax = value_int_q90),alpha  = 0.4)+
  geom_line()+ 
  theme_bw()+ 
  geom_line(aes(earliest_week_end_date, case), linetype = "dashed")+ 
  scale_y_continuous(name = expression(X[0*"-"*39 * "," ~ t]), breaks = scales::pretty_breaks(n=8))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank())


ggadult <- final_summary %>% 
  left_join(case_merged, by = c("earliest_week_end_date","age")) %>% 
  filter(age == "40 to 59") %>% 
  ggplot(aes(earliest_week_end_date, value_int_med))+ 
  geom_ribbon(aes(ymin = value_int_q025,ymax = value_int_q975),alpha  = 0.3)+
  # geom_ribbon(aes(ymin = value_int_q10,ymax = value_int_q90),alpha  = 0.4)+
  geom_line()+ 
  theme_bw()+ 
  geom_line(aes(earliest_week_end_date, case), linetype = "dashed")+ 
  scale_y_continuous(name = expression(X[40*"-"*59 * "," ~ t]), breaks = scales::pretty_breaks(n=8))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank())


ggold <- final_summary %>% 
  left_join(case_merged, by = c("earliest_week_end_date","age")) %>% 
  filter(age == "60 to 79") %>% 
  ggplot(aes(earliest_week_end_date, value_int_med))+ 
  geom_ribbon(aes(ymin = value_int_q025,ymax = value_int_q975),alpha  = 0.3)+
  # geom_ribbon(aes(ymin = value_int_q10,ymax = value_int_q90),alpha  = 0.4)+
  geom_line()+ 
  theme_bw()+ 
  geom_line(aes(earliest_week_end_date, case), linetype = "dashed")+ 
  scale_y_continuous(name = expression(X[60*"-"*79 * "," ~ t]), breaks = scales::pretty_breaks(n=8))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank())


ggolder <- final_summary %>% 
  left_join(case_merged, by = c("earliest_week_end_date","age")) %>% 
  filter(age == "80+") %>% 
  ggplot(aes(earliest_week_end_date, value_int_med))+ 
  geom_ribbon(aes(ymin = value_int_q025,ymax = value_int_q975),alpha  = 0.3)+
  # geom_ribbon(aes(ymin = value_int_q10,ymax = value_int_q90),alpha  = 0.4)+
  geom_line()+ 
  theme_bw()+ 
  geom_line(aes(earliest_week_end_date, case), linetype = "dashed")+ 
  scale_y_continuous(name = expression(X[80*"+"* "," ~ t]), breaks = scales::pretty_breaks(n=8))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank())

ggsave(filename = paste0("./age epidemic curve Toronto/plots/X_older.pdf"),
       plot = ggolder+ theme(axis.text.x.top = element_text(vjust = -78)), 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto/plots/X_young.pdf"),
       plot = ggyoung+ theme(axis.text.x.top = element_text(vjust = -78)), 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto/plots/X_adult.pdf"),
       plot = ggadult+ theme(axis.text.x.top = element_text(vjust = -78)), 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto/plots/X_old.pdf"),
       plot = ggold+theme(axis.text.x.top = element_text(vjust = -78)), 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)

#############
set.seed(2917)
ss <- sample(1:3000,500)
tmp = Reduce(rbind, lapply(as.list(ss),function(i){data_foranalysis_full$analysis_d[[i]]$version = i; data_foranalysis_full$analysis_d[[i]]}))
tmp <- tmp %>% 
  group_by(earliest_week_end_date,y) %>% 
  summarise(ratio_lwr = quantile(ratio,0.025),
            ratio_upr = quantile(ratio,0.975),
            ratio_med = quantile(ratio,0.5),
            ratio_v_u_lwr = quantile(ratio_v_u_fixed,0.025),
            ratio_v_u_upr = quantile(ratio_v_u_fixed,0.975),
            ratio_v_u_lwr2 = quantile(ratio_v_u_fixed,0.1),
            ratio_v_u_upr2 = quantile(ratio_v_u_fixed,0.9),
            ratio_v_u_med = quantile(ratio_v_u_fixed,0.5))


tmp <- tmp %>%ungroup() %>%  
  mutate(lagged_case = lag(y)) %>% 
  mutate(ratio_crude = y/lagged_case)

tmp2 <- tmp %>% 
  left_join(case_merged, by = "earliest_week_end_date")
tmp2 <- tmp2 %>% 
  group_by(age) %>% 
  mutate(lagged_case = lag(case)) %>% 
  mutate(ratio_crude = case/lagged_case)

tmp2 %>% filter(ratio_v_u_lwr >  1& year(earliest_week_end_date)>=2023)


gg_repro = tmp2 %>% 
  ungroup() %>% 
  mutate(age = factor(age, labels = c("0-39","40-59","60-79","80+"))) %>% 
  ggplot()+
  geom_ribbon(aes(earliest_week_end_date,ymax = ratio_v_u_upr, ymin = ratio_v_u_lwr), alpha = 0.3)+ 
  # geom_ribbon(aes(earliest_week_end_date,ymax = ratio_v_u_upr2, ymin = ratio_v_u_lwr2), alpha = 0.4)+  
  geom_line(aes(earliest_week_end_date,ratio_v_u_med))+  
  geom_line(aes(x=earliest_week_end_date, y=ratio_crude, col = age))+
  # geom_point(col="red")+
  theme_bw()+ 
  scale_y_continuous(name = expression(R[t]),
                     breaks = scales::pretty_breaks(n=10))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank(),
        legend.position = c(1, 1), # x and y coordinates, ranging from 0 to 1
        legend.justification = c(1.1, 1.1),
        legend.direction = "horizontal")+ 
  labs(col ="")

ggsave(filename = paste0("./age epidemic curve Toronto/plots/repro.pdf"),
       plot = gg_repro+theme(axis.text.x.top = element_text(vjust = -78)), 
       device = "pdf",
       width = 8/2*2, 
       height = 8/3,
       dpi = 300)


#############

load("./data/work_d_pho_testing_2024_10_07.RData") # Public health ontario testing

ggprob = final_summary %>% 
  group_by(earliest_week_end_date) %>% 
  slice(1) %>% 
  left_join(work_d_testing%>% 
              dplyr::select(week_end_date, total_number_of_tests) %>% 
              mutate(week_end_date = ymd(week_end_date)),
            by = c("earliest_week_end_date"="week_end_date")) %>% 
  ggplot(aes(earliest_week_end_date, Pi_med))+ 
  geom_ribbon(aes(ymin = Pi_q025,ymax = Pi_q975),alpha  = 0.3)+
  geom_line()+ 
  theme_bw()+ 
  scale_y_continuous(name = expression(pi[t]), breaks = scales::pretty_breaks(n=6))+
  scale_x_date(breaks=scales::pretty_breaks(n=6), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank())+
  geom_line(aes(earliest_week_end_date,total_number_of_tests/max(total_number_of_tests)),col = "black", linetype="dashed")

ggsave(filename = paste0("./age epidemic curve Toronto/plots/prob.pdf"),
       plot = ggprob+theme(axis.text.x.top = element_text(vjust = -78)), 
       device = "pdf",
       width = 8/2*2, 
       height = 8/3,
       dpi = 300)



