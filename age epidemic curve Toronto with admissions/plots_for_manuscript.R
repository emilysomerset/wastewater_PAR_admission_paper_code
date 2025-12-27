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

#### 
load("./age epidemic curve Toronto with admissions/results_nolag_toronto_thinned.RData")
load("./age epidemic curve Toronto with admissions/epidemic_model_nonlag_data.RData")
final_summary$age <- factor(final_summary$age, levels = 1:4, labels = levels(case_age_merged$age))
final_summary$earliest_week_end_date <- (final_summary$time-1)*7 + min(case_merged$earliest_week_end_date)

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

ggsave(filename = paste0("./age epidemic curve Toronto with admissions/plots/X_older_joint.pdf"),
       plot = ggolder+ theme(axis.text.x.top = element_text(vjust = -78)),
       device = "pdf",
       width = 8/2,
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto with admissions/plots/X_young_joint.pdf"),
       plot = ggyoung+ theme(axis.text.x.top = element_text(vjust = -78)),
       device = "pdf",
       width = 8/2,
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto with admissions/plots/X_adult_joint.pdf"),
       plot = ggadult+ theme(axis.text.x.top = element_text(vjust = -78)),
       device = "pdf",
       width = 8/2,
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto with admissions/plots/X_old_joint.pdf"),
       plot = ggold+theme(axis.text.x.top = element_text(vjust = -78)),
       device = "pdf",
       width = 8/2,
       height = 8/3,
       dpi = 300)



load("./data/work_d_pho_testing_2024_10_07.RData") # Public health ontario testing

ggprob = final_summary %>% 
  group_by(earliest_week_end_date) %>% 
  slice(1) %>% 
  left_join(work_d_testing%>% 
              dplyr::select(week_end_date, total_number_of_tests) %>% 
              mutate(week_end_date = ymd(week_end_date)),
            by = c("earliest_week_end_date"="week_end_date")) %>% 
  ggplot(aes(earliest_week_end_date, pY_med))+ 
  geom_ribbon(aes(ymin = pY_q025,ymax = pY_q975),alpha  = 0.3)+
  geom_line()+ 
  theme_bw()+ 
  scale_y_continuous(name =expression(pi*`'`[t]), breaks = scales::pretty_breaks(n=6))+
  scale_x_date(breaks=scales::pretty_breaks(n=6), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank())
# geom_line(aes(earliest_week_end_date,total_number_of_tests/max(total_number_of_tests)),col = "black", linetype="dashed")

ggsave(filename = paste0("./age epidemic curve Toronto with admissions/plots/prob_prime.pdf"),
       plot = ggprob+theme(axis.text.x.top = element_text(vjust = -78)),
       device = "pdf",
       width = 8/2*2,
       height = 8/3,
       dpi = 300)


ggtau_young = final_summary %>% 
  filter(age == "0 to 39") %>% 
  group_by(earliest_week_end_date) %>% 
  ggplot(aes(earliest_week_end_date, tau_med))+ 
  geom_ribbon(aes(ymin = tau_q025,ymax = tau_q975),alpha  = 0.3)+
  geom_line()+ 
  theme_bw()+ 
  scale_y_continuous(name = expression(tau[0*"-"*39 * "," ~ t]), breaks = scales::pretty_breaks(n=8))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank())

ggtau_adult = final_summary %>% 
  filter(age == "40 to 59") %>% 
  group_by(earliest_week_end_date) %>% 
  ggplot(aes(earliest_week_end_date, tau_med))+ 
  geom_ribbon(aes(ymin = tau_q025,ymax = tau_q975),alpha  = 0.3)+
  geom_line()+ 
  theme_bw()+ 
  scale_y_continuous(name = expression(tau[40*"-"*59 * "," ~ t]), breaks = scales::pretty_breaks(n=8))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank())

ggtau_old = final_summary %>% 
  filter(age == "60 to 79") %>% 
  group_by(earliest_week_end_date) %>% 
  ggplot(aes(earliest_week_end_date, tau_med))+ 
  geom_ribbon(aes(ymin = tau_q025,ymax = tau_q975),alpha  = 0.3)+
  geom_line()+ 
  theme_bw()+ 
  scale_y_continuous(name = expression(tau[60*"-"*79 * "," ~ t]), breaks = scales::pretty_breaks(n=8))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank())

ggtau_oldest = final_summary %>% 
  filter(age == "80+") %>% 
  group_by(earliest_week_end_date) %>% 
  ggplot(aes(earliest_week_end_date, tau_med))+ 
  geom_ribbon(aes(ymin = tau_q025,ymax = tau_q975),alpha  = 0.3)+
  geom_line()+ 
  theme_bw()+ 
  scale_y_continuous(name = expression(tau[80*"+"* "," ~ t]), breaks = scales::pretty_breaks(n=8))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x.top = element_text(vjust = -66),
        axis.ticks.x.top = element_blank())


ggsave(filename = paste0("./age epidemic curve Toronto with admissions/plots/tau_older_joint.pdf"),
       plot = ggtau_oldest+ theme(axis.text.x.top = element_text(vjust = -78)),
       device = "pdf",
       width = 8/2,
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto with admissions/plots/tau_young_joint.pdf"),
       plot = ggtau_young+ theme(axis.text.x.top = element_text(vjust = -78)),
       device = "pdf",
       width = 8/2,
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto with admissions/plots/tau_adult_joint.pdf"),
       plot = ggtau_adult+ theme(axis.text.x.top = element_text(vjust = -78)),
       device = "pdf",
       width = 8/2,
       height = 8/3,
       dpi = 300)

ggsave(filename = paste0("./age epidemic curve Toronto with admissions/plots/tau_old_joint.pdf"),
       plot = ggtau_old+theme(axis.text.x.top = element_text(vjust = -78)),
       device = "pdf",
       width = 8/2,
       height = 8/3,
       dpi = 300)

round(final_summary$tau_mean_med[1:4],2)
round(final_summary$tau_mean_q025[1:4],2)
round(final_summary$tau_mean_q975[1:4],2)

round(final_summary$rho_tau_med[1],2)
round(final_summary$rho_tau_q025[1],2)
round(final_summary$rho_tau_q975[1],2)
