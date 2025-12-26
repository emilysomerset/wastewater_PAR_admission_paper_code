rm(list=ls())

# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(aghq) # aghq_0.4.1
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(gridExtra) # gridExtra_2.3
library(cowplot) # cowplot_1.1.3
library(ggplot2)
library(lubridate)

load("./section_3/wastewater_model_CAN_phachist.RData")
source("./functions_general/process_results_original.R")

weights_datadic <- data.frame(site_id = df_full$site_id %>% unique()) %>% 
  left_join(
        data.frame(site_id = c("TAB","THC","THU","TNT"), 
                             weights = c(0.499/(0.177+0.232+0.064+0.499),
                                        0.177/(0.177+0.232+0.064+0.499),
                                        0.232/(0.177+0.232+0.064+0.499),
                                        0.064/(0.177+0.232+0.064+0.499))), by = "site_id")

weights_datadic <- weights_datadic %>% 
  mutate(weights = ifelse(is.na(weights),0, weights))


results <- process_results(df_full =df_full,
                           tmbdat = tmbdat,
                           samps1 = samps1,
                           polyOrder=3,
                           id_group = 1,
                           weights_df = weights_datadic)


##############################################
gg1 <- results$df_full %>% 
  ggplot(aes(x = sample_date, exp_v)) +
  geom_line(size=0.2) + 
  # geom_line(data=results$df_full, aes(x= sample_date, exp_v), col = "red")+
  geom_ribbon(aes(ymax = exp_v_upr, ymin = exp_v_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste(V,"(t)")), breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_blank(),
        axis.text.x.top = element_text(vjust = -65),
        axis.ticks.x.top = element_blank())


gg2 <- results$station_ave_df %>% 
  ggplot(aes(x = sample_date, ave_exp_v_u_fixed_med)) +
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = ave_exp_v_u_fixed_upr, ymin = ave_exp_v_u_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5), limits = c(0,650))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_blank(),
        axis.text.x.top = element_text(vjust = -65),
        axis.ticks.x.top = element_blank())

gg3 <- results$station_ave_df %>% 
  ggplot(aes(x = sample_date, ave_exp_v_u_fixed_deriv_med)) +
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = ave_exp_v_u_fixed_deriv_upr, ymin = ave_exp_v_u_fixed_deriv_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste(bar(mu),"'(t)")), breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_blank(),
        axis.text.x.top = element_text(vjust = -65),
        axis.ticks.x.top = element_blank())

gg4 <- results$df_full %>% 
  ggplot(aes(x = sample_date, exp_v_deriv)) +
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = exp_v_deriv_upr, ymin = exp_v_deriv_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste(V,"'(t)")), breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_blank(),
        axis.text.x.top = element_text(vjust = -65),
        axis.ticks.x.top = element_blank())



a = gg1+ theme(axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_3/plots/allsignals_CAN_wastewatermodel_a.pdf"),
       plot = a, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)
rstudioapi::viewer(paste0("./section_3/plots/allsignals_CAN_wastewatermodel_a.pdf"))

b = gg3+ theme(axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_3/plots/allsignals_CAN_wastewatermodel_b.pdf"),
       plot = b, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)

c = gg2+ theme(axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_3/plots/allsignals_CAN_wastewatermodel_c.pdf"),
       plot = c, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)

d = gg4+ theme(axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_3/plots/allsignals_CAN_wastewatermodel_d.pdf"),
       plot = d, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)


