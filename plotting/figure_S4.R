library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(fitdistrplus)
sapply(list.files(path = "functions/voi-master/R/", pattern = "\\.R$", full.names = TRUE), source)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(ggpattern)
library(patchwork)
source("functions/add_scenario_cols.R")
source("functions/order_scenario_cols.R")
source("functions/calc_cost_pc.R")
source("functions/find_best_strats.R")
source("functions/calc_cost_bdown.R")
source("functions/parse_inputs.R")
source("functions/voi_dec.R")
source("functions/voi_fit.R")
source("functions/table_formatting.R")

#extract fitted parameters from mmpx_dist dataframes
# print(rbind(mmp1_dist %>% group_by(igroup,variable) %>% slice_min(xvalue,n=1), 
#             mmp2_dist %>% group_by(igroup,variable) %>% slice_min(xvalue,n=1)) %>% 
#       filter(igroup == "LLMIC") %>% dplyr::select(-xmin,-xmax,-fit,-aic,-xvalue,-yvalue,-alpha) %>% 
#       mutate(across(where(is.numeric), ~ sprintf("%.4f", .x))), n = 30)

ctry_data <- read.csv("../input/country_data.csv") %>%
             mutate(igroup = factor(igroup, levels = c("LLMIC","UMIC","HIC")))

mmp1_data <- ctry_data %>% dplyr::select(igroup,Tres,sda,sdb,sdc,t_tit) %>%
             pivot_longer(cols = c(Tres,sda,sdb,sdc,t_tit), names_to = "variable", values_to = "value") %>%
             mutate(variable = case_when(variable == "Tres" ~ "Distancing: Response Time",
                                         variable == "sda" ~ "Distancing: Multiplier Intercept",
                                         variable == "sdb" ~ "Distancing: Death-Sensitivity",
                                         variable == "sdc" ~ "Distancing: Time-Relaxation",
                                         variable == "t_tit" ~ "Surveillance: Testing Start-Time")) %>%
             mutate(variable = factor(variable, levels = c("Distancing: Response Time", "Distancing: Multiplier Intercept", 
                                                           "Distancing: Death-Sensitivity", "Distancing: Time-Relaxation", 
                                                           "Surveillance: Testing Start-Time")))
mmp1_dist <- mmp1_data %>% filter(!is.na(value)) %>%
             mutate(candidate = case_when(variable == "Distancing: Response Time"        ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Distancing: Multiplier Intercept" ~ list(c("norm")), 
                                          variable == "Distancing: Death-Sensitivity"    ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Distancing: Time-Relaxation"      ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Surveillance: Testing Start-Time" ~ list(c("lnorm", "gamma", "weibull")))) %>%
             unnest_longer(candidate) %>%
             group_by(igroup, variable, candidate) %>%
             summarize(xmin     = min(value),
                       xmax     = max(value),
                       fit      = list(tryCatch(fitdist(value, unique(candidate), method = "mle"), error = function(e) NA)),
                       mean     = sapply(fit, function(f) tryCatch(f$estimate["mean"], error = function(e) NA)),
                       sd       = sapply(fit, function(f) tryCatch(f$estimate["sd"], error = function(e) NA)),
                       meanlog  = sapply(fit, function(f) tryCatch(f$estimate["meanlog"], error = function(e) NA)),
                       sdlog    = sapply(fit, function(f) tryCatch(f$estimate["sdlog"], error = function(e) NA)),
                       shape    = sapply(fit, function(f) tryCatch(f$estimate["shape"], error = function(e) NA)),
                       rate     = sapply(fit, function(f) tryCatch(f$estimate["rate"], error = function(e) NA)),
                       scale    = sapply(fit, function(f) tryCatch(f$estimate["scale"], error = function(e) NA)),
                       shape1   = sapply(fit, function(f) tryCatch(f$estimate["shape1"], error = function(e) NA)),
                       shape2   = sapply(fit, function(f) tryCatch(f$estimate["shape2"], error = function(e) NA)),
                       aic      = sapply(fit, function(f) tryCatch(f$aic, error = function(e) NA)), .groups = "drop") %>%
             group_by(igroup, variable) %>%
             mutate(xvalue = map2(unique(xmin), unique(xmax), ~seq(.x, .y, length.out = 500))) %>% unnest(cols = c(xvalue)) %>%
             ungroup() %>%
             mutate(yvalue = case_when(candidate == "norm"    ~ dnorm(xvalue, mean = mean, sd = sd),
                                       candidate == "lnorm"   ~ dlnorm(xvalue, meanlog = meanlog, sdlog = sdlog),
                                       candidate == "gamma"   ~ dgamma(xvalue, shape = shape, rate = rate),
                                       candidate == "weibull" ~ dweibull(xvalue, shape = shape, scale = scale),
                                       candidate == "beta"    ~ dbeta(xvalue, shape1 = shape1, shape2 = shape2))) %>%
             group_by(igroup, variable) %>%
             mutate(alpha = if_else(aic == min(aic, na.rm = TRUE), 1, 0.1)) %>%
             ungroup() %>%
             filter(alpha == 1)
### for plot legibility 
mmp1_dist <- mmp1_dist %>% mutate(yvalue = ifelse(variable == "Distancing: Death-Sensitivity" & xvalue < 0.02 |
                                                  igroup == "HIC" & variable == "Distancing: Time-Relaxation" & xvalue < 0.0001, 
                                                  NA, yvalue))
###
mmp1_labs <- data.frame(x      = rep(0.49, 5), 
                        y      = seq(0.79, 0, by = -0.195), 
                        xlabel = c("Seeding-to-Reponse Delay (doubling times)",
                                   "Transmission Multiplier Intercept Coefficient",
                                   "Transmission Multiplier Death-Sensitivity Coefficient",
                                   "Transmission Multiplier Time-Relaxation Coefficient",
                                   "Seeding-to-Testing-Start Delay (doubling times)"))
mmp1_leg  <- mmp1_dist %>% group_by(variable, igroup) %>%
             slice_min(xvalue, n = 1) %>%
             summarise(mu        = case_when(candidate == "norm"    ~ mean,
                                             candidate == "lnorm"   ~ exp(meanlog + (sdlog^2 / 2)),
                                             candidate == "gamma"   ~ shape / rate,
                                             candidate == "weibull" ~ scale * gamma(1 + (1/shape)),
                                             candidate == "beta"    ~ shape1 / (shape1 + shape2)),
                       candidate = case_when(candidate == "norm"    ~ "Normal",
                                             candidate == "lnorm"   ~ "Lognormal",
                                             candidate == "gamma"   ~ "Gamma",
                                             candidate == "weibull" ~ "Weibull",
                                             candidate == "beta"    ~ "Beta"), .groups = "drop") %>%
             mutate(x      = rep(c(0.26, 0.58, 0.90), 5), 
                    y      = rep(seq(0.94, 0.10, by = -0.195), each = 3),
                    legend = paste(candidate, "\n μ =", round(mu,3)))

mmp2_data <- ctry_data %>% dplyr::select(igroup,trate,Hmax,t_vax,arate,puptake) %>%
             pivot_longer(cols = c(trate,Hmax,t_vax,arate,puptake), names_to = "variable", values_to = "value") %>%
             mutate(variable = case_when(variable == "trate" ~ "Surveillance: Testing Rate",
                                         variable == "Hmax" ~ "Healthcare: Hospital Capacity",
                                         variable == "t_vax" ~ "Vaccination: Administration Start-Time",
                                         variable == "arate" ~ "Vaccination: Administration Rate",
                                         variable == "puptake" ~ "Vaccination: Coverage")) %>%
             mutate(variable = factor(variable, levels = c("Surveillance: Testing Rate", "Healthcare: Hospital Capacity", 
                                                           "Vaccination: Administration Start-Time", "Vaccination: Administration Rate", 
                                                           "Vaccination: Coverage")))
mmp2_dist <- mmp2_data %>% filter(!is.na(value)) %>%
             mutate(candidate = case_when(variable == "Surveillance: Testing Rate"             ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Healthcare: Hospital Capacity"          ~ list(c("lnorm", "gamma", "weibull")), 
                                          variable == "Vaccination: Administration Start-Time" ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Vaccination: Administration Rate"       ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Vaccination: Coverage"                  ~ list(c("lnorm", "gamma", "weibull")))) %>%
             unnest_longer(candidate) %>%
             group_by(igroup, variable, candidate) %>%
             summarize(xmin     = min(value),
                       xmax     = max(value),
                       fit      = list(tryCatch(fitdist(value, unique(candidate), method = "mle"), error = function(e) NA)),
                       mean     = sapply(fit, function(f) tryCatch(f$estimate["mean"], error = function(e) NA)),
                       sd       = sapply(fit, function(f) tryCatch(f$estimate["sd"], error = function(e) NA)),
                       meanlog  = sapply(fit, function(f) tryCatch(f$estimate["meanlog"], error = function(e) NA)),
                       sdlog    = sapply(fit, function(f) tryCatch(f$estimate["sdlog"], error = function(e) NA)),
                       shape    = sapply(fit, function(f) tryCatch(f$estimate["shape"], error = function(e) NA)),
                       rate     = sapply(fit, function(f) tryCatch(f$estimate["rate"], error = function(e) NA)),
                       scale    = sapply(fit, function(f) tryCatch(f$estimate["scale"], error = function(e) NA)),
                       shape1   = sapply(fit, function(f) tryCatch(f$estimate["shape1"], error = function(e) NA)),
                       shape2   = sapply(fit, function(f) tryCatch(f$estimate["shape2"], error = function(e) NA)),
                       aic      = sapply(fit, function(f) tryCatch(f$aic, error = function(e) NA)), .groups = "drop") %>%
             group_by(igroup, variable) %>%
             mutate(xvalue = map2(unique(xmin), unique(xmax), ~seq(.x, .y, length.out = 500))) %>% unnest(cols = c(xvalue)) %>%
             ungroup() %>%
             mutate(yvalue = case_when(candidate == "norm"    ~ dnorm(xvalue, mean = mean, sd = sd),
                                       candidate == "lnorm"   ~ dlnorm(xvalue, meanlog = meanlog, sdlog = sdlog),
                                       candidate == "gamma"   ~ dgamma(xvalue, shape = shape, rate = rate),
                                       candidate == "weibull" ~ dweibull(xvalue, shape = shape, scale = scale),
                                       candidate == "beta"    ~ dbeta(xvalue, shape1 = shape1, shape2 = shape2))) %>%
             group_by(igroup, variable) %>%
             mutate(alpha = if_else(aic == min(aic, na.rm = TRUE), 1, 0.1)) %>%
             ungroup() %>%
             filter(alpha == 1)
mmp2_labs <- data.frame(x      = rep(0.49, 5), 
                        y      = seq(0.79, 0, by = -0.195), 
                        xlabel = c("Tests Administered Daily (per 100k/day)",
                                   "Spare Hospital Beds (per 100k)",
                                   "Seeding-to-Vaccination-Start Delay (days)",
                                   "Vaccines Administered Daily (per 100k/day)",
                                   "Coverage Relative to Herd-Immunity Threshold (conditioned on IFR = 1%) (%)"))
mmp2_leg  <- mmp2_dist %>% group_by(variable, igroup) %>%
             slice_min(xvalue, n = 1) %>%
             summarise(mu        = case_when(candidate == "norm"    ~ mean,
                                             candidate == "lnorm"   ~ exp(meanlog + (sdlog^2 / 2)),
                                             candidate == "gamma"   ~ shape / rate,
                                             candidate == "weibull" ~ scale * gamma(1 + (1/shape)), 
                                             candidate == "beta"    ~ shape1 / (shape1 + shape2)),
                       candidate = case_when(candidate == "norm"    ~ "Normal",
                                             candidate == "lnorm"   ~ "Lognormal",
                                             candidate == "gamma"   ~ "Gamma",
                                             candidate == "weibull" ~ "Weibull",
                                             candidate == "beta"    ~ "Beta"), .groups = "drop") %>%
             mutate(mu     = ifelse(variable == "Vaccination: Coverage", 100*mu, mu),
                    x      = rep(c(0.25, 0.57, 0.89), 5), 
                    y      = rep(seq(0.94, 0.10, by = -0.195), each = 3),
                    legend = paste(candidate, "\n μ =", round(mu,3)))

p1 <- ggplot(mmp1_data, aes(x = value, fill = variable)) +
      facet_grid2(variable ~ igroup, switch = "y", scales = "free", independent = "all") +
      geom_histogram(aes(y = ..density..), color = "black") + 
      geom_line(data = mmp1_dist, aes(x = xvalue, y = yvalue), linewidth = 1, color = "black") +
      scale_fill_manual(values = c("Distancing: Response Time" = "slategray2", "Distancing: Multiplier Intercept" = "slategray2", 
                                   "Distancing: Death-Sensitivity" = "slategray2", "Distancing: Time-Relaxation" = "slategray2", 
                                   "Surveillance: Testing Start-Time"  = "palegreen")) + 
      theme_bw() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right") +
      theme(panel.spacing.y = unit(2.25, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "Relative Frequency")
p1 <- ggdraw(p1) + draw_text(mmp1_labs$xlabel, x = mmp1_labs$x, y = mmp1_labs$y, hjust = 0.5, vjust = 0.5, size = 9)
p1 <- ggdraw(p1) + draw_text(mmp1_leg$legend, x = mmp1_leg$x, y = mmp1_leg$y, hjust = 0.5, vjust = 0.5, size = 9)

p2 <- ggplot(mmp2_data, aes(x = value, fill = variable)) +
      facet_grid2(variable ~ igroup, switch = "y", scales = "free", independent = "all") +
      geom_histogram(aes(y = ..density..), color = "black") + 
      geom_line(data = mmp2_dist, aes(x = xvalue, y = yvalue), linewidth = 1, color = "black") +
      scale_fill_manual(values = c("Surveillance: Testing Rate" = "palegreen", "Healthcare: Hospital Capacity" = "purple", 
                                   "Vaccination: Administration Start-Time" = "lightsalmon", "Vaccination: Administration Rate" = "lightsalmon", 
                                   "Vaccination: Coverage" = "lightsalmon")) +
      theme_bw() +
      facetted_pos_scales(x = list(variable == "Vaccination: Coverage" ~ scale_x_continuous(labels = scales::label_number(scale = 100, suffix = "")))) +   
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right") +
      theme(panel.spacing.y = unit(2.25, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "Relative Frequency")
p2 <- ggdraw(p2) + draw_text(mmp2_labs$xlabel, x = mmp2_labs$x, y = mmp2_labs$y, hjust = 0.5, vjust = 0.5, size = 9)
p2 <- ggdraw(p2) + draw_text(mmp2_leg$legend, x = mmp2_leg$x, y = mmp2_leg$y, hjust = 0.5, vjust = 0.5, size = 9)

ggsave("figure_S4a.png", plot = p1, height = 14, width = 10)
ggsave("figure_S4b.png", plot = p2, height = 14, width = 10)