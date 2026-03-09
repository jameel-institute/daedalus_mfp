library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(fBasics)
library(fitdistrplus)
library(forecast)
library(scam)
library(ggplot2)
library(ggh4x)
library(ggdensity)
library(cowplot)
library(ggpattern)
library(patchwork)
source("functions/add_scenario_cols.R")
source("functions/order_scenario_cols.R")
source("functions/calc_loss_pc.R")
source("functions/parse_inputs.R")
source("functions/format_table.R")

#extract fitted parameters for marginal distributions from mmpx_dist dataframes
# print(rbind(mmp1_dist %>% group_by(igroup,variable) %>% slice_min(xvalue,n=1), 
#             mmp2_dist %>% group_by(igroup,variable) %>% slice_min(xvalue,n=1),
#             mmp3_dist %>% group_by(igroup,variable) %>% slice_min(xvalue,n=1)) %>% 
#       filter(igroup == "LLMIC") %>% dplyr::select(-xmin,-xmax,-fit,-aic,-xvalue,-yvalue,-alpha) %>% 
#       mutate(across(where(is.numeric), ~ sprintf("%.4f", .x))), n = 30)
#extract parameter correlations from mm_data dataframe
# test<- mm_data %>% filter(igroup == "LLMIC") %>% dplyr::select(-igroup,-country)
# SIG <- cor(test, use = "pairwise.complete.obs")
# write.csv(SIG, file = "SIG_correlation_matrix.csv", row.names = FALSE)

ctry_data <- read.csv("../input/country_data.csv") %>%
             mutate(igroup = factor(igroup, levels = c("LLMIC","UMIC","HIC")))
mmp1_data <- ctry_data %>% dplyr::select(igroup,country,Tres,sda,sdb,sdc) %>%
             pivot_longer(cols = c(Tres,sda,sdb,sdc), names_to = "variable", values_to = "value") %>%
             mutate(variable = case_when(variable == "Tres" ~ "Distancing: Response Time",
                                         variable == "sda" ~ "Distancing: Transmission Multiplier Intercept",
                                         variable == "sdb" ~ "Distancing: Transmission Multiplier Death-Sensitivity",
                                         variable == "sdc" ~ "Distancing: Transmission Multiplier Time-Decay")) %>%
             mutate(variable = factor(variable, levels = c("Distancing: Response Time", "Distancing: Transmission Multiplier Intercept", 
                                                           "Distancing: Transmission Multiplier Death-Sensitivity", "Distancing: Transmission Multiplier Time-Decay")))
mmp1_dist <- mmp1_data %>% filter(!is.na(value)) %>%
             mutate(candidate = case_when(variable == "Distancing: Response Time"        ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Distancing: Transmission Multiplier Intercept" ~ list(c("norm")), 
                                          variable == "Distancing: Transmission Multiplier Death-Sensitivity"    ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Distancing: Transmission Multiplier Time-Decay"           ~ list(c("lnorm", "gamma", "weibull")))) %>%
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
mmp1_dist <- mmp1_dist %>% mutate(yvalue = ifelse(variable == "Distancing: Transmission Multiplier Death-Sensitivity" & xvalue < 0.02 |
                                                  igroup == "HIC" & variable == "Distancing: Transmission Multiplier Time-Decay" & xvalue < 0.0001, 
                                                  NA, yvalue))
###
mmp1_labs <- data.frame(x      = rep(0.49, 4), 
                        y      = seq(0.742, 0.014, length.out = 4), 
                        xlabel = c("Response Time (doubling times)",
                                   "Transmission Multiplier Intercept Coefficient",
                                   "Transmission Multiplier Death-Sensitivity Coefficient",
                                   "Transmission Multiplier Time-Decay Coefficient"))
mmp1_leg  <- mmp1_dist %>% 
             group_by(variable, igroup) %>%
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
             mutate(x      = rep(c(0.255, 0.57, 0.885), 4), 
                    y      = rep(seq(0.93, 0.20, length.out = 4), each = 3),
                    legend = paste(candidate, "\n μ =", round(mu,3)))
mmp2_data <- ctry_data %>% dplyr::select(igroup,country,t_tit,trate,masc,Hmax) %>%
             pivot_longer(cols = c(t_tit,trate,masc,Hmax), names_to = "variable", values_to = "value") %>%
             mutate(variable = case_when(variable == "t_tit" ~ "Case Isolation: Testing Start-Time",
                                         variable == "trate" ~ "Case Isolation: Testing Rate",
                                         variable == "masc"  ~ "Case Isolation: Maximum Infection Ascertainment Ratio",
                                         variable == "Hmax"  ~ "Hospital Capacity: Spare Hospital Beds")) %>%
             mutate(variable = factor(variable, levels = c("Case Isolation: Testing Start-Time", "Case Isolation: Testing Rate", 
                                                           "Case Isolation: Maximum Infection Ascertainment Ratio", "Hospital Capacity: Spare Hospital Beds")))
mmp2_dist <- mmp2_data %>% filter(!is.na(value)) %>%
             mutate(candidate = case_when(variable == "Case Isolation: Testing Start-Time"  ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Case Isolation: Testing Rate"        ~ list(c("lnorm", "gamma", "weibull")), 
                                          variable == "Case Isolation: Maximum Infection Ascertainment Ratio" ~ list(c("beta")),
                                          variable == "Hospital Capacity: Spare Hospital Beds"     ~ list(c("lnorm", "gamma", "weibull")))) %>%
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
mmp2_labs <- data.frame(x      = rep(0.49, 4), 
                        y      = seq(0.742, 0.014, length.out = 4), 
                        xlabel = c("Testing Start-Time (doubling times)",
                                   "Testing Rate (per 100k/day)",
                                   "Maximum Infection Ascertainment Ratio (%)",
                                   "Spare Hospital Beds (per 100k)"))
mmp2_leg  <- mmp2_dist %>% 
             group_by(variable, igroup) %>%
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
             mutate(x      = rep(c(0.245, 0.56, 0.875), 4), 
                    y      = rep(seq(0.93, 0.20, length.out = 4), each = 3),
                    legend = paste(candidate, "\n μ =", round(mu,3)))
mmp3_data <- ctry_data %>% dplyr::select(igroup,country,t_vax,arate,puptake) %>%
             pivot_longer(cols = c(t_vax,arate,puptake), names_to = "variable", values_to = "value") %>%
             mutate(variable = case_when(variable == "t_vax" ~ "Vaccination: Administration Start-Time",
                                         variable == "arate" ~ "Vaccination: Administration Rate",
                                         variable == "puptake" ~ "Vaccination: Coverage")) %>%
             mutate(variable = factor(variable, levels = c("Vaccination: Administration Start-Time", 
                                                           "Vaccination: Administration Rate", 
                                                           "Vaccination: Coverage")))
mmp3_dist <- mmp3_data %>% filter(!is.na(value)) %>%
             mutate(candidate = case_when(variable == "Vaccination: Administration Start-Time" ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Vaccination: Administration Rate"       ~ list(c("lnorm", "gamma", "weibull")),
                                          variable == "Vaccination: Coverage"                  ~ list(c("beta")))) %>%
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
mmp3_labs <- data.frame(x      = rep(0.49, 3), 
                        y      = seq(0.66, 0.015, length.out = 3), 
                        xlabel = c("Administration Start-Time (days)",
                                   "Administration Rate (per 100k/day)",
                                   "Population-Level Vaccination Coverage (%)"))
mmp3_leg  <- mmp3_dist %>% 
             group_by(variable, igroup) %>%
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
             mutate(x      = rep(c(0.245, 0.56, 0.875), 3), 
                    y      = rep(seq(0.915, 0.27, length.out = 3), each = 3),
                    legend = paste(candidate, "\n μ =", round(mu,3)))

p1 <- ggplot(mmp1_data, aes(x = value, fill = variable)) +
      facet_grid2(variable ~ igroup, switch = "y", scales = "free", independent = "all",
                  labeller = labeller(variable = c("Distancing: Transmission Multiplier Death-Sensitivity" =
                                                   "Distancing:\nTransmission Multiplier Death-Sensitivity"))) +
      geom_histogram(aes(y = ..density..), color = "black") + 
      geom_line(data = mmp1_dist, aes(x = xvalue, y = yvalue), linewidth = 1, color = "black") +
      scale_fill_manual(values = c("Distancing: Response Time" = "slategray2", "Distancing: Transmission Multiplier Intercept" = "slategray2", 
                                   "Distancing: Transmission Multiplier Death-Sensitivity" = "slategray2", "Distancing: Transmission Multiplier Time-Decay" = "slategray2")) + 
      theme_bw() +
      scale_x_continuous(expand = c(0,0), position = "bottom", labels = scales::label_parse()) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right") +
      theme(panel.spacing.y = unit(2.25, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "Relative Frequency")
p1 <- ggdraw(p1) + draw_text(mmp1_labs$xlabel, x = mmp1_labs$x, y = mmp1_labs$y, hjust = 0.5, vjust = 0.5, size = 9)
p1 <- ggdraw(p1) + draw_text(mmp1_leg$legend, x = mmp1_leg$x, y = mmp1_leg$y, hjust = 0.5, vjust = 0.5, size = 9)

p2 <- ggplot(mmp2_data, aes(x = value, fill = variable)) +
      facet_grid2(variable ~ igroup, switch = "y", scales = "free", independent = "all",
                  labeller = labeller(variable = c("Case Isolation: Maximum Infection Ascertainment Ratio" =
                                                   "Case Isolation:\nMaximum Infection Ascertainment Ratio"))) +
      geom_histogram(aes(y = ..density..), color = "black") + 
      geom_line(data = mmp2_dist, aes(x = xvalue, y = yvalue), linewidth = 1, color = "black") +
      scale_fill_manual(values = c("Case Isolation: Testing Start-Time" = "palegreen", "Case Isolation: Testing Rate" = "palegreen", 
                                   "Case Isolation: Maximum Infection Ascertainment Ratio" = "palegreen", "Hospital Capacity: Spare Hospital Beds" = "purple")) +
      theme_bw() +
      scale_x_continuous(expand = c(0,0), position = "bottom", labels = scales::label_parse()) +
      facetted_pos_scales(x = list(variable == "Case Isolation: Maximum Infection Ascertainment Ratio" ~ scale_x_continuous(labels = scales::label_number(scale = 100, suffix = "")))) +   
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right") +
      theme(panel.spacing.y = unit(2.25, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "Relative Frequency")
p2 <- ggdraw(p2) + draw_text(mmp2_labs$xlabel, x = mmp2_labs$x, y = mmp2_labs$y, hjust = 0.5, vjust = 0.5, size = 9)
p2 <- ggdraw(p2) + draw_text(mmp2_leg$legend, x = mmp2_leg$x, y = mmp2_leg$y, hjust = 0.5, vjust = 0.5, size = 9)

p3 <- ggplot(mmp3_data, aes(x = value, fill = variable)) +
      facet_grid2(variable ~ igroup, switch = "y", scales = "free", independent = "all") +
      geom_histogram(aes(y = ..density..), color = "black") + 
      geom_line(data = mmp3_dist, aes(x = xvalue, y = yvalue), linewidth = 1, color = "black") +
      scale_fill_manual(values = c("Vaccination: Administration Start-Time" = "lightsalmon", 
                                   "Vaccination: Administration Rate" = "lightsalmon", 
                                   "Vaccination: Coverage" = "lightsalmon")) +
      theme_bw() +
      scale_x_continuous(expand = c(0,0), position = "bottom", labels = scales::label_parse()) +
      facetted_pos_scales(x = list(variable == "Vaccination: Coverage" ~ scale_x_continuous(labels = scales::label_number(scale = 100, suffix = "")))) +   
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right") +
      theme(panel.spacing.y = unit(2.25, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "Relative Frequency")
p3 <- ggdraw(p3) + draw_text(mmp3_labs$xlabel, x = mmp3_labs$x, y = mmp3_labs$y, hjust = 0.5, vjust = 0.5, size = 9)
p3 <- ggdraw(p3) + draw_text(mmp3_leg$legend, x = mmp3_leg$x, y = mmp3_leg$y, hjust = 0.5, vjust = 0.5, size = 9)

ggsave("figure_S2a.png", plot = p1, height = 14, width = 10)
ggsave("figure_S2b.png", plot = p2, height = 14, width = 10)
ggsave("figure_S2c.png", plot = p3, height = 11, width = 10)

###################################################################################################################################

mm_data <- rbind(mmp1_data,mmp2_data,mmp3_data)
mm_dist <- rbind(mmp1_dist,mmp2_dist,mmp3_dist) %>%
           group_by(igroup,variable) %>% slice_min(xvalue, n=1) %>% ungroup() %>%
           dplyr::select(-xmin,-xmax,-fit,-aic,-alpha,-xvalue,-yvalue)
mm_data <- mm_data %>%
           left_join(mm_dist, by = c("igroup", "variable")) %>%
           mutate(qu = case_when(candidate == "norm"    ~ pnorm(value, mean = mean, sd = sd),
                                 candidate == "lnorm"   ~ plnorm(value, meanlog = meanlog, sdlog = sdlog),
                                 candidate == "gamma"   ~ pgamma(value, shape = shape, rate = rate),
                                 candidate == "weibull" ~ pweibull(value, shape = shape, scale = scale),
                                 candidate == "beta"    ~ pbeta(value, shape1 = shape1, shape2 = shape2))) %>%
           mutate(zs = qnorm(qu, mean = 0, sd = 1)) %>%
           dplyr::select(igroup, country, variable, zs) %>%
           pivot_wider(names_from = variable, values_from = zs) 

p4 <- ggpairs(mm_data %>% filter(igroup == "LLMIC") %>% dplyr::select(-igroup,-country) %>% rename_with(~ gsub(": ", ":\n", .x)),
              lower = list(continuous = wrap("points", shape = 21, size = 1, stroke = 0.2, color = "black", fill = "orange")),
              diag  = list(continuous = wrap("blankDiag")),
              upper = list(continuous = wrap("cor", size = 3)),
              switch = "y") +
      theme_minimal() +
      scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 2), expand = c(0,0)) +
      scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 2), expand = c(0,0)) +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            strip.placement = "outside", 
            strip.text.x = element_text(angle = 90, hjust = 1),
            strip.text.y.left = element_text(angle = 0, hjust = 1),
            axis.title.y.left = element_blank()) + 
      labs(title = "", x = "Z-Score", y = "Z-Score")

p5 <- ggpairs(mm_data %>% filter(igroup == "UMIC") %>% dplyr::select(-igroup,-country) %>% rename_with(~ gsub(": ", ":\n", .x)),
              lower = list(continuous = wrap("points", shape = 21, size = 1, stroke = 0.2, color = "black", fill = "turquoise")),
              diag  = list(continuous = wrap("blankDiag")),
              upper = list(continuous = wrap("cor", size = 3)),
              switch = "y") +
      theme_minimal() +
      scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 2), expand = c(0,0)) +
      scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 2), expand = c(0,0)) +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            strip.placement = "outside", 
            strip.text.x = element_text(angle = 90, hjust = 1),
            strip.text.y.left = element_text(angle = 0, hjust = 1),
            axis.title.y.left = element_blank()) + 
      labs(title = "", x = "Z-Score", y = "Z-Score")

p6 <- ggpairs(mm_data %>% filter(igroup == "HIC") %>% dplyr::select(-igroup,-country) %>% rename_with(~ gsub(": ", ":\n", .x)),
              lower = list(continuous = wrap("points", shape = 21, size = 1, stroke = 0.2, color = "black", fill = "azure4")),
              diag  = list(continuous = wrap("blankDiag")),
              upper = list(continuous = wrap("cor", size = 3)),
              switch = "y") +
      theme_minimal() +
      scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 2), expand = c(0,0)) +
      scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 2), expand = c(0,0)) +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            strip.placement = "outside", 
            strip.text.x = element_text(angle = 90, hjust = 1),
            strip.text.y.left = element_text(angle = 0, hjust = 1),
            axis.title.y.left = element_blank()) + 
      labs(title = "", x = "Z-Score", y = "Z-Score")

ggsave("figure_S3a.png", plot = p4, height = 14, width = 10)
ggsave("figure_S3b.png", plot = p5, height = 14, width = 10)
ggsave("figure_S3c.png", plot = p6, height = 14, width = 10)

###################################################################################################################################

SIG_L <- mm_data %>% filter(igroup == "LLMIC") %>% dplyr::select(-igroup,-country) %>% 
         cor(, use = "pairwise.complete.obs")
SIG_M <- mm_data %>% filter(igroup == "UMIC") %>% dplyr::select(-igroup,-country) %>% 
         cor(, use = "pairwise.complete.obs")
SIG_H <- mm_data %>% filter(igroup == "HIC") %>% dplyr::select(-igroup,-country) %>% 
         cor(, use = "pairwise.complete.obs")
SIG   <- list(LLMIC = SIG_L, UMIC = SIG_M, HIC = SIG_H) 
mmd   <- mmp1_dist %>%
         group_by(igroup,variable) %>% slice_min(xvalue,n=1) %>% ungroup() %>%
         dplyr::select(-xmin,-xmax,-fit,-aic,-alpha,-xvalue,-yvalue) %>%
         group_by(igroup) %>%
         group_modify(~{
            vars <- .x$variable
            sig  <- SIG[[ .y$igroup ]][vars, vars]
            Z    <- MASS::mvrnorm(n  = 1000, mu = rep(0, length(vars)), Sigma = sig)
            ###
            as_tibble(Z) %>% setNames(vars) %>% mutate(sample = row_number()) %>%
            pivot_longer(-sample, names_to  = "variable", values_to = "Z") %>%
            left_join(.x, by = "variable")}) %>%
         ungroup() %>%
         mutate(U     = pnorm(Z),
                value = case_when(candidate == "norm"    ~ qnorm(U, mean = mean, sd = sd),
                                  candidate == "lnorm"   ~ qlnorm(U, meanlog = meanlog, sdlog = sdlog),
                                  candidate == "gamma"   ~ qgamma(U, shape = shape, rate = rate),
                                  candidate == "weibull" ~ qweibull(U, shape = shape, scale = scale),
                                  candidate == "beta"    ~ qbeta(U, shape1 = shape1, shape2 = shape2))) %>%
         dplyr::select(igroup,variable,sample,value) %>%
         pivot_wider(names_from = "variable", values_from = "value") %>%
         crossing(dvalue = 10^seq(-2,0), 
                  xvalue = seq(0, 1000, by = 250)) %>%
         mutate(yvalue = 1/(1 + exp(`Distancing: Transmission Multiplier Intercept` 
                                    + `Distancing: Transmission Multiplier Death-Sensitivity` * log10(dvalue) 
                                    - `Distancing: Transmission Multiplier Time-Decay` * (xvalue)))) %>%
         dplyr::select(igroup,sample,dvalue,xvalue,yvalue) %>%
         mutate(variable = "Distancing", 
                dvalue = factor(dvalue, levels = c("1","0.1","0.01")))
mms <-   mmp2_dist %>% filter(variable != "Hospital Capacity: Spare Hospital Beds") %>% 
         group_by(igroup,variable) %>% slice_min(xvalue,n=1) %>% ungroup() %>%
         dplyr::select(-xmin,-xmax,-fit,-aic,-alpha,-xvalue,-yvalue) %>%
         group_by(igroup) %>%
         group_modify(~{
           vars <- .x$variable
           sig  <- SIG[[ .y$igroup ]][vars, vars]
           Z    <- MASS::mvrnorm(n  = 100, mu = rep(0, length(vars)), Sigma = sig)
           ###
           as_tibble(Z) %>% setNames(vars) %>% mutate(sample = row_number()) %>%
           pivot_longer(-sample, names_to  = "variable", values_to = "Z") %>%
           left_join(.x, by = "variable")}) %>%
         ungroup() %>%
         mutate(U     = pnorm(Z),
                value = case_when(candidate == "norm"    ~ qnorm(U, mean = mean, sd = sd),
                                  candidate == "lnorm"   ~ qlnorm(U, meanlog = meanlog, sdlog = sdlog),
                                  candidate == "gamma"   ~ qgamma(U, shape = shape, rate = rate),
                                  candidate == "weibull" ~ qweibull(U, shape = shape, scale = scale),
                                  candidate == "beta"    ~ qbeta(U, shape1 = shape1, shape2 = shape2))) %>%  
         dplyr::select(igroup,variable,sample,value) %>%
         pivot_wider(names_from = "variable", values_from = "value") %>%
         crossing(xvalue = seq(0, 100)) %>%
         mutate(yvalue = `Case Isolation: Maximum Infection Ascertainment Ratio` * Heaviside(xvalue, a = `Case Isolation: Testing Start-Time`)) %>%
         dplyr::select(igroup,sample,xvalue,yvalue) %>%
         mutate(variable = "Case Isolation")
mmh <-   mmp2_dist %>% filter(variable == "Hospital Capacity: Spare Hospital Beds") %>% 
         group_by(igroup,variable) %>% slice_min(xvalue,n=1) %>% ungroup() %>%
         dplyr::select(-xmin,-xmax,-fit,-aic,-alpha,-xvalue,-yvalue) %>%
         crossing(sample = seq(1, 100)) %>%
         mutate(value = case_when(candidate == "norm"    ~ rnorm(n(), mean = mean, sd = sd),
                                  candidate == "lnorm"   ~ rlnorm(n(), meanlog = meanlog, sdlog = sdlog),
                                  candidate == "gamma"   ~ rgamma(n(), shape = shape, rate = rate),
                                  candidate == "weibull" ~ rweibull(n(), shape = shape, scale = scale),
                                  candidate == "beta"    ~ rbeta(n(), shape1 = shape1, shape2 = shape2))) %>%
         dplyr::select(igroup,variable,sample,value) %>%
         pivot_wider(names_from = "variable", values_from = "value") %>%
         crossing(xvalue = seq(0, 1000, by = 1000)) %>%
         mutate(yvalue = `Hospital Capacity: Spare Hospital Beds`) %>%
         dplyr::select(igroup,sample,xvalue,yvalue) %>%
         mutate(variable = "Hospital Capacity")
mmv <-   mmp3_dist %>% 
         group_by(igroup,variable) %>% slice_min(xvalue,n=1) %>% ungroup() %>%
         dplyr::select(-xmin,-xmax,-fit,-aic,-alpha,-xvalue,-yvalue) %>%
         group_by(igroup) %>%
         group_modify(~{
           vars <- .x$variable
           sig  <- SIG[[ .y$igroup ]][vars, vars]
           Z    <- MASS::mvrnorm(n  = 300, mu = rep(0, length(vars)), Sigma = sig)
           ###
           as_tibble(Z) %>% setNames(vars) %>% mutate(sample = row_number()) %>%
           pivot_longer(-sample, names_to  = "variable", values_to = "Z") %>%
           left_join(.x, by = "variable")}) %>%
         ungroup() %>%
         mutate(U     = pnorm(Z),
                value = case_when(candidate == "norm"    ~ qnorm(U, mean = mean, sd = sd),
                                  candidate == "lnorm"   ~ qlnorm(U, meanlog = meanlog, sdlog = sdlog),
                                  candidate == "gamma"   ~ qgamma(U, shape = shape, rate = rate),
                                  candidate == "weibull" ~ qweibull(U, shape = shape, scale = scale),
                                  candidate == "beta"    ~ qbeta(U, shape1 = shape1, shape2 = shape2))) %>%  
         dplyr::select(igroup,variable,sample,value) %>%
         pivot_wider(names_from = "variable", values_from = "value") %>%
         crossing(xvalue = seq(0, 1000)) %>%
         mutate(yvalue = pmin(`Vaccination: Coverage`, 
                         pmax(0, 10^-5*`Vaccination: Administration Rate` * (xvalue - `Vaccination: Administration Start-Time`)))) %>%
         dplyr::select(igroup,sample,xvalue,yvalue) %>%
         mutate(variable = "Vaccination")
mm  <-   bind_rows(mmd,mms,mmh,mmv) %>%
         mutate(variable =  factor(variable, levels = c("Distancing", "Case Isolation", "Hospital Capacity", "Vaccination")))
mm_labs <- data.frame(x      = rep(0.49, 4), 
                      y      = seq(0.742, 0.014, length.out = 4), 
                      xlabel = c("Time since Response Time (days)",
                                 "Time since Outbreak (doubling times)",
                                 "Time since Outbreak (days)",
                                 "Time since Outbreak (days)"))

p7 <- ggplot(mm, aes(x = xvalue, y = yvalue, color = variable, group = sample)) +
      facet_grid2(variable ~ igroup, switch = "y", scales = "free", independent = "x",
                  labeller = labeller(variable = c("Distancing" = "Distancing:\nTransmission Multiplier",
                                                   "Case Isolation" = "Case Isolation:\nMaximum Infection Ascertainment Ratio (%)",
                                                   "Hospital Capacity"  = "Hospital Capacity:\nSpare Hospital Beds (per 100k)", 
                                                   "Vaccination"   = "Vaccination:\nPopulation-Level Vaccination Coverage (%)"))) +
      geom_boxplot_pattern(data = subset(mm, variable == "Distancing"), 
                           aes(group = interaction(dvalue, xvalue), pattern = dvalue),   
                           position = position_dodge(width = 100), outlier.shape = NA, coef = 1.5, width = 80, 
                           linewidth = 0.2, color = "black", fill = "slategray2", 
                           pattern_density = 0.35, pattern_colour = "white", pattern_fill = "white") + 
      scale_pattern_manual(values = c("0.01" = "none", "0.1" = "stripe", "1" = "circle")) +
      geom_line(data = subset(mm, variable != "Distancing"), linewidth = 0.2) +
      scale_color_manual(values = c("Case Isolation" = "palegreen", "Hospital Capacity" = "purple", "Vaccination" = "lightsalmon")) +
      theme_bw() +
      scale_x_continuous(expand = c(0,0), position = "bottom", labels = scales::label_parse()) +
      facetted_pos_scales(y = list(variable == "Distancing" ~ scale_y_continuous(limits=c(0,1), expand = c(0,0), position = "right"),
                                   variable == "Case Isolation" ~ scale_y_continuous(limits=c(0,1), expand = c(0,0), labels = scales::label_number(scale = 100, suffix = ""), position = "right"),
                                   variable == "Hospital Capacity" ~ scale_y_log10(limits = c(3,300), expand = c(0,0), position = "right"),
                                   variable == "Vaccination" ~ scale_y_continuous(limits=c(0,1), expand = c(0,0), labels = scales::label_number(scale = 100, suffix = ""), position = "right"))) + 
      guides(color = "none", pattern = guide_legend(title = "Daily Deaths per 100k")) +
      theme(panel.spacing.x = unit(1.25, "lines"), panel.spacing.y = unit(2.25, "lines"), 
            legend.position = c(0.235, 0.838), legend.title = element_text(size = 9),  
            legend.key.width  = unit(0.6, "cm"), legend.key.height = unit(0.8, "cm")) + 
      labs(title = "", x = "", y = "")
p7 <- ggdraw(p7) + draw_text(mm_labs$xlabel, x = mm_labs$x, y = mm_labs$y, hjust = 0.5, vjust = 0.5, size = 9)

ggsave("figure_S4.png", plot = p7, height = 14, width = 10)