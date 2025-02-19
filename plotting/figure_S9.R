library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(fitdistrplus)
library(forecast)
sapply(list.files(path = "functions/voi-master/R/", pattern = "\\.R$", full.names = TRUE), source)
library(scam)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(ggpattern)
library(patchwork)
source("functions/add_scenario_cols.R")
source("functions/order_scenario_cols.R")
source("functions/calc_loss_pc.R")
source("functions/parse_inputs.R")
source("functions/voi_decision.R")
source("functions/voi_fit.R")
#source("functions/table_formatting.R")

list_files   <- list.files(path = "../output/archetypes/", pattern = "\\.csv$", full.names = TRUE)
input_files  <- list_files[grepl("_data\\.csv$", list_files)]
input_data   <- lapply(input_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>% parse_inputs() #slightly slow
output_files <- list_files[!grepl("_data\\.csv$", list_files)]
output_data  <- lapply(output_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>% 
                left_join(input_data %>% dplyr::select(location, country, mean_age, pr_le, AL_wavg, AHT_wavg, AS_wavg, workp,
                                                       epop, pr_workf, pr_gvapw, pr_wfh, Tres, sda, sdb, sdc, t_tit, trate, 
                                                       Hmax, t_vax, arate, puptake), by = c("location", "country")) %>%
                (function(x) calc_loss_pc(input_data,x)) %>% 
                dplyr::select(-gdp, -tdur, -starts_with("vlyl"), -starts_with("gdpl"), -vsyl, -VLYLpc, -GDPLpc, -VSYLpc)
evppi_values <- voi_decision(output_data, as.list(setdiff(names(output_data), c("location", "country", "disease", "strategy", "SLpc")))) %>%
                #filter(!parameter %in% c("sda","sdb","sdc")) %>% #slow              
                group_by(location, disease) %>% 
                slice_max(order_by = res, n = 1) %>% 
                ungroup()
output_data  <- output_data %>% left_join(evppi_values, by = c("location", "disease")) %>%
                rowwise() %>% 
                mutate(xaxis = get(parameter)) %>% 
                ungroup()
evppi_fits   <- voi_fit(output_data) %>% rename(xaxis = x, SLpc = y, lower = l, upper = u) %>% 
                left_join(evppi_values, by = c("location", "disease", "parameter"))
#prettification: vaccine coverage rescaling, gva per worker rescaling, alpha, filtering, xlabels
output_data  <- output_data %>% left_join(input_data %>% filter(country == 1) %>% 
                                          dplyr::select(location, ifr_1918, ifr_dlta, ifr_sars), by = c("location")) %>%
                mutate(xaxis = case_when(disease == "Influenza-1918-X" & parameter == "puptake" ~ xaxis*4^log10(100*ifr_1918),
                                         disease == "Covid-Delta-X" & parameter == "puptake" ~ xaxis*4^log10(100*ifr_dlta),
                                         disease == "SARS-X" & parameter == "puptake" ~ xaxis*4^log10(100*ifr_sars),
                                         TRUE ~ xaxis)) %>%
                mutate(xaxis = ifelse(parameter == "puptake", 100*xaxis, xaxis)) %>%
                mutate(xaxis = ifelse(parameter == "pr_gvapw", 10^6*xaxis, xaxis))
evppi_fits   <- evppi_fits %>% left_join(input_data %>% filter(country == 1) %>%
                                         dplyr::select(location, ifr_1918, ifr_dlta, ifr_sars), by = c("location")) %>%
                mutate(xaxis = case_when(disease == "Influenza-1918-X" & parameter == "puptake" ~ xaxis*4^log10(100*ifr_1918),
                                         disease == "Covid-Delta-X" & parameter == "puptake" ~ xaxis*4^log10(100*ifr_dlta),
                                         disease == "SARS-X" & parameter == "puptake" ~ xaxis*4^log10(100*ifr_sars),
                                         TRUE ~ xaxis)) %>%
                mutate(xaxis = ifelse(parameter == "puptake", 100*xaxis, xaxis)) %>%
                mutate(xaxis = ifelse(parameter == "pr_gvapw", 10^6*xaxis, xaxis))
evppi_fits   <- evppi_fits %>%
                group_by(location, disease, xaxis) %>% 
                mutate(alpha = ifelse(SLpc == min(SLpc), 1, 0.99)) %>% 
                ungroup() %>%
                group_by(location, disease, strategy) %>%
                mutate(alpha = ifelse(!any(alpha == 1) & alpha == 0.99, 0, alpha)) %>%
                ungroup() %>%
                group_by(location, disease) %>%
                mutate(alpha = ifelse(alpha > 0 & res < 0.25, 0.25, alpha)) %>% #threshold features in geom_ribbon below
                ungroup()
evppi_fitsX  <- evppi_fits %>%
                group_by(location, disease) %>% 
                summarise(nxaxis = unique(xaxis), .groups = "drop")
output_dataX <- output_data %>% 
                group_by(location, disease) %>% 
                summarise(xaxis = unique(xaxis), .groups = "drop") %>%
                left_join(evppi_fitsX, by = c("location", "disease")) %>%
                group_by(location, disease, xaxis) %>% 
                slice(which.min(abs(xaxis - nxaxis)))
output_data  <- output_data %>%
                left_join(output_dataX, by = c("location", "disease", "xaxis")) %>%
                left_join(evppi_fits %>% rename(nxaxis = xaxis) %>% dplyr::select(location, disease, strategy, nxaxis, alpha),
                          by = c("location", "disease", "strategy", "nxaxis")) %>%
                mutate(alpha = floor(alpha))
output_data  <- output_data %>%
                filter(!(location == "LLMIC" & disease == "Influenza-2009-X" & strategy == "Economic Closures") &
                       !(location == "LLMIC" & disease == "Covid-Omicron-X" & strategy == "Elimination") &
                       !(location == "LLMIC" & disease == "Covid-Delta-X" & strategy == "Elimination") &
                       !(location == "UMIC" & disease == "Influenza-2009-X" & strategy == "Economic Closures") &
                       !(location == "UMIC" & disease == "Influenza-1918-X" & strategy == "No Closures") &
                       !(location == "UMIC" & disease == "Covid-Omicron-X" & strategy == "Elimination") &
                       !(location == "UMIC" & disease == "SARS-X" & strategy == "School Closures") &
                       !(location == "UMIC" & disease == "SARS-X" & strategy == "Elimination") &
                       !(location == "HIC" & disease == "Influenza-2009-X" & strategy == "Economic Closures") &
                       !(location == "HIC" & disease == "Covid-Delta-X" & strategy == "School Closures") &
                       !(location == "HIC" & disease == "Covid-Delta-X" & strategy == "Economic Closures"))
evppi_fits   <- evppi_fits %>%
                filter(!(location == "LLMIC" & disease == "Influenza-2009-X" & strategy == "Economic Closures") &
                       !(location == "LLMIC" & disease == "Covid-Omicron-X" & strategy == "Elimination") &
                       !(location == "LLMIC" & disease == "Covid-Delta-X" & strategy == "Elimination") &
                       !(location == "UMIC" & disease == "Influenza-2009-X" & strategy == "Economic Closures") &
                       !(location == "UMIC" & disease == "Influenza-1918-X" & strategy == "No Closures") &
                       !(location == "UMIC" & disease == "Covid-Omicron-X" & strategy == "Elimination") &
                       !(location == "UMIC" & disease == "SARS-X" & strategy == "School Closures") &
                       !(location == "UMIC" & disease == "SARS-X" & strategy == "Elimination") &
                       !(location == "HIC" & disease == "Influenza-2009-X" & strategy == "Economic Closures") &
                       !(location == "HIC" & disease == "Covid-Delta-X" & strategy == "School Closures") &
                       !(location == "HIC" & disease == "Covid-Delta-X" & strategy == "Economic Closures")) 
evppi_labs   <- evppi_fits %>% 
                group_by(disease, location) %>%
                summarize(xlabel = unique(parameter), .groups = "drop") %>% #may vary by income group(@)
                mutate(xlabel = case_when(xlabel == "mean_age" ~ "Average Age (years)",
                                          #xlabel == "pr_le" ~ "@ Life Expectancy (years)",
                                          #xlabel == "AL_wavg" ~ "Household Contacts (#/person/day)",
                                          #xlabel == "AHT_wavg" ~ "Other-Location Contacts (#/person/day)",
                                          xlabel == "AS_wavg" ~ "School Contacts (#/child/day)",
                                          #xlabel == "workp" ~ "Workplace Contacts (#/adult/day)",
                                          #xlabel == "epop" ~ "Employment-Population Ratio (%)",
                                          #xlabel == "pr_workf" ~ "@ Sector Workforce (%)",
                                          xlabel == "pr_gvapw" ~ "GVA per Worker ($)",
                                          #xlabel == "pr_wfh" ~ "@ Home-Working Ratio (%)",
                                          xlabel == "Tres" ~ "Response Time (doubling times)",
                                          xlabel == "sda" ~ "Distancing Multiplier Intercept",
                                          xlabel == "sdb" ~ "Distancing Death-Sensitivity",
                                          xlabel == "sdc" ~ "Distancing Time-Decay",
                                          xlabel == "t_tit" ~ "Testing Start-Time (doubling times)",
                                          xlabel == "trate" ~ "Testing Rate (per 100k/day)",
                                          xlabel == "Hmax" ~ "Spare Hospital Beds (per 100k)",
                                          xlabel == "t_vax" ~ "Vaccination Start-Time (days)",
                                          #xlabel == "arate" ~ "Vaccination Rate (per 100k/day)",
                                          xlabel == "puptake" ~ "Vaccination Coverage (% of HIT)"),
                       xx = rep(c(0.17, 0.49, 0.81), 7),
                       xy = rep(seq(0.8372, 0.01, by = -0.1317), each = 3))

gg <- ggplot(output_data, aes(x = xaxis, y = SLpc, color = strategy, alpha = alpha)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_ribbon(data = evppi_fits, 
                  aes(ymin = lower, ymax = upper, fill = strategy, alpha = ceiling(alpha)*(0.05+0.45*(res > 0.25))), 
                  color = NA) +
      geom_point(shape = 19, size = 0.02, stroke = 0.2) +
      geom_line(data = evppi_fits, linewidth = 0.5) +
      scale_color_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy",
                                    "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_fill_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy",
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_alpha(range = c(0, 1)) +
      theme_bw() +
      facetted_pos_scales(
      x = list(
      location == "LLMIC" & disease == "Influenza-2009-X" ~ scale_x_log10(limits=c(1,300),breaks=c(1,3,10,30,100,300),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Influenza-2009-X" ~ scale_x_log10(limits=c(3,300),breaks=c(3,10,30,100,300),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Influenza-2009-X" ~ scale_x_log10(limits=c(10,1000),breaks=c(10,30,100,300,1000),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Influenza-1957-X" ~ scale_x_continuous(limits=c(-4,4),breaks=seq(-4,4,by=2),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Influenza-1957-X" ~ scale_x_continuous(limits=c(-5,5),breaks=seq(-5,5,by=2.5),labels=scales::label_parse(),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Influenza-1957-X" ~ scale_x_continuous(limits=c(-3,3),breaks=seq(-3,3,by=1.5),labels=scales::label_parse(),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Influenza-1918-X" ~ scale_x_log10(limits=c(1,300),breaks=c(1,3,10,30,100,300),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Influenza-1918-X" ~ scale_x_continuous(limits=c(-5,5),breaks=seq(-5,5,by=2.5),labels=scales::label_parse(),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Influenza-1918-X" ~ scale_x_continuous(limits=c(50,200),breaks=seq(50,200,by=50),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Covid-Omicron-X" ~ scale_x_continuous(limits=c(-4,4),breaks=seq(-4,4,by=2),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Covid-Omicron-X" ~ scale_x_continuous(limits=c(-5,5),breaks=seq(-5,5,by=2.5),labels=scales::label_parse(),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Covid-Omicron-X" ~ scale_x_log10(limits=c(1e-3,1e+1),breaks=c(1e-3,1e-2,1e-1,1e+0,1e+1),labels = scales::trans_format("log10", scales::math_format(10^.x)),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Covid-Delta-X" ~ scale_x_continuous(limits=c(-4,4),breaks=seq(-4,4,by=2),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Covid-Delta-X" ~ scale_x_log10(limits=c(1e-3,1e-1),breaks=c(1e-3,3e-3,1e-2,3e-2,1e-1),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Covid-Delta-X" ~ scale_x_log10(limits=c(1e-3,1e+1),breaks=c(1e-3,1e-2,1e-1,1e+0,1e+1),labels = scales::trans_format("log10", scales::math_format(10^.x)),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Covid-Wildtype-X" ~ scale_x_continuous(limits=c(-4,4),breaks=seq(-4,4,by=2),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Covid-Wildtype-X" ~ scale_x_log10(limits=c(1e-3,1e+1),breaks=c(1e-3,1e-2,1e-1,1e+0,1e+1),labels = scales::trans_format("log10", scales::math_format(10^.x)),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Covid-Wildtype-X" ~ scale_x_log10(limits=c(1e-3,1e+1),breaks=c(1e-3,1e-2,1e-1,1e+0,1e+1),labels = scales::trans_format("log10", scales::math_format(10^.x)),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "SARS-X" ~ scale_x_continuous(limits=c(0,250),breaks=seq(0,250,by=50),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "SARS-X" ~ scale_x_log10(limits=c(1,100),breaks=c(1,3,10,30,100),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "SARS-X" ~ scale_x_continuous(limits=c(4,14),breaks=seq(4,14,by=2),expand=c(0,0),position="bottom")),
      y = list(
      location == "LLMIC" & disease == "Influenza-2009-X" ~ scale_y_log10(limits=c(0.3,300),breaks=3*10^seq(-1,2),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Influenza-2009-X" ~ scale_y_log10(limits=c(0.3,300),breaks=3*10^seq(-1,2),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Influenza-2009-X" ~ scale_y_log10(limits=c(0.3,300),breaks=3*10^seq(-1,2),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Influenza-1957-X" ~ scale_y_log10(limits=c(0.3,300),breaks=3*10^seq(-1,2),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Influenza-1957-X" ~ scale_y_log10(limits=c(0.3,300),breaks=3*10^seq(-1,2),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Influenza-1957-X" ~ scale_y_log10(limits=c(0.3,300),breaks=3*10^seq(-1,2),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Influenza-1918-X" ~ scale_y_log10(limits=c(100,3000),breaks=c(100,300,1000,3000),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Influenza-1918-X" ~ scale_y_log10(limits=c(100,3000),breaks=c(100,300,1000,3000),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Influenza-1918-X" ~ scale_y_log10(limits=c(30,1000),breaks=c(30,100,300,1000),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Covid-Omicron-X" ~ scale_y_log10(limits=c(10,300),breaks=c(10,30,100,300),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Covid-Omicron-X" ~ scale_y_log10(limits=c(10,300),breaks=c(10,30,100,300),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Covid-Omicron-X" ~ scale_y_log10(limits=c(10,300),breaks=c(10,30,100,300),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Covid-Delta-X" ~ scale_y_log10(limits=c(30,1000),breaks=c(30,100,300,1000),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Covid-Delta-X" ~ scale_y_log10(limits=c(30,1000),breaks=c(30,100,300,1000),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Covid-Delta-X" ~ scale_y_log10(limits=c(30,1000),breaks=c(30,100,300,1000),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Covid-Wildtype-X" ~ scale_y_log10(limits=c(10,300),breaks=c(10,30,100,300),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Covid-Wildtype-X" ~ scale_y_log10(limits=c(10,300),breaks=c(10,30,100,300),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Covid-Wildtype-X" ~ scale_y_log10(limits=c(10,300),breaks=c(10,30,100,300),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "SARS-X" ~ scale_y_log10(limits=c(100,3000),breaks=c(100,300,1000,3000),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "SARS-X" ~ scale_y_log10(limits=c(30,1000),breaks=c(30,100,300,1000),labels=scales::label_parse(),expand=c(0,0),position="right"),
      location == "HIC" & disease == "SARS-X" ~ scale_y_log10(limits=c(30,1000),breaks=c(30,100,300,1000),labels=scales::label_parse(),expand=c(0,0),position="right"))) + 
      theme(panel.spacing = unit(0.75, "lines")) +
      labs(title = "", x = "", y = "Societal Loss (% of GDP)") +
      guides(color = guide_legend(title = NULL, override.aes = list(linewidth = 1.5)), fill = "none", alpha = "none") +
      theme(legend.position = "bottom", legend.box.just = "right", 
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9), 
            legend.margin = margin(0, 0, 0, 0))
gg <- ggdraw(gg) + draw_text(evppi_labs$xlabel, x = evppi_labs$xx, y = evppi_labs$xy, hjust = 0.5, vjust = 0.5, size = 9)

ggsave("figure_S9.png", plot = gg, height = 14, width = 10)