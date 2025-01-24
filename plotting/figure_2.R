library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(fitdistrplus)
library(forecast)
sapply(list.files(path = "functions/voi-master/R/", pattern = "\\.R$", full.names = TRUE), source)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(ggpattern)
library(patchwork)
source("functions/add_scenario_cols.R")
source("functions/order_scenario_cols.R")
source("functions/parse_inputs.R")
source("functions/calc_loss_pc.R")
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
                group_by(location, disease) %>% 
                slice_max(order_by = res, n = 1) %>% 
                ungroup()
output_data  <- output_data %>% left_join(evppi_values, by = c("location", "disease")) %>%
                rowwise() %>% 
                mutate(xaxis = get(parameter)) %>% 
                ungroup()
evppi_fits   <- voi_fit(output_data) %>% rename(xaxis = x, SLpc = y, lower = l, upper = u)
#prettification: bound x, define alpha, bound y, labels
output_data  <- output_data %>%
                group_by(location, disease) %>%
                mutate(bound_xl = quantile(xaxis, 0.05),
                       bound_xu = quantile(xaxis, 0.95)) %>%
                ungroup() %>%
                filter(xaxis >= bound_xl & xaxis <= bound_xu)
evppi_fits   <- evppi_fits %>%
                left_join(output_data %>% dplyr::select(location, disease, bound_xl, bound_xu) %>% 
                          distinct(location, disease, .keep_all = TRUE), by = c("location", "disease")) %>%
                filter(xaxis >= bound_xl & xaxis <= bound_xu) %>%
                group_by(location, disease, xaxis) %>% 
                mutate(alpha = ifelse(SLpc == min(SLpc), 1, 0.50)) %>% 
                ungroup() %>%
                mutate(group = cumsum(alpha == 0.50 | lag(alpha == 0.50, default = FALSE))) %>%
                group_by(location, disease, strategy) %>%
                mutate(alpha = ifelse(!any(alpha == 1) & alpha == 0.50, 0, alpha)) %>%
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
                mutate(alpha = 0.5*floor(alpha))
evppi_fits   <- evppi_fits %>%
                filter(alpha > 0)
evppi_fitsY  <- evppi_fits %>%
                group_by(location, disease) %>% 
                summarise(lower = min(lower),
                          upper = max(upper),, .groups = "drop")
output_data  <- output_data %>%
                left_join(evppi_fitsY, by = c("location", "disease")) %>%
                group_by(location, disease) %>%
                mutate(bound_yl = quantile(SLpc[alpha > 0], 0.05),
                       bound_yu = quantile(SLpc[alpha > 0], 0.95)) %>%
                mutate(bound_yl = min(bound_yl, lower),
                       bound_yu = max(bound_yu, upper)) %>%
                ungroup() %>%
                filter(SLpc >= bound_yl & SLpc <= bound_yu)
evppi_labs   <- evppi_fits %>% 
                group_by(disease, location) %>%
                summarize(xlabel = unique(parameter), .groups = "drop") %>% #may vary by income group(@)
                mutate(xlabel = case_when(xlabel == "mean_age" ~ "Average Age (years)",
                                          xlabel == "pr_le" ~ "@ Life Expectancy (years)",
                                          xlabel == "AL_wavg" ~ "Household Contacts (#/person/day)",
                                          xlabel == "AHT_wavg" ~ "Other-Location Contacts (#/person/day)",
                                          xlabel == "AS_wavg" ~ "School Contacts (#/child/day)",
                                          xlabel == "workp" ~ "Workplace Contacts (#/adult/day)",
                                          xlabel == "epop" ~ "Employment-Population Ratio (%)",
                                          xlabel == "pr_workf" ~ "@ Sector Workforce (%)",
                                          xlabel == "pr_gvapw" ~ "@ GVA per Worker ($)",
                                          xlabel == "pr_wfh" ~ "@ Home-Working Ratio (%)",
                                          xlabel == "Tres" ~ "Response Time (doubling times)",
                                          xlabel == "sda" ~ "Transmission Multiplier Intercept",
                                          xlabel == "sdb" ~ "Transmission Death-Sensitivity",
                                          xlabel == "sdc" ~ "Transmission Time-Relaxation",
                                          xlabel == "t_tit" ~ "Testing Start-Time (doubling times)",
                                          xlabel == "trate" ~ "Testing Rate (per 100k/day)",
                                          xlabel == "Hmax" ~ "Spare Hospital Capacity (per 100k)",
                                          xlabel == "t_vax" ~ "Vaccination Start-Time (days)",
                                          xlabel == "arate" ~ "Vaccination Rate (per 100k/day)",
                                          xlabel == "puptake" ~ "Vaccination Coverage (% of HIT)"),
                       xx = rep(c(0.16, 0.49, 0.81), 7),
                       xy = rep(seq(0.80, 0.01, by = -0.131), each = 3))

gg <- ggplot(output_data, aes(x = xaxis, y = SLpc, color = strategy, alpha = alpha)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_ribbon(data = evppi_fits, 
                  aes(ymin = lower, ymax = upper, fill = strategy, group = interaction(group, strategy)), 
                  color = NA, alpha = 0.50) +
      geom_point(shape = 19, size = 0.10, stroke = 0.25) +
      geom_line(data = evppi_fits, linewidth = 0.5) +
      scale_color_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy",
                                    "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_fill_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy",
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_alpha(range = c(0, 1)) +
      theme_bw() +
      scale_x_log10(expand = c(0, 0), position = "bottom") +
      scale_y_log10(expand = c(0, 0), position = "right") +
      #facetted_pos_scales(x = list(disease == "Covid-Omicron-X" & location == "HIC" ~ scale_x_continuous(labels = "Hospital Capacity"))) +
      theme(panel.spacing = unit(1.25, "lines")) +
      labs(title = "", x = "", y = "Societal Loss (% of GDP)") +
      guides(color = guide_legend(title = NULL), fill = "none", alpha = "none") +
      theme(legend.position = "top", legend.box.just = "right", legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))
gg <- ggdraw(gg) + draw_text(evppi_labs$xlabel, x = evppi_labs$xx, y = evppi_labs$xy, hjust = 0.5, vjust = 0.5, size = 9)

ggsave("figure_2n.png", plot = gg, height = 14, width = 10)