library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
sapply(list.files(path = "functions/voi-master/R/", pattern = "\\.R$", full.names = TRUE), source)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(ggpattern)
source("functions/add_scenario_cols.R")
source("functions/order_scenario_cols.R")
source("functions/calc_cost_pc.R")
source("functions/find_best_strats.R")
source("functions/calc_cost_bdown.R")
source("functions/parse_inputs.R")
source("functions/voi_dec.R")
source("functions/voi_fit.R")
source("functions/table_formatting.R")

file_list <- list.files(path = "../output/archetypes/", pattern = "\\.csv$", full.names = TRUE)
arch_data <- lapply(file_list, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>%
             calc_cost_pc() %>% parse_inputs()
arch_voi  <- voi_dec(arch_data, combn(c("mean_age", "sd_age", "skew_age", "le", "comm", "IGHC", "schoolA2", "workp", "EPOP", "TES", 
                                        "gvapw", "wfh", "Tres", "t_tit", "trate", "sdb", "sdl", "Hmax", "t_vax", "arate", "puptake"), 
                                      1, simplify = FALSE)) %>%
             group_by(location, disease) %>% slice_max(order_by = res, n = 1) %>% ungroup()
arch_data <- arch_data %>% left_join(arch_voi, by = c("location", "disease")) %>%
             rowwise() %>% mutate(xaxis = get(parameter)) %>% ungroup()
arch_fit  <- voi_fit(arch_data) %>% mutate(l = pmax(0,l)) %>% rename(xaxis = x, SECpc = y, lower = l, upper = u) %>%
             group_by(location, disease, xaxis) %>% mutate(alpha = ifelse(SECpc == min(SECpc), 1, 0)) %>% ungroup() %>%  
             mutate(group = cumsum(alpha == 0 | lag(alpha == 0, default = FALSE))) %>% 
             group_by(location, disease, strategy) %>%
             mutate(alpha = ifelse(any(alpha == 1) & alpha == 0, 0.25, alpha)) %>%
             ungroup()
arch_data <- arch_data %>% group_by(location, disease, strategy, xaxis) %>%
             mutate(alpha = {l <- location
                             d <- disease
                             s <- strategy
                             x <- xaxis
                             arch_fit %>% filter(location == l, disease == d, strategy == s) %>%
                                          slice(which.min(abs(xaxis - x))) %>% pull(alpha)}) %>% 
             mutate(alpha = floor(alpha)) %>% ungroup() 
arch_lab  <- arch_fit %>% group_by(disease, location) %>% 
             summarize(xlabel = unique(parameter), .groups = 'drop') %>%
             mutate(xlabel = case_when(xlabel == "mean_age" ~ "Mean Age (years)", 
                                       xlabel == "sd_age" ~ "SD Age (years)", 
                                       xlabel == "skew_age" ~ "Skewness Age", 
                                       xlabel == "le" ~ "Life Expectancy (years)",
                                       xlabel == "comm" ~ "Mean Household Contacts (per person/day)", 
                                       xlabel == "IGHC" ~ "Inter-Generation Household Contacts (%)",
                                       xlabel == "schoolA2" ~ "Mean School Contacts (per person/day)", 
                                       xlabel == "workp" ~ "Mean Workplace Contacts (per person/day)",
                                       xlabel == "EPOP" ~ "Employment-Population Ratio (%)", 
                                       xlabel == "TES" ~ "Teritary Sector Employment (%)",
                                       xlabel == "gvapw" ~ "GVA per Worker ($m)", 
                                       xlabel == "wfh" ~ "Home-Working Ratio (%)",
                                       xlabel == "Tres" ~ "Response Time (day)", 
                                       xlabel == "t_tit" ~ "Testing Start-Time (day)",
                                       xlabel == "trate" ~ "Testing Rate (per 100k/day)", 
                                       xlabel == "sdb" ~ "Distancing Sensitivity",
                                       xlabel == "sdl" ~ "Distancing Maximum", 
                                       xlabel == "Hmax" ~ "Hospital Capacity (per 100k)",
                                       xlabel == "t_vax" ~ "Vaccination Start-Time (day)", 
                                       xlabel == "arate" ~ "Vaccination Rate (per 100k/day)",
                                       xlabel == "puptake" ~ "Vaccination Coverage (%)"),  
                    xx = rep(c(0.16, 0.49, 0.81), 7),
                    xy = rep(seq(0.838, 0, by = -0.1375), each = 3))

gg <- ggplot(data = arch_data, aes(x = xaxis, y = SECpc, color = strategy, alpha = alpha)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_ribbon(data = arch_fit %>% filter(alpha == 1), 
                  aes(ymin = lower, ymax = upper, fill = strategy, group = interaction(group, strategy)), color = NA, alpha = 0.25) +
      geom_point(shape = 19, size = 0.2, stroke = 0.25) +
      geom_line(data = arch_fit, linewidth = 0.7) +
      scale_color_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                    "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_fill_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_alpha(range = c(0, 1)) +     
      theme_bw() +
      scale_x_log10(expand = c(0, 0), position = "bottom") +
      scale_y_continuous(expand = c(0, 0), position = "right") +
      #facetted_pos_scales(x = list(disease == "Covid-Omicron-X" & location == "HIC" ~ scale_x_continuous(labels = "Hospital Capacity"))) +   
      theme(panel.spacing = unit(1.25, "lines")) +
      labs(title = "", x = "", y = "Societal Loss (% of GDP)") +
      guides(color = guide_legend(title = NULL), fill = "none", alpha = "none") + 
      theme(legend.position = c(0.14, 0.9995), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 7))
gg <- ggdraw(gg) + draw_text(arch_lab$xlabel, x = arch_lab$xx, y = arch_lab$xy, hjust = 0.5, vjust = 0.5, size = 9)

ggsave("figure_3.png", plot = gg, height = 14, width = 10)