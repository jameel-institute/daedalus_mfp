library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
sapply(list.files(path = "functions/voi-master/R/", pattern = "\\.R$", full.names = TRUE), source)
library(ggplot2)
library(ggh4x)
library(cowplot)
source("functions/add_scenario_cols.R")
source("functions/order_scenario_cols.R")
source("functions/calc_cost_pc.R")
source("functions/find_best_strats.R")
source("functions/calc_cost_bdown.R")
source("functions/parse_inputs.R")
source("functions/voi_dec.R")
source("functions/voi_fit.R")

file_list <- list.files(path = "../output/archetypes/", pattern = "\\.csv$", full.names = TRUE)
arch_data <- lapply(file_list, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>%
             calc_cost_pc() %>% find_best_strats() %>% filter(min_mean == TRUE) %>% parse_inputs()
arch_voi  <- voi_est(arch_data, combn(c("Tres","t_tit","trate","sdb","sdl","Hmax","t_vax","arate","puptake"), 1, simplify = FALSE)) %>%
             group_by(location, disease, strategy) %>% slice_max(order_by = res, n = 1) %>% ungroup()
arch_data <- arch_data %>% left_join(arch_voi, by = c("location", "disease", "strategy")) %>%
             rowwise() %>% mutate(xaxis = get(parameter)) %>% ungroup()
arch_fit  <- voi_fit(arch_data) %>% mutate(l = pmax(0,l)) %>% rename(xaxis = x, SECpc = y, lower = l, upper = u)
arch_lab  <- arch_fit %>% group_by(disease, location, strategy) %>% 
             summarize(xlabel = unique(parameter), .groups = 'drop') %>%
             mutate(xlabel = case_when(xlabel == "Tres" ~ "Response Time (day)",
                                       xlabel == "t_tit" ~ "Testing Start-Time (day)",
                                       xlabel == "trate" ~ "Testing Rate (per 100k/day)",
                                       xlabel == "sdb" ~ "Distancing Sensitivity",
                                       xlabel == "sdl" ~ "Distancing Maximum",
                                       xlabel == "Hmax" ~ "Hospital Capacity (per 100k)",
                                       xlabel == "t_vax" ~ "Vaccination Start-Time (day)",
                                       xlabel == "arate" ~ "Vaccination Rate (per 100k/day)",
                                       xlabel == "puptake" ~ "Vaccination Coverage (%)"),  
                    xx = rep(c(0.16, 0.49, 0.82), 7),
                    xy = rep(seq(0.838, 0, by = -0.1375), each = 3))

gg <- ggplot(data = arch_data, aes(x = xaxis, y = SECpc, shape = strategy, color = parameter)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_ribbon(data = arch_fit, aes(ymin = lower, ymax = upper), color = NA, fill = "grey80", alpha = 0.5) +
      geom_point(size = 1.25, stroke = 0.35, fill = NA) +
      geom_line(data = arch_fit, linewidth = 0.75, color = "black") +
      scale_shape_manual(values = c("No Closures" = 21, "School Closures" = 23, "Economic Closures" = 22, "Elimination" = 25)) +
      scale_color_manual(values = c("t_tit" = "forestgreen", "trate" = "palegreen2", "sdb" = "slategray3", "Hmax" = "purple")) +
      theme_bw() +
      scale_x_log10(expand = c(0, 0), position = "bottom") +
      scale_y_continuous(expand = c(0, 0), position = "right") +
      #facetted_pos_scales(x = list(disease == "Covid-Omicron-X" & location == "HIC" ~ scale_x_continuous(labels = "Hospital Capacity"))) +   
      theme(panel.spacing = unit(1.25, "lines")) +
      labs(title = "", x = "", y = "Societal Loss (% of GDP)") +
      guides(shape = guide_legend(title = NULL), color = "none") + 
      theme(legend.position = c(0.29, 0.999), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 7))
gg <- ggdraw(gg) + draw_text(arch_lab$xlabel, x = arch_lab$xx, y = arch_lab$xy, hjust = 0.5, vjust = 0.5, size = 9)

ggsave("figure_3.png", plot = gg, height = 14, width = 10)