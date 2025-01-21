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
             calc_cost_pc() %>% parse_inputs()
arch_voi  <- voi_dec(arch_data, combn(c("mean_age", "sd_age", "skew_age", "le", "comm", "IGHC", "schoolA2", "workp", "EPOP", "TES", 
                                        "gvapw", "wfh", "Tres", "t_tit", "trate", "sdb", "sdl", "Hmax", "t_vax", "arate", "puptake"), 
                                      2, simplify = FALSE)) %>%
             group_by(location, disease) %>% slice_max(order_by = res, n = 1) %>% ungroup() #%>%
             #mutate(param1    = str_split(parameter, ",", simplify = TRUE)[, 1],
             #       param2    = str_split(parameter, ",", simplify = TRUE)[, 2],
             #       parameter = case_when(row_number() %in% c(1,2,5,6,8,9,10,11,12,13,16,18,21) ~ paste(param2, param1, sep = ","),
             #                             TRUE ~ parameter)) %>%
             #select(-param1,-param2)
arch_data <- arch_data %>% left_join(arch_voi, by = c("location", "disease")) 
arch_fit  <- voi_fit(arch_data) %>% rename(xaxis = x, yaxis = y, SECpc = z) %>% 
             group_by(location, disease, xaxis, yaxis) %>%
             mutate(SECpc = ifelse(SECpc == min(SECpc, na.rm = TRUE), SECpc, NA)) %>%
             ungroup()
arch_lab  <- arch_fit %>% group_by(disease, location) %>% 
             summarize(xlabel = str_split(unique(parameter), ",", simplify = TRUE)[, 1],
                       ylabel = str_split(unique(parameter), ",", simplify = TRUE)[, 2], .groups = 'drop') %>%
             mutate(across(c(xlabel, ylabel), ~ case_when(. == "mean_age" ~ "Mean Age (years)", 
                                                          . == "sd_age" ~ "SD Age (years)", 
                                                          . == "skew_age" ~ "Skewness Age", 
                                                          . == "le" ~ "Life Expectancy (years)",
                                                          . == "comm" ~ "Mean Household Contacts (per person/day)", 
                                                          . == "IGHC" ~ "Inter-Generation Household Contacts (%)",
                                                          . == "schoolA2" ~ "Mean School Contacts (per person/day)", 
                                                          . == "workp" ~ "Mean Workplace Contacts (per person/day)",
                                                          . == "EPOP" ~ "Employment-Population Ratio (%)", 
                                                          . == "TES" ~ "Teritary Sector Employment (%)",
                                                          . == "gvapw" ~ "GVA per Worker ($m)", 
                                                          . == "wfh" ~ "Home-Working Ratio (%)",
                                                          . == "Tres" ~ "Response Time (day)", 
                                                          . == "t_tit" ~ "Testing Start-Time (day)",
                                                          . == "trate" ~ "Testing Rate (per 100k/day)", 
                                                          . == "sdb" ~ "Distancing Sensitivity",
                                                          . == "sdl" ~ "Distancing Maximum", 
                                                          . == "Hmax" ~ "Hospital Capacity (per 100k)",
                                                          . == "t_vax" ~ "Vaccination Start-Time (day)", 
                                                          . == "arate" ~ "Vaccination Rate (per 100k/day)",
                                                          . == "puptake" ~ "Vaccination Coverage (%)")),
                    xx = rep(c(0.16, 0.49, 0.82), 7),
                    xy = rep(seq(0.838, 0, by = -0.1375), each = 3),
                    yx = rep(c(0.325, 0.655, 0.985), 7), 
                    yy = rep(seq(0.900, 0, by = -0.1375), each = 3))

gg <- ggplot(arch_fit, aes(x = xaxis, y = yaxis, z = SECpc, fill = strategy)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_contour_filled() +
      scale_fill_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      theme_bw() +
      scale_x_log10(expand = c(0, 0), position = "bottom") +
      scale_y_log10(expand = c(0, 0), position = "right") +
      #facetted_pos_scales(x = list(disease == "Covid-Omicron-X" & location == "HIC" ~ scale_x_continuous(labels = "Hospital Capacity"))) +   
      theme(panel.spacing = unit(1.25, "lines")) +
      labs(title = "", x = "", y = "") +
      guides(fill = "none") + 
      theme(legend.position = "right")
    #      legend.position = c(0.29, 0.999), legend.justification = c(1, 1), legend.box.just = "right", 
    #      legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 7))
gg <- ggdraw(gg) + draw_text(arch_lab$xlabel, x = arch_lab$xx, y = arch_lab$xy, hjust = 0.5, vjust = 0.5, size = 9) + 
      draw_text(arch_lab$ylabel, angle = 90, x = arch_lab$yx, y = arch_lab$yy, hjust = 0.5, vjust = 0.5, size = 9)

ggsave("figure_7.png", plot = gg, height = 14, width = 10)

#this script is not done!
#colorbar
#independent cost axes? possibly split figure by row and have colorbar per disease ...
#how to indicate optimal strategy? different colours/patterns?