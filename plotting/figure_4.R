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
#source("functions/calc_loss_pc.R")
#source("functions/voi_dec.R")
#source("functions/voi_fit.R")
#source("functions/table_formatting.R")

list_files   <- list.files(path = "../output/archetypes/", pattern = "\\.csv$", full.names = TRUE)
input_files  <- list_files[grepl("_data\\.csv$", list_files)]
input_data   <- lapply(input_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>% parse_inputs() #slightly slow
output_files <- list_files[!grepl("_data\\.csv$", list_files)]
output_data  <- lapply(output_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>%
                left_join(input_data %>% dplyr::select(location, country, gdp, mean_vly), by = c("location", "country")) %>%
                mutate(x = rowSums(across(starts_with("vlyl_"))), 
                       y = rowSums(across(starts_with("gdpl_"))) + vsyl) %>%
                dplyr::select(-tdur, -starts_with("vlyl"), -vsyl, -starts_with("gdpl")) %>%
                crossing(mean_vly_range = 10^seq(-3, log10(0.3), length.out = 100)) %>%
                mutate(x = x*(mean_vly_range / mean_vly)) %>%
                mutate(SLpc = 100*(x + y)/gdp) %>%
                group_by(location, disease, mean_vly_range, country) %>% 
                summarise(strategy = paste(strategy[SLpc == min(SLpc)], collapse = ", "), 
                          mean_vly = unique(mean_vly)) %>% #slightly slow
                group_by(location, disease, mean_vly_range, mean_vly) %>% #mean_vly added as a group for preservation
                count(strategy) %>% 
                mutate(proportion = n / sum(n)) %>%
                mutate(strategy = case_when(strategy %in% c("No Closures", "School Closures", "Economic Closures", "Elimination") ~ strategy,
                                            strategy == "School Closures, Economic Closures" ~ "Untriggered Closures",
                                            strategy == "School Closures, Economic Closures, Elimination" ~ "Untriggered Closures",
                                            TRUE ~ "Other")) %>%
                group_by(location, disease, mean_vly_range, strategy) %>%
                summarise(proportion = sum(proportion), 
                          mean_vly   = unique(mean_vly)) %>%
                mutate(strategy = factor(strategy, levels = c("No Closures", "Untriggered Closures", "School Closures", 
                                                              "Economic Closures", "Elimination", "Other")))

gg <- ggplot(output_data, aes(x = mean_vly_range, y = proportion, fill = strategy, pattern_density = strategy)) +
      facet_grid2(disease ~ location, switch = "y", scales = "fixed") +
      geom_area_pattern(pattern = "stripe", pattern_color = "darkgreen", pattern_fill = "darkgreen") +
      geom_vline(aes(xintercept = mean_vly), linetype = "dashed", color = "white") +
      scale_fill_manual(values = c("No Closures" = "magenta4", "Untriggered Closures" = "navy", "School Closures" = "navy",
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod",  "Other" = "grey"))+  
      scale_pattern_density_manual(values = c("No Closures" = 0, "Untriggered Closures" = 0.2, "School Closures" = 0,
                                              "Economic Closures" = 0, "Elimination" = 0, "Other" = 0),
                                   breaks = c("No Closures", "Untriggered Closures", "School Closures", "Economic Closures", "Elimination")) +
      theme_bw() +
      scale_x_log10(breaks = c(0.001,0.003,0.01,0.03,0.1,0.3), expand = c(0,0), position = "bottom", labels = scales::label_number(scale = 1000000, suffix = "")) +
      scale_y_continuous(expand = c(0,0), position = "right") +
      theme(panel.spacing = unit(1.00, "lines"), axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "", x = "Average VLY ($, nominal)", y = "Proportion Loss-Minimising") +
      guides(fill = "none", pattern_density = guide_legend(title = NULL, override.aes = list(fill = c("magenta4","navy","navy","darkgreen","goldenrod")))) +
      theme(legend.position = "top", legend.box.just = "right", legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

ggsave("figure_4n.png", plot = gg, height = 14, width = 10)