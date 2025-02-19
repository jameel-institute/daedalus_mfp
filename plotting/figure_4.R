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
library(ggdensity)
library(cowplot)
library(ggpattern)
library(patchwork)
source("functions/add_scenario_cols.R")
source("functions/order_scenario_cols.R")
#source("functions/calc_loss_pc.R")
source("functions/parse_inputs.R")
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
                crossing(mean_vly_range = seq(0, 10, length.out = 100)) %>%
                mutate(x = x*(mean_vly_range / mean_vly)) %>%
                mutate(SLpc = 100*(x + y)/gdp) %>%
                group_by(location, disease, mean_vly_range, country) %>% 
                summarise(strategy = paste(strategy[SLpc == min(SLpc)], collapse = ", "), 
                          mean_vly = unique(mean_vly)) %>% #slightly slow
                group_by(location, disease, mean_vly_range, mean_vly) %>% #mean_vly added as a group for preservation
                count(strategy) %>% 
                mutate(proportion = n / sum(n)) %>%
                mutate(strategy = case_when(strategy %in% c("No Closures", "School Closures", "Economic Closures", "Elimination") ~ strategy,
                                            strategy == "No Closures, School Closures, Economic Closures" ~ "Untriggered Closures",
                                            strategy == "School Closures, Economic Closures" ~ "Untriggered Closures",
                                            TRUE ~ "Other")) %>%
                group_by(location, disease, mean_vly_range, strategy) %>%
                summarise(proportion = sum(proportion), 
                          mean_vly   = unique(mean_vly)) %>%
                mutate(strategy = factor(strategy, levels = c("Untriggered Closures", "No Closures", "School Closures", 
                                                              "Economic Closures", "Elimination", "Other")))

gg <- ggplot(output_data, aes(x = mean_vly_range, y = proportion, fill = strategy, pattern_density = strategy)) +
      facet_grid2(disease ~ location, switch = "y", scales = "fixed") +
      geom_area_pattern(pattern = "stripe", pattern_size = 0.25, pattern_color = "navy", pattern_fill = "darkgreen") +
      geom_vline(aes(xintercept = mean_vly), linetype = "dashed", color = "white") +
      scale_fill_manual(values = c("Untriggered Closures" = "magenta4", "No Closures" = "magenta4", "School Closures" = "navy",
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod",  "Other" = "grey"))+  
      scale_pattern_density_manual(values = c("Untriggered Closures" = 0.2, "No Closures" = 0, "School Closures" = 0,
                                              "Economic Closures" = 0, "Elimination" = 0, "Other" = 0),
                                   breaks = c("Untriggered Closures", "No Closures", "School Closures", "Economic Closures", "Elimination")) +
      theme_bw() +
      scale_x_continuous(breaks=seq(0,10,by=2),   expand=c(0,0), position="bottom") + 
      scale_y_continuous(breaks=seq(0,1,by=0.25), expand=c(0,0), position="right") + 
      theme(panel.spacing = unit(0.75, "lines")) +
      labs(title = "", x = "Average VLY / GDP per capita", y = "Proportion Societal Loss-Minimising") +
      guides(fill = "none", 
             pattern_density = guide_legend(title = NULL, 
                                            override.aes = list(pattern_size = 0.75, 
                                                                fill = c("magenta4","magenta4","navy","darkgreen","goldenrod")))) +
      theme(legend.position = "bottom", legend.box.just = "right", 
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9), 
            legend.margin = margin(0, 0, 0, 0))

ggsave("figure_4.png", plot = gg, height = 14, width = 10)