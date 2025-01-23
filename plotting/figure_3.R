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
source("functions/calc_loss_pc.R")
#source("functions/parse_inputs.R")
#source("functions/voi_dec.R")
#source("functions/voi_fit.R")
#source("functions/table_formatting.R")

list_files   <- list.files(path = "../output/archetypes/", pattern = "\\.csv$", full.names = TRUE)
input_files  <- list_files[grepl("_data\\.csv$", list_files)]
input_data   <- lapply(input_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols()
output_files <- list_files[!grepl("_data\\.csv$", list_files)]
output_data  <- lapply(output_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>%
                (function(x) calc_loss_pc(input_data,x)) %>% mutate(x = VLYLpc, y = GDPLpc + VSYLpc) %>%
                group_by(location, disease, strategy) %>%
                mutate(mean_x    = mean(x),
                       q1_x      = quantile(x, 0.25),
                       q3_x      = quantile(x, 0.75),
                       bound_x   = quantile(x, 0.95),
                       mean_y    = mean(y),
                       q1_y      = quantile(y, 0.25),
                       q3_y      = quantile(y, 0.75),
                       bound_y   = quantile(y, 0.95),
                       mean_SLpc = mean(SLpc)) %>%
                group_by(location, disease) %>%
                mutate(bound_x = max(bound_x), 
                       bound_y = max(bound_y)) %>%
                filter(x <= bound_x & y <= bound_y) %>%
                group_by(location, country, disease) %>% 
                mutate(min_strat = (SLpc == min(SLpc))) %>% #includes strategies with equal losses
                group_by(location, disease, strategy) %>% 
                mutate(alpha = 1 - abs(cor(x, y))) %>% 
                #   mutate(alpha = case_when(
                #            disease == "Influenza-2009-X" & strategy == "Elimination" & (GDPLpc+VSYLpc) < 50 ~ alpha/6,
                #            disease == "Influenza-1957-X" & strategy != "No Closures" ~ alpha/3,
                #            disease == "Influenza-1918-X" | 
                #            disease == "Covid-Delta-X" | 
                #            disease == "Covid-Wildtype-X" & strategy == "No Closures" ~ alpha/3,
                #            location == "HIC" & disease == "SARS-X" & strategy == "No Closures" ~ alpha*20,
                #            TRUE ~ alpha)) %>%
                mutate(alpha = 0.05 + alpha/10) %>%
                mutate(alpha = ifelse(min_strat, 1, alpha)) %>%
                ungroup()
output_stats <- output_data %>% #for quicker plotting
                group_by(location, disease, strategy) %>% 
                summarise(mean_x    = unique(mean_x), 
                          q1_x      = unique(q1_x), 
                          q3_x      = unique(q3_x),
                          mean_y    = unique(mean_y), 
                          q1_y      = unique(q1_y), 
                          q3_y      = unique(q3_y),
                          mean_SLpc = unique(mean_SLpc)) %>%
                group_by(location, disease) %>%
                mutate(min_mean  = (mean_SLpc == min(mean_SLpc)),
                       min_mean2 = {
                         strategy_means <- mean_SLpc[!duplicated(strategy)]
                         sorted_means   <- sort(strategy_means)
                         (mean_SLpc == sorted_means[2])}) %>%
                group_by(location, disease, strategy) %>%
                mutate(min_any   = any(min_mean, min_mean2),
                       alpha     = ifelse(min_any, 1, 0.25))
output_grid  <- output_data %>%
                group_by(location, disease) %>%
                summarise(intercept = seq(0, unique(bound_x + bound_y), length.out = 20),
                          bound     = max(SLpc[min_strat == TRUE])) 

gg <- ggplot(output_stats, aes(x = x, y = y, fill = strategy, alpha = alpha)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_abline(data = output_grid, aes(slope = -1, intercept = intercept), linewidth = 0.1, color = "grey") +  
      geom_point(data = output_data %>% filter(min_strat == FALSE),
                 aes(x = x, y = y, color = strategy, fill = strategy, alpha = alpha),
                 shape = 19, size = 0.25, stroke = 0.25) +
      geom_point(data = output_data %>% filter(min_strat == TRUE),
                 aes(x = x, y = y, color = strategy, fill = strategy, alpha = alpha),
                 shape = 19, size = 0.25, stroke = 0.25) +
      # geom_point(data = arch_data %>% filter(min_strat == TRUE & strategy == "No Closures"),
      #            aes(x = x, y = y, color = strategy, fill = strategy, alpha = alpha),
      #            shape = 19, size = 0.2, stroke = 0.25) +
      geom_abline(data = output_grid, aes(slope = -1, intercept = bound), linewidth = 0.25, color = "black") +
      geom_linerange(aes(x = mean_x, y = mean_y, xmin = q1_x, xmax = q3_x), linewidth = 0.25, color = "black") +
      geom_linerange(aes(x = mean_x, y = mean_y, ymin = q1_y, ymax = q3_y), linewidth = 0.25, color = "black") +
      geom_point(aes(x = mean_x, y = mean_y),
                 shape = 21, size = 2.5, stroke = 0.25, color = "black") +
      scale_color_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                    "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_fill_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_alpha_continuous(range = c(0.05, 1)) +
      theme_bw() +
      scale_x_continuous(position = "bottom", expand = c(0.01, 0)) +
      scale_y_continuous(position = "right", expand = c(0.02, 0)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(1.00, "lines")) +
      labs(title = "", x = "VLYL (% of GDP)", y = "GDPL + VSYL (% of GDP)") +
      guides(color = "none", fill = guide_legend(title = NULL), alpha = "none") + 
      theme(legend.position = "top", legend.box.just = "right", legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

ggsave("figure_3n.png", plot = gg, height = 14, width = 10)