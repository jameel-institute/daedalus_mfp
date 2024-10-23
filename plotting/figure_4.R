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
             calc_cost_pc() %>% 
             mutate(x = VLYLpc, y = GDPLpc + VSYLpc) %>%  
             group_by(location, disease, country) %>% mutate(min_strat = (SECpc == min(SECpc))) %>% 
             group_by(location, disease, strategy) %>% 
               mutate(alpha = 1 - abs(cor(x, y))) %>% 
               mutate(alpha = case_when(
                        disease == "Influenza-2009-X" & strategy == "Elimination" & (GDPLpc+VSYLpc) < 50 ~ alpha/6,
                        disease == "Influenza-1957-X" & strategy != "No Closures" ~ alpha/3,
                        disease == "Influenza-1918-X" | 
                        disease == "Covid-Delta-X" | 
                        disease == "Covid-Wildtype-X" & strategy == "No Closures" ~ alpha/3,
                        location == "HIC" & disease == "SARS-X" & strategy == "No Closures" ~ alpha*20,
                        TRUE ~ alpha)) %>%
               mutate(alpha = 0.03 + alpha/10) %>%
               mutate(alpha = ifelse(min_strat, 1, alpha)) %>%
             ungroup()
arch_stat <- arch_data %>% 
             group_by(location, disease, strategy) %>%
             summarize(mean_x  = mean(x),
                       lower_x = quantile(x, 0.25),
                       upper_x = quantile(x, 0.75),
                       bound_x = quantile(x, 0.95),
                       mean_y  = mean(y),
                       lower_y = quantile(y, 0.25),
                       upper_y = quantile(y, 0.75),
                       bound_y = quantile(y, 0.95)) %>%
             group_by(location, disease) %>%
             mutate(bound_x = max(bound_x), 
                    bound_y = max(bound_y))
arch_data <- arch_data %>%
             left_join(arch_stat %>% select(location, disease, strategy, bound_x, bound_y), by = c("location", "disease", "strategy")) %>%
             group_by(location, disease) %>%
             filter(x <= bound_x & y <= bound_y) %>%
             ungroup()
arch_grid <- arch_data %>%
             group_by(location, disease) %>%
             summarise(intercept = seq(0, max(bound_x + bound_y), length.out = 20),
                       bound     = max(SECpc[min_strat == TRUE]),
                       .groups   = 'drop') 

gg <- ggplot(data = arch_stat, aes(x = mean_x, y = mean_y, color = strategy, fill = strategy)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_abline(data = arch_grid, aes(slope = -1, intercept = intercept), size = 0.1, color = "grey") +  
      geom_point(data = arch_data %>% filter(min_strat == FALSE),
                 aes(x = x, y = y, color = strategy, fill = strategy, alpha = alpha),
                 shape = 19, size = 0.2, stroke = 0.25) +
      geom_point(data = arch_data %>% filter(min_strat == TRUE),
                 aes(x = x, y = y, color = strategy, fill = strategy, alpha = alpha),
                 shape = 19, size = 0.2, stroke = 0.25) +
      geom_point(data = arch_data %>% filter(min_strat == TRUE & strategy == "No Closures"),
                 aes(x = x, y = y, color = strategy, fill = strategy, alpha = alpha),
                 shape = 19, size = 0.2, stroke = 0.25) +
      geom_abline(data = arch_grid, aes(slope = -1, intercept = bound), size = 0.2, color = "black") +
      geom_pointrange(aes(x = mean_x, y = mean_y, xmin = lower_x, xmax = upper_x),
                      size = 0.0, linewidth = 0.4, color = "black") +
      geom_pointrange(aes(x = mean_x, y = mean_y, ymin = lower_y, ymax = upper_y),
                      size = 0.0, linewidth = 0.4, color = "black") +
      geom_point(aes(x = mean_x, y = mean_y),
                 shape = 21, size = 2, stroke = 0.25, color = "black") +
      scale_color_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                    "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_fill_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_alpha_continuous(range = c(0.03, 1)) +
      theme_bw() +
      scale_x_continuous(position = "bottom", expand = c(0.01, 0)) +
      scale_y_continuous(position = "right", expand = c(0.02, 0)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(1.00, "lines")) +
      labs(title = "", x = "VLYL (% of GDP)", y = "GDPL + VSYL (% of GDP)") +
      guides(color = "none", fill = guide_legend(title = NULL), alpha = "none") + 
      theme(legend.position = "top", legend.box.just = "right", legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

ggsave("figure_4.png", plot = gg, height = 14, width = 10)