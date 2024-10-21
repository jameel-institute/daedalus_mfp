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
             select(location, country, starts_with("obj"), disease, strategy, VLYL, GDPL, VSYL) %>%
             crossing(vsl = seq(0, 200, by = 2)) %>% mutate(VLYL = VLYL*(vsl/160)) %>% mutate(SEC = VLYL + GDPL + VSYL) %>% 
             calc_cost_pc()
arch_prop <- arch_data %>% group_by(location, disease, vsl, country) %>%
             summarise(min_mean = if_else(sum(SECpc == min(SECpc)) > 1, 
                                  paste(strategy[SECpc == min(SECpc)], collapse = ", "), 
                                  strategy[which.min(SECpc)]), .groups = 'drop') %>%
             group_by(location, disease, vsl) %>%
             count(min_mean) %>% mutate(proportion = n / sum(n)) %>% ungroup() %>%
             mutate(min_mean = case_when(min_mean %in% c("No Closures", "School Closures", "Economic Closures", "Elimination") ~ min_mean,
                                         min_mean == "School Closures, Economic Closures" ~ "Untriggered Closures",
                                         TRUE ~ "Other")) %>%
             mutate(min_mean = factor(min_mean, levels = c("No Closures", "Untriggered Closures", "School Closures", 
                                                           "Economic Closures", "Elimination", "Other")))

gg <- ggplot(data = arch_prop , aes(x = vsl, y = proportion, fill = min_mean, pattern_density = min_mean)) +
      facet_grid2(disease ~ location, switch = "y", scales = "fixed") +
      geom_area_pattern(pattern = "stripe", pattern_color = "darkgreen", pattern_fill = "darkgreen") +
      scale_fill_manual(values = c("No Closures" = "magenta4", "Untriggered Closures" = "navy", "School Closures" = "navy",
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod",  "Other" = "grey"))+  
      scale_pattern_density_manual(values = c("No Closures" = 0, "Untriggered Closures" = 0.2, "School Closures" = 0,
                                              "Economic Closures" = 0, "Elimination" = 0, "Other" = 0),
                                   breaks = c("No Closures", "Untriggered Closures", "School Closures", "Economic Closures", "Elimination")) +
      theme_bw() +
      scale_x_continuous(expand = c(0,0), position = "bottom") +
      scale_y_continuous(expand = c(0,0), position = "right", labels = scales::label_number(scale = 100, suffix = "")) +
      theme(panel.spacing = unit(1.00, "lines")) +
      labs(title = "", x = "VSL / GDP per capita", y = "Loss-Minimising Probability (%)") +
      guides(fill = "none", pattern_density = guide_legend(title = NULL, override.aes = list(fill = c("magenta4","navy","navy","darkgreen","goldenrod")))) +
      theme(legend.position = "top", legend.box.just = "right", legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

ggsave("figure_5.png", plot = gg, height = 14, width = 10)

# #alternative WTP plot
#
# test2 <- arch_data %>%
#   pivot_longer(
#     cols = starts_with("SECpc_"),
#     names_to = "strategy",
#     values_to = "SECpc"
#   ) %>%
#   mutate(strategy = str_remove(strategy, "SECpc_")) %>%
#   group_by(location, disease, vsl, strategy) %>%
#   summarise(mean_SECpc = mean(SECpc, na.rm = TRUE), .groups = "drop") %>%
#   mutate(strategy = factor(strategy, levels = c("No Closures", "School Closures", "Economic Closures", "Elimination")))
# 
# 
# gg <- ggplot(data = test2 , aes(x = vsl, y = mean_SECpc, color = strategy)) +
#   facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
#   #geom_ribbon(data = arch_fit %>% filter(alpha == 1), 
#   #            aes(ymin = lower, ymax = upper, fill = strategy, group = interaction(group, strategy)), color = NA, alpha = 0.2) +
#   #geom_point(shape = 21, size = 0.1, stroke = 0.35, fill = NA) +
#   geom_line(linewidth = 0.75) +
#   scale_color_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
#                                 "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
#   theme_bw() +
#   scale_x_continuous(expand = c(0, 0), position = "bottom") +
#   scale_y_log10(limits = c(1, NA), expand = c(0, 0), position = "right") +
#   #facetted_pos_scales(x = list(disease == "Covid-Omicron-X" & location == "HIC" ~ scale_x_continuous(labels = "Hospital Capacity"))) +   
#   theme(panel.spacing = unit(1.25, "lines")) +
#   labs(title = "", x = "Value of Statistical Life (% of GDP per capita)", y = "Societal Loss (% of GDP)") +
#   guides(color = guide_legend(title = NULL)) + 
#   theme(legend.position = c(0.1395, 0.999), legend.justification = c(1, 1), legend.box.just = "right", 
#         legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 7))
# 
# ggsave("figure_S5.png", plot = gg, height = 14, width = 10)