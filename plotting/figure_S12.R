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
                (function(x) calc_loss_pc(input_data,x)) %>% 
                group_by(location, disease, strategy) %>%
                mutate(med_SLpc  = quantile(SLpc, 0.50),
                       mean_SLpc = mean(SLpc),
                       q3_SLpc   = quantile(SLpc, 0.75),
                       max_SLpc  = max(SLpc)) %>%
                group_by(location, disease) %>%
                mutate(min_med   = (med_SLpc  == min(med_SLpc)),
                       min_mean  = (mean_SLpc == min(mean_SLpc)),
                       min_mean2 = {
                         strategy_means <- mean_SLpc[!duplicated(strategy)]
                         sorted_means   <- sort(strategy_means)
                         (mean_SLpc == sorted_means[2])},
                       min_q3    = (q3_SLpc   == min(q3_SLpc))) %>%
                group_by(location, disease, strategy) %>%
                mutate(min_any   = any(min_mean, min_mean2)) %>%
                ungroup() %>%
                pivot_longer(cols = starts_with("vlyl_"), names_to = "vlyl", values_to = "value") %>%
                mutate(vlyl = factor(vlyl, levels = c("vlyl_1", "vlyl_2", "vlyl_3", "vlyl_4")))

gg <- ggplot(output_data, aes(x = strategy, y = value, fill = vlyl, alpha = min_any)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free_y") +
      geom_boxplot(outlier.shape = NA, coef = 1.5, width = 0.70, linewidth = 0.2) +
      scale_fill_manual(values = c("vlyl_1" = "magenta", "vlyl_2" = "green", "vlyl_3" = "brown", "vlyl_4" = "grey80"),
                        labels = c("Preschool-Age", "School-Age", "Working-Age", "Retired-Age")) +
      scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 1)) +
      theme_bw() + 
      facetted_pos_scales(y = list(
        scale_y_continuous(limits=c(0,4),    breaks=seq(0,4,1),      expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,16),   breaks=seq(0,16,4),     expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,1000), breaks=seq(0,1000,250), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,100),  breaks=seq(0,100,25),   expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,800),  breaks=seq(0,800,200),  expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,320),  breaks=seq(0,320,80),   expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,1000), breaks=seq(0,1000,250), expand=c(0,0), position="right"))) +
      theme(panel.spacing = unit(0.75, "lines"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "", y = "VLYL (% of GDP)") +
      guides(fill = guide_legend(title = NULL), alpha = "none") +
      theme(legend.position = "bottom", legend.box.just = "right", 
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9), 
            legend.margin = margin(0, 0, 0, 0))

ggsave("figure_S12.png", plot = gg, height = 14, width = 10)