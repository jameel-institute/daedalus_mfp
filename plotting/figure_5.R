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
                pivot_longer(cols = starts_with("gdpl_"), names_to = "gdpl", values_to = "value") %>%
                mutate(gdpl = factor(gdpl, levels = c("gdpl_1", "gdpl_2", "gdpl_3", "gdpl_4", "gdpl_5", 
                                                      "gdpl_6", "gdpl_7", "gdpl_8", "gdpl_9", "gdpl_10")))

gg <- ggplot(output_data, aes(x = strategy, y = value, fill = gdpl, alpha = min_any)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free_y") +
      geom_boxplot(outlier.shape = NA, coef = 1.5, width = 0.70, linewidth = 0.2) +
      scale_fill_manual(values = c("gdpl_1" = "green", "gdpl_2" = "grey50", "gdpl_3" = "red", "gdpl_4" = "cyan", "gdpl_5" = "orange",
                                   "gdpl_6" = "purple", "gdpl_7" = "black", "gdpl_8" = "white", "gdpl_9" = "blue", "gdpl_10" = "yellow"),
                        labels = c("Agriculture, Forestry, Fishing", "Mining", "Manufacturing", "Utilities", "Construction",
                                   "Retail, Hospitality", "Transport", "IT, Telecommunications", "Finance, Professional, Technical", "Public Administration, Other Services")) + 
      scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 1)) +
      theme_bw() + 
      facetted_pos_scales(y = list(
        scale_y_continuous(limits=c(0,12), breaks=seq(0,12,3), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,12), breaks=seq(0,12,3), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,24), breaks=seq(0,24,6), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,12), breaks=seq(0,12,3), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0), position="right"))) +
      theme(panel.spacing = unit(0.75, "lines"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "", y = "GDP Loss (% of GDP)") +
      guides(fill = guide_legend(title = NULL, byrow = TRUE), alpha = "none") +
      theme(legend.position = "bottom", legend.box.just = "right", 
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9), 
            legend.margin = margin(0, 0, 0, 0))

ggsave("figure_5.png", plot = gg, height = 14, width = 10)