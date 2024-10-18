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
             calc_cost_pc() %>% find_best_strats()
arch_tab  <- arch_data %>% group_by(location, disease, strategy) %>% 
             summarise(mean       = sprintf("%.1f", mean(SECpc)),
                       sd         = sprintf("%.1f", sd(SECpc)),
                       q1         = sprintf("%.1f", quantile(SECpc, 0.25)),
                       q2         = sprintf("%.1f", quantile(SECpc, 0.50)),
                       q3         = sprintf("%.1f", quantile(SECpc, 0.75)),
                       min_mean   = unique(min_mean),
                       min_med    = unique(min_med),
                       min_q3     = unique(min_q3)) %>% 
             group_by(location, disease) %>% 
             mutate(mean = if_else(min_mean, paste0("\\bfseries{",mean,"}"),             mean), 
                    q2   = if_else(min_med,  paste0(q2,"$^*$"),                          paste0(q2,"\\phantom{1}")),
                    q3   = if_else(min_q3,   paste0(q3,"\\textsuperscript\\textdagger"), paste0(q3,"\\phantom{1}"))) %>%
             ungroup() %>% 
             select(-starts_with("min")) %>%
             mutate(across(everything(), as.character)) %>%            
             insert_blank_rows("location")

write.table(arch_tab, file = "table_1.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)