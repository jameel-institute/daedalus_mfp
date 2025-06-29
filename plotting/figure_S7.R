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
#source("functions/voi_fit.R")
source("functions/format_table.R")

list_files   <- list.files(path = "../output/archetypes/", pattern = "\\.csv$", full.names = TRUE)
input_files  <- list_files[grepl("_data\\.csv$", list_files)]
input_data   <- lapply(input_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols()
output_files <- list_files[!grepl("_data\\.csv$", list_files)]
output_data  <- lapply(output_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>%
                (function(x) calc_loss_pc(input_data,x)) %>% 
                group_by(location, disease, strategy) %>%
                mutate(med_SLpc   = quantile(VSYLpc, 0.50),
                       mean_SLpc  = mean(VSYLpc),
                       q3_SLpc    = quantile(VSYLpc, 0.75),
                       max_VSYLpc = max(VSYLpc)) %>%
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
                ungroup()

gg <- ggplot(output_data, aes(x = strategy, y = VSYLpc, linewidth = min_any, alpha = min_any)) + 
      facet_grid2(disease ~ location, switch = "y", scales = "free_y") +
      geom_violin(fill = "blue") + 
      scale_linewidth_manual(values = c("FALSE" = 0.1, "TRUE" = 0.5)) +
      scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 1)) +
      geom_text(data = output_data %>% filter(min_med == TRUE), aes(x = strategy, y = max_VSYLpc),
                vjust = 0.4, label = "*", size = 6, color = "black", inherit.aes = FALSE) +
      geom_text(data = output_data %>% filter(min_q3  == TRUE), aes(x = strategy, y = max_VSYLpc),
                vjust = -1.6, label = "†", size = 3.5, color = "black", inherit.aes = FALSE) +
      theme_bw() + 
      scale_x_discrete(labels = c("School Closures" = "Reactive/Sustained-School Closures",
                                  "Economic Closures" = "Reactive/Reactive-School Closures")) +
      facetted_pos_scales(y = list(
        scale_y_continuous(limits=c(0,160),  breaks=seq(0,160,40),  expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,160),  breaks=seq(0,160,40),  expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,200), breaks=seq(0,200,50), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,160), breaks=seq(0,160,40), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,160), breaks=seq(0,160,40), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,160), breaks=seq(0,160,40), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,200), breaks=seq(0,200,50), expand=c(0,0), position="right"))) +
      theme(panel.spacing = unit(0.75, "lines"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "", y = "VSYL (% of GDP)") +
      guides(linewidth = "none", alpha = "none")

ggsave("figure_S7.png", plot = gg, height = 14, width = 10)

output_table <- output_data %>% 
                group_by(location, disease, strategy) %>% 
                summarise(mean    = sprintf("%.1f", unique(mean_SLpc)),
                          sd      = sprintf("%.1f", sd(VSYLpc)),
                          q1      = sprintf("%.1f", quantile(VSYLpc, 0.25)),
                          q2      = sprintf("%.1f", unique(med_SLpc)),
                          q3      = sprintf("%.1f", unique(q3_SLpc)),
                          min_any = unique(min_any),
                          min_med = unique(min_med),
                          min_q3  = unique(min_q3)) %>% 
                mutate(strategy = case_when(strategy == "School Closures" ~ "Reactive/Sustained-School Closures",
                                            strategy == "Economic Closures" ~ "Reactive/Reactive-School Closures",
                                            TRUE ~ strategy)) %>%
                mutate(strategy = if_else(min_any, paste0("\\bfseries{",strategy,"}"), strategy)) %>%
                mutate(strategy = if_else(min_med, paste0(strategy,"$^*$"), strategy)) %>%       
                mutate(strategy = if_else(min_q3,  paste0(strategy,"\\textsuperscript\\textdagger"), strategy)) %>%    
                mutate(q3       = if_else(str_detect(strategy, "Elimination"),  paste0(q3,"\\phantom{.}"), q3)) %>%
                mutate(q3       = if_else(str_detect(strategy, "Elimination") & disease == "SARS-X",  paste0(q3,"\\phantom{.}"), q3)) %>%
                dplyr::select(-starts_with("min")) %>%
                mutate(across(everything(), as.character)) %>%            
                format_table("location")

write.table(output_table, file = "table_S7.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)