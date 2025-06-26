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
                mutate(med_SLpc  = quantile(SLpc, 0.50),
                       mean_SLpc = mean(SLpc),
                       q3_SLpc   = quantile(SLpc, 0.75)) %>%
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
      scale_x_discrete(labels = c("School Closures" = "Reactive/Sustained-School Closures",
                                  "Economic Closures" = "Reactive/Reactive-School Closures")) +
      facetted_pos_scales(y = list(
        scale_y_continuous(limits=c(0,12), breaks=seq(0,12,3), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,16), breaks=seq(0,16,4), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,28), breaks=seq(0,28,7), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,12), breaks=seq(0,12,3), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0), position="right"))) +
      theme(panel.spacing = unit(0.75, "lines"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "", y = "GDPL (% of GDP)") +
      guides(fill = guide_legend(title = NULL, byrow = TRUE), alpha = "none") +
      theme(legend.position = "bottom", legend.box.just = "right", 
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9), 
            legend.margin = margin(0, 0, 0, 0))

ggsave("figure_S12.png", plot = gg, height = 14, width = 10)

output_table <- output_data %>%
                pivot_wider(names_from = gdpl, values_from = value) %>%
                group_by(location, disease, strategy) %>% 
                summarise(median_1  = sprintf("%.1f", quantile(gdpl_1,  0.50)),
                          q1_1      = sprintf("%.1f", quantile(gdpl_1,  0.25)),
                          q3_1      = sprintf("%.1f", quantile(gdpl_1,  0.75)),
                          median_2  = sprintf("%.1f", quantile(gdpl_2,  0.50)),
                          q1_2      = sprintf("%.1f", quantile(gdpl_2,  0.25)),
                          q3_2      = sprintf("%.1f", quantile(gdpl_2,  0.75)),
                          median_3  = sprintf("%.1f", quantile(gdpl_3,  0.50)),
                          q1_3      = sprintf("%.1f", quantile(gdpl_3,  0.25)),
                          q3_3      = sprintf("%.1f", quantile(gdpl_3,  0.75)),
                          median_4  = sprintf("%.1f", quantile(gdpl_4,  0.50)),
                          q1_4      = sprintf("%.1f", quantile(gdpl_4,  0.25)),
                          q3_4      = sprintf("%.1f", quantile(gdpl_4,  0.75)),
                          median_5  = sprintf("%.1f", quantile(gdpl_5,  0.50)),
                          q1_5      = sprintf("%.1f", quantile(gdpl_5,  0.25)),
                          q3_5      = sprintf("%.1f", quantile(gdpl_5,  0.75)),
                          median_6  = sprintf("%.1f", quantile(gdpl_6,  0.50)),
                          q1_6      = sprintf("%.1f", quantile(gdpl_6,  0.25)),
                          q3_6      = sprintf("%.1f", quantile(gdpl_6,  0.75)),
                          median_7  = sprintf("%.1f", quantile(gdpl_7,  0.50)),
                          q1_7      = sprintf("%.1f", quantile(gdpl_7,  0.25)),
                          q3_7      = sprintf("%.1f", quantile(gdpl_7,  0.75)),
                          median_8  = sprintf("%.1f", quantile(gdpl_8,  0.50)),
                          q1_8      = sprintf("%.1f", quantile(gdpl_8,  0.25)),
                          q3_8      = sprintf("%.1f", quantile(gdpl_8,  0.75)),
                          median_9  = sprintf("%.1f", quantile(gdpl_9,  0.50)),
                          q1_9      = sprintf("%.1f", quantile(gdpl_9,  0.25)),
                          q3_9      = sprintf("%.1f", quantile(gdpl_9,  0.75)),
                          median_10 = sprintf("%.1f", quantile(gdpl_10, 0.50)),
                          q1_10     = sprintf("%.1f", quantile(gdpl_10, 0.25)),
                          q3_10     = sprintf("%.1f", quantile(gdpl_10, 0.75)),
                          min_any   = unique(min_any),
                          min_med   = unique(min_med),
                          min_q3    = unique(min_q3)) %>%
                mutate(median_1  = paste0(median_1, "\\textcolor{white}{aaa} (",q1_1, "; ",q3_1, ")"),
                       median_2  = paste0(median_2, "\\textcolor{white}{aaa} (",q1_2, "; ",q3_2, ")"),
                       median_3  = paste0(median_3, "\\textcolor{white}{aaa} (",q1_3, "; ",q3_3, ")"),
                       median_4  = paste0(median_4, "\\textcolor{white}{aaa} (",q1_4, "; ",q3_4, ")"),
                       median_5  = paste0(median_5, "\\textcolor{white}{aaa} (",q1_5, "; ",q3_5, ")"),
                       median_6  = paste0(median_6, "\\textcolor{white}{aaa} (",q1_6, "; ",q3_6, ")"),
                       median_7  = paste0(median_7, "\\textcolor{white}{aaa} (",q1_7, "; ",q3_7, ")"),
                       median_8  = paste0(median_8, "\\textcolor{white}{aaa} (",q1_8, "; ",q3_8, ")"),
                       median_9  = paste0(median_9, "\\textcolor{white}{aaa} (",q1_9, "; ",q3_9, ")"),
                       median_10 = paste0(median_10,"\\textcolor{white}{aaa} (",q1_10,"; ",q3_10,")")) %>%
                # mutate(median_1  = if_else(min_any, paste0("\\bfseries{",median_1 ,"}"), median_1),
                #        median_2  = if_else(min_any, paste0("\\bfseries{",median_2 ,"}"), median_2),
                #        median_3  = if_else(min_any, paste0("\\bfseries{",median_3 ,"}"), median_3),
                #        median_4  = if_else(min_any, paste0("\\bfseries{",median_4 ,"}"), median_4),
                #        median_5  = if_else(min_any, paste0("\\bfseries{",median_5 ,"}"), median_5),
                #        median_6  = if_else(min_any, paste0("\\bfseries{",median_6 ,"}"), median_6),
                #        median_7  = if_else(min_any, paste0("\\bfseries{",median_7 ,"}"), median_7),
                #        median_8  = if_else(min_any, paste0("\\bfseries{",median_8 ,"}"), median_8),
                #        median_9  = if_else(min_any, paste0("\\bfseries{",median_9 ,"}"), median_9),
                #        median_10 = if_else(min_any, paste0("\\bfseries{",median_10,"}"), median_10)) %>%
                mutate(strategy = case_when(strategy == "School Closures" ~ "Reactive/Sustained-School Closures",
                                            strategy == "Economic Closures" ~ "Reactive/Reactive-School Closures",
                                            TRUE ~ strategy)) %>%
                mutate(strategy  = if_else(min_any, paste0("\\bfseries{",strategy,"}"), strategy)) %>%
                mutate(strategy  = if_else(min_med, paste0(strategy,"$^*$"), strategy)) %>%       
                mutate(strategy  = if_else(min_q3,  paste0(strategy,"\\textsuperscript\\textdagger"), strategy)) %>%    
                mutate(median_10 = if_else(str_detect(strategy, "Elimination"),  paste0(median_10,"\\phantom{.}"), median_10)) %>%
                mutate(median_10 = if_else(str_detect(strategy, "Elimination") & disease == "SARS-X",  paste0(median_10,"\\phantom{.}"), median_10)) %>%
                mutate(median_10 = if_else(str_detect(strategy, "Elimination") & disease == "Influenza-1918-X", str_remove_all(median_10, "\\\\phantom\\{\\.\\}"), median_10)) %>%
                dplyr::select(-starts_with("q"),-starts_with("min")) %>%
                mutate(across(everything(), as.character)) %>%            
                format_table("location")

write.table(output_table, file = "table_S10.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)