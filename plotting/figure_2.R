library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(fBasics)
library(fitdistrplus)
library(forecast)
library(scam)
library(ggplot2)
library(ggh4x)
library(ggdensity)
library(GGally)
library(cowplot)
library(ggpattern)
library(patchwork)
source("functions/add_scenario_cols.R")
source("functions/order_scenario_cols.R")
source("functions/calc_loss_pc.R")
source("functions/parse_inputs.R")
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
                mutate(pLYL  = mean(VLYLpc/SLpc),#mean of proportions
                       pGDPL = mean(GDPLpc/SLpc),
                       pSYL  = mean(VSYLpc/SLpc)) %>% 
                ungroup()

gg <- ggplot(output_data, aes(x = strategy, y = SLpc, fill = factor(..fill..), alpha = min_any), color = "black") + 
      facet_grid2(disease ~ location, switch = "y", scales = "free_y") +
      geom_violin(aes(width = pLYL+pGDPL+pSYL, linewidth = min_any, fill = "red")) + 
      geom_violin(aes(width = pGDPL+pSYL, fill = "white"), linewidth = 0, alpha = 1) +
      geom_violin(aes(width = pGDPL+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0, alpha = 1) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_linewidth_manual(values = c("FALSE" = 0.1, "TRUE" = 0.5)) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"), 
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 1)) +
      geom_text(data = output_data %>% filter(min_med == TRUE), aes(x = strategy, y = max_SLpc),
                vjust = 0.4, label = "*", size = 6, color = "black", inherit.aes = FALSE) +
      geom_text(data = output_data %>% filter(min_q3  == TRUE), aes(x = strategy, y = max_SLpc),
                vjust = -1.5, label = "†", size = 3.5, color = "black", inherit.aes = FALSE) +
      theme_bw() + 
      scale_x_discrete(labels = c("School Closures" = "Reactive-Business/Sustained-School Closures",
                                  "Economic Closures" = "Reactive-Business/Reactive-School Closures")) +
      facetted_pos_scales(y = list(
        scale_y_continuous(limits=c(0,200),  breaks=seq(0,200,50),    expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,320),  breaks=seq(0,320,80),    expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,2000), breaks=seq(0,2000,500),  expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,400),  breaks=seq(0,400,100),   expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,1600), breaks=seq(0,1600,400),  expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,800),  breaks=seq(0,800,200),   expand=c(0,0), position="right"))) +
        #scale_y_continuous(limits=c(0,2800), breaks=seq(0,2800,700),  expand=c(0,0), position="right"))) +
      theme(panel.spacing = unit(0.75, "lines"), axis.text.x = element_text(angle = 55, hjust = 1)) + 
      labs(title = "", x = "", y = "Socioeconomic Loss (% of GDP)") +
      guides(width = "none", linewidth = "none", color = "none", fill = guide_legend(title = NULL), alpha = "none") +
      theme(legend.position = "bottom", legend.box.just = "right", 
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9), 
            legend.margin = margin(0, 0, 0, 0))

ggsave("figure_2.png", plot = gg, height = 14, width = 10)

output_table <- output_data %>%
                group_by(location, disease, strategy) %>%
                summarise(mean    = unique(mean_SLpc),
                          sd      = sd(SLpc),
                          q1      = quantile(SLpc, 0.25),
                          q2      = unique(med_SLpc),
                          q3      = unique(q3_SLpc),
                          min_any = unique(min_any),
                          min_med = unique(min_med),
                          min_q3  = unique(min_q3)) %>%
                mutate(across(c(mean, sd, q1, q2, q3), ~ ifelse(.x < 10, sprintf("%.1f", .x), sprintf("%.0f", .x)))) %>%
                mutate(strategy = case_when(strategy == "School Closures" ~ "Reactive-Business/Sustained-School Closures",
                                            strategy == "Economic Closures" ~ "Reactive-Business/Reactive-School Closures",
                                            TRUE ~ strategy)) %>%
                mutate(strategy = if_else(min_any, paste0("\\bfseries{",strategy,"}"), strategy)) %>%
                mutate(strategy = if_else(min_med, paste0(strategy,"$^*$"), strategy)) %>%
                mutate(strategy = if_else(min_q3,  paste0(strategy,"\\textsuperscript\\textdagger"), strategy)) %>%
                mutate(q3       = if_else(str_detect(strategy, "Elimination"),  paste0(q3,"\\phantom{.}"), q3)) %>%
                mutate(q3       = if_else(str_detect(strategy, "Elimination") & disease == "Covid-Wildtype-X",  paste0(q3,"\\phantom{.}"), q3)) %>%
                dplyr::select(-starts_with("min")) %>%
                mutate(across(everything(), as.character)) %>%
                format_table("location")

write.table(output_table, file = "table_S5.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)