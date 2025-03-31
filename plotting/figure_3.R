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
                mutate(x = VLYLpc, y = GDPLpc + VSYLpc) %>%
                group_by(location, country, disease) %>%
                mutate(x = x[strategy == "No Closures"] - x,
                       y = y - y[strategy == "No Closures"]) %>%
                ungroup() %>%
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
                mutate(alpha = ifelse(min_any, 1, 0.05))
output_stats <- output_data %>% #for quicker plotting
                group_by(location, disease, strategy) %>% 
                summarise(mean_x  = mean(x), 
                          q1_x    = quantile(x, 0.25), 
                          q3_x    = quantile(x, 0.75),
                          mean_y  = mean(y), 
                          q1_y    = quantile(y, 0.25), 
                          q3_y    = quantile(y, 0.75),
                          min_any = unique(min_any),
                          min_med = unique(min_med),
                          min_q3  = unique(min_q3)) %>%
                mutate(alpha = ifelse(min_any, 1, 0.25))
                
gg <- ggplot(output_data, aes(x = x, y = y, color = strategy, fill = strategy, alpha = alpha)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_hline(yintercept = 0, linewidth = 0.10, color = "black") +  
      geom_vline(xintercept = 0, linewidth = 0.10, color = "black") +  
      geom_point(shape = 19, size = 0.02, stroke = 0.2) +
      geom_abline(slope = 1, linetype = "dashed", linewidth = 0.25, color = "black") +  
      geom_linerange(data = output_stats, aes(x = mean_x, y = mean_y, xmin = q1_x, xmax = q3_x), linewidth = 0.25, color = "black") +
      geom_linerange(data = output_stats, aes(x = mean_x, y = mean_y, ymin = q1_y, ymax = q3_y), linewidth = 0.25, color = "black") +
      geom_point(data = output_stats, aes(x = mean_x, y = mean_y), shape = 21, size = 2, stroke = 0.2, color = "black") +
      scale_color_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                    "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_fill_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      scale_alpha_continuous(range = c(0.05, 1)) +
      theme_bw() +
      facetted_pos_scales(
      x = list(
      location == "LLMIC" & disease == "Influenza-2009-X" ~ scale_x_continuous(limits=c(-3,9),breaks=seq(-3,9,3),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Influenza-2009-X" ~ scale_x_continuous(limits=c(-2,6),breaks=seq(-2,6,2),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Influenza-2009-X" ~ scale_x_continuous(limits=c(-1,3),breaks=seq(-1,3,1),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Influenza-1957-X" ~ scale_x_continuous(limits=c(-20,60),breaks=seq(-20,60,20),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Influenza-1957-X" ~ scale_x_continuous(limits=c(-20,60),breaks=seq(-20,60,20),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Influenza-1957-X" ~ scale_x_continuous(limits=c(-20,60),breaks=seq(-20,60,20),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Influenza-1918-X" ~ scale_x_continuous(limits=c(-400,1200),breaks=seq(-400,1200,400),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Influenza-1918-X" ~ scale_x_continuous(limits=c(-400,1200),breaks=seq(-400,1200,400),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Influenza-1918-X" ~ scale_x_continuous(limits=c(-400,1200),breaks=seq(-400,1200,400),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Covid-Omicron-X" ~ scale_x_continuous(limits=c(-30,90),breaks=seq(-30,90,30),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Covid-Omicron-X" ~ scale_x_continuous(limits=c(-30,90),breaks=seq(-30,90,30),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Covid-Omicron-X" ~ scale_x_continuous(limits=c(-30,90),breaks=seq(-30,90,30),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Covid-Delta-X" ~ scale_x_continuous(limits=c(-200,600),breaks=seq(-200,600,200),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Covid-Delta-X" ~ scale_x_continuous(limits=c(-200,600),breaks=seq(-200,600,200),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Covid-Delta-X" ~ scale_x_continuous(limits=c(-200,600),breaks=seq(-200,600,200),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Covid-Wildtype-X" ~ scale_x_continuous(limits=c(-100,300),breaks=seq(-100,300,100),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Covid-Wildtype-X" ~ scale_x_continuous(limits=c(-100,300),breaks=seq(-100,300,100),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Covid-Wildtype-X" ~ scale_x_continuous(limits=c(-100,300),breaks=seq(-100,300,100),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "SARS-X" ~ scale_x_continuous(limits=c(-500,1500),breaks=seq(-500,1500,500),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "SARS-X" ~ scale_x_continuous(limits=c(-500,1500),breaks=seq(-500,1500,500),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "SARS-X" ~ scale_x_continuous(limits=c(-500,1500),breaks=seq(-500,1500,500),expand=c(0,0),position="bottom")),
      y = list(
      location == "LLMIC" & disease == "Influenza-2009-X" ~ scale_y_continuous(limits=c(-50,150),breaks=seq(-50,150,50),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Influenza-2009-X" ~ scale_y_continuous(limits=c(-30,90),breaks=seq(-30,90,30),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Influenza-2009-X" ~ scale_y_continuous(limits=c(-20,60),breaks=seq(-20,60,20),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Influenza-1957-X" ~ scale_y_continuous(limits=c(-60,180),breaks=seq(-60,180,60),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Influenza-1957-X" ~ scale_y_continuous(limits=c(-30,90),breaks=seq(-30,90,30),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Influenza-1957-X" ~ scale_y_continuous(limits=c(-20,60),breaks=seq(-20,60,20),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Influenza-1918-X" ~ scale_y_continuous(limits=c(-80,240),breaks=seq(-80,240,80),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Influenza-1918-X" ~ scale_y_continuous(limits=c(-50,150),breaks=seq(-50,150,50),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Influenza-1918-X" ~ scale_y_continuous(limits=c(-50,150),breaks=seq(-50,150,50),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Covid-Omicron-X" ~ scale_y_continuous(limits=c(-60,180),breaks=seq(-60,180,60),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Covid-Omicron-X" ~ scale_y_continuous(limits=c(-40,120),breaks=seq(-40,120,40),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Covid-Omicron-X" ~ scale_y_continuous(limits=c(-25,75),breaks=seq(-25,75,25),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Covid-Delta-X" ~ scale_y_continuous(limits=c(-80,240),breaks=seq(-80,240,80),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Covid-Delta-X" ~ scale_y_continuous(limits=c(-60,180),breaks=seq(-60,180,60),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Covid-Delta-X" ~ scale_y_continuous(limits=c(-50,150),breaks=seq(-50,150,50),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Covid-Wildtype-X" ~ scale_y_continuous(limits=c(-80,240),breaks=seq(-80,240,80),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Covid-Wildtype-X" ~ scale_y_continuous(limits=c(-50,150),breaks=seq(-50,150,50),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Covid-Wildtype-X" ~ scale_y_continuous(limits=c(-30,90),breaks=seq(-30,90,30),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "SARS-X" ~ scale_y_continuous(limits=c(-200,200),breaks=seq(-200,200,100),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "SARS-X" ~ scale_y_continuous(limits=c(-100,100),breaks=seq(-100,100,50),expand=c(0,0),position="right"),
      location == "HIC" & disease == "SARS-X" ~ scale_y_continuous(limits=c(-100,100),breaks=seq(-100,100,50),expand=c(0,0),position="right"))) +
      theme(panel.spacing = unit(0.75, "lines")) +
      labs(title = "", x = "Decremental VLYL (% of GDP)", y = "Incremental GDPL + VSYL (% of GDP)") +
      guides(color = "none", 
             fill  = guide_legend(title = NULL),
             alpha = "none") + 
      theme(legend.position = "bottom", legend.box.just = "right", 
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9), 
            legend.margin = margin(0, 0, 0, 0))

ggsave("figure_3.png", plot = gg, height = 14, width = 10)

output_table <- output_stats %>% 
                mutate(mean_x = sprintf("%.1f", mean_x),
                       q1_x   = sprintf("%.1f", q1_x),
                       q3_x   = sprintf("%.1f", q3_x),
                       mean_y = sprintf("%.1f", mean_y),
                       q1_y   = sprintf("%.1f", q1_y),
                       q3_y   = sprintf("%.1f", q3_y)) %>%
                mutate(mean_x = paste0(mean_x," (",q1_x,"; ",q3_x,")"),
                       mean_y = paste0(mean_y," (",q1_y,"; ",q3_y,")")) %>%
                # mutate(mean_x = if_else(min_any, paste0("\\bfseries{",mean_x,"}"), mean_x),
                #        mean_y = if_else(min_any, paste0("\\bfseries{",mean_y,"}"), mean_y)) %>%
                mutate(strategy = if_else(min_any, paste0("\\bfseries{",strategy,"}"), strategy)) %>%
                mutate(strategy = if_else(min_med, paste0(strategy,"$^*$"), strategy)) %>%       
                mutate(strategy = if_else(min_q3,  paste0(strategy,"\\textsuperscript\\textdagger"), strategy)) %>%    
                mutate(mean_y   = if_else(str_detect(strategy, "Elimination"),  paste0(mean_y,"\\phantom{.}"), mean_y)) %>%
                mutate(mean_y   = if_else(str_detect(strategy, "Elimination") & disease == "SARS-X",  paste0(mean_y,"\\phantom{.}"), mean_y)) %>%
                dplyr::select(-starts_with("q"),-starts_with("min"),-alpha) %>%
                mutate(across(everything(), as.character)) %>%            
                format_table("location")

write.table(output_table, file = "table_S10.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)