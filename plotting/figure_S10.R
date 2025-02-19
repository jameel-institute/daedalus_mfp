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
                mutate(x = VLYLpc, y = GDPLpc + VSYLpc) %>%
                group_by(location, disease, strategy) %>%
                mutate(mean_SLpc = mean(SLpc)) %>%
                group_by(location, disease) %>%
                mutate(min_mean  = (mean_SLpc == min(mean_SLpc)),
                       min_mean2 = {
                         strategy_means <- mean_SLpc[!duplicated(strategy)]
                         sorted_means   <- sort(strategy_means)
                         (mean_SLpc == sorted_means[2])}) %>%
                group_by(location, disease, strategy) %>%
                mutate(min_any   = any(min_mean, min_mean2)) %>%
                ungroup()
output_stats <- output_data %>% #for quicker plotting
                group_by(location, disease, strategy) %>% 
                summarise(mean_x  = mean(x), 
                          q1_x    = quantile(x, 0.25), 
                          q3_x    = quantile(x, 0.75),
                          bound_x = quantile(x, 0.90),
                          mean_y  = mean(y), 
                          q1_y    = quantile(y, 0.25), 
                          q3_y    = quantile(y, 0.75),
                          bound_y = quantile(y, 0.90),
                          min_any = unique(min_any)) %>%
                mutate(alpha = ifelse(min_any, 1, 0.25))
output_grid  <- output_stats %>%
                group_by(location, disease) %>%
                summarise(intercept = seq(0, max(bound_x) + max(bound_y), length.out = 20)) 

gg <- ggplot(output_data, aes(x = x, y = y, fill = strategy)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_abline(data = output_grid, aes(slope = -1, intercept = intercept), linewidth = 0.10, color = "grey") +  
      geom_hdr(probs = c(0.001,0.20,0.50), color = NA) +
      geom_linerange(data = output_stats, aes(x = mean_x, y = mean_y, xmin = q1_x, xmax = q3_x), linewidth = 0.25, color = "black") +
      geom_linerange(data = output_stats, aes(x = mean_x, y = mean_y, ymin = q1_y, ymax = q3_y), linewidth = 0.25, color = "black") +
      geom_point(data = output_stats, aes(x = mean_x, y = mean_y), shape = 21, size = 2, stroke = 0.2, color = "black") +
      scale_fill_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod")) +
      theme_bw() +
      facetted_pos_scales(
      x = list(
      location == "LLMIC" & disease == "Influenza-2009-X" ~ scale_x_continuous(limits=c(0,4),breaks=seq(0,4,1),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Influenza-2009-X" ~ scale_x_continuous(limits=c(0,4),breaks=seq(0,4,1),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Influenza-2009-X" ~ scale_x_continuous(limits=c(0,4),breaks=seq(0,4,1),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Influenza-1957-X" ~ scale_x_continuous(limits=c(0,12),breaks=seq(0,12,3),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Influenza-1957-X" ~ scale_x_continuous(limits=c(0,12),breaks=seq(0,12,3),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Influenza-1957-X" ~ scale_x_continuous(limits=c(0,12),breaks=seq(0,12,3),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Influenza-1918-X" ~ scale_x_continuous(limits=c(0,1600),breaks=seq(0,1600,400),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Influenza-1918-X" ~ scale_x_continuous(limits=c(0,1200),breaks=seq(0,1200,300),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Influenza-1918-X" ~ scale_x_continuous(limits=c(0,400),breaks=seq(0,400,100),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Covid-Omicron-X" ~ scale_x_continuous(limits=c(0,60),breaks=seq(0,60,15),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Covid-Omicron-X" ~ scale_x_continuous(limits=c(0,60),breaks=seq(0,60,15),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Covid-Omicron-X" ~ scale_x_continuous(limits=c(0,60),breaks=seq(0,60,15),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Covid-Delta-X" ~ scale_x_continuous(limits=c(0,520),breaks=seq(0,520,130),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Covid-Delta-X" ~ scale_x_continuous(limits=c(0,720),breaks=seq(0,720,180),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Covid-Delta-X" ~ scale_x_continuous(limits=c(0,720),breaks=seq(0,720,180),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "Covid-Wildtype-X" ~ scale_x_continuous(limits=c(0,160),breaks=seq(0,160,40),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "Covid-Wildtype-X" ~ scale_x_continuous(limits=c(0,200),breaks=seq(0,200,50),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "Covid-Wildtype-X" ~ scale_x_continuous(limits=c(0,120),breaks=seq(0,120,30),expand=c(0,0),position="bottom"),
      location == "LLMIC" & disease == "SARS-X" ~ scale_x_continuous(limits=c(0,1200),breaks=seq(0,1200,300),expand=c(0,0),position="bottom"),
      location == "UMIC" & disease == "SARS-X" ~ scale_x_continuous(limits=c(0,1000),breaks=seq(0,1000,250),expand=c(0,0),position="bottom"),
      location == "HIC" & disease == "SARS-X" ~ scale_x_continuous(limits=c(0,400),breaks=seq(0,400,100),expand=c(0,0),position="bottom")),
      y = list(
      location == "LLMIC" & disease == "Influenza-2009-X" ~ scale_y_continuous(limits=c(0,160),breaks=seq(0,160,40),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Influenza-2009-X" ~ scale_y_continuous(limits=c(0,40),breaks=seq(0,40,10),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Influenza-2009-X" ~ scale_y_continuous(limits=c(0,20),breaks=seq(0,20,5),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Influenza-1957-X" ~ scale_y_continuous(limits=c(0,160),breaks=seq(0,160,40),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Influenza-1957-X" ~ scale_y_continuous(limits=c(0,60),breaks=seq(0,60,15),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Influenza-1957-X" ~ scale_y_continuous(limits=c(0,20),breaks=seq(0,20,5),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Influenza-1918-X" ~ scale_y_continuous(limits=c(0,320),breaks=seq(0,320,80),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Influenza-1918-X" ~ scale_y_continuous(limits=c(0,160),breaks=seq(0,160,40),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Influenza-1918-X" ~ scale_y_continuous(limits=c(0,60),breaks=seq(0,60,15),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Covid-Omicron-X" ~ scale_y_continuous(limits=c(0,240),breaks=seq(0,240,60),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Covid-Omicron-X" ~ scale_y_continuous(limits=c(0,120),breaks=seq(0,120,30),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Covid-Omicron-X" ~ scale_y_continuous(limits=c(0,60),breaks=seq(0,60,15),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Covid-Delta-X" ~ scale_y_continuous(limits=c(0,280),breaks=seq(0,280,70),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Covid-Delta-X" ~ scale_y_continuous(limits=c(0,200),breaks=seq(0,200,50),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Covid-Delta-X" ~ scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "Covid-Wildtype-X" ~ scale_y_continuous(limits=c(0,240),breaks=seq(0,240,80),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "Covid-Wildtype-X" ~ scale_y_continuous(limits=c(0,160),breaks=seq(0,160,40),expand=c(0,0),position="right"),
      location == "HIC" & disease == "Covid-Wildtype-X" ~ scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20),expand=c(0,0),position="right"),
      location == "LLMIC" & disease == "SARS-X" ~ scale_y_continuous(limits=c(0,280),breaks=seq(0,280,70),expand=c(0,0),position="right"),
      location == "UMIC" & disease == "SARS-X" ~ scale_y_continuous(limits=c(0,120),breaks=seq(0,120,30),expand=c(0,0),position="right"),
      location == "HIC" & disease == "SARS-X" ~ scale_y_continuous(limits=c(0,60),breaks=seq(0,60,15),expand=c(0,0),position="right"))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.75, "lines")) +
      labs(title = "", x = "VLYL (% of GDP)", y = "GDPL + VSYL (% of GDP)") +
      guides(fill  = guide_legend(title = NULL),
             alpha = "none") + 
      theme(legend.position = "bottom", legend.box.just = "right", 
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9), 
            legend.margin = margin(0, 0, 0, 0))

ggsave("figure_S10.png", plot = gg, height = 14, width = 10)