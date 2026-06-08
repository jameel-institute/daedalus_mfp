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

list_files   <- list.files(path = "../input/r0rt/", pattern = "\\.csv$", full.names = TRUE)
output_files <- list_files[!grepl("_data\\.csv$", list_files)]
output_data  <- lapply(output_files, function(path) {
                  title <- sub("\\.csv$", "", basename(path))
                  parts <- unlist(strsplit(title, "_"))
                  df    <- read.csv(path) %>% mutate(location = parts[1], disease = parts[2], strategy = parts[4])}) %>%
                bind_rows() %>% 
                mutate(location = factor(location, levels = c("LLMIC","UMIC","HIC")),
                       disease  = gsub(" ", "-", disease) %>% paste0("-X"),
                       disease  = factor(disease, levels = c("Influenza-2009-X","Influenza-1957-X","Influenza-1918-X",
                                                             "Covid-Omicron-X","Covid-Delta-X","Covid-Wildtype-X")),
                       strategy = factor(strategy)) %>%
                mutate(rt = 1 - (rt/r0))

gg <- ggplot(output_data, aes(x = strategy, y = rt, fill = strategy)) + 
      facet_grid2(disease ~ location, switch = "y", scales = "fixed") +
      geom_violin() + 
      scale_fill_manual(values = c("1" = "navy", "2" = "navy", "3" = "darkgreen", "4" = "darkgreen", "5" = "goldenrod")) +
      theme_bw() + 
      scale_x_discrete(labels = c("1" = "RBSS: Heavy EC, SC, WFH, Low ASC", 
                                  "2" = "RBSS: Light EC, SC, WFH, Low ASC", 
                                  "3" = "RBRS: Heavy EC, SC, WFH, Low ASC", 
                                  "4" = "RBRS: Light EC, WFH, Low ASC", 
                                  "5" = "ELIM: Light EC, WFH, High ASC")) +
      scale_y_continuous(limits = c(0,1), expand = c(0,0), position="right", labels = scales::label_number(scale = 100, suffix = "")) +
      theme(panel.spacing = unit(0.75, "lines"), axis.text.x = element_text(angle = 75, hjust = 1)) + 
      labs(title = "", x = "", y = expression("Relative Reduction in " * R[0] * " (%)")) +
      guides(fill = "none")

ggsave("figure_S5.png", plot = gg, height = 14, width = 10)