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
output_files <- list_files[!grepl("_data\\.csv$", list_files)]
output_data  <- lapply(output_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>% 
                filter(strategy %in% c("School Closures","Economic Closures")) %>%
                filter(!hthres == "[]") %>%
                mutate(hthres = gsub("\\[|\\]", "", hthres)) %>%
                separate_longer_delim(hthres, delim = ",") %>%
                mutate(hthres = as.numeric(hthres))

gg <- ggplot(output_data, aes(x = strategy, y = hthres, fill = strategy)) + 
      facet_grid2(disease ~ location, switch = "y", scales = "fixed") +
      geom_violin() + 
      scale_fill_manual(values = c("School Closures" = "navy", "Economic Closures" = "darkgreen")) + 
      theme_bw() +
      scale_x_discrete(labels = c("School Closures" = "Reactive-Business/Sustained-School Closures",
                                  "Economic Closures" = "Reactive-Business/Reactive-School Closures")) +
      scale_y_continuous(limits = c(0,1), expand = c(0,0), position="right", labels = scales::label_number(scale = 100, suffix = "")) +
      theme(panel.spacing = unit(0.75, "lines"), axis.text.x = element_text(angle = 65, hjust = 1)) + 
      labs(title = "", x = "", y = "Hospital Occupancy at Imposition of Heavy Closures (% of Spare Hospital Beds)") +
      guides(fill = "none")

ggsave("figure_S6.png", plot = gg, height = 14, width = 10)
ggsave("figure_S6.pdf", plot = gg, height = 14, width = 10)