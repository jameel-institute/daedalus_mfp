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
                pivot_longer(cols = starts_with("chad_"), names_to = "chad", values_to = "value") %>%
                mutate(chad = factor(chad, levels = c("chad_1", "chad_2", "chad_3", "chad_4")))

gg <- ggplot(output_data, aes(x = strategy, y = value, fill = chad)) + 
      facet_grid2(disease ~ location, switch = "y", scales = "free_y") +
      geom_boxplot(outlier.shape = NA, coef = 1.5, width = 0.70, linewidth = 0.2) +
      scale_fill_manual(values = c("chad_1" = "magenta", "chad_2" = "green", "chad_3" = "brown", "chad_4" = "grey80"),
                        labels = c("Preschool-Age", "School-Age", "Working-Age", "Retired-Age")) +
      theme_bw() + 
      scale_x_discrete(labels = c("School Closures" = "Reactive-Business/Sustained-School Closures",
                                  "Economic Closures" = "Reactive-Business/Reactive-School Closures")) +
      facetted_pos_scales(y = list(
        scale_y_continuous(limits=c(0,1000),  expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,8000),  expand=c(0,0), position="right", breaks = seq(0,8000,by=2000)),
        scale_y_continuous(limits=c(0,10000), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,10000), expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,28000), expand=c(0,0), position="right", breaks = seq(0,28000,by=7000)),
        scale_y_continuous(limits=c(0,10000), expand=c(0,0), position="right"))) +
          theme(panel.spacing = unit(0.75, "lines"), axis.text.x = element_text(angle = 55, hjust = 1)) + 
      labs(title = "", x = "", y = "Cumulative Hospital Admissions (per 100k)") +
      guides(fill = guide_legend(title = NULL)) +
      theme(legend.position = "bottom", legend.box.just = "right", 
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9), 
            legend.margin = margin(0, 0, 0, 0))

ggsave("figure_S13.png", plot = gg, height = 14, width = 10)