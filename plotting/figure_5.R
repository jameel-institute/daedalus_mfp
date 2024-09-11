library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
sapply(list.files(path = "functions/voi-master/R/", pattern = "\\.R$", full.names = TRUE), source)
library(ggplot2)
library(ggh4x)
library(cowplot)
source("functions/add_scenario_cols.R")
source("functions/order_scenario_cols.R")
source("functions/calc_cost_pc.R")
source("functions/find_best_strats.R")
source("functions/calc_cost_bdown.R")
source("functions/parse_inputs.R")
source("functions/voi_est.R")
source("functions/voi_est_fit.R")
source("functions/voi_dec.R")

file_list <- list.files(path = "../output/archetypes/", pattern = "\\.csv$", full.names = TRUE)
arch_data <- lapply(file_list, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>% 
             calc_cost_pc() %>% parse_inputs() 
arch_voi  <- voi_dec(arch_data, list(c("mean_age", "sd_age", "skew_age", "le"),
                                     c("comm", "IGHC", "schoolA2", "workp"),
                                     c("EPOP", "TES", "gvapw", "wfh"),
                                     c("Tres", "t_tit", "trate"),
                                     c("sdb", "sdl"),
                                     c("Hmax"),
                                     c("t_vax", "arate", "puptake"))) %>%
             mutate(parameter = case_when(parameter == "mean_age,sd_age,skew_age,le" ~ "Demography",
                                          parameter == "comm,IGHC,schoolA2,workp" ~ "Mixing",
                                          parameter == "EPOP,TES,gvapw,wfh" ~ "Economy",
                                          parameter == "Tres,t_tit,trate" ~ "Surveillance",
                                          parameter == "sdb,sdl" ~ "Distancing",
                                          parameter == "Hmax" ~ "Hospital Capacity",
                                          parameter == "t_vax,arate,puptake" ~ "Vaccination")) %>%
             mutate(parmeter = factor(parameter, levels = c("Surveillance", "Distancing", "Hospital Capacity", "Vaccination", 
                                                            "Demography", "Mixing", "Economy")))

gg <- ggplot(arch_voi, aes(x = rank, y = res, fill = parameter)) +  
      facet_grid(disease ~ location, switch = "y", scales = "free") +
      geom_bar(stat = "identity", position = "dodge", linewidth = 0.3, color = "black") +
      scale_fill_manual(values = c("Demography" = "black", "Mixing" = "white", "Economy" = "goldenrod",
                                   "Surveillance" = "palegreen", "Distancing" = "slategray2", 
                                   "Hospital Capacity" = "purple", "Vaccination" = "lightsalmon")) +
      theme_bw() + 
      scale_y_continuous(expand=c(0,0), position="right") + 
      theme(panel.spacing = unit(1, "lines")) + 
      labs(title = "", x = "", y = "EVPPI (% of GDP)") +
      guides(fill = guide_legend(title = NULL, ncol = 2, byrow = FALSE), alpha = "none") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),  
            legend.position = c(0.229, 0.9996), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 7))

ggsave("figure_5.png", plot = gg, height = 14, width = 10)