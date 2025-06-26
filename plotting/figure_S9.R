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
source("functions/voi_fit.R")
#source("functions/format_table.R")

list_files   <- list.files(path = "../output/archetypes/", pattern = "\\.csv$", full.names = TRUE)
input_files  <- list_files[grepl("_data\\.csv$", list_files)]
input_data   <- lapply(input_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols()
output_files <- list_files[!grepl("_data\\.csv$", list_files)]
output_data  <- lapply(output_files, add_scenario_cols) %>% bind_rows() %>% order_scenario_cols() %>%
                left_join(input_data %>% dplyr::select(location, country, Hmax), by = c("location", "country")) %>%
                (function(x) calc_loss_pc(input_data,x)) %>% 
                mutate(x = VLYLpc, y = GDPLpc + VSYLpc) %>%
                group_by(location, country, disease) %>%
                mutate(x = x[strategy == "No Closures"] - x,
                       y = y - y[strategy == "No Closures"]) %>%
                ungroup() %>%
                filter(strategy == "School Closures") %>%
                mutate(parameter = "Hmax", xaxis = Hmax)
output_fitX  <- voi_fit(output_data, x) %>% rename(xaxis = xval, yaxis = yval, lower = l, upper = u) %>% mutate(cost_eff = "x")
output_fitY  <- voi_fit(output_data, y) %>% rename(xaxis = xval, yaxis = yval, lower = l, upper = u) %>% mutate(cost_eff = "y")
output_data  <- output_data %>% pivot_longer(cols = c("x","y"), names_to = "cost_eff", values_to = "yaxis")

gg <- ggplot(output_data, aes(x = xaxis, y = yaxis, color = cost_eff, fill = cost_eff)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_ribbon(data = output_fitX, 
                  aes(ymin = lower, ymax = upper), color = NA, alpha = 0.25) +
      geom_ribbon(data = output_fitY, 
                  aes(ymin = lower, ymax = upper), color = NA, alpha = 1) +
      geom_point(shape = 19, size = 0.02, stroke = 0.2) +
      geom_line(data = output_fitX, linewidth = 0.5) +
      geom_line(data = output_fitY, linewidth = 0.5) +
      scale_color_manual(values = c("x" = "red", "y" = "blue"), 
                         labels = c("Decremental VLYL (% of GDP)", "Incremental GDPL + VSYL (% of GDP)")) +
      scale_fill_manual(values = c("x" = "red", "y" = "yellow"),
                        labels = c("Decremental VLYL (% of GDP)", "Incremental GDPL + VSYL (% of GDP)")) +
      theme_bw() +
      facetted_pos_scales(
        x = list(
        location == "LLMIC" ~ scale_x_log10(limits=c(1,300),breaks=c(1,3,10,30,100,300),expand=c(0,0),position="bottom"),
        location == "UMIC" ~ scale_x_log10(limits=c(3,300),breaks=c(3,10,30,100,300),expand=c(0,0),position="bottom"),
        location == "HIC" ~ scale_x_log10(limits=c(10,1000),breaks=c(10,30,100,300,1000),expand=c(0,0),position="bottom")),
        y = list (
        location == "LLMIC" & disease == "Influenza-2009-X" ~ scale_y_continuous(limits=c(-50,150),breaks=seq(-50,150,50),expand=c(0,0),position="right"),
        location == "UMIC" & disease == "Influenza-2009-X" ~ scale_y_continuous(limits=c(-20,60),breaks=seq(-20,60,20),expand=c(0,0),position="right"),
        location == "HIC" & disease == "Influenza-2009-X" ~ scale_y_continuous(limits=c(-10,30),breaks=seq(-10,30,10),expand=c(0,0),position="right"),
        location == "LLMIC" & disease == "Influenza-1957-X" ~ scale_y_continuous(limits=c(-50,150),breaks=seq(-50,150,50),expand=c(0,0),position="right"),
        location == "UMIC" & disease == "Influenza-1957-X" ~ scale_y_continuous(limits=c(-30,90),breaks=seq(-30,90,30),expand=c(0,0),position="right"),
        location == "HIC" & disease == "Influenza-1957-X" ~ scale_y_continuous(limits=c(-15,45),breaks=seq(-15,45,15),expand=c(0,0),position="right"),
        location == "LLMIC" & disease == "Influenza-1918-X" ~ scale_y_continuous(limits=c(-200,600),breaks=seq(-200,600,200),expand=c(0,0),position="right"),
        location == "UMIC" & disease == "Influenza-1918-X" ~ scale_y_continuous(limits=c(-200,600),breaks=seq(-200,600,200),expand=c(0,0),position="right"),
        location == "HIC" & disease == "Influenza-1918-X" ~ scale_y_continuous(limits=c(-200,600),breaks=seq(-200,600,200),expand=c(0,0),position="right"),
        location == "LLMIC" & disease == "Covid-Omicron-X" ~ scale_y_continuous(limits=c(-60,180),breaks=seq(-60,180,60),expand=c(0,0),position="right"),
        location == "UMIC" & disease == "Covid-Omicron-X" ~ scale_y_continuous(limits=c(-40,120),breaks=seq(-30,120,40),expand=c(0,0),position="right"),
        location == "HIC" & disease == "Covid-Omicron-X" ~ scale_y_continuous(limits=c(-20,60),breaks=seq(-20,60,20),expand=c(0,0),position="right"),
        location == "LLMIC" & disease == "Covid-Delta-X" ~ scale_y_continuous(limits=c(-100,300),breaks=seq(-100,300,100),expand=c(0,0),position="right"),
        location == "UMIC" & disease == "Covid-Delta-X" ~ scale_y_continuous(limits=c(-100,300),breaks=seq(-100,300,100),expand=c(0,0),position="right"),
        location == "HIC" & disease == "Covid-Delta-X" ~ scale_y_continuous(limits=c(-100,300),breaks=seq(-100,300,100),expand=c(0,0),position="right"),
        location == "LLMIC" & disease == "Covid-Wildtype-X" ~ scale_y_continuous(limits=c(-80,240),breaks=seq(-80,240,80),expand=c(0,0),position="right"),
        location == "UMIC" & disease == "Covid-Wildtype-X" ~ scale_y_continuous(limits=c(-80,240),breaks=seq(-80,240,80),expand=c(0,0),position="right"),
        location == "HIC" & disease == "Covid-Wildtype-X" ~ scale_y_continuous(limits=c(-80,240),breaks=seq(-80,240,80),expand=c(0,0),position="right"),
        location == "LLMIC" & disease == "SARS-X" ~ scale_y_continuous(limits=c(-400,1200),breaks=seq(-400,1200,400),expand=c(0,0),position="right"),
        location == "UMIC" & disease == "SARS-X" ~ scale_y_continuous(limits=c(-400,1200),breaks=seq(-400,1200,400),expand=c(0,0),position="right"),
        location == "HIC" & disease == "SARS-X" ~ scale_y_continuous(limits=c(-200,600),breaks=seq(-200,600,200),expand=c(0,0),position="right"))) +
      theme(panel.spacing = unit(0.75, "lines")) +
      labs(title = "", x = "Spare Hospital Beds (per 100k)", y = "") +
      guides(color = guide_legend(title = NULL),
             fill  = guide_legend(title = NULL)) +
      theme(legend.position = "bottom", legend.box.just = "right",
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9),
            legend.margin = margin(0, 0, 0, 0))

ggsave("figure_S9.png", plot = gg, height = 14, width = 10)