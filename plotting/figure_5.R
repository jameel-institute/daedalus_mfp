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
                (function(x) calc_loss_pc(input_data,x)) 

output_data1 <- output_data %>% 
                filter(disease == "Covid-Wildtype-X") %>%
                group_by(location, country) %>%
                filter(any(elimok == 1)) %>%
                ungroup() %>%
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
                ungroup()

g1 <- ggplot(output_data1, aes(x = strategy, y = SLpc, linewidth = min_any, fill = strategy, alpha = min_any)) + 
      facet_grid2(disease ~ location, switch = "y", scales = "fixed") +
      geom_violin(data = output_data %>% filter(disease == "Covid-Wildtype-X"),
                  aes(x = strategy, y = SLpc),
                  inherit.aes = FALSE, linewidth = 0.1, fill = NA) +
      geom_violin(color = NA, fill = "white", alpha = 1) +
      geom_violin(color = "black") + 
      scale_linewidth_manual(values = c("FALSE" = 0.1, "TRUE" = 0.5)) +
      scale_fill_manual(values = c("No Closures" = "magenta4", "School Closures" = "navy", 
                                   "Economic Closures" = "darkgreen", "Elimination" = "goldenrod"),
                        labels = c("No Closures" = "No Closures",
                                   "School Closures" = "Reactive-Business/Sustained-School Closures",
                                   "Economic Closures" = "Reactive-Business/Reactive-School Closures",
                                   "Elimination" = "Elimination")) + 
      scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 1)) +
      geom_text(data = output_data1 %>% filter(min_med == TRUE), aes(x = strategy, y = max_SLpc),
                vjust = 0.4, label = "*", size = 6, color = "black", inherit.aes = FALSE) +
      geom_text(data = output_data1 %>% filter(min_q3  == TRUE), aes(x = strategy, y = max_SLpc),
                vjust = -1.5, label = "†", size = 3.5, color = "black", inherit.aes = FALSE) +
      theme_bw() + 
      scale_x_discrete(limits = levels(output_data1$strategy), labels = NULL) + 
      scale_y_continuous(limits=c(0,800), breaks=seq(0,800,200), expand=c(0,0), position="right") +
      theme(panel.spacing = unit(2.00, "lines")) + 
      labs(title = "", x = "", y = "Socioeconomic Loss (% of GDP)") +
      guides(linewidth = "none", fill = guide_legend(title = NULL, override.aes = list(colour = NA)), alpha = "none") +
      theme(legend.position = "bottom", legend.box.just = "right",
            legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 9),
            legend.margin = margin(0, 0, 0, 0))

output_data2 <- output_data %>% filter(strategy %in% c("School Closures", "Economic Closures")) %>%
                filter(disease == "Influenza-1918-X") %>%
                mutate(x1 = wdunem, y1 = timldn) %>%
                mutate(x2 = GDPLpc, y2 = VSYLpc)
output_stats <- output_data2 %>% #for quicker plotting
                group_by(location, disease, strategy) %>%
                summarise(mean_x1  = mean(x1),
                          q1_x1    = quantile(x1, 0.25),
                          q3_x1    = quantile(x1, 0.75),
                          mean_y1  = mean(y1),
                          q1_y1    = quantile(y1, 0.25),
                          q3_y1    = quantile(y1, 0.75),
                          mean_x2  = mean(x2),
                          q1_x2    = quantile(x2, 0.25),
                          q3_x2    = quantile(x2, 0.75),
                          mean_y2  = mean(y2),
                          q1_y2    = quantile(y2, 0.25),
                          q3_y2    = quantile(y2, 0.75))

g2 <- ggplot(output_data2, aes(x = x1, y = y1, color = strategy, fill = strategy)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_point(shape = 19, size = 0.02, stroke = 0.2) +
      geom_linerange(data = output_stats, aes(x = mean_x1, y = mean_y1, xmin = q1_x1, xmax = q3_x1), linewidth = 0.25, color = "black") +
      geom_linerange(data = output_stats, aes(x = mean_x1, y = mean_y1, ymin = q1_y1, ymax = q3_y1), linewidth = 0.25, color = "black") +
      geom_point(data = output_stats, aes(x = mean_x1, y = mean_y1), shape = 21, size = 2, stroke = 0.2, color = "black") +
      scale_color_manual(values = c("School Closures" = "navy", "Economic Closures" = "darkgreen"),
                         labels = c("School Closures"   = "Reactive-Business/Sustained-School Closures",
                                    "Economic Closures" = "Reactive-Business/Reactive-School Closures")) +
      scale_fill_manual(values = c("School Closures" = "navy", "Economic Closures" = "darkgreen"),
                        labels = c("School Closures"   = "Reactive-Business/Sustained-School Closures",
                                   "Economic Closures" = "Reactive-Business/Reactive-School Closures")) +
      theme_bw() +
      facetted_pos_scales(
        x = list(location == "LLMIC" ~ scale_x_continuous(limits=c(0,280),breaks=seq(0,280,70),expand=c(0,0),position="bottom"),
                 location == "UMIC"  ~ scale_x_continuous(limits=c(0,200),breaks=seq(0,200,50),expand=c(0,0),position="bottom"),
                 location == "HIC"   ~ scale_x_continuous(limits=c(0,120),breaks=seq(0,120,30),expand=c(0,0),position="bottom")),
        y = list(location == "LLMIC" ~ scale_y_continuous(limits=c(0,800),breaks=seq(0,800,200),expand=c(0,0),position="right"),
                 location == "UMIC"  ~ scale_y_continuous(limits=c(0,600),breaks=seq(0,600,150),expand=c(0,0),position="right"),
                 location == "HIC"   ~ scale_y_continuous(limits=c(0,400),breaks=seq(0,400,100),expand=c(0,0),position="right"))) +
      theme(panel.spacing = unit(0.75, "lines")) +
      labs(title = "", x = "Furlough-Days per Worker", y = "Heavy-Closure Days") +
      theme(legend.position = "none")

g3 <- ggplot(output_data2, aes(x = x2, y = y2, color = strategy, fill = strategy)) +
      facet_grid2(disease ~ location, switch = "y", scales = "free", independent = "all") +
      geom_point(shape = 19, size = 0.02, stroke = 0.2) +
      geom_linerange(data = output_stats, aes(x = mean_x2, y = mean_y2, xmin = q1_x2, xmax = q3_x2), linewidth = 0.25, color = "black") +
      geom_linerange(data = output_stats, aes(x = mean_x2, y = mean_y2, ymin = q1_y2, ymax = q3_y2), linewidth = 0.25, color = "black") +
      geom_point(data = output_stats, aes(x = mean_x2, y = mean_y2), shape = 21, size = 2, stroke = 0.2, color = "black") +
      scale_color_manual(values = c("School Closures" = "navy", "Economic Closures" = "darkgreen"),
                         labels = c("School Closures"   = "Reactive-Business/Sustained-School Closures",
                                    "Economic Closures" = "Reactive-Business/Reactive-School Closures")) +
      scale_fill_manual(values = c("School Closures" = "navy", "Economic Closures" = "darkgreen"),
                        labels = c("School Closures"   = "Reactive-Business/Sustained-School Closures",
                                   "Economic Closures" = "Reactive-Business/Reactive-School Closures")) +
      theme_bw() +
      facetted_pos_scales(
        x = list(location == "LLMIC" ~ scale_x_continuous(limits=c(0,120),breaks=seq(0,120,30),expand=c(0,0),position="bottom"),
                 location == "UMIC"  ~ scale_x_continuous(limits=c(0,80),breaks=seq(0,80,20),expand=c(0,0),position="bottom"),
                 location == "HIC"   ~ scale_x_continuous(limits=c(0,40),breaks=seq(0,40,10),expand=c(0,0),position="bottom")),
        y = list(location == "LLMIC" ~ scale_y_continuous(limits=c(0,200),breaks=seq(0,200,50),expand=c(0,0),position="right"),
                 location == "UMIC"  ~ scale_y_continuous(limits=c(0,120),breaks=seq(0,120,30),expand=c(0,0),position="right"),
                 location == "HIC"   ~ scale_y_continuous(limits=c(0,40),breaks=seq(0,40,10),expand=c(0,0),position="right"))) +
      theme(panel.spacing = unit(0.75, "lines")) +
      labs(title = "", x = "GDPL (% of GDP)", y = "VSYL (% of GDP)") +
      theme(legend.position = "none")

output_data3 <- output_data %>% filter(strategy == "No Closures")
ctry_data    <- read.csv("../input/country_data.csv") %>%
                mutate(igroup = factor(igroup, levels = c("LLMIC","UMIC","HIC"))) %>%
                dplyr::select(igroup,Hmax) %>% filter(!is.na(Hmax)) %>%
                mutate(candidate = list(c("lnorm", "gamma", "weibull"))) %>% unnest_longer(candidate) %>%
                group_by(igroup, candidate) %>%
                summarize(fit      = list(tryCatch(fitdist(Hmax, unique(candidate), method = "mle"), error = function(e) NA)),
                          meanlog  = sapply(fit, function(f) tryCatch(f$estimate["meanlog"], error = function(e) NA)),
                          sdlog    = sapply(fit, function(f) tryCatch(f$estimate["sdlog"], error = function(e) NA)),
                          shape    = sapply(fit, function(f) tryCatch(f$estimate["shape"], error = function(e) NA)),
                          rate     = sapply(fit, function(f) tryCatch(f$estimate["rate"], error = function(e) NA)),
                          scale    = sapply(fit, function(f) tryCatch(f$estimate["scale"], error = function(e) NA)),
                          aic      = sapply(fit, function(f) tryCatch(f$aic, error = function(e) NA)), .groups = "drop") %>%
                group_by(igroup) %>%
                slice_min(aic, n = 1, with_ties = TRUE) %>%
                ungroup() %>%
                transmute(location = igroup,
                          ymin     = qlnorm(0.25, meanlog, sdlog),
                          ymax     = qlnorm(0.75, meanlog, sdlog))

g4 <- ggplot(output_data3, aes(x = disease, y = hpeak)) +
      facet_grid2(~ location, switch = "y", scales = "fixed") +
      geom_hline(data = ctry_data, aes(yintercept = ymin), linetype = "dashed", linewidth = 0.25, color = "black") +
      geom_hline(data = ctry_data, aes(yintercept = ymax), linetype = "dashed", linewidth = 0.25, color = "black") +
      geom_boxplot(outlier.shape = NA, coef = 1.5, width = 0.70, linewidth = 0.2, fill = "magenta4") +
      theme_bw() +
      scale_y_log10(limits = c(0.3,3000), breaks = c(0.3,3,30,300,3000), expand = c(0,0), position = "right",
                    labels = function(x) ifelse(x < 1, as.character(x), as.character(round(x)))) +
      theme(panel.spacing = unit(2.00, "lines"), axis.text.x = element_text(angle = 55, hjust = 1)) +
      labs(title = "", x = "", y = "Peak Hospital Occupancy (per 100k)") +
      theme(legend.position = "none")

lg <- get_legend(g1)
gg <- (g1 + theme(legend.position="none")) /
      (g2 + theme(legend.position="none")) /
      (g3 + theme(legend.position="none")) /
      (g4 + theme(legend.position="none")) +
      plot_annotation(tag_levels='A')
gg <- plot_grid(as_grob(gg), lg, ncol = 1, rel_heights = c(1,0.08))
ggsave("figure_5.png", plot = gg, height = 14, width = 10)