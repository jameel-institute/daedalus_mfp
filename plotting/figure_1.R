#blah
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
             calc_cost_pc() %>% #find_best_strats() %>% 
             calc_cost_bdown() 

gg <- ggplot(arch_data, aes(x = location, y = SECpc, fill = factor(..fill..)), alpha = 1) + 
      facet_grid(disease ~ strategy, switch = "y", scales = "free_y") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.3) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      #scale_linewidth_manual(values = c("FALSE" = 0.1, "TRUE" = 0.5)) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"), 
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      #scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 1)) +
      theme_bw() + 
      facetted_pos_scales(y = list(
        scale_y_continuous(limits=c(0,240),  breaks=seq(0,240,80),    expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,270),  breaks=seq(0,270,90),    expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,1800), breaks=seq(0,1800,600),  expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,360),  breaks=seq(0,360,120),   expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,1500), breaks=seq(0,1500,500),  expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,690),  breaks=seq(0,690,230),   expand=c(0,0), position="right"),
        scale_y_continuous(limits=c(0,3000), breaks=seq(0,3000,1000), expand=c(0,0), position="right"))) +
      theme(panel.spacing = unit(1, "lines"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "", y = "Societal Loss (% of GDP)") +
      guides(width = "none", fill = guide_legend(title = NULL)) +
      theme(legend.position = c(0.078, 0.999), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 7))

ggsave("figure_1.png", plot = gg, height = 14, width = 10)

# gg <- ggplot(arch_data, aes(x = strategy, y = SECpc, fill = factor(..fill..), alpha = min_mean)) +
#       facet_grid(disease ~ location, scales = "free_y", switch = "y") +
#       geom_violin(aes(width = pLYL+pGDP+pSYL, linewidth = min_mean, fill = "red")) +
#       geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0, alpha = 1) +
#       geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
#       geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0, alpha = 1) +
#       geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
#       scale_linewidth_manual(values = c("FALSE" = 0.1, "TRUE" = 0.5)) +
#       scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
#                         breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
#       scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 1)) +
#       geom_text(data = arch_data %>% filter(min_med == TRUE), aes(x = strategy, y = max_SECpc),
#                 vjust = 0.35, label = "*", size = 6, color = "black", inherit.aes = FALSE) +
#       geom_text(data = arch_data %>% filter(min_q3  == TRUE), aes(x = strategy, y = max_SECpc),
#                 vjust = -1.65, label = "†", size = 3.5, color = "black", inherit.aes = FALSE) +
#       theme_bw() +
#       facetted_pos_scales(y = list(
#         scale_y_continuous(position="right", limits=c(0,240), breaks=seq(0,240,80),   expand=c(0,0)),
#         scale_y_continuous(position="right", limits=c(0,270), breaks=seq(0,270,90),   expand=c(0,0)),
#         scale_y_continuous(position="right", limits=c(0,1800),breaks=seq(0,1800,600), expand=c(0,0)),
#         scale_y_continuous(position="right", limits=c(0,360), breaks=seq(0,360,120),  expand=c(0,0)),
#         scale_y_continuous(position="right", limits=c(0,1500),breaks=seq(0,1500,500), expand=c(0,0)),
#         scale_y_continuous(position="right", limits=c(0,690), breaks=seq(0,690,230),  expand=c(0,0)),
#         scale_y_continuous(position="right", limits=c(0,3000),breaks=seq(0,3000,1000),expand=c(0,0)))) +
#       theme(panel.spacing = unit(1, "lines"), axis.text.x = element_text(angle = 45, hjust = 1)) +
#       labs(title = "", x = "", y = "Societal Loss (% of GDP)") +
#       guides(width = "none", linewidth = "none", fill = guide_legend(title = NULL), alpha = "none") +
#       theme(legend.position = c(0.078, 0.999), legend.justification = c(1, 1), legend.box.just = "right", 
#             legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 7))
# 
# ggsave("figure_1alt.png", plot = gg, height = 14, width = 10)

# #alternative log-KDE calculation for violin plots
#
# log_density_df <- arch_data %>%
#   group_by(location, disease, strategy) %>%
#   filter(!is.na(SECpc)) %>% # Filter out missing values
#   do({
#     density_data <- logdensity(.$SECpc)
#     pdf_x <- density_data$x
#     pdf_y <- density_data$y
#     augmented_data <- data.frame(
#       pdf_x = c(pdf_x, rev(pdf_x)),
#       pdf_y = c(-pdf_y, rev(pdf_y)))
#     augmented_data
#   }) %>%
#   ungroup() %>%
#   group_by(location, disease) %>%
#   mutate(pdf_y = as.numeric(strategy) + 0.5*pdf_y/max(pdf_y, na.rm = TRUE)) %>%
#   ungroup()   #mutate(pdf_y = pdf_y + as.numeric(strategy))
# 
# gg <- ggplot(log_density_df, aes(x = pdf_y, y = pdf_x,  group = interaction(location, disease, strategy), alpha = 0.5)) + 
#     geom_polygon(linewidth = 0.5, color = "black", fill = "black") +
#     #coord_flip() +
#   facet_grid(rows = vars(disease), cols = vars(location), switch = "y", scales = "free_y") +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",  panel.spacing = unit(1, "lines")) + 
#   facetted_pos_scales(y = list(
#     scale_y_continuous(position="right", limits=c(0,300), expand=c(0,0), breaks=seq(0,300,100)),
#     scale_y_continuous(position="right", limits=c(0,300), expand=c(0,0), breaks=seq(0,300,100)),
#     scale_y_continuous(position="right", limits=c(0,1800), expand=c(0,0), breaks=seq(0,1800,600)),
#     scale_y_continuous(position="right", limits=c(0,750), expand=c(0,0), breaks=seq(0,750,250)),
#     scale_y_continuous(position="right", limits=c(0,1500), expand=c(0,0), breaks=seq(0,1500,500)),
#     scale_y_continuous(position="right", limits=c(0,450), expand=c(0,0), breaks=seq(0,450,150)),
#     scale_y_continuous(position="right", limits=c(0,3000), expand=c(0,0), breaks=seq(0,3000,1000)))) +
#   labs(title = "", x = "", y = "Societal Cost (% of Annual GDP)") + 
#   scale_x_continuous(breaks = 1:4, labels = c("No Closures", "School Closures", "Economic Closures", "Elimination"))