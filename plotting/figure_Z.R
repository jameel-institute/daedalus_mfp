library(dplyr)
library(whomapper)
library(stringr)
library(ggplot2)
library(ggh4x)
library(ggpattern)
library(patchwork)
source("functions/add_scenario_cols.R")
source("functions/calc_cost_pc.R")
source("functions/find_best_strat.R")
source("functions/calc_cost_bdown.R")

file_list  <- list.files(path = "../output/countries/", pattern = ".*_Influenza 2009_.*\\.csv$", full.names = TRUE)
ctry_data  <- lapply(file_list, add_scenario_cols) %>% bind_rows() %>%
              calc_cost_pc() %>% find_best_strat() %>% filter(best_strat == TRUE) %>% calc_cost_bdown()
ctry_dat1  <- ctry_data %>% distinct(location) %>% slice(1:15) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat2  <- ctry_data %>% distinct(location) %>% slice(16) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat3  <- ctry_data %>% distinct(location) %>% slice(17) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat4  <- ctry_data %>% distinct(location) %>% slice(18) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat5  <- ctry_data %>% distinct(location) %>% slice(19) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat6  <- ctry_data %>% distinct(location) %>% slice(20) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat7  <- ctry_data %>% distinct(location) %>% slice(21) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat8  <- ctry_data %>% distinct(location) %>% slice(22) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat9  <- ctry_data %>% distinct(location) %>% slice(23) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat10 <- ctry_data %>% distinct(location) %>% slice(24) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat11 <- ctry_data %>% distinct(location) %>% slice(25) %>% 
              inner_join(ctry_data, by = "location")
ctry_dat12 <- ctry_data %>% distinct(location) %>% slice(26:40) %>% 
              inner_join(ctry_data, by = "location")

who_adm0      <- pull_sfs(adm_level = 0, query_server = FALSE)
who_adm0$adm0 <- who_adm0$adm0 %>% rename(location = adm0_name) %>% 
                 mutate(location = str_to_title(location),
                        location = case_when(location == "Brunei Darussalam" ~ "Brunei",
                                             location == "Republic Of Korea" ~ "South Korea",
                                             location == "Lao People's Democratic Republic" ~ "Laos",
                                             location == "Netherlands (Kingdom Of The)" ~ "Netherlands",                           
                                             location == "Russian Federation" ~ "Russia",
                                             location == "Türki̇ye" ~ "Turkey",
                                             location == "The United Kingdom" ~ "United Kingdom",
                                             location == "United States Of America" ~ "United States",
                                             location == "Viet Nam" ~ "Vietnam",
                                             TRUE ~ location)) %>%
                 left_join(ctry_data %>% select(location,strategy,mean_SECpc) %>% unique(), by = "location")

gm <- ggplot() +
      geom_sf_pattern(data = who_adm0$adm0, 
                      linewidth = 0.1, col = "black", aes(pattern = strategy, fill = mean_SECpc),
                      pattern_density = 0.35, pattern_size = 0.05, pattern_spacing = 0.008,   
                      pattern_colour = "black", pattern_fill = "white") +
      scale_pattern_manual(values = c("No Closures" = "none", "School Closures" = "crosshatch", 
                                      "Economic Closures" = "circle", "Elimination" = "stripe", "NA" = "none"),
                           breaks = c("No Closures", "School Closures", "Economic Closures", "Elimination"),
                           labels = c("No Closures", "School Closures", "Economic Closures", "Elimination")) +
      scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey80") + #, limits = c(0,850)) +
      geom_sf(data = who_adm0$disp_border, linetype = "dashed", linewidth = 0.1, col = "black", fill = NA) +
      geom_sf(data = who_adm0$disp_area %>% filter(grepl("Lake|Sea",name)), linewidth = 0.1, col = "black", fill = "white") +
      theme_bw() +
      coord_sf(xlim = c(-180,180), ylim = c(-60,90), expand = FALSE) +
      scale_x_continuous(breaks = seq(-180, 180, by = 60)) +  
      scale_y_continuous(breaks = seq(-60, 90, by = 30)) + 
      theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
            axis.text.x = element_blank(), axis.text.y = element_blank()) +
      guides(pattern = guide_legend(title = NULL, order = 1), fill = guide_colorbar(title = "%", order = 2)) + 
      theme(legend.position = c(0.08, 0.28))

gt <- ggplot(ctry_dat1, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")
      #guides(fill = guide_legend(title = NULL)) +
      #theme(legend.position = c(0.05, 0.998), legend.justification = c(1, 1), legend.box.just = "right") 

g2 <- ggplot(ctry_dat2, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

g3 <- ggplot(ctry_dat3, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

g4 <- ggplot(ctry_dat4, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

g5 <- ggplot(ctry_dat5, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

g6 <- ggplot(ctry_dat6, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

g7 <- ggplot(ctry_dat7, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

g8 <- ggplot(ctry_dat8, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

g9 <- ggplot(ctry_dat9, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

g10<- ggplot(ctry_dat10, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

g11<- ggplot(ctry_dat11, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

gb <- ggplot(ctry_dat12, aes(x = strategy, y = SECpc, fill = factor(..fill..))) + 
      facet_grid(cols = vars(location), scales = "free") +
      geom_violin(aes(width = pLYL+pGDP+pSYL, fill = "red"), linewidth = 0.5) + 
      geom_violin(aes(width = pGDP+pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pGDP+pSYL, fill = "yellow"), linewidth = 0.1) +
      geom_violin(aes(width = pSYL, fill = "white"), linewidth = 0) +
      geom_violin(aes(width = pSYL, fill = "blue"), linewidth = 0.1) +
      scale_fill_manual(values = c("red" = "red", "white" = "white", "yellow" = "yellow", "blue" = "blue"),
                        breaks = c("red", "yellow", "blue"), labels = c("VLYL", "GDPL", "VSYL")) +
      theme_bw() + 
      facetted_pos_scales(y = list(scale_y_continuous(position="right", expand=c(0,0)))) +
      theme(panel.spacing = unit(1, "lines"), legend.position = "none") + 
      labs(title = "", x = "", y = "")

patch_left   <- (g2/g4/g6/g8/g10/g10)
patch_right  <- (g3/g5/g7/g9/g11/g11)
patch_middle <- (wrap_elements(patch_left)|wrap_elements(gm)|wrap_elements(patch_right)) + plot_layout(widths = c(1.6,10.8,1.6))
patchwork    <- (wrap_elements(gt)/wrap_elements(patch_middle)/wrap_elements(gb)) + plot_layout(heights = c(1.5,7,1.5))


patchwork    <- (wrap_elements(gt)/wrap_elements(gm)/wrap_elements(gb)) + 
                 plot_layout(heights = c(1.5,6,1.5))

ggsave("figure_4a.png", plot = patchwork, height = 10, width = 14)

###
# width and height
# how to arrange countries?
# area of densities ...
# caxis and ylim scalings for each disease ...