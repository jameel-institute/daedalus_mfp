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
#source("functions/add_scenario_cols.R")
#source("functions/order_scenario_cols.R")
#source("functions/calc_cost_pc.R")
#source("functions/parse_inputs.R")
#source("functions/voi_dec.R")
#source("functions/voi_fit.R")
#source("functions/table_formatting.R")

ctry_data <- read.csv("../input/country_data.csv") %>%
             mutate(igroup   = factor(igroup, levels = c("LLMIC","UMIC","HIC"))) %>%
             #demography
             mutate(pop      = rowSums(across(starts_with("Npop")))) %>%
             mutate(across(starts_with("Npop"), ~ .x/pop)) %>%
             #mixing
             mutate(pa_prop  = rowSums(across(Npop_1:Npop_1)),
                    sa_prop  = rowSums(across(Npop_2:Npop_4)),
                    wa_prop  = rowSums(across(Npop_5:Npop_13)),
                    ra_prop  = rowSums(across(Npop_14:Npop_21))) %>%
             rowwise() %>%
             mutate(matAL    = list(matrix(unlist(c_across(matAL_1:matAL_16)), nrow = 4, byrow = FALSE)),
                    AL_wavg  = sum(rowSums(matAL) * c_across(c(pa_prop, sa_prop, wa_prop, ra_prop))),
                    matAHT   = list(matrix(unlist(c_across(matAHT_1:matAHT_16)), nrow = 4, byrow = FALSE)),
                    AHT_wavg = sum(rowSums(matAHT) * c_across(c(pa_prop, sa_prop, wa_prop, ra_prop))),
                    matAS    = list(matrix(unlist(c_across(matAS_1:matAS_16)), nrow = 4, byrow = FALSE)),
                    AS_wavg  = sum(matAS[2,])) %>%
             ungroup() %>%
             #economy
             mutate(wa_pop   = pop*wa_prop) %>%
             mutate(across(starts_with("NNs"), ~ .x/wa_pop)) %>%
             mutate(NNs_46   = 1 - rowSums(across(starts_with("NNs")))) %>%
             mutate(across(obj_1:obj_45, ~ ifelse(wa_pop * get(sub("obj", "NNs", cur_column())) == 0, 0, 
                                           365*10^6* .x / (wa_pop * get(sub("obj", "NNs", cur_column()))))))
npop_data <- ctry_data %>% dplyr::select(igroup, country, starts_with("Npop_")) %>%
             pivot_longer(cols = starts_with("Npop_"), names_to = "Npop", values_to = "value") %>%
             mutate(Npop = parse_number(Npop))
lfex_data <- ctry_data %>% dplyr::select(igroup, country, starts_with("la_")) %>%
             pivot_longer(cols = starts_with("la_"), names_to = "la", values_to = "value") %>%
             mutate(la = parse_number(la))
nns_data  <- ctry_data %>% dplyr::select(igroup, country, starts_with("NNs_")) %>%
             pivot_longer(cols = starts_with("NNs_"), names_to = "NNs", values_to = "value") %>%
             mutate(NNs = parse_number(NNs)) %>%
             mutate(value = ifelse(value < 0.0001, 0.0001, value))
gvaw_data <- ctry_data %>% dplyr::select(igroup, country, starts_with("obj_")) %>%
             pivot_longer(cols = starts_with("obj_"), names_to = "obj", values_to = "value") %>%
             mutate(obj = parse_number(obj)) %>%
             filter(value != 0)#GVA is undefined if no workers in sector
wfh_data  <- ctry_data %>% dplyr::select(igroup, country, starts_with("wfh_")) %>%
             pivot_longer(cols = starts_with("wfh_"), names_to = "wfh", values_to = "value") %>%
             mutate(wfh = parse_number(wfh)) %>%
             mutate(value = ifelse(value < 0.0001, 0.0001, value))

p1 <- ggplot(npop_data, aes(x = Npop, y = value, group = country, color = igroup)) +
      facet_wrap(~ "Population by Age") +
      geom_line(linewidth = 0.5, alpha = 0.50) +
      scale_color_manual(values = c("LLMIC" = "orange", "UMIC" = "turquoise", "HIC" = "azure4")) +  
      theme_bw() +
      scale_x_continuous(breaks = 1:21, expand = c(0,0), position = "bottom",
                         labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", 
                                    "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99", "100+")) + 
      scale_y_continuous(limits = c(0,0.20), expand = c(0,0), position = "left", labels = scales::label_number(scale = 100, suffix = "")) +
      theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "Age Group", y = "Percentage of Population (%)") +
      guides(color = guide_legend(title = NULL)) +
      theme(legend.position = c(0.999, 0.995), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

p2 <- ggplot(lfex_data, aes(x = la, y = value, group = country, color = igroup)) +
      facet_wrap(~ "Life Expectancy by Age") +
      geom_line(linewidth = 0.5, alpha = 0.50) +
      scale_color_manual(values = c("LLMIC" = "orange", "UMIC" = "turquoise", "HIC" = "azure4")) +  
      theme_bw() +
      scale_x_continuous(breaks = 1:18, expand = c(0,0), position = "bottom",
                         labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", 
                                    "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+")) + 
      scale_y_continuous(limits = c(0,90), breaks = seq(0,90,30), expand = c(0,0), position = "left") +
      theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "Age Group", y = "Remaining Life Expectancy (years)") +
      guides(color = guide_legend(title = NULL)) +
      theme(legend.position = c(0.999, 0.995), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

p3 <- ggplot(ctry_data, aes(x = AL_wavg, fill = igroup)) +
      facet_wrap(~ "Household Contacts") +
      geom_histogram(aes(y = ..density..), position = "stack", color = "black") + 
      scale_fill_manual(values = c("LLMIC" = "orange", "UMIC" = "turquoise", "HIC" = "azure4")) +  
      theme_bw() +
      scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2), expand = c(0,0), position = "bottom") +
      scale_y_continuous(limits = c(0,3), breaks = seq(0,3,1), expand = c(0,0), position = "left", labels = scales::label_parse()) +
      # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "Population-Average Household Contacts (#/person/day)", y = "Relative Frequency") +
      guides(fill = guide_legend(title = NULL)) +
      theme(legend.position = c(0.202, 0.995), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

p4 <- ggplot(ctry_data, aes(x = AHT_wavg, fill = igroup)) +
      facet_wrap(~ "Other-Location Contacts") +
      geom_histogram(aes(y = ..density..), position = "stack", color = "black") + 
      scale_fill_manual(values = c("LLMIC" = "orange", "UMIC" = "turquoise", "HIC" = "azure4")) +  
      theme_bw() +
      scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2), expand = c(0,0), position = "bottom") +
      scale_y_continuous(limits = c(0,1.5), breaks = seq(0,1.5,0.5), expand = c(0,0), position = "left", labels = scales::label_parse()) +
      # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "Population-Average Other-Location Contacts (#/person/day)", y = "Relative Frequency") +
      guides(fill = guide_legend(title = NULL)) +
      theme(legend.position = c(0.202, 0.995), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

p5 <- ggplot(ctry_data, aes(x = AS_wavg, fill = igroup)) +
      facet_wrap(~ "School Contacts") +
      geom_histogram(aes(y = ..density..), position = "stack", color = "black") + 
      scale_fill_manual(values = c("LLMIC" = "orange", "UMIC" = "turquoise", "HIC" = "azure4")) +  
      theme_bw() +
      scale_x_continuous(limits = c(0,15), breaks = seq(0,15,3), expand = c(0,0), position = "bottom") +
      scale_y_continuous(limits = c(0,1.5), breaks = seq(0,1.5,0.5), expand = c(0,0), position = "left", labels = scales::label_parse()) +
      # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "School-Age School Contacts (#/person/day)", y = "Relative Frequency") +
      guides(fill = guide_legend(title = NULL)) +
      theme(legend.position = c(0.202, 0.995), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

p6 <- ggplot(ctry_data, aes(x = workp, fill = igroup)) +
      facet_wrap(~ "Workplace Contacts") +
      geom_histogram(aes(y = ..density..), position = "stack", color = "black") + 
      scale_fill_manual(values = c("LLMIC" = "orange", "UMIC" = "turquoise", "HIC" = "azure4")) +  
      theme_bw() +
      scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2), expand = c(0,0), position = "bottom") +
      scale_y_continuous(limits = c(0,1.5), breaks = seq(0,1.5,0.5), expand = c(0,0), position = "left", labels = scales::label_parse()) +
      # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(title = "", x = "Working-Age Workplace Contacts (#/person/day)", y = "Relative Frequency") +
      guides(fill = guide_legend(title = NULL)) +
      theme(legend.position = c(0.202, 0.995), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

p7 <- ggplot(nns_data, aes(x = NNs, y = value, group = country, color = igroup)) +
      facet_wrap(~ "Adult Population by Sector") +
      geom_line(linewidth = 0.5, alpha = 0.50) +
      scale_color_manual(values = c("LLMIC" = "orange", "UMIC" = "turquoise", "HIC" = "azure4")) +  
      theme_bw() +
      scale_x_continuous(breaks = 1:46, expand = c(0,0), position = "bottom",
                         labels = c("1: Agri, Hunting, Forestry", "2: Fishing & Aquaculture", "3: Mining & Energy Products",
                                    "4: Mining & Non-Energy Products", "5: Mining Support", "6: Food, Beverages & Tobacco",
                                    "7: Textiles, Leather & Footwear", "8: Wood & Cork Products", "9: Paper Products & Printing",
                                    "10: Coke & Petroleum Products", "11: Chemical Products", "12: Pharma & Medicinal Products",
                                    "13: Rubber & Plastics Products", "14: Non-Metallic Minerals", "15: Basic Metals",
                                    "16: Fabricated Metal Products", "17: Computer & Electronics", "18: Electrical Equipment",
                                    "19: Machinery & Equipment NEC", "20: Motor Vehicles & Trailers", "21: Transport Equipment",
                                    "22: Manufacturing & Machinery", "23: Electricity & Gas Supply", "24: Water & Waste Management",
                                    "25: Construction", "26: Wholesale & Retail", "27: Land Transport & Pipelines", "28: Water Transport",
                                    "29: Air Transport", "30: Warehousing & Support", "31: Postal & Courier Activities",
                                    "32: Accommodation & Food", "33: Publishing & Broadcasting", "34: Telecommunications",
                                    "35: IT & Info Services", "36: Financial & Insurance", "37: Real Estate", "38: Professional & Technical",
                                    "39: Admin & Support Services", "40: Public Admin & Social Sec", "41: Education", "42: Health & Social Work",
                                    "43: Arts & Recreation", "44: Other Services", "45: Households as Employers", "Not Working")) + 
      scale_y_log10(limits = c(0.0001,1), expand = c(0,0), position = "left", labels = scales::label_number(scale = 100, suffix = "")) + 
      theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1)) + 
      labs(title = "", x = "Economic Sector", y = "Percentage of Adult Population (%)") +
      guides(color = guide_legend(title = NULL, ncol = 3)) +
      theme(legend.position = c(0.263, 0.995), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

p8 <- ggplot(gvaw_data, aes(x = obj, y = value, group = country, color = igroup)) +
      facet_wrap(~ "GVA per Worker by Sector") +
      geom_line(linewidth = 0.5, alpha = 0.50) +
      scale_color_manual(values = c("LLMIC" = "orange", "UMIC" = "turquoise", "HIC" = "azure4")) +  
      theme_bw() +
      scale_x_continuous(breaks = 1:45, expand = c(0,0), position = "bottom",
                         labels = c("1: Agri, Hunting, Forestry", "2: Fishing & Aquaculture", "3: Mining & Energy Products",
                                    "4: Mining & Non-Energy Products", "5: Mining Support", "6: Food, Beverages & Tobacco",
                                    "7: Textiles, Leather & Footwear", "8: Wood & Cork Products", "9: Paper Products & Printing",
                                    "10: Coke & Petroleum Products", "11: Chemical Products", "12: Pharma & Medicinal Products",
                                    "13: Rubber & Plastics Products", "14: Non-Metallic Minerals", "15: Basic Metals",
                                    "16: Fabricated Metal Products", "17: Computer & Electronics", "18: Electrical Equipment",
                                    "19: Machinery & Equipment NEC", "20: Motor Vehicles & Trailers", "21: Transport Equipment",
                                    "22: Manufacturing & Machinery", "23: Electricity & Gas Supply", "24: Water & Waste Management",
                                    "25: Construction", "26: Wholesale & Retail", "27: Land Transport & Pipelines", "28: Water Transport",
                                    "29: Air Transport", "30: Warehousing & Support", "31: Postal & Courier Activities",
                                    "32: Accommodation & Food", "33: Publishing & Broadcasting", "34: Telecommunications",
                                    "35: IT & Info Services", "36: Financial & Insurance", "37: Real Estate", "38: Professional & Technical",
                                    "39: Admin & Support Services", "40: Public Admin & Social Sec", "41: Education", "42: Health & Social Work",
                                    "43: Arts & Recreation", "44: Other Services", "45: Households as Employers")) + 
      scale_y_log10(limits = c(10,10^7), expand = c(0,0), position = "left", labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
      theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1)) + 
      labs(title = "", x = "Economic Sector", y = "Annual GVA per Worker ($)") +
      guides(color = guide_legend(title = NULL, ncol = 3)) +
      theme(legend.position = c(0.263, 0.217), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

p9 <- ggplot(wfh_data, aes(x = wfh, y = value, group = country, color = igroup)) +
      facet_wrap(~ "Home-Working by Sector") +
      geom_line(linewidth = 0.5, alpha = 0.50) +
      scale_color_manual(values = c("LLMIC" = "orange", "UMIC" = "turquoise", "HIC" = "azure4")) +  
      theme_bw() +
      scale_x_continuous(breaks = 1:45, expand = c(0,0), position = "bottom",
                         labels = c("1: Agri, Hunting, Forestry", "2: Fishing & Aquaculture", "3: Mining & Energy Products",
                                    "4: Mining & Non-Energy Products", "5: Mining Support", "6: Food, Beverages & Tobacco",
                                    "7: Textiles, Leather & Footwear", "8: Wood & Cork Products", "9: Paper Products & Printing",
                                    "10: Coke & Petroleum Products", "11: Chemical Products", "12: Pharma & Medicinal Products",
                                    "13: Rubber & Plastics Products", "14: Non-Metallic Minerals", "15: Basic Metals",
                                    "16: Fabricated Metal Products", "17: Computer & Electronics", "18: Electrical Equipment",
                                    "19: Machinery & Equipment NEC", "20: Motor Vehicles & Trailers", "21: Transport Equipment",
                                    "22: Manufacturing & Machinery", "23: Electricity & Gas Supply", "24: Water & Waste Management",
                                    "25: Construction", "26: Wholesale & Retail", "27: Land Transport & Pipelines", "28: Water Transport",
                                    "29: Air Transport", "30: Warehousing & Support", "31: Postal & Courier Activities",
                                    "32: Accommodation & Food", "33: Publishing & Broadcasting", "34: Telecommunications",
                                    "35: IT & Info Services", "36: Financial & Insurance", "37: Real Estate", "38: Professional & Technical",
                                    "39: Admin & Support Services", "40: Public Admin & Social Sec", "41: Education", "42: Health & Social Work",
                                    "43: Arts & Recreation", "44: Other Services", "45: Households as Employers")) + 
      #scale_y_continuous(limits = c(0,0.50), expand = c(0,0), position = "left", labels = scales::label_number(scale = 100, suffix = "")) + 
      scale_y_log10(limits = c(0.0001,1), expand = c(0,0), position = "left", labels = scales::label_number(scale = 100, suffix = "")) + 
      theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1)) + 
      labs(title = "", x = "Economic Sector", y = "Percentage of Workers Capable of Home-Working (%)") +
      guides(color = guide_legend(title = NULL, ncol = 3)) +
      theme(legend.position = c(0.263, 0.995), legend.justification = c(1, 1), legend.box.just = "right", 
            legend.key.size = unit(0.80, "cm"), legend.text = element_text(size = 8))

patchwork <- p1 / p2 / (p3 + p4) / (p5 + p6)
ggsave("figure_S1a.png", plot = patchwork, height = 14, width = 10)
patchwork <- p7 / p8 / p9
ggsave("figure_S1b.png", plot = patchwork, height = 14, width = 10)