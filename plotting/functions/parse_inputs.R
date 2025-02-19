parse_inputs <- function(df_inps) {
  
  df <- df_inps %>% 
        rowwise() %>%
        mutate(#demography
               mean_age = sum(across(Npop_1:Npop_21) * seq(2.5, 102.5, by = 5)) / sum(across(Npop_1:Npop_21)),#assuming band uniformity
               #mixing
               pa_prop  = sum(across(Npop_1:Npop_1))   / sum(across(Npop_1:Npop_21)),
               sa_prop  = sum(across(Npop_2:Npop_4))   / sum(across(Npop_1:Npop_21)),
               wa_prop  = sum(across(Npop_5:Npop_13))  / sum(across(Npop_1:Npop_21)),
               ra_prop  = sum(across(Npop_14:Npop_21)) / sum(across(Npop_1:Npop_21)),
               matAL    = list(matrix(unlist(c_across(matAL_1:matAL_16)), nrow = 4, byrow = FALSE)),
               AL_wavg  = sum(rowSums(matAL) * c(pa_prop, sa_prop, wa_prop, ra_prop)),#depends on age but unavoidable
               matAHT   = list(matrix(unlist(c_across(matAHT_1:matAHT_16)), nrow = 4, byrow = FALSE)),
               AHT_wavg = sum(rowSums(matAHT) * c(pa_prop, sa_prop, wa_prop, ra_prop)),#depends on age but unavoidable
               matAS    = list(matrix(unlist(c_across(matAS_1:matAS_16)), nrow = 4, byrow = FALSE)),
               AS_wavg  = sum(matAS[2,]),#only contacts of school-age population 
               #economy
               epop     = sum(across(NNs_1:NNs_45)) / sum(across(Npop_5:Npop_13)),
               #disease
               ifr_1918 = sum(0.669*c(0.02284,0.00398,0.00478,0.00983,0.01700,0.02922,0.02470,0.02205,0.01647,0.01195,
                                      0.01647,0.01169,0.03081,0.04144,0.04941,0.04941,0.04941,0.04941,0.04941,0.04941,0.04941)
                              * across(Npop_1:Npop_21)) / sum(across(Npop_1:Npop_21)),
               ifr_dlta = sum(1.85*c(0.000016,0.000016,0.000070,0.000070,0.000309,0.000309,0.000844,0.000844,0.001610,0.001610,
                                     0.005950,0.005950,0.019300,0.019300,0.042800,0.042800,0.078000,0.078000,0.078000,0.078000,0.078000)
                              * across(Npop_1:Npop_21)) / sum(across(Npop_1:Npop_21)),
               ifr_sars = sum(0.867*c(0.017,0.017,0.017,0.017,0.024,0.024,0.024,0.024,0.089,0.089,
                                      0.089,0.089,0.255,0.255,0.255,0.255,0.177,0.177,0.177,0.177,0.177)
                              * across(Npop_1:Npop_21)) / sum(across(Npop_1:Npop_21)),
               #valuation
               na       = list(c(c_across(starts_with("Npop_"))[1:17], sum(across(starts_with("Npop_"))[18:21]))),
               la       = list(c_across(starts_with("la_"))),
               gdp      = 365*sum(across(starts_with("obj_"))),
               gdppc    = gdp/sum(na),
               vsl      = case_when(location == "LLMIC" || (location == "UMIC" && gdppc < 0.008809) ~ 
                                      10.9 * ((0.008809 / 0.060362)^0.85) * (gdppc / 0.008809),
                                    (location == "UMIC" && gdppc > 0.008809) || location == "HIC" ~ 
                                      10.9 * ((gdppc / 0.060362)^0.85),
                                    TRUE ~ 
                                      NA),
               vly      = vsl/(sum(la * na)/sum(na)),
               vly_gdpc = vly/gdppc,
               nstud    = sum(across(starts_with("Npop_"))[2:4]),
               llpc     = case_when(location == "LLMIC" ~ sum(c(0.62, 0.22) * c(27, 55)) / sum(c(27, 55)) / 0.33,
                                    location == "UMIC"  ~ 0.22 / 0.33,
                                    location == "HIC"   ~ 0.09 / 0.33),   
               vsy      = llpc*gdp/nstud,
               vsy_gdpc = vsy/gdppc) %>%
        group_by(location) %>%
        mutate(#demography
               pr_le    = {la_cols    <- pick(starts_with("la_"))
                           la_cols[[names(which.max(colSums(abs(cor(la_cols)))))]]},#mean le would not be independent of age distribution
               #economy
               pr_workf = {workf_cols <- across(starts_with("NNs_")) / rowSums(across(NNs_1:NNs_45))
                           workf_cols[[names(which.max(colSums(abs(cor(workf_cols)))))]]},
               pr_gvapw = {gvapw_cols <- across(starts_with("obj_")) / across(starts_with("NNs"))
                           pred_power <- cor(gvapw_cols, use = "pairwise.complete.obs")
                           noNA_power <- pred_power[,!is.na(colSums(gvapw_cols))] #don't consider predictors with any NA values
                           365*gvapw_cols[[names(which.max(colSums(noNA_power)))]]},
               pr_wfh   = {wfh_cols   <- pick(starts_with("wfh_")) %>% dplyr::select(., -wfh_41) 
                           wfh_cols[[names(which.max(colSums(abs(cor(wfh_cols)))))]]},
               #disease
               ifr_1918 = mean(ifr_1918),
               ifr_dlta = mean(ifr_dlta),
               ifr_sars = mean(ifr_sars),
               #valuation
               mean_vly = mean(vly_gdpc),
               mean_vsy = mean(vsy_gdpc)) %>% 
        ungroup()
  return(df)
  
}