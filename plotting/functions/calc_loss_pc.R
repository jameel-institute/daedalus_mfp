calc_loss_pc <- function(df_inps, df_outs) {
  
  df_gdp <- df_inps %>% mutate(pa_pop = rowSums(across(Npop_1:Npop_1)),
                               sa_pop = rowSums(across(Npop_2:Npop_4)),
                               wa_pop = rowSums(across(Npop_5:Npop_13)),
                               ra_pop = rowSums(across(Npop_14:Npop_21)),
                               wf     = rowSums(across(NNs_1:NNs_45)),
                               gdp    = 365*rowSums(across(starts_with("obj"))))
  df     <- df_outs %>% left_join(df_gdp %>% dplyr::select(location, country, pa_pop, sa_pop, wa_pop, ra_pop, wf, gdp), 
                                  by = c("location", "country")) %>%  
                        mutate(across(starts_with("fsiz") | starts_with("chad") | starts_with("mort"), ~ . * 10^5/
                                      case_when(grepl("_1$", cur_column()) ~ pa_pop,
                                                grepl("_2$", cur_column()) ~ sa_pop,
                                                grepl("_3$", cur_column()) ~ wa_pop,
                                                grepl("_4$", cur_column()) ~ ra_pop))) %>%
                        mutate(hpeak = 10^5 * hpeak /(pa_pop + sa_pop + wa_pop + ra_pop)) %>%
                        mutate(across(starts_with("vlyl") | starts_with("gdpl") | starts_with("vsyl"), ~ . * 100/gdp)) %>%
                        mutate(VLYLpc = rowSums(across(starts_with("vlyl"))),
                               GDPLpc = rowSums(across(starts_with("gdpl"))),
                               VSYLpc = vsyl) %>%
                        mutate(SLpc = VLYLpc + GDPLpc + VSYLpc) %>%
                        mutate(wdunem = wdunem / wf)
  return(df)
  
}