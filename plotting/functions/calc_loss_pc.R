calc_loss_pc <- function(df_inps, df_outs) {
  
  df_gdp <- df_inps %>% mutate(gdp = 365*rowSums(across(starts_with("obj"))))
  df     <- df_outs %>% left_join(df_gdp %>% dplyr::select(location, country, gdp), by = c("location", "country")) %>%  
                        mutate(across(starts_with("vlyl") | starts_with("gdpl") | starts_with("vsyl"), ~ . * 100/gdp)) %>%
                        mutate(VLYLpc = rowSums(across(starts_with("vlyl"))),
                               GDPLpc = rowSums(across(starts_with("gdpl"))),
                               VSYLpc = vsyl) %>%
                        mutate(SLpc = VLYLpc + GDPLpc + VSYLpc)
  return(df)
  
}