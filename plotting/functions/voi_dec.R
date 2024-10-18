voi_dec <- function(inp_df, p_list) {
  
  out_df <- inp_df %>% select(-SEC, -VLYL, -GDPL, -VSYL, -VLYLpc, -GDPLpc, -VSYLpc) %>%
                       pivot_wider(names_from = strategy, values_from = SECpc) %>% 
                       rename(t1 = "No Closures", t2 = "School Closures", t3 = "Economic Closures", t4 = "Elimination") %>%
                       drop_na(t1, t2, t3, t4) %>%
                       mutate(across(c(t1, t2, t3, t4), ~ - .)) %>%
                       mutate(trate = log10(trate),
                              sdb   = log10(sdb),
                              Hmax  = log10(Hmax),
                              arate = log10(arate)) %>% 
                       group_by(location, disease) %>%
                       summarise(parameter = sapply(p_list, function(x) paste(x, collapse = ",")),
                                 res       = evppi(outputs    = cur_data() %>% select(t1, t2, t3, t4) %>% as.data.frame(),
                                                   inputs     = pick(setdiff(names(inp_df), c("location", "disease", "strategy", 
                                                                                              "SEC", "VLYL", "GDPL", "VSYL",
                                                                                              "SECpc", "VLYLpc", "GDPLpc", "VSYLpc"))),
                                                   pars       = p_list,
                                                   verbose    = FALSE)$evppi / pmin(mean(-t1), mean(-t2), mean(-t3), mean(-t4))) %>%
                       arrange(res) %>% mutate(rank = letters[row_number(res)]) %>%
                       ungroup()
  return(out_df)
  
}