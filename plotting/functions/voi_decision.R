voi_decision <- function(df, p_list) {
  
  df <- df %>% pivot_wider(names_from = strategy, values_from = SLpc) %>% 
               rename(t1 = "No Closures", t2 = "School Closures", t3 = "Economic Closures", t4 = "Elimination") %>%
               mutate(across(c(t1, t2, t3, t4), ~ - .)) %>% #net benefit is negative loss
               group_by(location, disease) %>%
               mutate(across(all_of(unlist(p_list)), ~ if (all(.x > 0)) {predict(BoxCoxTrans(.x), .x)} else {.x})) %>% 
               summarise(parameter = sapply(p_list, function(x) paste(x, collapse = ",")),
                         res       = evppi(outputs = cur_data() %>% dplyr::select(t1, t2, t3, t4) %>% as.data.frame(),
                                           inputs  = pick(unlist(p_list)),
                                           pars    = p_list)$evppi / 
                                     pmin(mean(-t1), mean(-t2), mean(-t3), mean(-t4))) %>%
               arrange(-res) %>% mutate(rank = letters[row_number(-res)]) %>%
               ungroup()
  return(df)
  
}