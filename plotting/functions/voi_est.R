voi_est <- function(inp_df, p_list) {
  
  out_df <- inp_df %>% mutate(trate = log10(trate),
                              sdb   = log10(sdb),
                              Hmax  = log10(Hmax),
                              arate = log10(arate)) %>% 
                       group_by(location, disease, strategy) %>%
                       summarise(parameter = sapply(p_list, function(x) paste(x, collapse = ",")),
                                 res       = evppivar(outputs    = SECpc,
                                                      inputs     = pick(setdiff(names(inp_df), c("location", "disease", "strategy"))),
                                                      pars       = p_list,
                                                      verbose    = FALSE,
                                                      return_fit = FALSE)$evppi / var(SECpc)) %>%
                       arrange(res) %>% mutate(rank = letters[row_number(res)]) %>%
                       ungroup()
  return(out_df)
  
}