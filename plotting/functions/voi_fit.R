voi_fit <- function(inp_df) {
  
  out_df <- inp_df %>% mutate(trate = log10(trate),
                              sdb   = log10(sdb),
                              Hmax  = log10(Hmax),
                              arate = log10(arate)) %>% 
                       group_by(location, disease, strategy) %>%
                       summarise(parameter = unique(parameter),
                                 res       = evppivar(outputs    = SECpc,
                                                      inputs     = pick(setdiff(names(inp_df), c("location", "disease", "strategy"))),
                                                      pars       = unlist(strsplit(parameter, ",")),
                                                      verbose    = FALSE,
                                                      return_fit = TRUE)) %>%
                       unnest_wider(res) %>%                   
                       ungroup() %>%
                       mutate(parts = str_split(parameter, ",", simplify = TRUE),
                              p1    = parts[, 1],
                              p2    = if (ncol(parts) == 2) {parts[, 2]} else {NA_character_},
                              x     = if_else(p1 %in% c("trate", "sdb", "Hmax", "arate"), 10^x, x),
                              y     = if_else(p2 %in% c("trate", "sdb", "Hmax", "arate"), 10^y, y)) %>%
                       select(-parts, -p1, -p2)
  return(out_df)
  
}