voi_fit <- function(df, target_var) {
  
  df <- df %>% group_by(location, disease, strategy) %>%
               mutate(xaxis  = xaxis + 1e-10,
                      lambda = if (all(xaxis > 0)) {BoxCox.lambda(xaxis, lower = 0, upper = 1)} else {NA},
                      xaxis  = if (all(xaxis > 0)) {BoxCox(xaxis, unique(lambda))} else {xaxis}) %>%
               summarise(parameter = unique(parameter),
                         lambda    = unique(lambda),
                         res       = evppivar(outputs    = {{target_var}},
                                              inputs     = cur_data() %>% dplyr::select(xaxis) %>% rename_with(~parameter, xaxis),
                                              pars       = parameter,
                                              return_fit = TRUE)) %>%
               unnest_wider(res) %>%
               mutate(xval = if (!is.na(unique(lambda))) {InvBoxCox(xval, unique(lambda))} else {xval}) %>%
               ungroup() 
  return(df)
  
}