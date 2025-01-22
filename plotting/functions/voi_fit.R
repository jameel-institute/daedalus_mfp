voi_fit <- function(df) {
  
  df <- df %>% group_by(location, disease, strategy) %>%
               mutate(xaxis  = get(unique(parameter)),
                      lambda = if (all(xaxis > 0)) {BoxCoxTrans(xaxis)$lambda} else {NA},
                      xaxis  = if (all(xaxis > 0)) {predict(BoxCoxTrans(xaxis), xaxis)} else {xaxis}) %>%
               summarise(parameter = unique(parameter),
                         lambda    = unique(lambda),
                         res       = evppivar(outputs    = SLpc,
                                              inputs     = pick(xaxis),
                                              pars       = "xaxis",
                                              return_fit = TRUE)) %>%
               unnest_wider(res) %>%
               mutate(x = if (!is.na(unique(lambda))) {bc_inv(x, unique(lambda))} else {x}) %>%
               ungroup() 
  return(df)
  
}