order_scenario_cols <- function(df) {
  
  df <- df %>% mutate(location = factor(location, levels = c("LLMIC","UMIC","HIC")),
                      disease  = gsub(" ", "-", disease) %>% paste0("-X"),
                      disease  = factor(disease, levels = c("Influenza-2009-X","Influenza-1957-X","Influenza-1918-X",
                                                            "Covid-Omicron-X","Covid-Delta-X","Covid-Wildtype-X","SARS-X")),
                      strategy = factor(strategy, levels = c("No Closures","School Closures",
                                                              "Economic Closures","Elimination")))
  return(df)
  
}