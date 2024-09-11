calc_cost_bdown <- function(inp_df) {
  
  out_df <- inp_df %>% group_by(location, disease, strategy) %>%
                       mutate(pLYL = mean(VLYL/SEC),
                              pGDP = mean(GDPL/SEC),
                              pSYL = mean(VSYL/SEC)) %>% 
                       ungroup()
  return(out_df)
  
}