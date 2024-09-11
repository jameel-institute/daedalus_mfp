find_best_strats <- function(inp_df) {
  
  out_df <- inp_df %>% group_by(location, disease, strategy) %>%
            mutate(med_SECpc  = quantile(SECpc, 0.50),
                   mean_SECpc = mean(SECpc),
                   q3_SECpc   = quantile(SECpc, 0.75),
                   max_SECpc  = max(SECpc)) %>%
            group_by(location, disease) %>%
            mutate(min_med  = ifelse(med_SECpc  == min(med_SECpc), TRUE, FALSE),
                   min_mean = ifelse(mean_SECpc == min(mean_SECpc), TRUE, FALSE),
                   min_q3   = ifelse(q3_SECpc   == min(q3_SECpc), TRUE, FALSE),
                   min_max  = ifelse(max_SECpc  == min(max_SECpc), TRUE, FALSE)) %>%
            ungroup()
  return(out_df)
  
}