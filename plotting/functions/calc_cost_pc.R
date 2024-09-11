calc_cost_pc <- function(inp_df) {
  
  out_df <- inp_df %>% filter(!is.na(SEC))
  gva    <- out_df %>% select(starts_with("obj")) %>% rowSums(na.rm = TRUE) * 365
  out_df <- out_df %>% mutate(SECpc  = 100*SEC/gva,
                              VLYLpc = 100*VLYL/gva,
                              GDPLpc = 100*GDPL/gva,
                              VSYLpc = 100*VSYL/gva)
  return(out_df)
  
}