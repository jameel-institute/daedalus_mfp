add_scenario_cols <- function(inp_path) {
  
  title  <- sub("\\.csv$", "", basename(inp_path))
  parts  <- unlist(strsplit(title, "_"))
  out_df <- read.csv(inp_path) %>%
            mutate(location = parts[1], disease = parts[2], strategy = parts[3])
  return(out_df)
  
}