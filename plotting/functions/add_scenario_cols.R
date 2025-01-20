add_scenario_cols <- function(path) {
  
  title <- sub("\\.csv$", "", basename(path))
  parts <- unlist(strsplit(title, "_"))
  if (length(parts) < 3) {
    parts[2] <- NA_character_
    parts[3] <- NA_character_
  }
  df    <- read.csv(path) %>% mutate(location = parts[1], disease = parts[2], strategy = parts[3])
  return(df)
  
}