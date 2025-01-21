table_formatting <- function(dataframe, column) {
  
  tlabels <- as.character(unique(dataframe[[column]]))
  
  dataframe[[column]] <- factor(dataframe[[column]], levels = tlabels)
  
  df_split <- split(dataframe, dataframe[[column]])
  
  dataframe <- do.call(rbind, lapply(df_split, function(group) {
    # Create a data frame with NA values and the same number of columns as 'group'
    na_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(group)))
    colnames(na_row) <- colnames(group)  # Assign the column names
    rbind(group, na_row)  # Bind the group with the NA row
  }))
  
  na_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(dataframe)))
  colnames(na_row) <- colnames(dataframe)
  
  dataframe <- rbind(na_row,dataframe)#add NAs for top row header
  dataframe <- dataframe[-nrow(dataframe),]#remove NAs at bottom
  dataframe[[column]] <- NULL#remove column
  
  inds                 <- which(!complete.cases(dataframe))
  inds2                <- which(c(FALSE,dataframe[1:nrow(dataframe)-1,1]==dataframe[2:nrow(dataframe),1]))
  dataframe[inds,1]    <- tlabels
  dataframe[[1]][inds] <- paste0("\\bfseries{", dataframe[[1]][inds], "}")
  dataframe[inds2,1]   <- NA
  dataframe            <- dataframe %>% mutate_all(~ ifelse(is.na(.), "", .))
  
}