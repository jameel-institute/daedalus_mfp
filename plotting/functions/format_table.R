format_table <- function(df, column) {
  
  tlabels      <- as.character(unique(df[[column]]))
  df[[column]] <- factor(df[[column]], levels = tlabels)
  df_split     <- split(df, df[[column]])
  
  df <- do.call(rbind, lapply(df_split, function(group) {
    na_row           <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(group)))
    colnames(na_row) <- colnames(group)  # Assign the column names
    rbind(group, na_row)  # Bind the group with the NA row
    }))
  
  na_row           <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(df)))
  colnames(na_row) <- colnames(df)
  
  df           <- rbind(na_row,df)#add NAs for top row header
  df           <- df[-nrow(df),]#remove NAs at bottom
  df[[column]] <- NULL#remove column
  
  inds          <- which(!complete.cases(df))
  inds2         <- which(c(FALSE,df[1:nrow(df)-1,1]==df[2:nrow(df),1]))
  df[inds,1]    <- tlabels
  df[[1]][inds] <- paste0("\\textit{", df[[1]][inds], "}")
  df[inds2,1]   <- NA
  df            <- df %>% mutate_all(~ ifelse(is.na(.), "", .))
  
}