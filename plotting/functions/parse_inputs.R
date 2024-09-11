parse_inputs <- function(inp_df) {
  
  calculate_popw_av <- function(df, matrix) {
    
    Npop_mat <- as.matrix(df %>% select(starts_with("Npop")))
    Npop_sum <- rowSums(Npop_mat)# == 50m
    popw_av  <- rowSums(Npop_mat * matrix) / Npop_sum
    return(popw_av)
    
  }
  
  calculate_IGHC <- function(Npop_list, CM_list) {
    
    Npop <- matrix(unlist(Npop_list), nrow = 1, ncol = 21, byrow = TRUE)
    Npop <- c(Npop[1:15], sum(Npop[16:21]))
    CM   <- matrix(unlist(CM_list), nrow = 16, ncol = 16, byrow = FALSE)
    Ncar <- c(sum(Npop[1:4]), sum(Npop[5:13]), sum(Npop[14:16]))
    Ccar <- cbind(rowSums(CM[,1:4]), rowSums(CM[,5:13]), rowSums(CM[,14:16]))
    Ccar <- rbind((Npop[1:4]   %*% Ccar[1:4,])   / Ncar[1],
                  (Npop[5:13]  %*% Ccar[5:13,])  / Ncar[2],
                  (Npop[14:16] %*% Ccar[14:16,]) / Ncar[3])
    IGHC <- (Ncar %*% (rowSums(Ccar) - diag(Ccar))) / (Ncar %*% rowSums(Ccar))
    return(IGHC)
    
  }
  
  out_df <- inp_df %>% mutate(#demography
                              mean_age = {
                                age_midp_mat <- matrix(seq(2, 102, by = 5), nrow = nrow(.), ncol = 21, byrow = TRUE)
                                calculate_popw_av(., age_midp_mat)},
                              sd_age  = {
                                age_midp_mat <- matrix(seq(2, 102, by = 5), nrow = nrow(.), ncol = 21, byrow = TRUE)
                                mean_age_mat <- matrix(mean_age, nrow = nrow(.), ncol = 21, byrow = FALSE)
                                squared_dev  <- (age_midp_mat - mean_age_mat)^2
                                sqrt(calculate_popw_av(.,squared_dev))},
                              skew_age = {
                                age_midp_mat <- matrix(seq(2, 102, by = 5), nrow = nrow(.), ncol = 21, byrow = TRUE)
                                mean_age_mat <- matrix(mean_age, nrow = nrow(.), ncol = 21, byrow = FALSE)
                                cubed_dev    <- (age_midp_mat - mean_age_mat)^3
                                calculate_popw_av(.,cubed_dev) / sd_age^3},
                              le       = {
                                la_columns   <- select(., starts_with("la"))
                                la_columns[[names(which.max(colSums(cor(la_columns))))]]},#best predictor is 50-54 year-olds
                              #mixing
                              IGHC     = map2_dbl(select(., starts_with("Npop")) %>% split(1:nrow(.)), 
                                                  select(., starts_with("CM")) %>% split(1:nrow(.)), calculate_IGHC),
                              #economy
                              EPOP     = rowSums(select(., starts_with("NNs"))) / rowSums(select(., Npop5:Npop13)),
                              TES      = rowSums(select(., NNs26:NNs45)) / rowSums(select(., starts_with("NNs"))),
                              gvapw    = {
                                gvapw_cols   <- across(starts_with("obj")) / across(starts_with("NNs"))
                                pred_power   <- cor(gvapw_cols, use = "complete.obs")
                                noNA_power   <- pred_power[,!is.na(colSums(gvapw_cols))]
                                365*gvapw_cols[[names(which.max(colSums(noNA_power)))]]},#best predictor is retail sector (annual)
                              wfh      = {
                                wfh_columns  <- select(., starts_with("wfhu"))
                                wfh_columns[[names(which.max(colSums(cor(wfh_columns))))]]}) %>% #best predictor is IT-related sectors
                       select(-igroup,-country,-starts_with("Npop"),-starts_with("la"),-starts_with("CM"),-schoolA1,-travelA3,
                              -starts_with("NNs"),-starts_with("obj"),-gnipc,-starts_with("wfhl"),-starts_with("wfhu"),-vly)
  return(out_df)
  
}

# bin_quantile <- function(counts, quantile) {
#   
#   total_pop   <- sum(counts)
#   cumul_pop   <- cumsum(counts)
#   quant_ind   <- which(cumul_pop >= quantile*total_pop)[1]
#   cumul_pop_l <- ifelse(quant_ind == 1, 0, cumul_pop[quant_ind - 1])
#   bin_start   <- seq(0, 100, by = 5)[quant_ind]
#   bin_posit   <- (quantile*total_pop - cumul_pop_l)/counts[quant_ind]
#   bin_width   <- 5
#   quant_est   <- bin_start + bin_posit * bin_width
#   return(quant_est)
#   
# }
#
# calculate_FNR <- function(CM_columns) {
#   
#   CM  <- matrix(unlist(CM_columns), nrow = 16, ncol = 16, byrow = FALSE)
#   FNR <- norm(CM - diag(diag(CM)), type = "F") / norm(CM, type = "F")
#   return(FNR)
#   
# }
# 
# CDR      = rowSums(select(., Npop1:Npop3)) / rowSums(select(., Npop4:Npop13)),
# ODR      = rowSums(select(., Npop14:Npop21)) / rowSums(select(., Npop4:Npop13)),
# q1_age   = apply(select(., starts_with("Npop")), 1, bin_quantile, 0.25),
# q2_age   = apply(select(., starts_with("Npop")), 1, bin_quantile, 0.50),
# q3_age   = apply(select(., starts_with("Npop")), 1, bin_quantile, 0.75),
# mean_la  = {
#   la_midp_mat  <- as.matrix(select(., starts_with("la")))
#   la_midp_mat  <- cbind(la_midp_mat, replicate(3,la_midp_mat[,18]))
#   weighted_av(., la_midp_mat)},
# var_la   = {
#   la_midp_mat  <- as.matrix(select(., starts_with("la")))
#   la_midp_mat  <- cbind(la_midp_mat, replicate(3,la_midp_mat[,18]))
#   mean_la_mat  <- matrix(mean_la, nrow = nrow(.), ncol = 21, byrow = FALSE)
#   squared_dev  <- (la_midp_mat-mean_la_mat)^2
#   weighted_av(.,squared_dev)},
# skew_la  = {
#   la_midp_mat  <- as.matrix(select(., starts_with("la")))
#   la_midp_mat  <- cbind(la_midp_mat, replicate(3,la_midp_mat[,18]))
#   mean_la_mat  <- matrix(mean_la, nrow = nrow(.), ncol = 21, byrow = FALSE)
#   cubed_dev    <- (la_midp_mat-mean_la_mat)^3
#   weighted_av(.,cubed_dev) / var_la^(3/2)},
# FNR      = map_dbl(select(., starts_with("CM")) %>% split(1:nrow(.)), calculate_FNR),
# PES      = rowSums(select(., NNs1:NNs5)) / rowSums(select(., starts_with("NNs"))),
# SES      = rowSums(select(., NNs6:NNs25)) / rowSums(select(., starts_with("NNs"))),
# wfhu     = rowSums(across(starts_with("wfhu")) * across(starts_with("NNs"))) / 
#            rowSums(across(starts_with("NNs")))