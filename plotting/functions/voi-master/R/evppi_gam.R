fitted_gam <- function(y, inputs, pars, verbose=FALSE, ...){
    
    #new
    is_negative <- all(y <= 0)
    if (is_negative) 
      {y <- -y} 
    else {y}
    
    opts <- list(...)
    gam_formula <- opts$gam_formula
    pars <- clean_pars(pars)
    colnames(inputs) <- clean_pars(colnames(inputs))
    
    # if (is.null(gam_formula))
    #     gam_formula <- default_gam_formula(pars)
    # gam_formula <- formula(sprintf("y ~ %s", gam_formula))
    # model <- mgcv::gam(gam_formula,  family = Gamma(link = "log"), data = inputs)
    
    # if (pars %in% c("mean_age", "pr_le", "epop", "pr_workf", "pr_gvapw", "AL_wavg", "AHT_wavg", "AS_wavg", "workp", 
    #                 "Tres", "sdc", "t_tit", "t_vax")) {
    #   mdir <- 'mpi'
    # } else if (pars %in% c("pr_wfh", "sda", "sdb", "trate", "Hmax", "arate", "puptake")) {
    #   mdir <- 'mpd'
    # }
    # scam_formula <- formula(sprintf("y ~ s(%s, bs = '%s')", paste(pars, collapse=", "), mdir))
    # model        <- scam(scam_formula, family = Gamma(link = "log"), data = inputs)
    
    scam_formula <- formula(sprintf("y ~ s(%s, bs = 'mpi')", paste(pars, collapse=", ")))
    model1       <- scam(scam_formula, family = Gamma(link = "log"), data = inputs)
    scam_formula <- formula(sprintf("y ~ s(%s, bs = 'mpd')", paste(pars, collapse=", ")))
    model2       <- scam(scam_formula, family = Gamma(link = "log"), data = inputs)
    models       <- list(model1, model2)
    model        <- models[[which.min(c(AIC(model1), AIC(model2)))]]
    
    res <- (model$fitted) * ifelse(is_negative, -1, 1)
    attr(res, "model") <- model
    res
}

default_gam_formula <- function(pars){
    karg <- if (length(pars) >=4) ", k=4" else ""
    sprintf("te(%s, bs='cr'%s)", paste(pars, collapse=", "), karg)
}

fitted_rep_gam <- function(model, B) { 
    beta_rep <- mvtnorm::rmvnorm(B, coef(model), vcov(model))
    fitted_rep <- beta_rep %*% t(predict(model,type="lpmatrix"))
    fitted_rep
}

check_plot_gam <- function(mod){
  oldpar <- graphics::par(no.readonly=TRUE)
  on.exit(par(oldpar))
  graphics::par(mfrow=c(2,2))
  mgcv::gam.check(mod)
}

check_stats_gam <- function(mod){
    list(AIC = stats::AIC(mod)) 
}

#' Generate a string with all interactions of a certain degree, to be used in a GAM formula
#'
#' @param x Character vector of variable names
#'
#' @param degree Maximum interaction degree
#'
#' @return A string looking like the right hand side of a GAM formula with tensor product interactions.
#'
#' For example, if `x` is `c("x1","x2","x3")`, then `all_interactions(x, degree=2)` should return
#'
#' `"te(x1,x2) + te(x1,x3) + te(x1,x3)"`
#'
#' @examples
#' x <- c("x1","x2","x3")
#' all_interactions(x, 2)
#'
#' @export
#' 
all_interactions <- function(x, degree=2){
  c_mat <- utils::combn(x, degree)
  c_comma <- apply(c_mat, 2, function(y) paste(y, collapse = ","))
  c_tevec <- paste0("te(",c_comma,")")
  form_str <- paste(c_tevec, collapse = " + ")
  form_str
}
