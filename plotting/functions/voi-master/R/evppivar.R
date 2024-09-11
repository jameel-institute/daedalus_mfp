##' Calculate the expected value of partial perfect information for an estimation problem
##'
##' Calculate the expected value of partial perfect information for an estimation problem.  This computes the expected reduction in variance in some quantity of interest with perfect information about a parameter or parameters of interest.   
##'
##' @param outputs a vector of values for the quantity of interest, sampled from the uncertainty distribution of this quantity that is induced by the uncertainty about the parameters.   This can also be a data frame with one column. 
##'
##' Typically this will come from a Monte Carlo sample, where we first sample from the uncertainty distributions of the parameters, and then compute the quantity of interest as a function of the parameters.  It might also be produced by a Markov Chain Monte Carlo sample from the joint distribution of parameters and outputs. 
##'
##' @return A data frame with a column \code{pars}, indicating the parameter(s), and a column \code{evppi}, giving the corresponding EVPPI.
##' 
##' @inheritParams evppi
##' 
##' @references
##' Jackson, C., Presanis, A., Conti, S., & De Angelis, D. (2019). Value of information:
##' Sensitivity analysis and research design in Bayesian evidence synthesis.
##' Journal of the American Statistical Association, 114(528), 1436-1449.
##'
##' Jackson, C., Johnson, R., de Nazelle, A., Goel, R., de Sa, T. H.,
##' Tainio, M., & Woodcock, J. (2021). A guide to value of information
##' methods for prioritising research in health impact
##' modelling. Epidemiologic Methods, 10(1).
##'
##' Jackson, C. H., Baio, G., Heath, A., Strong, M., Welton, N. J., &
##' Wilson, E. C. (2022). Value of Information analysis in models to
##' inform health policy. Annual Review of Statistics and its
##' Application, 9, 95-118.
##' 
##' @export
evppivar <- function(outputs,
                     inputs,
                     pars=NULL,
                     method=NULL,
                     nsim=NULL,
                     verbose=TRUE,
                     return_fit=FALSE,
                     ...)
{
    inputs <- check_inputs(inputs, iname=deparse(substitute(inputs)))
    check_outputs_vector(outputs, inputs)
    if (is.list(pars)) {
        return(evppivar_list(outputs=outputs, inputs=inputs, pars=pars, 
                             method=method, nsim=nsim, verbose=verbose, return_fit=return_fit, ...))
    }
    pars <- check_pars(pars, inputs)
    opts <- list(...)
    if (is.null(method))
        method <- default_evppi_method(pars)

    if (is.null(nsim)) nsim <- nrow(inputs)
    outputs <- outputs[1:nsim]
    inputs <- inputs[1:nsim,,drop=FALSE]
    
    if (method %in% npreg_methods) {
        rese <- evppivar_npreg(outputs=outputs, inputs=inputs, 
                    pars=pars, method=method, verbose=verbose, ...)
    } else stop("Other methods not implemented yet")
    
    if (return_fit) {
      model  <- attr(rese[[2]], "model")
      x_var  <- inputs[[pars[1]]]
      x_seq  <- seq(min(x_var), max(x_var), length.out = 100)
      #x_seq  <- sort(unique(x_var)) 
      
      if (length(pars) == 1) {
        grid   <- expand.grid(x = x_seq)
        colnames(grid) <- pars
        pred   <- predict(model, newdata = grid, se.fit = TRUE)
        grid$y <- pred$fit
        grid$l <- pred$fit - pred$se.fit*qnorm(1-0.05/2)
        grid$u <- pred$fit + pred$se.fit*qnorm(1-0.05/2)
        #grid$l <- pred$fit - qt(1-0.05/2, df = model$df.residual)*sqrt(pred$se.fit^2 + sigma(model)^2)
        #grid$u <- pred$fit + qt(1-0.05/2, df = model$df.residual)*sqrt(pred$se.fit^2 + sigma(model)^2)
        colnames(grid) <- c("x", "y", "l", "u")} 
      else if (length(pars) == 2) {
        y_var  <- inputs[[pars[2]]]
        y_seq  <- seq(min(y_var), max(y_var), length.out = 100)
        #y_seq  <- sort(unique(y_var)) 
        grid   <- expand.grid(x = x_seq, y = y_seq)
        colnames(grid) <- pars  
        grid$z <- predict(model, newdata = grid)
        colnames(grid) <- c("x", "y", "z")} 
      else {
        stop("pars too long!")}
      
      res    <- grid}
    else {
      res    <- cbind(pars = paste(pars, collapse=","),
                      rese[[1]])}
    
    res
}

evppivar_list <- function(outputs, inputs, pars, method, nsim, verbose, return_fit, ...)
{
    npars <- length(pars)
    eres <- vector(npars, mode="list") 
    for (i in seq_len(npars)){
        eres[[i]] <- evppivar(outputs=outputs, inputs=inputs, pars=pars[[i]], 
                              method=method, nsim=nsim, verbose=verbose, return_fit=return_fit, ...)
    }
    do.call("rbind", eres)
}

check_outputs_vector <- function(outputs, inputs){
  if (is.data.frame(outputs)) {
    if (ncol(outputs) > 1)
      stop("if `outputs` is supplied as a data frame, it should have only one column")
    outputs <- unlist(outputs)
  }
  if (!is.numeric(outputs))
    stop("`outputs` should be a numeric vector or a data frame with one column")
  if (length(outputs) != nrow(inputs))
    stop(sprintf("Length of `outputs` (%s) should equal the number of rows of `inputs` (%s)",
                 length(outputs), nrow(inputs)))
}

evppivar_npreg <- function(outputs, inputs, pars, method, verbose, ...){
    fitted <- fitted_npreg_fn(method)(y=outputs, inputs=inputs, pars=pars, 
                                      verbose=verbose, ...)
    list(data.frame(evppi = var(fitted)), fitted)
} 