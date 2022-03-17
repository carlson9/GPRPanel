#' @description
#' The function to estimate the variance-covariance
#' returns the fitted model.
#'
#' @param form An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted
#' @param group Column name in data containing names of units under study (e.g., countries)
#' @param time Column name in data containing the temporal element (e.g., year)
#' @param data A data.frame containing everything in form and the group and time columns
#' @param loglik Logical indicating whether or not to estimate the log-likelihood. Defaults to FALSE
#' @param seed The seed for random number generation of the initial Stan parameter values. The default is 1919
#' @param n.cores The number of cores to use when executing the Markov chains in parallel. The default is 'max' (i.e., all available cores)
#' @param iter The number of iterations per chain. Default is 1000
#' @param chains A positive integer specifying the number of Markov chains. The default is 4
#' @param refresh Controls how often the progress of the sampling is reported
#' @param ... Additional arguments that are passed to sampling in rstan
#' 
#' @details
#' We need details
#' 
#' 
#' @return An rstan fit object
#' @author David Carlson \email{carlson.david@@wustl.edu} and Dilara Karaduman \email{karadumandilara3@@gmail.com}
#' @examples
#' 
#' \donttest{
#' we need an example
#' }
#'
#'
#' @seealso \code{\link{GPRPanelPred}}
#' @rdname GPRPanel
#' @aliases GPRPanel,ANY-method
#' @export
setGeneric(name="GPRPanel",
           def=function(form, group, time, data,
                        loglik = FALSE, seed = 1919,
                        n.cores = 'max',
                        iter = 1000, chains = 4, refresh = 100, ...)
           {standardGeneric("GPRPanel")}
)
#' @export
setMethod(f="GPRPanel",
          definition=function(form, group, time, data,
                              loglik = FALSE, seed = 1919,
                              n.cores = 'max',
                              iter = 1000, chains = 4, refresh = 100, ...){
    names = as.character(attr(terms(form), 'variables'))[-1]
    data = na.omit(data[, c(names, group, time)])
    g = model.matrix(~ as.factor(data[, group]) - 1)
    y = data[,names[1]]
    X = data[,names[2:length(names)]]
    X = apply(X, 2, scale)
    X = cbind(as.matrix(X), g)
    if(n.cores == 'max') rstan::options(mc.cores = parallel::detectCores())
    else rstan::options(mc.cores = n.cores)
    if(auto_write) rstan::rstan_options(auto_write = TRUE)
    t = scale(data[,time])
    X_corr = as.matrix(cbind(X, t))
    N=dim(X)[1]
    M=dim(X_corr)[2]
    K=dim(X)[2]
    if(loglik) fit = rstan::sampling(stanmodels$gp-fit-multivariate-loglik,
            data = list(X=X, N=N, K=K, y=y, M=M, X_corr=X_corr),
            seed = seed, iter=iter, chains=chains, refresh = refresh, ...)
    else fit = rstan::sampling(stanmodels$gp-fit,
            data = list(X=X, N=N, K=K, y=y, M=M, X_corr=X_corr),
            seed = seed, iter=iter, chains=chains, refresh = refresh, ...)
    return(fit)
}
