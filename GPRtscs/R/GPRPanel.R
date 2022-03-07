#' @description
#' The function to estimate the variance-covariance
#' returns the fitted model.
#'
#' @param form an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted
#' @param group A character vector/column of names of parameters.
#' @param time the time to evaluate the ODE system
#' @param data the data: a pdata.frame object or an ordinary data.frame
#' @param loglik the value of the log-likelihood. Defaults to FALSE.
#' @param seed The seed for random number generation. The default is generated from 1 to the maximum integer supported by R on the machine.
#' @param n.cores The number of cores to use when executing the Markov chains in parallel.
#' @param auto_write auto_write = TRUE, then a serialized version of the compiled model will be automatically saved to the hard disk in the same directory as the .stan file or in R’s temporary directory if the Stan program is expressed as a character string.
#' @param iter The number of iterations of IWLS default method ("glm.fit" uses iteratively reweighted least squares) used.
#' @param chains A positive integer specifying the number of Markov chains. The default is 4.
#' @param refresh to control how often the progress of the sampling is reported
#' @param ... Additional arguments that are passed to sampling.
#' @details
#'
#' @return as.array, as.matrix, and as.data.frame return an array, matrix, and data.frame, respectively.
#' @author David Carlson \email{carlson.david@@wustl.edu} and Dilara Karaduman \email{karadumandilara3@@gmail.com}
#' @examples
#'
#'
#' @seealso \code{\link{}} \code{\link{}} \code{\link{runMod}} \code{\link{}}
#' @rdname
#' @aliases
#' @export


GPRPanel = function(form, group, time, data,
        loglik = FALSE, seed = 1919,
        n.cores = 'max', auto_write = TRUE,
        iter = 1000, chains = 4, refresh = 100){
    g = model.matrix(~ as.factor(data[, group]) - 1)
    names = as.character(attr(terms(form), 'variables'))[-1]
    y = data[,names[1]]
    X = data[,names[2:length(names)]]
    X = apply(X, 2, scale)
    X = cbind(as.matrix(X), g)
    requireNamespace('rstan')
    if(n.cores == 'max') options(mc.cores = parallel::detectCores())
    else options(mc.cores = n.cores)
    if(auto_write) rstan_options(auto_write = TRUE)
    t = scale(data[,time])
    X_corr = as.matrix(cbind(X, t))
    N=dim(X)[1]
    M=dim(X_corr)[2]
    K=dim(X)[2]
    if(loglik) fit = rstan::stan(file = '../GPRtscs/inst/stan/gp-fit-multivariate-loglik.stan',
            data = list(X=X, N=N, K=K, y=y, M=M, X_corr=X_corr),
            seed = seed, iter=iter, chains=chains, refresh = refresh)
    else fit = rstan::stan(file = 'gp-fit.stan',
            data = list(X=X, N=N, K=K, y=y, M=M, X_corr=X_corr),
            seed = seed, iter=iter, chains=chains, refresh = refresh)
    return(fit)
}
?

