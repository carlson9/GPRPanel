#' @description
#' Fit a model defined in the Stan model and return the fitted
#' returns the fitted model.
#'
#' @param form An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted
#' @param group Column name in data containing names of units under study (e.g., countries)
#' @param time Column name in data containing the temporal element (e.g., year)
#' @param data A data.frame containing everything in form and the group and time columns
#' @param Z Z value for normal distribution
#' @param Z_corr z-scores for inference on parameters
#' @param loglik Logical indicating whether or not to estimate the log-likelihood. Defaults to FALSE
#' @param seed The seed for random number generation of the initial Stan parameter values. The default is 1919
#' @param n.cores The number of cores to use when executing the Markov chains in parallel. The default is 'max' (i.e., all available cores)
#' @param iter The number of iterations per chain. Default is 1000
#' @param chains A positive integer specifying the number of Markov chains. The default is 4
#' @param refresh Controls how often the progress of the sampling is reported
#' @param auto_write A logical scalar (defaulting to TRUE) that controls whether a compiled instance of a stanmodel-class is written to the hard disk in the same directory.
#' @param ... Additional arguments that are passed to sampling in rstan
#'
#' @details
#'GPRPAnelPred function sets run-time options used to fit the model.
#'
GPRPanelPred = function(form, group, time, data, Z, Z_corr,
        seed = 1919,
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

    #as of now, Z and Z_corr need to be supplied correctly by user - need to make more user friendly

    zN = dim(Z)[1]
    fit = rstan::stan(file = '../GPRtscs/inst/stan/gp-pred.stan',
            data = list(XZ=rbind(X,Z), N=N, zN=zN, K=K, y=y, M=M, XZ_corr=rbind(X_corr,Z_corr)),
            seed = seed, iter=iter, chains=chains, refresh = refresh)
    return(fit)
}



