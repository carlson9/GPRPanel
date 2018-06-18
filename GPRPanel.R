GPRPanel = function(form, group, time, data,
        loglik = FALSE, seed = 1919,
        cores = 'max', n.cores = NULL, auto_write = TRUE,
        iter = 1000, chains = 4, refresh = 100){
    g = model.matrix(~ group - 1, data)
    names = as.character(attr(terms(form), 'variables'))[-1]
    y = data[names[1]]
    X = data[names[2:length(names)]]
    X = cbind(as.matrix(X), g)
    requireNamespace(rstan)
    if(cores == 'max') options(mc.cores = parellel::detectCores())
    else options(mc.cores = n.cores)
    if(auto_write) rstan_options(auto_write = TRUE)
    t = scale(data[,time])
    X = apply(X, 2, scale)
    X_corr = as.matrix(cbind(X, t))
    N=dim(X)[1]
    M=dim(X_corr)[2]
    K=dim(X)[2]
    if(loglik) fit = stan(file = 'gp-fit-multivariate-loglik.stan',
            data = list(X=X, N=N, K=K, y=y, M=M, X_corr=X_corr),
            seed = seed, iter=iter, chains=chains, refresh = refresh)
    else fit = stan(file = 'gp-fit.stan',
            data = list(X=X, N=N, K=K, y=y, M=M, X_corr=X_corr),
            seed = seed, iter=iter, chains=chains, refresh = refresh)
    return(fit)
}
