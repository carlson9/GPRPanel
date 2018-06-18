GPRPanelPred = function(form, group, time, data, Z, Z_corr,
        seed = 1919,
        n.cores = 'max', auto_write = TRUE,
        iter = 1000, chains = 4, refresh = 100){
    g = model.matrix(~ group - 1, data)
    names = as.character(attr(terms(form), 'variables'))[-1]
    y = data[,names[1]]
    X = data[,names[2:length(names)]]
    X = cbind(as.matrix(X), g)
    requireNamespace(rstan)
    if(n.cores == 'max') options(mc.cores = parellel::detectCores())
    else options(mc.cores = n.cores)
    if(auto_write) rstan_options(auto_write = TRUE)
    t = scale(data[,time])
    X = apply(X, 2, scale)
    X_corr = as.matrix(cbind(X, t))
    N=dim(X)[1]
    M=dim(X_corr)[2]
    K=dim(X)[2]
    Z = apply(Z, 2, scale)
    Z_corr = apply(Z_corr, 2, scale)
    zN = dim(Z)[2]
    fit = stan(file = 'gp-pred.stan',
            data = list(XZ=rbind(X,Z), N=N, zN=zN, K=K, y=y, M=M, XZ_corr=rbind(X_corr,Z_corr)),
            seed = seed, iter=iter, chains=chains, refresh = refresh)
    return(fit)
}
