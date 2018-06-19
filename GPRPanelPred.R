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
    Z[,(1:(dim(Z)[2] - dim(g)[2]))] = apply(Z[,(1:(dim(Z)[2] - dim(g)[2]))], 2, scale)
    Z_corr = cbind(Z, scale(Z_corr[, dim(Z_corr)[2]))
    zN = dim(Z)[2]
    fit = rstan::stan(file = 'gp-pred.stan',
            data = list(XZ=rbind(X,Z), N=N, zN=zN, K=K, y=y, M=M, XZ_corr=rbind(X_corr,Z_corr)),
            seed = seed, iter=iter, chains=chains, refresh = refresh)
    return(fit)
    g = model.matrix(~ data[, group] - 1)
    g = model.matrix(~ data[, group] - 1)
    