data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N,K] X;
  vector[N] y;
//  real<lower=0> rho;
//  matrix[K,K] V;
  matrix[K,K] B;
  vector[K] nu;
}
parameters {
  real<lower=0> nug;
  real<lower=0> sig_sq;
  real<lower=0> tau_sq;
  vector<lower=0>[K] d1;
  vector<lower=0>[K] d2;
  vector[K] b;
  vector[K] b0;
  cholesky_factor_corr[K] Lcorr;
  vector<lower=0>[K] sigma;
}
model {
  matrix[N,N] Sigma;
  vector[N] mu;
  matrix[N,K] Mu;
  vector[K] d;
//  W ~ inv_wishart(rho, V);

  for(k in 1:K){
    d1[k] ~ gamma(1,20);
    d2[k] ~ gamma(10,10);
    d[k] = .5*(d1[k] + d2[k]);
  }
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      vector[K] summand;
      for(k in 1:K){
        summand[k] = -pow(X[i,k] - X[j,k],2)/d[k];
      }
      Sigma[i,j] = exp(sum(summand));
      Sigma[j,i] = Sigma[i,j];
    }
  }
  for (i in 1:N){
    for(k in 1:K){
      Mu[i,k] = X[i,k]*b[k];
    }
    mu[i]=sum(Mu[i,1:K]);
  }
  for (i in 1:N)
    Sigma[i,i] = 1 + nug; // + jitter

  sig_sq ~ inv_gamma(1,1);
  tau_sq ~ inv_gamma(1,1);
  nug ~ exponential(1);

  b0 ~ multi_normal(nu,B);
  sigma ~ cauchy(0,5);
  Lcorr ~ lkj_corr_cholesky(1);
  b ~ multi_normal_cholesky(b0, sig_sq*tau_sq*diag_pre_multiply(sigma,Lcorr));
  y ~ multi_normal(mu,sig_sq*Sigma);
}

