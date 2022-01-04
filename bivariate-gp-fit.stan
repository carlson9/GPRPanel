data {
int<lower=1> N;
int<lower=1> K1;
int<lower=1> K2;
int<lower=1> M1;
int<lower=1> M2;
matrix[N,K1] X1;
matrix[N,K2] X2;
matrix[N,M1] X_corr1;
matrix[N,M2] X_corr2;
vector[2] y[N];
}
parameters {
real<lower=0> nug1;
real<lower=0> sig_sq1;
real<lower=0> nug2;
real<lower=0> sig_sq2;
vector<lower=0>[M1] d11;
vector<lower=0>[M1] d21;
vector<lower=0>[M2] d12;
vector<lower=0>[M2] d22;
vector[K1] b1;
vector[K2] b2;
cholesky_factor_corr[2] L_Omega;
vector<lower=0>[2] L_sigma;
vector[N] ystar1;
vector[N] ystar2;
}
model {
matrix[N,N] Sigma1;
vector[N] mu1;
matrix[N,K1] Mu1;
vector[M1] d1;
matrix[N,N] Sigma2;
vector[N] mu2;
matrix[N,K2] Mu2;
vector[M2] d2;
row_vector[2] mu[N];
matrix[2,2] L_Sigma;

for(m in 1:M1){
d11[m] ~ gamma(1,20);
d21[m] ~ gamma(10,10);
d1[m] = .5*(d11[m] + d21[m]);
}

for(m in 1:M2){
d12[m] ~ gamma(1,20);
d22[m] ~ gamma(10,10);
d2[m] = .5*(d12[m] + d22[m]);
}



for (i in 1:(N-1)) {
for (j in (i+1):N) {
vector[M1] summand1;
for(m in 1:M1){
summand1[m] = -pow(X_corr1[i,m] - X_corr1[j,m],2)/d1[m];
}
Sigma1[i,j] = exp(sum(summand1));
Sigma1[j,i] = Sigma1[i,j];
}
}

for (i in 1:(N-1)) {
for (j in (i+1):N) {
vector[M2] summand2;
for(m in 1:M2){
summand2[m] = -pow(X_corr2[i,m] - X_corr2[j,m],2)/d2[m];
}
Sigma2[i,j] = exp(sum(summand2));
Sigma2[j,i] = Sigma2[i,j];
}
}


for (i in 1:N){
for(k in 1:K1){
Mu1[i,k] = X1[i,k]*b1[k];
}
mu1[i]=sum(Mu1[i,1:K1]);
}
for (i in 1:N)
Sigma1[i,i] = 1 + nug1; // + jitter

sig_sq1 ~ inv_gamma(1,1);
nug1 ~ exponential(1);

b1 ~ normal(0,10);
ystar1 ~ multi_normal(mu1,sig_sq1*Sigma1);


for (i in 1:N){
for(k in 1:K2){
Mu2[i,k] = X2[i,k]*b2[k];
}
mu2[i]=sum(Mu2[i,1:K2]);
}
for (i in 1:N)
Sigma2[i,i] = 1 + nug2; // + jitter

sig_sq2 ~ inv_gamma(1,1);
nug2 ~ exponential(1);

b2 ~ normal(0,10);
ystar2 ~ multi_normal(mu2,sig_sq2*Sigma2);

for(n in 1:N)
  mu[n] = [ystar1[n], ystar2[n]];

L_Sigma = diag_post_multiply(L_Omega, L_sigma);
  
L_Omega ~ lkj_corr_cholesky(4);
L_sigma ~ cauchy(0, 2.5);
  
y ~ multi_normal_cholesky(mu, L_Sigma);


}
