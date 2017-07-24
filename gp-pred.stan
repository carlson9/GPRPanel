data {
int<lower=1> N;
int<lower=1> zN;
int<lower=1> K;
int<lower=1> M;
//int<lower=1> zM;
matrix[N+zN,K] XZ;
matrix[N+zN,M] XZ_corr;
//matrix[zN,K] Z;
//matrix[zN,zM] Z_corr;
vector[N] y;
matrix[K,K] B;
vector[K] nu;
}
parameters {
real<lower=0> nug;
real<lower=0> sig_sq;
vector<lower=0>[M] d1;
vector<lower=0>[M] d2;
vector[K] b;
vector[zN] z;
}
model {
//matrix[N+zN,N+zN] Sigma;
matrix[N+zN,N+zN] Sigma;
vector[N+zN] mu;
matrix[N+zN,K] Mu;
vector[M] d;

vector[N+zN] yz;


//matrix[zN,zN] zSigma;
//vector[zN] zmu;
//matrix[zN,K] zMu;

for(m in 1:M){
d1[m] ~ gamma(1,20);
d2[m] ~ gamma(10,10);
d[m] = .5*(d1[m] + d2[m]);
}
for (i in 1:(N+zN-1)) {
for (j in (i+1):(N+zN)) {
vector[M] summand;
for(m in 1:M){
summand[m] = -pow(XZ_corr[i,m] - XZ_corr[j,m],2)/d[m];
}
Sigma[i,j] = exp(sum(summand));
Sigma[j,i] = Sigma[i,j];
}
}
for (i in 1:(N+zN)){
for(k in 1:K){
Mu[i,k] = XZ[i,k]*b[k];
}
mu[i]=sum(Mu[i,1:K]);
}
for (i in 1:(N+zN))
Sigma[i,i] = 1 + nug; // + jitter



//for (i in 1:(zN-1)) {
//  for (j in (i+1):zN) {
//   vector[zM] zsummand;
//  for(m in 1:zM){
//      zsummand[m] = -pow(Z_corr[i,m] - Z_corr[j,m],2)/d[m];
//    }
//    zSigma[i,j] = exp(sum(zsummand));
//    zSigma[j,i] = zSigma[i,j];
//  }
//}
//for (i in 1:zN){
//  for(k in 1:K){
//    zMu[i,k] = Z[i,k]*b[k];
//  }
//  zmu[i]=sum(zMu[i,1:K]);
//}
//for (i in 1:zN)
//  zSigma[i,i] = 1 + nug; // + jitter




sig_sq ~ inv_gamma(1,1);
nug ~ exponential(1);

b ~ normal(0,10);

for(n in 1:N) yz[n] = y[n];
for(n in 1:zN) yz[N+n] = z[n];

yz ~ multi_normal(mu,sig_sq*Sigma);

//z ~ multi_normal(zmu,sig_sq*zSigma);


}
