data {
  int<lower=0> n;
  vector[n] y;
  real x[n];
  matrix[n,n] dmat;
  real id[n,n];
}
parameters {
  real alpha;
  real beta;
  real<lower=0> tausq;
  real<lower=0> sigmasq;
  real<lower=0> phi;
}
transformed parameters {
  matrix[n,n] R;
  matrix[n,n] Sigma;
  vector[n] mu;

  for(i in 1:n){
    mu[i] = alpha+x[i]*beta;
    for(j in 1:n){
      R[i,j]=sigmasq*exp(-phi*dmat[i,j]);
    }
  }
  Sigma = R + diag_matrix(rep_vector(tausq,n));
}
model {
  y ~ multi_normal(mu,Sigma);
  beta ~ normal(0,10);
  sigmasq ~ uniform(0,10);
  tausq ~ uniform(0,10);
  phi ~ uniform(2,30);
}
