data {
  int<lower=0> n;
  vector[n] y;
  real x[n];

  int<lower=0> n_new;
  real x_new[n_new];
}
parameters {
  real alpha;
  real beta;
  real<lower=0> tausq;
  real<lower=0> phi;
}
model {
  for(i in 1:n){
    y[i] ~ normal(alpha+x[i]*beta,sqrt(tausq));
  }
  #y ~ multi_normal(mu+w,diag_matrix(rep_vector(tausq,n)));
  beta ~ normal(0,10);
  tausq ~ uniform(0,10);
  phi ~ uniform(2,30);
}
generated quantities {
  real y_new[n_new];
  vector[n] log_lik;

  {

  for(i in 1:n_new){
      y_new[i] = normal_rng(alpha+x_new[i]*beta,sqrt(tausq));
  }

  for (i in 1:n) {
    log_lik[i] = normal_lpdf(y[i] | alpha+x[i]*beta,sqrt(tausq));
  }
}
}
