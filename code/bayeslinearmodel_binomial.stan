data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0,upper=1> y[n];
  matrix[n,p] x;
  matrix[n,n] dmat;
}
parameters {
  real alpha;
  vector[p] beta;
}
model {
  for(i in 1:n){
    y[i] ~ binomial_logit(1,alpha+x[i,]*beta);
  }
  //y ~ multi_normal(mu+w,diag_matrix(rep_vector(tausq,n)));
  beta ~ multi_normal(rep_vector(0,p),10*diag_matrix(rep_vector(1,p)));
  //sigmasq ~ uniform(0,20);
  //phi ~ uniform(3/300,3/10);
}
generated quantities {
  vector[n] log_lik;

  for (i in 1:n) {
    log_lik[i] = binomial_logit_lpmf(y[i] | 1, alpha+x[i,]*beta);
  }
}
