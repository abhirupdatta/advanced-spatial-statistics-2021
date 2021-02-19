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
  vector[n] w;
  real<lower=0> tausq;
  real<lower=0> sigmasq;
  real<lower=0> phi;
}
transformed parameters {
  matrix[n,n] R;

  for(i in 1:n){
    for(j in 1:n){
      R[i,j]=sigmasq*exp(-phi*dmat[i,j]);
    }
  }
}
model {
  w ~ multi_normal(rep_vector(0,n),R);
  for(i in 1:n){
    y[i] ~ normal(alpha+x[i]*beta+w[i],sqrt(tausq));
  }
  #y ~ multi_normal(mu+w,diag_matrix(rep_vector(tausq,n)));
  beta ~ normal(0,10);
  sigmasq ~ uniform(0,10);
  tausq ~ uniform(0,10);
  phi ~ uniform(2,30);
}
