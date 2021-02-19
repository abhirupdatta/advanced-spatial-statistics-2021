data {
  int<lower=0> n;
  vector[n] y;
  real x[n];
  matrix[n,n] dmat;

  int<lower=0> n_new;
  real x_new[n_new];
  matrix[n_new,n] dmat_new;
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
generated quantities {
  real w_new[n_new];
  real y_new[n_new];
  vector[n] log_lik;

  {
  matrix[n,n] R_inv;
  matrix[n_new,n] R_cross;
  matrix[n_new,n] krig_weights;
  vector[n_new] krig_mean;
  vector[n_new] krig_var;

  R_inv=inverse(R);
  for(i in 1:n_new){
        for(j in 1:n){
          R_cross[i,j]=sigmasq*exp(-phi*dmat_new[i,j]);
          }
        }
  krig_weights=R_cross*R_inv;
  krig_mean=krig_weights*w;
  krig_var=rep_vector(sigmasq,n_new);

  for(i in 1:n_new){
       krig_var[i]=krig_var[i] - krig_weights[i,]*(R_cross[i,])';
    }

  for(i in 1:n_new){
      w_new[i] = normal_rng(krig_mean[i],sqrt(krig_var[i]));
      y_new[i] = normal_rng(alpha+x_new[i]*beta+w_new[i],sqrt(tausq));
  }

  for (i in 1:n) {
    log_lik[i] = normal_lpdf(y[i] | alpha+x[i]*beta+w[i],sqrt(tausq));
  }
}
}
