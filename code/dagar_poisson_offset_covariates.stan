data {
  int <lower=1> N;
  vector[N] N_nei; //number of directed neighbors for each node
  int <lower=1> N_edges;  // total numbere of edges
  int <lower=1, upper=N-1> nei[N_edges]; //vector stacking up the direected neighbors of each node
  int<lower=0, upper=N_edges> adjacency_ends[N]; //Where adjacency of each node ends in the nei vector

  int y[N];//observed counts
  real offset[N];//expected counts
  real X[N]; // single covariate (can be easily extended to multiple)
}

parameters {
  real <lower=0, upper=1> rho; // Spatial correlation parameter
  real <lower=0> std_dev_w; // Precision of the spatial effects
  //real <lower=0> tau; // Precision of the spatial effects
  vector[N] w; // Spatial Random Effect
  real alpha; // intercept
  real beta;

}

model {
  real logoff[N];
  vector[N] b; //
  vector[N] vec_var; //
  vector[N] t_rowsum; // only the rowsum of t is used
  vector[N] std_dev; // Rescaled std_dev by std_dev_w


  logoff=log(offset);
  // Construct w
  vec_var = (1 - rho * rho) ./ (1 + (N_nei - rep_vector(1,N)) * rho * rho);
  b = rho ./ (1 + (N_nei - rep_vector(1,N)) * rho * rho );
  // Linear in number of edges
  t_rowsum[1] = 0;
  for(i in 2:N){
    t_rowsum[i] = sum(w[nei[(adjacency_ends[i-1]+1):adjacency_ends[i]]]) * b[i];
  }

  std_dev = std_dev_w * sqrt(vec_var);

  for(i in 1:N){
    w[i] ~ normal(t_rowsum[i],std_dev[i]);
    y[i] ~ poisson_log(logoff[i] + alpha + X[i]*beta + w[i]);
  }

  //beta0 ~ normal(0, 100);
  rho ~ beta(1,1);

}

generated quantities {
  vector[N] log_lik;
  real logoff[N];

  logoff=log(offset);

  for (i in 1:N) {
    log_lik[i] = poisson_log_lpmf(y[i] | logoff[i] + alpha + X[i]*beta + w[i]);
  }
}
