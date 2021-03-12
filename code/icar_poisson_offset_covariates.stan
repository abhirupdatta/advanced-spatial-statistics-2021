### ICAR model for POisson data {
### part of the code taken from https://mc-stan.org/workshops/dec2017/spatial_smoothing_icar.html#15
}
data {
  int<lower=0> N;
    int<lower=0> N_edges;
    int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
    int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

    int<lower=0> y[N];
    real offset[N];//expected counts
    real X[N]; // single covariate (can be easily extended to multiple)

}

parameters {
  real <lower=0> std_dev_w; // std. dev. of the spatial effects
  vector[N-1] w_raw; // Spatial Random Effect
  real alpha; // intercept
  real beta;
}

transformed parameters {
  vector[N] w;
  w[1:(N-1)]=w_raw;
  w[N]=-sum(w_raw);
}

model {
  real logoff[N];

  logoff=log(offset);
  for(i in 1:N) {
    y[i] ~ poisson_log(logoff[i] + alpha + X[i]*beta + w[i]*std_dev_w);
  }
  target += -0.5 * dot_self(w[node1] - w[node2]);
}

generated quantities {
  vector[N] log_lik;
  real logoff[N];

  logoff=log(offset);

  for (i in 1:N) {
    log_lik[i] = poisson_log_lpmf(y[i] | logoff[i] + alpha + X[i]*beta + w[i]*std_dev_w);
  }
}
