data {
  int<lower=0> N;
    int<lower=0> N_edges;
    int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
    int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

    int<lower=0> y[N];
}

parameters {
  //real beta0; //the intercept
  real <lower=0> std_dev_w; // Precision of the spatial effects
  //real <lower=0> tau; // Precision of the spatial effects
  vector[N] w; // Spatial Random Effect

}

model {
  y ~ poisson_log(w * std_dev_w);
  target += -0.5 * dot_self(w[node1] - w[node2]);
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = poisson_log_lpmf(y[i] | std_dev_w*w[i]);
  }
}
