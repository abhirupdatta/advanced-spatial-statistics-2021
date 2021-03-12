data {
  int<lower=0> N;
    int<lower=0> y[N];
    real offset[N];//expected counts
    real X[N]; // single covariate (can be easily extended to multiple)

}

parameters {
  real alpha; // intercept
  real beta;
}

model {
  real logoff[N];

  logoff=log(offset);
  for(i in 1:N) {
    y[i] ~ poisson_log(logoff[i] + alpha + X[i]*beta);
  }
}

generated quantities {
  vector[N] log_lik;
  real logoff[N];

  logoff=log(offset);

  for (i in 1:N) {
    log_lik[i] = poisson_log_lpmf(y[i] | logoff[i] + alpha + X[i]*beta);
  }
}
