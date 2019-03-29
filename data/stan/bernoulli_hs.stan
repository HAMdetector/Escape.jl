data {
  int<lower=1> n_entries;
  int<lower=1> n_alleles;
  int<lower=0, upper=1> y[n_entries];
  matrix[n_entries, n_alleles] hla_matrix;
}

transformed data {
  real m0 = 1;
  real slab_scale = 0.5;
  real slab_scale2 = square(slab_scale);
  real slab_df = 25;
  real half_slab_df = 0.5 * slab_df;
}

parameters {
  real intercept;
  vector<lower=0>[n_alleles] lambda;
  vector[n_alleles] beta_tilde;
  real<lower=0> c2_tilde;
  real<lower=0> tau_tilde; 
}

transformed parameters {
  vector[n_alleles] beta_hla;
  {
    real tau0 = (m0 / (n_alleles - m0)) * (2 / sqrt(1.0 * n_entries)); 
    real tau = tau0 * tau_tilde;

    real c2 = slab_scale2 * c2_tilde;

    vector[n_alleles] lambda_tilde =
      sqrt( c2 * square(lambda) ./ (c2 + square(tau) * square(lambda)) );

    beta_hla = tau * lambda_tilde .* beta_tilde;
  }
}

model {
  intercept ~ normal(0, 2);

  lambda ~ cauchy(0, 1);
  tau_tilde ~ cauchy(0, 1);
  beta_tilde ~ normal(0, 1);
  c2_tilde ~ inv_gamma(half_slab_df, (slab_df * slab_scale2) / 2);

  y ~ bernoulli_logit(intercept + hla_matrix * beta_hla);
}

generated quantities {
  real y_rep[n_entries];
  real theta;
  vector[n_entries] log_lik;

  for (i in 1:n_entries) {
    log_lik[i] = bernoulli_logit_lpmf(y[i] | intercept + hla_matrix[i,] * beta_hla);
    y_rep[i] = bernoulli_rng(inv_logit(intercept + hla_matrix[i,] * beta_hla));
  }
  
  theta = sum(y_rep) / n_entries;
}