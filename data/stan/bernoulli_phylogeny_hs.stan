data {
  int<lower=1> n_entries;
  int<lower=1> n_alleles;
  int<lower=0, upper=1> y[n_entries];
  matrix[n_entries, n_alleles] hla_matrix;
  vector[n_entries] phylogeny_effect;
}

transformed data {
  real m0 = 1;
  real slab_scale = 3;
  real slab_scale2 = square(slab_scale);
  real slab_df = 3;
  real half_slab_df = 0.5 * slab_df;
  real col_sums[n_alleles];
  for (i in 1:n_alleles) {
    col_sums[i] = sum(hla_matrix[,i]);
  }
}

parameters {
  real intercept;
  vector[n_alleles] beta_tilde;
  vector<lower=0>[n_alleles] lambda;
  real<lower=0> c2_tilde;
  real<lower=0> tau_tilde;
  real phylogeny_coefficient;
}

transformed parameters {
  vector[n_alleles] beta_hla;
  {
    real sigma =  2;
    real tau0[n_alleles];
    real tau[n_alleles];
    vector[n_alleles] lambda_tilde;
    real c2 = slab_scale2 * c2_tilde;
    
    for (i in 1:n_alleles) {
      tau0[i] = (m0 / (n_alleles - m0)) * (sigma / sqrt(1.0 * fmax(1, col_sums[i])));
    }
    
    // tau ~ cauchy(0, tau0)
    for (i in 1:n_alleles) {
      tau[i] = tau0[i] * tau_tilde;
    }

    // c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scale2)
    // Implies that marginally beta ~ student_t(slab_df, 0, slab_scale)

    for (i in 1:n_alleles)
      lambda_tilde[i] = sqrt(c2 * square(lambda[i]) / (c2 + square(tau[i]) * square(lambda[i])));

    // beta ~ normal(0, tau * lambda_tilde);
    for (i in 1:n_alleles)
      beta_hla[i] = tau[i] * lambda_tilde[i] * beta_tilde[i];
  }
}

model {
  beta_tilde ~ normal(0, 1);
  lambda ~ cauchy(0, 1);
  tau_tilde ~ cauchy(0, 1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
  phylogeny_coefficient ~ student_t(7, 0, 3);

  intercept ~ student_t(3, 0, 5);

  y ~ bernoulli_logit(intercept + logit(phylogeny_effect) * phylogeny_coefficient +
                      hla_matrix * beta_hla);
}

generated quantities {
  real y_rep[n_entries];
  real theta;
  vector[n_entries] log_lik;

  for (i in 1:n_entries) {
    log_lik[i] = bernoulli_logit_lpmf(y[i] | intercept + 
      logit(phylogeny_effect[i]) * phylogeny_coefficient + hla_matrix[i,] * beta_hla);
    y_rep[i] = bernoulli_rng(inv_logit(intercept + 
      logit(phylogeny_effect[i]) * phylogeny_coefficient + hla_matrix[i,] * beta_hla));
  }
  
  theta = sum(y_rep)/n_entries;
}