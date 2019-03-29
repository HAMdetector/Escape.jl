data {
  int<lower=1> n_entries;
  int<lower=1> n_alleles;
  int<lower=0, upper=1> y[n_entries];
  matrix[n_entries, n_alleles] hla_matrix;
  vector[n_entries] phylogeny_effect;
}

parameters {
  real intercept;
  real<lower=0> phylogeny_coefficient;
  vector[n_alleles] beta_hla;
}

model {
  intercept ~ student_t(3, 0, 5);
  beta_hla ~ student_t(7, 0, 3);
  phylogeny_coefficient ~ student_t(7, 0, 3);
  
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