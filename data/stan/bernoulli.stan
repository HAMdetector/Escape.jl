data {
  int<lower=1> n_entries;
  int<lower=1> n_alleles;
  int<lower=0, upper=1> y[n_entries];
  matrix[n_entries, n_alleles] hla_matrix;
}

parameters {
  real intercept;
  vector[n_alleles] beta_hla;
}

model {
  intercept ~ student_t(3, 0, 5);
  beta_hla ~ student_t(7, 0, 3);

  y ~ bernoulli_logit(intercept + hla_matrix * beta_hla);
}

generated quantities {
  real y_rep[n_entries];
  real theta;
  
  for (i in 1:n_entries) {
    y_rep[i] = bernoulli_rng(inv_logit(intercept + hla_matrix[i,] * beta_hla));
  }
  
  theta = sum(y_rep)/n_entries;
}