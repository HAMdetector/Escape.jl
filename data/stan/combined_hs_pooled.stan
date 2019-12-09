functions {
    vector ll(vector beta, vector theta, real[] x, int[] y) {
        int N_obs = y[1];
        int N_alleles = num_elements(theta) - 1;
        matrix[N_obs, N_alleles] X;
        int y_obs[N_obs] = y[2:(N_obs + 1)];
        real lp;

        for (i in 1:N_obs) {
            X[i] = to_vector(x[(1 + (i - 1) * N_alleles) : (i * N_alleles)])';
        }
        lp = bernoulli_logit_glm_lpmf(y_obs | X, theta[1], theta[2:N_alleles + 1]);

        return [lp]';
    }
}

data {
    int N_obs;
    int N_alleles;
    int N_shards;
    int ys[N_shards, N_obs + 1];
    real xs[N_shards, N_obs * N_alleles];
}

transformed data {
    real slab_scale = 0.5;
    real slab_scale2 = square(slab_scale);
    real slab_df = 25;
    real half_slab_df = 0.5 * slab_df;
}

parameters {
    real m0;
    vector[N_shards] intercepts;
    vector<lower=0>[N_alleles] lambda[N_shards];
    vector[N_alleles] beta_tilde[N_shards];
    vector<lower=0>[N_shards] c2_tilde;
    vector<lower=0>[N_shards] tau_tilde;
}

transformed parameters {
    vector[N_alleles] beta_hla[N_shards];
    vector[N_alleles] lambda_tilde[N_shards];
    vector[N_alleles + 1] theta[N_shards];
    vector[0] beta;
    
    {
        for (i in 1:N_shards) {
            real tau0 = (m0 / (N_alleles - m0)) * (2 / sqrt(1.0 * N_obs));
            real tau = tau0 * tau_tilde[i];
            real c2 = slab_scale2 * c2_tilde[i];

            lambda_tilde[i] =
                sqrt(c2 * square(lambda[i]) ./ (c2 + square(tau) * square(lambda[i])));

            beta_hla[i] = tau * lambda_tilde[i] .* beta_tilde[i];

            theta[i][1] = intercepts[i];
            theta[i][2:N_alleles + 1] = beta_hla[i];
        }
    }
}


model {
    m0 ~ normal(0, 0.25);
    intercepts ~ normal(0, 5);
    for (i in 1:N_shards) {
        lambda[i] ~ cauchy(0, 1);
        tau_tilde[i] ~ cauchy(0, 1);
        beta_tilde[i] ~ normal(0, 1);
        c2_tilde[i] ~ inv_gamma(half_slab_df, (slab_df * slab_scale2) / 2);
    }

    target += map_rect(ll, beta, theta, xs, ys);
}