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

parameters {
    vector[N_alleles] b_alleles[N_shards];
    vector[N_shards] intercepts;
    real mu;
    real<lower=0> sigma;
}

transformed parameters {
    vector[N_alleles + 1] theta[N_shards];
    vector[2] beta;

    for (i in 1:N_shards) {
        theta[i][1] = intercepts[i];
        theta[i][2:N_alleles + 1] = b_alleles[i];

        beta[1] = mu;
        beta[2] = sigma;
    }
}

model {
    mu ~ normal(0, 1);
    sigma ~ normal(0, 3);

    for (i in 1:N_shards) {
        b_alleles[i] ~ normal(mu, sigma);
    }

    target += sum(map_rect(ll, beta, theta, xs, ys));
}