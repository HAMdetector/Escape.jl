functions {
    vector ll(vector beta, vector theta, real[] x, int[] y) {
        int N_obs = y[1];
        int N_alleles = num_elements(theta) - 1;
        matrix[N_obs, N_alleles] X;
        vector[N_obs] alpha;
        int y_obs[N_obs] = y[2:(N_obs + 1)];
        real lp;

        for (i in 1:N_obs) {
            X[i] = to_vector(x[(2 + (i - 1) * (N_alleles + 1)) : (i * (N_alleles + 1))])';
            alpha[i] = x[(1 + (i - 1)) * (N_alleles + 1)] * beta[1] + theta[1];
        }

        lp = bernoulli_logit_glm_lpmf(y_obs | X, alpha, theta[2:N_alleles + 1]);

        return [lp]';
    }
}

data {
    int N_obs;
    int N_alleles;
    int N_shards;
    int ys[N_shards, N_obs + 1];
    real xs[N_shards, N_obs * (N_alleles + 1)];
}

parameters {
    vector[N_alleles] b_alleles_raw[N_shards];
    vector[N_shards] intercepts;
    real mu;
    real<lower=0> sigma;
    real<lower=0> phylogeny;
}

transformed parameters {
    vector[N_alleles + 1] theta[N_shards];
    vector[N_alleles] b_alleles[N_shards];
    vector[1] beta;

    for (i in 1:N_shards) {
        b_alleles[i] = mu + sigma * b_alleles_raw[i];
        theta[i][1] = intercepts[i];
        theta[i][2:(N_alleles + 1)] = b_alleles[i];

        beta[1] = phylogeny;
    }
}

model {
    mu ~ normal(0, 1);
    sigma ~ normal(0, 3);
    phylogeny ~ normal(0, 1);

    for (i in 1:N_shards) {
        b_alleles_raw[i] ~ std_normal();
    }

    target += sum(map_rect(ll, beta, theta, xs, ys));
}