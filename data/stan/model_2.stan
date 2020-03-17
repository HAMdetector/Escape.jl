functions {
    vector ll(vector global_beta, vector local_beta, real[] xs, int[] ys) {
        int N_obs = ys[1];
        int y[N_obs] = segment(ys, 2, N_obs);
        int D = num_elements(local_beta) - 1;
        matrix[N_obs, D] X = to_matrix(segment(xs, 4, N_obs * D), N_obs, D);
        vector[N_obs] phylogeny = to_vector(segment(xs, 4 + N_obs * D + D, N_obs));

        real lp = bernoulli_logit_glm_lpmf(y | X, global_beta[1] * phylogeny + local_beta[1],
            local_beta[2:D + 1]
        );

        return [lp]';
    }
}

data {
    int N;
    int D;
    int R;
    real p0;
    int ys[R, N + 1];
    real xs[R, 3 + N * D + D + N];
}

transformed data {
    int N_obs[R];
    vector[R] y_means;
    vector[R] pseudo_variances;
    vector[R] pseudo_sigmas;
    vector[R] tau_0s;
    vector[D] y[R];

    for (i in 1:R) {
        N_obs[i] = ys[i, 1];
        y[i] = segment(to_vector(ys[i, ]), 2, N_obs[i]);
        y_means[i] = mean(y[i]);
        pseudo_variances[i] = (1 / y_means[i]) * (1 / (1 - y_means[i]));
        pseudo_sigmas[i] = sqrt(pseudo_variances[i]);
        tau_0s[i] = (p0 / (D - p0)) * (pseudo_sigmas[i] / sqrt(N_obs[i]));
    } 
}

parameters {
    vector[1] phylogeny;
    vector[R] b0;
    vector[D] beta_hla[R];
}

transformed parameters {
    vector[D + 1] beta[R];

    for (i in 1:R) {
        beta[i][1] = b0[i];
        beta[i][2:(D + 1)] = beta_hla[i];
    }
}

model {
    phylogeny ~ normal(0, 5);
    b0 ~ normal(0, 3);

    for (i in 1:R) {
        beta_hla[i] ~ std_normal();
    }

    target += sum(map_rect(ll, phylogeny, beta, xs, ys));
}

generated quantities {
    matrix[R, N] theta = rep_matrix(-1, R, N);

    for (i in 1:R) {
        matrix[N_obs[i], D] X = to_matrix(segment(xs[i, ], 4, N_obs[i] * D), N_obs[i], D);
        vector[N_obs[i]] phy = to_vector(segment(xs[i, ], 4 + N_obs[i] * D + D, N_obs[i]));

        for (j in 1:N_obs[i]) {
            theta[i, j] = inv_logit(phy[j] * phylogeny[1] + b0[i] + X[j] * beta_hla[i]);
        }
    }
}