functions {
    vector ll(vector global_beta, vector local_beta, real[] xs, int[] ys) {
        int N_obs = ys[1];
        int y[N_obs] = segment(ys, 2, N_obs);
        int D = num_elements(local_beta) - 1;
        matrix[N_obs, D] X = to_matrix(segment(xs, 4, N_obs * D), N_obs, D);

        real lp = bernoulli_logit_glm_lpmf(y | X, local_beta[1], local_beta[2:D + 1]);

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

parameters {
    vector[0] placeholder;
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
    b0 ~ normal(0, 3);

    for (i in 1:R) {
        beta_hla[i] ~ student_t(4, 0, 1);
    }

    target += sum(map_rect(ll, placeholder, beta, xs, ys));
}

generated quantities {
    matrix[R, N] theta = rep_matrix(-1, R, N);

    for (i in 1:R) {
        int N_obs = ys[i, 1];
        int y[N_obs] = segment(ys[i, ], 2, N_obs);
        matrix[N_obs, D] X = to_matrix(segment(xs[i, ], 4, N_obs * D), N_obs, D);

        for (j in 1:N_obs) {
            theta[i, j] = inv_logit(b0[i] + X[j] * beta_hla[i]);
        }
    }
}