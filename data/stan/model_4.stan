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
        pseudo_variances[i] = (1.0 / y_means[i]) * (1.0 / (1.0 - y_means[i]));
        pseudo_sigmas[i] = sqrt(pseudo_variances[i]);
        tau_0s[i] = (p0 / (D - p0)) * (pseudo_sigmas[i] / sqrt(N_obs[i]));
    } 
}

parameters {
    vector[1] phylogeny;
    vector[R] b0;
    vector[R] gamma0;
    real gamma;
    vector<lower=0>[R] aux1_tau;
    vector<lower=0>[R] aux2_tau;
    vector<lower=0>[D] aux1_lambda[R];
    vector<lower=0>[D] aux2_lambda[R];
    vector[D] z[R];
}

transformed parameters {
    vector[D + 1] beta[R];

    for (i in 1:R) {
        real tau;
        vector[D] lambda;

        tau = aux1_tau[i] * sqrt(aux2_tau[i]); // implies cauchy(0, tau_0^2) prior on tau
        lambda = aux1_lambda[i] .* sqrt(aux2_lambda[i]); // implies cauchy(0, 1) on lambda
        beta[i][1] = b0[i];
        beta[i][2:(D + 1)] = z[i] .* (lambda * tau); // implies normal(0, tau * lambda)
    }
}

model {
    phylogeny ~ normal(0, 5);
    b0 ~ normal(0, 3);
    gamma0 ~ std_normal();
    gamma ~ std_normal();

    for (i in 1:R) {
        real in_epitope[D] = segment(xs[i, ], 4 + N_obs[i] * D, D);

        z[i] ~ std_normal();
        aux1_tau[i] ~ std_normal();
        aux1_lambda[i] ~ std_normal();

        aux2_tau[i] ~ inv_gamma(0.5, square(tau_0s[i]) / 2);
        for (j in 1:D) {
            aux2_lambda[i][j] ~ inv_gamma(0.5, exp(gamma0[i] + gamma * in_epitope[j]) / 2);
        }
    }

    target += sum(map_rect(ll, phylogeny, beta, xs, ys));
}

generated quantities {
    vector[D] omega[R];
    vector[R] m_eff;
    matrix[R, N] theta = rep_matrix(-1, R, N);

    for (i in 1:R) {
        real tau = aux1_tau[i] * sqrt(aux2_tau[i]);
        m_eff[i] = 0;

        for (j in 1:D) {
            real lambda = aux1_lambda[i][j] * sqrt(aux2_lambda[i][j]);

            omega[i][j] = 1 - (1.0 / (1 + N_obs[i] * 1.0 / square(pseudo_sigmas[i]) * 
                square(lambda) * square(tau)));
            m_eff[i] += omega[i][j];
        }
    }

    for (i in 1:R) {
        matrix[N_obs[i], D] X = to_matrix(segment(xs[i, ], 4, N_obs[i] * D), N_obs[i], D);
        vector[N_obs[i]] phy = to_vector(segment(xs[i, ], 4 + N_obs[i] * D + D, N_obs[i]));

        for (j in 1:N_obs[i]) {
            theta[i, j] = inv_logit(phy[j] * phylogeny[1] + beta[i][1] + 
                X[j] * beta[i][2:(D + 1)]);
        }
    }
}