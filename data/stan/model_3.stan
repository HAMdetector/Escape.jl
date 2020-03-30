data {
    int N;
    int S;
    int D;
    int R;
    real p0;
    matrix[S, D] X;
    int y[N];
    int rs[N];
    int idx[N];
    matrix[R, S] phy;
    matrix[R, D] Z;
}

transformed data {
    vector[R] y_sums = rep_vector(0, R);
    vector[R] y_counts = rep_vector(0, R);
    vector[R] y_means;
    vector[R] pseudo_variances;
    vector[R] pseudo_sigmas;
    vector[R] tau_0s;
    vector[R] half_squared_tau_0s;
    
    for (i in 1:N) {
        y_sums[rs[i]] += y[i];
        y_counts[rs[i]] += 1; 
    }

    for (i in 1:R) {
        y_means[i] = y_sums[i] / y_counts[i];
        pseudo_variances[i] = (1.0 / y_means[i]) * (1.0 / (1.0 - y_means[i]));
        pseudo_sigmas[i] = sqrt(pseudo_variances[i]);
        tau_0s[i] = (p0 / (D - p0)) * (pseudo_sigmas[i] / sqrt(y_counts[i]));
        half_squared_tau_0s[i] = square(tau_0s[i]) / 2;
    }
}

parameters {
    real phylogeny;
    vector[R] b0_hla;
    vector<lower=0>[R] aux1_tau;
    vector<lower=0>[R] aux2_tau;
    vector<lower=0>[D] aux1_lambda[R];
    vector<lower=0>[D] aux2_lambda[R];
    vector[D] z[R];
}

transformed parameters {
    vector[D] beta_hla[R];
    for (i in 1:R) {
        real tau;
        vector[D] lambda;

        tau = aux1_tau[i] * sqrt(aux2_tau[i]);
        lambda = aux1_lambda[i] .* sqrt(aux2_lambda[i]);
        beta_hla[i] = z[i] .* (lambda * tau);
    }
}

model {
    b0_hla ~ normal(0, 10);
    phylogeny ~ normal(0, 3);

    for (i in 1:R) {
        z[i] ~ std_normal();
        aux1_lambda[i] ~ std_normal();
        aux1_tau[i] ~ std_normal();

        aux2_lambda[i] ~ inv_gamma(0.5, 0.5);
        aux2_tau[i] ~ inv_gamma(0.5, half_squared_tau_0s[i]);
    }

    {
        vector[N] theta_i;
        for (i in 1:N) {
            theta_i[i] = phylogeny * inv_logit(phy[rs[i], idx[i]]) +
                b0_hla[rs[i]] +
                X[idx[i]] * beta_hla[rs[i]];
        }

        y ~ bernoulli_logit(theta_i);
    }
}

generated quantities {
    matrix[R, S] theta = rep_matrix(-1, R, S);
    vector[D] omega[R];
    vector[R] m_eff;

    for (i in 1:N) {
        theta[rs[i], idx[i]] = inv_logit(b0_hla[rs[i]] + X[idx[i]] * beta_hla[rs[i]]);
    }

    for (i in 1:R) {
        real tau = aux1_tau[i] * sqrt(aux2_tau[i]);
        m_eff[i] = 0;

        for (j in 1:D) {
            real lambda = aux1_lambda[i][j] * sqrt(aux2_lambda[i][j]);
            omega[i][j] = 1.0 - (1.0 / (1.0 + y_counts[i] / square(pseudo_sigmas[i]) *
                square(lambda) * square(tau)));
            m_eff[i] += omega[i][j];
        }
    }
}