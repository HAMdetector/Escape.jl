// simple logistic regression model, without phylogeny and epitope prediction

functions {
    vector ll(vector global_pars, vector local_pars, real[] xs, int[] ys) {
        // extracting integer-valued data from ys
        int y_counts = ys[1];
        int S = ys[2];
        int D = ys[3];

        // extracting global pars
        real tau = global_pars[1];

        // extracting local parameters
        real b0_hla = local_pars[1];
        real c2 = local_pars[2];
        real b_phy = local_pars[3];
        
        vector[D] aux1_lambda = local_pars[4:(D + 3)];
        vector[D] aux2_lambda = local_pars[(4 + D):(3 + 2 * D)];
        vector[D] z_std = local_pars[(4 + 2 * D):(3 + 3 * D)];

        // model specification
        real lp = 0;

        lp += normal_lpdf(b0_hla | 0, 100);
        lp += inv_gamma_lpdf(c2 | 3.5, 3.5);
        lp += std_normal_lpdf(z_std | );
        lp += std_normal_lpdf(aux1_lambda | );
        lp += inv_gamma_lpdf(aux2_lambda | 0.5, 0.5);

        {
            vector[D] s_j_sq = to_vector(segment(xs, 1 + D + S + S * D, D));
            vector[D] lambda = (aux1_lambda .* sqrt(aux2_lambda)) .* s_j_sq;
            vector[D] lambda_tilde = sqrt((c2 * square(lambda)) ./
                (c2 + square(tau) * square(lambda)));
            vector[D] beta_hla = z_std .* (tau * lambda_tilde);
 
            lp += bernoulli_logit_glm_lpmf(segment(ys, 4 + S, S)[segment(ys, 4, y_counts)] | 
                to_matrix(segment(xs, 1 + D + S, S * D), S, D)[segment(ys, 4, y_counts)], 
                b0_hla + b_phy * logit(to_vector(segment(xs, 1 + D, S))[segment(ys, 4, y_counts)]), 
                beta_hla); 
        }

        return [lp]';
    }
}

data {
    int N;
    int S;
    int D;
    int R;
    int y[N];
    int rs[N];
    int idx[N];
    matrix[S, D] X;
    matrix[R, S] phy;
    matrix[R, D] Z;
}

transformed data {
    int y_counts_int[R] = rep_array(0, R);
    vector[D] s_j_sq;
    int ys[R, S] = rep_array(-1, R, S);
    int idxs[R, S] = rep_array(-1, R, S);

    int y_r[R, 3 + 2 * S]; // size(1); S(1); D(1); idx(S); y(S)
    real x_r[R, D + S + S * D + D]; // Z(D); phy(S); X(S, D, column major);  s_j_sq(D)

    // get size of y (y_counts) for each replacement, fill ys
    for (i in 1:N) {
        y_counts_int[rs[i]] += 1;
        ys[rs[i], idx[i]] = y[i]; 
    }

    // fill idxs
    for (i in 1:R) {
        int count = 0;
        for (j in 1:S) {
            if (ys[i, j] != -1) {
                count += 1;
                idxs[i, count] = j;
            }
        }
    }

    // variable variances
    for (i in 1:D) {
        s_j_sq[i] = variance(X[, i]);
    }

    // fill y_r and x_r
    for (i in 1:R) {
        y_r[i, 1] = y_counts_int[i];
        y_r[i, 2] = S;
        y_r[i, 3] = D;
        y_r[i, 4:(S + 3)] = idxs[i,];
        y_r[i, (4 + S):(3 + 2 * S)] = ys[i,];

        x_r[i, 1:D] = to_array_1d(Z[i,]);
        x_r[i, (1 + D):(D + S)] = to_array_1d(phy[i,]);
        x_r[i, (1 + D + S):(D + S + S * D)] = to_array_1d(X);
        x_r[i, (1 + D + S + S * D):(D + S + S * D + D)] = to_array_1d(s_j_sq);
    }
}

parameters {
    vector[R] b0_hla;
    vector<lower=0>[R] c2;
    vector<lower=0>[D] aux1_lambda[R];
    vector<lower=0>[D] aux2_lambda[R];
    vector[D] z_std[R];

    real<lower=0> aux1_tau;
    real<lower=0> aux2_tau;

    vector[R] b_phy;
    real mu_phy;
    real<lower=0> sigma_phy;
}

model {
    aux1_tau ~ std_normal();
    aux2_tau ~ inv_gamma(0.5, 0.5);

    sigma_phy ~ normal(0, 0.5);
    mu_phy ~ normal(1, 1);
    b_phy ~ normal(mu_phy, sigma_phy);

    {
        vector[3 + 3 * D] local_pars[R];
        vector[1] global_pars;
        
        global_pars[1] = aux1_tau * sqrt(aux2_tau) * ((2.0 / (D - 2.0)) * (2.0 / sqrt(S)));

        for (r in 1:R) {
            local_pars[r][1] = b0_hla[r];
            local_pars[r][2] = c2[r];
            local_pars[r][3] = b_phy[r];
            local_pars[r][4:(D + 3)] = aux1_lambda[r];
            local_pars[r][(4 + D):(3 + 2 * D)] = aux2_lambda[r];
            local_pars[r][(4 + 2 * D):(3 + 3 * D)] = z_std[r];
        }

        target += sum(map_rect(ll, global_pars, local_pars, x_r, y_r));
    }
}

generated quantities {
    vector[D] beta_hla[R];

    real tau = aux1_tau * sqrt(aux2_tau) * ((2.0 / (D - 2.0)) * (2.0 / sqrt(S)));

    for (r in 1:R) {
        vector[D] lambda = (aux1_lambda[r] .* sqrt(aux2_lambda[r])) .* s_j_sq;
        vector[D] lambda_tilde = sqrt(c2[r] * square(lambda) ./ 
                (c2[r] + square(tau) * square(lambda)));
        beta_hla[r] = z_std[r] .* (tau * lambda_tilde);
    }
}
