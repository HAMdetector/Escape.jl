// full model with HS prior, phylogeny and epitope prediction

functions {
    vector ll(vector global_pars, vector local_pars, real[] xs, int[] ys) {
        // extracting integer-valued data from ys
        int D = (num_elements(local_pars) - 4) / 4;
        int y_counts = ys[1];
        int S = ys[2];
        
        // extracting global pars
        real b_epi = global_pars[1];
        // real b_var = global_pars[2];

        // extracting local parameters
        real b0_hla = local_pars[1];
        real aux1_tau = local_pars[2];
        real aux2_tau = local_pars[3];
        real b_phy = local_pars[4];
        vector[D] aux1_lambda = local_pars[5:(D + 4)];
        vector[D] aux2_lambda = local_pars[(D + 5):(2*D + 4)];
        vector[D] z_std = local_pars[(2*D + 5):(3*D + 4)];
        vector[D] c2 = local_pars[(3*D + 5):(4*D + 4)];

        // model specification
        real lp = 0;
        lp += normal_lpdf(b0_hla | 0, 5);
        lp += std_normal_lpdf(z_std | );
        lp += std_normal_lpdf(aux1_lambda | );
        lp += std_normal_lpdf(aux1_tau | );
        lp += inv_gamma_lpdf(c2 | 2.5, 2.5);
        lp += inv_gamma_lpdf(aux2_tau | 0.5, 0.5 * square(xs[1]));
        lp += inv_gamma_lpdf(aux2_lambda | 0.5, (1 + b_epi * to_vector(segment(xs, 2, D))) * 0.5);
            
        {
            real tau = aux1_tau * sqrt(aux2_tau);
            vector[D] lambda = aux1_lambda .* sqrt(aux2_lambda);
            vector[D] lambda_tilde = sqrt((c2 .* square(lambda)) ./ 
                (c2 + square(tau) * square(lambda)));
            vector[D] beta_hla = z_std .* (tau * lambda_tilde);

            lp += bernoulli_logit_glm_lpmf(segment(ys, 3 + S, S)[segment(ys, 3, y_counts)] | 
                to_matrix(segment(xs, 2 + D + S, S * D), S, D)[segment(ys, 3, y_counts)], 
                b0_hla + b_phy * logit(to_vector(segment(xs, 2 + D, S))[segment(ys, 3, y_counts)]), beta_hla); //b_phy * logit(to_vector(segment(xs, 2 + D, S))[segment(ys, 3, y_counts)]) + b0_hla, 
        }

        return [lp]';
    }
}

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
    int y_counts_int[R] = rep_array(0, R);
    vector[R] y_means;
    vector[R] pseudo_variances;
    vector[R] pseudo_sigmas;
    vector[R] tau_0s;
    vector[D] s_j_sq;
    int ys[R, S] = rep_array(-1, R, S);
    int idxs[R, S] = rep_array(-1, R, S);

    int y_r[R, 2 + S + S]; // size(1); S(1); idx(S); y(S)
    real x_r[R, 1 + D + S + S * D + D]; // tau_0(1); Z(D); phy(S); X(S, D, column major); s_j_sq(D)

    // get size of y (y_counts) for each replacement, fill ys
    for (i in 1:N) {
        y_sums[rs[i]] += y[i];
        y_counts[rs[i]] += 1;
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

    // calculate tau_0 for each replacement
    for (i in 1:R) {
        y_means[i] = y_sums[i] / y_counts[i];
        pseudo_variances[i] = (1.0 / y_means[i]) * (1.0 / (1.0 - y_means[i]));
        pseudo_sigmas[i] = sqrt(pseudo_variances[i]);
        tau_0s[i] = (p0 / (D - p0)) * (pseudo_sigmas[i] / sqrt(y_counts[i]));
    }

    // variable variances
    for (i in 1:D) {
        s_j_sq[i] = variance(X[, i]);
    }

    // fill y_r and x_r, this is indexing madness
    for (i in 1:R) {
        y_r[i, 1] = y_counts_int[i];
        y_r[i, 2] = S;
        y_r[i, 3:(S + 2)] = idxs[i,];
        y_r[i, (S + 3):(2 + S * 2)] = ys[i,];

        x_r[i, 1] = tau_0s[i];
        x_r[i, 2:(D + 1)] = to_array_1d(Z[i,]);
        x_r[i, (D + 2):(1 + D + S)] = to_array_1d(phy[i,]);
        x_r[i, (2 + D + S):(1 + D + S + S * D)] = to_array_1d(X);
        x_r[i, (2 + D + S + S * D):(1 + D + S + S * D + D)] = to_array_1d(s_j_sq);
    }
}

parameters {
    vector[R] b0_hla;
    vector<lower=0>[R] aux1_tau;
    vector<lower=0>[R] aux2_tau;
    vector<lower=0>[D] aux1_lambda[R];
    vector<lower=0>[D] aux2_lambda[R];
    vector<lower=0>[D] c2[R];
    vector[D] z_std[R];

    real<lower=0> b_epi;
    real mu_phy;
    real<lower=0> sigma_phy;
    vector[R] b_phy;
}

model {
    mu_phy ~ normal(1, 1);
    sigma_phy ~ student_t(4, 0, 0.5);
    b_phy ~ normal(mu_phy, sigma_phy);
    b_epi ~ normal(0, 3);

    {
        vector[4 + 4 * D] local_pars[R];
        vector[1] global_pars;

        global_pars[1] = b_epi;

        for (i in 1:R) {
            local_pars[i][1] = b0_hla[i];
            local_pars[i][2] = aux1_tau[i];
            local_pars[i][3] = aux2_tau[i];
            local_pars[i][4] = b_phy[i];
            local_pars[i][5:(D + 4)] = aux1_lambda[i];
            local_pars[i][(D + 5):(D + D + 4)] = aux2_lambda[i];
            local_pars[i][(D + D + 5):(D + D + D + 4)] = z_std[i];
            local_pars[i][(D + D + D + 5):(D + D + D + D + 4)] = c2[i];
        }

        target += sum(map_rect(ll, global_pars, local_pars, x_r, y_r));
    }
}

generated quantities {
    vector[D] beta_hla[R];

    for (r in 1:R) {
        real tau = aux1_tau[r] * sqrt(aux2_tau[r]);
        vector[D] lambda = aux1_lambda[r] .* sqrt(aux2_lambda[r]);
        vector[D] lambda_tilde = sqrt(c2[r] .* square(lambda) ./ 
                (c2[r] + square(tau) * square(lambda)));
        beta_hla[r] = z_std[r] .* (tau * lambda_tilde);
    }
}