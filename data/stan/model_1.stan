functions {
    vector ll(vector global_pars, vector local_pars, real[] xs, int[] ys) {
        // extracting integer-valued data from ys
        int D = (num_elements(local_pars) - 1) / 2;
        int y_counts = ys[1];
        int S = ys[2];
        int idxs[y_counts] = segment(ys, 3, y_counts);
        int y[S] = segment(ys, 3 + S, S);

        // extracting local parameters
        real b0_hla = local_pars[1];
        vector[D] z_std = local_pars[2:(D + 1)];
        vector[D] lambda = local_pars[(2 + D):(2*D + 1)];
        vector[D] beta_hla = z_std .* sqrt(lambda);

        // model specification
        real lp = 0;
        lp += std_normal_lpdf(z_std | );
        lp += inv_gamma_lpdf(lambda | 2, 2);
        lp += normal_lpdf(b0_hla | 0, 10);

        lp += bernoulli_logit_glm_lpmf(y[idxs] | 
            to_matrix(segment(xs, 2 + D + S, S * D), S, D)[idxs], 
            b0_hla, 
            beta_hla);

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
    vector[R] half_squared_tau_0s;
    int ys[R, S] = rep_array(-1, R, S);
    int idxs[R, S] = rep_array(-1, R, S);

    int y_r[R, 2 + S + S]; // size(1); S(1); idx(S); y(S)
    real x_r[R, 1 + D + S + S * D]; // tau_0(1); Z(D); phy(S); X(S, D, column major) 

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
        half_squared_tau_0s[i] = square(tau_0s[i]) / 2;
    }

    // fill y_r and x_r, this is indexing madness
    for (i in 1:R) {
        y_r[i, 1] = y_counts_int[i];
        y_r[i, 2] = S;
        y_r[i, 3:(S + 2)] = idxs[i,];
        y_r[i, (S + 3):(2 + S * 2)] = ys[i,];

        x_r[i, 1] = half_squared_tau_0s[i];
        x_r[i, 2:(D + 1)] = to_array_1d(Z[i,]);
        x_r[i, (D + 2):(1 + D + S)] = to_array_1d(phy[i,]);
        x_r[i, (2 + D + S):(1 + D + S + S * D)] = to_array_1d(X);
    }
}

parameters {
    vector[R] b0_hla;
    vector<lower=0>[D] lambda[R];
    vector[D] z_std[R];
    vector[0] global_pars;
}

transformed parameters {
    vector[2*D + 1] local_pars[R];
    for (i in 1:R) {
        local_pars[i][1] = b0_hla[i];
        local_pars[i][2:(D + 1)] = z_std[i];
        local_pars[i][(2 + D):(2*D + 1)] = lambda[i]; 
    }
}

model {
    target += sum(map_rect(ll, global_pars, local_pars, x_r, y_r));
}

generated quantities {
    vector[N] theta;
    vector[D] beta_hla[R];

    for (i in 1:R) {
        beta_hla[i] = z_std[i] .* sqrt(lambda[i]);
    }

    for (i in 1:N) {
        theta[i] = inv_logit(
            b0_hla[rs[i]] + 
            X[idx[i]] * beta_hla[rs[i]]
        );
    }
}