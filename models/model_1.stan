// simple logistic regression model, without phylogeny and epitope prediction

functions {
    vector ll(vector global_pars, vector local_pars, array[] real xs, array[] int ys) {
        // extracting integer-valued data from ys
        int y_counts = ys[1];
        int S = ys[2];
        int D = ys[3];

        // extracting local parameters
        real b0_hla = local_pars[1];
        vector[D] beta_hla = to_vector(local_pars[2:(D + 1)]);

        // model specification
        real lp = 0;
        lp += normal_lpdf(b0_hla | 0, 100);
        lp += student_t_lpdf(beta_hla | 7, 0, 1);

        lp += bernoulli_logit_glm_lpmf(segment(ys, 4 + S, S)[segment(ys, 4, y_counts)] | 
            to_matrix(segment(xs, 2 + D + S, S * D), S, D)[segment(ys, 4, y_counts)], 
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
    matrix[S, D] X;
    array[N] int y;
    array[N] int rs;
    array[N] int idx;
    matrix[R, S] phy;
    matrix[R, D] Z;
}

transformed data {
    array[R] int y_counts_int = rep_array(0, R);
    vector[D] s_j_sq;
    array[R, S] int ys = rep_array(-1, R, S);
    array[R, S] int idxs = rep_array(-1, R, S);

    array[R, 3 + S + S] int y_r; // size(1); S(1); D(1); idx(S); y(S)
    array[R, 1 + D + S + S * D + D] real x_r; // tau_0(1); Z(D); phy(S); X(S, D, column major); 
        // s_j_sq(D)

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
        y_r[i, (S + 4):(3 + S * 2)] = ys[i,];

        x_r[i, 1] = 2;//tau_0s[i];
        x_r[i, 2:(D + 1)] = to_array_1d(Z[i,]);
        x_r[i, (D + 2):(1 + D + S)] = to_array_1d(phy[i,]);
        x_r[i, (2 + D + S):(1 + D + S + S * D)] = to_array_1d(X);
        x_r[i, (2 + D + S + S * D):(1 + D + S + S * D + D)] = to_array_1d(s_j_sq);
    }
}

parameters {
    vector[R] b0_hla;
    array[R] vector[D] beta_hla;
}

model {
    {
        array[R] vector[D + 1] local_pars;
        vector[0] global_pars;

        for (i in 1:R) {
            local_pars[i][1] = b0_hla[i];
            local_pars[i][2:(D + 1)] = beta_hla[i];
        }

        target += sum(map_rect(ll, global_pars, local_pars, x_r, y_r));
    }
}

generated quantities {

}
