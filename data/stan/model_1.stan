// basic logistic regression model without HS prior

functions {
    vector ll(vector global_pars, vector local_pars, real[] xs, int[] ys) {
        // extracting integer-valued data from ys
        int D = num_elements(local_pars) - 1;
        int y_counts = ys[1];
        int S = ys[2];

        // extracting local parameters
        real b0_hla = local_pars[1];
        vector[D] beta_hla = local_pars[2:(D + 1)];

        // model specification
        real lp = 0;
        lp += normal_lpdf(b0_hla | 0, 5);
        lp += student_t_lpdf(beta_hla | 7, 0, 1);

        lp += bernoulli_logit_glm_lpmf(segment(ys, 3 + S, S)[segment(ys, 3, y_counts)] | 
            to_matrix(segment(xs, 1, S * D), S, D)[segment(ys, 3, y_counts)], 
            b0_hla, beta_hla);

        return [lp]';
    }
}

data {
    int N;
    int S;
    int D;
    int R;
    real p0; // unused in model 1
    matrix[S, D] X;
    int y[N];
    int rs[N];
    int idx[N];
    matrix[R, S] phy; // unused in model 1
    matrix[R, D] Z; // unused in model 1
}

transformed data {
    vector[R] y_counts = rep_vector(0, R);
    int y_counts_int[R] = rep_array(0, R);
    int ys[R, S] = rep_array(-1, R, S);
    int idxs[R, S] = rep_array(-1, R, S);

    int y_r[R, 2 + S + S]; // size(1); S(1); idx(S); y(S)
    real x_r[R, S * D]; // X(S, D, column major) 

    // get size of y (y_counts) for each replacement, fill ys
    for (i in 1:N) {
        y_counts[rs[i]] += 1;
        y_counts_int[rs[i]] += 1;
        ys[rs[i], idx[i]] = y[i]; 
    }

    // fill idxs
    for (i in 1:R) {
        int count = 1;
        for (j in 1:S) {
            if (ys[i, j] != -1) {
                idxs[i, count] = j;
                count += 1;
            }
        }
    }

    // fill y_r and x_r, this is indexing madness
    for (i in 1:R) {
        y_r[i, 1] = y_counts_int[i];
        y_r[i, 2] = S;
        y_r[i, 3:(S + 2)] = idxs[i,];
        y_r[i, (S + 3):(2 + S * 2)] = ys[i,];
        x_r[i, :] = to_array_1d(X);
    }
}

parameters {
    vector[R] b0_hla;
    vector<lower=0>[D] beta_hla[R];
}

model {
    {
        vector[1 + D] local_pars[R];
        vector[0] global_pars;

        for (i in 1:R) {
            local_pars[i][1] = b0_hla[i];
            local_pars[i][2:(D + 1)] = beta_hla[i];
        }

        target += sum(map_rect(ll, global_pars, local_pars, x_r, y_r));
    }
}