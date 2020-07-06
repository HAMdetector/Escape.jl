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
    int y_counts_int[R];
    vector[R] y_means;
    vector[R] pseudo_variances;
    vector[R] tau_0s;

    for (i in 1:N) {
        y_sums[rs[i]] += y[i];
        y_counts[rs[i]] += 1;
        y_counts_int[rs[i]] += 1;
    }

    for (i in 1:R) {
        y_means[i] = y_sums[i] / y_counts[i];
        pseudo_variances[i] = (1.0 / y_means[i]) * (1.0 / (1.0 - y_means[i]));
        tau_0s[i] = (p0 / (D - p0)) * (sqrt(pseudo_variances[i]) / sqrt(y_counts[i]));
    }
}

parameters {
    vector[R] b0;
    vector[D] beta_hla[R];
}

model {
    b0 ~ normal(0, 10);

    for (i in 1:1) {
        int ii = 1;
        int y_[y_counts_int[i]];
        matrix[y_counts_int[i], D] X_;

        for (n in 1:N) {
            if (rs[n] == i) {
                X_[ii] = X[idx[n]];
                y_[ii] = y[n];
                ii += 1;
            }
        }

        beta_hla[i] ~ student_t(7, 0, 1);
        y_ ~ bernoulli_logit_glm(X_, b0[i], beta_hla[i]);
    }
}