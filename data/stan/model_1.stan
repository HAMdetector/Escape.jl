data {
    int N;
    int S;
    int D;
    int R;
    matrix[S, D] X;
    int y[N];
    int rs[N];
    int idx[N];
    matrix[R, S] phy;
    matrix[R, D] Z;
}

parameters {
    vector[R] b0_hla;
    vector[D] beta_hla[R];
}

model {
    b0_hla ~ normal(0, 10);
    for (i in 1:R) {
        beta_hla[i] ~ student_t(4, 0, 1);
    }

    {
        vector[N] theta_i;
        for (i in 1:N) {
            theta_i[i] = b0_hla[rs[i]] + X[idx[i]] * beta_hla[rs[i]];
        }

        y ~ bernoulli_logit(theta_i);
    }
}

generated quantities {
    matrix[R, S] theta = rep_matrix(-1, R, S);

    for (i in 1:N) {
        theta[rs[i], idx[i]] = inv_logit(b0_hla[rs[i]] + X[idx[i]] * beta_hla[rs[i]]);
    }
}