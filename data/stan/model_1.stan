data {
    int<lower=1> N;
    int<lower=3> D;
    int<lower=1> R;
    int ys[R, N + 1];
    matrix[R, N * D] xs;
}

parameters {
    vector[R] intercepts;
    vector[D] beta_hla[R];
}

model {
    intercepts ~ normal(0, 20);

    for (i in 1:R) {
        beta_hla[i] ~ normal(0, 3);

        {
            matrix[ys[i, 1], D] X;
            int y[ys[i, 1]];

            for (j in 1:ys[i, 1]) {
                X[j,] = xs[i, ][(1 + (j - 1) * D):(j * D)];
                y[j] = ys[i, 1 + j];
            }

            y ~ bernoulli_logit(intercepts[i] + X * beta_hla[i]);
        }
    }
}