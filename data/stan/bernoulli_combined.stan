data {
    int N_y;
    int y[N_y];
    int ii[N_y];
    int jj[N_y];
    int N;
    int M;
    int R;
    matrix[N, M] H;
}

parameters {
    real intercept[R];
    real mu;
    real<lower=0> sigma;
    matrix<offset=mu, multiplier=sigma>[M, R] beta;
}

model {
    for (m in 1:M)
        beta[m,] ~ normal(mu, sigma);

    mu ~ normal(0, 2);
    sigma ~ normal(0, 2);
    intercept ~ normal(0, 3);

    for (i in 1:N_y) {
        y[i] ~ bernoulli_logit(intercept[jj[i]] + H[ii[i],] * beta[,jj[i]]);
    }
}