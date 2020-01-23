functions {
    vector ll(vector beta, vector theta, real[] x, int[] y) {
        int N = y[1];
        int D = num_elements(theta) - 1;
        matrix[N, D] X;
        real lp;

        for (i in 1:N) {
            X[i] = to_vector(x[(1 + (i - 1) * D) : (i * D)])';
        }

        // used for error checking
        // for (i in 1:N) {
        //     for (j in 1:D) {
        //         if (X[i, j] < 0) {
        //             reject("assertion error");
        //         }
        //     }

        //     if (y[2:(N + 1)][i] < 0) {
        //         reject("assertion error")
        //     }
        // }

        lp = bernoulli_logit_glm_lpmf(y[2:(N + 1)] | X, theta[1], theta[2:D + 1]);

        return [lp]';
    }
}

data {
    int N;
    int D;
    int R;
    int ys[R, N + 1];
    real y_mean[R];
    real hla_mean[D];
    real xs[R, N * D];
    real p0;
}

transformed data {
    real slab_scale = 3;
    real slab_df = 4;
}

parameters {
    vector[R] intercepts;
    vector[D] z[R];
    real<lower=0> aux1_global[R];
    real<lower=0> aux2_global[R];
    vector<lower=0>[D] aux1_local[R];
    vector<lower=0>[D] aux2_local[R];
    real<lower=0> caux[R];
    vector[0] beta;
}

transformed parameters {
    real<lower=0> tau[R];
    vector<lower=0>[D] lambda[R];
    vector<lower=0>[D] lambda_tilde[R];
    real<lower=0> c[R];
    vector[D] beta_hla[R];
    vector[D + 1] theta[R];

    {
        for (i in 1:R) {
            real scale_global = (p0 / (D - p0)) * 
                (1.0 / sqrt(((y_mean[i] * (1 - y_mean[i])))) / sqrt(1.0 * N));

            lambda[i] = aux1_local[i] .* sqrt(aux2_local[i]); // implies cauchy distribution
            tau[i] = aux1_global[i] * sqrt(aux2_global[i]) * scale_global;
            c[i] = slab_scale * sqrt(caux[i]);
            lambda_tilde[i] = sqrt(c[i]^2 * square(lambda[i]) ./ 
                (c[i]^2 + tau[i]^2 * square(lambda[i])));
            beta_hla[i] = z[i] .* lambda_tilde[i] * tau[i];

            theta[i][1] = intercepts[i];
            theta[i][2:D + 1] = beta_hla[i];
        }
    }
}


model {
    for (i in 1:R) {
        z[i] ~ normal(0, 1);
        aux1_local[i] ~ normal(0, 1);
        aux1_global[i] ~ normal(0, 1);
        aux2_global[i] ~ inv_gamma(0.5, 0.5);
        caux[i] ~ inv_gamma(0.5* slab_df, 0.5 * slab_df);
        intercepts[i] ~ normal(0, 20);

        for (j in 1:D) {
            aux2_local[i][j] ~ inv_gamma(0.5, 0.5 * ((1.0 / ((hla_mean[j]) * (1 - hla_mean[j])))^2));
        }
    }

    target += sum(map_rect(ll, beta, theta, xs, ys));
}

generated quantities {
    vector[N] log_lik[R];

    for (i in 1:R) {
        int samples = ys[i, 1];
        for (n in 1:samples) {
            int y;
            row_vector[D] x_i;

            y = ys[i, n + 1];
            x_i = to_vector(xs[i, (1 + (n - 1) * D) : (n * D)])';

            log_lik[i][n] = bernoulli_logit_lpmf(y| intercepts[i] + x_i * beta_hla[i]);
        }
    }
}