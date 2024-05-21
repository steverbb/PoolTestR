// generated with brms 2.21.2
functions {
  /* compute a latent Gaussian process
   * Args:
   *   x: array of continuous predictor values
   *   sdgp: marginal SD parameter
   *   lscale: length-scale parameter
   *   zgp: vector of independent standard normal variables
   * Returns:
   *   a vector to be added to the linear predictor
   */
  vector gp(data array[] vector x, real sdgp, vector lscale, vector zgp) {
    int Dls = rows(lscale);
    int N = size(x);
    matrix[N, N] cov;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      cov = gp_exponential_cov(x, sdgp, lscale[1]);
    } else {
      // multi-dimensional non-isotropic GP
      cov = gp_exponential_cov(x[, 1], sdgp, lscale[1]);
      for (d in 2:Dls) {
        cov = cov .* gp_exponential_cov(x[, d], 1, lscale[d]); // .* is element-wise multiplication
      }
    }
    for (n in 1:N) {
      // deal with numerical non-positive-definiteness
      cov[n, n] += 1e-12;
    }
    return cholesky_decompose(cov) * zgp;
  }
  /* Spectral density function of a Gaussian process
   * with exponential covariance kernel
   * Args:
   *   x: array of numeric values of dimension NB x D
   *   sdgp: marginal SD parameter
   *   lscale: vector of length-scale parameters
   * Returns:
   *   numeric values of the function evaluated at 'x'
   */
  vector spd_cov_exp_quad(data array[] vector x, real sdgp, vector lscale) {
    int NB = dims(x)[1];
    int D = dims(x)[2];
    int Dls = rows(lscale);
    vector[NB] out;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      real constant = square(sdgp) * pow(2, D) * pow(pi(), (D - 1) / 2.0) * tgamma((D + 1) / 2.0) / lscale[1];
      real inverse_l2 = 1 / square(lscale[1]);
      real four_pi2 = 4 * square(pi());
      real power = -(D + 1) / 2.0;
      for (m in 1:NB) {
        out[m] = constant * pow(inverse_l2 + four_pi2 * dot_self(x[m]), power);
      }
    } else {
      // multi-dimensional non-isotropic GP
      real norm_l = sqrt(dot_self(lscale));
      real constant = square(sdgp) * pow(2, D) * pow(pi(), (D - 1) / 2.0) * tgamma((D + 1) / 2.0) / norm_l;
      real inverse_l2 = 1 / square(norm_l);
      real four_pi2 = 4 * square(pi());
      real power = -(D + 1) / 2.0;
      for (m in 1:NB) {
        out[m] = constant * pow(inverse_l2 + four_pi2 * dot_self(x[m]), power);
      }
    }
    return out;
  }
  /* how many elements are in a range of integers?
   * Args:
   *   x: an integer array
   *   start: start of the range (inclusive)
   *   end: end of the range (inclusive)
   * Returns:
   *   a scalar integer
   */
  int size_range(array[] int x, int start, int end) {
    int out = 0;
    for (i in 1:size(x)) {
      out += (x[i] >= start && x[i] <= end);
    }
    return out;
  }
  /* which elements are in a range of integers?
   * Args:
   *   x: an integer array
   *   start: start of the range (inclusive)
   *   end: end of the range (inclusive)
   * Returns:
   *   an integer array
   */
  array[] int which_range(array[] int x, int start, int end) {
    array[size_range(x, start, end)] int out;
    int j = 1;
    for (i in 1:size(x)) {
      if (x[i] >= start && x[i] <= end) {
        out[j] = i;
        j += 1;
      }
    }
    return out;
  }
  /* adjust array values to x - start + 1
   * Args:
   *   x: an integer array
   *   start: start of the range of values in x (inclusive)
   * Returns:
   *   an integer array
   */
  array[] int start_at_one(array[] int x, int start) {
    array[size(x)] int out;
    for (i in 1:size(x)) {
      out[i] = x[i] - start + 1;
    }
    return out;
  }

    real pool_bernoulli_logit_lpmf(int y, real alpha, int N) {
      return bernoulli_lpmf(1-y | exp(log1m_inv_logit(alpha) * N));
    }
    int pool_bernoulli_logit_rng(real alpha, int N) {
      return 1-bernoulli_rng(exp(log1m_inv_logit(alpha) * N));
    }

}
data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  // data for custom integer vectors
  array[N] int vint1;
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  // data related to GPs
  int<lower=1> Kgp_1;  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Dgp_1;  // GP dimension
  // number of latent GP groups
  int<lower=1> Nsubgp_1;
  // indices of latent GP groups per observation
  array[N] int<lower=1> Jgp_1;
  array[Nsubgp_1] vector[Dgp_1] Xgp_1;  // covariates of the GP
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  vector<lower=0>[Kgp_1] sdgp_1;  // GP standard deviation parameters
  array[Kgp_1] vector<lower=0>[1] lscale_1;  // GP length-scale parameters
  vector[Nsubgp_1] zgp_1;  // latent variables of the GP
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(b | 6, 0, 1.5);
  lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
  lprior += student_t_lpdf(sdgp_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += inv_gamma_lpdf(lscale_1[1][1] | 33.492838, 21.513469);
}
model {
  // likelihood including constants
  if (!prior_only) {
    vector[Nsubgp_1] gp_pred_1 = gp(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Xc * b + gp_pred_1[Jgp_1];
    for (n in 1:N) {
      target += pool_bernoulli_logit_lpmf(Y[n] | mu[n], vint1[n]);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zgp_1);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
