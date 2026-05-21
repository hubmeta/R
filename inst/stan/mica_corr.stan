// =============================================================
// MICA — Bayesian correlation-matrix completion
//
// Non-centered hierarchical model on Fisher-z observed cells with
// regression-informed priors for missing cells and PD-by-construction
// correlation draws via a Cholesky factor.
// =============================================================

data {
  int<lower=2> d;

  int<lower=0> n_obs;
  array[n_obs, 2] int<lower=1, upper=d> obs_idx;
  vector[n_obs] z_obs;
  vector<lower=0>[n_obs] within_prec_obs;

  int<lower=0> n_miss;
  array[n_miss, 2] int<lower=1, upper=d> miss_idx;
  vector[n_miss] prior_mean_z;
  vector<lower=0>[n_miss] prior_prec_z;

  real<lower=0> tau2_scale;
  real<lower=0> lkj_eta;
}

parameters {
  cholesky_factor_corr[d] L;
  real<lower=0> tau;
  vector[n_obs] true_z_raw;
}

transformed parameters {
  matrix[d, d] R = multiply_lower_tri_self_transpose(L);
  vector[n_obs] true_z;

  for (o in 1:n_obs) {
    int i = obs_idx[o, 1];
    int j = obs_idx[o, 2];
    real r_clip = fmax(fmin(R[i, j], 0.9999), -0.9999);
    true_z[o] = atanh(r_clip) + tau * true_z_raw[o];
  }
}

model {
  L ~ lkj_corr_cholesky(lkj_eta);
  tau ~ cauchy(0, tau2_scale);
  true_z_raw ~ std_normal();

  for (m in 1:n_miss) {
    int i = miss_idx[m, 1];
    int j = miss_idx[m, 2];
    real r_clip = fmax(fmin(R[i, j], 0.9999), -0.9999);
    real z_ij = atanh(r_clip);
    target += normal_lpdf(z_ij | prior_mean_z[m], inv_sqrt(prior_prec_z[m]));
  }

  for (o in 1:n_obs) {
    z_obs[o] ~ normal(true_z[o], inv_sqrt(within_prec_obs[o]));
  }
}

generated quantities {
  matrix[d, d] R_out = R;
}

