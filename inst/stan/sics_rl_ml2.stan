data {
  int<lower=0> N; // number of subjects
  int<lower=0> T; // total number of trials (per subject)
  array[T,N] int choice;
  array[T,N] int reward;
  array[T,N] int reset;
  array[T,N] int<lower=0> k; // number of times the current choice has been chosen, including the current trial
  int<lower=0> n;            // dimension of learning rate parameter >= max(k)


  // Parameters for hyperprior for population parameters for alpha and beta
  //// NOTE: Same hyperpriors for all alpha_j, 1 =< j <= n.
  // alpha
  real hyper_mu_mu_alpha;
  real hyper_sigma_mu_alpha;
  real hyper_mu_sigma_alpha;
  real hyper_sigma_sigma_alpha;
  // beta
  real hyper_mu_mu_beta;
  real hyper_sigma_mu_beta;
  real hyper_mu_sigma_beta;
  real hyper_sigma_sigma_beta;
  // Correlation
  real hyper_eta;

  int sample_posterior;  // Toggle for prior/posterior sampling

  /*vvv SICS vvv*/
  int<lower=0> SICS_n;
  int<lower=0> SICS_N;
  // A: equality constraints
  int<lower=0> SICS_mA;
  matrix[SICS_n,SICS_mA] SICS_MA;
  vector[SICS_mA] SICS_a;
  // B: inequality constraints
  int<lower=0> SICS_mB;
  matrix[SICS_n,SICS_mB] SICS_MB;
  vector[SICS_mB] SICS_b;
  // C: interval constraints
  int<lower=0> SICS_mC;
  matrix[SICS_n,SICS_mC] SICS_MC;
  vector[SICS_mC] SICS_c0;
  vector[SICS_mC] SICS_c1;
  // U: unconstrained part
  int<lower=0> SICS_mU;
  matrix[SICS_n,SICS_mU] SICS_MU;
  /*^^^ SICS ^^^*/
}

transformed data{
  /*vvv SICS vvv*/
  vector[SICS_n]  SICS_MA_x_a;
  if (SICS_mA > 0) {
    SICS_MA_x_a = SICS_MA * SICS_a;
  }
  /*^^^ SICS ^^^*/
}

parameters {

  vector[n] mu_alpha;
  vector[n] log_sigma_alpha;
  real mu_beta;
  real log_sigma_beta;
  cholesky_factor_corr[n+1] L;

  // array[N] vector[n] alpha; //SICS-ized
  /*vvv SICS vvv*/
  array[SICS_N] vector<lower=rep_array(SICS_b, SICS_N)>[SICS_mB] SICS_wB;
  array[SICS_N] vector<lower=rep_array(SICS_c0, SICS_N), upper=rep_array(SICS_c1,N)>[SICS_mC] SICS_wC;
  array[SICS_N] vector[SICS_mU] SICS_wU;
  /*^^^ SICS ^^^*/
  array[N] real log_beta;
}

transformed parameters {
  array[N] vector[n] alpha;

  /*vvv SICS vvv*/
  for (i in 1:N) {
    vector[n] SICS_wsum = rep_vector(0,SICS_n);
    if (SICS_mA > 0) SICS_wsum += SICS_MA_x_a;
    if (SICS_mB > 0) SICS_wsum += SICS_MB * SICS_wB[i];
    if (SICS_mC > 0) SICS_wsum += SICS_MC * SICS_wC[i];
    if (SICS_mU > 0) SICS_wsum += SICS_MU * SICS_wU[i];
    alpha[i] = SICS_wsum;
  }
  /*^^^ SICS ^^^*/
}

model {

  /*** Population level Priors ***/

   target += normal_lpdf(mu_alpha|hyper_mu_mu_alpha, hyper_sigma_mu_alpha);
   target += normal_lpdf(log_sigma_alpha|hyper_mu_sigma_alpha, hyper_sigma_sigma_alpha);
   target += normal_lpdf(mu_beta|hyper_mu_mu_beta, hyper_sigma_mu_beta);
   target += normal_lpdf(log_sigma_beta|hyper_mu_sigma_beta, hyper_sigma_sigma_beta);
   target += lkj_corr_cholesky_lpdf(L|hyper_eta);

  for (i in 1:N) {
    /*** Participant level distribution ***/
    /* alphas AND beta are sampled from a multivariate normal governed by pop parameters*/
    target += multi_normal_cholesky_lpdf(
      append_row(alpha[i], log_beta[i])|
      append_row(mu_alpha, mu_beta),
      diag_pre_multiply(exp(append_row(log_sigma_alpha,log_sigma_beta)), L)
      );

    if (sample_posterior == 1){
      vector[n] alpha_t = inv_logit(alpha[i]);
      real beta = exp(log_beta[i]);
      real v[2];
      row_vector[T] v_diff;

      for (t in 1:T) {  // Loop over trials

      if (reset[t,i] == 1) {v[1] = 0.5; v[2] = 0.5;} // Reset values at start of new block

      if (reward[t,i] == -9) { // Missing Value handling
      v_diff[t] = 0;           // Missing value is handled like a completely random choice
      // nothing else happens (v0,v1 not updated)
      }
      else {
        v_diff[t] = (v[2] - v[1]);

        // Value update ---
        v[choice[t,i]+1] += alpha_t[k[t,i]] * (reward[t,i] - v[choice[t,i]+1]);
      } //END Missing value handling

      } //END loop over trials
      target += bernoulli_logit_lpmf(choice[:,i]|v_diff * beta);
    }
  }
}
