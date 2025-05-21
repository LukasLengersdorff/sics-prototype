data {
  int<lower=0> T; // total number of trials
  int choice[T];
  int reward[T];
  int reset[T];
  int<lower=0> k[T];     // number of times the current choice has been chosen, including the current trial
  int<lower=0> n;        // dimension of learning rate parameter >= max(k)
  real exp_lambda;       // Parameter of exponential prior for beta
  real beta_shape[2];    // Shape parameters of beta prior for alpha
  int sample_posterior;  // Toggle for prior/posterior sampling

  /*vvv SICS vvv*/
  int<lower=0> SICS_n;
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
  real log_beta[sample_posterior];

  /*vvv SICS vvv*/
  vector<lower=SICS_b>[SICS_mB] SICS_wB;
  vector<lower=SICS_c0, upper=SICS_c1>[SICS_mC] SICS_wC;
  vector[SICS_mU] SICS_wU;
  /*^^^ SICS ^^^*/
}

transformed parameters {
  vector[n] alpha;

  /*vvv SICS vvv*/
  {
    vector[SICS_n] SICS_wsum = rep_vector(0,SICS_n);
    if (SICS_mA > 0) SICS_wsum += SICS_MA_x_a;
    if (SICS_mB > 0) SICS_wsum += SICS_MB * SICS_wB;
    if (SICS_mC > 0) SICS_wsum += SICS_MC * SICS_wC;
    if (SICS_mU > 0) SICS_wsum += SICS_MU * SICS_wU;
    alpha = SICS_wsum;
  }
  /*^^^ SICS ^^^*/
}

model {
  vector[n] alpha_t = inv_logit(alpha);
  real beta[sample_posterior] = exp(log_beta);
  /*Priors*/
  target += beta_lpdf(alpha_t|beta_shape[1],beta_shape[2]) + logistic_lpdf(alpha|0,1);

  if (sample_posterior == 1){
    vector[T] LL;
    real v[2];
    row_vector[T] v_diff;

    target += exponential_lpdf(beta| exp_lambda) + log_beta[1];

    for (t in 1:T) {  // Loop over trials

    if (reset[t] == 1) {v[1] = 0.5; v[2] = 0.5;} // Reset values at start of new block

    if (reward[t] == -9) { // Missing Value handling
    v_diff[t] = 0;         // Missing value is handled like a completely random choice
    // nothing else happens (v0,v1 not updated)
    }
    else {
      v_diff[t] = (v[2] - v[1]);

      // Value update ---
      v[choice[t]+1] += alpha_t[k[t]] * (reward[t] - v[choice[t]+1]);
    } //END Missing value handling

    } //END loop over trials
    target += bernoulli_logit_lpmf(choice|v_diff * beta[1]);
  }
}

generated quantities {
  vector[n] alpha_t = inv_logit(alpha);
}
