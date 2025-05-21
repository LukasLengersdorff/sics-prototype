data {
  /*Multinomial model*/
  int<lower=0> n;             // Number of categories
  int<lower=0> y[n];          // Counts per category
  vector<lower=0>[n] alpha;   // Dirichlet parameter alpha
  int sample_posterior;       // Toggle for prior/posterior sampling

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
  /*vvv SICS vvv*/
  vector[SICS_mB] SICS_wB;
  vector[SICS_mC] SICS_wC;
  vector[SICS_mU] SICS_wU;
  /*^^^ SICS ^^^*/
}

transformed parameters {
  vector[n] eta;

  /*vvv SICS vvv*/
  {
    vector[SICS_n] SICS_wsum = rep_vector(0,SICS_n);
    if (SICS_mA > 0) SICS_wsum += SICS_MA_x_a;
    if (SICS_mB > 0) SICS_wsum += SICS_MB * SICS_wB;
    if (SICS_mC > 0) SICS_wsum += SICS_MC * SICS_wC;
    if (SICS_mU > 0) SICS_wsum += SICS_MU * SICS_wU;
    eta = SICS_wsum;
  }
  /*^^^ SICS ^^^*/

}

model {
  vector[n] p = softmax(eta);

  // Prior
  target += dirichlet_lpdf(p|alpha) + sum(log(p)); // sum(log(p)) is the log determinant of the Jacobian of the softmax
  // Likelihood
  if (sample_posterior == 1) {
    target += multinomial_lpmf(y|p);
  }
}

generated quantities{
  vector[n] p = softmax(eta);
}
