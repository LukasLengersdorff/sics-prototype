## Multinomial model ----

#' SICS Frame - Multinomial model
#'
#' Collect data and prior information for a multinomial model
#'
#' @param y Vector of counts
#' @param alpha ...
#' @param alpha_default ...
#' @returns A SICSFrame_Multinomial object
#'
#' @export
SICSFrame_Multinomial = function(y,
                                 alpha = NULL,
                                 alpha_default = 1) {
  n = length(y)
  if (is.null(alpha)) {
    alpha = rep(alpha_default, n)
  } else if (length(alpha) == 1) {
    alpha = rep(alpha, n)
  } else if (length(alpha) != n) {
    stop("Length of alpha is not equal to length of y")
  }

  structure(list(y = y,
                 n = n,
                 alpha = alpha,
                 SICS_n = n),
            class = c("SICSFrame_Multinomial", "SICSFrame"))
}

#'
#'
#' @export
model.SICSFrame_Multinomial = function(SF,
                                       constraint = NULL,
                                          add_A0 = TRUE,
                                          index_A0 = 1,
                                       fit_prior = TRUE,
                                       fit_posterior = TRUE,
                                       no_fit = FALSE,
                                       alternative_stanmodel = NULL, # NULL to use precompiled STAN model. For experimental use of alternative stanmodel, pass it here.
                                       ...) {


  if (add_A0) {
    A0 = ConstantConstraint(indices = index_A0, constant = 0)
    constraint = constraint + A0
  }
  if (length(constraint$a) != constraint$mA) {constraint$a = numeric(constraint$mA)}
  # Solves a bug that happens when a = NULL but A is not empty,
  # which should normally result in the default of a = numeric(mA),
  # but doesn't work with the current A0 solution

  no_params = ifelse(constraint$mA == constraint$n, "both", FALSE)

  if (is.null(alternative_stanmodel)) {
    stanmod = stanmodels$sics_multinomial
  } else {
    stanmod = alternative_stanmodel
  }

  Modelobj = model.SICSFrame(SF = SF,
                             constraint = constraint,
                             fit_prior = fit_prior,
                             fit_posterior = fit_posterior,
                             no_fit = no_fit,
                             no_params = no_params,
                             stanmod = stanmod,
                             ...)
  class(Modelobj) = c("SICSModel_Multinomial", "SICSModel")
  return(Modelobj)
}


## RL ----

#'
#'
#' @export
SICSFrame_RL = function(choice,
                        reward,
                        reset,
                        k,
                        max_k = max(k),
                        exp_lambda = 1,
                        beta_shape1 = 1,
                        beta_shape2 = 1) {
  T = length(choice)
  k_base = k
  k = pmin(k, max_k)
  structure(list(T = T,
                 choice = choice,
                 reward = reward,
                 reset = reset,
                 k = k,
                 k_base = k_base,
                 n = max_k,
                 SICS_n = max_k,
                 exp_lambda = exp_lambda,
                 beta_shape = c(beta_shape1, beta_shape2)),
            class = c("SICSFrame_RL","SICSFrame"))
}

#'
#'
#' @export
model.SICSFrame_RL = function(SF,
                              constraint = NULL,
                              fit_prior = TRUE,
                              fit_posterior = TRUE,
                              no_fit = FALSE,
                              alternative_stanmodel = NULL,
                              ...) {

  no_params = ifelse(constraint$mA == constraint$n, "prior", FALSE)

  if (is.null(alternative_stanmodel)) {
    stanmod = stanmodels$sics_rl
  } else {
    stanmod = alternative_stanmodel
  }

  Modelobj = model.SICSFrame(SF = SF,
                             constraint = constraint,
                             fit_prior = fit_prior,
                             fit_posterior = fit_posterior,
                             no_fit = no_fit,
                             no_params = no_params,
                             stanmod = stanmod,
                             ...)
  class(Modelobj) = c("SICSModel_RL", "SICSModel")
  return(Modelobj)
}

## RL multilevel 1: constraint on population parameters----

#'
#'
#' @export
SICSFrame_RL_ml1 = function(choice,
                            reward,
                            reset,
                            k,
                            max_k = max(k),
                            hyperprior_mu_alpha = c(0,2),
                            hyperprior_sigma_alpha = c(1,1),
                            hyperprior_mu_beta = c(0,1),
                            hyperprior_sigma_beta = c(1,1),
                            hyperprior_lkj_eta = 1) {
  T = nrow(choice)
  N = ncol(choice)
  k_base = k
  k = pmin(k, max_k)
  structure(list(N = N,
                 T = T,
                 choice = choice,
                 reward = reward,
                 reset = reset,
                 k = k,
                 k_base = k_base,
                 n = max_k,
                 SICS_n = max_k,
                 hyper_mu_mu_alpha = hyperprior_mu_alpha[1],
                 hyper_sigma_mu_alpha = hyperprior_mu_alpha[2],
                 hyper_mu_sigma_alpha = hyperprior_sigma_alpha[1],
                 hyper_sigma_sigma_alpha = hyperprior_sigma_alpha[2],
                 hyper_mu_mu_beta = hyperprior_mu_beta[1],
                 hyper_sigma_mu_beta = hyperprior_mu_beta[2],
                 hyper_mu_sigma_beta = hyperprior_sigma_beta[1],
                 hyper_sigma_sigma_beta = hyperprior_sigma_beta[2],
                 hyper_eta = hyperprior_lkj_eta),
            class = c("SICSFrame_RL_ml1", "SICSFrame_RL_ml", "SICSFrame"))
}

#'
#'
#' @export
#'
model.SICSFrame_RL_ml1 = function(SF,
                                  constraint = NULL,
                                  fit_prior = TRUE,
                                  fit_posterior = TRUE,
                                  no_fit = FALSE,
                                  alternative_stanmodel = NULL,
                                  ...) {

  ####################
  no_params = FALSE ##
  ####################
  # TBD: The Stan model could be made more efficient by removing all non-SICSized parameters (mu_beta, log_sigma_beta, L...)
  # if sample_posterior = 0. Then we would have
  # no_params = ifelse(constraint$mA == constraint$n, "prior", FALSE)
  ####################


  if (is.null(alternative_stanmodel)) {
    stanmod = stanmodels$sics_rl_ml1
  } else {
    stanmod = alternative_stanmodel
  }

  Modelobj = model.SICSFrame(SF = SF,
                             constraint = constraint,
                             fit_prior = fit_prior,
                             fit_posterior = fit_posterior,
                             no_fit = no_fit,
                             no_params = no_params,
                             stanmod = stanmod,
                             ...)
  class(Modelobj) = c("SICSModel_RL_ml1", "SICSModel_RL_ml", "SICSModel")
  return(Modelobj)
}


## RL multilevel 2: constraint on every subject parameter----

#'
#'
#' @export
SICSFrame_RL_ml2 = function(choice,
                            reward,
                            reset,
                            k,
                            max_k = max(k),
                            hyperprior_mu_alpha = c(0,2),
                            hyperprior_sigma_alpha = c(1,1),
                            hyperprior_mu_beta = c(0,1),
                            hyperprior_sigma_beta = c(1,1),
                            hyperprior_lkj_eta = 1) {
  T = nrow(choice)
  N = ncol(choice)
  k_base = k
  k = pmin(k, max_k)
  structure(list(N = N,
                 SICS_N = N,
                 T = T,
                 choice = choice,
                 reward = reward,
                 reset = reset,
                 k = k,
                 k_base = k_base,
                 n = max_k,
                 SICS_n = max_k,
                 hyper_mu_mu_alpha = hyperprior_mu_alpha[1],
                 hyper_sigma_mu_alpha = hyperprior_mu_alpha[2],
                 hyper_mu_sigma_alpha = hyperprior_sigma_alpha[1],
                 hyper_sigma_sigma_alpha = hyperprior_sigma_alpha[2],
                 hyper_mu_mu_beta = hyperprior_mu_beta[1],
                 hyper_sigma_mu_beta = hyperprior_mu_beta[2],
                 hyper_mu_sigma_beta = hyperprior_sigma_beta[1],
                 hyper_sigma_sigma_beta = hyperprior_sigma_beta[2],
                 hyper_eta = hyperprior_lkj_eta),
            class = c("SICSFrame_RL_ml2", "SICSFrame_RL_ml", "SICSFrame"))
}

#'
#'
#' @export
model.SICSFrame_RL_ml2 = function(SF,
                                  constraint = NULL,
                                  fit_prior = TRUE,
                                  fit_posterior = TRUE,
                                  no_fit = FALSE,
                                  alternative_stanmodel = NULL,
                                  ...) {

  fit0 = NULL
  fit1 = NULL

  ####################
  no_params = FALSE ##
  ####################
  # In RL_ml2, there are always parameters to be sampled

  if (is.null(alternative_stanmodel)) {
    stanmod = stanmodels$sics_rl_ml2
  } else {
    stanmod = alternative_stanmodel
  }

  Modelobj = model.SICSFrame(SF = SF,
                             constraint = constraint,
                             fit_prior = fit_prior,
                             fit_posterior = fit_posterior,
                             no_fit = no_fit,
                             no_params = no_params,
                             stanmod = stanmod,
                             ...)
  class(Modelobj) = c("SICSModel_RL_ml2", "SICSModel_RL_ml", "SICSModel")
  return(Modelobj)
}
