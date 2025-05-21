### --------- Generic functions --------

#'
#'
#' @export
#'
#'
model <- function (x, ...) {
  UseMethod("model", x)
}

#'
#'
#' @export
#'
#'
fit <- function (x, ...) {
  UseMethod("fit", x)
}

#'
#'
#' @export
#'
#'
logml <- function (x, ...) {
  UseMethod("logml", x)
}



### --------- Common methods ----------

#' Prepare and fit a SICS model (SICSFrame + Constraint)
#'
#' This method is called internally by all submethods (model.SICSFrame_*).
#'
#' Use case for manual use: Development of new SICS Models without overhead.
#' @export
model.SICSFrame = function(SF,
                           constraint = NULL,
                           fit_prior = TRUE,
                           fit_posterior = TRUE,
                           no_fit = FALSE,
                           no_params = FALSE,
                           stanmod = stop("model.SICSFrame called without specifying stanmod"),
                           ...) {

  # The following things should happen in the subclass method model.SICSFrame_*:
  # - Determination of no_params
  # - Necessary additions or changes to the constraint (e.g., obligatory equality constraints)
  # - Specification of correct stan_model object stanmod

  standata = c(SF,
               list(SICS_mA = constraint$mA,
                    SICS_mB = constraint$mB,
                    SICS_mC = constraint$mC,
                    SICS_mU = constraint$mU,
                    SICS_MA = constraint$MA,
                    SICS_MB = constraint$MB,
                    SICS_MC = constraint$MC,
                    SICS_MU = constraint$MU,
                    SICS_a  = as.array(constraint$a),
                    SICS_b  = as.array(constraint$b),
                    SICS_c0 = as.array(constraint$c0),
                    SICS_c1 = as.array(constraint$c1)))

  Modelobj = structure(list(fit0 = NULL, fit1 = NULL, constraint = constraint, standata = standata, no_params = no_params, stanmod = stanmod), class = "SICSModel")
  if (!no_fit) {
    Modelobj = fit.SICSModel(Modelobj, fit_prior = fit_prior, fit_posterior = fit_posterior, stanmod = stanmod, ...)
  }
  return(Modelobj)
}

#' Fits a SICS Model (SICSFrame + Constraint)
#'
#' In most cases, you will NOT need to manually call this function,
#' as it is called internally by all model.SICSFrame_* methods.
#' @export
fit.SICSModel = function(Modelobj,
                         fit_prior = TRUE,
                         fit_posterior = TRUE,
                         iter_posterior = 2000,
                         chains_posterior = 4,
                         iter_prior = 5*iter_posterior,
                         chains_prior = 4,
                         force_prior_sampling = FALSE,
                         stanmod = Modelobj$stanmod, ...) {

  if (fit_prior) {
    if (Modelobj$no_params %in% c("prior","both")) {
      .algo = "Fixed_param"
      .chains = 1
      .iter = 1
    } else {
      .algo = "NUTS"
      .chains = chains_prior
      .iter = iter_prior
    }

    if (!attr(Modelobj$constraint, "empty") || force_prior_sampling) {
      Modelobj$fit0 = rstan::sampling(stanmod, data = c(Modelobj$standata, sample_posterior = 0),
                               algorithm = .algo,
                               chains = .chains,
                               iter = .iter,
                               ...)
    } else {
      Modelobj$fit0 = "Skipped"
      message("\nSkipped sampling from the prior because constraint is empty.\n")
    }
  }
  if (fit_posterior) {
    if (Modelobj$no_params %in% c("both")) {
      .algo = "Fixed_param"
      .chains = 1
      .iter = 1
    } else {
      .algo = "NUTS"
      .chains = chains_posterior
      .iter = iter_posterior
      Modelobj$fit1 = rstan::sampling(stanmod, data = c(Modelobj$standata, sample_posterior = 1),
                               algorithm = .algo,
                               chains = .chains,
                               iter = .iter,
                               ...)
    }
  }
  Modelobj
}


#'
#'
#' @export
logml.SICSModel = function(modelobj, add_logdet_jacobian = TRUE, silent = TRUE,  ...) {
  if (is.null(modelobj$fit0) || is.null(modelobj$fit1)) stop("Both prior and posterior need to be fitted to compute Log Marginal Likelihood")
  if (modelobj$no_params == FALSE) {
    ml0 = ifelse(attr(modelobj$constraint, "empty"),
                 0,
                 bridgesampling::bridge_sampler(modelobj$fit0, silent = silent, ...)$logml)
    ml1 = bridgesampling::bridge_sampler(modelobj$fit1, silent = silent, ...)$logml
  } else {
    ml0 = rstan::log_prob(object = modelobj$fit0, upars = numeric(0))
    ml1 = ifelse(modelobj$no_params == "both",
                 rstan::log_prob(object = modelobj$fit1, upars = numeric(0)),
                 bridgesampling::bridge_sampler(modelobj$fit1, silent = silent, ...)$logml)
  }
  if (add_logdet_jacobian) {
    N = modelobj$standata$SICS_N
    if (is.null(N)) N = 1
    D = compute_logdet_jacobian(modelobj$constraint)
    ml0 = ml0 + N*D
    ml1 = ml1 + N*D
  }
  return(c("logml" = ml1-ml0, "logpC" = ml0, "logpCD" = ml1))
}
