#' Prior distributions
#'
#' @name priors
#'
#' @description These functions are used to specify priors on selected
#'   (hyper)parameters in **salmonIPM** models.
#'
#'   The default priors used in the various models are intended to be weakly
#'   informative, in that they provide moderate regularization and help
#'   stabilize sampling. Priors on scaling parameters, e.g. `Rmax` or `mu_Rmax`,
#'   are automatically adjusted to be weakly informative but consistent with the
#'   observed marginal distribution of population density. For many applications
#'   these defaults will perform well, but if external information not included
#'   in `fish_data` is available, it can be incorporated via user-specified
#'   priors on key parameters. See **Details** for a table of available prior
#'   options.
#'
#' @param mean Prior mean for normal or generalized normal distribution.
#' @param sd Prior standard deviation for normal distribution.
#' @param scale Prior scale for generalized normal distribution. Equivalent to
#'   `alpha` in [gnorm], but renamed to avoid confusion with the spawner-recruit
#'   intrinsic productivity parameter.
#' @param shape Prior shape for generalized normal distribution. Equivalent to
#'   `beta` in [gnorm] but renamed to avoid confusion with covariate slopes.
#' @param meanlog,sdlog Prior log-scale mean and standard deviation,
#'   respectively, for lognormal distribution. See [Lognormal].
#' @param a,b Prior shape parameters for the beta distribution. Equivalent to
#'   `shape1` and `shape2`, respectively, in [Beta].
#' @param concentration Vector of shape parameters for the Dirichlet
#'   distribution. Equivalent to `alpha` in [gtools::dirichlet], but renamed to avoid
#'   confusion with the spawner-recruit intrinsic productivity parameter.
#' @param lb,ub Lower and upper bounds for the uniform distribution.
#' @param eta Prior shape parameter for the [LKJ
#'   distribution](https://mc-stan.org/docs/functions-reference/correlation_matrix_distributions.html)
#'   over correlation matrices.
#'
#' @details The table below shows the parameters in each model that can be given
#'   user-specified priors and the corresponding distributions. Note that users
#'   can modify the prior parameters but not the distribution families;
#'   attempting to do the latter will result in an error.
#'
#'   Priors for parameters that are bounded on the positive real line (e.g.
#'   `tau`, `tau_S` and `tau_M`) are automatically left-truncated at zero.
#'
#'   For parameters that are modeled as functions of covariates using the
#'   `par_models` argument to [salmonIPM()], the specified prior applies when
#'   all predictors are at their sample means.
#'
#'   If `RRS != "none"`, the global spawner-recruit parameters must be replaced
#'   with their `W` and `H` counterparts; e.g. if `RRS == "alpha"` then instead
#'   of a prior on `alpha` one would specify priors on `alpha_W` and `alpha_H`.
#'   If the former is provided, it will have no effect. See [salmonIPM()] for
#'   details of the `RRS` argument.
#'
#'   The generalized normal density with `shape >> 1` is useful as a platykurtic
#'   "soft-uniform" prior to regularize the posterior away from regions of
#'   parameter space that may cause computational or sampling problems. In the
#'   case of spawner and smolt observation error log-SDs, the default prior
#'   bounds them &#8819; 0.1.
#'
#'   The uniform distribution and the LKJ distribution are  included for
#'   internal use; currently no correlation matrices have user-specified priors.
#' 
#' |                    |                         |                         |                     |                        |                        |                        | **Parameter (PDF)**     |                    |                        |                    |                     |                                   |
#' |:-------------------|:-----------------------:|:-----------------------:|:-----------------------:|:------------------:|:----------------------:|:----------------------:|:-----------------------:|:------------------:|:----------------------:|:------------------:|:-------------------:|:---------------------------------:|
#' | **Model**          | `alpha` \cr `lognormal` | `mu_alpha` \cr `normal` | `mu_psi` \cr `beta` | `Rmax` \cr `lognormal` | `mu_Rmax` \cr `normal` | `Mmax` \cr `lognormal` | `mu_Mmax` \cr `normal`  | `mu_MS` \cr `beta` | `mu_p` \cr `dirichlet` | `mu_SS` \cr `beta` | `tau` \cr `gnormal` | `tau_S` \cr `tau_M` \cr `gnormal` |
#' | `IPM_SS_np`        | &#x2611;                | &#x2610;                | &#x2610;            | &#x2611;               | &#x2610;               | &#x2610;               | &#x2610;                | &#x2610;           | &#x2611;               | &#x2610;           | &#x2611;            | &#x2610;                          |
#' | `IPM_SSiter_np`    | &#x2611;                | &#x2610;                | &#x2610;            | &#x2611;               | &#x2610;               | &#x2610;               | &#x2610;                | &#x2610;           | &#x2611;               | &#x2611;           | &#x2611;            | &#x2610;                          |
#' | `IPM_SS_pp`        | &#x2610;                | &#x2611;                | &#x2610;            | &#x2610;               | &#x2611;               | &#x2610;               | &#x2610;                | &#x2610;           | &#x2611;               | &#x2610;           | &#x2611;            | &#x2610;                          |
#' | `IPM_SSiter_pp`    | &#x2610;                | &#x2611;                | &#x2610;            | &#x2610;               | &#x2611;               | &#x2610;               | &#x2610;                | &#x2610;           | &#x2611;               | &#x2611;           | &#x2611;            | &#x2610;                          |
#' | `IPM_SMS_np`       | &#x2611;                | &#x2610;                | &#x2610;            | &#x2610;               | &#x2610;               | &#x2611;               | &#x2610;                | &#x2611;           | &#x2611;               | &#x2610;           | &#x2610;            | &#x2611;                          |
#' | `IPM_SMS_pp`       | &#x2610;                | &#x2611;                | &#x2610;            | &#x2610;               | &#x2610;               | &#x2610;               | &#x2611;                | &#x2611;           | &#x2611;               | &#x2610;           | &#x2610;            | &#x2611;                          |
#' | `IPM_SMaS_np`      | &#x2611;                | &#x2610;                | &#x2610;            | &#x2610;               | &#x2610;               | &#x2611;               | &#x2610;                | &#x2610;           | &#x2610;               | &#x2610;           | &#x2610;            | &#x2611;                          |
#' | `IPM_LCRchum_pp`   | &#x2610;                | &#x2610;                | &#x2611;            | &#x2610;               | &#x2610;               | &#x2610;               | &#x2611;                | &#x2611;           | &#x2611;               | &#x2610;           | &#x2610;            | &#x2610;                          |
#' 
#' @return A named list to be used internally by the **salmonIPM** model-fitting
#'   and summary functions.
#' 
NULL

#' @rdname priors
#' @export
normal <- function(mean = 0, sd = 1) 
{
  stopifnot(sd > 0)
  list(dist = "normal", mean = mean, sd = sd)
}

#' @rdname priors
#' @export
gnormal <- function(mean = 0, scale = 1, shape = 1)
{
  stopifnot(scale > 0 && shape > 0)
  list(dist = "gnormal", mean = mean, scale = scale, shape = shape)
}

#' @rdname priors
#' @export
lognormal <- function(meanlog = 0, sdlog = 1) 
{
  stopifnot(sdlog > 0)
  list(dist = "lognormal", meanlog = meanlog, sdlog = sdlog)
}

#' @rdname priors
#' @export
beta <- function(a = 1, b = 1)
{
  stopifnot(a > 0 && b > 0)
  list(dist = "beta", a = a, b = b)
}

#' @rdname priors
#' @export
dirichlet <- function(concentration = 1) 
{
  stopifnot(all(concentration > 0))
  list(dist = "dirichlet", concentration = concentration)
}

#' @rdname priors
#' @export
uniform <- function(lb = 0, ub = 1) {
  stopifnot(ub > lb)
  list(dist = "uniform", lb = lb, ub = ub)
}

#' @rdname priors
#' @export
lkj_corr <- function(eta = 1) 
{
  stopifnot(eta > 0)
  list(dist = "lkj_corr", eta = eta)
}

