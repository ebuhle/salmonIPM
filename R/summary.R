#' Summary method for `salmonIPMfit` objects
#'
#' Summarize the posterior distributions of estimated parameters and derived
#' quantities using the MCMC draws in a [salmonIPMfit] object.
#'
#' @name summary.salmonIPMfit
#' @method summary salmonIPMfit
#'
#' @param object An object of class [salmonIPMfit].
#' @param pars A character vector specifying (hyper)parameters, states, and/or
#'   quantities of interest ("parameters") to summarize. Parameters can be
#'   explicitly named or one or more shortcuts can be used to specify
#'   hierarchical levels of parameters; see [stan_pars()] for details. The
#'   default is `"hyper"`, i.e. the top-level hyperparameters that are given
#'   priors.
#' @param include Logical scalar defaulting to `TRUE` indicating whether to
#'   include or exclude the parameters given by `pars`. If `FALSE`, only entire
#'   multidimensional parameters can be excluded, rather than particular
#'   elements of them.
#' @param probs A numeric vector of posterior quantiles to return. The default
#'   is `c(0.05, 0.5, 0.95)`, i.e. the median and 90% credible interval.
#' @param funs Named list of summary or diagnostic functions. The provided names
#'   will be used as the names of the columns in the result unless a function
#'   returns a named vector, in which case the latter names are used for the
#'   corresponding columns. The functions can be specified in any format
#'   supported by [rlang::as_function()]. Passed to the `...` argument of
#'   [summarize_draws()].
#' @param .cores The number of cores to use for computing summaries for
#'   different variables in parallel. The default is `.cores = 1`, in which case
#'   no parallelization is used.
#'
#' @return A data frame whose first column contains the variable names and the
#'   remaining columns contain summary statistics and diagnostics.
#'
#' @details Internally the posterior samples are converted to a [draws] object
#'   and [summarize_draws()] is called to compute the summary statistics and
#'   diagnostics. See that function's documentation for available diagnostic
#'   functions, syntax and examples. The default summaries are `mean`, `sd`, and
#'   `~ quantile(.x, probs)`, and the default [diagnostics] are the bulk
#'   effective sample size [ess_bulk()] and the Gelman-Rubin potential scale
#'   reduction factor [rhat()]. If other summary or diagnostic functions are
#'   specified via `funs`, they override the defaults rather than augmenting
#'   them.
#'
#'   If `pars` includes correlation matrices, only the lower triangular elements
#'   are returned. This avoids redundant summary output as well as false
#'   positive diagnostic results such as `bulk-ESS` and `Rhat` being `NaN` for
#'   the diagonal elements.
#'
#' @examples
#' # <under construction>
#'
#' @seealso [draws], [draws_summary()], [diagnostics], [salmonIPMfit]
#' @importFrom posterior summarize_draws ess_bulk rhat
#' @export

summary.salmonIPMfit <- function(object, pars = "hyper", include = TRUE, 
                                 probs = c(0.05, 0.50, 0.95), funs = list(),
                                 .cores = 1) 
{
  default_funs <- list(mean = mean, sd = sd, qnt = ~ quantile(.x, probs = probs), 
                       ess_bulk = ess_bulk, rhat = rhat)
  funs <- if(length(funs)) funs else default_funs
  drf <- as_draws_df(object, pars = pars, include = include)
  drf <- draws_df_lower_tri(drf)  # only lower triangle of correlation matrices
  sdf <- do.call(summarize_draws, list(.x = drf, funs, .cores = .cores))
  return(as.data.frame(sdf))
}
