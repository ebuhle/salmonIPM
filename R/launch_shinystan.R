#' Launch the ShinyStan GUI with `salmonIPM` models
#'
#' The ShinyStan interface provides visual and numerical summaries of 
#' model parameters and convergence diagnostics.
#'
#' @aliases launch_shinystan
#' @param object A [salmonIPMfit] object.
#' @param ... Additional arguments passed to [launch_shinystan.default].
#' @seealso [shinystan::launch_shinystan()], [salmonIPM()]
#' @importFrom shinystan launch_shinystan
#' @exportS3Method shinystan::launch_shinystan
#' @export launch_shinystan
launch_shinystan.salmonIPMfit <- function(object, ...) {
  launch_shinystan(object$stanfit, ...)
}