#' Validate prior specification for an integrated spawner-recruit model.
#' 
#' Takes an optional user-specified prior function and otherwise returns default
#' prior parameters for the given `stan_model`. This function is for internal use. 
#' 
#' @param user A user-specified prior function as described in [`priors`].
#' @param default The default prior function as described in [`priors`].
#'
#' @return A named array of user-specifiable prior parameters whose elements are part of
#' the `data` argument passed to [rstan::sampling()] when fitting **salmonIPM** models.
#' The `dist` attribute is a character string giving the name of the distribution.
#' 
#' @importFrom utils as.relistable
stan_prior <- function(user, default)
{
  pd <- eval(default)
  if(is.null(user)) {
    p <- pd
  } else {
    p <- eval(user)
    stopifnot(all(names(p) == names(pd)))
    stopifnot(user$dist == pd$dist)
  }
  
  pdist <- p$dist
  ppars <- p[setdiff(names(p), "dist")]
  pout <- unlist(as.relistable(ppars))
  return(structure(as.array(pout), dist = pdist))
}
