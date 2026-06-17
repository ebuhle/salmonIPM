#' Spawner-recruit functions
#' 
#' Compute recruitment given spawner abundance, a spawner-recruit function and parameters.
#' 
#' @param alpha Numeric vector, matrix, data frame, [posterior::rvar()],
#'   or `draws` of intrinsic productivity (i.e., recruits per spawner at zero spawner 
#'   density; slope of the spawner-recruit function at the origin).
#' @param alpha_W,alpha_H Numeric vectors, matrices, data frames, [posterior::rvar()] 
#'   or `draws` of intrinsic productivity for wild- and hatchery-origin spawners, 
#'   respectively. If they differ, both must be specified and `alpha` must not be used
#'   (and conversely). 
#' @param Rmax Numeric vector, matrix, data frame, [posterior::rvar()], 
#'   or `draws` of maximum recruitment per unit of habitat (length or area). This corresponds to 
#'   the asymptote of the Beverton-Holt or the mode of the Ricker.
#' @param Rmax_W,Rmax_H Numeric vectors, matrices, data frames, [posterior::rvar()] 
#'   or `draws` of maximum recruitment per unit of habitat for wild- and hatchery-origin spawners, 
#'   respectively. If they differ, both must be specified and `Rmax` must not be used
#'   (and conversely).
#' @param S Numeric vector, matrix, data frame, [posterior::rvar()], 
#'   or `draws` of spawner abundance.
#' @param S_W,S_H Numeric vectors, matrices, data frames, [posterior::rvar()] 
#'   or `draws` of wild- and hatchery-origin spawner abundance, respectively. 
#'   Must be specified if either `alpha` or `Rmax` differ by rearing type, 
#'   in which case `S` must not be used (and conversely).
#' @param A Numeric vector, matrix, data frame, [posterior::rvar()], or `draws` of 
#'   spawning habitat size (either stream length or area), used to standardize `Rmax`. 
#'   The default is 1, in which case `Rmax` is in units of abundance (which is also density).
#' @param R_per_S Logical indicating whether to return recruits per spawner rather than
#'   recruits (the default).
#' @inheritParams salmonIPM
#'   
#' @return A vector, matrix, data frame, [posterior::rvar()] or `draws`, depending on the 
#'   argument types, containing either recruits or recruits per spawner. Calculations are 
#'   vectorized and elements of shorter arguments are recycled as necessary.
#'   
#' @details The **salmonIPM** package uses a nonstandard parameterization of the Ricker
#'   model by the maximum recruitment `Rmax`. This is typically better identified by
#'   data than per capita density dependence, and it facilitates a common interpretation 
#'   and priors with the Beverton-Holt. Here \eqn{e} is the base of the natural logarithm.
#'   
#'   Note that the functions for the `RRS != "none"` case are written below in their most general
#'   form, with both `alpha` and `Rmax` differing between wild and hatchery spawners. If
#'   only one parameter is specified in `RRS`, then the `_W` and `_H` values of the other
#'   parameter are equal and the expression can be further simplified.
#'
#' `RRS == "none"`    
#' 
#' \eqn{
#' R = 
#' \begin{cases}
#' \alpha S & \text{exponential} 
#' \\\\
#' \dfrac{\alpha S}{1 + \dfrac{\alpha S}{A R_\text{max}}} & \text{Beverton-Holt} 
#' \\\\
#' \alpha S \text{exp} {\left(- \dfrac{\alpha S}{e A R_\text{max}} \right)} & \text{Ricker}
#' \end{cases}
#' }
#' 
#' `RRS != "none"`
#' 
#' \eqn{
#' R = 
#' \begin{cases}
#' \alpha_\text{W} S_\text{W} + \alpha_\text{H} S_\text{H} & \text{exponential} 
#' \\\\
#' \dfrac{\alpha_\text{W} S_\text{W} + \alpha_\text{H} S_\text{H}}{1 + \dfrac{\alpha_\text{W} S_\text{W}}{A R_\text{max,W}} + \dfrac{\alpha_\text{H} S_\text{H}}{A R_\text{max,H}}} & \text{Beverton-Holt (Leslie-Gower)} 
#' \\\\
#' \left(\alpha_\text{W} S_\text{W} + \alpha_\text{H} S_\text{H} \right)  \text{exp}\left(-\dfrac{\alpha_\text{W} S_\text{W}}{e A R_\text{max,W}} - \dfrac{\alpha_\text{H} S_\text{H}}{e A R_\text{max,H}} \right) & \text{Ricker}
#' \end{cases}
#' }
#' 
#' @examples 
#' alpha <- 3
#' Rmax <- 1000
#' S <- 500
#' 
#' # default is Beverton-Holt
#' SR(alpha = alpha, Rmax = Rmax, S = S) 
#'
#' # approximately Rmax
#' SR(alpha = alpha, Rmax = Rmax, S = 1e6) 
#'
#' # scale Rmax by habitat area
#' SR(alpha = alpha, Rmax = Rmax, S = S, A = 0.1) 
#'
#' # discrete exponential ignores Rmax
#' SR(SR_fun = "exp", alpha = alpha, Rmax = Rmax, S = S) 
#'
#' # vectorization with recycling
#' SR(alpha = rep(alpha, 10), Rmax = rep(Rmax, 10), S = matrix(S, 10, 4)) 
#'
#' # return recruits per spawner
#' SR(alpha = alpha, Rmax = Rmax, S = S, R_per_S = TRUE) 
#'
#' # plot Ricker and show Rmax
#' curve(SR(SR_fun = "Ricker", alpha = alpha, Rmax = Rmax, S = x), from = 0, to = 2000,
#'       xlab = "Spawners", ylab = "Recruits", main = "Ricker")
#' abline(h = Rmax, lty = 2)
#'
#' # differential hatchery / wild relative reproductive success
#' alpha_W <- 3
#' alpha_H <- 2
#' S_W <- 400
#' S_H <- 100
#' # compare to BH result above
#' SR(alpha_W = alpha_W, alpha_H = alpha_H, Rmax = Rmax, S_W = S_W, S_H = S_H)
#'
#' @seealso [salmonIPM()] for fitting models, [simIPM()] for simulating data
#' 
#' @importFrom posterior is_rvar
#' @export

SR <- function(SR_fun = c("BH","B-H","bh","b-h","Ricker","ricker","exp"),
               alpha = NULL, alpha_W = NULL, alpha_H = NULL, 
               Rmax = NULL, Rmax_W = NULL, Rmax_H = NULL,
               S = NULL, S_W = NULL, S_H = NULL, A = 1, R_per_S = FALSE) 
{
  SR_fun <- match.arg(SR_fun)
  if(SR_fun == "DI") SR_fun <- "exp"
  if(SR_fun %in% c("B-H","bh","b-h")) SR_fun <- "BH"
  if(SR_fun == "ricker") SR_fun <- "Ricker"
  
  if((!is.null(alpha) & (!is.null(alpha_W) | !is.null(alpha_H))) | 
     is.null(alpha_W) != is.null(alpha_H))
    stop("If alpha_W and alpha_H differ, both must be specified and `alpha`",
         "\n  must not be used (and conversely)")
  if((!is.null(Rmax) & (!is.null(Rmax_W) | !is.null(Rmax_H))) | 
     is.null(Rmax_W) != is.null(Rmax_H))
    stop("If Rmax_W and Rmax_H differ, both must be specified and `Rmax`",
         "\n  must not be used (and conversely)")
  if((!is.null(alpha_W) | !is.null(alpha_H) | !is.null(Rmax_W) | !is.null(Rmax_H)) & 
     (is.null(S_W) | is.null(S_H)))
    stop("S_W and S_H must be specified if either `alpha` or `Rmax` differ by rearing type",
         "\n  in which case `S` must not be used (and conversely)")

  # Use most general form of unequal-RRS model to avoid repetition
  if(!is.null(alpha_W) & is.null(Rmax_W)) {
    Rmax_W <- Rmax
    Rmax_H <- Rmax
  }
  if(is.null(alpha_W) & !is.null(Rmax_W)) {
    alpha_W <- alpha
    alpha_H <- alpha
  }
  
  for(i in c("alpha","alpha_W","alpha_H","Rmax","Rmax_W","Rmax_H","S","S_W","S_H","A")) {
    if(is.data.frame(get(i))) assign(i, as.matrix(get(i)))
  }

  if(is.null(alpha_W) & is.null(Rmax_W)) { # relative reproductive success == 1
    R <- switch(
      SR_fun,
      exp = alpha*S,
      BH = alpha*S / (1 + alpha*S/(A*Rmax)),
      Ricker = alpha*S*exp(-alpha*S/(exp(1)*A*Rmax))
    ) 
  } else {  # relative reproductive success != 1
    S <- S_W + S_H
    R <- switch(
      SR_fun,
      exp = alpha_W*S_W + alpha_H*S_H,
      BH = (alpha_W*S_W + alpha_H*S_H) / (1 + alpha_W*S_W/(A*Rmax_W) + alpha_H*S_H/(A*Rmax_H)),
      Ricker = (alpha_W*S_W + alpha_H*S_H) * exp(-alpha_W*S_W/(exp(1)*A*Rmax_W) - alpha_H*S_H/(exp(1)*A*Rmax_H))
    )
  }
  
  return(if(R_per_S) R/S else R) 
}

