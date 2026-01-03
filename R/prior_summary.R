#' Summarize the priors used in a `salmonIPM` model
#'
#' Print a summary of the priors on hyperparameters in a fitted [salmonIPMfit]
#' model object.
#'
#' @name prior_summary.salmonIPMfit
#' @aliases prior_summary
#'
#' @param object An object of class [salmonIPMfit] with prior information stored
#'   in `prior.info`.
#' @param digits Number of decimal places to print.
#' @param ... Currently ignored.
#'
#' @return A character vector containing the prior summary information.
#'
#' @details This function is called mainly for its side effect of printing the
#'   prior summary. The character vector to be printed is returned invisibly,
#'   but the same information can be found in more usable form in
#'   `object$prior.info`. Priors are identified as user-specifiable or
#'   hard-coded in the printed output.
#'
#'   Note that correlation matrices are given
#'   [`lkj_corr`](https://mc-stan.org/docs/functions-reference/correlation_matrix_distributions.html)
#'   priors, but the actual Stan implementation uses the more efficient lower
#'   Cholesky factor version of the LKJ distribution.
#'
#' @seealso [prior_summary()], [priors], [salmonIPM()], [salmonIPMfit]
#'
#' @importFrom rstantools prior_summary
#' @exportS3Method rstantools::prior_summary
#' @export prior_summary

prior_summary.salmonIPMfit <- function(object, digits = 2, ...) {
  prior.info <- object$prior.info
  pars <- names(prior.info)
  types <- sapply(prior.info, attr, "type")
  user_priors <- as.list(object$call)$priors
  user_pars <- pars %in% sapply(user_priors, function(f) all.vars(f)[1])
  user <- ifelse(types == "user", paste0("[", ifelse(user_pars, "x", " "), "]  "), "     ")
  names(user) <- pars
  pad <- sapply(max(nchar(pars)) - nchar(pars), function(n) paste(rep(" ", n), collapse = ""))
  names(pad) <- pars
  
  prsum <- sapply(pars, function(.par) {
    prinfo <- prior.info[[.par]]
    dist <- prinfo$dist
    args <- lapply(prinfo[setdiff(names(prinfo), "dist")], round, digits = digits)
    args <- paste(paste0(names(args), " = ", args), collapse = ", ")
    bounds <- attr(prinfo, "bounds")
    if(!is.null(bounds)) 
      bounds <- paste0(" T[", gsub("NA", "", paste0(bounds, collapse = ", ")), "]")
    paste0(user[.par], pad[.par], .par, " ~ ", dist, "(", args, ")", bounds, "\n") 
  })
  prsum <- unlist(prsum)
  
  cat("Priors for model ", object$stan_model, "\n", 
      "[ ] user-specifiable (default) \n",
      "[x] set by user \n\n",
      prsum,
      sep = "")
  
  invisible(prsum)
}

#' Extract prior info from Stan-formatted data and model code
#'
#' Helper function to extract both user-specifiable and hard-coded prior data
#' and assemble it in the appropriate form for the `prior.info` element of a
#' [salmonIPMfit] object.
#'
#' @param stan_data Named list of input data passed to [rstan::sampling()] for
#'   fitting an IPM, as returned by [stan_data()].
#' @param stanmodel An object of S4 class `stanmodel` as constructed by
#'   [stan_model()] or by **rstantools** during installation of
#'   [salmonIPM-package] and stored internally in `salmonIPM:::stanmodels`.
#' @param pars Character vector of hyperparameters for which priors are
#'   specified.
#'
#' @return A named list of priors on `pars`. Each element is itself a named list
#'   representing the prior in the format returned by the functions in [priors]
#'   with an additional attribute `"type"`.
#'
#' @details Priors include (1) Those that are user-specifiable through the
#'   `salmonIPM(priors)` argument (whether modified or defaults) and returned in
#'   Stan format by [stan_data()], (2) those that are explicitly hard-coded in
#'   the `model` block of `stanmodel`, and (3) those that are implicitly defined
#'   by bounds on `parameter` declarations.
#'
#'   To extract type (1), `get_prior_info()` calls the [priors] functions using
#'   the contents of `stan_data`.
#'
#'   To extract type (2), `get_prior_info()` parses `stanmodel` into lines of
#'   text and then uses regex to pull out the sampling statement for each
#'   parameter as a string, which corresponds directly to a [priors] function
#'   call.
#'
#'   Type (3) is extracted similarly to type (2), but bound declarations are
#'   also used to set the `bounds` attribute on otherwise unbounded priors (e.g.
#'   normal) of types (1) and (2). A special case is simplex parameters,
#'   which unless otherwise specified are assigned a `dirichlet(1)` prior
#'   (ignoring the actual dimension for now).
#'
#' @importFrom utils relist

get_prior_info <- function(stan_data, stanmodel, pars) 
{
  pars <- pars[!grepl("delta_.*alpha|delta_.*max", pars)]  # H/W discounts don't have priors
  code_pars <- gsub("R_", "L_", pars)  # correlation matrices -> Cholesky factors
  pars <- setNames(code_pars, pars)
  smtext <- strsplit(stanmodel@model_code, "\\n")[[1]]
  
  # User-specifiable priors from stan_data
  prior_data <- stan_data[names(stan_data) %in% paste0("prior_", pars)]
  names(prior_data) <- gsub("prior_", "", names(prior_data))
  user_priors <- lapply(prior_data, function(x) {
    structure(do.call(attr(x, "dist"), relist(x)), type = "user")
  })
  
  # Hard-coded priors from stanmodel
  #
  # .*      any character, 0 or more (e.g. indent spaces, any transformation function)
  # .par    hyperparameter name
  # .*      any character, 0 or more (e.g. closing paren of transformation, space) 
  # ~ *     this line contains a sampling statement (may be 0 or more spaces after "~")
  # (.+);   pdf capture group: any character, 1 or more, then end of sampling statement
  # .*      any character, 0 or more (e.g. comment)
  hard_pars <- pars[!pars %in% names(user_priors)] 
  hard_priors <- lapply(hard_pars, function(.par) {
    regex <- paste0(".*", .par, ".*~ *(.+);.*")
    prtext <- gsub(regex, "\\1", grep(regex, smtext, value = TRUE))
    prtext <- gsub("_cholesky", "", prtext)  # Cholesky factors -> correlation matrices
    if(length(prtext)) structure(as.list(eval(parse(text = prtext))), type = "hard")
  })
  hard_priors <- hard_priors[!sapply(hard_priors, is.null)]

  # Implicit priors declared by bounds in stanmodel
  #
  # .+          any character, 1 or more (e.g. indent spaces, variable type, <)
  # lower *= *  beginning of bounds (may or may not have spaces around "=")
  # (-?\\d+)    lb capture group: optional negative sign, then 1 or more digits 
  # upper *= *  upper bound declaration
  # (-?\\d+)    ub capture group 
  # .*          any character, 0 or more (e.g. >, dims)
  #  .par       space before hyperparameter name
  # ;.*         end of declaration followed by 0 or more characters (e.g. comment)
  bounded_priors <- sapply(pars, as.null, 0)
  for(.par in pars) {
    if(grepl("rho_alpha.max", .par)) {  # hack b/c declared parameter is L_alpha[R|M]max
      bounded_priors[[.par]] <- structure(uniform(-1,1), type = "bounded")
    } else {
      lregex <- paste0(".+lower *= *(-?\\d+).* ", .par, ";.*")
      lb <- gsub(lregex, "\\1", grep(lregex, smtext, value = TRUE))
      lb <- ifelse(length(lb), as.numeric(lb), NA)
      uregex <- paste0(".+upper *= *(-?\\d+).* ", .par, ";.*")
      ub <- gsub(uregex, "\\1", grep(uregex, smtext, value = TRUE))
      ub <- ifelse(length(ub), as.numeric(ub), NA)
      bounds <- c(lb = lb, ub = ub)
      not_NA <- sum(!is.na(bounds))
      if(!is.null(user_priors[[.par]]) && not_NA && user_priors[[.par]]$dist != "lognormal") {
        attr(user_priors[[.par]], "bounds") <- bounds
      } else if(!is.null(hard_priors[[.par]]) && not_NA && hard_priors[[.par]]$dist != "lognormal") {
        attr(hard_priors[[.par]], "bounds") <- bounds
      } else if(not_NA == 2) {
        bounded_priors[[.par]] <- structure(do.call(uniform, as.list(bounds)), type = "bounded")
      }
    }
  }
  bounded_priors <- bounded_priors[!sapply(bounded_priors, is.null)]
  
  # Implicit bounded priors declared by parameter type: simplex
  #
  # .+        any character, 1 or more (e.g. indent spaces, array declaration) 
  # simplex   parameter type
  # \\[.+\\]  simplex dimension (don't use this info for now)
  #  .par     space before hyperparameter name
  # ;.*       end of declaration followed by 0 or more characters (e.g. comment)
  simplex_pars <- pars[!pars %in% c(names(user_priors), names(hard_priors))] 
  simplex_priors <- lapply(simplex_pars, function(.par) {
    sregex <- paste0(".+simplex\\[.+\\] ", .par, ";.*")
    smatch <- any(grepl(sregex, smtext))
    if(smatch) structure(dirichlet(1), type = "bounded")
  })
  simplex_priors <- simplex_priors[!sapply(simplex_priors, is.null)]
  
  all_priors <- c(user_priors, hard_priors, bounded_priors, simplex_priors)
  return(all_priors[match(names(pars), names(all_priors))])
}
