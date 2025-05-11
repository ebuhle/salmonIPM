#' Return the (hyper)parameters and states in a specified **salmonIPM** model
#'
#' @param stan_model Character string giving the name of the model. See
#'   [salmonIPM()] for details.
#' @param pars An optional character vector specifying one or more hierarchical
#'   levels of parameters. Options are `"all"` (the default), `"hyper"`
#'   (top-level hyperparameters that are given priors), `"group"` (`pop`- or
#'   `year`-level parameters shared by multiple states), `"states"` (the lowest
#'   level, corresponding to unique rows in `fish_data`), and `"ppd"` (only if
#'   `model == "RR"`, observation-level predictions drawn from the posterior
#'   predictive distribution).
#' @param object A [salmonIPMfit] object. If this is provided then `SR_fun`,
#'   `RRS` and `par_models` are not needed and will be ignored; their values
#'   are extracted from `object`.
#' @inheritParams salmonIPM
#'
#' @return Character vector with names of selected parameters and states
#'
#' @export

stan_pars <- function(stan_model = c("IPM_SS_np","IPM_SSiter_np","IPM_SS_pp","IPM_SSiter_pp",
                                     "IPM_SMS_np","IPM_SMS_pp","IPM_SMaS_np",
                                     "IPM_LCRchum_pp","RR_SS_np","RR_SS_pp"), 
                      pars = c("all","hyper","group","states","ppd"), 
                      SR_fun = "BH", RRS = "none", par_models = NULL, object = NULL)
{
  if(!is.null(object)) {
    stopifnot("salmonIPMfit" %in% class(object))
    stan_model <- object$stan_model
    SR_fun <- object$SR_fun
    RRS <- object$RRS
    par_models <- object$par_models
  } else {
    stan_model <- match.arg(stan_model)
  }
  stanmodel <- gsub("iter", "", stan_model)  # same Stan code for iteroparity
  pars <- match.arg(pars, several.ok = TRUE)

  par_list <- list( 
    IPM_SS_np = list(
      hyper = c("alpha","alpha_W","alpha_H","delta_alpha","beta_alpha",
                "Rmax","Rmax_W","Rmax_H","delta_Rmax","beta_Rmax",
                "beta_R","rho_R","sigma_R","mu_p","sigma_p","R_p",
                "mu_SS","beta_SS","rho_SS","sigma_SS","tau"),
      states = c("R","p","s_SS","S","q","p_HOS")
    ),
    
    IPM_SS_pp = list(
      hyper = c("mu_alpha","mu_alpha_W","mu_alpha_H","delta_mu_alpha","beta_alpha","sigma_alpha",
                "mu_Rmax","mu_Rmax_W","mu_Rmax_H","delta_mu_Rmax","beta_Rmax","sigma_Rmax",
                "rho_alphaRmax","R_alphaRmax","beta_R","sigma_year_R","rho_R","sigma_R",
                "mu_p","sigma_pop_p","R_pop_p","sigma_p","R_p",
                "mu_SS","beta_SS","rho_SS","sigma_year_SS","sigma_SS","tau"),
      group = c("alpha","alpha_W","alpha_H","delta_alpha","Rmax","Rmax_W","Rmax_H","delta_Rmax",
                "eta_year_R","mu_pop_alr_p","eta_year_SS"),
      states = c("R","p","s_SS","S","p_HOS","q")
    ),
    
    IPM_SMS_np = list(
      hyper = c("alpha","alpha_W","alpha_H","delta_alpha","beta_alpha",
                "Mmax","Mmax_W","Mmax_H","delta_Mmax","beta_Mmax",
                "beta_M","rho_M","sigma_M","mu_MS","beta_MS","rho_MS","sigma_MS",
                "mu_p","sigma_p","R_p","tau_M","tau_S"),
      states = c("M","s_MS","p","S","q","p_HOS")
    ),
    
    IPM_SMS_pp = list(
      hyper = c("mu_alpha","mu_alpha_W","mu_alpha_H","delta_mu_alpha","beta_alpha","sigma_alpha",
                "mu_Mmax","mu_Mmax_W","mu_Mmax_H","delta_mu_Mmax","beta_Mmax","sigma_Mmax",
                "rho_alphaMmax","R_alphaMmax","beta_M","rho_M","sigma_year_M","sigma_M","tau_M",
                "mu_MS","beta_MS","rho_MS","sigma_year_MS","sigma_MS",
                "mu_p","sigma_pop_p","R_pop_p","sigma_p","R_p","tau_S"),
      group = c("alpha","alpha_W","alpha_H","delta_alpha","Mmax","Mmax_W","Mmax_H","delta_Mmax",
                "eta_year_M","eta_year_MS","mu_pop_alr_p"),
      states = c("M","s_MS","p","S","q","p_HOS")
    ),
    
    IPM_SMaS_np = list(
      hyper = c("alpha","alpha_W","alpha_H","delta_alpha","beta_alpha",
                "Mmax","Mmax_W","Mmax_H","delta_Mmax","beta_Mmax",
                "beta_M","rho_M","sigma_M","mu_p_M","sigma_p_M","R_p_M","tau_M",
                "mu_MS","beta_MS","rho_MS","sigma_MS","R_MS",
                "mu_p_MS","sigma_p_MS","R_p_MS","tau_S"),
      states = c("p_M","M","q_M","s_MS","p_MS","S","q_MS","q_GR","p_HOS")
    ),
    
    # IPM_ICchinook_pp = list(
    #   hyper = c(if("mu_alpha" %in% RRS) c("mu_alpha_W","mu_alpha_H","delta_mu_alpha") else "mu_alpha",
    #             "beta_alpha","sigma_alpha",
    #             switch(SR_fun, exp = NULL,
    #                    c(if("mu_Mmax" %in% RRS) c("mu_Mmax_W","mu_Mmax_H","delta_mu_Mmax") else "mu_Mmax",
    #                      "beta_Mmax","sigma_Mmax")),
    #             ifelse(identical(RRS, "none"), "rho_alphaMmax", "R_alphaMmax"),
    #             "beta_M","rho_M","sigma_M","mu_D","beta_D","rho_D","sigma_D",
    #             "mu_SAR","beta_SAR","rho_SAR","sigma_SAR","mu_U","beta_U","rho_U","sigma_U",
    #             "mu_p","sigma_pop_p","R_pop_p","sigma_p","R_p","tau_S"),
    #   group = c("alpha","alpha_W","alpha_H","delta_alpha","Mmax","Mmax_W","Mmax_H","delta_Mmax",
    #             "s_D","SAR","s_U","mu_pop_alr_p"),
    #   states = c("M","p","p_HOS","S","q")
    # ),
    
    IPM_LCRchum_pp = list(
      hyper = c("mu_E","sigma_E","delta_NG",
                "mu_psi","mu_psi_W","mu_psi_H","delta_mu_psi","beta_psi","sigma_psi",
                "mu_Mmax","mu_Mmax_W","mu_Mmax_H","delta_mu_Mmax","beta_Mmax","sigma_Mmax",
                "beta_M","rho_M","sigma_year_M","sigma_M","mu_tau_M","sigma_tau_M",
                "mu_MS","beta_MS","rho_MS","sigma_year_MS","sigma_MS",
                "mu_p","sigma_pop_p","R_pop_p","sigma_p","R_p",
                "mu_F","sigma_pop_F","sigma_F","p_D","mu_tau_S","sigma_tau_S"),
      group = c("psi","psi_W","psi_H","delta_psi","Mmax","Mmax_W","Mmax_H","delta_Mmax",
                "eta_year_M","eta_year_MS","mu_pop_alr_p"),
      states = c("M","tau_M","s_MS","p","p_F","S","tau_S","q","q_F","q_O","p_HOS")
    ),
    
    RR_SS_np = list(
      hyper = c("alpha","Rmax","rho_R","sigma_R","R_hat"),
      ppd = c("S_sim","R_sim")
    ),
    
    RR_SS_pp = list(
      hyper = c("mu_alpha","sigma_alpha","mu_Rmax","sigma_Rmax","rho_alphaRmax",
                "rho_R","sigma_year_R","sigma_R"),
      group = c("alpha","Rmax","eta_year_R"),
      ppd = c("R_hat","S_sim","R_sim")
    )
  )
  
  pars_out <- if(identical(pars, "all")) {
    unlist(par_list[[stanmodel]])
  } else {
    unlist(par_list[[stanmodel]][pars])
  }
  
  # drop unused pars based on RRS
  RRS_opts <- unique(gsub('(.+)_[WH]$', '\\1', grep('_[WH]$', pars_out, value = TRUE)))
  not_RRS <- RRS_opts[!grepl(paste(RRS, collapse = "|"), RRS_opts)]
  RRS_hyper_opts <- intersect(RRS_opts, par_list[[stanmodel]]$hyper)
  RRS_hyper <- grep(paste(RRS, collapse = "|"), RRS_hyper_opts, value = TRUE)
  drop_pars <- c(RRS, RRS_hyper, # if base param is in RRS, drop it and keep H/W versions
                 paste0(not_RRS, "_W"), paste0(not_RRS, "_H"), paste0("delta_", not_RRS),
                 paste0(ifelse(identical(RRS, "none"), "R_", "rho_"), 
                        paste(gsub("mu_", "", RRS_hyper_opts), collapse = "")))

  # drop maximum recruitment for density-independent S-R function
  if(SR_fun %in% c("exp","DI")) 
    drop_pars <- c(drop_pars, grep("max", pars_out, value = TRUE))
  
  # drop unused iteroparity pars
  if(!grepl("iter", stan_model))
    drop_pars <- c(drop_pars, grep("_SS", pars_out, value = TRUE))
  
  # drop unused betas based on par_models 
  betas <- grep("beta_", pars_out, value = TRUE)
  modeled_pars <- sapply(par_models, function(f) all.vars(f)[1])
  modeled_pars <- gsub("s_", "", modeled_pars)  # hack to deal with s_SS, s_MS
  drop_pars <- c(drop_pars, setdiff(betas, paste0("beta_", modeled_pars)))
  
  pars_out <- setdiff(pars_out, drop_pars)
  return(pars_out)
}
