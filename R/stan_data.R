#' Prepare input data for IPMs or run-reconstruction models
#' 
#' This function is mostly used internally, but may occasionally be useful for diagnosing
#' problems (e.g., checking the numeric coding of populations and years or the 
#' replacement values for `NA`s) or for plotting.
#' 
#' @param stan_model Character string specifying the **salmonIPM** model to be fit.
#' @inheritParams salmonIPM
#'
#' @return A named list that is passed to [rstan::sampling()] as the `data`
#'   argument used when fitting **salmonIPM** models.
#'
#' @examples
#' # Simulate data for a multi-population spawner-to-spawner model
#' set.seed(1234)
#' N_pop <- 10
#' N_year <- 20
#' N <- N_pop*N_year
#' 
#' pars <- list(mu_alpha = 2, sigma_alpha = 0.5, mu_Rmax = 5, sigma_Rmax = 0.5, 
#'              rho_alphaRmax = 0.3, rho_R = 0.7, sigma_year_R = 0.5, sigma_R = 0.3, 
#'              tau = 0.5, mu_p = c(0.05, 0.55, 0.4), sigma_pop_p = c(0.1, 0.2), 
#'              R_pop_p = diag(2), sigma_p = c(0.5, 0.5), R_p = diag(2), S_init_K = 0.7)
#' 
#' fd <- data.frame(pop = rep(1:N_pop, each = N_year), year = rep(1:N_year, N_pop),
#'                  A = 1, p_HOS = 0, F_rate = rbeta(N,7,3), B_rate = 0,
#'                  n_age_obs = 50, n_HW_obs = 0)
#' 
#' sim_out <- simIPM(pars = pars, fish_data = fd, N_age = 3, max_age = 5)
#' 
#' # Prepare simulated data for Stan
#' dat <- stan_data("IPM_SS_pp", fish_data = sim_out$sim_dat)
#'
#' @importFrom tidyr pivot_longer starts_with
#' @export

stan_data <- function(stan_model = c("IPM_SS_np","IPM_SSiter_np","IPM_SS_pp","IPM_SSiter_pp",
                                     "IPM_SMS_np","IPM_SMS_pp","IPM_SMaS_np",
                                     "IPM_LCRchum_pp","RR_SS_np","RR_SS_pp"), 
                      SR_fun = "BH", RRS = "none", ages = NULL, 
                      par_models = NULL, center = TRUE, scale = TRUE, 
                      priors = NULL, fish_data, age_F = NULL, age_B = NULL, 
                      age_S_obs = NULL, age_S_eff = NULL, conditionGRonMS = FALSE, 
                      fecundity_data = NULL)
{
  stan_model <- match.arg(stan_model)
  
  # Helper function to simplify capturing regex from string
  cregex <- function(string, regex) {
    rgxpr <- regexpr(regex, string)
    return(regmatches(string, rgxpr))
  }
  
  # Basic objects and data dimensions
  # General error-checking
  fish_data <- as.data.frame(fish_data)
  for(i in c("pop","year","A","fit_p_HOS","B_take_obs"))
    if(!is.null(fish_data[[i]]) & any(is.na(fish_data[[i]])))
      stop("Missing values not allowed in fish_data$", i)
  for(i in names(fish_data)) assign(i, fish_data[[i]])
  N <- nrow(fish_data)
  pop <- as.numeric(factor(pop))
  year <- as.numeric(factor(year))
  if(any(unlist(tapply(year, pop, diff)) != 1))
    stop("Non-consecutive years not allowed in fish_data")
  pool_pops <- switch(strsplit(stan_model, "_")[[1]][3], np = FALSE, pp = TRUE)
  life_cycle <- strsplit(stan_model, "_")[[1]][2]
  iter <- grepl("iter", life_cycle)   # iteroparous life history?
  model_life_cycle <- paste(strsplit(stan_model, "_")[[1]][1], 
                            gsub("iter", "", life_cycle), # same Stan code for iteroparity
                            sep = "_")
  validate_RRS(stan_model = stan_model, SR_fun = SR_fun, RRS = RRS)
  validate_par_models(stan_model = stan_model, SR_fun = SR_fun, RRS = RRS,
                      par_models = par_models)
  X <- par_model_matrix(par_models = par_models, fish_data = fish_data, 
                        center = center, scale = scale)

  # Age of fixed subadult stages
  if(!life_cycle %in% c("SS","SSiter","SMaS") & (any(is.na(ages)) | is.null(ages)))
    stop("Multi-stage models must specify age for all fixed-age subadult stages")
  
  # Age-frequency data
  if(life_cycle == "SMaS") {
    if(!is.logical(conditionGRonMS))
      stop("conditionGRonMS must be TRUE or FALSE for model 'IPM_SMaS_np'")
    
    n_Mage_obs <- as.matrix(fish_data[, grep("n_Mage", names(fish_data))])
    if(any(is.na(n_Mage_obs))) 
      stop("NA found in compositional data (n_Mage_obs). If an age class was not seen ",
           "in the sample or no sample was taken, the observed frequency is 0.")
    N_Mage <- ncol(n_Mage_obs)
    if(all(n_Mage_obs == 0) & conditionGRonMS)
      stop("ocean age-composition (n_MSage_obs) must be included if conditionGRonMS == TRUE")
    max_Mage <- max(as.numeric(cregex(colnames(n_Mage_obs), "\\d+")))
    
    n_MSage_obs <- as.matrix(fish_data[, grep("n_MSage", names(fish_data))])
    if(any(is.na(n_MSage_obs))) 
      stop("NA found in compositional data (n_MSage_obs). If an age class was not seen ",
           "in the sample or no sample was taken, the observed frequency is 0.")
    N_MSage <- ncol(n_MSage_obs)
    max_MSage <- max(as.numeric(cregex(colnames(n_MSage_obs), "\\d+")))
    max_age <- max_Mage + max_MSage
    
    n_GRage_obs <- as.matrix(fish_data[, grep("n_GRage", names(fish_data))])
    if(any(is.na(n_GRage_obs)))
      stop("NA found in compositional data (n_GRage_obs). If an age class was not seen ",
           "in the sample or no sample was taken, the observed frequency is 0.")
    N_GRage <- sum(grepl("n_GRage", names(fish_data)))
    if(N_GRage != N_Mage*N_MSage)
      stop("ncol(n_GRage_obs) does not equal ncol(n_Mage_obs) * ncol(n_MSage_obs).
           Make sure Gilbert-Rich age frequencies include all combinations of smolt
           and (implicitly) ocean age. See ?salmonIPM for details.")
  } else {
    n_age_obs <- as.matrix(fish_data[, grep("n_age", names(fish_data))])
    if(any(is.na(n_age_obs)))
      stop("NA found in compositional data (n_age_obs). If an age class was not seen ", 
           "in the sample or no sample was taken, the observed frequency is 0.")
    if(!life_cycle == "SSiter" & any(colSums(n_age_obs) == 0))
      stop("n_age_obs contains at least one age class that was never observed")
    N_age <- ncol(n_age_obs) / ifelse(life_cycle == "SSiter", 2, 1)
    max_age <- max(as.numeric(cregex(colnames(n_age_obs), "\\d+"))) - 
      ifelse(life_cycle == "SSiter", 1, 0)
  }

  # Origin-frequency, sex-frequency, and fecundity data
  # Conditional broodstock transfer probabilities
  if(life_cycle == "LCRchum") {
    if(any(is.na(n_M_obs + n_F_obs))) 
      stop("NA found in compositional data (n_[M|F]_obs). If a sex was not seen ",
           "in the sample or no sample was taken, the observed frequency is 0.")
    if(is.null(fecundity_data))
      stop("Fecundity data must be provided for Lower Columbia chum IPM") 
    if(any(is.na(fecundity_data)))
      stop("Missing observations are not allowed in fecundity data")
    
    which_H_pop <- grep("Hatchery", levels(factor(fish_data$pop)))
    which_D_pop <- which(!1:max(pop) %in% which_H_pop)
    use_M_obs <- !is.na(M_obs) & !(M_obs == 0 & pop %in% which_H_pop)
    use_S_obs <- !is.na(S_obs) & !(S_obs == 0 & pop %in% which_H_pop)

    n_O_obs <- as.matrix(fish_data[, grep("n_O", names(fish_data))])
    if(any(is.na(n_O_obs))) 
      stop("NA found in compositional data (n_O_obs). If an origin was not seen ",
           "in the sample or no sample was taken, the observed frequency is 0.")
    which_O_pop <- sapply(strsplit(colnames(n_O_obs[,-1]), "_"), 
                          function(x) as.numeric(substring(x[2], 2)))
    which_U_pop <- which(!1:max(pop) %in% which_O_pop)
    # identify cells of n_O_obs that are structural zeros because
    # (1) a given known-origin pop did not exist in any parental brood year
    # (2) a known-origin pop did not receive any broodstock in a given year
    N_O_pop <- length(which_O_pop)
    use_n_O_obs <- matrix(NA, N, 1 + N_O_pop)
    use_n_O_obs[,1] <- 
      ifelse(pop %in% which_O_pop, !is.na(S_obs[i]) & S_obs[i] > 0, TRUE)
    for(i in 1:N) {
      for(j in 1:N_O_pop) {
        use_n_O_obs[i,1+j] <- 
          any(pop == which_O_pop[j] & year == year[i] - (max_age - N_age + 1) & 
                !is.na(S_obs) & S_obs > 0) &
          ifelse(pop[i] %in% which_O_pop, !is.na(S_obs[i]) & S_obs[i] > 0, TRUE)
      }
    }
    
    n_B_obs <- fish_data[, grep("n_B", names(fish_data))]
    if(any(is.na(n_B_obs))) 
      stop("NA found in compositional data (n_B_obs). If an origin was not seen ",
           "in the broodstock or no broodstock was taken, the observed frequency is 0.")
    p_B_obs <- cbind(location = factor(fish_data$pop), year = year, 
                     sweep(n_B_obs, 1, pmax(B_take_obs, 1), "/"))
    p_B_obs <- pivot_longer(p_B_obs, cols = starts_with("n_B"), values_to = "p",
                            names_to = "disposition", names_pattern = "n_B(.*)_obs", 
                            names_transform = as.integer)
    p_B_obs$disposition <- factor(levels(p_B_obs$location)[p_B_obs$disposition],
                                  levels = levels(p_B_obs$location))
    P_B <- tapply(p_B_obs$p, list(p_B_obs$year, p_B_obs$disposition, p_B_obs$location),
                  FUN = identity, default = 0)
  } else {
    fit_p_HOS <- as.logical(fit_p_HOS)
    if(any(is.na(n_W_obs + n_H_obs)))
      stop("NA found in compositional data (n_[W|H]_obs). If a rearing type was not seen ",
           "in the sample or no sample was taken, the observed frequency is 0.")
    if(any(!fit_p_HOS & n_H_obs > 0))
      warning("Hatchery spawners present in rows ", 
              paste(which(!fit_p_HOS & n_H_obs > 0), collapse = ", "),
              " but fit_p_HOS == FALSE. See ?salmonIPM to learn why this is problematic.")
    if(any(fit_p_HOS & n_H_obs + n_W_obs == 0))
      warning("fit_p_HOS == TRUE in rows ", which(fit_p_HOS & n_H_obs + n_W_obs == 0),
              " but origin-frequency samples are missing. ",
              "See ?salmonIPM to learn why this is problematic.")
    n_W_obs <- n_W_obs[fit_p_HOS]
    n_H_obs <- n_H_obs[fit_p_HOS]
  }
  
  # F_rate error checking
  F_rate_NA <- tapply(F_rate, pop, function(x) any(is.na(x[-c(1:max_age)])))
  if(any(F_rate_NA))
    stop("Missing values not allowed in fish_data$F_rate except in years 1:max_age")
  
  # Fishery and broodstock age-selectivity
  if(is.null(age_F)) 
    age_F <- rep(1, switch(life_cycle, SSiter = N_age + 1, SMaS = N_MSage, N_age))
  age_F <- as.numeric(age_F)
  
  if(is.null(age_B)) 
    age_B <- rep(1, switch(life_cycle, SSiter = N_age + 1, SMaS = N_MSage, N_age))
  age_B <- as.numeric(age_B)
  
  # Partially observed and effective age indicators
  if(model_life_cycle == "IPM_SS") {
    if(is.null(age_S_obs)) 
      age_S_obs <- rep(1, switch(life_cycle, SSiter = N_age + 1, N_age))
    age_S_obs <- as.numeric(age_S_obs)
    
    if(is.null(age_S_eff)) 
      age_S_eff <- rep(1, switch(life_cycle, SSiter = N_age + 1, N_age))
    age_S_eff <- as.numeric(age_S_eff)
  }
  
  # Priors
  if(is.null(priors)) {
    priors <- list()
  } else {   # formulas must be 2-sided
    stopifnot(all(sapply(priors, function(f) attr(terms(f), "response"))))
  }
  pars <- unlist(sapply(priors, function(f) as.character(f[[2]])))
  pdfs <- lapply(priors, function(f) f[[3]])
  names(pdfs) <- pars
  stopifnot(all(sapply(pdfs, is.call)))
  
  if(stan_model %in% c("IPM_SS_np","IPM_SSiter_np","IPM_SS_pp","IPM_SSiter_pp")) {
    # autoscale Rmax to data
    log_RA <- log((S_obs + B_take_obs)/((1 - F_rate)*A)) 
    log_RA <- log_RA[is.finite(log_RA)]
    prior_Rmax_mean <- quantile(log_RA, 0.8, na.rm = TRUE, names = FALSE)
    prior_Rmax_sd <- 2*sd(log_RA, na.rm = TRUE)
    
    if(grepl("_np", stan_model)) {
      pars <- match.arg(pars, c("alpha","alpha_W","alpha_H","Rmax","Rmax_W","Rmax_H","mu_p","mu_SS","tau"), several.ok = TRUE)
      prior_alpha <- stan_prior(pdfs$alpha, lognormal(2,2))
      prior_alpha_W <- stan_prior(pdfs$alpha_W, lognormal(2,2)) 
      prior_alpha_H <- stan_prior(pdfs$alpha_H, lognormal(2,2))
      prior_Rmax <- stan_prior(pdfs$Rmax, lognormal(prior_Rmax_mean, prior_Rmax_sd))
      prior_Rmax_W <- stan_prior(pdfs$Rmax_W, lognormal(prior_Rmax_mean, prior_Rmax_sd))
      prior_Rmax_H <- stan_prior(pdfs$Rmax_H, lognormal(prior_Rmax_mean, prior_Rmax_sd))
      prior_tau <- stan_prior(pdfs$tau, gnormal(1, 0.85, 30)) # squash tau < 0.1 to avoid divergences
    }
    
    if(grepl("_pp", stan_model)) {
      pars <- match.arg(pars, c("mu_alpha","mu_alpha_W","mu_alpha_H",
                                "mu_Rmax","mu_Rmax_W","mu_Rmax_H","mu_p","mu_SS","tau"), 
                        several.ok = TRUE)
      prior_mu_alpha <- stan_prior(pdfs$mu_alpha, normal(2,5))
      prior_mu_alpha_W <- stan_prior(pdfs$mu_alpha_W, normal(2,5))
      prior_mu_alpha_H <- stan_prior(pdfs$mu_alpha_H, normal(2,5))
      prior_mu_Rmax <- stan_prior(pdfs$mu_Rmax, normal(prior_Rmax_mean, prior_Rmax_sd))
      prior_mu_Rmax_W <- stan_prior(pdfs$mu_Rmax_W, normal(prior_Rmax_mean, prior_Rmax_sd))
      prior_mu_Rmax_H <- stan_prior(pdfs$mu_Rmax_H, normal(prior_Rmax_mean, prior_Rmax_sd))
      prior_tau <- stan_prior(pdfs$tau, gnormal(0,1,2)) # default normal
    }
    
    prior_mu_p <- stan_prior(pdfs$mu_p, dirichlet(rep(1, N_age)))
    prior_mu_SS <- stan_prior(pdfs$mu_SS, beta(1,1))
  }
  
  if(stan_model %in% c("IPM_SMS_np","IPM_SMS_pp","IPM_SMaS_np","IPM_ICchinook_pp")) {
    if(stan_model == "IPM_ICchinook_pp") {  
      # autoscale Mmax to data with ballpark SAR
      log_MA <- log((S_obs + B_take_obs)/((1 - F_rate)*0.01*A))
    } else {  
      # autoscale Mmax to smolt data
      log_MA <- log(M_obs/A)
    }
    log_MA <- log_MA[is.finite(log_MA)]
    prior_Mmax_mean <- quantile(log_MA, 0.8, na.rm = TRUE, names = FALSE)
    prior_Mmax_sd <- 2*sd(log_MA, na.rm = TRUE)
    
    if(grepl("_np", stan_model)) {
      pars <- match.arg(pars, c("alpha","Mmax","mu_MS","mu_p","tau_M","tau_S"), several.ok = TRUE)
      prior_alpha <- stan_prior(pdfs$alpha, lognormal(5,5))
      prior_alpha_W <- stan_prior(pdfs$alpha_W, lognormal(5,5))
      prior_alpha_H <- stan_prior(pdfs$alpha_H, lognormal(5,5))
      prior_Mmax <- stan_prior(pdfs$Mmax, lognormal(prior_Mmax_mean, prior_Mmax_sd))
      prior_Mmax_W <- stan_prior(pdfs$Mmax_W, lognormal(prior_Mmax_mean, prior_Mmax_sd))
      prior_Mmax_H <- stan_prior(pdfs$Mmax_H, lognormal(prior_Mmax_mean, prior_Mmax_sd))
      prior_tau_M <- stan_prior(pdfs$tau_M, gnormal(1, 0.85, 30)) # squash tau < 0.1 to avoid divergences
      prior_tau_S <- stan_prior(pdfs$tau_S, gnormal(1, 0.85, 30)) # squash tau < 0.1 to avoid divergences
    }
    
    if(grepl("_pp", stan_model)) {
      pars <- match.arg(pars, c("mu_alpha","mu_Mmax","mu_MS","mu_p","tau_M","tau_S"), several.ok = TRUE)
      prior_mu_alpha <- stan_prior(pdfs$mu_alpha, normal(5,5))
      prior_mu_alpha_W <- stan_prior(pdfs$mu_alpha_W, normal(5,5))
      prior_mu_alpha_H <- stan_prior(pdfs$mu_alpha_H, normal(5,5))
      prior_mu_Mmax <- stan_prior(pdfs$mu_Mmax, normal(prior_Mmax_mean, prior_Mmax_sd))
      prior_mu_Mmax_W <- stan_prior(pdfs$mu_Mmax_W, normal(prior_Mmax_mean, prior_Mmax_sd))
      prior_mu_Mmax_H <- stan_prior(pdfs$mu_Mmax_H, normal(prior_Mmax_mean, prior_Mmax_sd))
      prior_tau_M <- stan_prior(pdfs$tau_M, gnormal(0,1,2)) # normal
      prior_tau_S <- stan_prior(pdfs$tau_S, gnormal(0,1,2)) # normal
    }
    
    if(grepl("SMS", stan_model)) {
      prior_mu_p <- stan_prior(pdfs$mu_p, dirichlet(rep(1, N_age)))
      prior_mu_MS <- stan_prior(pdfs$mu_MS, beta(1,1))
    }
  }
  
  if(stan_model == "IPM_LCRchum_pp") {
    # autoscale Mmax to smolt data (can == 0 for hatcheries)
    log_MA <- log(replace(M_obs, M_obs == 0, NA)/A) 
    log_MA <- log_MA[is.finite(log_MA)]
    prior_Mmax_mean <- quantile(log_MA, 0.8, na.rm = TRUE, names = FALSE)
    prior_Mmax_sd <- 2*sd(log_MA, na.rm = TRUE)
    
    pars <- match.arg(pars, c("mu_psi","mu_Mmax","mu_MS","mu_p"), several.ok = TRUE)
    prior_mu_psi <- stan_prior(pdfs$mu_psi, beta(1,1))
    prior_mu_psi_W <- stan_prior(pdfs$mu_psi_W, beta(1,1))
    prior_mu_psi_H <- stan_prior(pdfs$mu_psi_H, beta(1,1))
    prior_mu_Mmax <- stan_prior(pdfs$mu_Mmax, normal(prior_Mmax_mean, prior_Mmax_sd))
    prior_mu_Mmax_W <- stan_prior(pdfs$mu_Mmax_W, normal(prior_Mmax_mean, prior_Mmax_sd))
    prior_mu_Mmax_H <- stan_prior(pdfs$mu_Mmax_H, normal(prior_Mmax_mean, prior_Mmax_sd))
    prior_mu_MS <- stan_prior(pdfs$mu_MS, beta(1,1))
    prior_mu_p <- stan_prior(pdfs$mu_p, dirichlet(rep(1, N_age)))
  }
  
  # # Prior data
  # ## Deprecated: move this into fish_data columns
  # if(life_cycle == "ICchinook") {
  #   if(is.null(prior_data)) {
  #     stop("Priors for survival must be specified for model IPM_ICchinook_pp")
  #   } else {
  #     s <- prior_data$s[prior_data$s$year %in% fish_data$year,]
  #     
  #     prior_D <- na.omit(s[,c("year","mu_prior_D","sigma_prior_D")])
  #     N_prior_D <- nrow(prior_D)
  #     which_prior_D <- match(prior_D$year, sort(unique(fish_data$year)))
  #     mu_prior_D <- prior_D$mu_prior_D
  #     sigma_prior_D <- prior_D$sigma_prior_D
  #     
  #     prior_SAR <- na.omit(s[,c("year","mu_prior_SAR","sigma_prior_SAR")])
  #     N_prior_SAR <- nrow(prior_SAR)
  #     which_prior_SAR <- match(prior_SAR$year, sort(unique(fish_data$year)))
  #     mu_prior_SAR <- prior_SAR$mu_prior_SAR
  #     sigma_prior_SAR <- prior_SAR$sigma_prior_SAR
  #     
  #     prior_U <- na.omit(s[,c("year","mu_prior_U","sigma_prior_U")])
  #     N_prior_U <- nrow(prior_U)
  #     which_prior_U <- match(prior_U$year, sort(unique(fish_data$year)))
  #     mu_prior_U <- prior_U$mu_prior_U
  #     sigma_prior_U <- prior_U$sigma_prior_U
  #   }
  # }
  
  # Run reconstruction
  if(model_life_cycle == "RR_SS") {
    rr_dat <- run_recon(fish_data, age_F = age_F, age_B = age_B)
    which_fit <- which(!is.na(rr_dat$R_obs) & !is.na(rr_dat$S_obs))
    N_fit <- length(which_fit)
    R_obs <- rr_dat$R_obs
    p_pop_obs <- aggregate(rr_dat[,grep("p_age", names(rr_dat))], list(pop = rr_dat$pop), mean, na.rm = TRUE)[,-1]
  }
  
  ## Data to return
  
  out <- switch(model_life_cycle,
                IPM_SS = list(
                  # info for observed data
                  N = N,
                  pop = pop, 
                  year = year,
                  # recruitment
                  SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                  RRS = array(c(any(grepl("alpha", RRS)), any(grepl("Rmax", RRS)))),
                  A = A,
                  K_alpha = ifelse(is.null(X$alpha), 0, ncol(X$alpha)), 
                  X_alpha = if(is.null(X$alpha)) matrix(0,N,0) else X$alpha,
                  prior_alpha = if(pool_pops) NULL else prior_alpha,
                  prior_alpha_W = if(pool_pops) NULL else prior_alpha_W,
                  prior_alpha_H = if(pool_pops) NULL else prior_alpha_H,
                  prior_mu_alpha = if(pool_pops) prior_mu_alpha else NULL,
                  prior_mu_alpha_W = if(pool_pops) prior_mu_alpha_W else NULL,
                  prior_mu_alpha_H = if(pool_pops) prior_mu_alpha_H else NULL,
                  K_Rmax = ifelse(is.null(X$Rmax), 0, ncol(X$Rmax)), 
                  X_Rmax = if(is.null(X$Rmax)) matrix(0,N,0) else X$Rmax,
                  prior_Rmax = if(pool_pops) NULL else prior_Rmax,
                  prior_Rmax_W = if(pool_pops) NULL else prior_Rmax_W,
                  prior_Rmax_H = if(pool_pops) NULL else prior_Rmax_H,
                  prior_mu_Rmax = if(pool_pops) prior_mu_Rmax else NULL,
                  prior_mu_Rmax_W = if(pool_pops) prior_mu_Rmax_W else NULL,
                  prior_mu_Rmax_H = if(pool_pops) prior_mu_Rmax_H else NULL,
                  K_R = ifelse(is.null(X$R), 0, ncol(X$R)), 
                  X_R = if(is.null(X$R)) matrix(0,N,0) else X$R,
                  # kelt survival
                  iter = as.integer(iter),
                  K_SS = ifelse(is.null(X$s_SS), 0, ncol(X$s_SS)), 
                  X_SS = if(is.null(X$s_SS)) matrix(0,N,0) else X$s_SS,
                  prior_mu_SS = prior_mu_SS,
                  # spawner abundance
                  N_S_obs = sum(!is.na(S_obs)),
                  which_S_obs = as.vector(which(!is.na(S_obs))),
                  S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                  prior_tau = prior_tau,
                  # spawner or [maiden | repeat] age structure
                  N_age = N_age,
                  max_age = max_age,
                  n_age_obs = as.matrix(fish_data[,grep("n_age", names(fish_data))]),
                  age_S_obs = age_S_obs,
                  age_S_eff = age_S_eff,
                  prior_mu_p = prior_mu_p,
                  # H/W composition
                  N_H = sum(fit_p_HOS),
                  which_H = as.vector(which(fit_p_HOS)),
                  n_W_obs = n_W_obs,
                  n_H_obs = n_H_obs,
                  # fishery and hatchery removals
                  F_rate = replace(F_rate, is.na(F_rate), 0),
                  age_F = age_F,
                  N_B = sum(B_take_obs > 0),
                  which_B = as.vector(which(B_take_obs > 0)),
                  B_take_obs = B_take_obs[B_take_obs > 0],
                  age_B = age_B
                ),
                
                IPM_SMS = list(
                  # info for observed data
                  N = N,
                  pop = pop, 
                  year = year,
                  # smolt production
                  SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                  RRS = array(c(any(grepl("alpha", RRS)), any(grepl("Mmax", RRS)))),
                  A = A,
                  K_alpha = ifelse(is.null(X$alpha), 0, ncol(X$alpha)), 
                  X_alpha = if(is.null(X$alpha)) matrix(0,N,0) else X$alpha,
                  prior_alpha = if(pool_pops) NULL else prior_alpha,
                  prior_alpha_W = if(pool_pops) NULL else prior_alpha_W,
                  prior_alpha_H = if(pool_pops) NULL else prior_alpha_H,
                  prior_mu_alpha = if(pool_pops) prior_mu_alpha else NULL,
                  prior_mu_alpha_W = if(pool_pops) prior_mu_alpha_W else NULL,
                  prior_mu_alpha_H = if(pool_pops) prior_mu_alpha_H else NULL,
                  K_Mmax = ifelse(is.null(X$Mmax), 0, ncol(X$Mmax)), 
                  X_Mmax = if(is.null(X$Mmax)) matrix(0,N,0) else X$Mmax,
                  prior_Mmax = if(pool_pops) NULL else prior_Mmax,
                  prior_Mmax_W = if(pool_pops) NULL else prior_Mmax_W,
                  prior_Mmax_H = if(pool_pops) NULL else prior_Mmax_H,
                  prior_mu_Mmax = if(pool_pops) prior_mu_Mmax else NULL,
                  prior_mu_Mmax_W = if(pool_pops) prior_mu_Mmax_W else NULL,
                  prior_mu_Mmax_H = if(pool_pops) prior_mu_Mmax_H else NULL,
                  K_M = ifelse(is.null(X$M), 0, ncol(X$M)), 
                  X_M = if(is.null(X$M)) matrix(0,N,0) else X$M,
                  # smolt abundance
                  N_M_obs = sum(!is.na(M_obs)),
                  which_M_obs = as.vector(which(!is.na(M_obs))),
                  M_obs = replace(M_obs, is.na(M_obs) | M_obs==0, 1),
                  prior_tau_M = prior_tau_M,
                  smolt_age = ages$M,
                  # SAR
                  K_MS = ifelse(is.null(X$s_MS), 0, ncol(X$s_MS)), 
                  X_MS = if(is.null(X$s_MS)) matrix(0,N,0) else X$s_MS,
                  prior_mu_MS = prior_mu_MS,
                  # spawner abundance
                  N_S_obs = sum(!is.na(S_obs)),
                  which_S_obs = as.vector(which(!is.na(S_obs))),
                  S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                  prior_tau_S = prior_tau_S,
                  # spawner age structure
                  N_age = N_age, 
                  max_age = max_age,
                  n_age_obs = n_age_obs,
                  prior_mu_p = prior_mu_p,
                  # H/W composition
                  N_H = sum(fit_p_HOS),
                  which_H = as.vector(which(fit_p_HOS)),
                  n_W_obs = n_W_obs,
                  n_H_obs = n_H_obs,
                  # fishery and hatchery removals
                  F_rate = replace(F_rate, is.na(F_rate), 0),
                  age_F = age_F,
                  N_B = sum(B_take_obs > 0),
                  which_B = as.vector(which(B_take_obs > 0)),
                  B_take_obs = B_take_obs[B_take_obs > 0],
                  age_B = age_B
                ),
                
                IPM_SMaS = list(
                  # info for observed data
                  N = N,
                  pop = pop, 
                  year = year,
                  # smolt production
                  SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                  RRS = array(c(any(grepl("alpha", RRS)), any(grepl("Mmax", RRS)))),
                  A = A,
                  K_alpha = ifelse(is.null(X$alpha), 0, ncol(X$alpha)), 
                  X_alpha = if(is.null(X$alpha)) matrix(0,N,0) else X$alpha,
                  prior_alpha = prior_alpha,
                  prior_alpha_W = if(pool_pops) NULL else prior_alpha_W,
                  prior_alpha_H = if(pool_pops) NULL else prior_alpha_H,
                  K_Mmax = ifelse(is.null(X$Mmax), 0, ncol(X$Mmax)), 
                  X_Mmax = if(is.null(X$Mmax)) matrix(0,N,0) else X$Mmax,
                  prior_Mmax = prior_Mmax,
                  prior_Mmax_W = if(pool_pops) NULL else prior_Mmax_W,
                  prior_Mmax_H = if(pool_pops) NULL else prior_Mmax_H,
                  K_M = ifelse(is.null(X$M), 0, ncol(X$M)), 
                  X_M = if(is.null(X$M)) matrix(0,N,0) else X$M,
                  # smolt abundance
                  N_M_obs = sum(!is.na(M_obs)),
                  which_M_obs = as.vector(which(!is.na(M_obs))),
                  M_obs = replace(M_obs, is.na(M_obs) | M_obs==0, 1),
                  prior_tau_M = prior_tau_M,
                  # smolt age structure
                  N_Mage = N_Mage,
                  max_Mage = max_Mage,
                  n_Mage_obs = n_Mage_obs,
                  # SAR
                  K_MS = ifelse(is.null(X$s_MS), 0, ncol(X$s_MS)), 
                  X_MS = if(is.null(X$s_MS)) matrix(0,N,0) else X$s_MS,
                  # spawner abundance
                  N_S_obs = sum(!is.na(S_obs)),
                  which_S_obs = as.vector(which(!is.na(S_obs))),
                  S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                  prior_tau_S = prior_tau_S,
                  # spawner ocean and Gilbert-Rich age structure
                  N_MSage =  N_MSage, 
                  max_MSage = max_MSage,
                  n_MSage_obs = n_MSage_obs,
                  n_GRage_obs = n_GRage_obs,
                  conditionGRonMS = as.numeric(conditionGRonMS),
                  # H/W composition
                  N_H = sum(fit_p_HOS),
                  which_H = as.vector(which(fit_p_HOS)),
                  n_W_obs = n_W_obs,
                  n_H_obs = n_H_obs,
                  # fishery and hatchery removals
                  F_rate = replace(F_rate, is.na(F_rate), 0),
                  MSage_F = age_F,
                  N_B = sum(B_take_obs > 0),
                  which_B = as.vector(which(B_take_obs > 0)),
                  B_take_obs = B_take_obs[B_take_obs > 0],
                  MSage_B = age_B
                ),
                
                IPM_LCRchum = list(
                  # info for observed data
                  N = N,
                  pop = pop, 
                  year = year,
                  # egg deposition
                  SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                  RRS = array(c(any(grepl("psi", RRS)), any(grepl("Mmax", RRS)))),
                  A = A,
                  N_E = nrow(fecundity_data),
                  age_E = fecundity_data$age_E,
                  E_obs = fecundity_data$E_obs,
                  # spawner-smolt productivity
                  smolt_age = ages$M,
                  K_psi = ifelse(is.null(X$psi), 0, ncol(X$psi)), 
                  X_psi = if(is.null(X$psi)) matrix(0,N,0) else X$psi,
                  prior_mu_psi = prior_mu_psi,
                  prior_mu_psi_W = prior_mu_psi_W,
                  prior_mu_psi_H = prior_mu_psi_H,
                  K_Mmax = ifelse(is.null(X$Mmax), 0, ncol(X$Mmax)), 
                  X_Mmax = if(is.null(X$Mmax)) matrix(0,N,0) else X$Mmax,
                  prior_mu_Mmax = prior_mu_Mmax,
                  prior_mu_Mmax_W = prior_mu_Mmax_W,
                  prior_mu_Mmax_H = prior_mu_Mmax_H,
                  K_M = ifelse(is.null(X$M), 0, ncol(X$M)), 
                  X_M = if(is.null(X$M)) matrix(0,N,0) else X$M,
                  # smolt abundance and observation error
                  N_M_obs = sum(use_M_obs),
                  which_M_obs = as.vector(which(use_M_obs)),
                  M_obs = replace(M_obs, is.na(M_obs) | M_obs==0, 1),
                  N_tau_M_obs = sum(!is.na(tau_M_obs)),
                  which_tau_M_obs = as.vector(which(!is.na(tau_M_obs))),
                  N_tau_M_fit = sum(!is.na(tau_M_obs) & !pop %in% which_H_pop),
                  which_tau_M_fit = as.vector(which(!is.na(tau_M_obs) & !pop %in% which_H_pop)),
                  tau_M_obs = replace(tau_M_obs, is.na(tau_M_obs), 0),
                  N_upstream = sum(!is.na(downstream_trap)),
                  which_upstream = as.vector(which(!is.na(downstream_trap))),
                  downstream_trap = na.omit(downstream_trap),
                  # SAR
                  K_MS = ifelse(is.null(X$s_MS), 0, ncol(X$s_MS)), 
                  X_MS = if(is.null(X$s_MS)) matrix(0,N,0) else X$s_MS,
                  prior_mu_MS = prior_mu_MS,
                  # spawner abundance and observation error
                  N_S_obs = sum(use_S_obs),
                  which_S_obs = as.vector(which(use_S_obs)),
                  S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                  N_tau_S_obs = sum(!is.na(tau_S_obs)),
                  which_tau_S_obs = as.vector(which(!is.na(tau_S_obs))),
                  N_tau_S_fit = sum(!is.na(tau_S_obs) & !pop %in% which_O_pop),
                  which_tau_S_fit = as.vector(which(!is.na(tau_S_obs) & !pop %in% which_O_pop)),
                  tau_S_obs = replace(tau_S_obs, is.na(tau_S_obs), 0),
                  # spawner age structure and sex ratio
                  N_age = N_age,
                  max_age = max_age,
                  n_age_obs = n_age_obs,
                  prior_mu_p = prior_mu_p,
                  n_MF_obs = as.vector(n_M_obs + n_F_obs),
                  n_F_obs = as.vector(n_F_obs),
                  p_G_obs = as.vector(p_G_obs),
                  # origin composition
                  N_O_pop = N_O_pop,
                  which_O_pop = which_O_pop,
                  N_U_pop = length(which_U_pop),
                  which_U_pop = which_U_pop,
                  N_n_O_obs = sum(use_n_O_obs),
                  which_n_O_obs = which(as.vector(use_n_O_obs)),
                  n_O_obs = n_O_obs,
                  # fishery and hatchery removals and translocations
                  F_rate = replace(F_rate, is.na(F_rate), 0),
                  age_F = age_F,
                  N_H_pop = length(which_H_pop),
                  which_H_pop = which_H_pop,
                  N_D_pop = length(which_D_pop),
                  which_D_pop = which_D_pop,
                  N_B = sum(B_take_obs > 0),
                  which_B = as.vector(which(B_take_obs > 0)),
                  B_take_obs = B_take_obs[B_take_obs > 0],
                  age_B = age_B,
                  P_B = P_B
                ),
                
                # IPM_ICchinook = list(
                #   # info for observed data
                #   N = N,
                #   pop = pop, 
                #   year = year,
                #   # info for forward simulations
                #   N_fwd = N_fwd,
                #   pop_fwd = as.vector(fish_data_fwd$pop),
                #   year_fwd = as.vector(fish_data_fwd$year),
                #   A_fwd = as.vector(fish_data_fwd$A),
                #   B_rate_fwd = as.vector(fish_data_fwd$B_rate),
                #   F_rate_fwd = as.vector(fish_data_fwd$F_rate),
                #   p_HOS_fwd = as.vector(fish_data_fwd$p_HOS),
                #   # recruitment
                #   SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                #   smolt_age = ages$M,
                #   A = A,
                #   K_alpha = ifelse(is.null(X$alpha), 0, ncol(X$alpha)), 
                #   X_alpha = if(is.null(X$alpha)) matrix(0,N,0) else X$alpha,
                #   prior_mu_alpha = prior_mu_alpha,
                #   K_Mmax = ifelse(is.null(X$Mmax), 0, ncol(X$Mmax)), 
                #   X_Mmax = if(is.null(X$Mmax)) matrix(0,N,0) else X$Mmax,
                #   prior_mu_Mmax = prior_mu_Mmax,
                #   K_M = ifelse(is.null(X$M), 0, ncol(X$M)), 
                #   X_M = if(is.null(X$M)) matrix(0,N,0) else X$M,
                #   # downstream, SAR, upstream survival
                #   # K_D = ncol(env_data$s_D),
                #   # X_D = as.matrix(env_data$s_D),
                #   N_prior_D = N_prior_D,
                #   which_prior_D = as.vector(which_prior_D),
                #   mu_prior_D = as.vector(mu_prior_D),
                #   sigma_prior_D = as.vector(sigma_prior_D),
                #   # K_SAR = ncol(env_data$s_SAR),
                #   # X_SAR = as.matrix(env_data$s_SAR),
                #   N_prior_SAR = N_prior_SAR,
                #   which_prior_SAR = as.vector(which_prior_SAR),
                #   mu_prior_SAR = as.vector(mu_prior_SAR),
                #   sigma_prior_SAR = as.vector(sigma_prior_SAR),
                #   # K_U = ncol(env_data$s_U),
                #   # X_U = as.matrix(env_data$s_U),
                #   N_prior_U = N_prior_U,
                #   which_prior_U = as.vector(which_prior_U),
                #   mu_prior_U = as.vector(mu_prior_U),
                #   sigma_prior_U = as.vector(sigma_prior_U),
                #   # spawner abundance
                #   N_S_obs = sum(!is.na(S_obs)),
                #   which_S_obs = as.vector(which(!is.na(S_obs))),
                #   S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                #   prior_tau_S = prior_tau_S,
                #   # spawner age structure
                #   N_age = N_age, 
                #   max_age = max_age,
                #   n_age_obs = n_age_obs,
                #   prior_mu_p = prior_mu_p,
                #   # H/W composition
                #   N_H = sum(fit_p_HOS),
                #   which_H = as.vector(which(fit_p_HOS)),
                #   n_W_obs = n_W_obs,
                #   n_H_obs = n_H_obs,
                #   # fishery and hatchery removals
                #   F_rate = replace(F_rate, is.na(F_rate), 0),
                #   age_F = age_F,
                #   N_B = sum(B_take_obs > 0),
                #   which_B = as.vector(which(B_take_obs > 0)),
                #   B_take_obs = B_take_obs[B_take_obs > 0],
                #   age_B = age_B
                # ),
                
                RR_SS = list(
                  SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                  N = length(S_obs),
                  pop = pop, 
                  year = year,
                  N_fit = N_fit,
                  which_fit = as.vector(which_fit),
                  S_obs = replace(S_obs, S_obs == 0 | is.na(S_obs), 1),
                  R_obs = replace(R_obs, R_obs == 0 | is.na(R_obs), 1),
                  A = A,
                  S_NA = as.vector(as.integer(is.na(S_obs))),
                  R_NA = as.vector(as.integer(is.na(R_obs))),
                  N_age = N_age,
                  max_age = max(as.numeric(substring(grep("p_age", names(rr_dat), value = TRUE), 6, 6))),
                  p_pop_obs = as.matrix(p_pop_obs)
                )
  )  # end switch()
  
  return(out)
}
