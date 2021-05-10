#' Generate initial values for fitting IPMs or run-reconstruction spawner-recruit models.
#'
#' @inheritParams salmonIPM
#' @param data Named list of input data for fitting either an IPM or
#'   run-reconstruction spawner-recruit model in [Stan](http://mc-stan.org), 
#'   as returned by [stan_data()].
#' @param chains A positive integer specifying the number of Markov chains.
#'
#' @importFrom stats aggregate na.omit
#'
#' @return A list with initial starting values for the parameters and states in
#'   the model.
#'
#' @export
stan_init <- function(stan_model, data, chains = 1) 
{
  if(!stan_model %in% c("IPM_SS_np","IPM_SS_pp","IPM_SMS_np","IPM_SMS_pp",
                        "IPM_SMaS_np","IPM_LCRchum_pp","IPM_ICchinook_pp",
                        "RR_SS_np","RR_SS_pp"))
    stop("Stan model ", stan_model, " does not exist")
  
  model <- strsplit(stan_model, "_")[[1]][1]
  for(i in names(data)) assign(i, data[[i]])
  
  if(model == "IPM") {
    S_obs_noNA <- S_obs
    S_obs[-which_S_obs] <- NA
    p_HOS_obs <- pmin(pmax(n_H_obs/(n_H_obs + n_W_obs), 0.1), 0.9)
    p_HOS_obs[n_H_obs + n_W_obs == 0] <- 0.5
    p_HOS_all <- rep(0,N)
    p_HOS_all[which_H] <- p_HOS_obs
    min_age <- max_age - N_age + 1
    adult_ages <- min_age:max_age
    q_obs <- sweep(n_age_obs, 1, rowSums(n_age_obs), "/")
    q_obs_NA <- apply(is.na(q_obs), 1, any)
    q_obs[q_obs_NA,] <- rep(colMeans(na.omit(q_obs)), each = sum(q_obs_NA))
    R_a <- matrix(NA, N, N_age)
    S_W_obs <- S_obs*(1 - p_HOS_all)
    B_rate_all <- rep(0,N)
    B_rate <- pmin(pmax(B_take_obs/(S_W_obs[which_B]*(1 - q_obs[which_B,1]) + B_take_obs), 0.01), 0.99)
    B_rate[is.na(B_rate)] <- 0.1
    B_rate_all[which_B] <- B_rate
    
    # Maybe figure out a way to do this with a call to run_recon?
    for(i in 1:N)
      for(j in 1:N_age)
      {
        if(year[i] + adult_ages[j] <= max(year[pop==pop[i]]))
        {
          b <- ifelse(j==1, 0, B_rate_all[i+adult_ages[j]])
          f <- ifelse(j==1, 0, F_rate[i + adult_ages[j]])
          R_a[i,j] <- S_W_obs[i + adult_ages[j]]*q_obs[i + adult_ages[j],j]/((1 - b)*(1 - f))
        }
      }
    
    R_a <- pmax(R_a, min(1, R_a[R_a > 0], na.rm = T))
    R <- rowSums(R_a)
    R[is.na(R)] <- max(R, na.rm = T)
    p <- sweep(R_a, 1, R, "/")
    p_NA <- apply(is.na(p), 1, any)
    p[p_NA, ] <- rep(colMeans(na.omit(p)), each = sum(p_NA))
    alr_p <- sweep(log(p[, 1:(N_age-1), drop = FALSE]), 1, log(p[,N_age]), "-")
    zeta_p <- apply(alr_p, 2, scale)
    mu_p <- aggregate(p, list(pop), mean)
    zeta_pop_p <- aggregate(alr_p, list(pop), mean)[,-1, drop = FALSE]
    zeta_pop_p <- apply(zeta_pop_p, 2, scale)
  }
  
  if(stan_model %in% c("IPM_SMS_np","IPM_SMS_pp","IPM_LCRchum_pp")) {
    if(N_M_obs < N)
      M_obs[-which_M_obs] <- median(M_obs[which_M_obs])
    s_MS <- pmin(S_obs_noNA/M_obs, 0.9)
  }
  
  if(stan_model == "IPM_SMaS_np") {
    # This is a bit of a hack to avoid tedious age-structured run reconstruction
    N_GRage <- N_Mage*N_MSage
    max_age <- max_Mage + max_MSage
    q_M_obs <- sweep(n_Mage_obs, 1, rowSums(n_Mage_obs), "/")
    s_MS <- mean(pmin(S_obs/M_obs, 0.99), na.rm = TRUE)
    q_MS_obs <- sweep(n_MSage_obs, 1, rowSums(n_MSage_obs), "/")
    q_MS_obs_NA <- apply(is.na(q_MS_obs), 1, any)
    q_MS_obs[q_MS_obs_NA,] <- rep(colMeans(na.omit(q_MS_obs)), each = sum(q_MS_obs_NA))
    q_MS_pop <- as.matrix(aggregate(q_MS_obs, list(pop), mean))[,-1,drop=FALSE]
    mu_q_MS <- array(0, c(N_pop,N_Mage,N_MSage))
    for(i in 1:N_pop)
      for(j in 1:N_Mage)
        mu_q_MS[i,j,] <- as.vector(q_MS_pop[i,])
    q_GR_obs <- sweep(n_GRage_obs, 1, rowSums(n_GRage_obs), "/")
    p_HOS_obs <- pmin(pmax(n_H_obs/(n_H_obs + n_W_obs), 0.01), 0.99)
    p_HOS_obs[n_H_obs + n_W_obs == 0] <- 0.5
    p_HOS_all <- rep(0,N)
    p_HOS_all[which_H] <- p_HOS_obs
    S_W_obs <- S_obs*(1 - p_HOS_all)
    B_rate <- pmin(pmax(B_take_obs/(S_W_obs[which_B] + B_take_obs), 0.01), 0.99)
    B_rate[is.na(B_rate)] <- 0.1
  }
  
  if(stan_model == "IPM_LCRchum_pp") {
    E <- S_obs_noNA*0.5*mean(fecundity_data$E_obs)
    s_EM <- pmin(M_obs/E, 0.9)
  }
  
  out <- lapply(1:chains, function(i) {
    switch(stan_model,
           IPM_SS_np = list(
             # recruitment
             alpha = array(exp(runif(N_pop, 1, 3)), dim = N_pop),
             beta_alpha = matrix(rnorm(K_alpha*N_pop, 0, 1), N_pop, K_alpha),
             Rmax = array(rlnorm(N_pop, log(tapply(R/A, pop, quantile, 0.9)), 0.5), dim = N_pop),
             beta_Rmax = matrix(rnorm(K_Rmax*N_pop, 0, 1), N_pop, K_Rmax),
             beta_R = matrix(rnorm(K_R*N_pop, 0, 1), N_pop, K_R),
             rho_R = array(runif(N_pop, 0.1, 0.7), dim = N_pop),
             sigma_R = array(runif(N_pop, 0.05, 2), dim = N_pop), 
             zeta_R = as.vector(scale(log(R)))*0.1,
             # spawner age structure
             mu_p = mu_p,
             sigma_p = matrix(runif(N_pop*(N_age-1), 0.5, 1), N_pop, N_age-1),
             zeta_p = zeta_p,
             # H/W composition, removals
             p_HOS = p_HOS_obs,
             B_rate = B_rate,
             # initial spawners, observation error
             S_init = rep(median(S_obs_noNA), max_age*N_pop),
             q_init = matrix(colMeans(q_obs), max_age*N_pop, N_age, byrow = T),
             tau = array(runif(N_pop, 0.5, 1), dim = N_pop)
           ),
           
           IPM_SS_pp = list(
             # recruitment
             mu_alpha = runif(1, 1, 3),
             beta_alpha = array(rnorm(K_alpha, 0, 1), dim = K_alpha),
             sigma_alpha = runif(1, 0.1, 0.5),
             zeta_alpha = as.vector(runif(N_pop, -1, 1)),
             mu_Rmax = rnorm(1, log(quantile(R/A,0.9)), 0.5),
             beta_Rmax = array(rnorm(K_Rmax, 0, 1), dim = K_Rmax),
             sigma_Rmax = runif(1, 0.1, 0.5),
             zeta_Rmax = as.vector(runif(N_pop,-1,1)),
             rho_alphaRmax = runif(1, -0.5, 0.5),
             beta_R = array(rnorm(K_R, 0, 1), dim = K_R),
             rho_R = runif(1, 0.1, 0.7),
             sigma_year_R = runif(1, 0.1, 0.5),
             zeta_year_R = as.vector(rnorm(max(year, year_fwd), 0, 0.1)),
             sigma_R = runif(1, 0.5, 1),
             zeta_R = as.vector(scale(log(R)))*0.1,
             # spawner age structure
             mu_p = colMeans(p),
             sigma_pop_p = array(runif(N_age - 1, 0.5, 1), dim = N_age - 1),
             zeta_pop_p = zeta_pop_p,
             sigma_p = array(runif(N_age-1, 0.5, 1), dim = N_age - 1),
             zeta_p = zeta_p,
             # H/W composition, removals
             p_HOS = p_HOS_obs,
             B_rate = B_rate,
             # initial spawners, observation error
             S_init = rep(median(S_obs_noNA), max_age*N_pop),
             q_init = matrix(colMeans(q_obs), max_age*N_pop, N_age, byrow = T),
             tau = runif(1, 0.5, 1)
           ),
           
           IPM_SMS_np = list(
             # smolt recruitment
             alpha = array(exp(runif(N_pop,1,3)), dim = N_pop),
             beta_alpha = matrix(rnorm(K_alpha*N_pop, 0, 1), N_pop, K_alpha),
             Mmax = array(rlnorm(N_pop, log(tapply(R/A, pop, quantile, 0.9)), 0.5), dim = N_pop),
             beta_Mmax = matrix(rnorm(K_Mmax*N_pop, 0, 1), N_pop, K_Mmax),
             beta_M = matrix(rnorm(K_M*N_pop,0,1), N_pop, K_M),
             rho_M = array(runif(N_pop, 0.1, 0.7), dim = N_pop),
             sigma_M = array(runif(N_pop, 0.05, 2), dim = N_pop), 
             zeta_M = as.vector(scale(log(M_obs)))*0.1,
             # SAR
             mu_MS = array(plogis(rnorm(N_pop, mean(qlogis(s_MS)), 0.5)), dim = N_pop),
             beta_MS = matrix(rnorm(K_MS*N_pop,0,1), N_pop, K_MS),
             rho_MS = array(runif(N_pop, 0.1, 0.7), dim = N_pop),
             sigma_MS = array(runif(N_pop, 0.05, 2), dim = N_pop), 
             zeta_MS = as.vector(scale(qlogis(s_MS))),
             # spawner age structure
             mu_p = mu_p,
             sigma_p = matrix(runif(N_pop*(N_age-1),0.5,1), N_pop, N_age-1),
             zeta_p = zeta_p,
             # H/W composition, removals
             p_HOS = p_HOS_obs,
             B_rate = B_rate,
             # initial states, observation error
             M_init = array(rep(median(M_obs), smolt_age*N_pop), dim = smolt_age*N_pop),
             S_init = array(rep(median(S_obs_noNA), (max_age - smolt_age)*N_pop), dim = (max_age - smolt_age)*N_pop),
             q_init = matrix(colMeans(q_obs), (max_age - smolt_age)*N_pop, N_age, byrow = T),
             tau_M = array(runif(N_pop, 0.5, 1), dim = N_pop),
             tau_S = array(runif(N_pop, 0.5, 1), dim = N_pop)
           ),
           
           IPM_SMS_pp = list(
             # smolt recruitment
             mu_alpha = runif(1, 1, 3),
             beta_alpha = array(rnorm(K_alpha, 0, 1), dim = K_alpha),
             sigma_alpha = runif(1, 0.1, 0.5),
             zeta_alpha = as.vector(rnorm(N_pop, 0, 1)),
             mu_Mmax = rnorm(1, log(quantile(R/A,0.9)), 0.5),
             beta_Mmax = array(rnorm(K_Mmax, 0, 1), dim = K_Mmax),
             sigma_Mmax = runif(1, 0.1, 0.5),
             zeta_Mmax = as.vector(rnorm(N_pop, 0, 1)),
             rho_alphaMmax = runif(1, -0.5, 0.5),
             beta_M = array(rnorm(K_M, 0, 1), dim = K_M),
             rho_M = runif(1, 0.1, 0.7),
             sigma_year__M = runif(1, 0.1, 0.5),
             zeta_year__M = as.vector(rnorm(max(year), 0, 0.1)),
             sigma_M = runif(1, 0.5, 1),
             zeta_M = as.vector(scale(log(M_obs)))*0.1,
             # SAR
             mu_MS = plogis(rnorm(1, mean(qlogis(s_MS)), 0.5)),
             beta_MS = array(rnorm(K_MS,0,1), dim = K_MS),
             rho_MS = runif(1, 0.1, 0.7),
             sigma_year__MS = runif(1, 0.05, 2), 
             sigma_MS = runif(1, 0.5, 1),
             zeta_MS = as.vector(scale(qlogis(s_MS))),
             # spawner age structure
             mu_p = colMeans(p),
             sigma_pop_p = array(runif(N_age - 1, 0.5, 1), dim = N_age - 1),
             zeta_pop_p = zeta_pop_p,
             sigma_p = array(runif(N_age-1, 0.5, 1), dim = N_age - 1),
             zeta_p = zeta_p,
             # H/W composition, removals
             p_HOS = p_HOS_obs,
             B_rate = B_rate,
             # initial states, observation error
             M_init = array(rep(median(M_obs), smolt_age*N_pop), dim = smolt_age*N_pop),
             S_init = array(rep(median(S_obs_noNA), (max_age - smolt_age)*N_pop), dim = (max_age - smolt_age)*N_pop),
             q_init = matrix(colMeans(q_obs), (max_age - smolt_age)*N_pop, N_age, byrow = T),
             tau_M = runif(1, 0.5, 1),
             tau_S = runif(1, 0.5, 1)
           ),
           
           IPM_LCRchum_pp = list(
             # egg deposition
             mu_E = rlnorm(N_age, tapply(log(E_obs), age_E, mean), 0.5),
             sigma_E = rlnorm(N_age, log(tapply(E_obs, age_E, sd)), 0.5), 
             delta_NG = runif(1, 0.7, 1),
             # egg-smolt survival
             mu_psi = plogis(rnorm(1, mean(qlogis(s_EM)), 0.3)),
             beta_psi = array(rnorm(K_psi, 0, 0.3), dim = K_psi),
             sigma_psi = runif(1, 0.1, 0.5),
             zeta_psi = rnorm(N_pop, 0, 1),
             mu_Mmax = rnorm(1, mean(log(S_obs[which_S_obs])), 3),
             beta_Mmax = array(rnorm(K_Mmax, 0, 1), dim = K_Mmax),
             sigma_Mmax = runif(1, 0.5, 2),
             zeta_Mmax = rnorm(N_pop, 0, 1),
             beta_M = array(rnorm(K_M, 0, 1), dim = K_M),
             rho_M = runif(1, 0.1, 0.7),
             sigma_year_M = runif(1, 0.1, 0.5),
             zeta_year_M = rnorm(max(year), 0, 0.1),
             sigma_M = runif(1, 0.1, 0.5),
             zeta_M = as.vector(scale(log(M_obs)))*0.1,
             # SAR
             mu_MS = plogis(rnorm(1, mean(qlogis(s_MS)), 0.5)),
             beta_MS = array(rnorm(K_MS,0,1), dim = K_MS),
             rho_MS = runif(1, 0.1, 0.7),
             sigma_year_MS = runif(1, 0.05, 2), 
             zeta_year_MS = as.vector(tapply(scale(qlogis(s_MS)), year, mean)),
             sigma_MS = runif(1, 0.5, 1),
             zeta_MS = as.vector(scale(qlogis(s_MS))),
             # spawner age structure and sex ratio
             mu_p = colMeans(p),
             sigma_pop_p = runif(N_age - 1, 0.5, 1),
             zeta_pop_p = zeta_pop_p,
             sigma_p = runif(N_age-1, 0.5, 1),
             zeta_p = zeta_p,
             mu_F = runif(1,0.4,0.6),
             sigma_pop_F = runif(1, 0.5, 0.3),
             zeta_pop_F = rnorm(N_pop, 0, 0.3),
             sigma_F = runif(1, 0.5, 0.3),
             zeta_F = rnorm(N, 0, 0.3),
             # H/W composition, removals
             p_HOS = p_HOS_obs,
             B_rate = B_rate,
             # initial states, observation error
             M_init = rep(median(M_obs), smolt_age*N_pop),
             S_init = rep(median(S_obs_noNA), (max_age - smolt_age)*N_pop),
             q_init = matrix(colMeans(q_obs), (max_age - smolt_age)*N_pop, N_age, byrow = T),
             mu_tau_M = runif(1, 0, 0.5),
             sigma_tau_M = runif(1, 0, 0.5),
             mu_tau_S = runif(1, 0, 0.5),
             sigma_tau_S = runif(1, 0, 0.5)
           ),
           
           IPM_ICchinook_pp = list(
             # smolt recruitment
             mu_alpha = runif(1, 1, 3),
             beta_alpha = array(rnorm(K_alpha, 0, 1), dim = K_alpha),
             sigma_alpha = runif(1, 0.1, 0.5),
             zeta_alpha = runif(N_pop, -1, 1),
             mu_Mmax = rnorm(1, log(quantile(R/A,0.9)), 0.5),
             beta_Mmax = array(rnorm(K_Mmax, 0, 1), dim = K_Mmax),
             sigma_Mmax = runif(1, 0.1, 0.5),
             zeta_Mmax = runif(N_pop, -1, 1),
             rho_alphaMmax = runif(1, -0.5, 0.5),
             beta_M = array(rnorm(K_M, 0, 1), dim = K_M),
             rho_M = runif(1, 0.1, 0.7),
             sigma_M = runif(1, 0.05, 2), 
             zeta_M = as.vector(scale(log(R)))*0.01,
             M_init = rep(median(S_obs_noNA)*100, smolt_age*N_pop),
             # downstream, SAR, upstream survival
             mu_D = qlogis(0.8),
             beta_D = array(rnorm(K_D, 0, 1), dim = K_D),
             rho_D = runif(1, 0.1, 0.7),
             sigma_D = runif(1, 0.05, 2),
             zeta_D = rnorm(max(year,year_fwd), 0, 0.1),
             mu_SAR = qlogis(0.01),
             beta_SAR = array(rnorm(K_SAR, 0, 1), dim = K_SAR),
             rho_SAR = runif(1, 0.1, 0.7),
             sigma_SAR = runif(1, 0.05, 2),
             zeta_SAR = rnorm(max(year,year_fwd), 0, 0.1),
             mu_U = qlogis(0.8),
             beta_U = array(rnorm(K_U, 0, 1), dim = K_U),
             rho_U = runif(1, 0.1, 0.7),
             sigma_U = runif(1, 0.05, 2),
             zeta_U = rnorm(max(year,year_fwd), 0, 0.1),
             # spawner age structure
             mu_p = colMeans(p),
             sigma_pop_p = runif(N_age - 1, 0.5, 1),
             zeta_pop_p = zeta_pop_p,
             sigma_p = runif(N_age-1, 0.5, 1),
             zeta_p = zeta_p,
             # H/W composition, removals
             p_HOS = p_HOS_obs,
             B_rate = B_rate,
             # initial spawners, observation error
             S_init = rep(median(S_obs_noNA), (max_age - smolt_age)*N_pop),
             q_init = matrix(colMeans(q_obs), (max_age - smolt_age)*N_pop, N_age, byrow = TRUE),
             tau_S = runif(1, 0.5, 1)
           ),
           
           IPM_SMaS_np = list(
             # smolt recruitment
             alpha = array(rlnorm(N_pop, max(log(M_obs/S_obs), na.rm = TRUE), 1), dim = N_pop),
             beta_alpha = matrix(rnorm(K_alpha*N_pop, 0, 1), N_pop, K_alpha),
             Mmax = array(rlnorm(N_pop, log(tapply(M_obs/A, pop, quantile, 0.9, na.rm = TRUE)), 0.5), dim = N_pop),
             beta_Mmax = matrix(rnorm(K_Mmax*N_pop, 0, 1), N_pop, K_Mmax),
             beta_M = matrix(rnorm(K_M*N_pop,0,1), N_pop, K_M),
             rho_M = array(runif(N_pop, 0.1, 0.7), dim = N_pop),
             sigma_M = array(runif(N_pop, 0.05, 2), dim = N_pop), 
             zeta_M = rnorm(N,0,0.1), 
             # smolt age structure
             mu_p_M = aggregate(q_M_obs, list(pop), mean, na.rm = TRUE),
             sigma_p_M = matrix(runif(N_pop*(N_Mage - 1), 0.05, 2), N_pop, N_Mage - 1),
             zeta_p_M = matrix(rnorm(N*(N_Mage - 1), 0, 0.1), N, N_Mage - 1),
             # SAR
             mu_MS = matrix(plogis(rnorm(N_pop*N_Mage, qlogis(s_MS), 0.5)), N_pop, N_Mage),
             beta_MS = matrix(rnorm(K_MS*N_pop,0,1), N_pop, K_MS),
             rho_MS = matrix(runif(N_pop, 0.1, 0.7), N_pop, N_Mage),
             sigma_MS = matrix(runif(N_pop, 0.05, 2), N_pop, N_Mage), 
             zeta_MS = matrix(rnorm(N*N_Mage, 0, 0.1), N, N_Mage),
             # ocean age structure
             mu_p_MS = mu_q_MS,
             sigma_p_MS = array(runif(N_pop*N_Mage*(N_MSage - 1), 0.05, 2), 
                                c(N_pop, N_Mage, N_MSage - 1)), 
             zeta_p_MS = matrix(rnorm(N*N_Mage*(N_MSage - 1), 0, 0.5), N, N_Mage*(N_MSage - 1)),
             # H/W composition, removals
             p_HOS = p_HOS_obs,
             B_rate = B_rate,
             # initial states, observation error
             M_init = array(rep(median(M_obs), max_Mage*N_pop), dim = max_Mage*N_pop),
             q_M_init = matrix(colMeans(q_M_obs, na.rm = TRUE), max_Mage*N_pop, N_Mage, byrow = T),
             S_init = array(rep(median(S_obs, na.rm = TRUE), N_pop*max_MSage), dim = N_pop*max_MSage),
             q_GR_init = matrix(colMeans(q_GR_obs, na.rm = TRUE), max_MSage*N_pop, N_GRage, byrow = T),
             tau_S = array(runif(N_pop, 0.01, 0.05), dim = N_pop),
             tau_M = array(runif(N_pop, 0.01, 0.05), dim = N_pop)
           ),
           
           RR_SS_pp = list(
             # This is currently not based on the input data
             mu_alpha = runif(1, 3, 6), 
             sigma_alpha = runif(1, 0.1, 0.5),
             zeta_alpha = array(runif(N_pop, -1, 1), dim = N_pop), 
             mu_Rmax = rnorm(1, log(quantile(S/A, 0.9, na.rm = T)), 0.5),
             sigma_Rmax = runif(1, 0.1, 0.5),
             zeta_Rmax = array(runif(N_pop, -1, 1), dim = N_pop), 
             rho_alphaRmax = runif(1, -0.5, 0.5),
             rho_R = runif(1, 0.1, 0.7),
             sigma_year_R = runif(1, 0.1, 0.5), 
             sigma_R = runif(1, 0.1, 2)
           ),
           
           RR_SS_np = list(
             # This is currently not based on the input data
             alpha = array(exp(runif(N_pop, 1, 3)), dim = N_pop),
             Rmax = array(exp(runif(N_pop, -1, 0)), dim = N_pop),
             rho_R = array(runif(N_pop, 0.1, 0.7), dim = N_pop),
             sigma_R = array(runif(N_pop, 0.5, 1), dim = N_pop)
           )
    )  # end switch()
  })  # end lappply()
  
  return(out)
}
