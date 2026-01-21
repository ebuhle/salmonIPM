functions {
  #include /include/SR.stan
  #include /include/gnormal_lpdf.stan
  #include /include/mat_lmult.stan
  #include /include/row_sums.stan
}

data {
  // info for observed data
  int<lower=1> N;                      // total number of cases in all pops and years
  array[N] int<lower=1,upper=N> pop;   // population index
  array[N] int<lower=1,upper=N> year;  // calendar year index
  // smolt production
  int<lower=1> SR_fun;                 // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  array[2] int<lower=0,upper=1> RRS;   // fit W vs. H {alpha, Mmax} (1) or not (0)?
  int<lower=1> smolt_age;              // smolt age
  vector[N] A;                         // habitat area associated with each spawner abundance obs
  int<lower=0> K_alpha;                // number of intrinsic productivity covariates
  matrix[N,K_alpha] X_alpha;           // intrinsic productivity covariates
  array[2] real prior_mu_alpha;        // prior mean, sd for hyper-mean log intrinsic productivity
  array[2] real prior_mu_alpha_W;      // prior mean, sd for hyper-mean log W intrinsic productivity
  array[2] real prior_mu_alpha_H;      // prior mean, sd for hyper-mean log H intrinsic productivity
  int<lower=0> K_Mmax;                 // number of maximum smolt recruitment covariates
  matrix[N,K_Mmax] X_Mmax;             // maximum smolt recruitment covariates
  array[2] real prior_mu_Mmax;         // prior mean, sd for hyper-mean log maximum smolt recruitment
  array[2] real prior_mu_Mmax_W;       // prior mean, sd for hyper-mean log W maximum smolt recruitment
  array[2] real prior_mu_Mmax_H;       // prior mean, sd for hyper-mean log H maximum smolt recruitment
  int<lower=0> K_M;                    // number of smolt recruitment covariates
  array[N] row_vector[K_M] X_M;        // smolt recruitment covariates
  // smolt abundance
  int<lower=1,upper=N> N_M_obs;        // number of cases with non-missing smolt abundance obs 
  array[N_M_obs] int<lower=1,upper=N> which_M_obs; // cases with non-missing smolt abundance obs
  vector<lower=0>[N] M_obs;            // observed annual smolt abundance (not density)
  array[3] real prior_tau_M;           // prior mean, scale, shape for smolt observation error SD
  // SAR (sMolt-Spawner survival)
  int<lower=0> K_MS;                   // number of SAR covariates
  matrix[N,K_MS] X_MS;                 // SAR covariates
  array[2] real prior_mu_MS;           // prior a, b for mean SAR
  // spawner abundance
  int<lower=1,upper=N> N_S_obs;        // number of cases with non-missing spawner abundance obs 
  array[N_S_obs] int<lower=1,upper=N> which_S_obs; // cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;            // observed annual total spawner abundance (not density)
  array[3] real prior_tau_S;           // prior mean, scale, shape for spawner observation error SD
  // spawner age structure
  int<lower=2> N_age;                  // number of adult age classes
  int<lower=2> max_age;                // maximum adult age
  matrix<lower=0>[N,N_age] n_age_obs;  // observed wild spawner age frequencies (all zero row = NA)  
  vector<lower=0>[N_age] prior_mu_p;   // prior concentration for mean age distribution
  // H/W composition
  int<lower=0,upper=N> N_H;            // number of years with p_HOS > 0
  array[N_H] int<lower=1,upper=N> which_H; // years with p_HOS > 0
  array[N_H] int<lower=0> n_W_obs;     // count of wild spawners in samples
  array[N_H] int<lower=0> n_H_obs;     // count of hatchery spawners in samples
  // fishery and hatchery removals
  vector[N] F_rate;                     // fishing mortality rate of wild adults
  vector<lower=0,upper=1>[N_age] age_F; // is age a (non)selected (0/1) by fishery?
  int<lower=0,upper=N> N_B;             // number of years with B_take > 0
  array[N_B] int<lower=1,upper=N> which_B; // years with B_take > 0
  vector[N_B] B_take_obs;               // observed broodstock take of wild adults
  vector<lower=0,upper=1>[N_age] age_B; // is age a (non)selected (0/1) in broodstock?
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years
  array[N] int<lower=1> pop_year;          // index of years within each pop, starting at 1
  array[N_age] int<lower=0> ocean_ages;    // ocean ages
  int<lower=1> max_ocean_age = max_age - smolt_age; // maximum ocean age
  int<lower=1> min_ocean_age = max_ocean_age - N_age + 1; // minimum ocean age
  array[N_H] int<lower=0> n_HW_obs;        // total sample sizes for H/W frequencies
  real mu_M_init = mean(log(M_obs[which_M_obs])); // prior log-mean of smolt abundance in years 1:smolt_age
  real sigma_M_init = sd(log(M_obs[which_M_obs])); // prior log-SD of smolt abundance in years 1:smolt_age
  vector[max_ocean_age*N_pop] mu_S_init;   // prior mean of total spawner abundance in years 1:max_ocean_age
  real sigma_S_init = sd(log(S_obs[which_S_obs])); // prior log-SD of spawner abundance in years 1:max_ocean_age
  matrix[N_age,max_ocean_age*N_pop] mu_q_init; // prior counts of wild spawner age distns in years 1:max_ocean_age
  int<lower=0,upper=1> any_RRS = max(RRS); // does either S-R parameter differ by W vs. H?
  int<lower=2,upper=4> N_SR = 2 + sum(RRS); // number of S-R parameter hyper-means
  row_vector[3] prior_mu_alpha_mean;       // prior means for hyper-mean log intrinsic productivity (all/W/H)
  row_vector[3] prior_mu_alpha_sd;         // prior SDs for hyper-mean log intrinsic productivity (all/W/H)
  row_vector[3] prior_mu_Mmax_mean;        // prior means for hyper-mean log maximum recruitment (all/W/H)
  row_vector[3] prior_mu_Mmax_sd;          // prior SDs for hyper-mean log maximum recruitment (all/W/H)
  
  // Indices and array-based calculations
  for(a in 1:N_age) ocean_ages[a] = min_ocean_age - 1 + a;
  for(i in 1:N_H) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];

  // Years within each pop starting at 1 (determines when to use initial values)
  pop_year[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year[i] = 1;
    else
      pop_year[i] = pop_year[i-1] + 1;
  }
  
  // Priors on S_init, q_init accounting for amalgamation of orphan age classes
  for(i in 1:max_ocean_age)
  {
    int N_orphan_age = N_age - max(i - min_ocean_age, 0); // number of orphan age classes
    int N_amalg_age = N_age - N_orphan_age + 1;           // number of amalgamated age classes
    
    for(j in 1:N_pop)
    {
      int ii = (j - 1)*max_ocean_age + i; // index into S_init, q_init

      // S_init prior mean scales observed log-mean by fraction of orphan age classes
      mu_S_init[ii] = mean(log(S_obs[which_S_obs])) + log(N_orphan_age) - log(N_age);
      
      // prior on q_init implies q_orphan ~ Dir(1)
      mu_q_init[,ii] = append_row(rep_vector(1.0/N_amalg_age, N_amalg_age), 
                                  rep_vector(1, N_orphan_age - 1));
    }
  }
  
  // Extract scalar prior parameters into vectors for efficient sampling statements
  prior_mu_alpha_mean = [prior_mu_alpha[1], prior_mu_alpha_W[1], prior_mu_alpha_H[1]];
  prior_mu_alpha_sd = [prior_mu_alpha[2], prior_mu_alpha_W[2], prior_mu_alpha_H[2]];
  prior_mu_Mmax_mean = [prior_mu_Mmax[1], prior_mu_Mmax_W[1], prior_mu_Mmax_H[1]];
  prior_mu_Mmax_sd = [prior_mu_Mmax[2], prior_mu_Mmax_W[2], prior_mu_Mmax_H[2]];
}

parameters {
  // smolt recruitment
  real mu_alpha;                         // hyper-mean log intrinsic spawner-smolt productivity
  real mu_alpha_W;                       // hyper-mean log W spawner-smolt intrinsic productivity
  real mu_alpha_H;                       // hyper-mean log H spawner-smolt intrinsic productivity
  vector[K_alpha] beta_alpha;            // regression coefs for log alpha
  real<lower=0> sigma_alpha;             // hyper-SD log intrinsic spawner-smolt productivity
  vector[!RRS[1]*N_pop] zeta_alpha;      // log spawner-smolt intrinsic productivity (Z-scores)
  vector[RRS[1]*N_pop] zeta_alpha_W;     // log W spawner-smolt intrinsic productivity (Z-scores)
  vector[RRS[1]*N_pop] zeta_alpha_H;     // log H spawner-smolt intrinsic productivity (Z-scores)
  real mu_Mmax;                          // hyper-mean log maximum smolt recruitment
  real mu_Mmax_W;                        // hyper-mean log W maximum smolt recruitment
  real mu_Mmax_H;                        // hyper-mean log H maximum smolt recruitment
  vector[K_Mmax] beta_Mmax;              // regression coefs for log Mmax
  real<lower=0> sigma_Mmax;              // hyper-SD log maximum smolt recruitment
  vector[!RRS[2]*N_pop] zeta_Mmax;       // log maximum smolt recruitment (Z-scores)
  vector[RRS[2]*N_pop] zeta_Mmax_W;      // log W maximum smolt recruitment (Z-scores)
  vector[RRS[2]*N_pop] zeta_Mmax_H;      // log H maximum smolt recruitment (Z-scores)
  cholesky_factor_corr[N_SR] L_alphaMmax; // Cholesky factor of correlation among log alpha, Mmax (all/W/H)
  vector[K_M] beta_M;                    // regression coefs for smolt recruitment
  real<lower=-1,upper=1> rho_M;          // AR(1) coef for log smolt productivity anomalies
  real<lower=0> sigma_year_M;            // process error SD of log smolt productivity anomalies
  vector[N_year] zeta_year_M;            // log smolt productivity anomalies (Z-scores)
  real<lower=0> sigma_M;                 // unique smolt recruitment process error SD
  vector[N] zeta_M;                      // unique smolt recruitment process errors (Z-scores)
  // SAR
  real<lower=0,upper=1> mu_MS;           // mean SAR
  vector[K_MS] beta_MS;                  // regression coefs for logit SAR anomalies
  real<lower=-1,upper=1> rho_MS;         // AR(1) coef for logit SAR anomalies
  real<lower=0> sigma_year_MS;           // process error SD of logit SAR anomalies
  vector[N_year] zeta_year_MS;           // logit SAR anomalies (Z-scores)
  real<lower=0> sigma_MS;                // unique logit SAR process error SD
  vector[N] zeta_MS;                     // unique logit SAR process errors (Z-scores)
  // spawner age structure
  simplex[N_age] mu_p;                   // among-pop mean of age distributions
  vector<lower=0>[N_age-1] sigma_pop_p;  // among-pop SD of mean log-ratio age distributions
  cholesky_factor_corr[N_age-1] L_pop_p; // Cholesky factor of among-pop correlation matrix of mean log-ratio age distns
  matrix[N_pop,N_age-1] zeta_pop_p;      // population mean log-ratio age distributions (Z-scores)
  vector<lower=0>[N_age-1] sigma_p;      // SD of log-ratio cohort age distributions
  cholesky_factor_corr[N_age-1] L_p;     // Cholesky factor of correlation matrix of cohort log-ratio age distributions
  matrix[N,N_age-1] zeta_p;              // log-ratio cohort age distributions (Z-scores)
  // H/W composition, removals
  vector<lower=0,upper=1>[N_H] p_HOS;     // true p_HOS in years which_H
  vector<lower=0,upper=1>[N_B] B_rate;    // true broodstock take rate when B_take > 0
  // initial states, observation error
  vector<lower=0>[smolt_age*N_pop] M_init; // true smolt abundance in years 1:smolt_age
  vector<lower=0>[max_ocean_age*N_pop] S_init; // true total spawner abundance in years 1:max_ocean_age
  array[max_ocean_age*N_pop] simplex[N_age] q_init; // true wild spawner age distributions in years 1:max_ocean_age
  real<lower=0> tau_M;                   // smolt observation error SDs
  real<lower=0> tau_S;                   // spawner observation error SDs
}

transformed parameters {
  // smolt recruitment
  vector<lower=0>[!RRS[1]*N_pop] alpha;  // intrinsic spawner-smolt productivity
  vector<lower=0>[RRS[1]*N_pop] alpha_W; // W intrinsic spawner-smolt productivity
  vector<lower=0>[RRS[1]*N_pop] alpha_H; // H intrinsic spawner-smolt productivity
  vector[N] alpha_Xbeta;                 // intrinsic productivity including covariate effects
  vector[N] alpha_W_Xbeta;               // W intrinsic productivity including covariate effects
  vector[N] alpha_H_Xbeta;               // H intrinsic productivity including covariate effects
  vector<lower=0>[!RRS[2]*N_pop] Mmax;   // maximum smolt recruitment 
  vector<lower=0>[RRS[2]*N_pop] Mmax_W;  // W maximum smolt recruitment 
  vector<lower=0>[RRS[2]*N_pop] Mmax_H;  // H maximum smolt recruitment 
  vector[N] Mmax_Xbeta;                  // maximum recruitment including covariate effects
  vector[N] Mmax_W_Xbeta;                // W maximum recruitment including covariate effects
  vector[N] Mmax_H_Xbeta;                // H maximum recruitment including covariate effects
  vector[N_year] eta_year_M;             // log brood year spawner-smolt productivity anomalies
  vector<lower=0>[N] M_hat;              // expected smolt abundance (not density) by brood year
  vector<lower=0>[N] M0;                 // true smolt abundance (not density) by brood year
  vector<lower=0>[N] M;                  // true smolt abundance (not density) by outmigration year
  // SAR
  vector[N_year] eta_year_MS;            // logit SAR anomalies by outmigration year
  vector<lower=0,upper=1>[N] s_MS;       // true SAR by outmigration year
  // H/W spawner abundance, removals
  vector[N] p_HOS_all;                   // true p_HOS in all years (can == 0)
  vector<lower=0>[N] S_W;                // true total wild spawner abundance
  vector[N] S_H;                         // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                  // true total spawner abundance
  vector<lower=0,upper=1>[N] B_rate_all; // true broodstock take rate in all years
  // spawner age structure
  row_vector[N_age-1] mu_alr_p;          // mean of log-ratio cohort age distributions
  matrix[N_pop,N_age-1] mu_pop_alr_p;    // population mean log-ratio age distributions
  matrix<lower=0,upper=1>[N,N_age] p;    // true adult age distributions by outmigration year
  matrix<lower=0,upper=1>[N,N_age] q;    // true spawner age distributions
  
  // Multivariate Matt trick for log([alpha, alpha_W, alpha_H, Mmax, Mmax_W, Mmax_H])
  {
    matrix[N_pop,N_SR] zeta_alphaMmax; // random effects (z-scored)
    matrix[N_pop,N_SR] eta_alphaMmax;  // random effects
    vector[N_SR] sigma_alphaMmax;      // SD vector
    vector[N] Xbeta_alpha = exp(mat_lmult(X_alpha, beta_alpha)); // covariate effects on alpha 
    vector[N] Xbeta_Mmax = exp(mat_lmult(X_Mmax, beta_Mmax));    // covariate effects on Mmax
    
    sigma_alphaMmax = append_row(rep_vector(sigma_alpha, 1 + RRS[1]), 
                                 rep_vector(sigma_Mmax, 1 + RRS[2]));
    zeta_alphaMmax = append_col(RRS[1] ? append_col(zeta_alpha_W, zeta_alpha_H) : to_matrix(zeta_alpha),
                                RRS[2] ? append_col(zeta_Mmax_W, zeta_Mmax_H) : to_matrix(zeta_Mmax));            
    eta_alphaMmax = diag_pre_multiply(sigma_alphaMmax, L_alphaMmax * zeta_alphaMmax')';
    
    if(RRS[1]) {
      alpha_W = exp(mu_alpha_W + eta_alphaMmax[,1]);
      alpha_W_Xbeta = alpha_W[pop] .* Xbeta_alpha;
      alpha_H = exp(mu_alpha_H + eta_alphaMmax[,2]);
      alpha_H_Xbeta = alpha_H[pop] .* Xbeta_alpha;
    } else {
      alpha = exp(mu_alpha + eta_alphaMmax[,1]);
      alpha_Xbeta = alpha[pop] .* Xbeta_alpha;
    }
    
    if(RRS[2]) {
      Mmax_W = exp(mu_Mmax_W + eta_alphaMmax[,RRS[1]+2]);
      Mmax_W_Xbeta = Mmax_W[pop] .* Xbeta_Mmax;
      Mmax_H = exp(mu_Mmax_H + eta_alphaMmax[,RRS[1]+3]);
      Mmax_H_Xbeta = Mmax_H[pop] .* Xbeta_Mmax;
    } else {
      Mmax = exp(mu_Mmax + eta_alphaMmax[,RRS[1]+2]);
      Mmax_Xbeta = Mmax[pop] .* Xbeta_Mmax;
    }
  }
  
  // AR(1) model for spawner-smolt productivity and SAR anomalies
  eta_year_M[1] = zeta_year_M[1] * sigma_year_M / sqrt(1 - rho_M^2);     // initial anomaly
  eta_year_MS[1] = zeta_year_MS[1] * sigma_year_MS / sqrt(1 - rho_MS^2); // initial anomaly
  for(t in 2:N_year)
  {
    eta_year_M[t] = rho_M * eta_year_M[t-1] + zeta_year_M[t] * sigma_year_M;
    eta_year_MS[t] = rho_MS * eta_year_MS[t-1] + zeta_year_MS[t] * sigma_year_MS;
  }
  // constrain "fitted" log or logit anomalies to sum to 0
  eta_year_M = eta_year_M - mean(eta_year_M);
  eta_year_MS = eta_year_MS - mean(eta_year_MS);
  // annual population-specific SAR
  s_MS = inv_logit(logit(mu_MS) + mat_lmult(X_MS, beta_MS) + eta_year_MS[year] + zeta_MS*sigma_MS);
  
  // Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;
  
  // Multivariate Matt trick for age vectors (pop-specific mean and within-pop, time-varying)
  mu_alr_p = to_row_vector(log(head(mu_p, N_age-1)) - log(mu_p[N_age]));
  mu_pop_alr_p = rep_matrix(mu_alr_p,N_pop) + diag_pre_multiply(sigma_pop_p, L_pop_p * zeta_pop_p')';
  // Inverse log-ratio (softmax) transform of cohort age distn
  {
    matrix[N,N_age-1] alr_p = mu_pop_alr_p[pop,] + diag_pre_multiply(sigma_p, L_p * zeta_p')';
    matrix[N,N_age] exp_alr_p = append_col(exp(alr_p), rep_vector(1,N));
    p = diag_pre_multiply(inv(row_sums(exp_alr_p)), exp_alr_p);
  }

  // Calculate true total wild and hatchery spawners, spawner age distribution, and smolts,
  // and predict smolt recruitment from brood year i
  for(i in 1:N)
  {
    int ii;                  // index into S_init and q_init
    // number of orphan age classes <lower=0,upper=N_age>
    int N_orphan_age = max(N_age - to_int(fdim(pop_year[i], min_ocean_age)), N_age); 
    vector[N_orphan_age] q_orphan; // orphan age distribution (amalgamated simplex)
    row_vector[N_age] S_W_a; // true wild spawners by age

    // Smolt recruitment
    if(pop_year[i] <= smolt_age) // use initial values
      M[i] = M_init[(pop[i]-1)*smolt_age + pop_year[i]];  
    else  // smolts from appropriate brood year
      M[i] = M0[i-smolt_age];  
    
    // Spawners and age structure
    // Use initial values for orphan age classes, otherwise use process model
    if(pop_year[i] <= max_ocean_age)
    {
      ii = (pop[i] - 1)*max_ocean_age + pop_year[i];
      q_orphan = append_row(sum(head(q_init[ii], N_age - N_orphan_age + 1)), 
                            tail(q_init[ii], N_orphan_age - 1));
    }
    
    for(a in 1:N_age)
    {
      if(pop_year[i] <= ocean_ages[a]) // use initial values
        S_W_a[a] = S_init[ii]*(1 - p_HOS_all[i])*q_orphan[a - (N_age - N_orphan_age)];
      else // use recruitment process model
        S_W_a[a] = M[i-ocean_ages[a]]*s_MS[i-ocean_ages[a]]*p[i-ocean_ages[a],a] * 
                   (1 - age_F[a]*F_rate[i])*(1 - age_B[a]*B_rate_all[i]);
    }
    
    // Total spawners and age structure
    S_W[i] = sum(S_W_a);
    S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
    S[i] = S_W[i] + S_H[i];
    q[i,] = S_W_a/S_W[i];

    // Smolt production from brood year i
    M_hat[i] = SR(SR_fun, RRS, alpha_Xbeta[i], alpha_W_Xbeta[i], alpha_H_Xbeta[i],
                  Mmax_Xbeta[i], Mmax_W_Xbeta[i], Mmax_H_Xbeta[i], S[i], S_W[i], S_H[i], A[i]);
    M0[i] = M_hat[i] * exp(eta_year_M[year[i]] + dot_product(X_M[i], beta_M) + sigma_M*zeta_M[i]);
  }
}

model {
  vector[N_B] log_B_take; // log of true broodstock take when B_take_obs > 0
  
  // Priors
  
  // smolt recruitment
  // [log(alpha), log(Mmax)] ~ MVN([mu_alpha, mu_Mmax], D*R_alphaMmax*D)
  // where D = diag_matrix(sigma_alpha, sigma_Mmax)
  [mu_alpha, mu_alpha_W, mu_alpha_H] ~ normal(prior_mu_alpha_mean, prior_mu_alpha_sd);
  beta_alpha ~ normal(0,5);
  sigma_alpha ~ normal(0,3);
  append_row(zeta_alpha, append_row(zeta_alpha_W, zeta_alpha_H)) ~ std_normal();  
  [mu_Mmax, mu_Mmax_W, mu_Mmax_H] ~ normal(prior_mu_Mmax_mean, prior_mu_Mmax_sd);
  beta_Mmax ~ normal(0,5);
  sigma_Mmax ~ normal(0,3);
  append_row(zeta_Mmax, append_row(zeta_Mmax_W, zeta_Mmax_H)) ~ std_normal();
  L_alphaMmax ~ lkj_corr_cholesky(1);
  beta_M ~ normal(0,5);
  rho_M ~ gnormal(0,0.85,30);  // mildly regularize to ensure stationarity
  sigma_year_M ~ normal(0,3);
  zeta_year_M ~ std_normal();  // eta_year_M[i] ~ N(rho_M*eta_year_M[i-1], sigma_year_M)
  sigma_M ~ normal(0,3);
  zeta_M ~ std_normal();       // total recruits: M ~ lognormal(log(M_hat), sigma)

  // SAR
  mu_MS ~ beta(prior_mu_MS[1], prior_mu_MS[2]);
  beta_MS ~ normal(0,3);
  rho_MS ~ gnormal(0,0.85,20);    // mildly regularize rho to ensure stationarity
  sigma_year_MS ~ normal(0,3);
  zeta_year_MS ~ std_normal(); // eta_year_MS[i] ~ N(rho_MS*eta_year_MS[i-1], sigma_year_MS)
  sigma_MS ~ normal(0,3);
  zeta_MS ~ std_normal();      // SAR: logit(s_MS) ~ normal(logit(s_MS_hat), sigma_MS)

  // spawner age structure
  mu_p ~ dirichlet(prior_mu_p);
  to_vector(sigma_pop_p) ~ normal(0,3);
  to_vector(sigma_p) ~ normal(0,3);
  L_pop_p ~ lkj_corr_cholesky(1); // R_pop_p ~ lkj_corr(1)
  L_p ~ lkj_corr_cholesky(1);     // R_p ~ lkj_corr(1)
  // pop mean age probs logistic MVN: 
  // mu_pop_alr_p[i,] ~ MVN(mu_alr_p,D*R_pop_p*D), where D = diag_matrix(sigma_pop_p)
  to_vector(zeta_pop_p) ~ std_normal();
  // age probs logistic MVN: 
  // alr_p[i,] ~ MVN(mu_pop_alr_p[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
  to_vector(zeta_p) ~ std_normal();

  // removals
  log_B_take = log(S_W[which_B]) + log(q[which_B,]*age_B) + logit(B_rate); 
  // implies B_take[i] = S_W[i] * (q[i,] * ageB) * B_rate[i] / (1 - B_rate[i])
  B_take_obs ~ lognormal(log_B_take, 0.05); // penalty to force pred and obs broodstock take to match 

  // initial states
  // (accounting for amalgamation of q_init to q_orphan)
  M_init ~ lognormal(mu_M_init, sigma_M_init);
  S_init ~ lognormal(mu_S_init, sigma_S_init);
  {
    matrix[N_age,max_ocean_age*N_pop] q_init_mat;
    
    for(j in 1:size(q_init)) q_init_mat[,j] = q_init[j];
    target += sum((mu_q_init - 1) .* log(q_init_mat)); // q_init[i] ~ Dir(mu_q_init[,i])
  }

  // observation error
  tau_M ~ gnormal(prior_tau_M[1], prior_tau_M[2], prior_tau_M[3]);
  tau_S ~ gnormal(prior_tau_S[1], prior_tau_S[2], prior_tau_S[3]);

  // Observation model
  M_obs[which_M_obs] ~ lognormal(log(M[which_M_obs]), tau_M);  // observed smolts
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau_S);  // observed spawners
  n_H_obs ~ binomial(n_HW_obs, p_HOS); // observed counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q));  // obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
}

generated quantities {
  real delta_mu_alpha;              // H vs. W discount in hyper-mean log intrinsic productivity
  vector[RRS[1]*N_pop] delta_alpha; // H vs. W discount in log intrinsic productivity
  real delta_mu_Mmax;               // H vs. W discount in hyper-mean log maximum smolt recruitment
  vector[RRS[2]*N_pop] delta_Mmax;  // H vs. W discount in log maximum smolt recruitment
  real rho_alphaMmax;               // correlation between log(alpha) and log(Mmax)
  corr_matrix[N_SR] R_alphaMmax;    // correlation among log alpha, Mmax (all/W/H)
  matrix[N_pop,N_age] mu_pop_p; // population mean age distributions
  corr_matrix[N_age-1] R_pop_p; // among-pop correlation matrix of mean log-ratio age distns 
  corr_matrix[N_age-1] R_p;     // correlation matrix of within-pop cohort log-ratio age distns 
  vector[N] LL_M_obs;           // pointwise log-likelihood of smolts
  vector[N] LL_S_obs;           // pointwise log-likelihood of spawners
  vector[N_H] LL_n_H_obs;       // pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL_n_age_obs;       // pointwise log-likelihood of wild age frequencies
  vector[N] LL;                 // total pointwise log-likelihood                              
  
  // H vs. W discounts
  delta_mu_alpha = mu_alpha_H - mu_alpha_W;
  delta_alpha = log(alpha_H) - log(alpha_W);
  delta_mu_Mmax = mu_Mmax_H - mu_Mmax_W;
  delta_Mmax = log(Mmax_H) - log(Mmax_W);
  
  // Correlation matrix of log([alpha, alpha_W, alpha_H, Rmax, Rmax_W, Rmax_H])
  R_alphaMmax = multiply_lower_tri_self_transpose(L_alphaMmax);
  if(!max(RRS)) rho_alphaMmax = R_alphaMmax[2,1];
  
  // Population mean age distributions
  {
    matrix[N_pop,N_age] exp_pop_alr_p = append_col(exp(mu_pop_alr_p), rep_vector(1,N_pop));
    mu_pop_p = diag_pre_multiply(inv(row_sums(exp_pop_alr_p)), exp_pop_alr_p);
  }

  // Correlation matrices of log-ratio conditional age-at-return distributions
  R_pop_p = multiply_lower_tri_self_transpose(L_pop_p);
  R_p = multiply_lower_tri_self_transpose(L_p);
  
  // Pointwise observation log-likelihoods
  LL_M_obs = rep_vector(0,N);
  for(i in 1:N_M_obs)
    LL_M_obs[which_M_obs[i]] = lognormal_lpdf(M_obs[which_M_obs[i]] | log(M[which_M_obs[i]]), tau_M); 
  LL_S_obs = rep_vector(0,N);
  for(i in 1:N_S_obs)
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau_S); 
  LL_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1,N_age);
  LL_n_H_obs = rep_vector(0,N_H);
  for(i in 1:N_H)
    LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
  LL = LL_M_obs + LL_S_obs + LL_n_age_obs;
  LL[which_H] = LL[which_H] + LL_n_H_obs;
}
