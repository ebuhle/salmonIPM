functions {
  #include /include/SR.stan
  #include /include/gnormal_lpdf.stan
  #include /include/mat_lmult.stan
  #include /include/row_sums.stan
  #include /include/veq.stan
}

data {
  // info for observed data
  int<lower=1> N;                      // total number of cases in all natural pops and years
  array[N] int<lower=1,upper=N> pop;   // natural population ID
  array[N] int<lower=1,upper=N> year;  // natural calendar year ID
  // fecundity
  int<lower=0> N_E;                    // number of fecundity observations
  array[N_E] int<lower=1> age_E;       // female age (years) for fecundity observations
  vector<lower=0>[N_E] E_obs;          // observed fecundity
  // spawner-smolt productivity
  int<lower=1> SR_fun;                 // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  array[2] int<lower=0,upper=1> RRS;   // fit W vs. H {alpha, Mmax} (1) or not (0)?
  vector[N] A;                         // habitat area associated with each spawner abundance obs
  int<lower=1> smolt_age;              // smolt age
  int<lower=0> K_psi;                  // number of egg-to-smolt survival covariates
  matrix[N,K_psi] X_psi;               // egg-to-smolt survival covariates
  array[2] real prior_mu_psi;          // prior a, b for hyper-mean egg-to-smolt survival
  array[2] real prior_mu_psi_W;        // prior a, b for hyper-mean W egg-to-smolt survival
  array[2] real prior_mu_psi_H;        // prior a, b for hyper-mean H egg-to-smolt survival
  int<lower=0> K_Mmax;                 // number of maximum smolt recruitment covariates
  matrix[N,K_Mmax] X_Mmax;             // maximum smolt recruitment covariates
  array[2] real prior_mu_Mmax;         // prior mean, sd for hyper-mean log maximum smolt recruitment
  array[2] real prior_mu_Mmax_W;       // prior mean, sd for hyper-mean W log maximum smolt recruitment
  array[2] real prior_mu_Mmax_H;       // prior mean, sd for hyper-mean H log maximum smolt recruitment
  int<lower=0> K_M;                    // number of smolt recruitment covariates
  array[N] row_vector[K_M] X_M;        // smolt recruitment covariates
  // smolt abundance and observation error
  int<lower=1,upper=N> N_M_obs;        // number of cases with non-missing smolt abundance obs 
  array[N_M_obs] int<lower=1,upper=N> which_M_obs; // cases with non-missing smolt abundance obs
  vector<lower=0>[N] M_obs;            // observed annual natural smolt abundance (not density)
  int<lower=1,upper=N> N_tau_M_obs;    // number of cases with known smolt observation error SD
  array[N_tau_M_obs] int<lower=1,upper=N> which_tau_M_obs; // cases with known smolt observation error SD
  vector<lower=0>[N] tau_M_obs;        // known smolt observation error SDs
  int<lower=0,upper=N> N_upstream;     // number of smolt obs upstream of a trap in another pop
  array[N_upstream] int<lower=1,upper=N> which_upstream;  // smolt obs upstream of a trap in another pop
  array[N_upstream] int<lower=1,upper=N> downstream_trap; // row index for the downstream trap
  // SAR (smolt-spawner survival)
  int<lower=0> K_MS;                   // number of SAR covariates
  matrix[N,K_MS] X_MS;                 // SAR covariates
  array[2] real prior_mu_MS;           // prior a, b for mean SAR
  // spawner abundance and observation error
  int<lower=1,upper=N> N_S_obs;        // number of cases with non-missing spawner abundance obs 
  array[N_S_obs] int<lower=1,upper=N> which_S_obs; // cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;            // observed annual total spawner abundance (not density)
  int<lower=1,upper=N> N_tau_S_obs;    // number of cases with known spawner observation error SD
  array[N_tau_S_obs] int<lower=1,upper=N> which_tau_S_obs; // cases with known spawner observation error SD
  vector<lower=0>[N] tau_S_obs;        // known spawner observation error SDs
  // spawner age structure and sex ratio
  int<lower=2> N_age;                  // number of adult age classes
  int<lower=2> max_age;                // maximum adult age
  matrix<lower=0>[N,N_age] n_age_obs;  // observed wild spawner age frequencies (all zero row = NA)  
  vector<lower=0>[N_age] prior_mu_p;   // prior concentration for mean age distribution
  array[N] int<lower=0> n_M_obs;       // observed count of male spawners
  array[N] int<lower=0> n_F_obs;       // observed count of female spawners
  vector<lower=0,upper=1>[N] p_G_obs;  // proportion green (fully fecund) females
  // origin composition
  int<lower=0> N_H_pop;                // number of hatchery pops
  array[N_H_pop] int<lower=1,upper=max(pop)> which_H_pop; // hatchery pop IDs
  matrix<lower=0>[N,1+N_H_pop] n_O_obs; // observed spawner origin frequencies (1st col is NOR)
  // fishery and hatchery removals and translocations
  vector[N] F_rate;                    // fishing mortality rate of wild adults
  vector<lower=0,upper=1>[N_age] age_F; // is age a (non)selected (0/1) by fishery?
  int<lower=0,upper=N> N_B;            // number of years with B_take > 0
  array[N_B] int<lower=1,upper=N> which_B; // years with B_take > 0
  vector[N_B] B_take_obs;              // observed broodstock take of wild adults
  vector<lower=0,upper=1>[N_age] age_B; // is age a (non)selected (0/1) in broodstock?
  vector<lower=0>[N] S_add_obs;        // number of translocated spawners added to population
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years
  array[N] int<lower=1> pop_year;          // index of years within each pop, starting at 1
  int<lower=0,upper=N_pop> N_W_pop = N_pop - N_H_pop;  // number of natural (wild) pops
  array[N_W_pop] int<lower=1,upper=N_pop> which_W_pop; // wild pop IDs
  vector<lower=0,upper=1>[N] indx_H;       // is case from a hatchery (1) or wild (0) pop?
  array[N_age] int<lower=0> ocean_ages;    // ocean ages
  vector[N_age] ones_N_age = rep_vector(1,N_age); // for rowsums of p matrix 
  vector[N] ones_N = rep_vector(1,N);      // for elementwise inverse of rowsums 
  int<lower=1> max_ocean_age = max_age - smolt_age; // maximum ocean age
  int<lower=1> min_ocean_age = max_ocean_age - N_age + 1; // minimum ocean age
  int<lower=1> min_age = max_age - N_age + 1; // minimum adult age
  array[N_E] int<lower=1> a_E;             // female age index for fecundity observations
  array[N] int<lower=0> n_MF_obs;          // total sample sizes for sex frequencies
  vector<lower=0,upper=1>[N] p_NG_obs = 1 - p_G_obs; // proportion non-green females
  real mu_mu_Mmax = max(log(M_obs[which_M_obs] ./ A[which_M_obs])); // prior mean of mu_Mmax
  real sigma_mu_Mmax = sd(log(M_obs[which_M_obs])); // prior SD of mu_Mmax
  real mu_M_init = mean(log(M_obs[which_M_obs]));  // prior log-mean of smolt abundance in years 1:smolt_age
  real sigma_M_init = 2*sd(log(M_obs[which_M_obs])); // prior log-SD of smolt abundance in years 1:smolt_age
  vector[max_ocean_age*N_pop] mu_S_init;   // prior mean of total spawner abundance in years 1:max_ocean_age
  real sigma_S_init = sd(log(S_obs[which_S_obs])); // prior log-SD of spawner abundance in years 1:max_ocean_age
  matrix[N_age,max_ocean_age*N_pop] mu_q_init; // prior counts of wild spawner age distns in years 1:max_ocean_age
  int<lower=0,upper=1> any_RRS = max(RRS); // does either S-R parameter differ by W vs. H?
  int<lower=2,upper=4> N_SR = 2 + sum(RRS); // number of S-R parameter hyper-means
  row_vector[3] prior_mu_psi_a;            // prior a for hyper-mean egg-to-smolt survival (all/W/H)
  row_vector[3] prior_mu_psi_b;            // prior b for hyper-mean egg-to-smolt survival (all/W/H)
  row_vector[3] prior_mu_Mmax_mean;        // prior means for hyper-mean log maximum recruitment (all/W/H)
  row_vector[3] prior_mu_Mmax_sd;          // prior SDs for hyper-mean log maximum recruitment (all/W/H)
  
  {
    int jj = 1;
    for(j in 1:N_pop) 
    {
      if(max(veq(which_H_pop, j)) == 0)
      {
        which_W_pop[jj] = j;
        jj = jj + 1;
      }
    }
  }

  for(i in 1:N) indx_H[i] = max(veq(which_H_pop, pop[i]));
  for(a in 1:N_age) ocean_ages[a] = min_ocean_age - 1 + a;
  for(i in 1:N) n_MF_obs[i] = n_M_obs[i] + n_F_obs[i];
  for(i in 1:N_E) a_E[i] = age_E[i] - min_age + 1;

  pop_year[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year[i] = 1;
    else
      pop_year[i] = pop_year[i-1] + 1;
  }
  
  for(i in 1:max_ocean_age)
  {
    int N_orphan_age = N_age - max(i - min_ocean_age, 0); // number of orphan age classes
    int N_amalg_age = N_age - N_orphan_age + 1;           // number of amalgamated age classes
    
    for(j in 1:N_pop)
    {
      int ii = (j - 1)*max_ocean_age + i; // index into S_init, q_init
      
      // S_init prior mean that scales observed log-mean by fraction of orphan age classes
      mu_S_init[ii] = mean(log(S_obs[which_S_obs])) + log(N_orphan_age) - log(N_age);
      
      // prior on q_init that implies q_orphan ~ Dir(1)
      mu_q_init[,ii] = append_row(rep_vector(1.0/N_amalg_age, N_amalg_age), 
                                  rep_vector(1, N_orphan_age - 1));
    }
  }
  
  prior_mu_psi_a = [prior_mu_psi[1], prior_mu_psi_W[1], prior_mu_psi_H[1]];
  prior_mu_psi_b = [prior_mu_psi[2], prior_mu_psi_W[2], prior_mu_psi_H[2]];
  prior_mu_Mmax_mean = [prior_mu_Mmax[1], prior_mu_Mmax_W[1], prior_mu_Mmax_H[1]];
  prior_mu_Mmax_sd = [prior_mu_Mmax[2], prior_mu_Mmax_W[2], prior_mu_Mmax_H[2]];
}

parameters {
  // fecundity
  vector<lower=0>[N_age] mu_E;           // mean fecundity at each female age
  vector<lower=0>[N_age] sigma_E;        // among-female SD of fecundity at each age
  real<lower=0,upper=1> delta_NG;        // fecundity of non-green females as proportion of green
  // spawner-smolt productivity
  real<lower=0,upper=1> mu_psi;          // hyper-mean density-independent egg-smolt survival
  real<lower=0,upper=1> mu_psi_W;        // hyper-mean W density-independent egg-smolt survival
  real<lower=0,upper=1> mu_psi_H;        // hyper-mean H density-independent egg-smolt survival
  vector[K_psi] beta_psi;                // regression coefs for logit psi
  real<lower=0> sigma_psi;               // hyper-SD logit density-independent egg-smolt survival
  vector[!RRS[1]*N_pop] zeta_psi;        // logit density-independent egg-smolt survival (Z-scores)
  vector[RRS[1]*N_pop] zeta_psi_W;       // logit W density-independent egg-smolt survival (Z-scores)
  vector[RRS[1]*N_pop] zeta_psi_H;       // logit H density-independent egg-smolt survival (Z-scores)
  real mu_Mmax;                          // hyper-mean log maximum smolt recruitment
  real mu_Mmax_W;                        // hyper-mean log W maximum smolt recruitment
  real mu_Mmax_H;                        // hyper-mean log H maximum smolt recruitment
  vector[K_Mmax] beta_Mmax;              // regression coefs for log Mmax
  real<lower=0> sigma_Mmax;              // among-pop hyper-SD of log maximum smolt recruitment
  vector[!RRS[2]*N_pop] zeta_Mmax;       // log maximum smolt recruitment (Z-scores)
  vector[RRS[2]*N_pop] zeta_Mmax_W;      // log W maximum smolt recruitment (Z-scores)
  vector[RRS[2]*N_pop] zeta_Mmax_H;      // log H maximum smolt recruitment (Z-scores)
  vector[K_M] beta_M;                    // regression coefs for annual spawner-smolt productivity anomalies
  real<lower=-1,upper=1> rho_M;          // AR(1) coef of annual spawner-smolt productivity anomalies
  real<lower=0> sigma_year_M;            // innovation SD of annual spawner-smolt productivity anomalies
  vector[N_year] zeta_year_M;            // annual spawner-smolt productivity anomalies (Z-scores)
  real<lower=0> sigma_M;                 // SD of unique spawner-smolt productivity process errors
  vector[N] zeta_M;                      // unique spawner-smolt productivity process errors (Z-scores)
  // SAR
  real<lower=0,upper=1> mu_MS;           // hyper-mean SAR
  vector[K_MS] beta_MS;                  // regression coefs for annual logit SAR anomalies
  real<lower=-1,upper=1> rho_MS;         // AR(1) coef of annual logit SAR anomalies
  real<lower=0> sigma_year_MS;           // innovation SD of annual logit SAR anomalies
  vector[N_year] zeta_year_MS;           // annual logit SAR anomalies (Z-scores)
  real<lower=0> sigma_MS;                // SD of unique logit SAR process errors
  vector[N] zeta_MS;                     // unique logit SAR process errors (Z-scores)
  // spawner age structure
  simplex[N_age] mu_p;                   // among-pop mean of age distributions
  vector<lower=0>[N_age-1] sigma_pop_p;  // among-pop SD of mean log-ratio age distributions
  cholesky_factor_corr[N_age-1] L_pop_p; // Cholesky factor of among-pop correlation matrix of mean log-ratio age distns
  matrix[N_pop,N_age-1] zeta_pop_p;      // population mean log-ratio age distributions (Z-scores)
  vector<lower=0>[N_age-1] sigma_p;      // SD of log-ratio cohort age distributions
  cholesky_factor_corr[N_age-1] L_p;     // Cholesky factor of correlation matrix of cohort log-ratio age distributions
  matrix[N,N_age-1] zeta_p;              // log-ratio cohort age distributions (Z-scores)
  // spawner sex ratio
  real<lower=0,upper=1> mu_F;            // hyper-mean proportion female spawners 
  real<lower=0> sigma_pop_F;             // among-pop SD of mean logit proportion female
  vector[N_pop] zeta_pop_F;              // population mean logit proportion female (Z-scores)
  real<lower=0> sigma_F;                 // annual SD of logit proportion female
  vector[N] zeta_F;                      // logit proportion females by outmigration year (Z-scores)
  // origin composition, removals
  array[N_H_pop] simplex[N_W_pop] P_D;   // P(dispersal from hatchery h to wild pop j) 
  vector<lower=0,upper=1>[N_B] B_rate;   // true broodstock take rate when B_take > 0
  // initial states and observation error
  vector<lower=0>[smolt_age*N_pop] M_init; // true smolt abundance in years 1:smolt_age
  vector<lower=0>[max_ocean_age*N_pop] S_init; // true total spawner abundance in years 1:max_ocean_age
  array[max_ocean_age*N_pop] simplex[N_age] q_init; // true wild spawner age distributions in years 1:max_ocean_age
  vector<lower=0,upper=1>[max_ocean_age*N_pop] q_F_init; // true proportion females in years 1:max_ocean_age
  real<lower=0> mu_tau_M;                // median smolt observation error SD
  real<lower=0> sigma_tau_M;             // log-SD of smolt observation error SDs
  vector[N] tau_M_z;                     // log of smolt observation error SDs (Z-scores)
  real<lower=0> mu_tau_S;                // median spawner observation error SD
  real<lower=0> sigma_tau_S;             // log-SD of spawner observation error SDs
  vector[N] tau_S_z;                     // log of spawner observation error SDs (Z-scores)
}

transformed parameters {
  // fecundity
  vector[N] E_hat = rep_vector(0,N);     // expected egg deposition
  vector[N] E_hat_W = rep_vector(0,N);   // expected W egg deposition
  vector[N] E_hat_H = rep_vector(0,N);   // expected H egg deposition
  // spawner-smolt productivity
  vector<lower=0,upper=1>[!RRS[1]*N_pop] psi;  // density-independent egg-smolt survival 
  vector<lower=0,upper=1>[RRS[1]*N_pop] psi_W; // W density-independent egg-smolt survival 
  vector<lower=0,upper=1>[RRS[1]*N_pop] psi_H; // H density-independent egg-smolt survival 
  vector[N] psi_Xbeta;                         // egg-smolt survival including covariate effects
  vector[N] psi_W_Xbeta;                       // W egg-smolt survival including covariate effects
  vector[N] psi_H_Xbeta;                       // H egg-smolt survival including covariate effects
  vector<lower=0>[!RRS[2]*N_pop] Mmax;   // maximum smolt recruitment 
  vector<lower=0>[RRS[2]*N_pop] Mmax_W;  // W maximum smolt recruitment 
  vector<lower=0>[RRS[2]*N_pop] Mmax_H;  // H maximum smolt recruitment 
  vector[N] Mmax_Xbeta;                  // maximum recruitment including covariate effects
  vector[N] Mmax_W_Xbeta;                // W maximum recruitment including covariate effects
  vector[N] Mmax_H_Xbeta;                // H maximum recruitment including covariate effects
  vector[N_year] eta_year_M;             // annual spawner-smolt productivity anomalies
  vector<lower=0>[N] M_hat = rep_vector(0,N); // expected smolt abundance by brood year
  vector<lower=0>[N] M0 = rep_vector(0,N);    // true smolt abundance by brood year
  vector<lower=0>[N] M = rep_vector(0,N);     // true smolt abundance by outmigration year
  // SAR
  vector[N_year] eta_year_MS;            // annual logit SAR anomalies
  vector<lower=0,upper=1>[N] s_MS;       // true SAR
  // H/W spawner abundance, removals
  vector<lower=0>[N] S_W = rep_vector(0,N); // true total wild spawner abundance
  vector<lower=0>[N] S_H = rep_vector(0,N); // true total hatchery-origin spawner abundance
  vector<lower=0>[N] S = rep_vector(0,N);   // true total spawner abundance
  array[N_pop] vector[N_W_pop] P_D_all;     // dispersal matrix padded with zeros 
  array[N_year] matrix[N_pop,N_pop] S_O = 
    rep_array(rep_matrix(0,N_pop,N_pop), N_year); // true spawners by origin and return location
  matrix[N,1+N_H_pop] q_O = rep_matrix(1,N,1+N_H_pop); // true origin distns (col 1: unknown / natural)
  vector<lower=0,upper=1>[N] B_rate_all; // broodstock take rate in all years
  // spawner age structure
  row_vector[N_age-1] mu_alr_p;          // mean of log-ratio cohort age distributions
  matrix[N_pop,N_age-1] mu_pop_alr_p;    // population mean log-ratio age distributions
  matrix<lower=0,upper=1>[N,N_age] p;    // true adult age distributions by outmigration year
  matrix<lower=0,upper=1>[N,N_age] q = rep_matrix(1,N,N_age); // true spawner age distributions
  // spawner sex ratio
  vector<lower=0,upper=1>[N] p_F;        // proportion females by outmigration year
  vector<lower=0,upper=1>[N] q_F = rep_vector(0,N); // proportion females by return year
  // observation error SDs
  vector<lower=0>[N] tau_M;              // smolt observation error SDs
  vector<lower=0>[N] tau_S;              // spawner observation error SDs

  // Pad B_rate
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;
  
  // Matt trick for population-specific egg-smolt survival and capacity (uncorrelated)
  {
    vector[N] Xbeta_psi = mat_lmult(X_psi, beta_psi);          // covariate effects on psi 
    vector[N] Xbeta_Mmax = exp(mat_lmult(X_Mmax, beta_Mmax));  // covariate effects on Mmax

    if(RRS[1]) {
      vector[N_pop] logit_psi_W = logit(mu_psi_W) + sigma_psi*zeta_psi_W;
      vector[N_pop] logit_psi_H = logit(mu_psi_H) + sigma_psi*zeta_psi_H;
      psi_W = rep_vector(0,N_pop);
      psi_W[which_W_pop] = inv_logit(logit_psi_W[which_W_pop]);
      psi_W_Xbeta = inv_logit(logit_psi_W[pop] + Xbeta_psi) .* (1 - indx_H);
      psi_H = rep_vector(0,N_pop);
      psi_H[which_W_pop] = inv_logit(logit_psi_H[which_W_pop]);
      psi_H_Xbeta = inv_logit(logit_psi_H[pop] + Xbeta_psi) .* (1 - indx_H);
    } else {
      vector[N_pop] logit_psi = logit(mu_psi) + sigma_psi*zeta_psi;
      psi = rep_vector(0,N_pop);
      psi[which_W_pop] = inv_logit(logit_psi[which_W_pop]);
      psi_Xbeta = inv_logit(logit_psi[pop] + Xbeta_psi) .* (1 - indx_H);
    }
    
    if(RRS[2]) {
      Mmax_W = rep_vector(0,N_pop);
      Mmax_W[which_W_pop] = exp(mu_Mmax_W + sigma_Mmax*zeta_Mmax_W[which_W_pop]);
      Mmax_W_Xbeta = Mmax_W[pop] .* Xbeta_Mmax;
      Mmax_H = rep_vector(0,N_pop);
      Mmax_H[which_W_pop] = exp(mu_Mmax_H + sigma_Mmax*zeta_Mmax_H[which_W_pop]);
      Mmax_H_Xbeta = Mmax_H[pop] .* Xbeta_Mmax;
    } else {
      Mmax = rep_vector(0,N_pop);
      Mmax[which_W_pop] = exp(mu_Mmax + sigma_Mmax*zeta_Mmax[which_W_pop]);
      Mmax_Xbeta = Mmax[pop] .* Xbeta_Mmax;
    }
  }
  
  // AR(1) model for spawner-smolt productivity and SAR anomalies
  eta_year_M[1] = zeta_year_M[1] * sigma_year_M / sqrt(1 - rho_M^2);     // initial anomaly
  eta_year_MS[1] = zeta_year_MS[1] * sigma_year_MS / sqrt(1 - rho_MS^2); // initial anomaly
  for(i in 2:N_year)
  {
    eta_year_M[i] = rho_M * eta_year_M[i-1] + zeta_year_M[i] * sigma_year_M;
    eta_year_MS[i] = rho_MS * eta_year_MS[i-1] + zeta_year_MS[i] * sigma_year_MS;
  }
  // constrain "fitted" anomalies to sum to 0 (columns of X should be centered)
  eta_year_M = eta_year_M - mean(eta_year_M);
  eta_year_MS = eta_year_MS - mean(eta_year_MS);
  // annual population-specific SAR
  s_MS = inv_logit(logit(mu_MS) + mat_lmult(X_MS, beta_MS) + eta_year_MS[year] + zeta_MS*sigma_MS);
  
  // Multivariate Matt trick for age vectors (pop-specific mean and within-pop, time-varying)
  mu_alr_p = to_row_vector(log(head(mu_p, N_age-1)) - log(mu_p[N_age]));
  mu_pop_alr_p = rep_matrix(mu_alr_p,N_pop) + diag_pre_multiply(sigma_pop_p, L_pop_p * zeta_pop_p')';
  // Inverse log-ratio (softmax) transform of cohort age distn
  {
    matrix[N,N_age-1] alr_p = mu_pop_alr_p[pop,] + diag_pre_multiply(sigma_p, L_p * zeta_p')';
    matrix[N,N_age] exp_alr_p = append_col(exp(alr_p), ones_N);
    p = diag_pre_multiply(ones_N ./ (exp_alr_p * ones_N_age), exp_alr_p);
  }
  
  // Annual population-specific proportion females
  p_F = inv_logit(logit(mu_F) + sigma_pop_F*zeta_pop_F[pop] + sigma_F*zeta_F);
  
  // Pad straying matrix
  P_D_all = rep_array(rep_vector(0, N_W_pop), N_pop);
  P_D_all[which_H_pop] = P_D;

  // Calculate true hatchery spawners by origin and return location in return year i
  // so they are available for the wild pop loop 
  for(i in 1:N)
  {
    // This loop is only for hatchery populations
    if(indx_H[i])
    {
      row_vector[N_age] S_H_a;  // hatchery spawners by age
      
      //// NOTE: should use orphan formulation 
      if(pop_year[i] <= max_ocean_age)  // use initial values
      {
        int ii = (pop[i] - 1)*max_ocean_age + pop_year[i]; // index into S_init
        S_O[year[i]][which_W_pop,pop[i]] = S_init[ii]*P_D_all[pop[i]];
      }
      else // use hatchery smolt-to-adult process model 
      {
        for(a in 1:N_age)
          S_H_a[a] = M_obs[i-ocean_ages[a]] * s_MS[i-ocean_ages[a]] * p[i-ocean_ages[a],a];
        S_O[year[i]][which_W_pop,pop[i]] = S_H_a * (1 - age_F*F_rate[i]) * P_D_all[pop[i]];
      }
    }
  }

  // Calculate true total wild and hatchery spawners, spawner age distribution, 
  // eggs and smolts, and predict smolt recruitment from brood year i
  for(i in 1:N)
  {
    // This loop is only for wild populations
    if(!indx_H[i])
    {
      int ii;                        // index into S_init and q_init
      // number of orphan age classes <lower=0,upper=N_age>
      int N_orphan_age = max(N_age - to_int(fdim(pop_year[i], min_ocean_age)), N_age); 
      vector[N_orphan_age] q_orphan; // orphan age distribution
      row_vector[N_age] S_W_a;       // wild spawners by age
      vector[N_age] q_F_a;           // proportion of each age that are female
      
      // Smolt recruitment
      if(pop_year[i] <= smolt_age)
        M[i] = M_init[(pop[i]-1)*smolt_age + pop_year[i]];  // use initial values
      else
        M[i] = M0[i-smolt_age];  // smolts from appropriate brood year
      
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
        {
          S_W_a[a] = S_init[ii]*q_orphan[a - (N_age - N_orphan_age)];
          q_F_a[a] = q_F_init[ii];
        }
        else // use recruitment process model 
        {
          S_W_a[a] = M[i-ocean_ages[a]] * s_MS[i-ocean_ages[a]] * p[i-ocean_ages[a],a] *
                     (1 - age_F[a]*F_rate[i]) * (1 - age_B[a]*B_rate_all[i]);
          q_F_a[a] = p_F[i-ocean_ages[a]];
        }
      }
      
      // Total spawners including translocations and strays
      // age and sex structure
      S_W[i] = sum(S_W_a);
      q[i,] = S_W_a/S_W[i];
      q_F[i] = q[i,]*q_F_a;
      S_H[i] = sum(S_O[year[i]][pop[i],which_H_pop]);
      S[i] = S_W[i] + S_H[i] + S_add_obs[i];
      q_O[i,1] = S_W[i] + S_add_obs[i];             // unknown natural origin
      q_O[i,2:] = S_O[year[i]][pop[i],which_H_pop]; // known (hatchery) origin
      q_O[i,] = q_O[i,]/sum(q_O[i,]);

      // Expected egg production from brood year i
      // weighted by age structure and sex ratio 
      // discounted for proportion of non-green (not fully fecund) females
      {
        // weighted per capita fecundity
        real wt_E = q[i,] * mu_E * q_F[i] * (p_G_obs[i] + delta_NG * p_NG_obs[i]);
        
        if(any_RRS) {
          E_hat_W[i] = wt_E * S_W[i];
          E_hat_H[i] = wt_E * S_H[i];
        } else {
          E_hat[i] = wt_E * S[i];
        }
      }
      
      // Smolt production from brood year i
      M_hat[i] = SR(SR_fun, RRS, psi_Xbeta[i], psi_W_Xbeta[i], psi_H_Xbeta[i],
                    Mmax_Xbeta[i], Mmax_W_Xbeta[i], Mmax_H_Xbeta[i], 
                    E_hat[i], E_hat_W[i], E_hat_H[i], A[i]);
      M0[i] = M_hat[i] * exp(eta_year_M[year[i]] + dot_product(X_M[i], beta_M) + sigma_M*zeta_M[i]);
    }
  }
  
  // Impute missing observation error SDs from lognormal hyperdistribution
  // (use Matt trick)
  tau_M = exp(log(mu_tau_M) + sigma_tau_M*tau_M_z);
  tau_M[which_tau_M_obs] = tau_M_obs[which_tau_M_obs];
  tau_S = exp(log(mu_tau_S) + sigma_tau_S*tau_S_z);
  tau_S[which_tau_S_obs] = tau_S_obs[which_tau_S_obs];
}

model {
  vector[N_B] log_B_take; // log of true broodstock take when B_take_obs > 0
  vector[N] M_downstream; // total smolts including upstream populations 

  // Priors
  
  // fecundity and sex ratio
  mu_E ~ normal(2500,500);
  sigma_E ~ normal(500,1000);

  // spawner-smolt productivity
  [mu_psi, mu_psi_W, mu_psi_H] ~ beta(prior_mu_psi_a, prior_mu_psi_b);
  beta_psi ~ normal(0,5);
  sigma_psi ~ normal(0,1);     // mildly regularize long right tail
  // logit(psi) ~ N(logit(mu_psi), sigma_psi)
  append_row(zeta_psi, append_row(zeta_psi_W, zeta_psi_H)) ~ std_normal();  
  [mu_Mmax, mu_Mmax_W, mu_Mmax_H] ~ normal(prior_mu_Mmax_mean, prior_mu_Mmax_sd);
  beta_Mmax ~ normal(0,5);
  sigma_Mmax ~ normal(0,3);
  // log(Mmax) ~ N(mu_Mmax, sigma_Mmax)
  append_row(zeta_Mmax, append_row(zeta_Mmax_W, zeta_Mmax_H)) ~ std_normal();
  beta_M ~ normal(0,5);
  rho_M ~ gnormal(0,0.85,20);  // regularize away from 1 to ensure stationarity
  sigma_year_M ~ normal(0,2);
  zeta_year_M ~ std_normal();  // eta_year_M[i] ~ N(rho_M*eta_year_M[i-1], sigma_year_M)
  sigma_M ~ normal(0,2);
  zeta_M ~ std_normal();       // M0 ~ lognormal(log(M_hat) + eta_year_M, sigma_M)
  
  // SAR
  mu_MS ~ beta(prior_mu_MS[1], prior_mu_MS[2]);
  beta_MS ~ normal(0,3);
  rho_MS ~ gnormal(0,0.85,20); // regularize away from 1 to ensure stationarity
  sigma_year_MS ~ normal(0,2);
  zeta_year_MS ~ std_normal(); // eta_year_MS[i] ~ N(rho_MS*eta_year_MS[i-1], sigma_year_MS)
  sigma_MS ~ normal(0,2);
  zeta_MS ~ std_normal();      // SAR: logit(s_MS) ~ N(logit(s_MS_hat), sigma_MS)

  // spawner age structure
  mu_p ~ dirichlet(prior_mu_p);
  to_vector(sigma_pop_p) ~ normal(0,3);
  to_vector(sigma_p) ~ normal(0,3);
  L_pop_p ~ lkj_corr_cholesky(1);
  L_p ~ lkj_corr_cholesky(1);
  // pop mean age probs logistic MVN: 
  // mu_pop_alr_p[i,] ~ MVN(mu_alr_p, D*R_pop_p*D), where D = diag_matrix(sigma_pop_p)
  to_vector(zeta_pop_p) ~ std_normal();
  // age probs logistic MVN: 
  // alr_p[i,] ~ MVN(mu_pop_alr_p[pop[i],], D*R_p*D), where D = diag_matrix(sigma_p)
  to_vector(zeta_p) ~ std_normal();
  
  // spawner sex ratio
  sigma_pop_F ~ normal(0,2);
  zeta_pop_F ~ std_normal();   // eta_pop_F ~ N(0, sigma_pop_F)
  sigma_F ~ normal(0,2);
  zeta_F ~ std_normal();       // logit(p_F[i]) ~ N(eta_pop_F[pop[i]], sigma_F)
  q_F_init ~ beta(3,3);        // mildly regularize initial states toward 0.5

  // removals
  log_B_take = log(S_W[which_B]) + log(q[which_B,]*age_B) + logit(B_rate); 
  // implies B_take[i] = S_W[i] * (q[i,] * ageB) * B_rate[i] / (1 - B_rate[i])
  B_take_obs ~ lognormal(log_B_take, 0.05); // penalty to force pred and obs broodstock to match 

  // initial states
  // (accounting for amalgamation of q_init to q_orphan)
  M_init ~ lognormal(mu_M_init, sigma_M_init);
  S_init ~ lognormal(mu_S_init, sigma_S_init);
  {
    matrix[N_age,max_ocean_age*N_pop] q_init_mat;
    
    for(j in 1:size(q_init)) q_init_mat[,j] = q_init[j];
    target += sum((mu_q_init - 1) .* log(q_init_mat)); // q_init[i] ~ Dir(mu_q_init[,i])
  }

  // observation error SDs
  // hyperparameters
  mu_tau_M ~ normal(0,1); // half-normal
  sigma_tau_M ~ normal(0,5);
  mu_tau_S ~ normal(0,1); // half-normal
  sigma_tau_S ~ normal(0,5);
  // unknown SDs
  tau_M_z ~ std_normal(); // tau_M[-which_tau_M_obs] ~ lognormal(log(mu_tau_M), sigma_tau_M);
  tau_S_z ~ std_normal(); // tau_S[-which_tau_S_obs] ~ lognormal(log(mu_tau_S), sigma_tau_S);
  // known SDs
  tau_M_obs[which_tau_M_obs] ~ lognormal(log(mu_tau_M), sigma_tau_M);
  tau_S_obs[which_tau_S_obs] ~ lognormal(log(mu_tau_S), sigma_tau_S);

  // Observation model
  E_obs ~ normal(mu_E[a_E], sigma_E[a_E]);              // observed fecundity
  target += -normal_lccdf(0 | mu_E[a_E], sigma_E[a_E]); // zero-truncated normal: T[0, ] not vectorized
  M_downstream = M;
  for(i in 1:N_upstream) M_downstream[downstream_trap[i]] += M[which_upstream[i]];
  // M_downstream[which_M_obs] ~ lognormal(log(M_obs[which_M_obs]), tau_M[which_M_obs]); // prior on smolts
  M_obs[which_M_obs] ~ lognormal(log(M_downstream[which_M_obs]), tau_M[which_M_obs]); // likelihood of smolts
  // S[which_S_obs] ~ lognormal(log(S_obs[which_S_obs]), tau_S[which_S_obs]); // prior on spawners
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau_S[which_S_obs]); // likelihood of spawners
  target += sum(n_age_obs .* log(q)); // obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
  target += sum(n_O_obs .* log(q_O)); // obs origin freq: n_O_obs[i] ~ multinomial(q_O[i])
  n_F_obs ~ binomial(n_MF_obs, q_F);  // observed counts of female vs. male spawners
}

generated quantities {
  real delta_mu_psi;               // H vs. W discount in hyper-mean density-independent egg-smolt survival
  vector[RRS[1]*N_pop] delta_psi;  // H vs. W discount in logit density-independent egg-smolt survival
  real delta_mu_Mmax;              // H vs. W discount in hyper-mean log maximum smolt recruitment
  vector[RRS[2]*N_pop] delta_Mmax; // H vs. W discount in log maximum smolt recruitment
  corr_matrix[N_age-1] R_pop_p; // among-pop correlation matrix of mean log-ratio age distns 
  corr_matrix[N_age-1] R_p;     // correlation matrix of within-pop cohort log-ratio age distns 
  vector[N] p_HOS = row_sums(q_O[,2:]); // true proportion hatchery-origin spawners
  // ?? vector[N_E] LL_E_obs; ??
  vector[N] M_downstream;       // total smolts including upstream populations 
  vector[N] LL_M_obs;           // pointwise log-likelihood of smolts
  vector[N] LL_S_obs;           // pointwise log-likelihood of spawners
  vector[N] LL_n_age_obs;       // pointwise log-likelihood of wild age frequencies
  vector[N] LL_n_O_obs;         // pointwise log-likelihood of origin-frequencies
  vector[N] LL_n_F_obs;         // pointwise log-likelihood of female vs. male frequencies
  vector[N] LL;                 // total pointwise log-likelihood     

  delta_mu_psi = mu_psi_H - mu_psi_W;
  delta_psi = logit(psi_H) - logit(psi_W);
  delta_mu_Mmax = mu_Mmax_H - mu_Mmax_W;
  delta_Mmax = log(Mmax_H) - log(Mmax_W);

  R_pop_p = multiply_lower_tri_self_transpose(L_pop_p);
  R_p = multiply_lower_tri_self_transpose(L_p);
  
  M_downstream = M;
  for(i in 1:N_upstream) M_downstream[downstream_trap[i]] += M[which_upstream[i]];
  LL_M_obs = rep_vector(0,N);
  for(i in 1:N_M_obs)
    LL_M_obs[which_M_obs[i]] = lognormal_lpdf(M[which_M_obs[i]] | log(M_obs[which_M_obs[i]]), tau_M[which_M_obs[i]]); 
    // LL_M_obs[which_M_obs[i]] = lognormal_lpdf(M_obs[which_M_obs[i]] | log(M[which_M_obs[i]]), tau_M[which_M_obs[i]]); 
  LL_S_obs = rep_vector(0,N);
  for(i in 1:N_S_obs)
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S[which_S_obs[i]] | log(S_obs[which_S_obs[i]]), tau_S[which_S_obs[i]]); 
    // LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau_S[which_S_obs[i]]); 
  LL_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1,N_age);
  LL_n_O_obs = (n_O_obs .* log(q_O)) * rep_vector(1, 1 + N_H_pop);
  for(i in 1:N)
    LL_n_F_obs[i] = binomial_lpmf(n_F_obs[i] | n_MF_obs[i], q_F[i]);
  LL = LL_M_obs + LL_S_obs + LL_n_age_obs + LL_n_O_obs + LL_n_F_obs;
}
