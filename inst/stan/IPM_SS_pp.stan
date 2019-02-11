functions {
  # spawner-recruit functions
  real SR(int SR_fun, real alpha, real Rmax, real S, real A) {
    real R;
    
    if(SR_fun == 1)      # discrete exponential
      R = alpha*S/A;
    else if(SR_fun == 2) # Beverton-Holt
      R = alpha*S/(A + alpha*S/Rmax);
    else if(SR_fun == 3) # Ricker
      R = alpha*(S/A)*exp(-alpha*S/(A*e()*Rmax));
    
    return(R);
  }
  
  # Generalized normal (aka power-exponential) unnormalized log-probability
  real pexp_lpdf(real y, real mu, real sigma, real shape) {
    return(-(fabs(y - mu)/sigma)^shape);
  }
  
  # convert matrix to array of column vectors
  vector[] matrix_to_array(matrix m) {
    vector[2] arr[cols(m)];
    
    for(i in 1:cols(m))
      arr[i] = col(m,i);
    return(arr);
  }

  # Vectorized logical equality
  int[] veq(int[] x, int y) {
    int xeqy[size(x)];
    for(i in 1:size(x))
      xeqy[i] = x[i] == y;
    return(xeqy);
  }

  # Vectorized logical &&
  int[] vand(int[] cond1, int[] cond2) {
    int cond1_and_cond2[size(cond1)];
    for(i in 1:size(cond1))
      cond1_and_cond2[i] = cond1[i] && cond2[i];
    return(cond1_and_cond2);
  }

  # R-style conditional subsetting
  int[] rsub(int[] x, int[] cond) {
    int xsub[sum(cond)];
    int pos;
    pos = 1;
    for (i in 1:size(x))
      if (cond[i])
      {
        xsub[pos] = x[i];
        pos = pos + 1;
      }
    return(xsub);
  }

  # Equivalent of R: which(cond), where sum(cond) == 1
  int which(int[] cond) {
    int which_cond;
    for(i in 1:size(cond))
      if(cond[i])
        which_cond = i;
      return(which_cond);
  }
}

data {
  int<lower=1> SR_fun;                 # S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  int<lower=1> N;                      # total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];         # population identifier
  int<lower=1,upper=N> year[N];        # brood year identifier
  int<lower=0,upper=max(pop)> N_pop_H; # number of populations with hatchery input
  int<lower=1,upper=max(pop)> which_pop_H[max(N_pop_H,1)]; # populations with hatchery input
  int<lower=1,upper=N> N_S_obs;        # number of cases with non-missing spawner abundance obs 
  int<lower=1,upper=N> which_S_obs[N_S_obs]; # cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;        # observed annual total spawner abundance (not density)
  int<lower=2> N_age;                  # number of adult age classes
  int<lower=2> max_age;                # maximum adult age
  matrix<lower=0>[N,N_age] n_age_obs;  # observed wild spawner age frequencies (all zero row = NA)  
  int<lower=0,upper=N> N_H;            # number of years with p_HOS > 0
  int<lower=1,upper=N> which_H[max(N_H,1)]; # years with p_HOS > 0
  int<lower=0> n_W_obs[max(N_H,1)];    # count of wild spawners in samples (assumes no NAs)
  int<lower=0> n_H_obs[max(N_H,1)];    # count of hatchery spawners in samples (assumes no NAs)
  vector<lower=0>[N] A;                # habitat area associated with each spawner abundance obs
  vector<lower=0,upper=1>[N] F_rate;   # fishing mortality of wild adults
  int<lower=0,upper=N> N_B;            # number of years with B_take > 0
  int<lower=1,upper=N> which_B[max(N_B,1)]; # years with B_take > 0
  vector[max(N_B,1)] B_take_obs;       # observed broodstock take of wild adults
  int<lower=0> N_fwd;                  # total number of cases in forward simulations
  int<lower=1,upper=N> pop_fwd[max(N_fwd,1)]; # population identifier for forward simulations
  int<lower=1,upper=N+N_fwd> year_fwd[max(N_fwd,1)]; # brood year identifier for forward simulations
  vector<lower=0>[max(N_fwd,1)] A_fwd; # habitat area for each forward simulation
  vector<lower=0,upper=1>[max(N_fwd,1)] F_rate_fwd; # fishing mortality for forward simulations
  vector<lower=0,upper=1>[max(N_fwd,1)] B_rate_fwd; # broodstock take rate for forward simulations
  vector<lower=0,upper=1>[max(N_fwd,1)] p_HOS_fwd; # p_HOS for forward simulations
  int<lower=1> N_X;                    # number of productivity covariates
  matrix[max(max(year),max(year_fwd)),N_X] X; # brood-year productivity covariates (if none, use vector of zeros)
}

transformed data {
  int<lower=1,upper=N> N_pop;        # number of populations
  int<lower=1,upper=N> N_year;       # number of years, not including forward simulations
  int<lower=1,upper=N> N_year_all;   # total number of years, including forward simulations
  int<lower=2> ages[N_age];          # adult ages
  int<lower=0> n_HW_obs[max(N_H,1)]; # total sample sizes for H/W frequencies
  int<lower=1> pop_year_indx[N];     # index of years within each pop, starting at 1
  int<lower=0,upper=N> fwd_init_indx[max(N_fwd,1),N_age]; # links "fitted" brood years to recruits in forward sims
  
  N_pop = max(pop);
  N_year = max(year);
  N_year_all = max(max(year), max(year_fwd));
  for(a in 1:N_age)
    ages[a] = max_age - N_age + a;
  for(i in 1:max(N_H,1)) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];
  
  pop_year_indx[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year_indx[i] = 1;
    else
      pop_year_indx[i] = pop_year_indx[i-1] + 1;
  }
  
  fwd_init_indx = rep_array(0, max(N_fwd,1), N_age);
  if(N_fwd > 0)
  {
    for(i in 1:N_fwd)
    {
      for(a in 1:N_age)
      {
        if(year_fwd[i] - ages[a] < min(rsub(year_fwd, veq(pop_fwd, pop_fwd[i]))))
          fwd_init_indx[i,a] = which(vand(veq(pop, pop_fwd[i]), veq(year, year_fwd[i] - ages[a])));
      }
    }
  }
}

parameters {
  real mu_alpha;                         # hyper-mean log intrinsic productivity
  real<lower=0> sigma_alpha;             # hyper-SD log intrinsic productivity
  vector[N_pop] epsilon_alpha_z;         # log intrinsic prod (Z-scores)
  real mu_Rmax;                          # hyper-mean log asymptotic recruitment
  real<lower=0> sigma_Rmax;              # hyper-SD log asymptotic recruitment
  vector[N_pop] epsilon_Rmax_z;          # log asymptotic recruitment (Z-scores)
  real<lower=-1,upper=1> rho_alphaRmax;  # correlation between log(alpha) and log(Rmax)
  vector[N_X] beta_phi;                  # regression coefs for log productivity anomalies
  real<lower=-1,upper=1> rho_phi;        # AR(1) coef for log productivity anomalies
  real<lower=0> sigma_phi;               # hyper-SD of brood year log productivity anomalies
  vector[N_year_all] epsilon_phi_z;      # log brood year productivity anomalies (Z-scores)
  real<lower=0> sigma;                   # unique process error SD
  simplex[N_age] mu_p;                   # among-pop mean of age distributions
  vector<lower=0>[N_age-1] sigma_gamma;  # among-pop SD of mean log-ratio age distributions
  cholesky_factor_corr[N_age-1] L_gamma; # Cholesky factor of among-pop correlation matrix of mean log-ratio age distns
  matrix[N_pop,N_age-1] epsilon_gamma_z; # population mean log-ratio age distributions (Z-scores)
  vector<lower=0>[N_age-1] sigma_p;      # SD of log-ratio cohort age distributions
  cholesky_factor_corr[N_age-1] L_p;     # Cholesky factor of correlation matrix of cohort log-ratio age distributions
  matrix[N,N_age-1] epsilon_p_z;         # log-ratio cohort age distributions (Z-scores)
  vector<lower=0>[max_age*N_pop] S_init; # true total spawner abundance in years 1-max_age
  simplex[N_age] q_init[max_age*N_pop];  # true wild spawner age distributions in years 1-max_age
  vector<lower=0,upper=1>[max(N_H,1)] p_HOS; # true p_HOS in years which_H
  vector[N] epsilon_R_z;                 # log true recruit abundance (not density) by brood year (z-scores)
  vector<lower=0,upper=1>[max(N_B,1)] B_rate; # true broodstock take rate when B_take > 0
  real<lower=0> tau;                     # observation error SD of total spawners
}

transformed parameters {
  vector<lower=0>[N_pop] alpha;          # intrinsic productivity 
  vector<lower=0>[N_pop] Rmax;           # asymptotic recruitment 
  vector[N_year_all] phi;                # log brood year productivity anomalies
  vector<lower=0>[N] S_W;                # true total wild spawner abundance
  vector[N] S_H;                         # true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                  # true total spawner abundance
  row_vector[N_age-1] mu_gamma;          # mean of log-ratio cohort age distributions
  matrix[N_pop,N_age-1] gamma;           # population mean log-ratio age distributions
  matrix<lower=0,upper=1>[N,N_age] p;    # cohort age distributions
  matrix<lower=0,upper=1>[N,N_age] q;    # true spawner age distributions
  vector[N] p_HOS_all;                   # true p_HOS in all years (can == 0)
  vector<lower=0>[N] R_hat;              # expected recruit abundance (not density) by brood year
  vector<lower=0>[N] R;                  # true recruit abundance (not density) by brood year
  vector<lower=0,upper=1>[N] B_rate_all; # true broodstock take rate in all years
  
  # Multivariate Matt trick for [log(alpha), log(Rmax)]
  {
    matrix[2,2] L_alphaRmax;       # temp variable: Cholesky factor of corr matrix of log(alpha), log(Rmax)
    matrix[N_pop,2] epsilon_alphaRmax_z; # temp variable [log(alpha), log(Rmax)] random effects (z-scored)
    matrix[N_pop,2] epsilon_alphaRmax; # temp variable: [log(alpha), log(Rmax)] random effects
    vector[2] sigma_alphaRmax;     # temp variable: SD vector of [log(alpha), log(Rmax)]
    
    L_alphaRmax[1,1] = 1;
    L_alphaRmax[2,1] = rho_alphaRmax;
    L_alphaRmax[1,2] = 0;
    L_alphaRmax[2,2] = sqrt(1 - rho_alphaRmax^2);
    sigma_alphaRmax[1] = sigma_alpha;
    sigma_alphaRmax[2] = sigma_Rmax;
    epsilon_alphaRmax_z = append_col(epsilon_alpha_z, epsilon_Rmax_z);
    epsilon_alphaRmax = (diag_matrix(sigma_alphaRmax) * L_alphaRmax * epsilon_alphaRmax_z')';
    alpha = exp(mu_alpha + col(epsilon_alphaRmax,1));
    Rmax = exp(mu_Rmax + col(epsilon_alphaRmax,2));
  }
  
  # AR(1) model for phi
  phi[1] = epsilon_phi_z[1]*sigma_phi/sqrt(1 - rho_phi^2); # initial anomaly
  for(i in 2:N_year_all)
    phi[i] = rho_phi*phi[i-1] + epsilon_phi_z[i]*sigma_phi;
  # constrain "fitted" log anomalies to sum to 0 (X should be centered)
  phi = phi - mean(phi[1:N_year]) + X*beta_phi;
  
  # Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  if(N_H > 0)
    p_HOS_all[which_H] = p_HOS;
  
  B_rate_all = rep_vector(0,N);
  if(N_B > 0)
    B_rate_all[which_B] = B_rate;
  
  # Multivariate Matt trick for age vectors (pop-specific mean and within-pop, time-varying)
  mu_gamma = to_row_vector(log(mu_p[1:(N_age-1)]) - log(mu_p[N_age]));
  gamma = rep_matrix(mu_gamma,N_pop) + (diag_matrix(sigma_gamma) * L_gamma * epsilon_gamma_z')';
  p = append_col(gamma[pop,] + (diag_matrix(sigma_p) * L_p * epsilon_p_z')', rep_vector(0,N));
  
  # Calculate true total wild and hatchery spawners and spawner age distribution
  # and predict recruitment from brood year i
  for(i in 1:N)
  {
    row_vector[N_age] exp_p; # temp variable: exp(p[i,])
    row_vector[N_age] S_W_a; # temp variable: true wild spawners by age
    int ii;                  # temp variable: index into S_init and q_init
    
    # Inverse log-ratio transform of cohort age distn
    # (built-in softmax function doesn't accept row vectors)
    exp_p = exp(p[i,]);
    p[i,] = exp_p/sum(exp_p);
    ii = (pop[i] - 1)*max_age + pop_year_indx[i];
    
    if(pop_year_indx[i] <= max_age)
    {
      # Use initial values
      S_W[i] = S_init[ii]*(1 - p_HOS_all[i]);        
      S_H[i] = S_init[ii]*p_HOS_all[i];
      q[i,1:N_age] = to_row_vector(q_init[ii, 1:N_age]);
      S_W_a = S_W[i]*q[i,];
    }
    else
    {
      # Use recruitment process model
      for(a in 1:N_age)
        S_W_a[a] = R[i-ages[a]]*p[i-ages[a],a];
      for(a in 2:N_age)  # catch and broodstock removal (assumes no take of age 1)
        S_W_a[a] = S_W_a[a]*(1 - F_rate[i])*(1 - B_rate_all[i]);
      S_W[i] = sum(S_W_a);
      S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
      q[i,] = S_W_a/S_W[i];
    }
    
    S[i] = S_W[i] + S_H[i];
    R_hat[i] = A[i] * SR(SR_fun, alpha[pop[i]], Rmax[pop[i]], S[i], A[i]);
    R[i] = R_hat[i]*exp(phi[year[i]] + sigma*epsilon_R_z[i]);
  }
}

model {
  vector[max(N_B,1)] B_take; # true broodstock take when B_take_obs > 0
  
  # Priors
  mu_alpha ~ normal(2,5);
  sigma_alpha ~ pexp(0,3,10);
  mu_Rmax ~ normal(0,10);
  sigma_Rmax ~ pexp(0,3,10);
  beta_phi ~ normal(0,5);
  rho_phi ~ pexp(0,0.85,50);  # mildly regularize to ensure stationarity
  sigma_phi ~ pexp(0,2,10);
  sigma ~ pexp(0,1,10);
  for(i in 1:(N_age-1))
  {
    sigma_gamma[i] ~ pexp(0,2,5);
    sigma_p[i] ~ pexp(0,2,5); 
  }
  L_gamma ~ lkj_corr_cholesky(1);
  L_p ~ lkj_corr_cholesky(1);
  tau ~ pexp(0,1,10);
  S_init ~ lognormal(0,10);
  if(N_B > 0)
  {
    B_take = B_rate .* S_W[which_B] .* (1 - q[which_B,1]) ./ (1 - B_rate);
    B_take_obs ~ lognormal(log(B_take), 0.1); # penalty to force pred and obs broodstock take to match 
  }
  
  # Hierarchical priors
  # [log(alpha), log(Rmax)] ~ MVN([mu_alpha, mu_Rmax], D*R_aRmax*D), 
  # where D = diag_matrix(sigma_alpha, sigma_Rmax)
  epsilon_alpha_z ~ normal(0,1);
  epsilon_Rmax_z ~ normal(0,1);
  epsilon_phi_z ~ normal(0,1);   # phi[i] ~ N(rho_phi*phi[i-1], sigma_phi)
  # pop mean age probs logistic MVN: 
  # gamma[i,] ~ MVN(mu_gamma,D*R_gamma*D), 
  # where D = diag_matrix(sigma_gamma)
  to_vector(epsilon_gamma_z) ~ normal(0,1);
  
  # Process model
  # age probs logistic MVN: 
  # alr_p[i,] ~ MVN(gamma[pop[i],], D*R_p*D), 
  # where D = diag_matrix(sigma_p)
  to_vector(epsilon_p_z) ~ normal(0,1);
  epsilon_R_z ~ normal(0,1); # total recruits: R ~ lognormal(log(R_hat), sigma)
  
  # Observation model
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau);  # observed total spawners
  if(N_H > 0) n_H_obs ~ binomial(n_HW_obs, p_HOS); # observed counts of hatchery vs. wild spawners
  target += sum(n_age_obs .* log(q)); # obs wild age freq: n_age_obs[i] ~ multinomial(q[i])
}

generated quantities {
  corr_matrix[N_age-1] R_gamma;     # among-pop correlation matrix of mean log-ratio age distns 
  corr_matrix[N_age-1] R_p;         # correlation matrix of within-pop cohort log-ratio age distns 
  vector<lower=0>[N_fwd] S_W_fwd;   # true total wild spawner abundance in forward simulations
  vector[N_fwd] S_H_fwd;            # true total hatchery spawner abundance in forward simulations
  vector<lower=0>[N_fwd] S_fwd;     # true total spawner abundance in forward simulations
  matrix<lower=0,upper=1>[N_fwd,N_age] p_fwd; # cohort age distributions in forward simulations
  matrix<lower=0,upper=1>[N_fwd,N_age] q_fwd; # spawner age distributions in forward simulations
  vector<lower=0>[N_fwd] R_hat_fwd; # expected recruit abundance by brood year in forward simulations
  vector<lower=0>[N_fwd] R_fwd;     # true recruit abundance by brood year in forward simulations
  vector[N] LL_S_obs;               # pointwise log-likelihood of total spawners
  vector[max(N_H,1)] LL_n_H_obs;    # pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL_n_age_obs;           # pointwise log-likelihood of wild age frequencies
  vector[N] LL;                     # total pointwise log-likelihood                              
  
  R_gamma = multiply_lower_tri_self_transpose(L_gamma);
  R_p = multiply_lower_tri_self_transpose(L_p);
  
  # Calculate true total wild and hatchery spawners and spawner age distribution
  # and simulate recruitment from brood year i
  # (Note that if N_fwd == 0, this block will not execute)
  for(i in 1:N_fwd)
  {
    vector[N_age-1] alr_p_fwd;   # temp variable: alr(p_fwd[i,])'
    row_vector[N_age] S_W_a_fwd;   # temp variable: true wild spawners by age

    # Inverse log-ratio transform of cohort age distn
    alr_p_fwd = multi_normal_cholesky_rng(to_vector(gamma[pop_fwd[i],]), L_p);
    p_fwd[i,] = to_row_vector(softmax(append_row(alr_p_fwd,0)));

    for(a in 1:N_age)
    {
      if(fwd_init_indx[i,a] != 0)
      {
        # Use estimated values from previous cohorts
        S_W_a_fwd[a] = R[fwd_init_indx[i,a]]*p[fwd_init_indx[i,a],a];
      }
      else
      {
        S_W_a_fwd[a] = R_fwd[i-ages[a]]*p_fwd[i-ages[a],a];
      }
    }

    for(a in 2:N_age)  # catch and broodstock removal (assumes no take of age 1)
      S_W_a_fwd[a] = S_W_a_fwd[a]*(1 - F_rate_fwd[i])*(1 - B_rate_fwd[i]);

    S_W_fwd[i] = sum(S_W_a_fwd);
    S_H_fwd[i] = S_W_fwd[i]*p_HOS_fwd[i]/(1 - p_HOS_fwd[i]);
    q_fwd[i,] = S_W_a_fwd/S_W_fwd[i];
    S_fwd[i] = S_W_fwd[i] + S_H_fwd[i];
    R_hat_fwd[i] = A_fwd[i] * SR(SR_fun, alpha[pop_fwd[i]], Rmax[pop_fwd[i]], S_fwd[i], A_fwd[i]);
    R_fwd[i] = lognormal_rng(log(R_hat_fwd[i]) + phi[year_fwd[i]], sigma);
  }
  
  LL_S_obs = rep_vector(0,N);
  for(i in 1:N_S_obs)
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau); 
  LL_n_age_obs = (n_age_obs .* log(q)) * rep_vector(1,N_age);
  LL_n_H_obs = rep_vector(0,max(N_H,1));
  if(N_H > 0)
  {
    for(i in 1:N_H)
      LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
  }
  LL = LL_S_obs + LL_n_age_obs;
  LL[which_H] = LL[which_H] + LL_n_H_obs;
}
