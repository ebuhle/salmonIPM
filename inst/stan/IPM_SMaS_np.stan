functions {
  #include /include/SR.stan
  #include /include/gnormal_lpdf_vec.stan
  #include /include/col_sums.stan
  #include /include/rep_vec.stan
  #include /include/to_row_vector_row_major.stan
  #include /include/seq.stan
}

data {
  // info for observed data
  int<lower=1> N;                       // total number of cases in all pops and years
  array[N] int<lower=1,upper=N> pop;    // population index
  array[N] int<lower=1,upper=N> year;   // calendar year index
  // smolt production
  int<lower=1> SR_fun;                  // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  array[2] int<lower=0,upper=1> RRS;    // fit H vs. W {alpha, Rmax} (1) or not (0)?
  vector[N] A;                          // habitat area associated with each spawner abundance obs
  int<lower=0> K_alpha;                 // number of intrinsic productivity covariates
  matrix[N,K_alpha] X_alpha;            // intrinsic productivity covariates
  array[2] real prior_alpha;            // prior meanlog, sdlog for intrinsic productivity
  array[2] real prior_alpha_W;          // prior meanlog, sdlog for W intrinsic productivity
  array[2] real prior_alpha_H;          // prior meanlog, sdlog for H intrinsic productivity
  int<lower=0> K_Mmax;                  // number of maximum smolt recruitment covariates
  matrix[N,K_Mmax] X_Mmax;              // maximum smolt recruitment covariates
  array[2] real prior_Mmax;             // prior meanlog, sdlog for maximum smolt recruitment
  array[2] real prior_Mmax_W;           // prior meanlog, sdlog for W maximum smolt recruitment
  array[2] real prior_Mmax_H;           // prior meanlog, sdlog for H maximum smolt recruitment
  int<lower=0> K_M;                     // number of smolt recruitment covariates
  array[N] row_vector[K_M] X_M;         // smolt recruitment covariates
  // smolt abundance
  int<lower=1,upper=N> N_M_obs;         // number of cases with non-missing smolt abundance obs 
  array[N_M_obs] int<lower=1,upper=N> which_M_obs; // cases with non-missing smolt abundance obs
  vector<lower=0>[N] M_obs;             // observed annual smolt abundance (not density)
  array[3] real prior_tau_M;            // prior mean, scale, shape for smolt observation error SD
  // smolt age structure
  int<lower=2> N_Mage;                  // number of smolt age classes
  int<lower=2> max_Mage;                // maximum smolt age
  matrix<lower=0>[N,N_Mage] n_Mage_obs; // observed smolt age frequencies (all zero row = NA)  
  // SAR (sMolt-Spawner survival)
  int<lower=0> K_MS;                    // number of SAR covariates
  array[N] row_vector[K_MS] X_MS;       // SAR covariates
  // spawner abundance
  int<lower=1,upper=N> N_S_obs;         // number of cases with non-missing spawner abundance obs 
  array[N_S_obs] int<lower=1,upper=N> which_S_obs; // cases with non-missing spawner abundance obs
  vector<lower=0>[N] S_obs;             // observed annual total spawner abundance (not density)
  array[3] real prior_tau_S;            // prior mean, scale, shape for spawner observation error SD
  // spawner ocean age and Gilbert-Rich age structure
  int<lower=2> N_MSage;                 // number of ocean age classes
  int<lower=1> max_MSage;               // maximum ocean age
  matrix<lower=0>[N,N_MSage] n_MSage_obs; // observed ocean age frequencies (all zero row = NA)  
  matrix<lower=0>[N,N_Mage*N_MSage] n_GRage_obs; // obs W spawner Gilbert-Rich age freqs (all zero row = NA)  
  int<lower=0,upper=1> conditionGRonMS; // is n_GRage_obs unconditional (0) or conditioned on ocean age totals (1)?
  // H/W composition
  int<lower=0,upper=N> N_H;             // number of years with p_HOS > 0
  array[N_H] int<lower=1,upper=N> which_H; // years with p_HOS > 0
  array[N_H] int<lower=0> n_W_obs;     // count of wild spawners in samples
  array[N_H] int<lower=0> n_H_obs;     // count of hatchery spawners in samples
  // fishery and hatchery removals
  vector[N] F_rate;                    // fishing mortality rate of wild adults
  vector<lower=0,upper=1>[N_MSage] MSage_F; // is ocean age a (non)selected (0/1) by fishery?
  int<lower=0,upper=N> N_B;            // number of years with B_take > 0
  array[N_B] int<lower=1,upper=N> which_B; // years with B_take > 0
  vector[N_B] B_take_obs;              // observed broodstock take of wild adults
  vector<lower=0,upper=1>[N_MSage] MSage_B; // is ocean age a (non)selected (0/1) in broodstock?
}

transformed data {
  int<lower=1,upper=N> N_pop = max(pop);   // number of populations
  int<lower=1,upper=N> N_year = max(year); // number of years
  int<lower=1> min_Mage = max_Mage - N_Mage + 1; // minimum smolt age
  array[N_Mage] int<lower=0> Mages = seq(min_Mage, max_Mage); // smolt ages
  int<lower=0> min_MSage = max_MSage - N_MSage + 1; // minimum ocean age
  array[N_MSage] int<lower=0> MSages = seq(min_MSage, max_MSage); // ocean ages
  int<lower=2> max_age = max_Mage + max_MSage; // maximum adult age
  int<lower=2> N_GRage = N_Mage*N_MSage;   // number of Gilbert-Rich age classes
  array[N] int<lower=1> pop_year;          // index of years within each pop, starting at 1
  array[N_H] int<lower=0> n_HW_obs;        // total sample sizes for H/W frequencies
  vector[max_Mage*N_pop] mu_M_init;        // prior mean of total smolt abundance in years 1:max_Mage
  real sigma_M_init = sd(log(M_obs[which_M_obs])); // prior log-SD of smolt abundance in years 1:max_Mage
  matrix[N_Mage,max_Mage*N_pop] mu_q_M_init; // prior counts of smolt age distns in years 1:max_Mage
  vector[max_MSage*N_pop] mu_S_init;       // prior mean of total spawner abundance in years 1:max_MSage
  real sigma_S_init = sd(log(S_obs[which_S_obs])); // prior log-SD of spawner abundance in years 1:max_MSage
  matrix[N_GRage,max_MSage*N_pop] mu_q_GR_init; // prior counts of wild spawner age distns in years 1:max_MSage
  
  // Years within each pop starting at 1 (determines when to use initial values)
  pop_year[1] = 1;
  for(i in 1:N)
  {
    if(i == 1 || pop[i-1] != pop[i])
      pop_year[i] = 1;
    else
      pop_year[i] = pop_year[i-1] + 1;
  }
  
  // Array-based calculations 
  for(i in 1:N_H) n_HW_obs[i] = n_H_obs[i] + n_W_obs[i];

  // Priors on M_init, q_M_init accounting for amalgamation of orphan smolt age classes
  for(i in 1:max_Mage)
  {
    int N_orphan_Mage = N_Mage - max(i - min_Mage, 0); // number of orphan smolt age classes
    int N_amalg_Mage = N_Mage - N_orphan_Mage + 1; // number of amalgamated smolt age classes
    
    for(j in 1:N_pop)
    {
      int ii = (j - 1)*max_Mage + i; // index into M_init, q_M_init

      // M_init prior mean scales observed log-mean by fraction of orphan age classes
      mu_M_init[ii] = mean(log(M_obs[which_M_obs])) + log(N_orphan_Mage) - log(N_Mage);

      // prior on q_M_init implies q_M_orphan ~ Dir(1)
      mu_q_M_init[,ii] = append_row(rep_vector(1.0/N_amalg_Mage, N_amalg_Mage), 
                                    rep_vector(1, N_orphan_Mage - 1));
    }
  }
  
  // Priors on S_init, q_GR_init accounting for amalgamation of orphan ocean-age classes
  for(i in 1:max_MSage)
  {
    int N_orphan_MSage = N_MSage - max(i - min_MSage, 0); // number of orphan ocean-age classes
    int N_amalg_MSage = N_MSage - N_orphan_MSage + 1; // number of amalgamated ocean-age classes
    
    for(j in 1:N_pop)
    {
      int ii = (j - 1)*max_MSage + i; // index into S_init, q_GR_init
      
      // S_init prior mean scales observed log-mean by fraction of orphan age classes
      mu_S_init[ii] = mean(log(S_obs[which_S_obs])) + log(N_orphan_MSage) - log(N_MSage);

      // prior on q_GR_init implies q_GR_orphan ~ Dir(1)
      mu_q_GR_init[,ii] = rep_vec(append_row(rep_vector(1.0/N_amalg_MSage, N_amalg_MSage), 
                                                        rep_vector(1, N_orphan_MSage - 1)), N_Mage);
    }
  }
}

parameters {
  //?// indicates params that could be arrays instead of matrices
  // smolt recruitment
  vector<lower=0>[!RRS[1]*N_pop] alpha;       // intrinsic spawner-smolt productivity
  vector<lower=0>[RRS[1]*N_pop] alpha_W;      // intrinsic spawner-smolt productivity of wild spawners
  vector<lower=0>[RRS[1]*N_pop] alpha_H;      // intrinsic spawner-smolt productivity of hatchery spawners
  matrix[N_pop,K_alpha] beta_alpha;           // regression coefs for log alpha
  vector<lower=0>[!RRS[2]*N_pop] Mmax;        // maximum smolt recruitment
  vector<lower=0>[RRS[2]*N_pop] Mmax_W;       // maximum smolt recruitment of wild spawners
  vector<lower=0>[RRS[2]*N_pop] Mmax_H;       // maximum smolt recruitment of hatchery spawners
  matrix[N_pop,K_Mmax] beta_Mmax;             // regression coefs for log Mmax
  matrix[N_pop,K_M] beta_M;                   //?// regression coefs for spawner-smolt productivity
  vector<lower=-1,upper=1>[N_pop] rho_M;      // AR(1) coefs for spawner-smolt productivity
  vector<lower=0>[N_pop] sigma_M;             // spawner-smolt process error SDs
  vector[N] zeta_M;                           // smolt recruitment process errors (z-scored)
  // smolt age structure
  array[N_pop] simplex[N_Mage] mu_p_M;        // population mean smolt age distributions
  matrix<lower=0>[N_pop,N_Mage-1] sigma_p_M;  //?// log-ratio cohort smolt age distribution SDs
  array[N_pop] cholesky_factor_corr[N_Mage-1] L_p_M; // Cholesky factors of corr matrices of log-ratio smolt age distns
  matrix[N,N_Mage-1] zeta_p_M;                //?// log-ratio cohort smolt age distn errors (Z-scored)
  // SAR
  matrix<lower=0,upper=1>[N_pop,N_Mage] mu_MS; //?// mean SAR for each smolt age
  matrix[N_pop,K_MS] beta_MS;                 //?// regression coefs for SAR (independent of smolt age)  
  matrix<lower=-1,upper=1>[N_pop,N_Mage] rho_MS; //?// AR(1) coefs of SAR for each smolt age
  matrix<lower=0>[N_pop,N_Mage] sigma_MS;     //?// SAR process error SDs for each smolt age
  array[N_pop] cholesky_factor_corr[N_Mage] L_MS; // Cholesky factors of corr matrices of SAR across smolt ages
  matrix[N,N_Mage] zeta_MS;                   //?// SAR process errors for each smolt age (z-scored)
  // ocean age structure
  array[N_pop,N_Mage] simplex[N_MSage] mu_p_MS; // pop mean ocean age distributions for each smolt age
  array[N_pop,N_Mage] vector<lower=0>[N_MSage-1] sigma_p_MS; // log-ratio ocean age SDs for each smolt age
  array[N_pop] cholesky_factor_corr[N_Mage*(N_MSage-1)] L_p_MS; // Cholesky factors of corr matrices of log-ratio ocean age
  matrix[N,N_Mage*(N_MSage-1)] zeta_p_MS;     //?// log-ratio ocean age errors (Z-scored)
  // H/W composition, removals
  vector<lower=0,upper=1>[N_H] p_HOS;         // true p_HOS in years which_H
  vector<lower=0,upper=1>[N_B] B_rate;        // true broodstock take rate when B_take > 0
  // initial states, observation error
  vector<lower=0>[max_Mage*N_pop] M_init;     // true smolt abundance in years 1:max_Mage
  array[max_Mage*N_pop] simplex[N_Mage] q_M_init; // true smolt age distns in years 1:max_Mage
  vector<lower=0>[max_MSage*N_pop] S_init;    // true total spawner abundance in years 1:max_MSage
  array[max_MSage*N_pop] simplex[N_GRage] q_GR_init; // true wild spawner age distns in years 1:max_MSage
  vector<lower=0>[N_pop] tau_M;               // smolt observation error SDs
  vector<lower=0>[N_pop] tau_S;               // spawner observation error SDs
}

transformed parameters {
  //?// indicates transformed params that could be arrays instead of matrices
  // smolt recruitment
  vector[N] alpha_Xbeta;                 // intrinsic productivity including covariate effects
  vector[N] alpha_W_Xbeta;               // wild intrinsic productivity including covariate effects
  vector[N] alpha_H_Xbeta;               // hatchery intrinsic productivity including covariate effects
  vector[N] Mmax_Xbeta;                  // maximum recruitment including covariate effects
  vector[N] Mmax_W_Xbeta;                // wild maximum recruitment including covariate effects
  vector[N] Mmax_H_Xbeta;                // hatchery maximum recruitment including covariate effects
  vector<lower=0>[N] M_hat;              // expected smolt abundance (not density) by brood year
  vector[N] epsilon_M;                   // process error in smolt abundance by brood year
  vector<lower=0>[N] M0;                 // true smolt abundance (not density) by brood year
  vector<lower=0>[N] M;                  // true smolt abundance (not density) by outmigration year
  // smolt age structure
  array[N_pop] vector[N_Mage-1] mu_alr_p_M; // population mean log ratio smolt age distributions
  array[N] simplex[N_Mage] p_M;          // true smolt age distributions by brood year
  matrix<lower=0,upper=1>[N,N_Mage] q_M; // true smolt age distributions by calendar year
  // SAR
  matrix[N,N_Mage] epsilon_MS;           //?// SAR process errors by smolt age and outmigration year
  matrix<lower=0,upper=1>[N,N_Mage] s_MS; //?// true SAR by smolt age and outmigration year
  array[N_pop,N_Mage] vector[N_MSage-1] mu_alr_p_MS; // population mean log ratio age distributions
  array[N,N_Mage] simplex[N_MSage] p_MS; // true ocean age distns by outmigration year
  // H/W spawner abundance, removals
  vector[N] p_HOS_all;                   // true p_HOS in all years (can == 0)
  vector<lower=0>[N] S_W;                // true total wild spawner abundance
  vector[N] S_H;                         // true total hatchery spawner abundance (can == 0)
  vector<lower=0>[N] S;                  // true total spawner abundance
  vector<lower=0,upper=1>[N] B_rate_all; // true broodstock take rate in all years
  // spawner age structure
  matrix<lower=0,upper=1>[N,N_GRage] q_GR; // true Gilbert-Rich age distns of spawners
  matrix<lower=0,upper=1>[N,N_MSage] q_MS; // true ocean age distns of spawners
  
  // Pad p_HOS and B_rate
  p_HOS_all = rep_vector(0,N);
  p_HOS_all[which_H] = p_HOS;
  B_rate_all = rep_vector(0,N);
  B_rate_all[which_B] = B_rate;
  
  // S-R parameters including covariate effects
  if(RRS[1]) {
    alpha_W_Xbeta = alpha_W[pop] .* exp(rows_dot_product(X_alpha, beta_alpha[pop,]));
    alpha_H_Xbeta = alpha_H[pop] .* exp(rows_dot_product(X_alpha, beta_alpha[pop,]));
  } else {
    alpha_Xbeta = alpha[pop] .* exp(rows_dot_product(X_alpha, beta_alpha[pop,]));
  }
  if(RRS[2]) {
    Mmax_W_Xbeta = Mmax_W[pop] .* exp(rows_dot_product(X_Mmax, beta_Mmax[pop,]));
    Mmax_H_Xbeta = Mmax_H[pop] .* exp(rows_dot_product(X_Mmax, beta_Mmax[pop,]));
  } else {
    Mmax_Xbeta = Mmax[pop] .* exp(rows_dot_product(X_Mmax, beta_Mmax[pop,]));
  }

  // Log-ratio transform of pop-specific mean cohort age distributions
  for(j in 1:N_pop)
  {
    // Smolt age 
    mu_alr_p_M[j] = log(head(mu_p_M[j], N_Mage-1)) - log(mu_p_M[j,N_Mage]);
    
    // Ocean age
    for(a in 1:N_Mage)
      mu_alr_p_MS[j,a] = log(head(mu_p_MS[j,a,], N_MSage-1)) - log(mu_p_MS[j,a,N_MSage]);
  }
  
  // Calculate true smolts and total wild and hatchery spawners by age,
  // and predict smolt recruitment from brood year i
  for(i in 1:N)
  {
    int mm;                                   // index into M_init, q_M_init
    // number of orphan smolt age classes <lower=0,upper=N_Mage>
    int N_orphan_Mage = max(N_Mage - to_int(fmin(pop_year[i], min_Mage)), N_Mage); 
    // orphan smolt age distn (amalgamated simplex)
    vector[N_orphan_Mage] q_M_orphan;
    vector[N_Mage] alr_p_M;                   // alr(p_M[i,])
    row_vector[N_Mage] M_a;                   // true smolts by age
    vector[N_Mage*(N_MSage-1)] mu_alr_p_MS_i; // mu_alr_p_MS[i,,] 
    vector[N_Mage*(N_MSage-1)] sigma_p_MS_i;  // sigma_p_MS[i,,]
    vector[N_Mage*(N_MSage-1)] alr_p_MS;      // alr(p_MS[i,]) 
    int ss;                                   // index into S_init, q_GR_init
    // number of orphan ocean-age classes <lower=0,upper=N_MSage>
    int N_orphan_MSage = max(N_MSage - to_int(fmin(pop_year[i], min_MSage)), N_MSage); 
    // slice of orphan G-R age distn for a given smolt age (amalgamated simplex)
    vector[N_orphan_MSage] q_GR_orphan;
    matrix[N_Mage,N_MSage] S_W_a;             // true W spawners by smolt and ocean age
    
    // Time-varying IID age vectors (multivariate Matt trick)
    // Smolt age
    alr_p_M = rep_vector(0,N_Mage);
    alr_p_M[1:(N_Mage-1)] = mu_alr_p_M[pop[i]] + sigma_p_M[pop[i],]' .* (L_p_M[pop[i]] * zeta_p_M[i,]');
    alr_p_M = exp(alr_p_M);
    p_M[i] = alr_p_M / sum(alr_p_M);
    
    // Ocean age
    // assemble mean and SD vectors by flattening across smolt age
    for(a in 1:N_Mage)
    {
      mu_alr_p_MS_i[((a-1)*(N_MSage-1) + 1):(a*(N_MSage-1))] = mu_alr_p_MS[pop[i],a];
      sigma_p_MS_i[((a-1)*(N_MSage-1) + 1):(a*(N_MSage-1))] = sigma_p_MS[pop[i],a];
    }
    alr_p_MS = mu_alr_p_MS_i + sigma_p_MS_i .* (L_p_MS[pop[i]] * zeta_p_MS[i,]');
    // inverse log-ratio transform and assign back to array
    for(a in 1:N_Mage)
    {
      p_MS[i,a] = exp(append_row(alr_p_MS[((a-1)*(N_MSage-1) + 1):(a*(N_MSage-1))], 0));
      p_MS[i,a] = p_MS[i,a] / sum(p_MS[i,a]);
    }
    
    // AR(1) smolt recruitment process errors and MAR(1) SAR process errors  
    if(pop_year[i] == 1) // initial process error
    {
      epsilon_M[i] = zeta_M[i] * sigma_M[pop[i]] / sqrt(1 - rho_M[pop[i]]^2);
      // cheat: doesn't use MAR(1) stationary covariance
      epsilon_MS[i,] = zeta_MS[i,] .* sigma_MS[pop[i],] ./ sqrt(1 - square(rho_MS[pop[i],]));
    }
    else
    {
      epsilon_M[i] = rho_M[pop[i]] * epsilon_M[i-1] + zeta_M[i] * sigma_M[pop[i]];
      epsilon_MS[i,] = rho_MS[pop[i],] .* epsilon_MS[i-1,] + sigma_MS[pop[i],] .* (L_MS[pop[i]] * zeta_MS[i,]')';
    }
    // SAR for outmigration year i
    s_MS[i,] = inv_logit(logit(mu_MS[pop[i],]) + dot_product(X_MS[i], beta_MS[pop[i],]) + epsilon_MS[i,]);
    
    // Smolt recruitment
    // Use initial values for orphan age classes, otherwise use process model
    if(pop_year[i] <= max_Mage)
    {
      mm = (pop[i] - 1)*max_Mage + pop_year[i];
      q_M_orphan = append_row(sum(head(q_M_init[mm], N_Mage - N_orphan_Mage + 1)), 
                              tail(q_M_init[mm], N_orphan_Mage - 1));
    }
    
    for(a in 1:N_Mage)
    {
      if(pop_year[i] <= Mages[a]) // use initial values
        M_a[a] = M_init[mm]*q_M_orphan[a - (N_Mage - N_orphan_Mage)];
      else // age-a smolts from appropriate brood year
        M_a[a] = M0[i-Mages[a]]*p_M[i-Mages[a],a]; 
    }
    
    M[i] = sum(M_a);
    q_M[i,] = M_a / M[i];
    
    // Spawners and age structure
    // Use initial values for orphan age classes, otherwise use process model
    for(sa in 1:N_Mage)
    {
      if(pop_year[i] <= max_MSage)
      {
        int aa[N_MSage] = seq((sa - 1)*N_MSage + 1, sa*N_MSage); // slice of q_GR_init for smolt age sa
        ss = (pop[i] - 1)*max_MSage + pop_year[i];
        q_GR_orphan = append_row(sum(head(q_GR_init[ss][aa], N_MSage - N_orphan_MSage + 1)), 
                                 tail(q_GR_init[ss][aa], N_orphan_MSage - 1));
      }
      
      for(oa in 1:N_MSage)
      {
        if(pop_year[i] <= MSages[oa]) // use initial values
          S_W_a[sa,oa] = S_init[ss]*(1 - p_HOS_all[i])*q_GR_orphan[oa - (N_MSage - N_orphan_MSage)];
        else // use recruitment process model
          S_W_a[sa,oa] = M[i-MSages[oa]]*q_M[i-MSages[oa],sa]*s_MS[i-MSages[oa],sa]*p_MS[i-MSages[oa],sa][oa] * 
                         (1 - MSage_F[oa]*F_rate[i])*(1 - MSage_B[oa]*B_rate_all[i]);
      }
    }
    
    // Total spawners and age structure
    S_W[i] = sum(S_W_a);
    S_H[i] = S_W[i]*p_HOS_all[i]/(1 - p_HOS_all[i]);
    S[i] = S_W[i] + S_H[i];
    q_MS[i,] = col_sums(S_W_a)/S_W[i];
    if(conditionGRonMS)  // q_GR := probabilities of smolt age conditioned on ocean age
      q_GR[i,] = to_row_vector_row_major(diag_post_multiply(S_W_a/S_W[i], 1 ./ q_MS[i,]));
    else  // q_GR := unconditional probabilities of Gilbert-Rich age 
      q_GR[i,] = to_row_vector_row_major(S_W_a/S_W[i]);

    // Smolt production from brood year i
    M_hat[i] = SR(SR_fun, RRS, alpha_Xbeta[i], alpha_W_Xbeta[i], alpha_H_Xbeta[i],
                  Mmax_Xbeta[i], Mmax_W_Xbeta[i], Mmax_H_Xbeta[i], S[i], S_W[i], S_H[i], A[i]);
    M0[i] = M_hat[i] * exp(dot_product(X_M[i], beta_M[pop[i],]) + epsilon_M[i]);
  }
}

model {
  vector[N_B] log_B_take; // log of true broodstock take when B_take_obs > 0
  
  // Priors
  
  // smolt recruitment
  if(RRS[1]) {
    alpha_W ~ lognormal(prior_alpha_W[1], prior_alpha_W[2]);
    alpha_H ~ lognormal(prior_alpha_H[1], prior_alpha_H[2]);
  } else {
    alpha ~ lognormal(prior_alpha[1], prior_alpha[2]);
  }
  to_vector(beta_alpha) ~ normal(0,5);
  if(RRS[2]) {
    Mmax_W ~ lognormal(prior_Mmax_W[1], prior_Mmax_W[2]);
    Mmax_H ~ lognormal(prior_Mmax_H[1], prior_Mmax_H[2]);
  } else {
    Mmax ~ lognormal(prior_Mmax[1], prior_Mmax[2]);
  }
  to_vector(beta_Mmax) ~ normal(0,5);
  to_vector(beta_M) ~ normal(0,5);
  rho_M ~ gnormal(0,0.85,20); // mildly regularize rho to ensure stationarity
  sigma_M ~ normal(0,3);
  zeta_M ~ std_normal();       // total smolts: log(M) ~ normal(log(M_hat), sigma_M)

  // smolt age structure
  to_vector(sigma_p_M) ~ normal(0,3);
  for(j in 1:N_pop)
    L_p_M[j] ~ lkj_corr_cholesky(1); // R_p_M[j] ~ lkj_corr(1)
  // smolt age probs logistic MVN: 
  // alr(p_M[i,]) ~ MVN(mu_alr_p_M[pop[i],], D*R_p_M*D), where D = diag_matrix(sigma_p_M[pop[i],])
  to_vector(zeta_p_M) ~ std_normal();

  // SAR
  to_vector(beta_MS) ~ normal(0,3);
  to_vector(sigma_MS) ~ normal(0,3);
  to_vector(rho_MS) ~ gnormal(0,0.85,20); // mildly regularize rho to ensure stationarity
  for(j in 1:N_pop)
    L_MS[j] ~ lkj_corr_cholesky(1);  // R_MS[j] ~ lkj_corr(1)
  to_vector(zeta_MS) ~ std_normal(); // SAR: logit(s_MS) ~ normal(logit(s_MS_hat), sigma_MS)

  // ocean age structure
  for(j in 1:N_pop)
  {
    for(a in 1:N_Mage)
      sigma_p_MS[j,a] ~ normal(0,3);
    L_p_MS[j] ~ lkj_corr_cholesky(1); // R_p_MS[j] ~ lkj_corr(1)
  }
  // ocean age probs logistic MVN: 
  // alr(p_MS[i,]) ~ MVN(mu_alr_p_MS[pop[i],,], D*R_p_MS*D), where D = diag_matrix(sigma_p_MS[pop[i],,])
  to_vector(zeta_p_MS) ~ std_normal();

  // removals
  log_B_take = log(S_W[which_B]) + log(q_MS[which_B,]*MSage_B) + logit(B_rate); 
  // implies B_take[i] = S_W[i] * (q_MS[i,] * MSage_B) * B_rate[i] / (1 - B_rate[i])
  B_take_obs ~ lognormal(log_B_take, 0.05); // penalty to force pred and obs broodstock take to match 

  // initial states (accounting for amalgamation of q_init to q_orphan)
  // smolt abundance and age structure
  M_init ~ lognormal(mu_M_init, sigma_M_init);
  {
    matrix[N_Mage,max_Mage*N_pop] q_M_init_mat;
    
    for(j in 1:size(q_M_init)) q_M_init_mat[,j] = q_M_init[j];
    target += sum((mu_q_M_init - 1) .* log(q_M_init_mat)); // q_M_init[i] ~ Dir(mu_q_M_init[,i])
  }
  // spawner abundance and age structure
  S_init ~ lognormal(mu_S_init, sigma_S_init);
  {
    matrix[N_GRage,max_MSage*N_pop] q_GR_init_mat;
    
    for(j in 1:size(q_GR_init)) q_GR_init_mat[,j] = q_GR_init[j];
    target += sum((mu_q_GR_init - 1) .* log(q_GR_init_mat)); // q_GR_init[i] ~ Dir(mu_q_GR_init[,i])
  }

  // observation error
  tau_M ~ gnormal(prior_tau_M[1], prior_tau_M[2], prior_tau_M[3]);
  tau_S ~ gnormal(prior_tau_S[1], prior_tau_S[2], prior_tau_S[3]);

  // Observation model
  M_obs[which_M_obs] ~ lognormal(log(M[which_M_obs]), tau_M[pop[which_M_obs]]);  // observed smolts
  target += sum(n_Mage_obs .* log(q_M));   // obs smolt age freq: n_Mage_obs[i] ~ multinomial(q_M[i,])
  S_obs[which_S_obs] ~ lognormal(log(S[which_S_obs]), tau_S[pop[which_S_obs]]);  // observed spawners
  // obs wild ocean age freq: n_MSage_obs[i,] ~ multinomial(q_MS[i,])
  target += sum(n_MSage_obs .* log(q_MS));
  // obs wild Gilbert-Rich age freq
  // if conditionGRonMS == 0: n_GRage_obs[i,] ~ multinomial(q_GR[i,])
  //                          P(freshwater age & ocean age)
  // if conditionGRonMS == 1: n_GRage_obs[i,oa] ~ multinomial(q_GR[i,oa]) where oa indexes ocean age
  //                          P(freshwater age | ocean age)
  target += sum(n_GRage_obs .* log(q_GR));
  n_H_obs ~ binomial(n_HW_obs, p_HOS);  // observed counts of hatchery vs. wild spawners
}

generated quantities {
  vector[RRS[1]*N_pop] delta_alpha;   // H vs. W discount in log intrinsic productivity
  vector[RRS[2]*N_pop] delta_Mmax;    // H vs. W discount in log maximum recruitment
  array[N_pop] corr_matrix[N_Mage-1] R_p_M; // correlation matrices of log-ratio smolt age distns
  array[N_pop] corr_matrix[N_Mage] R_MS; // correlation matrices of logit SAR by smolt age
  corr_matrix[N_Mage*(N_MSage-1)] R_p_MS[N_pop]; // correlation matrices of log-ratio ocean age distns
  vector[N] LL_M_obs;        // pointwise log-likelihood of smolts
  vector[N] LL_n_Mage_obs;   // pointwise log-likelihood of smolt age frequencies
  vector[N] LL_S_obs;        // pointwise log-likelihood of spawners
  vector[N] LL_n_MSage_obs;  // pointwise log-likelihood of ocean age frequencies
  vector[N] LL_n_GR_age_obs; // pointwise log-likelihood of Gilbert-Rich age frequencies
  vector[N_H] LL_n_H_obs;    // pointwise log-likelihood of hatchery vs. wild frequencies
  vector[N] LL;              // total pointwise log-likelihood                              

  // H vs. W discounts
  delta_alpha = log(alpha_H) - log(alpha_W);
  delta_Mmax = log(Mmax_H) - log(Mmax_W);
  
  // Correlation matrices 
  for(j in 1:N_pop)
  {
    R_p_M[j] = multiply_lower_tri_self_transpose(L_p_M[j]);
    R_MS[j] = multiply_lower_tri_self_transpose(L_MS[j]);
    R_p_MS[j] = multiply_lower_tri_self_transpose(L_p_MS[j]);
  }
  
  // Pointwise observation log-likelihoods
  LL_M_obs = rep_vector(0,N);
  for(i in 1:N_M_obs)
    LL_M_obs[which_M_obs[i]] = lognormal_lpdf(M_obs[which_M_obs[i]] | log(M[which_M_obs[i]]), tau_M); 
  LL_n_Mage_obs = (n_Mage_obs .* log(q_M)) * rep_vector(1,N_Mage);
  LL_S_obs = rep_vector(0,N);
  for(i in 1:N_S_obs)
    LL_S_obs[which_S_obs[i]] = lognormal_lpdf(S_obs[which_S_obs[i]] | log(S[which_S_obs[i]]), tau_S); 
  LL_n_MSage_obs = (n_MSage_obs .* log(q_MS)) * rep_vector(1,N_MSage);
  LL_n_GR_age_obs = (n_GRage_obs .* log(q_GR)) * rep_vector(1,N_GRage);
  LL_n_H_obs = rep_vector(0,N_H);
  for(i in 1:N_H)
    LL_n_H_obs[i] = binomial_lpmf(n_H_obs[i] | n_HW_obs[i], p_HOS[i]);
  LL = LL_M_obs + LL_n_Mage_obs + LL_S_obs + LL_n_MSage_obs + LL_n_GR_age_obs;
  LL[which_H] = LL[which_H] + LL_n_H_obs;
}
