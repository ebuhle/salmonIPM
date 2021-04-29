functions {
  // spawner-recruit functions
  real SR(int SR_fun, real alpha, real Rmax, real S, real A) {
    real R;
    
    if(SR_fun == 1)      // discrete exponential
      R = alpha*S;
    else if(SR_fun == 2) // Beverton-Holt
      R = alpha*S/(1 + alpha*S/(A*Rmax));
    else if(SR_fun == 3) // Ricker
      R = alpha*S*exp(-alpha*S/(A*e()*Rmax));
    
    return(R);
  }
  
  // Generalized normal (aka power-exponential) unnormalized log-probability
  real pexp_lpdf(real y, real mu, real sigma, real shape) {
    return(-(fabs(y - mu)/sigma)^shape);
  }
    
  // Quantiles of a vector
  real quantile(vector v, real p) {
    int N = num_elements(v);
    real Np = round(N*p);
    real q;
    
    for(i in 1:N) {
      if(i - Np == 0.0) q = v[i];
    }
    return(q);
  }
}

data {
  int<lower=1> SR_fun;          // S-R model: 1 = exponential, 2 = BH, 3 = Ricker
  int<lower=1> N;               // total number of cases in all pops and years
  int<lower=1,upper=N> pop[N];  // population identifier
  int<lower=1,upper=N> year[N]; // brood year identifier
  int<lower=1,upper=N> N_fit;   // number of cases used in fitting (non-missing S and R)
  int<lower=1,upper=N> which_fit[N_fit]; // cases used in fitting
  vector<lower=0>[N] S;         // observed annual total spawner abundance (not density)
  vector<lower=0>[N] R;         // total natural recruit abundance (not density), including harvest and broodstock removals
  vector[N] A;                  // habitat area associated with each spawner abundance obs
  int<lower=0,upper=1> S_NA[N]; // logical indicating whether S is missing and should be simulated
  int<lower=0,upper=1> R_NA[N]; // logical indicating whether R is missing and should be simulated
  int<lower=2> N_age;           // number of adult age classes
  int<lower=2> max_age;         // maximum adult age
  matrix<lower=0,upper=1>[max(pop),N_age] p;  // average recruit age distributions for each pop 
}

transformed data {
  int<lower=1,upper=N> N_pop;   // number of populations
  int<lower=1,upper=N> N_year;  // number of years
  int<lower=2> ages[N_age];     // adult ages
  real mu_mu_Rmax = quantile(log(R[which_fit]), 0.9);  // prior mean of mu_Rmax
  real sigma_mu_Rmax = sd(log(R[which_fit]));  // prior SD of mu_Rmax
  
  N_pop = max(pop);
  N_year = max(year);
  for(a in 1:N_age)
    ages[a] = max_age - N_age + a;
}

parameters {
  real mu_alpha;                        // hyper-mean log intrinsic productivity of wild spawners
  real<lower=0> sigma_alpha;            // hyper-SD log intrinsic productivity
  vector[N_pop] zeta_alpha;             // log intrinsic prod of wild spawners (Z-scores)
  real mu_Rmax;                         // hyper-mean log asymptotic recruitment
  real<lower=0> sigma_Rmax;             // hyper-SD log asymptotic recruitment
  vector[N_pop] zeta_Rmax;              // log asymptotic recruitment (Z-scores)
  real<lower=-1,upper=1> rho_alphaRmax; // correlation between log(alpha) and log(Rmax)
  real<lower=-1,upper=1> rho_phi;       // AR(1) coef for brood year log productivity anomalies
  real<lower=0> sigma_phi;              // hyper-SD of brood year log productivity anomalies
  vector[max(year)] zeta_phi;           // log brood year productivity anomalies (Z-scores)
  real<lower=0> sigma;                  // residual error SD
}

transformed parameters {
  vector<lower=0>[N_pop] alpha; // intrinsic productivity 
  vector<lower=0>[N_pop] Rmax;  // asymptotic recruitment 
  vector<lower=0>[N_year] phi;  // log brood year productivity anomalies
  vector<lower=0>[N] R_hat;     // expected recruit abundance (not density) by brood year
  
  // Multivariate Matt trick for [log(alpha), log(Rmax)]
  {
    matrix[2,2] L_alphaRmax;           // Cholesky factor of corr matrix of log(alpha), log(Rmax)
    matrix[N_pop,2] zeta_alphaRmax;    // [log(alpha), log(Rmax)] random effects (z-scored)
    matrix[N_pop,2] epsilon_alphaRmax; // [log(alpha), log(Rmax)] random effects
    vector[2] sigma_alphaRmax;         // SD vector of [log(alpha), log(Rmax)]
    
    L_alphaRmax[1,1] = 1;
    L_alphaRmax[2,1] = rho_alphaRmax;
    L_alphaRmax[1,2] = 0;
    L_alphaRmax[2,2] = sqrt(1 - rho_alphaRmax^2);
    sigma_alphaRmax[1] = sigma_alpha;
    sigma_alphaRmax[2] = sigma_Rmax;
    zeta_alphaRmax = append_col(zeta_alpha, zeta_Rmax);
    epsilon_alphaRmax = diag_pre_multiply(sigma_alphaRmax, L_alphaRmax * zeta_alphaRmax')';
    alpha = exp(mu_alpha + epsilon_alphaRmax[,1]);
    Rmax = exp(mu_Rmax + epsilon_alphaRmax[,2]);
  }
  
  // AR(1) model for phi
  phi[1] = zeta_phi[1]*sigma_phi/sqrt(1 - rho_phi^2); // initial anomaly
  for(i in 2:N_year)
    phi[i] = rho_phi*phi[i-1] + zeta_phi[i]*sigma_phi;
  phi = phi - mean(phi);  // constrain log anomalies to sum to zero

  // Predict recruitment
  R_hat = rep_vector(0,N);
  for(i in 1:N_fit)
    R_hat[which_fit[i]] = SR(SR_fun, alpha[pop[which_fit[i]]], Rmax[pop[which_fit[i]]], 
                             S[which_fit[i]], A[which_fit[i]]);
}

model {
  // Priors
  mu_alpha ~ normal(0,5);
  sigma_alpha ~ pexp(0,3,10);
  mu_Rmax ~ normal(mu_mu_Rmax, sigma_mu_Rmax);
  sigma_Rmax ~ normal(0,3);
  rho_phi ~ pexp(0,0.85,50);  // mildly regularize to ensure stationarity
  sigma_phi ~ normal(0,3);
  sigma ~ normal(0,2);
  
  // Hierarchical priors
  // [log(alpha), log(Rmax)] ~ MVN(0, D*R_log_aRmax*D), where D = diag_matrix(sigma_alpha, sigma_Rmax)
  zeta_alpha ~ std_normal();
  zeta_Rmax ~ std_normal();
  zeta_phi ~ std_normal();    // phi ~ N(0, sigma_phi)
  
  // Likelihood
  R[which_fit] ~ lognormal(log(R_hat[which_fit]) + phi[year[which_fit]], sigma);
}

generated quantities {
  vector[N] S_sim;    // simulated spawners
  vector[N] R_sim;    // simulated recruits

  S_sim = S;
  R_sim = R;

  for(i in 1:N)
  {
    if(S_NA[i] == 1)
    {
      if(i >= max_age && pop[i-max_age] == pop[i])
      {
        S_sim[i] = 0;
        for(a in 1:N_age)
          S_sim[i] = S_sim[i] + R_sim[i-ages[a]]*p[pop[i],a];
      }
    }

    if(R_NA[i] == 1)
      R_sim[i] = SR(SR_fun,alpha[pop[i]], Rmax[pop[i]], S_sim[i], A[i]) * lognormal_rng(phi[year[i]], sigma);
  }
}
