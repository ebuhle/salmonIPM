#' Fit an IPM or run-reconstruction regression model
#'
#' Fits integrated population models (IPMs) or run-reconstruction regressions for anadromous 
#' Pacific salmonid population dynamics using a hierarchical Bayesian framework 
#' implemented in Stan. Various models are available, representing alternative 
#' life-history structures and independent or hierarchically pooled populations.
#' Users can specify stage-specific covariate effects and hyper-priors using formula syntax.
#'
#' @param stan_model Character string specifying the **salmonIPM** model to be
#'   fit. A more concise alternative to specifying `model`, `life_cycle`, and
#'   `pool_pops` and will override those arguments.
#' @param model Either `"IPM"` or `"RR"`, indicating whether to fit an IPM or
#'   run-reconstruction regression model.
#' @param life_cycle Character string indicating which life-cycle model to fit.
#'   One of the following options (must be `"SS"` if `model == "RR"`):
#'   * `"SS"`  Spawner-to-spawner (the default)
#'   * `"SSiter"` Spawner-to-spawner with iteroparity (an alias for `"SS"` with
#'   `stan_data("IPM_SS_[x]p", fish_data = fish_data)$iter` set to 1)
#'   * `"SMS"`  Spawner-smolt-spawner with a fixed smolt age
#'   * `"SMaS"`  Spawner-smolt-spawner with multiple smolt age classes (currently only
#'   available for `pool_pops == FALSE`)
#'   * `"LCRchum"`  Customized spawner-smolt-spawner model for Lower Columbia River chum
#'   (`pool_pops == TRUE`)
#' @param pool_pops Logical defaulting to `TRUE` if multiple populations are
#'   present in `fish_data` and `FALSE` otherwise, indicating whether to model
#'   multiple populations hierarchically rather than as independent "fixed
#'   effects". It is possible to fit a model to multiple populations
#'   simultaneously even though they share no parameters; indeed this is more
#'   efficient than fitting them one at a time because calculations are
#'   vectorized and warmup is shared.
#' @param ages For multi-stage models, a named list giving the ages in years of
#'   all fixed-age subadult life stages. This is not needed for `IPM_SMaS_np`
#'   because in that case smolt age structure is provided in `fish_data`.
#' @param SR_fun One of `"exp"` (density-independent discrete exponential),
#'   `"BH"` (Beverton-Holt, the default), or `"Ricker"`, indicating which
#'   spawner-recruit function to fit. Synonyms `"DI"`, `"B-H"`, `"bh"`, `"b-h"`
#'   and `"ricker"` are accepted.
#' @param RRS A character string or vector of strings naming parameters of the
#'   function specified by `SR_fun` that differ between wild- and
#'   hatchery-origin spawners, such that the relative reproductive success of
#'   hatchery spawners is not equal to 1. If `pool_pops == TRUE`, these should
#'   be the names of the population-specific parameters, **not** their
#'   hyper-means. For example, if `life_cycle %in% c("SS","SSiter")`, the
#'   options are `"none"` (the default), `"alpha"`, `"Rmax"`, or
#'   `c("alpha","Rmax")`.
#' @param par_models  Optional list of two-sided formulas of the form `theta ~
#'   x1 + ... + xK`, where `theta` is a (hyper)parameter or state in the model
#'   specified by `stan_model` that accepts covariates (see **Details** for
#'   available model-parameter combinations) and `x1 + ... + xK` are terms
#'   involving variables in `fish_data`. Standard formula syntax such as `:` and
#'   `*` may be used; see [stats::formula()].
#' @param center Logical indicating whether the terms in model matrices
#'   constructed from `fish_data` using the formulas in `par_models` should be
#'   centered. It is usually recommended to use the default (`TRUE`) so the
#'   baseline parameter estimate applies when predictors are at their sample
#'   means, but in some cases such as factor predictors `center = FALSE` may be
#'   appropriate. If combining categorical and numeric predictors, the latter
#'   can be centered and scaled prior to modeling.
#' @param scale  Logical indicating whether the model matrices constructed from
#'   `fish_data` using the formulas in `par_models` should be scaled to have
#'   column SDs of 1. Unit-scaling predictors is less critical than centering,
#'   but is advisable if variables have scales far from 1.
#' @param age_F Logical or 0/1 vector of length `N_age` indicating whether each
#'   adult age is fully (non)selected by the fishery. The default is all
#'   selected. If `life_cycle == "SSiter"`, `N_age` refers to the total number
#'   of maiden and repeat age classes (counting the repeat plus group as 1).
#' @param age_B Logical or 0/1 vector of length `N_age` indicating whether each
#'   adult age is fully (non)selected in broodstock collection. The default is
#'   all selected. If `life_cycle == "SSiter"`, `N_age` refers to the total
#'   number of maiden and repeat age classes (counting the repeat plus group as
#'   1).
#' @param age_S_obs If `stan_model == "IPM_SS_pp"`, a logical or 0/1 integer
#'   vector indicating, for each adult age, whether the observed total spawner
#'   data includes that age. The default is to treat `S_obs` as including
#'   spawners of all ages. This option may be useful if certain age classes are
#'   not counted. If `life_cycle == "SSiter"`, `N_age` refers to the total
#'   number of maiden and repeat age classes (counting the repeat plus group as
#'   1).
#' @param age_S_eff If `stan_model == "IPM_SS_pp"`, a logical or 0/1 vector
#'   indicating, for each adult age, whether spawners of that age contribute to
#'   reproduction. This can be used, e.g., to exclude jacks from the effective
#'   breeding population. The default is to include spawners of all ages. If
#'   `life_cycle == "SSiter"`, `N_age` refers to the total number of maiden and
#'   repeat age classes (counting the repeat plus group as 1).
#' @param conditionGRonMS If `life_cycle == "SMaS"`, logical indicating whether
#'   the Gilbert-Rich age frequencies `n_GRage_obs` in `fish_data` are
#'   conditioned on ocean age. If `FALSE` (the default) the counts are assumed
#'   to be sampled randomly from the population. If `TRUE`, it is assumed that
#'   the number of spawners of each ocean age (e.g., jacks vs 1-ocean) is
#'   arbitrary, but smolt (FW) age is randomly sampled within each ocean age;
#'   i.e., in a `smolt age x ocean age` contingency table, the cell frequencies
#'   are conditioned on the column totals.
#' @param priors Optional list of two-sided formulas of the form `theta ~
#'   distribution(params)`, where `theta` is a hyperparameter that can take a
#'   user-specified prior and `distribution()` is its canonical prior family.
#'   See [`priors`] for details on the available parameters in each model and
#'   their corresponding families. Any hyperparameters not given explicit priors
#'   will use the default values of the prior `params`.
#' @param fish_data Data frame where each row corresponds to a unique population
#'   `x` year, that includes the following `colnames` in no particular order
#'   except where noted (namely age- and origin-composition data):
#'   * `pop`  Required factor, numeric, or character population name or ID. Will be coerced to a
#'   factor, but it is recommended that this be a factor with concise,
#'   informative levels, e.g. `"Johnson Cr"`. This is especially true if there
#'   are multiple populations, in which case `levels(factor(pop))` can be used
#'   in interpreting and plotting the posterior draws.
#'   * `year`  Required numeric or integer giving the calendar year corresponding to each
#'   observation. Note that `fish_data` is not indexed by brood year. For a
#'   brood table run reconstruction, see [run_recon()].
#'   * `A`  Required spawning habitat size (either stream length or area). Will often be
#'   time-invariant within a population, but need not be. Habitat size is used
#'   internally to convert population size scaling parameters (e.g., `Rmax`)
#'   from density to abundance, so if `A == 1` no rescaling is done and these
#'   parameters are in units of fish. This is fine if `pool_pops == FALSE`, but
#'   in hierarchical multi-population models it is advisable to provide habitat
#'   size so that the assumption of hierarchical exchangeability is more valid.
#'   * `M_obs`  If `life_cycle %in% c("SMS","SMaS","LCRchum")`,
#'   the observed number of wild-origin smolts (`integer` or `numeric`). Missing
#'   / unknown observations are coded as `NA`.
#'   * `tau_M_obs`  If `life_cycle == "LCRchum"`,  known lognormal observation error SDs
#'   for smolt abundance. Missing values (`NA`) will be imputed.
#'   * `downstream_trap`  If `life_cycle == "LCRchum"`, row indices
#'   corresponding to a downstream smolt trap in a different population whose
#'   catch additionally includes the smolts produced in one or more upstream
#'   populations, assuming no extra mortality en route. Each upstream population
#'   can have at most one downstream trap (in addition to its own, if any) but a
#'   trap can have multiple upstream populations. If `downstream_trap[i] == j`,
#'   then `M_downstream[j] = M[j] + M[i]`. If `is.na(downstream_trap[i])` then
#'   `M[i]` is not double-counted.
#'   * `n_Mage[min_Mage]_obs...n_Mage[max_Mage]_obs`  If `life_cycle == "SMaS"`,
#'   multiple columns of observed smolt age sample frequencies (counts), where
#'   `[min_Mage]` and `[max_Mage]` are the numeral age in years of the youngest
#'   and oldest smolts, respectively. Note that age is measured in calendar
#'   years from the brood year (i.e., the Gilbert-Rich system).
#'   * `S_obs`  Required observed total escapement of all wild and hatchery-origin spawners
#'   (`integer` or `numeric`). Missing / unknown observations are coded as `NA`.
#'   * `tau_S_obs`  If `life_cycle == "LCRchum"`, known lognormal observation error SDs
#'   for spawner abundance. Missing values (`NA`) will be imputed.
#'   * `n_age[min_age]_obs ... n_age[max_age]_obs`  Integer columns of
#'   observed spawner age sample frequencies (counts), where `[min_age]` and
#'   `[max_age]` are the numeral age in years (total, not ocean age) of the
#'   youngest and oldest spawners. Required unless `life_cycle %in% c("SSiter","SMaS")`, 
#'   in which case it is ignored. Note that `n_age_obs` and all other compositional 
#'   data types must not contain `NA`. If the sample included no individuals of a 
#'   given category or if no samples were collected, the observed frequency is 0.
#'   * `n_age[min_age]M_obs ... n_age[max_age]M_obs n_age[min_age + 1]K_obs ... n_age[max_age + 1]K_obs`,
#'   If `life_cycle == "SSiter"`, integer columns of observed first-time
#'   (maiden) and repeat (former kelt) spawner age frequencies where `[min_age]`
#'   and `[max_age]` are the total age in years of the youngest and oldest
#'   **maiden** spawners, respectively. Contiguous maiden age columns denoted by
#'   `M` are followed by an equal number of contiguous repeat age columns
#'   denoted by `K`, where each repeat age is 1 year greater than the
#'   corresponding maiden age. The maximum repeat age class is a plus-group,
#'   i.e. it includes all repeat spawners age `max_age + 1` or older.
#'   * `n_MSage[min_MSage]_obs ... n_MSage[max_MSage]_obs`  If `life_cycle == "SMaS"`,
#'   integer columns of observed spawner ocean age sample frequencies, where
#'   `[min_MSage]` and `[max_MSage]` are the youngest and oldest ocean age in
#'   years, respectively. Nonzero ocean age frequencies are only required if
#'   `conditionGRonMS == TRUE` (the columns must be present in any case). If
#'   `conditionGRonMS == FALSE`, then `n_MSage_obs` represents **independent**
#'   samples, not simply the (implicit) ocean-age marginal totals of
#'   `n_GRage_obs`.
#'   * `n_GRage[min_age]_[min_Mage]_obs ... n_GRage[max_age]_[max_Mage]_obs`  If
#'   `life_cycle == "SMaS"`, integer columns of observed Gilbert-Rich age sample
#'   frequencies, varying fastest by smolt age (`min_Mage:max_Mage`) and then by
#'   total age (`min_age:max_age`). For example, a life history with subyearling
#'   or yearling smolts and ocean ages 2:3 would have column names
#'   `c("n_GRage_3_1_obs", "n_GRage_4_1_obs", "n_GRage_4_2_obs",
#'   "n_GRage_5_2_obs")`. All combinations of smolt age and (implicitly) ocean
#'   age must be represented, even if some were never observed.
#'   * `n_W_obs`  Required integer observed sample frequencies of natural-origin spawners.
#'   * `n_H_obs`  Required integer observed sample frequencies of hatchery-origin spawners.
#'   * `fit_p_HOS`  Logical or 0/1 indicating for each row `i` in `fish_data` whether the
#'   model should estimate `p_HOS[i] > 0`. Required if `model == "IPM"` unless
#'   `life_cycle == "LCRchum"`. [stan_data()] will give a warning if any row `i`
#'   meets either of two conditions: `as.logical(fit_p_HOS[i]) == FALSE` but
#'   `n_W_obs[i] + n_H_obs[i] > 0`, or `as.logical(fit_p_HOS[i]) == TRUE` but
#'   `n_W_obs[i] + n_H_obs[i] == 0`. The first means HOR were observed, so not
#'   accounting for them risks biasing the estimated parameters and states (aka
#'   "masking"). The second means the model is being asked to estimate
#'   `p_HOS[i]` with no case-specific hatchery / wild origin-frequency data.
#'   Because `p_HOS[i]` is an *a priori* independent parameter (a "fixed
#'   effect"), this is a difficult task. There may be some shared information
#'   via the process model to indirectly inform it, but this is likely to lead
#'   to poor estimates and sampling problems.
#'   * `n_O0_obs n_O[which_O_pop[1]]_obs ... n_O[which_O_pop[N_O_pop]]_obs`
#'   If `life_cycle = "LCRchum"`, multiple columns of observed origin sample
#'   frequencies. The first column, named "O" for origin and "0" for null /
#'   naught, refers to unknown natural origin, i.e. unmarked spawners presumed
#'   to be NOR. The next `N_O_pop` columns are numbered by the levels of
#'   `factor(fish_data$pop)` corresponding to the set of known-origin
#'   populations. Typically these are hatcheries, but NOR may be identified by
#'   PIT tags, parentage-based tagging, or other means. The `LCRchum` model uses
#'   origin-composition observations to infer the dispersal rates of hatchery
#'   (or other known-origin) fish, so `n_W_obs` (the same as `n_O0_obs` assuming
#'   all known origins are hatcheries) and `n_H_obs` (equal to
#'   `sum(n_O_obs[-1])` in that case) are not needed, although they can be
#'   included in `fish_data` for informational purposes or for post-processing
#'   draws. Likewise `fit_p_HOS` is not needed and will be ignored.
#'   * `n_M_obs`  If `life_cycle == "LCRchum"`, integer observed frequencies
#'   of male spawners.
#'   * `n_F_obs`  If `life_cycle == "LCRchum"`, integer observed frequencies of
#'   female spawners.
#'   * `p_G_obs`  If `life_cycle == "LCRchum"`, observed proportion (assumed known
#'   without error) of female spawners that are "green", i.e. fully fecund.
#'   * `F_rate`  Total harvest rate (proportion) of natural-origin fish. Required for 
#'   all models, even if all values are 0.
#'   * `B_take_obs`  Number of adults taken for hatchery broodstock. Required
#'   for all models, even if all values are 0.
#'   * `S_add_obs`  If `stan_model == "IPM_LCRchum_pp"`, number of adults translocated into
#'   population.
#'   * `...`  Additional variables to be used as covariates. These can vary spatially and/or
#'   temporally.
#' @param fecundity_data If `life_cycle == "LCRchum"`, data frame with the
#'   following columns, representing observations of fecundity with each row
#'   corresponding to a female:
#' * `age_E`  Female age in years.
#' * `E_obs`  Observed fecundity.
#' @param pars A character vector specifying (hyper)parameters, states, and/or
#'   quantities of interest ("parameters") to be saved. The default is to save
#'   all parameters. Parameters can be explicitly named or one or more shortcuts
#'   can be used to specify hierarchical levels of parameters; see [stan_pars()]
#'   for details. If parameters are explicitly named, the user should make sure
#'   they exist in `stan_model`, e.g. by calling `stan_pars(stan_model)`.
#' @param include Logical scalar defaulting to `TRUE` indicating whether to
#'   include or exclude the parameters given by `pars`. If `FALSE`, only entire
#'   multidimensional parameters can be excluded, rather than particular
#'   elements of them.
#' @param log_lik Logical scalar indicating whether the pointwise log-likelihood
#'   should be saved, e.g. for later use with [loo::loo()].
#' @param init A list of named lists of initial values to be passed to
#'   [rstan::sampling()]. The default (and recommended) `NULL` randomly
#'   generates initial values for each chain given the data and model using
#'   [stan_init()].
#' @param chains Positive integer specifying the number of HMC chains; see
#'   [rstan::sampling()].
#' @param iter Positive integer specifying the number of iterations for each
#'   chain (including warmup); see [rstan::sampling()].
#' @param warmup Positive integer specifying the number of warmup iterations per
#'   chain. If step-size adaptation is enabled (the default), this also controls
#'   the number of iterations for which adaptation is run; hence these warmup
#'   samples should not be used for inference. The number of warmup iterations
#'   should not be larger than `iter`. See [rstan::sampling()].
#' @param thin Positive integer specifying the period for saving samples. The
#'   default is 1, which is usually the recommended value. See
#'   [rstan::sampling()].
#' @param cores Number of cores to use when running chains in parallel. Defaults
#'   to the number of physical cores available. See [rstan::sampling()].
#' @param control A named list of options to control sampler behavior. See
#'   [rstan::stan()] for details and available options. In contrast to
#'   **rstan**, the default value of `adapt_delta` in **salmonIPM** is increased
#'   to 0.95 as we have found this necessary to minimize divergences in most
#'   cases.
#' @param save_data Logical scalar defaulting to `FALSE` indicating whether to
#'   save the data passed to Stan by [stan_data()] in the [salmonIPMfit] object.
#'   Can be useful for posterior predictive checking and reproducibility, at the
#'   cost of increased object size.
#' @param ... Additional arguments to pass to [rstan::sampling()].
#'
#' @details The following parameters and states can be modeled as functions of
#'   covariates using the argument `par_models`. The response distribution
#'   families and link functions are automatically implemented so there is no
#'   need to `log()`- or `qlogis()`-transform the left-hand side of the formula,
#'   although such syntax will also work (because the LHS is not found in
#'   `fish_data` and is just syntactic sugar to determine the parameter name).
#'   The design matrices passed to the Stan model cannot include an intercept,
#'   but it is not necessary to manually remove it in the RHS; if present by
#'   default, [par_model_matrix()] will automatically remove it.
#'
#'   As with any regression model, the user must ensure the effects specified
#'   are estimable given the design matrix. For example, the effect of a
#'   spatially varying but time-invariant predictor would not be identifiable in
#'   a `np` model because populations are modeled independently.
#' 
#' |                    |                         |                             |                        | **Response (family)**  |                    |                     |                              |                              |
#' |:-------------------|:-----------------------:|:---------------------------:|:----------------------:|:----------------------:|:------------------:|:-------------------:|:----------------------------:|:----------------------------:|    
#' | Model              | `alpha` \cr (lognormal) | `psi` \cr (logistic normal) | `Rmax` \cr (lognormal) | `Mmax` \cr (lognormal) | `R` \cr (lognormal)| `M` \cr (lognormal) | `s_MS` \cr (logistic normal) | `s_SS` \cr (logistic normal) |
#' | `IPM_SS_np`        | &#x2611;                | &#x2610;                    | &#x2611;               | &#x2610;               | &#x2611;           | &#x2610;            | &#x2610;                     | &#x2610;                     | 
#' | `IPM_SSiter_np`    | &#x2611;                | &#x2610;                    | &#x2611;               | &#x2610;               | &#x2611;           | &#x2610;            | &#x2610;                     | &#x2611;                     | 
#' | `IPM_SS_pp`        | &#x2611;                | &#x2610;                    | &#x2611;               | &#x2610;               | &#x2611;           | &#x2610;            | &#x2610;                     | &#x2610;                     | 
#' | `IPM_SSiter_pp`    | &#x2611;                | &#x2610;                    | &#x2611;               | &#x2610;               | &#x2611;           | &#x2610;            | &#x2610;                     | &#x2611;                     | 
#' | `IPM_SMS_np`       | &#x2611;                | &#x2610;                    | &#x2610;               | &#x2611;               | &#x2610;           | &#x2611;            | &#x2611;                     | &#x2610;                     | 
#' | `IPM_SMS_pp`       | &#x2611;                | &#x2610;                    | &#x2610;               | &#x2611;               | &#x2610;           | &#x2611;            | &#x2611;                     | &#x2610;                     | 
#' | `IPM_SMaS_np`      | &#x2611;                | &#x2610;                    | &#x2610;               | &#x2611;               | &#x2610;           | &#x2611;            | &#x2611;                     | &#x2610;                     | 
#' | `IPM_LCRchum_pp`   | &#x2610;                | &#x2611;                    | &#x2610;               | &#x2611;               | &#x2610;           | &#x2611;            | &#x2611;                     | &#x2610;                     | 
#'
#' @return An object of class `salmonIPMfit` representing the fitted model. See
#'   [`salmonIPMfit-class`] for details.
#'
#' @importFrom rstan sampling get_elapsed_time
#' @importFrom utils packageVersion 
#' @export
#' @encoding UTF-8

salmonIPM <- function(stan_model = paste(model, life_cycle, ifelse(pool_pops, "pp", "np"), sep = "_"), 
                      model = c("IPM","RR"), life_cycle = c("SS","SSiter","SMS","SMaS","LCRchum"), 
                      ages = NULL, pool_pops = nlevels(factor(fish_data$pop)) > 1, 
                      SR_fun = c("BH","B-H","bh","b-h","Ricker","ricker","exp"), RRS = "none", 
                      par_models = NULL, center = TRUE, scale = TRUE, 
                      age_F = NULL, age_B = NULL, age_S_obs = NULL, age_S_eff = NULL, 
                      conditionGRonMS = FALSE, priors = NULL, 
                      fish_data, fecundity_data = NULL, 
                      pars = "all", include = TRUE, log_lik = FALSE, 
                      init = NULL, chains = 4, iter = 2000, warmup = floor(iter/2), thin = 1, 
                      cores = parallel::detectCores(logical = FALSE), 
                      control = NULL, save_data = FALSE, ...)
{
  model <- match.arg(model)
  life_cycle <- match.arg(life_cycle)
  force(stan_model)
  mlp <- strsplit(stan_model, "_")[[1]]
  model <- mlp[1]             # override these args if stan_model is specified
  life_cycle <- mlp[2]
  pool_pops <- mlp[3] == "pp"
  stanmodel <- gsub("iter", "", stan_model) # the same Stan code handles semel/iteroparity 
  SR_fun <- match.arg(SR_fun)
  if(SR_fun == "DI") SR_fun <- "exp"
  if(SR_fun %in% c("B-H","bh","b-h")) SR_fun <- "BH"
  if(SR_fun == "ricker") SR_fun <- "Ricker"

  .call <- as.list(match.call(expand.dots = TRUE))
  for(n in grep("_data", names(.call), value = TRUE, invert = TRUE)) 
    .call[[n]] <- eval(.call[[n]])
  .call <- as.call(.call)
  
  names(par_models) <- lapply(par_models, function(f) all.vars(f)[1]) # names are responses

  dat <- stan_data(stan_model = stan_model, SR_fun = SR_fun, RRS = RRS, ages = ages, 
                   par_models = par_models, center = center, scale = scale, priors = priors, 
                   fish_data = fish_data, age_F = age_F, age_B = age_B, 
                   age_S_obs = age_S_obs, age_S_eff = age_S_eff, 
                   conditionGRonMS = conditionGRonMS, fecundity_data = fecundity_data)
  
  pars <- include_pars(pars = pars, stan_model = stan_model, SR_fun = SR_fun, RRS = RRS, 
                       par_models = par_models, include = include, log_lik = log_lik)
  hyper <- stan_pars(stan_model = stan_model, pars = "hyper", SR_fun = SR_fun, 
                     RRS = RRS, par_models = par_models)
  prior.info <- get_prior_info(stan_data = dat, stanmodel = stanmodels[[stanmodel]], pars = hyper)

  if(is.null(init)) 
    init <- stan_init(stan_model = stan_model, stan_data = dat, chains = chains)
  if(is.null(control$adapt_delta)) 
    control$adapt_delta <- 0.95
  
  stanfit <- rstan::sampling(stanmodels[[stanmodel]], 
                             data = dat, init = init, pars = pars,
                             chains = chains, cores = cores, 
                             iter = iter, warmup = warmup, thin = thin, 
                             control = control, ...)
  
  out <- salmonIPMfit(stanfit = stanfit, call = .call, stan_model = stan_model,
                      model = model, life_cycle = life_cycle, ages = ages, 
                      pool_pops = pool_pops, SR_fun = SR_fun, RRS = RRS, 
                      par_models = par_models, center = center, scale = scale, 
                      prior.info = prior.info, age_S_obs = age_S_obs, 
                      age_S_eff = age_S_eff, conditionGRonMS = conditionGRonMS,
                      stan_data = if(save_data) dat else NULL,
                      dims = list(N = dat$N, N_pop = max(dat$pop), N_year = max(dat$year)),
                      pops = levels(factor(fish_data$pop)),
                      elapsed_time = get_elapsed_time(stanfit),
                      salmonIPM_version = packageVersion("salmonIPM"))
  
  return(out)
}
