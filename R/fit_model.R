#' Fit Bayesian model for estimating effects on the *growth rate* of BANS burrow
#' counts between years (assuming burrow counts are accurate)
#'
#' @param inputdata Inputdata created by running \code{\link{setup_BBS_model}}
#' @param n.adapt Number of iterations for adaptation, defaults to 1000
#' @param n.burnin Number of iterations for burn-in, defaults to 5000
#' @param n.sample Number of iterations for sampling, defauts to 10000
#' @param n.chains Number of MCMC chains, defaults to 3
#' @param holdout Number of years to hold out of model fitting, and generate
#'   separate predictions for from the model. The selected number will be held
#'   out from the end of the time series. Defaults to 5.
#' @param vars Character string listing the variables in the model to be
#'   monitored.
#' @param type Model type; one of: 
#'   - simple_pois_od (exponential growth, poisson error, overdispersion)
#'   - simple_negbin_od (exponential growth, negative binomial error, overdispersion)
#'   - simple_negbin (exponential growth, negative binomial error, no overdispersion term)
#'   - thetalog_pois_varK (density-dependence with uncertain theta, poisson error, variable carrying capacity depending on predictor 1),
#'   - thetalog_negbin_varK (density-dependence with uncertain theta, negative binomial error, variable carrying capacity depending on predictor 1),
#'   - thetalog_negbin (density-dependence with uncertain theta, negative binomial error)
#'   - ricker (density-dependence with theta = 1)
#' @param ... Additional arguments passed to \code{\link[rjags]{jags.model}}.
#'
#' @return Returns a coda.samples object from rjags
#'
#' @export
#' @import rjags
#' @importFrom stats update
#' @importFrom MCMCvis MCMCsummary MCMCtrace
#'   
fit_BANSpopgrowth_model <- function(inputdata, n.adapt = 1000, n.burnin = 5000, 
                             n.sample = 10000, n.chains = 3, holdout = 5, 
                             type = 'growth_simple',
                             vars = c('N', 'beta', 'ypred'), 
                             thin = 1,
                             ...) {
  if (type == 'growth_logistic_simple') {
    # constant effect of density-dependence
    
    idat = inputdata
    idat$fityears = length(idat$N)
    idat$nBcov = dim(idat$Bpredictors)[2]
    idat$nKcov = dim(idat$Kpredictors)[2]
    idat$npred = dim(idat$predmatrix)[1]
    
    if (holdout > 0) {
      startpred = idat$fityears - holdout + 1
      idat$y[startpred:idat$fityears] = NA #replace hold out data with NA
    }
    
    modelstring = "
  model {

    ## PRIORS
    var.o ~ dunif(0, 2)
    tau.o = 1/var.o
    
    r0 ~ dnorm(0, 1/100)
    k ~ dnorm(0, 1/100)
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }
    var.p ~ dunif(0, 2)
    tau.p = 1/var.p

    ## LIKELIHOOD: 
    # y = observed growth rate from t to t+1
    # r0 = intrinsic growth rate (when predictors are at zero/mean values)
    # N = # burrows in year t (assumed measured without error)
    # k = strength of density dependence (~ -r0/K), where K =
    # carrying capacity; assumed constant across all years
    
    for (t in 1:fityears){
    
      y[t] ~ dnorm(mu[t], tau.o) #error in observation
      
      mu[t] = r0 + k * N[t] +
        beta[1] * Bpredictors[t, 1] + 
        beta[2] * Bpredictors[t, 2] + 
        beta[3] * Bpredictors[t, 3] + 
        beta[4] * Bpredictors[t, 4] + 
        e[t]
        
      # overdisperson by observation/year
      e[t] ~ dnorm(0, tau.p) #error in process model

      # simulated growth rates for the same 21 years as in inputdat$y, including
      # for years where y[t] wasn't observed (2005-06, 2020-21), but not where N
      # wasn't observed (2006-07)
      ysim[t] ~ dnorm(mu[t], tau.o)
      
      # squared error (in growth rates)
      sq[t] = (y[t] - mu[t])^2
      sqsim[t] = (ysim[t] - mu[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (for growth rates)
    mean.y = mean(y[])
    mean.ysim = mean(ysim[])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise

    sd.y = sd(y[])
    sd.ysim = sd(ysim[])
    pvalue.sd = step(sd.ysim - sd.y)

    fit = sum(sq[])
    fitsim = sum(sqsim[])
    pvalue.fit = step(fitsim - fit)
    
    ## PREDICTED BURROW NUMBERS:
    Nsim[1] = N_obs[1] #start with observed value
    
    # 2000 through 2006 (i=2 to i=8):
    # --> include prediction of unobserved burrow numbers for 2006 from
    # observations in 2005 and estimated growth rate:
    #    Nsim[2006] = N[2005] * (1 + ysim[2005-06])
    
    for (i in 2:8) {
      Nsim[i] = N_obs[i-1] * (1 + ysim[i-1])
    }
    
    # # 2007 (i=9) use observed value to fill in
    # Nsim[9] = N_obs[9]
    
    # for 2007, prediction of burrows depends on unknown number of
    # burrows in 2006 and unknown growth rate from 2006 to 2007:
    #    Nsim[2007] = N[2006] * (1 + ysim[2006-07])
    # --> so use predicted value for burrows in 2006 instead, and estimate
    # unobserved growth from 2006 to 2007 from predicted burrows in 2006 and
    # observed burrows in 2007 (can't use model above because no N observation
    # in 2006 to account for density dependence)
    #    ysim[2006-07] = (N[2007] - N[2006]) / N[2006]

    ysim.2006 = (N_obs[9] - Nsim[8]) / Nsim[8]

    Nsim[9] = Nsim[8] * (1 + ysim.2006) #will end up being the same as N[2007] (confounded)

    # 2008 through 2021:
    # use i-2 for ysim because skipping over 2006 in ysim vector
    for (i in 10:(fityears + 2)) {
      Nsim[i] = N_obs[i-1] * (1 + ysim[i-2])
    }
    
    # ## GROWTH RATE VS. PREDICTORS
    # # all other predictors held at mean value/0, including # of burrows, thus no
    # # additional effect of density dependence
    # 
    # for (i in 1:npred) {
    #   # e.pred[i] ~ dnorm(0, tau.p) #error in process model (treat as same across predictors)
    #   # k.pred[i] ~ dnorm(mu.k.pred[i], tau.k)
    #   # mu.k.pred[i] = k0
    #   
    #   ypred[i, 1] ~ dnorm(mu.y.pdsi[i], tau.o) #error in observation
    #   mu.y.aflow[i] = r0 + beta[1] * predmatrix[i,1] 
    #   
    #   ypred[i, 2] ~ dnorm(mu.y.pdsi[i], tau.o) #error in observation
    #   mu.y.pdsi[i] = r0 + beta[2] * predmatrix[i,2] 
    #   
    #   ypred[i, 3] ~ dnorm(mu.y.t[i], tau.o) #error in observation
    #   mu.y.t[i] = r0 + beta[3] * predmatrix[i,3] 
    #   
    #   # add effect on carrying capacity/strength of density dependence
    #   ypred[i, 4] ~ dnorm(mu.dd.pred[i], tau.k)
    #   mu.dd.pred[i] = k0 + kappa[1] * predmatrix[i,4]
    # }

    
  }"
  } else if (type == 'growth_logistic_vark') {
    # density-dependence, but let strength of effect vary randomly across years
    
    idat = inputdata
    idat$fityears = length(idat$N)
    idat$nBcov = dim(idat$Bpredictors)[2]
    idat$nKcov = dim(idat$Kpredictors)[2]
    idat$npred = dim(idat$predmatrix)[1]
    
    if (holdout > 0) {
      startpred = idat$fityears - holdout + 1
      idat$y[startpred:idat$fityears] = NA #replace hold out data with NA
    }
    
    modelstring = "
  model {

    ## PRIORS
    var.o ~ dunif(0, 2)
    tau.o = 1/var.o
    
    r0 ~ dnorm(0, 1/100)
    k0 ~ dnorm(0, 1/100)
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }
    var.p ~ dunif(0, 2)
    tau.p = 1/var.p
    
    var.k ~ dunif(0, 2)
    tau.k = 1/var.k

    ## LIKELIHOOD: 
    # y = observed growth rate from t to t+1
    # r0 = intrinsic growth rate (when predictors are at zero/mean values)
    # N = # burrows in year t (assumed measured without error)
    # k = strength of density dependence in year t (~ -r0/K), where K =
    # carrying capacity, and k varies randomly
    
    for (t in 1:fityears){
    
      y[t] ~ dnorm(mu[t], tau.o) #error in observation
      
      mu[t] = r0 + k[t] * N[t] +
        beta[1] * Bpredictors[t, 1] + 
        beta[2] * Bpredictors[t, 2] + 
        beta[3] * Bpredictors[t, 3] + 
        beta[4] * Bpredictors[t, 4] + 
        e[t]
      
      # annual random variation in k
      k[t] ~ dnorm(k0, tau.k)
      
      # overdisperson by observation/year
      e[t] ~ dnorm(0, tau.p) #error in process model

      # simulated growth rates for the same 21 years as in inputdat$y, including
      # for years where y[t] wasn't observed (2005-06, 2020-21), but not where N
      # wasn't observed (2006-07)
      ysim[t] ~ dnorm(mu[t], tau.o)
      
      # squared error (in growth rates)
      sq[t] = (y[t] - mu[t])^2
      sqsim[t] = (ysim[t] - mu[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (for growth rates)
    mean.y = mean(y[])
    mean.ysim = mean(ysim[])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise

    sd.y = sd(y[])
    sd.ysim = sd(ysim[])
    pvalue.sd = step(sd.ysim - sd.y)

    fit = sum(sq[])
    fitsim = sum(sqsim[])
    pvalue.fit = step(fitsim - fit)
    
    ## PREDICTED BURROW NUMBERS:
    Nsim[1] = N_obs[1] #start with observed value
    
    # 2000 through 2006 (i=2 to i=8):
    # --> include prediction of unobserved burrow numbers for 2006 from
    # observations in 2005 and estimated growth rate:
    #    Nsim[2006] = N[2005] * (1 + ysim[2005-06])
    
    for (i in 2:8) {
      Nsim[i] = N_obs[i-1] * (1 + ysim[i-1])
    }
    
    # # 2007 (i=9) use observed value to fill in
    # Nsim[9] = N_obs[9]
    
    # for 2007, prediction of burrows depends on unknown number of
    # burrows in 2006 and unknown growth rate from 2006 to 2007:
    #    Nsim[2007] = N[2006] * (1 + ysim[2006-07])
    # --> so use predicted value for burrows in 2006 instead, and estimate
    # unobserved growth from 2006 to 2007 from predicted burrows in 2006 and
    # observed burrows in 2007 (can't use model above because no N observation
    # in 2006 to account for density dependence)
    #    ysim[2006-07] = (N[2007] - N[2006]) / N[2006]

    ysim.2006 = (N_obs[9] - Nsim[8]) / Nsim[8]

    Nsim[9] = Nsim[8] * (1 + ysim.2006) #will end up being the same as N[2007] (confounded)

    # 2008 through 2021:
    # use i-2 for ysim because skipping over 2006 in ysim vector
    for (i in 10:(fityears + 2)) {
      Nsim[i] = N_obs[i-1] * (1 + ysim[i-2])
    }
    
    ## GROWTH RATE VS. PREDICTORS
    # all other predictors held at mean value/0, including # of burrows, thus no
    # additional effect of density dependence

    for (i in 1:npred) {
      e.pred[i] ~ dnorm(0, tau.p) #error in process model (treat as same across predictors)
      k.pred[i] ~ dnorm(k0, tau.k)

      ypred[i, 1] ~ dnorm(mu.y.aflow[i], tau.o) #error in observation
      mu.y.aflow[i] = r0 + beta[1] * predmatrix[i,1] + k.pred[i]*mean(N[]) + e.pred[i]
      
      ypred[i, 2] ~ dnorm(mu.y.cflow[i], tau.o) #error in observation
      mu.y.cflow[i] = r0 + beta[2] * predmatrix[i,2] + k.pred[i]*mean(N[]) + e.pred[i]

      ypred[i, 3] ~ dnorm(mu.y.pdsi[i], tau.o) #error in observation
      mu.y.pdsi[i] = r0 + beta[3] * predmatrix[i,3] + k.pred[i]*mean(N[]) + e.pred[i]

      ypred[i, 4] ~ dnorm(mu.y.t[i], tau.o) #error in observation
      mu.y.t[i] = r0 + beta[4] * predmatrix[i,4] + k.pred[i]*mean(N[]) + e.pred[i]
    }

    
  }"
  } else if (type == 'growth_logistic_predk') {
    # density-dependence, with strength of effect predicted by an environmental covariate

    
    idat = inputdata
    idat$fityears = length(idat$N)
    idat$nBcov = dim(idat$Bpredictors)[2]
    idat$nKcov = dim(idat$Kpredictors)[2]
    idat$npred = dim(idat$predmatrix)[1]
    
    if (holdout > 0) {
      startpred = idat$fityears - holdout + 1
      idat$y[startpred:idat$fityears] = NA #replace hold out data with NA
    }
    
    modelstring = "
  model {

    ## PRIORS
    var.o ~ dunif(0, 2)
    tau.o = 1/var.o
    
    r0 ~ dnorm(0, 1/100)
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }
    var.p ~ dunif(0, 2)
    tau.p = 1/var.p
    
    k0 ~ dunif(-5, 5) #in tens of thousands
    for (i in 1:nKcov) {
      kappa[i] ~ dnorm(0, 1/100)
    }
    var.k ~ dunif(0, 2)
    tau.k = 1/var.k

    ## LIKELIHOOD: 
    # y = observed growth rate from t to t+1
    # r0 = intrinsic growth rate (when predictors are at zero/mean values)
    # N = # burrows in year t (assumed measured without error)
    # k = strength of density dependence in year t (~ -r0/K), where K =
    # carrying capacity, and k is also affected by environmental factors
    
    for (t in 1:fityears){
    
      y[t] ~ dnorm(mu[t], tau.o) #error in observation
      
      mu[t] = r0 + k[t] * N[t] +
        beta[1] * Bpredictors[t, 1] + 
        beta[2] * Bpredictors[t, 2] + 
        beta[3] * Bpredictors[t, 3] + 
        e[t]
        
      # overdisperson by observation/year
      e[t] ~ dnorm(0, tau.p) #error in process model

      # carrying capacity (strength of density dependence) varies with predictor
      # kappa = slope of linear trend effect on carrying capacity
      # k0 = intercept; strength of density dependence when predictor = mean value
      k[t] ~ dnorm(mu.k[t], tau.k)
      mu.k[t] = k0 + kappa[1] * Kpredictors[t, 1]

      # simulated growth rates for the same 21 years as in inputdat$y, including
      # for years where y[t] wasn't observed (2005-06, 2020-21), but not where N
      # wasn't observed (2006-07)
      ysim[t] ~ dnorm(mu[t], tau.o)
      
      # squared error (in growth rates)
      sq[t] = (y[t] - mu[t])^2
      sqsim[t] = (ysim[t] - mu[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (for growth rates)
    mean.y = mean(y[])
    mean.ysim = mean(ysim[])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise

    sd.y = sd(y[])
    sd.ysim = sd(ysim[])
    pvalue.sd = step(sd.ysim - sd.y)

    fit = sum(sq[])
    fitsim = sum(sqsim[])
    pvalue.fit = step(fitsim - fit)
    
    ## PREDICTED BURROW NUMBERS:
    Nsim[1] = N_obs[1] #start with observed value
    
    # 2000 through 2006 (i=2 to i=8):
    # --> include prediction of unobserved burrow numbers for 2006 from
    # observations in 2005 and estimated growth rate:
    #    Nsim[2006] = N[2005] * (1 + ysim[2005-06])
    
    for (i in 2:8) {
      Nsim[i] = N_obs[i-1] * (1 + ysim[i-1])
    }
    
    # # 2007 (i=9) use observed value to fill in
    # Nsim[9] = N_obs[9]
    
    # for 2007, prediction of burrows depends on unknown number of
    # burrows in 2006 and unknown growth rate from 2006 to 2007:
    #    Nsim[2007] = N[2006] * (1 + ysim[2006-07])
    # --> so use predicted value for burrows in 2006 instead, and estimate
    # unobserved growth from 2006 to 2007 from predicted burrows in 2006 and
    # observed burrows in 2007 (can't use model above because no N observation
    # in 2006 to account for density dependence)
    #    ysim[2006-07] = (N[2007] - N[2006]) / N[2006]

    ysim.2006 = (N_obs[9] - Nsim[8]) / Nsim[8]

    Nsim[9] = Nsim[8] * (1 + ysim.2006) #will end up being the same as N[2007] (confounded)

    # 2008 through 2021:
    # use i-2 for ysim because skipping over 2006 in ysim vector
    for (i in 10:(fityears + 2)) {
      Nsim[i] = N_obs[i-1] * (1 + ysim[i-2])
    }
    
    ## GROWTH RATE VS. PREDICTORS
    # all other predictors held at mean value/0, including # of burrows, thus no
    # additional effect of density dependence
    
    for (i in 1:npred) {
      # e.pred[i] ~ dnorm(0, tau.p) #error in process model (treat as same across predictors)
      # k.pred[i] ~ dnorm(mu.k.pred[i], tau.k)
      # mu.k.pred[i] = k0
      
      ypred[i, 1] ~ dnorm(mu.y.pdsi[i], tau.o) #error in observation
      mu.y.aflow[i] = r0 + beta[1] * predmatrix[i,1] 
      
      ypred[i, 2] ~ dnorm(mu.y.pdsi[i], tau.o) #error in observation
      mu.y.pdsi[i] = r0 + beta[2] * predmatrix[i,2] 
      
      ypred[i, 3] ~ dnorm(mu.y.t[i], tau.o) #error in observation
      mu.y.t[i] = r0 + beta[3] * predmatrix[i,3] 
      
      # add effect on carrying capacity/strength of density dependence
      ypred[i, 4] ~ dnorm(mu.dd.pred[i], tau.k)
      mu.dd.pred[i] = k0 + kappa[1] * predmatrix[i,4]
    }

    
  }"
  }
  
  
  jm = rjags::jags.model(file = textConnection(modelstring),
                         data = idat,
                         n.adapt = n.adapt,
                         n.chains = n.chains,
                         ...)
  update(jm, n.iter = n.burnin)
  results = rjags::coda.samples(jm, 
                                variable.name = vars,
                                n.iter = n.sample,
                                thin = thin)
  return(results)
}



plot_model_predictions = function(mod, year, obsdat, scale = 1, center = 0) {
  # pull other stats:
  predvalues = full_join(MCMCvis::MCMCpstr(mod, params = c('ysim', 'ysim.2006', 'Nsim'), 
                                           func = function(x) median(x)) %>%
                           purrr::map_df(as_tibble, .id = 'metric') %>%
                           mutate(year = year),
                         MCMCvis::MCMCpstr(mod, params = c('ysim', 'ysim.2006', 'Nsim'),
                                           func = function(x) HDInterval::hdi(x, 0.95)) %>%
                           purrr::map_df(as_tibble, .id = 'metric') %>%
                           mutate(year = year),
                         by = c('metric', 'year')) %>% 
    # backtransform
    mutate_at(vars(value, lower, upper),
              ~case_when(metric == 'Nsim' ~ . * scale + center,
                         TRUE ~ .)) %>% 
    mutate(metric = recode(metric, ysim.2006 = 'ysim')) %>% 
    complete(year = 1999:2021, nesting(metric)) %>% 
    mutate(type = 'predicted',
           name = if_else(metric %in% c('ysim', 'ysim.2006'), 
                          'per-capita growth rate', 
                          'number of burrows'))

  p = bind_rows(predvalues %>% select(year, name, value:upper, type), 
                obsdat %>% select(year, pgr, burrows) %>% 
                  pivot_longer(-year, names_to = 'metric') %>% 
                  mutate(name = recode(metric, 'pgr' = 'per-capita growth rate',
                                       'burrows' = 'number of burrows'),
                         type = 'observed') %>% 
                  select(year, name, value, type)) %>% 
    ggplot(aes(year, value)) + 
    facet_wrap(~name, scales = 'free_y', nrow = 2) +
    geom_ribbon(data = . %>% filter(type == 'predicted'),
                aes(ymin = lower, ymax = upper), fill = 'gray80') +
    geom_line(aes(color = type)) + geom_point(aes(color = type)) +
    scale_color_grey(end = 0.5) + 
    labs(x = NULL, y = NULL)

  return(p)
}

plot_cov_relationships = function(mod, xvals, nm, 
                                  obsdat, scale, center) {
  # pull other stats:
  predvalues = bind_cols(MCMCvis::MCMCpstr(mod, params = c('ypred'), 
                                           func = function(x) mean(x)) %>%
                           purrr::map_df(as_tibble) %>% 
                           set_names(paste0(nm, '.median')),
                         MCMCvis::MCMCpstr(mod, params = c('ypred'),
                                           func = function(x) HDInterval::hdi(x, 0.95)) %>%
                           purrr::map_df(as_tibble) %>% 
                           set_names(paste(rep(nm, 2), 
                                           rep(c('lower', 'upper'), each = 4),
                                           sep = '.'))
  ) %>% 
    mutate(index = c(1:nrow(.))) %>% 
    pivot_longer(-index) %>% 
    separate(name, into = c('cov', 'metric'), sep = '\\.') %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    arrange(cov, index) %>% 
    left_join(xvals, by = c('index', 'cov' = 'name')) %>% 
    mutate(cov = recode(cov, 'pdsi' = 'Annual mean Palmer Drought Severity Index, Apr-Aug',
                        'flow14' = 'Annual stream flow > 14,000 cfs (millions)',
                        'flow14_3' = '3-year cumulative stream flow > 14,000 cfs (millions)',
                        't' = 'Long-term trend'))

  p = ggplot(predvalues, aes(value, median)) + 
    facet_wrap(~cov, scales = 'free_x', nrow = 2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') +
    geom_line() + geom_point() + 
    geom_point(data = obsdat %>% 
                 mutate(cov = recode(cov, 
                                     'pdsi' = 'Annual mean Palmer Drought Severity Index, Apr-Aug',
                                     'flow14' = 'Annual stream flow > 14,000 cfs (millions)',
                                     'flow14_3' = '3-year cumulative stream flow > 14,000 cfs (millions)',
                                     't' = 'Long-term trend')), 
               aes(value, pgr), color = 'blue') +
    labs(y = 'per-capita growth rate', x = NULL)
  
  return(p)
}

plot_K_predictions = function(mod, predictor, obsdat, scale = 1, center = 0) {
  # # PREDICTED CARRYING CAPACITY (in same scale/units as N)
  # # for a range of the predictor values
  # for (j in 1:npred) {
  #   dsim[j] = dd_intercept + beta[3] * Kpred[j]
  #   K[j] = -1 * rate/dsim[j]
  # }
  dat = MCMCvis::MCMCchains(mod, params = c('dd_intercept', 'kappa', 'rate'))
  xvals = (log(predictor + 1) - 14.64198)/(15.96866-12.91454)
  
  predvalues = purrr::map(xvals,
                          function(x) -dat[,'rate'] / 
                            (dat[,'dd_intercept'] + x*dat[,'kappa']) * scale) %>% 
    purrr::map_df(function(x) {
      c('value' = median(x), HDInterval::hdi(x, 0.95))
    }) %>% 
    mutate(predictor = predictor)
  
  p = ggplot(predvalues, aes(predictor/1000, value)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') +
    geom_point() + geom_line() +
    geom_point(data = obsdat, aes(flow14_cum3/1000, burrows, 
                                  color = as.factor(year))) +
    scale_color_viridis_d() +
    xlab('Cumulative flow >14cfs (3 years) (thousands)') +
    ylab('carrying capacity')

  return(p)
}