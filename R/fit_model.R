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
  if (type == 'growth_simple') {
    # growth rate depends only on density-independent environmental predictors
    
    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(idat$Bpredictors)[2] + 1
    
    modelstring = "
  model {

    ## PRIORS
    rate ~ dnorm(0, 1/100)
    sigma ~ dunif(0, 2)
    tau = 1/sigma^2

    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }

    ## LIKELIHOOD: 
    # y = observed growth rate from t to t+1

    for (t in 1:fityears){
    
      y[t] ~ dnorm(mu[t], tau)
      
      mu[t] = rate +
        beta[1] * Bpredictors[t, 1] +
        beta[2] * Bpredictors[t, 2] + 
        beta[3] * Bpredictors[t, 3] + 
        beta[4] * year[t]

      # simulated growth rates
      ysim[t] ~ dnorm(mu[t], tau)
      # Nsim[t+1] = N[t] * (1 + ysim[t])

      # squared error
      sq[t] = (y[t] - mu[t])^2
      sqsim[t] = (ysim[t] - mu[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (separate for fitted & predicted portions of time series)
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
    Nsim[1] = N[1]
    for (i in 2:7) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
    Nsim[8] = Nsim[7] * (1 + ysim[7])
    for (i in 9:(fityears+1)) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
    

  }"
  } else if (type == 'growth_logistic') {
    # growth rate depends on density (relative to carrying capacity) and
    # density-independent environmental predicotrs; carrying capacity does not
    # vary over time
    
    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(idat$Bpredictors)[2] + 1
    
    modelstring = "
  model {

    ## PRIORS
    K ~ dunif(0, 10) # scaled in tens of thousands
    rate ~ dnorm(0, 1/100)
    sigma ~ dunif(0, 2)
    tau = 1/sigma^2

    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }

    ## LIKELIHOOD: 
    # y = observed growth rate from t to t+1
    # N = # burrows in year t (assumed measured without error)
    # K = carrying capacity (assumed not to vary?)
    
    for (t in 1:fityears){
    
      y[t] ~ dnorm(mu[t], tau)
      
      mu[t] = rate * (1 - N[t]/K) +
        beta[1] * Bpredictors[t, 1] +
        beta[2] * Bpredictors[t, 2] + 
        beta[3] * Bpredictors[t, 3] + 
        beta[4] * year[t]
        
      # simulated growth rates
      ysim[t] ~ dnorm(mu[t], tau)
      
      # squared error
      sq[t] = (y[t] - mu[t])^2
      sqsim[t] = (ysim[t] - mu[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (separate for fitted & predicted portions of time series)
    mean.y = mean(y[])
    mean.ysim = mean(ysim[])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise

    sd.y = sd(y[])
    sd.ysim = sd(ysim[])
    pvalue.sd = step(sd.ysim - sd.y)

    # cv.y = sd(y[]) / mean.y
    # cv.ysim = sd(ysim[]) / mean.ysim
    # pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[])
    fitsim = sum(sqsim[])
    pvalue.fit = step(fitsim - fit)
    
    ## PREDICTED BURROW NUMBERS:
    Nsim[1] = N[1]
    for (i in 2:7) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
    Nsim[8] = Nsim[7] * (1 + ysim[7])
    for (i in 9:(fityears+1)) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
  }"
  } else if (type == 'growth_logistic_alt') {
    # alternative parameterization of above model, allowing strength of density
    # dependence to be estimated independent of r and K
    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(idat$Bpredictors)[2] + 1
    
    modelstring = "
  model {

    ## PRIORS
    dd ~ dnorm(0, 1/100)
    rate ~ dnorm(0, 1/100)
    sigma ~ dunif(0, 2)
    tau = 1/sigma^2

    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }

    ## LIKELIHOOD: 
    # y = observed growth rate from t to t+1 (observed with error)
    # rate = intrinsic growth rate (when all predictors at mean values)
    # dd = strength of density dependence ~ -rate/K (carrying capacity);
    #   assumed not to vary
    # N = # burrows in year t (assumed measured without error)
    
    for (t in 1:fityears){
    
      y[t] ~ dnorm(mu[t], tau)
      
      mu[t] = rate + dd * N[t] +
        beta[1] * Bpredictors[t, 1] +
        beta[2] * Bpredictors[t, 2] + 
        beta[3] * Bpredictors[t, 3] + 
        beta[4] * year[t]

      # simulated growth rates
      ysim[t] ~ dnorm(mu[t], tau)
      # Nsim[t+1] = ysim[t] * N[t] + N[t]

      # squared error
      sq[t] = (y[t] - mu[t])^2
      sqsim[t] = (ysim[t] - mu[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (separate for fitted & predicted portions of time series)
    mean.y = mean(y[])
    mean.ysim = mean(ysim[])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise

    sd.y = sd(y[])
    sd.ysim = sd(ysim[])
    pvalue.sd = step(sd.ysim - sd.y)

    # cv.y = sd(y[]) / mean.y
    # cv.ysim = sd(ysim[]) / mean.ysim
    # pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[])
    fitsim = sum(sqsim[])
    pvalue.fit = step(fitsim - fit)
    
    ## PREDICTED BURROW NUMBERS:
    Nsim[1] = N[1]
    for (i in 2:7) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
    Nsim[8] = Nsim[7] * (1 + ysim[7])
    for (i in 9:(fityears+1)) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
  }"
  } else if (type == 'growth_logistic_variabledd') {
    # growth rate depends on density (relative to carrying capacity) and
    # density-independent environmental predictors; carrying capacity also
    # varies with a long-term trend
    
    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(idat$Bpredictors)[2] + 1
    
    modelstring = "
  model {

    ## PRIORS
    dd_intercept ~ dunif(-5, 5) #in tens of thousands
    rate ~ dnorm(0, 1/100)
    sigma ~ dunif(0, 2)
    tau = 1/sigma^2
    
    mu.d ~ dnorm(0, 1/100)
    sigma.d ~ dunif(0, 2)
    tau.d = 1/sigma.d^2
    
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }

    ## LIKELIHOOD: 
    # y = observed growth rate from t to t+1
    # N = # burrows in year t (assumed measured without error)
    # k_intercept = carrying capacity in year 1 (total # of burrows possible)
    # kappa = slope of linear trend effect on carrying capacity
    
    for (t in 1:fityears){
    
      y[t] ~ dnorm(mu[t], tau)
      
      mu[t] = rate + dd[t] * N[t] +
        beta[1] * Bpredictors[t, 1] + 
        beta[2] * Bpredictors[t, 2] + 
        beta[3] * Bpredictors[t, 3] + 
        beta[4] * year[t]

      # carrying capacity (strength of density dependence) varies annually
      dd[t] ~ dnorm(mu.d, tau.d)

      # simulated growth rates
      ysim[t] ~ dnorm(mu[t], tau)

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
    Nsim[1] = N[1]
    for (i in 2:7) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
    Nsim[8] = Nsim[7] * (1 + ysim[7])
    for (i in 9:(fityears+1)) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
  }"
  } else if (type == 'growth_logistic_trenddd') {
    # growth rate depends on density (relative to carrying capacity) and
    # density-independent environmental predictors; carrying capacity also
    # varies with a long-term trend
    
    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(idat$Bpredictors)[2] + 1
    
    modelstring = "
  model {

    ## PRIORS
    dd_intercept ~ dunif(-5, 5) 
    rate ~ dnorm(0, 1/100)
    sigma ~ dunif(0, 2)
    tau = 1/sigma^2
    
    sigma.d ~ dunif(0, 2)
    tau.d = 1/sigma.d^2
    kappa ~ dnorm(0, 1/100)
    
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }

    ## LIKELIHOOD: 
    # y = observed growth rate from t to t+1
    # N = # burrows in year t (assumed measured without error)
    
    for (t in 1:fityears){
    
      y[t] ~ dnorm(mu[t], tau)
      
      mu[t] = rate + dd[t] * N[t] +
        beta[1] * Bpredictors[t, 1] + 
        beta[2] * Bpredictors[t, 2] + 
        beta[3] * Bpredictors[t, 3] + 
        beta[4] * year[t]

      # carrying capacity (strength of density dependence) varies with long-term trend
      # kappa = slope of linear trend effect on carrying capacity
      dd[t] ~ dnorm(mu.d[t], tau.d)
      mu.d[t] = dd_intercept + kappa * year[t]
      
      # simulated growth rates
      ysim[t] ~ dnorm(mu[t], tau)

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
    Nsim[1] = N[1]
    for (i in 2:7) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
    Nsim[8] = Nsim[7] * (1 + ysim[7])
    for (i in 9:(fityears+1)) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
  }"
  } else if (type == 'growth_logistic_varK') {
    # growth rate depends on density (relative to carrying capacity) and
    # density-independent environmental predictors; carrying capacity also
    # varies with predictors
    
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
    # N = # burrows in year t (assumed measured without error)
    # dd = strength of density dependence in year t (~ -rate/K), where K =
    # carrying capacity, and dd is affected by environmental factors (stream)
    
    for (t in 1:fityears){
    
      y[t] ~ dnorm(mu[t], tau.o) #error in observation
      
      mu[t] = r0 + k[t] * N[t] +
        beta[1] * Bpredictors[t, 1] + 
        beta[2] * Bpredictors[t, 2] + 
        beta[3] * Bpredictors[t, 3] + 
        error[t]
        
      # overdisperson by observation/year
      error[t] ~ dnorm(0, tau.p) #error in process model

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
    Nsim[1] = N[1] #start with observed value
    
    # 2000 through 2006 (i=2 to i=8):
    # --> include prediction of unobserved burrow numbers for 2006 from
    # observations in 2005 and estimated growth rate:
    #    Nsim[2006] = N[2005] * (1 + ysim[2005-06])
    
    for (i in 2:8) {
      Nsim[i] = N[i-1] * (1 + ysim[i-1])
    }
    
    # 2007 (i=9) use observed value to fill in
    Nsim[9] = N[8]
    
    # # for 2007, prediction of burrows depends on unknown number of
    # # burrows in 2006 and unknown growth rate from 2006 to 2007:
    # #    Nsim[2007] = N[2006] * (1 + ysim[2006-07])
    # # --> so use predicted value for burrows in 2006 instead, and estimate
    # # unobserved growth from 2006 to 2007 from predicted burrows in 2006 and
    # # observed burrows in 2007 (can't use model above because no N observation
    # # in 2006 to account for density dependence)
    # #    ysim[2006-07] = (N[2007] - N[2006]) / N[2006]
    # # Note: N[8] refers to 2007 because there is no entry in N for 2006
    # 
    # ysim.2006 = (N[8] - Nsim[8]) / Nsim[8]
    # 
    # Nsim[9] = Nsim[8] * (1 + ysim.2006) #will end up being the same as N[2007] (confounded)

    # 2008 through 2021:
    for (i in 10:(fityears + 2)) {
      Nsim[i] = N[i-2] * (1 + ysim[i-2])
    }
    
    ## GROWTH RATE VS. PREDICTORS
    # all other predictors held at mean value/0, including # of burrows, thus no
    # additional effect of density dependence
    
    for (i in 1:npred) {
      e[i] ~ dnorm(0, tau.p) #error in process model (treat as same across predictors)
      k.pred[i] ~ dnorm(mu.k.pred[i], tau.k)
      mu.k.pred[i] = k0
      
      ypred[i, 1] ~ dnorm(mu.y.pdsi[i], tau.o) #error in observation
      mu.y.aflow[i] = r0 + k.pred[i] * mean(N[]) + beta[1] * predmatrix[i,1] + e[i]
      
      
      ypred[i, 2] ~ dnorm(mu.y.pdsi[i], tau.o) #error in observation
      mu.y.pdsi[i] = r0 + k.pred[i] * mean(N[]) + beta[2] * predmatrix[i,2] + e[i]
      
      ypred[i, 3] ~ dnorm(mu.y.t[i], tau.o) #error in observation
      mu.y.t[i] = r0 + k.pred[i] * mean(N[]) + beta[3] * predmatrix[i,3] + e[i]
      
      # add effect on carrying capacity/strength of density dependence
      ypred[i, 4] ~ dnorm(mu.dd.pred[i], tau.k)
      mu.dd.pred[i] = k0 + kappa[1] * predmatrix[i,4]
    }

    
  }"
  } else if (type == 'N_simple_poisson') {
    # burrow counts are imperfect estimates of true N (with poisson error);
    # model change between true Ns directly; no density dependence
    
    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(inputdat$Bpredictors)[2] + 1
    
    modelstring = "
  model {
    ## PRIORS
    N[1] ~ dpois(y[1])
    Nsim[1] ~ dpois(y[1])
    
    rate ~ dnorm(0, 1/100)
    
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }
    
    sigma ~ dunif(0, 2)
    tau = 1/sigma^2

    ## LIKELIHOOD: y = # burrows counted
    for (t in 2:fityears){
      y[t] ~ dpois(N[t])

      log(N[t]) = log(N[t-1]) + rate +
        beta[1] * Bpredictors[t-1, 1] + 
        beta[2] * Bpredictors[t-1, 2] + 
        beta[3] * Bpredictors[t-1, 3] + 
        beta[4] * year[t] +
        error[t]

      # additional random error by observation/year (overdispersion)
      error[t] ~ dnorm(0, tau)
      
      # simulated burrow count data & per capita growth rate
      Nsim[t] ~ dpois(N[t])
      ysim[t-1] = (Nsim[t] - Nsim[t-1])/Nsim[t-1]
      
      # squared error
      sq[t] = (y[t] - N[t])^2
      sqsim[t] = (Nsim[t] - N[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: 
    mean.y = mean(y[2:fityears])
    mean.ysim = mean(Nsim[2:fityears])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    sd.y = sd(y[2:fityears])
    sd.ysim = sd(Nsim[2:fityears])
    pvalue.sd = step(sd.ysim - sd.y)
    
    cv.y = sd(y[2:fityears]) / mean.y
    cv.ysim = sd(Nsim[2:fityears]) / mean.ysim
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[2:fityears])
    fitsim = sum(sqsim[2:fityears])
    pvalue.fit = step(fitsim - fit) 

  }"
  } else if (type == 'N_logistic_poisson') {
    # burrow counts are imperfect estimates of true N (with poisson error);
    # model change between true Ns directly; constant density dependent effect
    
    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(inputdat$Bpredictors)[2] + 1
    
    modelstring = "
  model {
    ## PRIORS
    N[1] ~ dpois(y[1])
    Nsim[1] ~ dpois(y[1])
    
    rate ~ dnorm(0, 1/100)
    dd ~ dnorm(0, 1/100)
    
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }
    
    sigma ~ dunif(0, 2)
    tau = 1/sigma^2

    ## LIKELIHOOD: y = # burrows counted
    for (t in 2:fityears){
      y[t] ~ dpois(N[t])

      log(N[t]) = log(N[t-1]) + rate + dd * (N[t-1]) +
        beta[1] * Bpredictors[t-1, 1] + 
        beta[2] * Bpredictors[t-1, 2] + 
        beta[3] * Bpredictors[t-1, 3] + 
        beta[4] * year[t] +
        error[t]

      # additional random error by observation/year (overdispersion)
      error[t] ~ dnorm(0, tau)
      
      # simulated burrow count data & per-capita growth rate
      Nsim[t] ~ dpois(N[t])
      ysim[t-1] = (Nsim[t] - Nsim[t-1])/Nsim[t-1]
      
      # squared error
      sq[t] = (y[t] - N[t])^2
      sqsim[t] = (Nsim[t] - N[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: 
    mean.y = mean(y[2:fityears])
    mean.ysim = mean(Nsim[2:fityears])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    sd.y = sd(y[2:fityears])
    sd.ysim = sd(Nsim[2:fityears])
    pvalue.sd = step(sd.ysim - sd.y)
    
    cv.y = sd(y[2:fityears]) / mean.y
    cv.ysim = sd(Nsim[2:fityears]) / mean.ysim
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[2:fityears])
    fitsim = sum(sqsim[2:fityears])
    pvalue.fit = step(fitsim - fit) 

  }"
  } else if (type == 'N_simple_negbin') {
    # no density dependence
    
    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(inputdat$Bpredictors)[2] + 1
    
    modelstring = "
  model {

    ## PRIORS
    N[1] ~ dnegbin(k / (k + y[1]), k) 
    Nsim[1] ~ dnegbin(k / (k + y[1]), k)
    k ~ dgamma(0.001,0.001) #dispersion parameter for negative binomial

    rate ~ dnorm(0, 1/100)

    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }
    
    sigma ~ dunif(0, 2)
    tau = 1/sigma^2
    
    # var.error ~ dgamma(0.001, 0.001)
    # tau.error = 1/var.error
    
    # var.error ~ dunif(0, 10)
    # tau.error = 1/var.error
    
    ## LIKELIHOOD: y = # burrows counted
    for (t in 2:fityears){
      y[t] ~ dnegbin(p[t], k) 
      p[t] = k / (k + N[t])

      log(N[t]) = log(N[t-1]) + rate +
        beta[1] * Bpredictors[t-1, 1] + 
        beta[2] * Bpredictors[t-1, 2] + 
        beta[3] * Bpredictors[t-1, 3] + 
        beta[4] * year[t] +
        error[t]
        
      # overdisperson by observation/year
      error[t] ~ dnorm(0, tau)

      # simulated burrow count data & per-capita growth
      Nsim[t] ~ dnegbin(p[t], k) #if negative binomial
      ysim[t-1] = (Nsim[t] - Nsim[t-1])/Nsim[t-1]
      
      # squared error
      sq[t] = (y[t] - N[t])^2
      sqsim[t] = (Nsim[t] - N[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (separate for fitted & predicted portions of time series)
    mean.y = mean(y[2:fityears])
    mean.ysim = mean(Nsim[2:fityears])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    sd.y = sd(y[2:fityears])
    sd.ysim = sd(Nsim[2:fityears])
    pvalue.sd = step(sd.ysim - sd.y)
    
    cv.y = sd(y[2:fityears]) / mean.y
    cv.ysim = sd(Nsim[2:fityears]) / mean.ysim
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[2:fityears])
    fitsim = sum(sqsim[2:fityears])
    pvalue.fit = step(fitsim - fit) 

  }"
  } else if (type == 'N_logistic_negbin') {
    # no density dependence
    
    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(inputdat$Bpredictors)[2] + 1
    
    modelstring = "
  model {

    ## PRIORS
    N[1] ~ dnegbin(k / (k + y[1]), k) 
    Nsim[1] ~ dnegbin(k / (k + y[1]), k)
    k ~ dgamma(0.001,0.001) #dispersion parameter for negative binomial

    rate ~ dnorm(0, 1/100)
    dd ~ dnorm(0, 1/100)
    
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }
    
    sigma ~ dunif(0, 2)
    tau = 1/sigma
    
    # var.error ~ dgamma(0.001, 0.001)
    # tau.error = 1/var.error
    
    # var.error ~ dunif(0, 10)
    # tau.error = 1/var.error
    
    ## LIKELIHOOD: y = # burrows counted
    for (t in 2:fityears){
      y[t] ~ dnegbin(p[t], k) 
      p[t] = k / (k + N[t])

      log(N[t]) = log(N[t-1]) + rate + dd * N[t-1] +
        beta[1] * Bpredictors[t-1, 1] + 
        beta[2] * Bpredictors[t-1, 2] + 
        beta[3] * Bpredictors[t-1, 3] + 
        beta[4] * year[t] +
        error[t]
        
      # overdisperson by observation/year
      error[t] ~ dnorm(0, tau)

      # simulated burrow count data & per-capita growth
      Nsim[t] ~ dnegbin(p[t], k) #if negative binomial
      ysim[t-1] = (Nsim[t] - Nsim[t-1])/Nsim[t-1]
      
      # squared error
      sq[t] = (y[t] - N[t])^2
      sqsim[t] = (Nsim[t] - N[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (separate for fitted & predicted portions of time series)
    mean.y = mean(y[2:fityears])
    mean.ysim = mean(Nsim[2:fityears])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    sd.y = sd(y[2:fityears])
    sd.ysim = sd(Nsim[2:fityears])
    pvalue.sd = step(sd.ysim - sd.y)
    
    cv.y = sd(y[2:fityears]) / mean.y
    cv.ysim = sd(Nsim[2:fityears]) / mean.ysim
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[2:fityears])
    fitsim = sum(sqsim[2:fityears])
    pvalue.fit = step(fitsim - fit) 

  }"
  } else if (type == 'N_logistic_negbin_varK') {

    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(idat$Bpredictors)[2] + 1
    idat$nKcov = dim(idat$Kpredictors)[2]
    
    modelstring = "
  model {

    ## PRIORS
    N[1] ~ dnegbin(k / (k + y[1]), k) 
    Nsim[1] ~ dnegbin(k / (k + y[1]), k)
    k ~ dgamma(0.001,0.001) #dispersion parameter for negative binomial

    rate ~ dnorm(0, 1/100)
    sigma ~ dunif(0, 10)
    tau = 1/sigma
    # tau ~ dgamma(0.01, 0.01)
    
    dd_intercept ~ dunif(-5, 5) #in tens of thousands
    sigma.d ~ dunif(0, 10)
    tau.d = 1/sigma.d
    # tau.d ~ dgamma(0.01, 0.01)
    
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }
    for (i in 1:nKcov) {
      kappa[i] ~ dnorm(0, 1/100)
    }

    ## LIKELIHOOD: y = # burrows counted
    for (t in 2:fityears){
      y[t] ~ dnegbin(p[t], k) 
      p[t] = k / (k + N[t])

      log(N[t]) = log(N[t-1]) + rate + dd[t-1] * N[t-1] +
        beta[1] * Bpredictors[t-1, 1] + 
        beta[2] * Bpredictors[t-1, 2] + 
        beta[3] * year[t] +
        error[t]
        
      # overdisperson by observation/year
      error[t] ~ dnorm(0, tau)
      
      # dd affected by predictors
      dd[t-1] ~ dnorm(mu.d[t-1], tau.d)
      mu.d[t-1] = dd_intercept + kappa[1] * Kpredictors[t-1, 1]

      # simulated burrow count data & per-capita growth
      Nsim[t] ~ dnegbin(p[t], k) #if negative binomial
      ysim[t-1] = (Nsim[t] - Nsim[t-1])/Nsim[t-1]
      
      # squared error
      sq[t] = (y[t] - N[t])^2
      sqsim[t] = (Nsim[t] - N[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (separate for fitted & predicted portions of time series)
    mean.y = mean(y[2:fityears])
    mean.ysim = mean(Nsim[2:fityears])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    sd.y = sd(y[2:fityears])
    sd.ysim = sd(Nsim[2:fityears])
    pvalue.sd = step(sd.ysim - sd.y)
    
    cv.y = sd(y[2:fityears]) / mean.y
    cv.ysim = sd(Nsim[2:fityears]) / mean.ysim
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[2:fityears])
    fitsim = sum(sqsim[2:fityears])
    pvalue.fit = step(fitsim - fit) 

  }"
  } else if (type == 'simple_negbin') {
    # no density dependence
    
    idat = inputdata
    idat$fityears = length(idat$year)
    idat$nBcov = dim(inputdat$Bpredictors)[2] + 1
    
    modelstring = "
  model {

    ## PRIORS
    N[1] ~ dnegbin(k / (k + y[1]), k)
    k ~ dgamma(0.001, 0.001) #dispersion parameter for negative binomial
    
    rate ~ dnorm(0, 1/100)
    
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }
    
    ## LIKELIHOOD: y = # burrows counted
    for (t in 2:fityears){
      y[t] ~ dnegbin(p[t], k) 
      p[t] = k / (k + N[t])

      log(N[t]) = log(N[t-1]) + rate +
        beta[1] * predictors[t, 1] + 
        beta[2] * predictors[t, 2] + 
        beta[3] * predictors[t, 3] + 
        beta[4] * predictors[t, 4] + 
        beta[5] * predictors[t, 5] +
        beta[6] * predictors[t, 6] +
        beta[7] * year[t]

      # simulated burrow count data
      ysim[t] ~ dnegbin(p[t], k)
      
      # squared error
      sq[t] = (y[t] - N[t])^2
      sqsim[t] = (ysim[t] - N[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (separate for fitted & predicted portions of time series)
    mean.y = mean(y[2:fityears])
    mean.ysim = mean(ysim[2:fityears])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    sd.y = sd(y[2:fityears])
    sd.ysim = sd(ysim[2:fityears])
    pvalue.sd = step(sd.ysim - sd.y)
    
    cv.y = sd(y[2:fityears]) / mean.y
    cv.ysim = sd(ysim[2:fityears]) / mean.ysim
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[2:fityears])
    fitsim = sum(sqsim[2:fityears])
    pvalue.fit = step(fitsim - fit) 

  }"
  } else if (type == 'thetalog_negbin_varK') {
    idat = inputdata
    idat$fityears = length(idat$year)
    
    modelstring = "
  model {

    ## PRIORS
    N[1] ~ dnegbin(k / (k + y[1]), k)
    k ~ dgamma(0.01,0.01)
    
    theta ~ dunif(0, 5)
    Kappa ~ dunif(10000, 50000)
    rate ~ dunif(0, 1)

    for (i in 1:5) {
      beta[i] ~ dnorm(0, 1/100)
    }
    
    ## LIKELIHOOD: y = # burrows counted
    for (t in 2:fityears){
      y[t] ~ dnegbin(p[t], k) 
      p[t] = k / (k + N[t])

      # theta-logistic model of population growth for N ('true' number of burrows)
      #  rate = intrinsic growth rate
      #  K[t] = year-specific carrying capacity
      #  theta = form/strength of density dependence
      #  beta = effects of other density-independent predictors on K or rate;
      #    beta1 affects K, others affect rate
      #  Kappa = average/baseline carrying capacity
      
      log(N[t]) = log(N[t-1]) + rate * (1 - (N[t-1] / K[t])^theta) + 
        beta[2] * predictors[t, 2] + 
        beta[3] * predictors[t, 3] + 
        beta[4] * predictors[t, 4] + 
        beta[5] * year[t]

      K[t] = Kappa + beta[1] * predictors[t, 1]
      
      # simulated burrow count data
      ysim[t] ~ dnegbin(p[t], k) 
      
      # squared error
      sq[t] = (y[t] - N[t])^2
      sqsim[t] = (ysim[t] - N[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS: (separate for fitted & predicted portions of time series)
    mean.y = mean(y[2:fityears])
    mean.ysim = mean(ysim[2:fityears])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    sd.y = sd(y[2:fityears])
    sd.ysim = sd(ysim[2:fityears])
    pvalue.sd = step(sd.ysim - sd.y)
    
    cv.y = sd(y[2:fityears]) / mean.y
    cv.ysim = sd(ysim[2:fityears]) / mean.ysim
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[2:fityears])
    fitsim = sum(sqsim[2:fityears])
    pvalue.fit = step(fitsim - fit) 

  }"
    
  } else if (type == 'thetalog_negbin') {
    idat = inputdata
    idat$fityears = length(idat$year)
    
    modelstring = "
  model {

    ## PRIORS
    N[1] ~ dnegbin(k / (k + y[1]), k)
    k ~ dgamma(0.01,0.01)
    
    theta ~ dunif(0, 5)
    K ~ dunif(10000, 50000)
    
    rate ~ dnorm(0, 1/100)
    error ~ dnorm(0, 1/100)
    
    for (i in 1:5) {
      beta[i] ~ dnorm(0, 1/100)
    }
    
    ## LIKELIHOOD: y = # burrows counted
    for (t in 2:fityears){
      y[t] ~ dnegbin(p[t], k)
      p[t] = k / (k + N[t])

      log(N[t]) = log(N[t-1]) + rate * (1 - (N[t-1] / K)^theta) + 
        beta[1] * predictors[t, 1] + 
        beta[2] * predictors[t, 2] +
        beta[3] * predictors[t, 3] + 
        beta[4] * predictors[t, 4] +
        beta[5] * year[t] + 
        error

      # simulated burrow count data
      ysim[t] ~ dnegbin(p[t], k)
      
      # squared error
      sq[t] = (y[t] - N[t])^2
      sqsim[t] = (ysim[t] - N[t])^2
    }

    ## POSTERIOR PREDICTIVE CHECKS:
    mean.y = mean(y[2:fityears])
    mean.ysim = mean(ysim[2:fityears])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    sd.y = sd(y[2:fityears])
    sd.ysim = sd(ysim[2:fityears])
    pvalue.sd = step(sd.ysim - sd.y)
    
    cv.y = sd(y[2:fityears]) / mean.y
    cv.ysim = sd(ysim[2:fityears]) / mean.ysim
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[2:fityears])
    fitsim = sum(sqsim[2:fityears])
    pvalue.fit = step(fitsim - fit) 

  }"
  } else if (type == 'ricker_nb') {
    if (holdout > 0) {
    idat = inputdata
    idat$fityears = length(idat$year) - holdout #number of years to include in model fitting
    idat$startpred = idat$fityears + 1 #start year for predictions
    idat$endpred = length(idat$year) #end year for predictions
    idat$nBcov = dim(idat$Bpredictors)[2] + 1 # +1 for trend
    
    modelstring = "
  model {

    ## PRIORS
    for (i in 1:nBcov) {
      beta[i] ~ dnorm(0, 1/100)
    }
    #error ~ dnorm(0, 1/100)
    r ~ dgamma(0.01,0.01) #for negative binomial
    # N[1] ~ dpois(y[1]) #base N1 on initial observation
    N[1] ~ dnegbin(r / (r + y[1]), r)
    
    ## LIKELIHOOD: y = # burrows counted
    for (t in 2:fityears){
      # assume y[t] ~ N[t] (true # of burrows), with some error:
      # y[t] ~ dpois(N[t])
      # try negative binomial, since not enough noise with poisson
      # --> based on: http://doingbayesiandataanalysis.blogspot.com/2012/04/negative-binomial-reparameterization.html
      y[t] ~ dnegbin(p[t], r)
      p[t] = r / (r + N[t])
      # v[t] = r * (1 - p[t]) / (p[t] * p[t])
      
      # true number of burrows varies with population growth rate, density
      # dependence, and density independent environmental factors:
      # beta1 = intrinsic growth rate at mean values for all predictors
      # beta2 = strength of density dependence
      # beta3-8 = change in growth rate attributable to environmental factors
      # beta9 = long-term trend in growth rate not explained by anything else
      # error = extra, unmodeled variation in growth rate
      log(N[t]) = log(N[t-1]) + beta[1] + beta[2] * N[t-1] + 
        beta[3] * predictors[t, 1] + beta[4] * predictors[t, 2] +
        beta[5] * predictors[t, 3] + beta[6] * predictors[t, 4] +
        beta[7] * predictors[t, 5] + beta[8] * predictors[t, 6] +
        beta[9] * year[t]
        #+ error

      # simulated burrow count data
      ysim[t] ~ dnegbin(p[t], r)
      
      # squared error
      sq[t] = (y[t] - N[t])^2
      sqsim[t] = (ysim[t] - N[t])^2
    }

    ## PREDICTIONS
    for (t in startpred:endpred) {
        log(N[t]) = log(N[t-1]) + beta[1] + beta[2] * N[t-1] + 
          beta[3] * predictors[t, 1] + beta[4] * predictors[t, 2] +
          beta[5] * predictors[t, 3] + beta[6] * predictors[t, 4] +
          beta[7] * predictors[t, 5] + beta[8] * predictors[t, 6] +
          beta[9] * year[t]
        #+ error
        
        p[t] = r / (r + N[t])
        ypred[t] ~ dnegbin(p[t], r)
        
        # squared error
        sq[t] = (y[t] - N[t])^2
        sqpred[t] = (ypred[t] - N[t])^2
      }
    
    ## POSTERIOR PREDICTIVE CHECKS: (separate for fitted & predicted portions of time series)
    mean.y = mean(y[2:fityears])
    mean.ysim = mean(ysim[2:fityears])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    sd.y = sd(y[2:fityears])
    sd.ysim = sd(ysim[2:fityears])
    pvalue.sd = step(sd.ysim - sd.y)
    
    cv.y = sd(y[2:fityears]) / mean.y
    cv.ysim = sd(ysim[2:fityears]) / mean.ysim
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[2:fityears])
    fitsim = sum(sqsim[2:fityears])
    pvalue.fit = step(fitsim - fit) 
      
    # for holdout:
    hmean.y = mean(y[startpred:endpred])
    hmean.ypred = mean(ypred[startpred:endpred])
    pvalue.hmean = step(hmean.ypred - hmean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    hsd.y = sd(y[startpred:endpred])
    hsd.ypred = sd(ypred[startpred:endpred])
    pvalue.hsd = step(hsd.ypred - hsd.y)
     
    hcv.y = sd(y[startpred:endpred]) / hmean.y
    hcv.ypred = sd(ypred[startpred:endpred]) / hmean.ypred
    pvalue.hcv = step(hcv.ypred - hcv.y)

    hfit = sum(sq[startpred:endpred])
    hfitpred = sum(sqpred[startpred:endpred])
    pvalue.hfit = step(hfitpred - hfit) 

  }"
  } else {
    idat = inputdata
    idat$fityears = length(idat$year)
    
    modelstring = "
  model {

    ## PRIORS
    for (i in 1:9) {
      beta[i] ~ dnorm(0, 1/100)
    }
    #error ~ dnorm(0, 1/100)
    r ~ dgamma(0.01,0.01) #for negative binomial
    # N[1] ~ dpois(y[1]) #base N1 on initial observation
    N[1] ~ dnegbin(r / (r + y[1]), r)
    
    ## LIKELIHOOD: y = # burrows counted
    for (t in 2:fityears){
      
      # # if we had more data about variation in occupancy:
      # y[t] ~ dpois(lambda[t]) # burrow count in year t is proportional to 'true' burrow count
      # theta[t] ~ dbeta(5, 5) # proportion of burrows occupied in year t; pull toward 50%
      # lambda[t] = N[t] / theta[t] #true burrow count is a function of pop size and proportion occupied
      
      # for now, assume y[t] ~ N[t], with some error:
      # y[t] ~ dpois(N[t])
      # try negative binomial, since not enough noise with poisson
      # --> based on: http://doingbayesiandataanalysis.blogspot.com/2012/04/negative-binomial-reparameterization.html
      y[t] ~ dnegbin(p[t], r)
      p[t] = r / (r + N[t])
      # v[t] = r * (1 - p[t]) / (p[t] * p[t])
      
      # beta1 = intrinsic growth rate at mean values for all predictors
      # beta2 = strength of density dependence
      # beta3-8 = change in growth rate attributable to environmental factors
      # beta9 = long-term trend in growth rate not explained by anything else
      # error = extra, unmodeled variation in growth rate
      log(N[t]) = log(N[t-1]) + beta[1] + beta[2] * N[t-1] + 
        beta[3] * predictors[t, 1] + beta[4] * predictors[t, 2] +
        beta[5] * predictors[t, 3] + beta[6] * predictors[t, 4] +
        beta[7] * predictors[t, 5] + beta[8] * predictors[t, 6] +
        beta[9] * year[t]
        #+ error

      # simulated burrow count data
      ysim[t] ~ dnegbin(p[t], r)
      
      # squared error
      sq[t] = (y[t] - N[t])^2
      sqsim[t] = (ysim[t] - N[t])^2
    }
    
    ## POSTERIOR PREDICTIVE CHECKS: (separate for fitted & predicted portions of time series)
    mean.y = mean(y[2:fityears])
    mean.ysim = mean(ysim[2:fityears])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise
    
    sd.y = sd(y[2:fityears])
    sd.ysim = sd(ysim[2:fityears])
    pvalue.sd = step(sd.ysim - sd.y)
    
    cv.y = sd(y[2:fityears]) / mean.y
    cv.ysim = sd(ysim[2:fityears]) / mean.ysim
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[2:fityears])
    fitsim = sum(sqsim[2:fityears])
    pvalue.fit = step(fitsim - fit) 

  }"
  }
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


fit_model <- function(inputdata, n.adapt = 1000, n.burnin = 5000,
                       n.sample = 10000, n.chains = 3, type = 'growth', ...) {
  vars = c('slope', 'intercept', 'beta', 'pvalue.mean', 'pvalue.fit', 'pvalue.sd')

  if (type == 'growth') {

    modelstring = "
  model {

    ## LIKELIHOOD: y = observed rate of increase from t-1 to t
    for (i in 1:length(y)){
      y[i] ~ dnorm(mu[i], tau)

      mu[i] = intercept + slope[1] * year[i] + slope[2] * year[i]^2 +
        beta[1] * predictors[i,1] + beta[2] * predictors[i,2] + 
        beta[3] * predictors[i,3] + beta[4] * predictors[i,4] + 
        beta[5] * predictors[i,5] + beta[6] * predictors[i,6] +
        beta[7] * predictors[i,7]

      # simulated data with same mean and error:
      ysim[i] ~ dnorm(mu[i], tau)

      sq[i] = (y[i] - mu[i])^2
      sqsim[i] = (ysim[i] - mu[i])^2
    }

    ## PRIORS
    for (i in 1:2) {
      slope[i] ~ dnorm(0, 1/10000)
    }
    intercept ~ dnorm(0, 1/10000)
    
    for (i in 1:npred) {
      beta[i] ~ dnorm(0, 1/10000)
    }
    sigma ~ dunif(0, 10)
    tau = 1/sigma^2

    ## POSTERIOR PREDICTIVE CHECKS:
    mean.y = mean(y[])
    mean.ysim = mean(ysim[])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 and 1 otherwise
    
    sd.y = sd(y[])
    sd.ysim = sd(ysim[])
    pvalue.sd = step(sd.ysim - sd.y)
    
    # cv.y = sd(y[]) / mean.y
    # cv.ysim = sd(ysim[]) / mean.ysim
    # pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[])
    fitsim = sum(sqsim[])
    pvalue.fit = step(fitsim - fit)

  }"
} else if (type == 'count') {
  vars = c('beta', 'pvalue.mean', 'pvalue.cv', 'pvalue.fit', 'ypred', 'sigma.overdispersion')

  modelstring = "
  model {

    ## LIKELIHOOD: y = count in year t
    for (i in 1:length(y)){
      y[i] ~ dpois(lambda[i])

      log(lambda[i]) = beta[1] + beta[2] * year[i] + overdispersion[i]

      overdispersion[i] ~ dnorm(0, tau.overdispersion)
      
      # simulated data with same mean and error:
      ysim[i] ~ dpois(lambda[i])

      sq[i] = (y[i] - lambda[i])^2
      sqsim[i] = (ysim[i] - lambda[i])^2
    }

    ## PRIORS
    beta[1] ~ dnorm(10, 1)
    beta[2] ~ dnorm(0, 1/10000)

    sigma.overdispersion ~ dunif(0, 0.5)
    tau.overdispersion = 1/sigma.overdispersion^2
    
    ## POSTERIOR PREDICTIVE CHECKS:
    mean.y = mean(y[])
    mean.ysim = mean(ysim[])
    pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 and 1 otherwise

    cv.y = sd(y[]) / mean(y[])
    cv.ysim = sd(ysim[]) / mean(ysim[])
    pvalue.cv = step(cv.ysim - cv.y)

    fit = sum(sq[])
    fitsim = sum(sqsim[])
    pvalue.fit = step(fitsim - fit)
    
    ## PREDICTED VALUES: (without overdispersion)
    for (t in 1:length(year)) {
      log(ypred[t]) <- beta[1] + beta[2] * year[t]
    }

  }"
  }

  # results = runjags::run.jags(model = modelstring,
  #                             monitor = vars,
  #                             data = inputdata[-which(names(inputdata) %in% c('year.pred', 'dat'))],
  #                             n.chains = n.chains,
  #                             adapt = n.adapt,
  #                             burnin = n.burnin,
  #                             sample = n.sample,
  #                             modules = c('glm'),
  #                             method = method,
  #                             ...)

  load.module('glm')

  jm = rjags::jags.model(file = textConnection(modelstring),
                         data = inputdata,
                         n.adapt = n.adapt,
                         n.chains = n.chains,
                         ...)
  update(jm, n.iter = n.burnin)
  results = rjags::coda.samples(jm, variable.names = vars, n.iter = n.sample,
                                thin = 1)

  # MCMCvis::MCMCsummary(results, 
  #                      params = vars[-which(grepl('ypred', vars))]) %>% 
  #   print()
  # MCMCvis::MCMCtrace(results, 
  #                    params = vars[-which(grepl('pvalue|ypred', vars))], 
  #                    pdf = FALSE)
  return(results)
}

plot_model_predictions = function(mod, year, obsdat, scale, center) {
  # pull other stats:
  predvalues = full_join(MCMCvis::MCMCpstr(mod, params = c('ysim', 'ysim.2006', 'Nsim'), 
                                           func = function(x) mean(x)) %>%
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

  p = full_join(predvalues %>% select(year, name, value:upper), 
                obsdat %>% select(year, pgr, burrows) %>% 
                  pivot_longer(-year, names_to = 'metric', 
                               values_to = 'observed') %>% 
                  mutate(name = recode(metric, 'pgr' = 'per-capita growth rate',
                                       'burrows' = 'number of burrows')) %>% 
                  select(year, name, observed)) %>% 
    ggplot(aes(year, value)) + 
    facet_wrap(~name, scales = 'free_y', nrow = 2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') +
    geom_line() + geom_point() +
    geom_line(aes(y = observed), color = 'blue') + 
    geom_point(aes(y = observed), color = 'blue')
    xlab(NULL) 

  return(p)
}

plot_cov_relationships = function(mod, xvals, nm, r0,
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
    left_join(xvals, by = c('index', 'cov' = 'name'))

  p = ggplot(predvalues, aes(value, median)) + 
    facet_wrap(~cov, scales = 'free', nrow = 2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') +
    geom_line() + geom_point() + 
    geom_hline(aes(yintercept = r0), color = 'red') +
    geom_point(data = obsdat, aes(value, pgr), color = 'blue')
  
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