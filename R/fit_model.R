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
                             n.sample = 10000, n.chains = 3, holdout = 0, 
                             vars = c('beta'), thin = 1,
                             ...) {
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
      #e.pred[i] ~ dnorm(0, tau.p) #error in process model (treat as same across predictors)

      ypred[i, 1] ~ dnorm(mu.ypred[i, 1], tau.o) #error in observation
      mu.ypred[i, 1] = r0 + beta[1] * predmatrix[i, 1]
      # + e.pred[i]
      
      ypred[i, 2] ~ dnorm(mu.ypred[i, 2], tau.o) #error in observation
      mu.ypred[i, 2] = r0 + beta[2] * predmatrix[i, 2]
      # + e.pred[i]

      ypred[i, 3] ~ dnorm(mu.ypred[i, 3], tau.o) #error in observation
      mu.ypred[i, 3] = r0 + beta[3] * predmatrix[i, 3]
      #+ e.pred[i]

      ypred[i, 4] ~ dnorm(mu.ypred[i, 4], tau.o) #error in observation
      mu.ypred[i, 4] = r0 + beta[4] * predmatrix[i, 4]
      # + e.pred[i]
      
      # effect of density dependence for diff pop sizes:
      ypred[i, 5] ~ dnorm(mu.ypred[i, 5], tau.o) #error in observation
      mu.ypred[i, 5] = r0 + k0 * predmatrix[i, 5] 
      #+ e.pred[i]
      # k.pred[i] ~ dnorm(k0, tau.k) #variance in k
    }
  }"
  
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