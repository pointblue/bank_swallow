model {

  ## PRIORS
  # intrinsic growth rate (positive by definition)
  r ~ dunif (0, 2)

  # carrying capacity when all vars at mean (positive by definition); max
  # observed value 19170
  k0 ~ dunif(0, 40000)
  
  # trend in carrying capacity:
  k1 ~ dnorm(0, 1/100)

  # # other potential influences on growth rate:  
  # for (i in 1:nBcov) {
  #   beta[i] ~ dnorm(0, 1/100)
  # }

  # observation error:
  sigma.o ~ dunif(0, 2)
  tau.o <- 1/sigma^2

  # # process error:
  # sigma.p ~ dunif(0, 2)
  # tau.p <- 1/sigma^2

  ## LIKELIHOOD: 
  # y = observed growth rate from t-1 to t (with observation error)
  #     y[t] = (N[t] - N[t-1])/N[t-1]
  #     y[t] = R(1 - N[t-1]/K)
  #     y[t] = r0 - r0/K * N[t-1]
  # r0 = intrinsic growth rate (when predictors are at zero/mean values)
  # N = # burrows in year t-1 (assumed measured without error)
  # k0 = -r0/K = strength of density dependence in year t-1; assumed constant
  
  for (t in 1:length(N[])) {
    
    # true growth rate
    mu[t] = r - (r/K[t]) * N[t] 
    # +
    #   beta[1] * Bpredictors[t, 1] + 
    #   beta[2] * Bpredictors[t, 2] + 
    #   beta[3] * Bpredictors[t, 3] + 
    #   beta[4] * Bpredictors[t, 4] + 
    #   beta[5] * Bpredictors[t, 5] +
    #   e[t]

    # # overdisperson by observation/year
    # e[t] ~ dnorm(0, tau.p) #error in process model

    # annual trend in K (carrying capacity)
    K[t] = k0 + k1 * Kpredictors[t, 1]
    
    # observed growth rate (with some observation error presumed)
    y[t] ~ dnorm(mu[t], tau.o)
    
    # SIMULATED GROWTH RATES
    # for the same 21 years as in inputdat$y, including for years where y[t]
    # wasn't observed (2005-06, 2020-21), but not where N wasn't observed
    # (2006-07)
    ysim[t] ~ dnorm(mu[t], tau.o)

    # squared error (in growth rates)
    sq[t] = (y[t] - mu[t])^2
    sqsim[t] = (ysim[t] - mu[t])^2
  }
}