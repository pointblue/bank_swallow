model {

  ## PRIORS
  # max growth rate when N is small (on log.scale)
  rmax ~ dunif (-1, 1)

  # effect sizes
  beta[1] ~ dnorm(0, 1/100) # effect of prior year flow
  beta[2] ~ dnorm(0, 1/100) # effect of density dependence

  # observation error in estimates of growth rate:
  sigma.o ~ dunif(0, 0.5)
  tau.o <- 1/sigma.o^2
  
  # annual process error in growth rate (R)
  sigma.r ~ dunif(0, 0.5)
  tau.r <- 1/sigma.r^2

  # possible "true" values of N1 (number of prior burrows, t-1)
  for (i in 1:length(y1[])) {
    N1[i] ~ dunif(10000, 25000)
  }

  for (t in 1:length(r[])) {
    # observed prior burrow count (y1[t]) depends on true prior burrow count
    # (N1[t]) -- required to be modeled because missing data
    y1[t] ~ dpois(N1[t])
    
    # observed annual growth rate from t-1 to t (r[t]) depends on "true" annual
    # growth rate (R[t])
    r[t] ~ dnorm(R[t], tau.o)
    
    # true population growth rate: effects of density dependence, flow, and
    # annual trend in the effect on density-dependence
    log(R[t]) = rmax + 
      beta[1] * flowt1[t] + #direct effect of flowt1 on population growth
      beta[2] * (N1[t]) + 
      er[t]
    er[t] ~ dnorm(0, tau.r)

    # simulated
    r.sim[t] ~ dnorm(R[t], tau.r)
    y.sim[t] ~ dpois(N1[t])
  }
}