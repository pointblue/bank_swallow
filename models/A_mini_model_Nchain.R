# mini Nchain model has one effect of flow and one effect of density dependence

# >> model inspired by Bayesian Models: A Statistical Primer for Ecologists;
# wildebeest example (sections 1.2.1 and 8.5)

model {

  ## PRIORS
  # max growth rate when N is small (on log.scale)
  rmax ~ dunif (-1, 1)

  # effect sizes
  beta[1] ~ dnorm(0, 1/100) # effect of prior year flow
  beta[2] ~ dnorm(0, 1/100) # effect of density dependence

  # error
  sigma.o ~ dunif(0, 5)
  tau.o = 1/sigma.o
  
  sigma.p ~ dunif(0, 5)
  tau.p = 1/sigma.p

  # estimate N1, N.sim1 and R1 from observed data
  N[1] ~ dlnorm(log(y[1]), tau.o) 
  N.sim[1] ~ dlnorm(log(y[1]), tau.o)
  R[1] = 1.495942
  
  # MODEL: start t at 2 because we need N[t-1] (a.k.a. 1999)
  for (t in 2:length(y[])) {
    # observation model:
    # observed current burrow count (y[t]) depends on mean count N with
    # observation error
    y[t] ~ dlnorm(log(N[t]), tau.o)
    
    # process model:
    # true current burrow count (N[t]) depends on true prior burrow count
    # (N[t-1]) and annual growth rate agr(t) from t-1 to t, with process error
    N[t] ~ dlnorm(log(mean.N[t]), tau.p)
    mean.N[t] = max(.0001, R[t] * N[t-1])
    
    # growth rate: effects of density dependence, flow, and
    # annual trend in the effect on density-dependence
    log(R[t]) = rmax + 
      beta[1] * flowt1[t] + #direct effect of flowt1 on population growth
      beta[2] * N[t-1] #effect of density dependence

  }
  
  # simulated values for observed growth and burrow counts
  for (t in 2:length(y[])) {
    N.sim[t] ~ dlnorm(log(R[t] * N.sim[t-1]), tau.p)
    y.sim[t] ~ dlnorm(log(N[t]), tau.o)
    R.sim[t] = N.sim[t]/N.sim[t-1]
  }
}