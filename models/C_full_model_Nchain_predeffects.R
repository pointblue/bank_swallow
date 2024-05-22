# mini Nchain model has one effect of flow and one effect on density dependence

# >> model inspired by Bayesian Models: A Statistical Primer for Ecologists;
# wildebeest example (sections 1.2.1 and 8.5)

model {

  ## PRIORS
  # on log scale:
  rmax ~ dunif (-2, 2) # max growth rate when N is small (on log.scale)
  beta[1] ~ dnorm(0, 1/100) # flowt
  beta[2] ~ dnorm(0, 1/100) # flowt1
  beta[3] ~ dnorm(0, 1/100) # drought1
  beta[4] ~ dnorm(0, 1/100) # trend
  beta[5] ~ dnorm(0, 1/100) # effect of density dependence (=-rmax/K)
  
  # error
  sigma.o ~ dunif(0, 5)
  tau.o = 1/sigma.o
  
  sigma.p ~ dunif(0, 5)
  tau.p = 1/sigma.p

  # estimate N1 and R1 from observed data
  N[1] ~ dlnorm(log(y[1]), tau.o) 
  N.sim[1] ~ dlnorm(log(y[1]), tau.o)
  R[1] = 1.495942
  
  # MODEL: start t at 2 because we need N[t-1] (a.k.a. 1999)
  for (t in 2:length(y[])) {
    # observation model:
    # observed current burrow count (y[t]) depends on mean count N with
    # observation error
    y[t] ~ dlnorm(log(N[t]), tau.o) # observed burrows
    
    # process model:
    N[t] ~ dlnorm(log(mean.N[t]), tau.p) #expected number of burrows
    mean.N[t] = max(.0001, R[t] * N[t-1])
    
    # growth rate: effects of density dependence, flow, and
    # annual trend in the effect on density-dependence
    log(R[t]) = rmax + 
      beta[1] * predictors[t, 1] + #flowt
      beta[2] * predictors[t, 2] + #flowt1
      beta[3] * predictors[t, 3] + #drought1
      beta[4] * predictors[t, 4] + #long-term trend
      beta[5] * N[t-1]/1000  # effect of density dependence (~ -rmax/K) 
    # (divide by 1000 to put the Beta on a more similar scale to the others)
  }
  
  # model checking for every year of modeled data:
  for (t in 2:length(y[])) {
    # simulated values for observed growth and burrow counts
    N.sim[t] ~ dlnorm(log(R[t] * N.sim[t-1]), tau.p)
    y.sim[t] ~ dlnorm(log(N[t]), tau.o)
    R.sim[t] = N.sim[t]/N.sim[t-1]
    
    # squared error
    sq[t] = (y[t] - N[t])^2 #observed and expected burrow counts
    sqsim[t] = (y.sim[t] - N[t])^2 #simulated and expected burrow counts
    
    # residuals for checking autocorrelation
    e[t] = y[t] - N[t] 
  }

  # POSTERIOR PREDICTIVE CHECKS--------
  # for every year of modeled data
  mean.y = mean(y[2:length(y[])])
  mean.ysim = mean(y.sim[2:length(y.sim[])])
  pvalue.mean = step(mean.ysim - mean.y) # step(x) returns 0 if x<0 (y > ysim) and 1 otherwise

  sd.y = sd(y[2:length(y.sim[])])
  sd.ysim = sd(y.sim[2:length(y.sim[])])
  pvalue.sd = step(sd.ysim - sd.y)

  # calculate Bayesian P value
  fit <- sum(sq[2:length(sq[])]) #sum of squares
  fit.new <- sum(sqsim[2:length(sqsim[])]) #sum of squares
  pvalue <- step(fit.new - fit) #difference: step(x) returns 0 if x<0 and 1 otherwise

  # effect of each predictor on estimated growth rate, with all other predictors
  # held at mean value/0
  for (i in 1:npred) {
    
    for (j in 1:4) { # include effect of density dependence at median pop size
      log(R.pred[i,j]) = rmax + beta[j] * predmatrix[i,j] + beta[5]*15000/1000
    }
    
    log(R.pred[i,5]) = rmax + beta[5] * predmatrix[i,5] #effect of density dependence
  }

}