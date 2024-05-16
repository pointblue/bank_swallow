# FIT FULL MODEL: dynamic population model
# >> include additional covariates for effects on growth rate
# >> include additional model checks (e.g., autocorrelation)
#
# HYPOTHESES---------
# H1. population growth (t-1 to t) is driven by river flow in the immediate
# previous winter (year t) creating more quality habitat, and driving
# immigration
#
# H2. population growth is driven by river flow in the previous winter (year
# t-1) creating more high quality habitat, resulting in higher reproductive
# success and recruitment
#
# H3. population growth is affected by drought/food availability in the previous
# breeding season, resulting in reduced reproductive success (and maybe adult
# survival?)
#
# H4. population growth is influenced by density-dependent reproductive success,
# and the population size in the previous breeding season
# 
# consider also long-term trend in growth rate

# load packages and functions
source('R/packages.R')
source('R/functions_figs.R')

modeldat = read_csv('data/modeldat.csv') |> filter(WY >= 1999)

# same as for mini model
inits = list(
  list(sigma.p = 0.05, sigma.o = 0.01),
  list(sigma.p = 0.01, sigma.o = 0.10),
  list(sigma.p = 0.10, sigma.o = 0.05))


# INPUT DATA----------
# include 1999 so that we can access the burrow count in 1999 to inform model
# from 2000 onward (n = 25); all other vectors need to be the same length
inputdat1 = list(
  # y = all observed burrow counts, 1999 through 2023; use original (unrounded)
  # data
  y = modeldat |> pull(burrows_orig), #|> scale(center = FALSE) |> as.vector(),

  # effects on R driven by current flow, prior flow, and prior drought
  predictors = modeldat |> 
    select(log.flowt, log.flowt1, drought1, WY) |>
    mutate_all(~scale(.))
)

# FIT MODEL-----------
mod = rjags::jags.model(file = 'models/B_full_model_Nchain.R',
                        data = inputdat1,
                        n.adapt = 20000,
                        n.chains = 3,
                        inits = inits)
update(mod, n.iter = 20000)
mod_res = rjags::coda.samples(
  mod, 
  variable.name = c('R', 'N', 'beta', 'rmax', 'e', 'fit', 'fit.new',
                    'sigma.o', 'sigma.p', 'y.sim', 'N.sim', 'R.sim',
                    'pvalue', 'pvalue.mean', 'pvalue.sd', 'sq', 'sqsim'), 
  n.iter = 100000, thin = 25)

# INSPECT OUTPUT------------
# key model parameters:
MCMCvis::MCMCsummary(mod_res, 
                     params = c('beta', 'rmax', 'sigma.o', 'sigma.p'), 
                     HPD = TRUE, pg0 = TRUE)
# convergence good 
# n.eff good
# beta1 (flowt) may be negative (p=0.77)
# beta2 (flowt1) positive (p=0.99)
# beta3 (drought1) neutral (p=0.55 to be positive)
# beta4 (density dependence) negative (p=0.99)
# beta5 (trend in density dependence) neutral
# rmax 0.47 (on the log scale) = 1.6

MCMCvis::MCMCtrace(mod_res,
                   params = c('beta', 'rmax', 'sigma.o', 'sigma.p'),
                   pdf = FALSE)
# trace plots look good, and rmax no longer bumping up against 1 after extending
# range for prior to -2 to 2 (from -1 to 1) and beta[1] no longer bumping up
# against lower limit

# model checking:
MCMCvis::MCMCsummary(mod_res, c('pvalue', 'pvalue.mean', 'pvalue.sd'), 
                     HPD = TRUE, pg0 = TRUE)
# mod1d: pvalue.mean and pvalue close to 0.5 (good); pvalue.sd a bit high
# (meaning sd of burrow counts tends higher for observed than simulated)
MCMCvis::MCMCplot(mod_res, params = c('sq', 'sqsim'), HPD = TRUE)
MCMCvis::MCMCplot(mod_res, params = c('pvalue', 'pvalue.mean', 'pvalue.sd'), HPD = TRUE)

# check for autocorrelation in the residuals
e = MCMCpstr(mod_res, params = "e", func = mean)
acf(e$e)
# no evidence of autocorrelation

# observed vs. predicted:
plot_observed_predicted(mod_res, obsdat = modeldat, HPD = TRUE, pg0 = TRUE)
