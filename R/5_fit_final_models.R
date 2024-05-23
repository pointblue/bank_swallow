# FINAL MODELS: dynamic population model with effects on annual growth rate
# >> add structure to estimate partial effects plots

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


# load packages and functions
source('R/packages.R')
source('R/functions_figs.R')

modeldat = read_csv('data/modeldat.csv') |> filter(WY >= 1999)

# same as for mini model
inits = list(
  list(sigma.p = 0.05, sigma.o = 0.01),
  list(sigma.p = 0.01, sigma.o = 0.10),
  list(sigma.p = 0.10, sigma.o = 0.05))


# INPUT DATA-----
# include 1999 so that we can access the burrow count in 1999 to
# inform model from 2000 onward (n = 25); all other vectors need to be the same
# length
inputdat_final = list(
  # y = all observed burrow counts, 1999 through 2023; use original (unrounded)
  # data
  y = modeldat |> pull(burrows_orig), 
  
  # effects on R driven by current flow, prior flow, prior drought, density,
  # long-term trend
  predictors = modeldat |> 
    select(log.flowt, log.flowt1, drought1, WY) |>
    mutate_all(~scale(.)),
  
  # range of scaled values for each predictor
  predmatrix = tibble(log.flowt = seq(-1.68, 1.99, length.out = 11),
                      log.flowt1 = seq(-1.60, 1.73, length.out = 11),
                      drought1 = seq(-1.70, 1.65, length.out = 11),
                      WY = seq(-1.63, 1.63, length.out = 11),
                      N = seq(10000, 21000, length.out = 11)/1000),
  npred = 11
)

# FIT MODEL-----
mod = rjags::jags.model(file = 'models/C_full_model_Nchain_predeffects.R',
                        data = inputdat_final,
                        n.adapt = 20000,
                        n.chains = 3,
                        inits = inits)
update(mod, n.iter = 20000)
mod_res_final = rjags::coda.samples(
  mod, 
  variable.name = c('R', 'R.pred', 'N', 'beta', 'rmax', 'e', 'fit', 'fit.new',
                    'sigma.o', 'sigma.p', 'y.sim', 'R.sim', 'N.sim',
                    'pvalue', 'pvalue.mean', 'pvalue.sd', 'sq', 'sqsim'), 
  n.iter = 100000, thin = 25)
save(inputdat_final, mod_res_final, file = 'output/model_results.RData')

## INSPECT OUTPUT-------
# model parameters:
MCMCvis::MCMCsummary(mod_res_final, 
                     params = c('beta', 'rmax', 'sigma.o', 'sigma.p'), 
                     HPD = TRUE, pg0 = TRUE)
# convergence good, n.eff good
# beta1 (flowt) may still be negative (p=0.23)
# beta2 (flowt1) still positive (p=0.99)
# beta3 (drought1) still neutral (p=0.56 to be positive)
# beta4 (long-term trend) still neutral (p=0.51)
# beta5 (density dependence) negative (p=0.99)
# rmax 0.48 (on the log scale)

MCMCvis::MCMCtrace(mod_res_final,
                   params = c('beta', 'rmax', 'sigma.o', 'sigma.p'),
                   pdf = FALSE)
# trace plots could be a bit better for beta[1] and rmax, but rmax no longer
# bumping up against 1 after extending range for prior to -2 to 2 (from -1 to 1)
# and beta[1] no longer bumping up against lower limit

# model checking:
MCMCvis::MCMCsummary(mod_res_final, 
                     params = c('pvalue', 'pvalue.mean', 'pvalue.sd'), 
                     HPD = TRUE, pg0 = TRUE)
# same as before
MCMCvis::MCMCplot(mod_res_final, params = c('sq', 'sqsim'), HPD = TRUE)
MCMCvis::MCMCplot(mod_res_final, 
                  params = c('pvalue', 'pvalue.mean', 'pvalue.sd'), HPD = TRUE)

# check for autocorrelation in the residuals
e = MCMCpstr(mod_res_final, params = "e", func = mean)
acf(e$e)
# still no evidence of autocorrelation

## ESTIMATE K?--------
# K ~ -rmax/beta[5]

Kdat = MCMCpstr(mod_res_final, c('rmax', 'beta[5]'), ISB = FALSE, type = 'chains')
Kdat$K = -Kdat$rmax/Kdat$`beta[5]`
Kdat_coda = as.mcmc(Kdat[['K']]) |> as.mcmc.list()

# PARAMETER ESTIMATES--------
model_stats = MCMCvis::MCMCsummary(mod_res_final, 
                                   params = c('beta', 'rmax'), 
                                   HPD = TRUE, pg0 = TRUE) |> 
  as_tibble(rownames = 'var') |> 
  mutate(varname = c(names(inputdat_final$predmatrix), 'rmax')) |> 
  bind_rows(
    MCMCvis::MCMCsummary(Kdat_coda, HPD = TRUE, pg0 = TRUE) |> 
      as_tibble() |> mutate(varname = 'K')) |> 
  select(varname, everything())

write_csv(model_stats, 'output/model_stats.csv')

# Estimates of N and R
N_est = MCMCvis::MCMCsummary(mod_res_final,  params = c('N'), HPD = TRUE) |> 
  as_tibble(rownames = 'parameter') |> 
  mutate(WY = c(1999:2023))
N_est |> filter(WY == 2006)

R_est = MCMCvis::MCMCsummary(mod_res_final,  params = c('R'), HPD = TRUE) |> 
  as_tibble(rownames = 'parameter') |> 
  mutate(WY = c(1999:2023))
R_est |> filter(WY %in% c(2006, 2007))
