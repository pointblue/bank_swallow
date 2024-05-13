# FULL MODEL: dynamic population model with effects on annual growth rate
#
# HYPOTHESES---------
# H1. population growth (t-1 to t) is influenced by recruitment via an increase
# in carrying capacity during the current breeding season (t), influenced by one
# or more of:
#  1a) extent of rip rap in the current water year (t)
#  1b) river flow during the immediate prior winter/same water year (t)
#
# H2. population growth is influenced by recruitment via density-dependent nest
# success in the previous year, due to varying carrying capacity (the number of
# available nest sites), which may be influenced by one or more of:
#   1a) extent of rip rap in the previous winter/water year (t-1)
#   1b) river flow in the previous one or more years (t-1)
#
# H3. population growth is influenced by recruitment via nest success in the
# previous year due to varying habitat quality (freshness of nest sites),
# which may be influenced by river flow in the previous one or more years
#
# H4. population growth is influenced via nest success, juvenile and adult
# survival in the previous year, which may be influenced by drought during the
# previous breeding season/fall prior to migration

# load packages and functions
source('R/packages.R')
#source('R/functions.R')

modeldat = read_csv('data/modeldat.csv') |> filter(WY >= 1999)

# same as for mini model
inits = list(
  list(sigma.r = 0.40, sigma.p = 0.05, sigma.o = 0.01),
  list(sigma.r = 0.25, sigma.p = 0.01, sigma.o = 0.10),
  list(sigma.r = 0.10, sigma.p = 0.10, sigma.o = 0.05))


# INPUT DATA----------
# include 1999 so that we can access the burrow count in 1999 to inform model
# from 2000 onward (n = 25); all other vectors need to be the same length
inputdat1 = list(
  # y = all observed burrow counts, 1999 through 2023; use original (unrounded)
  # data
  y = modeldat |> pull(burrows_orig), #|> scale(center = FALSE) |> as.vector(),

  # effects on R driven by current flow, prior flow, and prior drought
  predictors = modeldat |> 
    select(log.flowt, log.flowt1, log.flowt2, log.flowt3, drought1, WY) |>
    mutate_all(~scale(.))
)

# FIT MODEL-----------
mod1d = rjags::jags.model(file = 'models/1_full_model_Nchain.R',
                         data = inputdat1,
                         n.adapt = 20000,
                         n.chains = 3,
                         inits = inits)
update(mod1d, n.iter = 20000)
mod1d_res = rjags::coda.samples(
  mod1d, variable.name = c('R', 'N', 'beta', 'rmax', 'e', 'fit', 'fit.new',
                          'sigma.r', 'sigma.o', 'sigma.p', 'y.sim', 'r.sim', 
                          'pvalue', 'pvalue.mean', 'pvalue.sd', 'sq', 'sqsim'), 
  n.iter = 100000, thin = 25)
save(inputdat1, mod1d_res, file = 'output/model1d_20240509.RData')

# INSPECT OUTPUT------------
# key model parameters:
MCMCvis::MCMCsummary(mod1d_res, 
                     params = c('beta', 'rmax', 'sigma.r', 'sigma.o', 'sigma.p'), 
                     HPD = TRUE, pg0 = TRUE)
# convergence good 
# n.eff pretty good (could be higher for beta[1] and rmax
# beta1 (density dependence) negative (p=0.99)
# beta2 (flowt) may be negative (p=0.66)
# beta3 (flowt1) positive (p=0.96)
# beta4 (flowt2) positive (p=0.86)
# beta5 (flowt3) may be positive (p=0.78)
# beta6 (drought1) neutral (p=0.61 to be positive)
# beta7 (trend in density dependence) neutral
# rmax 0.79 (on the log scale) = 2.20

MCMCvis::MCMCtrace(mod1d_res,
                   params = c('beta', 'rmax', 'sigma.r', 'sigma.o', 'sigma.p'),
                   pdf = FALSE)
# trace plots could be a bit better for beta[1] and rmax, but rmax no longer
# bumping up against 1 after extending range for prior to -2 to 2 (from -1 to 1)
# and beta[1] no longer bumping up against lower limit

# model checking:
MCMCvis::MCMCsummary(mod1d_res, c('pvalue', 'pvalue.mean', 'pvalue.sd',
                                  'sq', 'sqsim'), HPD = TRUE, pg0 = TRUE)
# mod1d: pvalue.mean and pvalue close to 0.5 (good); pvalue.sd a bit high
# (meaning sd is a bit higher in the original data than the simulated)
MCMCvis::MCMCplot(mod1d_res, params = c('sq', 'sqsim'), HPD = TRUE)
MCMCvis::MCMCplot(mod1d_res, params = c('pvalue', 'pvalue.mean', 'pvalue.sd'), HPD = TRUE)

# check for autocorrelation in the residuals
e = MCMCpstr(mod1d_res, params = "e", func = mean)
acf(e$e)
# no evidence of autocorrelation

# PLOT OBSERVED VS. PREDICTED-----
# plot observed N with estimated n and K
p1_N = MCMCvis::MCMCsummary(mod1c_res, c('N', 'y.sim'), HPD = TRUE, pg0 = TRUE) %>% 
  as_tibble(rownames = 'var', .name_repair = 'universal') %>%
  rename(lci = `..95._HPDL`, uci = `..95._HPDU`) %>% 
  mutate(var = gsub('\\[.*$', '', var),
         WY = c(1999:2023, #N from t-1 through last t
                1999:2023)) %>% 
  # add observed data
  bind_rows(modeldat %>% select(WY, mean = burrowst) %>% mutate(var = 'observed')) %>%
  #filter(var != 'K') %>%
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = var), alpha = 0.5) +
  geom_line(aes(color = var)) + 
  geom_point(aes(color = var)) + 
  scale_color_manual(values = c('observed' = 'red', 'N' = 'blue4', 'y.sim' = 'gray40')) +
  scale_fill_manual(values = c('N' = 'lightblue4', 'y.sim' = 'gray70')) +
  theme_bw() + ylim(0, NA) 

# plot observed agr with estimated R and r.sim
p1_agr = MCMCvis::MCMCsummary(mod1c_res, c('R', 'r.sim'), HPD = TRUE, pg0 = TRUE) %>% 
  as_tibble(rownames = 'var', .name_repair = 'universal') %>%
  rename(lci = `..95._HPDL`, uci = `..95._HPDU`) %>% 
  mutate(var = gsub('\\[.*$', '', var),
         WY = c(1999:2023, #N from t-1 through last t
                1999:2023)) %>% 
  # add observed data
  bind_rows(modeldat %>% select(WY, mean = agr) %>% mutate(var = 'observed')) %>%
  #filter(var != 'K') %>%
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = var), alpha = 0.5) +
  geom_line(aes(color = var)) + 
  geom_point(aes(color = var)) + 
  scale_color_manual(values = c('observed' = 'red', 'R' = 'blue4', 'r.sim' = 'gray40')) +
  scale_fill_manual(values = c('R' = 'lightblue4', 'r.sim' = 'gray70')) +
  theme_bw() + ylim(0, NA)

p1_agr/p1_N
