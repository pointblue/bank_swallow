# test alternate, relatively simple mini models to build core model structure,
# with effects only of prior flow on growth rate and annual trend in K

# load packages
source('R/packages.R')

modeldat = read_csv('data/modeldat.csv') |> filter(WY >= 1999)

# explore variation in response vars

mean(modeldat$log.agr, na.rm = TRUE) # mean = -0.025
sd(modeldat$log.agr, na.rm = TRUE) # sd = 0.255
range(modeldat$log.agr, na.rm = TRUE) # -0.548 to 0.421

# simulate = function(sigma) {
#   draws = rnorm(10000, mean = 0, sd = sigma)
#   hist(draws, freq = F, main = paste('random draws for sd =', sigma),
#        breaks = 20)
#   #hist(plogis(draws), freq = F, main = 'inverse logit of draws', breaks = 20)
# }
hist(modeldat$log.agr)
simulate(0.255)
simulate(0.5)
simulate(1)

mean(modeldat$burrows_orig, na.rm = TRUE) #14925.08
var(modeldat$burrows_orig, na.rm = TRUE) #11685030 (much greater than the mean)
range(modeldat$burrows_orig, na.rm = TRUE) #7836 to 20299
# compare to Poisson
hist(modeldat$burrows_orig, xlim = c(5000, 25000))
hist(rpois(25, lambda = 14925.08), col = 'red', alpha = 0.5, add = TRUE)
# compare to Normal approx
hist(rnorm(25, mean = 14925.08, sd = sqrt(14925.08)), 
     col = 'blue', alpha = 0.5, add = TRUE)
hist(rnorm(25, mean = 14925.08, sd = 3418), 
     col = 'green', alpha = 0.5, add = TRUE)

qqplot(modeldat$burrows_orig, 
       qpois(ppoints(modeldat$burrows_orig), lambda = 14925.08),
       xlim = c(7000, 21000), ylim = c(7000, 21000))
abline(0, 1, col = 2, lty = 2)
# compare to normal approx
qqplot(modeldat$burrows_orig, 
       qnorm(ppoints(modeldat$burrows_orig), mean = 14925.08, sd = sqrt(14925.08)),
       xlim = c(7000, 21000), ylim = c(7000, 21000))
abline(0, 1, col = 2, lty = 2)

# R----------
#
# use estimated growth rate as the response variable (assuming error in the
# estimation), but include density dependent reproduction in the prior year,
# based on previous burrow count (estimated with Poisson error)
#
# >> Note: this model doesn't explicitly model the count process (Also the
# counts are actually the average of multiple observers, and I don't have access
# to the original counts of each observer.) Is modeling/accounting for the count
# process beyond the scope of this study?
#

## INPUT DATA-----
# excluding 1999 because prior burrow count in 1998 unreliable, and
# therefore so is the growth rate from 1998 to 1999 (n = 24)
inputdat_r = list(
  # y1 = observed prior burrow count (t-1) [treated as estimate of N1 (true count)]
  y1 = modeldat |> filter(WY > 1999) |> pull(burrowst1),
  
  # r = observed/estimated annual growth rate from t-1 to t [treated as estimate of R (true growth)]
  r = modeldat |> filter(WY > 1999) |> pull(agr),
  
  # prior year flow (use log-transformed version)
  flowt1 = modeldat |> filter(WY > 1999) |> pull(log.flowt1) |> scale() |> 
    as.vector(),

  # annual trend
  year = modeldat |> filter(WY > 1999) |> pull(WY) |> scale() |> as.vector()
)

## simulation---------
# start by simulating some data to gauge initial conditions

# create vectors of n, mu with length nyears + 1
nyr <- length(inputdat_r$r)
r <- numeric(nyr)
R <- numeric(nyr)

# estimate priors (edit as needed)
sigma.r <- modeldat |> drop_na() |> pull(agr) |> log() |> sd() # ~ 0.25
sigma.r <- 0.5
sigma.r <- 0.1

# simulate values
for (t in 1:length(inputdat_r$r)) {
  R[t] <- rlnorm(1, log(1.0005), sigma.r)
  r[t] <- rpois(1, lambda = R[t])
}
# compare - is it reasonable?
ggplot(modeldat |> filter(WY>1999), aes(WY, agr)) + geom_point() + geom_line() +
  geom_line(aes(y = R), col = 'red')

## fit model-------
inits_r = list(
  list(sigma.r = 0.40),
  list(sigma.r = 0.25),
  list(sigma.r = 0.10))

mod_r = rjags::jags.model(file = 'models/0_mini_model_R.R',
                          data = inputdat_r,
                          n.adapt = 20000,
                          n.chains = 3,
                          inits = inits_r)
update(mod_r, n.iter = 20000)
mod_r_res = rjags::coda.samples(
  mod_r, variable.name = c('R', 'r', 'N1', 'beta', 'rmax', 'er', 'y.sim', 'r.sim', 
                           'sigma.o', 'sigma.r'), 
  n.iter = 100000, thin = 25)

MCMCvis::MCMCsummary(mod_r_res, 
                     params = c('beta', 'rmax', 'sigma.o', 'sigma.r'), 
                     HPD = TRUE, pg0 = TRUE)
MCMCvis::MCMCtrace(mod_r_res,
                   params = c('beta', 'rmax', 'sigma.o', 'sigma.r'),
                   pdf = FALSE)

# plot observed N with estimated N1, K1
mod_r_N = MCMCvis::MCMCsummary(mod_r_res, c('N1', 'y.sim'), HPD = TRUE, pg0 = TRUE) %>% 
  as_tibble(rownames = 'var', .name_repair = 'universal') %>%
  rename(lci = `..95._HPDL`, uci = `..95._HPDU`) %>% 
  mutate(var = gsub('\\[.*$', '', var),
         WY = c(1999:2022,
                1999:2022)) |>
  # add observed data
  bind_rows(modeldat %>% select(WY, mean = burrowst) %>% mutate(var = 'observed')) %>%
  #filter(var != 'K') %>%
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = var), alpha = 0.5) +
  geom_line(aes(color = var)) + 
  geom_point(aes(color = var)) + 
  scale_color_manual(values = c('observed' = 'red', 'N1' = 'blue4', 'y.sim' = 'gray40')) +
  scale_fill_manual(values = c('N1' = 'lightblue4', 'y.sim' = 'gray70')) +
  theme_bw() + ylim(0, NA) 

# plot observed and estimated agr
mod_r_agr = MCMCvis::MCMCsummary(mod_r_res, c('R', 'r.sim'), HPD = TRUE, pg0 = TRUE) %>% 
  as_tibble(rownames = 'var', .name_repair = 'universal') %>%
  rename(lci = `..95._HPDL`, uci = `..95._HPDU`) %>% 
  mutate(var = gsub('\\[.*$', '', var),
         WY = c(2000:2023, 2000:2023)) %>% 
  # add observed data
  bind_rows(modeldat %>% select(WY, mean = agr) %>% mutate(var = 'observed')) %>%
  #filter(var != 'K') %>%
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = var), alpha = 0.5) +
  geom_line(aes(color = var)) + 
  geom_point(aes(color = var)) + 
  scale_color_manual(values = c('observed' = 'red', 'R' = 'blue4', 'r.sim' = 'dodgerblue')) +
  scale_fill_manual(values = c('R' = 'lightblue4', 'r.sim' = 'lightblue')) +
  theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = 1), linetype = 'dashed')


# N chain-------

# attempt an alternative dynamic population model that more explicitly accounts
# for error in the annual counts, and constrains growth as the relationship
# between the "true" counts in subsequent years (rather than an independent
# quantity)
#
# >> Note: model requires assuming burrow count in 1999 is "true" without error

## INPUT DATA----
# include 1999 so that we can access the burrow count in 1999 to
# inform model from 2000 onward (n = 25); all other vectors need to be the same
# length
inputdat_chain = list(
  # y = all observed burrow counts, 1999 through 2023; use original (unrounded)
  # data
  y = modeldat |> pull(burrows_orig), #|> scale(center = FALSE) |> as.vector(),
  
  # prior year flow (use log-transformed version)
  flowt1 = modeldat |> pull(log.flowt1) |> scale() |> as.vector(),
  
  # annual trend
  year = modeldat |> pull(WY) |> scale() |> as.vector()
)

## simulation---------
# start by simulating some data to gauge initial conditions

# create vectors of n, mu with length nyears + 1
nyr <- length(inputdat_chain$y)
n <- numeric(nyr)
mu <- numeric(nyr)
y <- numeric(nyr)
R <- numeric(nyr)

# estimate priors (edit as needed)
sigma.r <- modeldat |> drop_na() |> pull(agr) |> log() |> sd() # ~ 0.25
sigma.p <- 0.05
sigma.o <- 0.05
n[1] <- inputdat_chain$y[1]
mu[1] <- n[1] # mean from deterministic model to simulate
y[1] <- n[1]
R[1] <- modeldat$agr[1]

# simulate values
for (t in 2:length(inputdat_chain$y)) {
  R[t] <- rlnorm(1, log(1.0005), sigma.r)
  mu[t] <- R[t] * n[t - 1]
  n[t] <- rlnorm(1, log(mu[t]), sigma.p)
  y[t] <- rlnorm(1, log(n[t]), sigma.o)
}
# compare - is it reasonable?
ggplot(modeldat, aes(WY, burrows_orig)) + geom_point() + geom_line() +
  geom_line(aes(y = n), col = 'red') + 
  geom_line(aes(y = y), col = 'blue')

## fit model-----
inits = list(
  list(sigma.r = 0.40, sigma.p = 0.05, sigma.o = 0.01),
  list(sigma.r = 0.25, sigma.p = 0.01, sigma.o = 0.10),
  list(sigma.r = 0.10, sigma.p = 0.10, sigma.o = 0.05))

mod_Nchain = rjags::jags.model(file = 'models/0_mini_model_Nchain_v2.R',
                               data = inputdat_chain,
                               n.adapt = 20000,
                               n.chains = 3,
                               inits = inits)
update(mod_Nchain, n.iter = 20000)
mod_Nchain_res = rjags::coda.samples(
  mod_Nchain, variable.name = c('R', 'N', 'beta', 'rmax', 'sigma.r', 'sigma.o', 
                                'sigma.p', 'y.sim', 'r.sim'), 
  n.iter = 100000, thin = 25)

MCMCvis::MCMCsummary(mod_Nchain_res, 
                     params = c('beta', 'rmax', 'sigma.r', 'sigma.o', 'sigma.p'), 
                     HPD = TRUE, pg0 = TRUE)
MCMCvis::MCMCtrace(mod_Nchain_res,
                   params = c('beta', 'rmax', 'sigma.r', 'sigma.o', 'sigma.p'),
                   pdf = FALSE)

# plot observed N with estimated n and K
mod_Nchain_N = MCMCvis::MCMCsummary(mod_Nchain_res, c('N', 'y.sim'), HPD = TRUE, pg0 = TRUE) %>% 
  as_tibble(rownames = 'var', .name_repair = 'universal') %>%
  rename(lci = `..95._HPDL`, uci = `..95._HPDU`) %>% 
  mutate(var = gsub('\\[.*$', '', var),
         WY = c(1999:2023,
                #1999:2022,
                2000:2023)) %>% 
  # add observed data
  bind_rows(modeldat %>% select(WY, mean = burrowst) %>% mutate(var = 'observed')) %>%
  #filter(var != 'K') %>%
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = var), alpha = 0.5) +
  geom_line(aes(color = var)) + 
  geom_point(aes(color = var)) + 
  scale_color_manual(values = c('observed' = 'red', 'N' = 'blue4', 'y.sim' = 'dodgerblue')) +
  scale_fill_manual(values = c('N' = 'lightblue4', 'y.sim' = 'lightblue')) +
  theme_bw() + ylim(0, NA) #+
  #facet_wrap(~var, ncol = 1)
# now N1 and N directly depend on each other

# >> N-chain option seems a lot cleaner/simpler, but n.eff is tiny; might
# require excessive iterations

# plot observed and estimated agr
mod_Nchain_agr = MCMCvis::MCMCsummary(mod_Nchain_res, c('R', 'r.sim'), HPD = TRUE, pg0 = TRUE) %>% 
  as_tibble(rownames = 'var', .name_repair = 'universal') %>%
  rename(lci = `..95._HPDL`, uci = `..95._HPDU`) %>% 
  mutate(var = gsub('\\[.*$', '', var),
         WY = c(2000:2023, 2000:2023)) %>% 
  # add observed data
  bind_rows(modeldat %>% select(WY, mean = agr) %>% mutate(var = 'observed')) %>%
  #filter(var != 'K') %>%
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = var), alpha = 0.5) +
  geom_line(aes(color = var)) + 
  geom_point(aes(color = var)) + 
  scale_color_manual(values = c('observed' = 'red', 'R' = 'blue4', 'r.sim' = 'dodgerblue')) +
  scale_fill_manual(values = c('R' = 'lightblue4', 'r.sim' = 'lightblue')) +
  theme_bw() + ylim(0, NA) +
  geom_hline(aes(yintercept = 1), linetype = 'dashed') #+
#facet_wrap(~var, ncol = 1)
#


# compare----------
library(patchwork)
mod_r_agr/mod_Nchain_agr # these look extremely (suspiciously?) similar!
mod_r_N/mod_Nchain_N # variance in Nchain model looks more reasonable (suspiciously narrow for agr)

MCMCvis::MCMCsummary(mod_r_res, 
                     params = c('beta', 'rmax', 'sigma.o', 'sigma.r'), 
                     HPD = TRUE, pg0 = TRUE)
MCMCvis::MCMCsummary(mod_Nchain_res, 
                     params = c('beta', 'rmax', 'sigma.r', 'sigma.o', 'sigma.p'), 
                     HPD = TRUE, pg0 = TRUE)
