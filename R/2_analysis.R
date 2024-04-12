# MODEL---------
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
source('R/functions.R')
#source('R/fit_model.R')
#source('R/fit_BANS_model.R')

birddat = read_csv('data/birddat.csv')
flowdat = read_csv('data/flowdat_20220930.csv')
droughtdat = read_csv('data/droughtdat.csv')
riprapdat = read_csv('data/riprapdat.csv')


# Generate model predictors----------

## flow-------

flowdat_all = flowdat %>% 
  select(WY, flowt = flowtotal, flow14t = flowtotal_threshold) %>% 
  # extend to 2023 and 2024
  bind_rows(tibble(WY = c(2023:2025))) %>% 
  # find lag values
  mutate(flowt1 = lag(flowt),
         flowt2 = lag(flowt, n = 2),
         flowt3 = lag(flowt, n = 3),
         # repeat for total flow above threshold
         flow14t1 = lag(flow14t),
         flow14t2 = lag(flow14t, n = 2),
         flow14t3 = lag(flow14t, n = 3)) %>% 
  select(WY, starts_with('flowt'), starts_with('flow14t'))

# check correlations:
corrplot::corrplot.mixed(
  cor(flowdat_all %>% filter(WY >= 2000) %>% drop_na(), 
      method = 'spearman'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# total and threshold totals are all highly correlated (>0.9)
# lagged values not correlated with others
# >> all flow data moderately negatively correlated with WY, most strongly with
# flowt (-0.46)

## drought--------
# Note: monthly PDSI values do not depict drought on timescales < ~12 months;
# it's designed to be highly correlated month to month to account for "land
# memory" (https://climatedataguide.ucar.edu/climate-data/palmer-drought-severity-index-pdsi)
# --> the Z index is an alternative that is more responsive to short term drought conditions

drought_sum = bind_rows(
  droughtdat %>% 
    filter(WY >= 1995) %>% 
    filter(!is.na(season)) %>% 
    group_by(index, WY, season) %>%
    summarize(value = mean(value),
              .groups = 'drop'),
  droughtdat %>% 
    filter(WY >= 1995) %>% 
    group_by(index, WY) %>%
    summarize(value = mean(value),
              .groups = 'drop') %>% 
    mutate(season = 'WY')) %>% 
  pivot_wider(names_from = c('index', 'season'))

# correlations:
corrplot::corrplot.mixed(
  cor(drought_sum %>% filter(WY >= 2000 & WY <= 2022) %>% drop_na(), method = 'pearson'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# zndx and pdsi strongly correlated within season (breeding = 0.90; rainy =
# 0.87; WY = 0.94)
# zndx breeding and rainy season metrics not strongly correlated (0.51), but
# each is strongly correlated with WY (0.87, 0.86)
# pdsi breeding and rainy strongly correlated (0.73), and both strongly 
# correlated with WY (0.90, 0.95)
# no correlation with WY for any

# focus on PDSI_breeding
drought_final = drought_sum %>% select(WY, drought = pdsidv_breeding) %>% 
  mutate(drought1 = lag(drought))
corrplot::corrplot.mixed(
  cor(drought_final %>% filter(WY >= 2000 & WY <= 2022) %>% drop_na(), method = 'pearson'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# no strong correlations

## riprap---------


## compile predictors----------------
# check for correlations among predictors

modeldat = list(birddat %>% select(WY = year, pgr, burrowst = burrows) %>% 
                  mutate(burrowst = burrowst,
                         burrowst1 = lag(burrowst)),
                flowdat_all %>% 
                  # convert to cfs (millions)
                  mutate(across(c(starts_with('flowt'), starts_with('flow14t')), 
                                ~./1000000)),
                drought_final,
                riprap_proj_spline %>% filter(WY >= 1999) %>% 
                  mutate(riprap1 = lag(riprap))) %>% 
  purrr::reduce(full_join) %>% 
  filter(WY >= 1999 & WY <= 2023) #exclude data prior to 1999
write_csv(modeldat, 'data/modeldat.csv')

# check correlations with birddat & across metrics
# no transformation
corrplot::corrplot.mixed(
  cor(modeldat %>% select(!starts_with('flow14t')) %>% drop_na(), method = 'spearman'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')
corrplot::corrplot.mixed(
  cor(modeldat %>% select(!starts_with('flowt')) %>% drop_na(), method = 'spearman'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')

# strongest: 
# - pgr & flowt1 (total or threshold) (0.69, 0.71)
# - pgr & burrows_prior (-0.62)
# - drought & flowt (0.64, 0.59); drought1 & flowt1 (0.65, 0.63)
# >> WY and riprap perfectly correlated (of course)

ggplot(modeldat, aes(flow14t1)) + geom_density()
ggplot(modeldat, aes(log(flow14t1))) + geom_density()
ggplot(modeldat, aes(flowt1)) + geom_density()
ggplot(modeldat, aes(log(flowt1))) + geom_density()

mean(modeldat$pgr, na.rm = TRUE) # mean = 0.041
var(modeldat$pgr, na.rm = TRUE) # var = 0.057

simulate = function(sigma) {
  draws = rnorm(10000, mean = 0, sd = sigma)
  hist(draws, freq = F, main = paste('random draws for var =', sigma),
       breaks = 20)
  #hist(plogis(draws), freq = F, main = 'inverse logit of draws', breaks = 20)
}


# test analysis: N mini model-----------
# model should have effects of current flow, prior flow, and drought on R, plus
# effects of riprap and prior flow on K

testdat = list(
  # y = observed burrow count (t); exclude 1999 because prior burrow count in 1998 unreliable
  y = modeldat %>% filter(WY != 1999 & WY != 2023) %>% pull(burrowst) %>% round(),
  
  # N = observed prior burrow count (t-1)
  N = modeldat %>% filter(WY != 1999 & WY != 2023) %>% pull(burrowst1) %>% round(),
  
  # effects on K driven by prior flow and rip rap
  Kpredictors = modeldat %>% filter(WY != 1999 & WY != 2023) %>% 
    select(riprap_t1 = riprap1,
           flow_t1 = flowt1,
           flow_t2 = flowt2,
           flow_t3 = flowt3) %>% 
    # transform, scale, center
    mutate(across(starts_with('flow'), ~log(. + 1)),
           across(starts_with('flow'), ~(. - 1.5)/0.3),
           riprap_t1 = (riprap_t1 - 83000)/3000),
  
  # effects on R driven by current flow, prior flow, and prior drought
  Rpredictors = modeldat %>% filter(WY != 1999 & WY != 2023) %>% 
    select(flow_t = flowt,
           flow_t1 = flowt1,
           flow_t2 = flowt2,
           flow_t3 = flowt3,
           drought_t1 = drought1) %>% 
    # transform, scale, center
    mutate(across(starts_with('flow'), ~log(. + 1)),
           across(starts_with('flow'), ~(. - 1.5)/0.3),
           drought_t1 = (drought_t1 + 1)/2),
           
  nKcov = 4, #number of covariates in Kpredictors
  nRcov = 5 #number of covariates in Bpredictors
)

test1 = rjags::jags.model(file = 'models/0_mini_model_N.R',
                          data = testdat,
                          n.adapt = 5000,
                          n.chains = 1)
update(test1, n.iter = 5000)
test1_coda = rjags::coda.samples(
  test1, variable.name = c('r', 'R', 'K', 'mu', 'n', 'pgr', 'k0', 'k1',
                           'sigma.r', 'sigma.k'), 
  n.iter = 5000, thin = 10)

MCMCvis::MCMCsummary(test1_coda, 
                     params = c('r', 'k0', 'k1', 'sigma.r', 'sigma.k'), 
                     HPD = TRUE, pg0 = TRUE)
MCMCvis::MCMCtrace(test1_coda,
                   params = c('r', 'k0', 'k1', 'sigma.r', 'sigma.k'),
                   pdf = FALSE)
MCMCvis::MCMCplot(test1_coda, params = c('k1'), HPD = TRUE)
MCMCvis::MCMCplot(test1_coda, params = c('K'), HPD = TRUE)
MCMCvis::MCMCplot(test1_coda, params = c('n'), HPD = TRUE)
MCMCvis::MCMCplot(test1_coda, params = c('mu'), HPD = TRUE)
MCMCvis::MCMCplot(test1_coda, params = c('pgr'), HPD = TRUE)
MCMCvis::MCMCplot(test1_coda, params = c('R'), HPD = TRUE)

# plot observed N with estimated n and K
MCMCvis::MCMCsummary(test1_coda, c('n', 'K', 'mu'), HPD = TRUE, pg0 = TRUE) %>% 
  as_tibble(rownames = 'var', .name_repair = 'universal') %>%
  rename(lci = `..95._HPDL`, uci = `..95._HPDU`) %>% 
  mutate(var = gsub('\\[.*$', '', var),
         WY = c(1999:2021, #n (burrows prior)
                1999:2021, #K (prior year)
                2000:2022)) %>% #mu (current year)
  # add observed data
  bind_rows(modeldat %>% select(WY, mean = burrowst) %>% mutate(var = 'observed')) %>%
  #filter(var != 'K') %>%
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = var), alpha = 0.5) +
  geom_line(aes(color = var)) + 
  geom_point(aes(color = var)) + 
  theme_bw() + ylim(0, NA)

# test analysis: N chain-------
testdat_nchain = list(
  
  # N = observed prior burrow count (t-1)
  N = modeldat %>% pull(burrowst) %>% round(),
  
  # effects on K driven by prior flow and rip rap;
  # they must be in this order!
  Kpredictors = modeldat %>% 
    select(riprap_t1 = riprap1,
           flow_t1 = flowt1,
           flow_t2 = flowt2,
           flow_t3 = flowt3) %>% 
    # transform, scale, center
    mutate(across(starts_with('flow'), ~log(. + 1)),
           across(starts_with('flow'), ~(. - 1.5)/0.3),
           riprap_t1 = (riprap_t1 - 83000)/3000),
  
  # effects on R driven by current flow, prior flow, and drought
  Rpredictors = modeldat %>% 
    select(flow_t = flowt,
           flow_t1 = flowt1,
           flow_t2 = flowt2,
           flow_t3 = flowt3,
           drought_t1 = drought1) %>% 
    # transform, scale, center
    mutate(across(starts_with('flow'), ~log(. + 1)),
           across(starts_with('flow'), ~(. - 1.5)/0.3),
           drought_t1 = (drought_t1 + 1)/2),
  
  nKcov = 4, #number of covariates in Kpredictors
  nRcov = 5 #number of covariates in Bpredictors
)

test1_chain = rjags::jags.model(file = 'models/0_mini_model_Nchain.R',
                          data = testdat_nchain,
                          n.adapt = 5000,
                          n.chains = 1)
update(test1_chain, n.iter = 5000)
test1_coda_chain = rjags::coda.samples(
  test1_chain, variable.name = c('r', 'R', 'K', 'n', 'k0', 'k1',
                           'sigma.r', 'sigma.k'), 
  n.iter = 5000, thin = 10)

MCMCvis::MCMCsummary(test1_coda_chain, 
                     params = c('r', 'k0', 'k1', 'sigma.r', 'sigma.k'), 
                     HPD = TRUE, pg0 = TRUE)
MCMCvis::MCMCtrace(test1_coda_chain,
                   params = c('r', 'k0', 'k1', 'sigma.r', 'sigma.k'),
                   pdf = FALSE)
MCMCvis::MCMCplot(test1_coda_chain, params = c('k1'), HPD = TRUE)
MCMCvis::MCMCplot(test1_coda_chain, params = c('K'), HPD = TRUE)
MCMCvis::MCMCplot(test1_coda_chain, params = c('n'), HPD = TRUE)
MCMCvis::MCMCplot(test1_coda_chain, params = c('R'), HPD = TRUE)

# plot observed N with estimated n and K
MCMCvis::MCMCsummary(test1_coda_chain, c('n', 'K'), HPD = TRUE, pg0 = TRUE) %>% 
  as_tibble(rownames = 'var', .name_repair = 'universal') %>%
  rename(lci = `..95._HPDL`, uci = `..95._HPDU`) %>% 
  mutate(var = gsub('\\[.*$', '', var),
         WY = c(1999:2023, 
                1999:2022)) %>% 
  # add observed data
  bind_rows(modeldat %>% select(WY, mean = burrowst) %>% mutate(var = 'observed')) %>%
  #filter(var != 'K') %>%
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = var), alpha = 0.5) +
  geom_line(aes(color = var)) + 
  geom_point(aes(color = var)) + 
  scale_color_manual(values = c('observed' = 'red', 'n' = 'blue4', 'K' = 'gray40')) +
  scale_fill_manual(values = c('n' = 'lightblue4', 'K' = 'gray70')) +
  theme_bw() + ylim(0, NA)

# ANALYSIS 1-----------
# using flow totals (vs threshold)

inputdat_nchain_totals = list(
  
  # N = observed prior burrow count (t-1)
  N = modeldat %>% filter(WY != 2023) %>%pull(burrowst) %>% round(),
  
  # effects on K driven by rip rap
  Kpredictors = modeldat %>% filter(WY != 2023) %>% 
    select(riprap_t1 = riprap1) %>% 
    # transform, scale, center
    mutate(riprap_t1 = (riprap_t1 - 83000)/3000),
  
  # effects on R driven by current flow, prior flow, and drought
  Rpredictors = modeldat %>% filter(WY != 2023) %>%
    select(flow_t = flowt,
           flow_t1 = flowt1,
           #flow_t2 = flowt2,
           #flow_t3 = flowt3,
           #drought_t1 = drought1
           ) %>% 
    # transform, scale, center
    mutate(across(starts_with('flow'), ~log(. + 1)),
           across(starts_with('flow'), ~(. - 1.5)/0.3),
           #drought_t1 = (drought_t1 + 1)/2
           ),
  
  nKcov = 1, #number of covariates in Kpredictors
  nRcov = 2 #number of covariates in Bpredictors
)

jm1 = rjags::jags.model(file = 'models/1_Nchain_covariates.R',
                        data = inputdat_nchain_totals,
                        n.adapt = 50000,
                        n.chains = 3)
update(jm1, n.iter = 100000)
res1 = rjags::coda.samples(
  jm1, variable.name = c('r', 'R', 'n', 'k0', 'beta',
                         'sigma.r'#, 'k', 'K', 'sigma.k'), 
                        ),
  n.iter = 300000, thin = 100)

MCMCvis::MCMCsummary(
  res1, c('r', 'k0', 'beta', 'sigma.r'), #, 'k', 'sigma.k'),
  HPD = TRUE, pg0 = TRUE)
MCMCvis::MCMCtrace(
  res1,
  params = c('r', 'k0', 'beta', 'sigma.r'), #, 'k', 'sigma.k'),
  pdf = FALSE)

MCMCvis::MCMCplot(res1, params = c('k'), HPD = TRUE)
MCMCvis::MCMCplot(res1, params = c('beta'), HPD = TRUE)
MCMCvis::MCMCplot(res1, params = c('K'), HPD = TRUE)
MCMCvis::MCMCplot(res1, params = c('n'), HPD = TRUE)
MCMCvis::MCMCplot(res1, params = c('R'), HPD = TRUE)

# plot observed N with estimated n and K
bind_rows(
  MCMCvis::MCMCsummary(res1, 'n', HPD = TRUE, pg0 = TRUE) %>% 
    as_tibble(rownames = 'var', .name_repair = 'universal') %>% 
    mutate(WY = c(1999:2022)),
  MCMCvis::MCMCsummary(res1, 'k0', HPD = TRUE, pg0 = TRUE) %>% 
    as_tibble(rownames = 'var', .name_repair = 'universal') %>% 
    expand_grid(WY = c(1999:2022))) %>%
  rename(lci = `...4`, uci = `...5`) %>% 
  mutate(var = gsub('\\[.*$', '', var)) %>% 
  # add observed data
  bind_rows(modeldat %>% filter(WY != 2023) %>% 
              select(WY, mean = burrowst) %>% mutate(var = 'observed')) %>%
  #filter(var != 'K') %>%
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = var), alpha = 0.5) +
  geom_line(aes(color = var)) + 
  geom_point(aes(color = var)) + 
  scale_color_manual(values = c('observed' = 'red', 'n' = 'blue4', 'k0' = 'gray40')) +
  scale_fill_manual(values = c('n' = 'lightblue4', 'k0' = 'gray70')) +
  theme_bw() + ylim(0, NA)



# simulated vs. observed pgr & burrows-----------
sim1b = list(
  'pgr' = MCMCvis::MCMCsummary(res1b, c('ysim', 'ysim.2007'),
                               HPD = TRUE, pg0 = TRUE) %>% 
    mutate(WY = c(2000:2006, 2008:2022, 2007)) %>%
    left_join(modeldat %>% select(WY, pgr)) %>% 
    arrange(WY),
  'N' = MCMCvis::MCMCsummary(res1b, c('Nsim'),
                             HPD = TRUE, pg0 = TRUE) %>% 
    mutate(WY = c(1999:2022)) %>%
    full_join(modeldat %>% select(WY, burrowst)) %>% 
    arrange(WY))

plot_simulations(sim = sim1b, ylim2 = c(0, 28))
ggsave(filename = 'fig/simulations.png', width = 6, height = 6, units = 'in')

## predicted pgr over time ---------
# (holding betas constant/at mean but assuming K changes)
MCMCvis::MCMCsummary(res1c, c('ypred2'), HPD = TRUE) %>% 
  as_tibble(.name_repair = 'universal') %>% 
  rename(lci = ..95._HPDL, uci = ..95._HPDU) %>% 
  mutate(value = inputdat1$predmatrix2$t,
         value_orig = value * 5 + 55) %>% 
  ggplot(aes(value_orig, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = 'gray85') + 
  geom_line() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ylim(-1, 1) +
  geom_point(data = modeldat, aes(riprap, pgr)) +
  theme_bw()

## predicted carrying capacity over time---------
MCMCvis::MCMCsummary(res1c, c('K'), HPD = TRUE) %>% 
  as_tibble(.name_repair = 'universal') %>% 
  rename(lci = ..95._HPDL, uci = ..95._HPDU) %>% 
  mutate(N = inputdat1$N,
         WY = modeldat %>% filter(!WY %in% c(1999, 2007)) %>% pull(WY)) %>% 
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci),  fill = 'gray80') +
  geom_line() + geom_point(aes(y = N)) +
  #ylim(0, 3) +
  theme_bw()

MCMCvis::MCMCsummary(res1c, c('k.pred'), HPD = TRUE) %>% 
  as_tibble(.name_repair = 'universal') %>% 
  rename(lci = ..95._HPDL, uci = ..95._HPDU) %>% 
  mutate(value = inputdat1$predmatrix2$t,
         value_orig = value * 5 + 55) %>% 
  ggplot(aes(value_orig, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = 'gray85') + 
  geom_line() + theme_bw() +
  labs(x = 'riprap (miles)', y = 'carrying capacity (burrows)') +
  #ylim(0, 4) +
  geom_point(data = modeldat, aes(riprap, burrowst)) + 
  theme_bw()

# ANALYSIS 2-----------
# using total flow data

inputdat2 = list(
  # y = per-capita growth rate t-1 to t
  # >> exclude 1999 because data from 1998 less reliable, so estimate of pgr
  # from 98 to 99 is also unreliable
  # >> exclude 2007 because burrows_prior (from 2006) is missing to estimate pgr
  # >> exclude 2021
  y = modeldat %>% filter(!WY %in% c(1999, 2007)) %>% pull(pgr),
  
  # density-dependent reproduction driven by burrows_prior
  N = modeldat %>% filter(!WY %in% c(1999, 2007)) %>% pull(burrowst1),
  
  # consider trend in density dependence (due to changes in rip rap)
  Kpredictors = modeldat %>% filter(!WY %in% c(1999, 2007)) %>%
    # center on 55 with sd of 5
    select(riprap) %>% mutate(riprap = (riprap - 55)/5),
  # center on 2010 with sd of 5
  #select(t = WY) %>% mutate(t = (t - 2010)/5),
  
  # Betas = density-independent predictors of growth rate from t-1 to t
  # >> exclude 1999 and 2007 as above
  Bpredictors = modeldat %>% filter(!WY %in% c(1999, 2007)) %>%  
    select(flow_t = flowt,
           flow_t1 = flowt1,
           flow_t2 = flowt2,
           flow_t3 = flowt3) %>%
    mutate(
      across(starts_with('flow'), ~log(. + 1)) # flowdata: log transform
      #across(everything(), scale) # everything: center and scale
    ) %>%
    as.matrix() %>% scale(),
  
  nBcov = 4,
  
  # range of scaled values for beta predictors
  predmatrix = tibble(flow_t = seq(-1.4, 2.1, length.out = 11),
                      flow_t1 = seq(-1.6, 2.4, length.out = 11),
                      flow_t2 = seq(-1.4, 2.1, length.out = 11),
                      flow_t3 = seq(-1.6, 2.1, length.out = 11)),
  npred = 11,
  predmatrix2 = tibble(t = seq(-2, 2, length.out = 11)),
  npred2 = 11
)

## riprap-------
jm2 = rjags::jags.model(file = 'R/full_model_sim.R',
                        data = inputdat2,
                        n.adapt = 5000,
                        n.chains = 3)
update(jm2, n.iter = 5000)
res2 = rjags::coda.samples(
  jm2, variable.name = c('r', 'k', 'K', 'k0', 'k1', 'sigma', 'beta', 
                         'ysim', 'ypred', 'k.pred', 'ypred2',
                         'pvalue.mean', 'pvalue.sd', 'pvalue.fit'), 
  n.iter = 100000, thin = 10)

MCMCvis::MCMCsummary(
  res2, c('r', 'k0', 'k1', 'sigma', 'beta', 
          'pvalue.mean', 'pvalue.sd', 'pvalue.fit'), 
  HPD = TRUE, pg0 = TRUE)
MCMCvis::MCMCplot(res2, params = c('beta'), HPD = TRUE)

## constant-------
jm2b = rjags::jags.model(file = 'R/full_model_simple.R',
                        data = inputdat2,
                        n.adapt = 5000,
                        n.chains = 3)
update(jm2b, n.iter = 5000)
res2b = rjags::coda.samples(
  jm2b, variable.name = c('r', 'k0', 'sigma', 'beta', 
                         'ysim', 'ypred', 
                         'pvalue.mean', 'pvalue.sd', 'pvalue.fit'), 
  n.iter = 100000, thin = 10)

MCMCvis::MCMCsummary(
  res2b, c('r', 'k0', 'sigma', 'beta', 
          'pvalue.mean', 'pvalue.sd', 'pvalue.fit'), 
  HPD = TRUE, pg0 = TRUE)
MCMCvis::MCMCplot(res2b, params = c('beta'), HPD = TRUE)

## effects plots------
full_join(MCMCvis::MCMCsummary(res2b, c('ypred'), HPD = TRUE) %>% 
            as_tibble(.name_repair = 'universal') %>% 
            mutate(var = rep(dimnames(inputdat2$Bpredictors)[[2]], each = inputdat2$npred),
                   index = rep(c(1:11), inputdat2$nBcov)),
          inputdat2$predmatrix %>% mutate(index = c(1:inputdat2$npred)) %>% 
            pivot_longer(-index, names_to = 'var') %>% 
            left_join(
              attributes(inputdat2$Bpredictors)$`scaled:center` %>% 
                as_tibble(rownames = 'var') %>% set_names(c('var', 'center'))) %>% 
            left_join(
              attributes(inputdat2$Bpredictors)$`scaled:scale` %>% 
                as_tibble(rownames = 'var') %>% set_names(c('var', 'scale'))) %>% 
            mutate(value_orig = value * scale + center,
                   value_orig = if_else(var %in% c('flow_t', 'flow_t1', 'flow_t2', 'flow_t3'),
                                        exp(value_orig)-1,
                                        value_orig)),
          by = c('index', 'var')) %>% 
  rename(pgr = mean, lci = ..95._HPDL, uci = ..95._HPDU) %>% 
  #filter(!var %in% 'drought1') %>%
  #filter(value_orig > 1) %>% 
  ggplot(aes(value_orig, pgr)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci),
              fill = 'gray85') + geom_line() + facet_wrap(~var, ncol = 1) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point(data = modeldat %>% filter(!WY %in% c(1999, 2006, 2007)) %>% 
               select(pgr, starts_with('flowt')) %>% 
               pivot_longer(-pgr, names_to = 'var', values_to = 'value_orig') %>% 
               mutate(var = gsub('flowt', 'flow_t', var))) +
  #ylim(-1, 1) +
  xlim(0, 10) +
  theme_bw()

## predicted pgr over time ---------
# (holding betas constant/at mean but assuming K changes)
MCMCvis::MCMCsummary(res2, c('ypred2'), HPD = TRUE) %>% 
  as_tibble(.name_repair = 'universal') %>% 
  rename(lci = ..95._HPDL, uci = ..95._HPDU) %>% 
  mutate(value = inputdat2$predmatrix2$t,
         value_orig = value * 5 + 55) %>% 
  ggplot(aes(value_orig, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = 'gray85') + 
  geom_line() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ylim(-1, 1) +
  geom_point(data = modeldat, aes(riprap, pgr)) +
  theme_bw()

## predicted carrying capacity over time---------
MCMCvis::MCMCsummary(res2, c('K'), HPD = TRUE) %>% 
  as_tibble(.name_repair = 'universal') %>% 
  rename(lci = ..95._HPDL, uci = ..95._HPDU) %>% 
  mutate(N = inputdat2$N,
         WY = modeldat %>% filter(!WY %in% c(1999, 2007)) %>% pull(WY)) %>% 
  ggplot(aes(WY, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci),  fill = 'gray80') +
  geom_line() + geom_point(aes(y = N)) + 
  theme_bw()

MCMCvis::MCMCsummary(res2, c('k.pred'), HPD = TRUE) %>% 
  as_tibble(.name_repair = 'universal') %>% 
  rename(lci = ..95._HPDL, uci = ..95._HPDU) %>% 
  mutate(value = inputdat2$predmatrix2$t,
         value_orig = value * 5 + 55) %>% 
  ggplot(aes(value_orig, mean)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = 'gray85') + 
  geom_line() + theme_bw() +
  labs(x = 'riprap (miles)', y = 'carrying capacity (burrows)') +
  #ylim(0, 1) +
  geom_point(data = modeldat, aes(riprap, burrowst)) +
  theme_bw()



# MODEL COMPARISONS---------
dic1 = dic.samples(model = jm1, n.iter = 100000, thin = 10, type = 'pD')
dic1b = dic.samples(model = jm1b, n.iter = 100000, thin = 10, type = 'pD')
dic2 = dic.samples(model = jm2, n.iter = 100000, thin = 10, type = 'pD')
dic2b = dic.samples(model = jm2b, n.iter = 100000, thin = 10, type = 'pD')

dic1 - dic2 # 4; slight preference for mod2 over mod1
dic1 - dic1b # 1.8; slight pref for mod1b over mod1
dic2 - dic2b # -0.03
dic1b - dic2 # 2; pref for mod2 over all


# FINALIZE & INTERPRET--------
# mod2c
mod_final = mod2c

MCMCvis::MCMCtrace(mod_final$results,
                   params = c('r0', 'beta', 'k0', 'var.o'),
                   pdf = FALSE)
# fits well; n.eff sufficient; trace plots well-mixed;

effect_stats = MCMCvis::MCMCsummary(mod_final$results,
                     c('r0', 'beta', 'k0', 'k1', 'var.o', #'var.p', 
                       'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
                     HPD = TRUE, 
                     pg0 = TRUE) %>% 
  as_tibble(rownames = 'var') %>% 
  mutate(predictor = c('growth', names(inputdat2$Bpredictors), 
                       'density_int', 'density_trend',
                       'var.obs', #'var.process', 
                       'pvalue.mean', 'pvalue.sd', 'pvalue.fit')) %>% 
  select(var, predictor, everything())
write_csv(effect_stats, 'output/effect_stats.csv')

# plot effect sizes
effect_stats %>% 
  filter(grepl('beta|k', var)) %>% 
  select(predictor, mean, lcl = `95%_HPDL`, ucl = `95%_HPDU`, pg0 = `p>0`) %>% 
  mutate(label = if_else(round(pg0, 2) >= 0.9 | round(pg0, 2) <= 0.1, 
                         '*', ''),
         predictor = factor(predictor, 
                            levels = rev(c('flow_concurrent', 'flow_prior', 'flow_older2',
                                       'drought_prior', 
                                       'density_int', 'density_trend')))) %>% 
  ggplot(aes(predictor, mean)) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed') +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.25) +
  geom_point() + 
  geom_text(aes(y = ucl + 0.01, label = label), size = 4) +
  labs(x = NULL, y = 'effect size') + 
  theme_classic() +
  coord_flip()

# simulated vs. observed pgr & burrows:
simulations = list(
  'pgr' = MCMCvis::MCMCsummary(mod_final$results, c('ysim', 'ysim.2007'),
                     HPD = TRUE, 
                     pg0 = TRUE) %>% 
    mutate(WY = c(1999:2006, 2008:2021, 2007)) %>%
    left_join(modeldat %>% select(WY, pgr)) %>% 
    arrange(WY),
  'N' = MCMCvis::MCMCsummary(mod_final$results, c('Nsim'),
                             HPD = TRUE, 
                             pg0 = TRUE) %>% 
    mutate(WY = c(1998:2021)) %>%
    full_join(modeldat %>% select(WY, burrows)) %>% 
    arrange(WY))

plot_simulations(sim = simulations, maxyear = 2020, ylim2 = c(0, 28))
ggsave(filename = 'fig/simulations.png', width = 6, height = 6, units = 'in')

## effects plots------
pred = bind_rows(
  generate_predictions(mod_final$results, inputdat = inputdat2, obsdat = modeldat),
  generate_predictions(mod_final$results, inputdat = inputdat2, 
                       density = TRUE, obsdat = modeldat))
  
flowlevels = c('flow_concurrent', 'flow_prior', 'flow_older2', 
               'drought_prior', 'burrows_prior')
flowlabels = c('A. Concurrent flow (t)', 
               'B. Prior year flow (t-1)', 
               'C. Cumulative prior flow (t-2 + t-3)', 
               'D. Prior year drought (t-1)', 
               'E. Prior year burrow count (t-1)')

obsdat = modeldat %>% 
  filter(!is.na(pgr)) %>% 
  pivot_longer(burrows_prior:last_col(), names_to = 'predictor', 
               values_to = 'value') %>% 
  mutate(predictor = gsub('flowtotal$', 'flow_concurrent', predictor),
         predictor = gsub('flowtotal_', 'flow_', predictor)) %>% 
  filter(predictor %in% pred$predictor) %>% 
  #filter(grepl('flow', predictor)) %>% 
  mutate(predictor = factor(predictor, levels = flowlevels,
                            labels = flowlabels),
         # convert cfs to acre-feet
         value = if_else(grepl('flow', predictor), 
                         value * 1.9834591996927,
                         value))

betadat = effect_stats %>% select(predictor, pg0 = `p>0`) %>% 
  filter(!grepl('var|pvalue|trend|growth', predictor)) %>% 
  mutate(predictor = if_else(predictor == 'density_int', 'burrows_prior', predictor),
         predictor = factor(predictor, levels = flowlevels, labels = flowlabels),
         problab = paste0('P{\u03b2 > 0} = ', pg0),
         problab = if_else(pg0 >= 0.9 | pg0 <= 0.1, paste0(problab, '*'), problab))

pred %>% 
  mutate(predictor = factor(predictor, levels = flowlevels, labels = flowlabels),
         # convert cfs to acre-feet
         value = if_else(grepl('flow', predictor), 
                         value * 1.9834591996927,
                         value)) %>%
  ggplot(aes(value, mean)) + 
  facet_wrap(~predictor, ncol = 2) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray85') + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed', linewidth = 0.2) +
  geom_line(linewidth = 0.6) + 
  #scale_x_log10(limits = c(1, NA)) + 
  geom_point(data = obsdat, aes(value, pgr)) +
  geom_text(data = obsdat %>% 
              filter(WY %in% c(1999, 2000, 2018, 2017, 2019, 2010, 2015)), 
            aes(value, pgr, label = WY), 
            color = 'black', size = 3, hjust = 0, nudge_x = 0.2, 
            vjust = 'outward', check_overlap = FALSE) +
  geom_text(data = betadat, aes(label = problab, x = 0, y = Inf),
            hjust = 0, nudge_x = 0.2, vjust = 1, size = 3) +
  ylim(-0.6, 0.6) + 
  scale_x_continuous(limits = c(0, 33), expand = c(0, 0)) +
  #xlim(0, NA) +
  labs(x = 'acre-feet (millions)', y = 'per-capita growth rate (t-1 to t') +
  theme_classic() + 
  theme(strip.text = element_text(hjust = 0),
        strip.background = element_rect(color = NA))
ggsave(filename = 'fig/effects_flow.png',
       width = 6, height = 6, units = 'in')


pred2 %>% 
  ggplot(aes(value, mean)) + 
  #facet_wrap(~predictor, ncol = 1) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray85') + 
  geom_hline(aes(yintercept = 0), linetype = 'dashed', linewidth = 0.2) +
  geom_line(linewidth = 0.6) + 
  #scale_x_log10(limits = c(1, NA)) + 
  geom_point(data = obsdat, aes(value, pgr)) +
  geom_text(data = obsdat %>% 
              filter(WY %in% c(1999, 2000, 2018, 2017, 2019, 2010, 2015)), 
            aes(value, pgr, label = WY), 
            color = 'black', size = 3, hjust = 0, nudge_x = 0.2, 
            vjust = 'outward', check_overlap = FALSE) +
  geom_text(data = betadat, aes(label = problab, x = 0, y = Inf),
            hjust = 0, nudge_x = 0.2, vjust = 1, size = 3) +
  ylim(-0.6, 0.6) + 
  scale_x_continuous(limits = c(0, 33), expand = c(0, 0)) +
  #xlim(0, NA) +
  labs(x = 'acre-feet (millions)', y = 'per-capita growth rate (t-1 to t') +
  theme_classic() + 
  theme(strip.text = element_text(hjust = 0),
        strip.background = element_rect(color = NA))

# # examine variability in k estimates over the years:
# MCMCvis::MCMCplot(mod, params = 'k')
# inputdat2$Bpredictors %>%
#   mutate(k = MCMCvis::MCMCpstr(mod_growth_vark, 'k', func = median)$k) %>%
#   pivot_longer(-k) %>%
#   ggplot(aes(value, k)) + geom_point() + geom_smooth(method = 'gam') +
#   facet_wrap(~name, scales = 'free')
# # possibly stronger density dependence with highest flow14 values

## Summary stats & figure------------
cov_stats = modeldat %>% select(WY, currentflow, laggedflow, pdsi) %>%
  mutate(t = WY) %>% pivot_longer(-WY) %>%
  group_by(name) %>%
  summarize(mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            min = min(value, na.rm = TRUE),
            .groups = 'drop') %>% 
  write_csv('output/cov_stats.csv')

fig_covariates = ggsave(filename = 'fig/covariates.png',
                        plot = plot_covariates(modeldat),
                        width = 6, height = 6, units = 'in')


# ANALYSIS 2--------
# drop threshold concept and use total flow

modeldat2 = list(birddat %>% select(WY = year, burrows, pgr) %>% 
                  mutate(burrows_prior = lag(burrows)) %>% 
                  # add burrows_prior for 2021
                  bind_rows(
                    tibble(WY = 2021, burrows = NA, burrows_prior = 19035)),
                flowdat_all %>% 
                  select(WY, flowtotal, flowtotal_prior, flowtotal_older2),
                drought_sum %>% 
                  mutate(pdsi_prior = lag(pdsi_breeding)) %>% 
                  select(WY, pdsi_prior)) %>% 
  purrr::reduce(full_join) %>% 
  filter(WY >= 1999 & WY <= 2021) 

inputdat2 = list(
  # per-capita growth rate t-1 to t
  # >> exclude 2007 because burrows_prior (from 2006) is missing to estimate pgr
  y = modeldat2 %>% filter(WY != 2007) %>% pull(pgr),
  # density-independent predictors of growth rate from t-1 to t:
  # >> exclude 2007 because burrows_prior (from 2006) is missing to estimate pgr
  Bpredictors = modeldat2 %>% 
    select(WY, 
           flow_concurrent = flowtotal, 
           flow_prior = flowtotal_prior, 
           flow_older2 = flowtotal_older2,
           #flow_prior3 = flow_threshold_prior3,
           pdsi_prior) %>%
    #mutate(t = WY) %>% 
    # log transform flow data
    mutate_at(vars(starts_with('flow')), ~log(. + 1)) %>%
    # scale
    mutate_at(vars(flow_concurrent:pdsi_prior), 
              ~scale(., center = TRUE, scale = TRUE)) %>% 
    filter(WY != 2007) %>% select(-WY),
  # Kpredictors = modeldat2 %>% 
  #   select(t = WY) %>% 
  #   mutate(t = scale(t, center = TRUE, scale = TRUE)),
  # density-dependent reproduction driven by burrows_prior
  N = modeldat2 %>% filter(WY != 2007) %>% pull(burrows_prior) %>% 
    scale(center = TRUE, scale = TRUE) %>% as.vector(),
  # keep full list of unscaled N_obs for prediction purposes (including 2007)
  N_obs = modeldat2 %>% pull(burrows_prior) %>% round(digits = 0),
  # range of scaled values for each predictor
  predmatrix = tibble(flow_concurrent = seq(-1.6, 2.1, length.out = 11),
                      flow_prior = seq(-1.6, 2.1, length.out = 11),
                      #flow_prior3 = seq(-2, 1.35, length.out = 11),
                      flow_older2 = seq(-1.6, 2.1, length.out = 11),
                      pdsi_prior = seq(-1.3, 2.6, length.out = 11),
                      t = (seq(1999, 2020, 2) - attr(scale(modeldat$WY), 'scaled:center'))/attr(scale(modeldat$WY), 'scaled:scale'),
                      N = (seq(-1.6, 1.4, length.out = 11))))

mod2 = fit_test(
  inputdat2, n.adapt = 5000, n.burnin = 5000, n.sample = 100000, n.chains = 3,
  vars = c('r0', 'beta', 'k0', 'var.o', 'var.p', 
           'ysim', 'ysim.2007', 'Nsim', 'ypred', 'mu.ypred',
           'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  thin = 10)

MCMCvis::MCMCtrace(mod2$results,
                   params = c('r0', 'beta', 'k0', 'var.o', 'var.p'),
                   pdf = FALSE)
MCMCvis::MCMCsummary(mod2$results,
                     c('r0', 'beta', 'k0', 'var.o', 'var.p',
                       'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
                     HPD = TRUE, 
                     pg0 = TRUE) %>% 
  as_tibble(rownames = 'var') %>% 
  mutate(predictor = c('growth', names(inputdat2$Bpredictors), 
                       'density_dep', 
                       'var.obs', 'var.process', 
                       'pvalue.mean', 'pvalue.sd', 'pvalue.fit')) %>% 
  select(var, predictor, everything())
# fits well; n.eff sufficient; trace plots well-mixed;

dic2 = dic.samples(model = mod2$jm, n.iter = 100000, thin = 10, type = 'pD')

# effect_stats2 = get_sumstats(mod = mod2, 
#                             param = c('beta', 'k0', 'var.k', 'r0'),
#                             id = c(inputdat2$Bpredictors %>% names(), 
#                                    'k0', 'var.k', 'r0'))
# # high probability of positive effect of prior flow and older 2 flow, negative effect of k0
# write_csv('output/effect_stats2.csv')


# simulated vs. observed pgr & burrows:
simulations2 = extract_simulations(mod2$results, obsdat = modeldat2)
plot_simulations(sim = simulations2, maxyear = 2020)
ggsave(filename = 'fig/simulations2.png', width = 6, height = 6, units = 'in')

# strength of effects
plot_sumstats(mod2$results, param = c('beta', 'k0'),
              id = c(inputdat2$Bpredictors %>% names(), 'k0')) 
ggsave(filename = 'fig/effect_size2.png', width = 6, height = 3, units = 'in')

# effects plots
pred2 = generate_predictions(mod2$results, inputdat = inputdat2, obsdat = modeldat2)
plot_cov_relationships(
  predvalues = pred2, 
  obsdat = modeldat2 %>% 
    filter(!is.na(pgr)) %>% 
    select(N = burrows_prior,
           pgr, 
           flow_concurrent = flowtotal,
           flow_prior = flowtotal_prior,
           flow_older2 = flowtotal_older2,
           #N = burrows_prior, 
           pdsi_prior,
           #t = WY
           ) %>% 
    mutate(across(starts_with('flow'), ~./1000000)),
  ylim = c(-0.7, 0.8))
ggsave(filename = 'fig/covariate_effects2.png',
       width = 6, height = 6, units = 'in')

# ANALYSIS 3--------
# use ndays above threshold instead of flow data

modeldat3 = list(birddat %>% select(WY = year, burrows, pgr) %>% 
                   mutate(burrows_prior = lag(burrows)) %>% 
                   # add burrows_prior for 2021
                   bind_rows(
                     tibble(WY = 2021, burrows = NA, burrows_prior = 19035)),
                 flowdat_all %>% 
                   select(WY, ndays_threshold, ndays_threshold_prior, ndays_threshold_older2),
                 drought_sum %>% 
                   mutate(pdsi_prior = lag(pdsi_breeding)) %>% 
                   select(WY, pdsi_prior)) %>% 
  purrr::reduce(full_join) %>% 
  filter(WY >= 1999 & WY <= 2021) 

inputdat3 = list(
  # per-capita growth rate t-1 to t
  # >> exclude 2007 because burrows_prior (from 2006) is missing to estimate pgr
  y = modeldat3 %>% filter(WY != 2007) %>% pull(pgr),
  # density-independent predictors of growth rate from t-1 to t:
  # >> exclude 2007 because burrows_prior (from 2006) is missing to estimate pgr
  Bpredictors = modeldat3 %>% 
    select(WY, 
           flow_concurrent = ndays_threshold, 
           flow_prior = ndays_threshold_prior, 
           flow_older2 = ndays_threshold_older2,
           #flow_prior3 = flow_threshold_prior3,
           pdsi_prior) %>%
    mutate(t = WY) %>% 
    # scale
    mutate_at(vars(flow_concurrent:t), ~scale(., center = TRUE, scale = TRUE)) %>% 
    filter(WY != 2007) %>% select(-WY),
  # density-dependent reproduction driven by burrows_prior
  N = modeldat2 %>% filter(WY != 2007) %>% pull(burrows_prior) %>% 
    scale(center = TRUE, scale = TRUE) %>% as.vector(),
  # keep full list of unscaled N_obs for prediction purposes (including 2007)
  N_obs = modeldat2 %>% pull(burrows_prior) %>% round(digits = 0),
  # range of scaled values for each predictor
  predmatrix = tibble(flow_concurrent = seq(-1, 2, length.out = 11),
                      flow_prior = seq(-1, 2.7, length.out = 11),
                      #flow_prior3 = seq(-2, 1.35, length.out = 11),
                      flow_older2 = seq(-1.2, 2.2, length.out = 11),
                      pdsi_prior = seq(-1.25, 2.5, length.out = 11),
                      t = (seq(1999, 2020, 2) - attr(scale(modeldat$WY), 'scaled:center'))/attr(scale(modeldat$WY), 'scaled:scale'),
                      N = (seq(-1.6, 1.4, length.out = 11))))

mod3 = fit_test(
  inputdat, n.adapt = 5000, n.burnin = 5000, n.sample = 100000, n.chains = 3,
  vars = c('r0', 'beta', 'k0', 'k', 'var.o', 'var.p', 'var.k',
           'ysim', 'ysim.2007', 'Nsim', 'ypred', 'mu.ypred',
           'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  thin = 10)

MCMCvis::MCMCtrace(mod3$results,
                   params = c('r0', 'beta', 'k0', 'var.o', 'var.p', 'var.k'),
                   pdf = FALSE)
MCMCvis::MCMCsummary(mod3$results,
                     c('r0', 'beta', 'k0', 'var.o', 'var.p', 'var.k',
                       'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
                     HPD = TRUE, 
                     pg0 = TRUE) %>% 
  as_tibble(rownames = 'var') %>% 
  mutate(predictor = c('growth', names(inputdat$Bpredictors), 'densitydep', 
                       'var.obs', 'var.process', 'var.densitydep', 
                       'pvalue.mean', 'pvalue.sd', 'pvalue.fit')) %>% 
  select(var, predictor, everything())
# fits well; n.eff sufficient; trace plots well-mixed;

dic3 = dic.samples(model = mod3$jm, n.iter = 100000, thin = 10, type = 'pD')



# effect_stats3 = get_sumstats(mod = mod3, 
#                              param = c('beta', 'k0', 'var.k', 'r0'),
#                              id = c(inputdat3$Bpredictors %>% names(), 
#                                     'k0', 'var.k', 'r0'))
# # high probability of positive effect of prior flow and older 2 flow, negative effect of k0
# write_csv('output/effect_stats3.csv')


# simulated vs. observed pgr & burrows:
simulations3 = extract_simulations(mod3, obsdat = modeldat3)
plot_simulations(sim = simulations3, maxyear = 2020, ylim2 = c(0, 26))
ggsave(filename = 'fig/simulations3.png', width = 6, height = 6, units = 'in')

# strength of effects
plot_sumstats(mod3, param = c('beta', 'k0', 'var.k'),
              id = c(inputdat3$Bpredictors %>% names(), 'k0', 'var.k')) + 
  facet_wrap(~id, scales = 'free')
ggsave(filename = 'fig/effect_size3.png', width = 6, height = 3, units = 'in')

# effects plots
pred3 = generate_predictions(mod3, inputdat = inputdat3, obsdat = modeldat3, logflow = FALSE)
plot_cov_relationships(
  predvalues = pred3, 
  obsdat = modeldat3 %>% 
    filter(!is.na(pgr)) %>% 
    select(pgr, 
           flow_concurrent = ndays_threshold, 
           flow_prior = ndays_threshold_prior, 
           flow_older2 = ndays_threshold_older2,
           N = burrows_prior, 
           pdsi_prior,
           t = WY))
ggsave(filename = 'fig/covariate_effects2.png',
       width = 6, height = 6, units = 'in')

# ANALYSIS 4--------
# drop threshold concept and use total flow

modeldat4 = list(birddat %>% select(WY = year, burrows, pgr) %>% 
                   mutate(burrows_prior = lag(burrows)) %>% 
                   # add burrows_prior for 2021
                   bind_rows(
                     tibble(WY = 2021, burrows = NA, burrows_prior = 19035)),
                 flowdat_all %>% 
                   select(WY, flowtotal, flowtotal_prior, flowtotal_older2),
                 drought_sum %>% 
                   mutate(pdsi_prior = lag(pdsi_breeding)) %>% 
                   select(WY, pdsi_prior)) %>% 
  purrr::reduce(full_join) %>% 
  filter(WY >= 1999 & WY <= 2021) 

mean(modeldat4$pgr, na.rm = TRUE) # mean = 0.041
var(modeldat4$pgr, na.rm = TRUE) # var = 0.057


inputdat2 = list(
  # per-capita growth rate t-1 to t
  # >> exclude 2007 because burrows_prior (from 2006) is missing to estimate pgr
  y = modeldat2 %>% filter(WY != 2007) %>% pull(pgr),
  # density-independent predictors of growth rate from t-1 to t:
  # >> exclude 2007 because burrows_prior (from 2006) is missing to estimate pgr
  Bpredictors = modeldat2 %>% 
    select(WY, 
           flow_concurrent = flowtotal, 
           flow_prior = flowtotal_prior, 
           flow_older2 = flowtotal_older2,
           #flow_prior3 = flow_threshold_prior3,
           pdsi_prior) %>%
    mutate(t = WY) %>% 
    # log transform flow data
    mutate_at(vars(starts_with('flow')), ~log(. + 1)) %>%
    # scale
    mutate_at(vars(flow_concurrent:t), ~scale(., center = TRUE, scale = TRUE)) %>% 
    filter(WY != 2007) %>% select(-WY),
  # density-dependent reproduction driven by burrows_prior
  N = modeldat2 %>% filter(WY != 2007) %>% pull(burrows_prior) %>% 
    scale(center = TRUE, scale = TRUE) %>% as.vector(),
  # keep full list of unscaled N_obs for prediction purposes (including 2007)
  N_obs = modeldat2 %>% pull(burrows_prior) %>% round(digits = 0),
  # range of scaled values for each predictor
  predmatrix = tibble(flow_concurrent = seq(-1.6, 2, length.out = 11),
                      flow_prior = seq(-1.6, 2, length.out = 11),
                      #flow_prior3 = seq(-2, 1.35, length.out = 11),
                      flow_older2 = seq(-1.8, 1.8, length.out = 11),
                      pdsi_prior = seq(-1.25, 2.5, length.out = 11),
                      t = (seq(1999, 2020, 2) - attr(scale(modeldat$WY), 'scaled:center'))/attr(scale(modeldat$WY), 'scaled:scale'),
                      N = (seq(-1.6, 1.4, length.out = 11))))

mod2 = fit_test(
  inputdat, n.adapt = 5000, n.burnin = 5000, n.sample = 100000, n.chains = 3,
  vars = c('r0', 'beta', 'k0', 'k', 'var.o', 'var.p', 'var.k',
           'ysim', 'ysim.2007', 'Nsim', 'ypred', 'mu.ypred',
           'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  thin = 10)

MCMCvis::MCMCtrace(mod2$results,
                   params = c('r0', 'beta', 'k0', 'var.o', 'var.p', 'var.k'),
                   pdf = FALSE)
MCMCvis::MCMCsummary(mod2$results,
                     c('r0', 'beta', 'k0', 'var.o', 'var.p', 'var.k',
                       'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
                     HPD = TRUE, 
                     pg0 = TRUE) %>% 
  as_tibble(rownames = 'var') %>% 
  mutate(predictor = c('growth', names(inputdat$Bpredictors), 'densitydep', 
                       'var.obs', 'var.process', 'var.densitydep', 
                       'pvalue.mean', 'pvalue.sd', 'pvalue.fit')) %>% 
  select(var, predictor, everything())
# fits well; n.eff sufficient; trace plots well-mixed;

dic2 = dic.samples(model = mod2$jm, n.iter = 100000, thin = 10, type = 'pD')

# effect_stats2 = get_sumstats(mod = mod2, 
#                             param = c('beta', 'k0', 'var.k', 'r0'),
#                             id = c(inputdat2$Bpredictors %>% names(), 
#                                    'k0', 'var.k', 'r0'))
# # high probability of positive effect of prior flow and older 2 flow, negative effect of k0
# write_csv('output/effect_stats2.csv')


# simulated vs. observed pgr & burrows:
simulations2 = extract_simulations(mod2$results, obsdat = modeldat2)
plot_simulations(sim = simulations2, maxyear = 2021, ylim1 = c(-1.2,1.2), ylim2 = c(-5, 40))
ggsave(filename = 'fig/simulations2.png', width = 6, height = 6, units = 'in')

# strength of effects
plot_sumstats(mod2, param = c('beta', 'k0', 'var.k'),
              id = c(inputdat2$Bpredictors %>% names(), 'k0', 'var.k')) + 
  facet_wrap(~id, scales = 'free')
ggsave(filename = 'fig/effect_size2.png', width = 6, height = 3, units = 'in')

# effects plots
pred2 = generate_predictions(mod2, inputdat = inputdat2, obsdat = modeldat2)
plot_cov_relationships(
  predvalues = pred2, 
  obsdat = modeldat2 %>% 
    filter(!is.na(pgr)) %>% 
    select(pgr, 
           flow_concurrent = flowtotal,
           flow_prior = flowtotal_prior,
           flow_older2 = flowtotal_older2,
           N = burrows_prior, 
           pdsi_prior,
           t = WY) %>% 
    mutate(across(starts_with('flow'), ~./1000000)))
ggsave(filename = 'fig/covariate_effects2.png',
       width = 6, height = 6, units = 'in')

# OLD STUFF-------
# # Compiled predictors----------------
# # check for correlations among predictors
# 
# modeldat = list(birddat %>% 
#                   select(WY = year, burrows, pgr),
#                 # %>% 
#                 #   mutate(burrows_tp1 = lead(burrows)),
#                 # line up current water year's flow data (starting the previous October) with the same water year's #burrows
#                 mflow_sum %>% 
#                   select(WY, currentflow = flowtotal, laggedflow = flow2_lag1),
#                   # # add leading effects (from year t+1, just before the next
#                   # # breeding season)
#                   # mutate(flow14_tp1 = lead(flow14),
#                   #        flow14_3_tp1 = lead(flow14_3)),
#                 drought_sum %>% 
#                   select(WY, pdsi = pdsi_breeding) 
#                 # %>% 
#                 #   mutate(pdsi_tp1 = lead(pdsi))
#                 ) %>% 
#   purrr::reduce(full_join) %>% 
#   filter(WY >= 1999) %>% 
#   write_csv('data/modeldat.csv')
# # missing only _tp1 data for 2020
#   
# # check correlations with birddat & across metrics (log-transform flow metrics)
# corrplot::corrplot.mixed(
#   cor(modeldat %>% mutate_at(vars(ends_with('flow')), ~log(.+1)) %>%
#         drop_na(),
#       method = 'pearson'),
#   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# # strong positive correlations:
# # - laggedflow & burrows_tp1 (0.66)
# # - pgr & currentflow (0.69)
# # strongest negative correlation:
# # - burrows & WY (-0.63)
#   
# # Summary stats & figure------------
# cov_stats = modeldat %>% select(WY, currentflow, laggedflow, pdsi) %>%
#   mutate(t = WY) %>% pivot_longer(-WY) %>%
#   group_by(name) %>%
#   summarize(mean = mean(value, na.rm = TRUE),
#             median = median(value, na.rm = TRUE),
#             sd = sd(value, na.rm = TRUE),
#             max = max(value, na.rm = TRUE),
#             min = min(value, na.rm = TRUE),
#             .groups = 'drop') %>% 
#   write_csv('output/cov_stats.csv')
# 
# fig_covariates = ggsave(filename = 'fig/covariates.png',
#                         plot = plot_covariates(modeldat),
#                         width = 6, height = 6, units = 'in')
# 
# # ANALYSIS-----------
# 
# # input data: drop 2006 because can't having missing values in
# # N_obs since it is a predictor
# inputdat = list(
#   # per-capita growth rate t to t+1
#   y = modeldat %>% filter(WY != 2006) %>% pull(pgr),
#   # density-independent predictors of growth rate from t to t+1:
#   Bpredictors = modeldat %>% 
#     select(WY, currentflow, laggedflow, pdsi) %>%
#     mutate(t = WY) %>% 
#     # log transform flow data
#     mutate_at(vars(ends_with('flow')), ~log(.+1)) %>%
#     # scale
#     mutate_at(vars(currentflow:t), ~scale(., center = TRUE, scale = TRUE)) %>% 
#     filter(WY != 2006) %>% select(-WY),
#   N = modeldat %>% filter(WY != 2006) %>% 
#     pull(burrows) %>% scale(center = TRUE, scale = TRUE) %>% as.vector(),
#   # keep full list of unscaled N_obs for prediction purposes
#   N_obs = modeldat %>% pull(burrows) %>% round(digits = 0),
#   # range of values for each predictor
#   predmatrix = tibble(currentflow = seq(-2.15, 1.45, length.out = 11),
#                       laggedflow = seq(-1.75, 1.45, length.out = 11),
#                       pdsi = seq(-1.35, 2.05, length.out = 11),
#                       t = (seq(1999, 2020, 2) - attr(scale(modeldat$WY), 'scaled:center'))/attr(scale(modeldat$WY), 'scaled:scale'),
#                       N = (seq(-1.7, 1.4, length.out = 11)))
# )
# 
# # fit logistic model: allows strength of density dependence to vary annually
# # from a normal distribution with mean k0 and error var.k
# mod = fit_BANSpopgrowth_model(
#   inputdat, n.adapt = 5000, n.burnin = 15000, n.sample = 50000, n.chains = 3,
#   vars = c('r0', 'beta', 'k0', 'k', 'var.o', 'var.p', 'var.k',
#            'ysim', 'ysim.2006', 'Nsim', 'ypred', 'mu.ypred',
#            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
#   thin = 25)
# save(mod, file = 'output/model.RData')
# 
# # MCMCvis::MCMCsummary(mod,
# #                      c('r0', 'beta', 'k0', 'var.o', 'var.p', 'var.k',
# #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
# # MCMCvis::MCMCtrace(mod,
# #                    params = c('r0', 'beta', 'k0', 'var.o', 'var.p', 'var.k'),
# #                    pdf = FALSE)
# # fits well; n.eff sufficient; trace plots well-mixed;
#   
# # Summary stats & figures---------
# effect_stats = get_sumstats(mod, param = c('beta', 'k0', 'var.k', 'r0'),
#                             id = c(inputdat$Bpredictors %>% names(), 
#                                    'k0', 'var.k', 'r0')) %>% 
#   write_csv('output/effect_stats.csv')
# # high probability of positive effect of flow14_3, negative effect of k0
#   
# # simulated vs. observed pgr & burrows:
# fig_simulations = ggsave(filename = 'fig/simulations.png',
#                          plot = plot_simulations(mod, obsdat = modeldat),
#                          width = 6, height = 6, units = 'in')
# 
# # strength of effects
# fig_effectsizes = ggsave(filename = 'fig/effect_size.png',
#                          plot = plot_sumstats(mod, c('beta', 'k0', 'var.k'),
#                                               id = c(inputdat$Bpredictors %>% 
#                                                        names(), 'k0', 'var.k')),
#                          width = 6, height = 3, units = 'in')
#   
# fig_covpred = ggsave(filename = 'fig/covariate_effects.png',
#                      plot = plot_cov_relationships(mod, inputdat = inputdat,
#                                                    obsdat = modeldat),
#                      width = 6, height = 6, units = 'in')
#   
#   # # examine variability in k estimates over the years:
#   # MCMCvis::MCMCplot(mod, params = 'k')
#   # inputdat2$Bpredictors %>%
#   #   mutate(k = MCMCvis::MCMCpstr(mod_growth_vark, 'k', func = median)$k) %>%
#   #   pivot_longer(-k) %>%
#   #   ggplot(aes(value, k)) + geom_point() + geom_smooth(method = 'gam') +
#   #   facet_wrap(~name, scales = 'free')
#   # # possibly stronger density dependence with highest flow14 values
# 
