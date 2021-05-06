# load packages and functions
source('R/packages.R')
source('R/functions.R')
source('R/fit_model.R')

# update readme for github repo https://github.com/pointblue/bank_swallow
readme_page = render_Rmd("Rmd/README.Rmd", "README.md")
  
# RESPONSE DATA--------
# BANS survey data: annual burrow counts between river miles 144-243 (roughly
# Colusa to Red Bluff); counts are the average of two separate counters'
# totals, so sometimes not an integer; survey protocols changed ~1999, so
# analysis focuses on data from 1999 onward

birddat = read_csv('data/BANS_burrow_counts.csv', col_types = cols()) %>% 
  # check for / fill in any missing years
  complete(year = seq(min(year), max(year))) %>%
  # calculate annual per-capita growth rate as the change in burrow count to
  # the following year
  mutate(agr = lead(burrows)/burrows,
         pgr = (lead(burrows) - burrows)/burrows) %>% 
  write_csv('data/birddat.csv')
# Note: no count data from 2006, so no agr/pgr calculated for 2005-06
# interval, or 2006-07 interval

# COVARIATES---------
# Focus on effects of environmental conditions on reproductive success, and
# subsequent growth to the following year; consider also density-dependence

# 1. Nest habitat quantity & quality-----------
# total daily mean flow above a threshold is proportional to stream power and
# extent of erosion, increasing the quantity and quality of nest habitat;
# thresholds considered: 14000 cfs and 50000 cfs; 14-15,000 cited as a
# threshold above which erosion occurs; 50,000 cfs cited as causing widespread
# erosion throughout the system (Larsen et al. 2006; BANS-TAC 2013)

# > manually download mean daily flow data for each gauge from:
# https://cdec.water.ca.gov/dynamicapp/wsSensorData

# relevant gauges (N to S): 
# - Vina Bridge (VIN): data available from 1993-01-01
# - Hamilton City (HMC): data available from 2001-01-01
# - Ord Ferry (ORD): data available from 1993-01-01
# - Butte City Gauge (BTC): data available from 1993-01-01

mflowdat = compile_waterdat(dir = 'data/CDEC', type = 'MFLOW',
                            col_types = cols()) %>% 
  cleanup_waterdat(mindays = 365) %>% # filter out incomplete annual data
  calculate_annual_total(threshold = 14000, #water-year totals above threshold
                         mindays = 0.9 * 365) %>% 
  write_csv('data/mflowdat.csv')
# Note: originally explored seasonal flow totals, but annual totals highly
# correlated with rainy season totals, and breeding season totals very low

# check correlations across stations & metrics:
# - each metric highly correlated across 4 stations (>0.99)
# - strong correlations across metrics also (>0.85)
corrplot::corrplot.mixed(
  cor(mflowdat %>% filter(WY >= 1999) %>% select(-ndays) %>%
        pivot_wider(names_from = STATION_ID,
                    values_from = flowtotal) %>%
        select(-WY) %>% drop_na(),
      method = 'pearson'),
  tl.col = 'black', tl.pos = 'lt')

# summarize average across stations, and calculate rolling total
mflow_sum = mflowdat %>% 
  group_by(WY) %>% 
  summarize(nstations = length(flowtotal[!is.na(flowtotal)]), #note sometimes only 3 stations available per time step
            flowtotal = mean(flowtotal, na.rm = TRUE),
            .groups = 'drop') %>% 
  select(-nstations) %>% #ranges 3-4 stations contributing data each water year
  # calculate running totals for current and 2 prior water years
  mutate(flowtotal_3 = zoo::rollsum(flowtotal, 3, fill = NA, align = 'right')) %>% 
  # filter after computing lag data
  filter(WY >= 1999)

# check correlations:
corrplot::corrplot.mixed(
  cor(mflow_sum %>% drop_na(), method = 'pearson'),
    order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# no strong correlations among flow metrics or with year

# 2. Breeding season conditions---------
# loss of foraging habitat for BANS (and more broadly, aerial insectivores) is
# likely an important factor in pop declines; drought/primary productivity may
# play an important role in nestling survival

# drought indices from NOAA for Sacramento Valley area:

# code to pull straight from NCDC FTP website:
# droughtdat = purrr::map_df(c('pdsidv', 'zndxdv'),
#                            ~get_drought_indices(datname = .x)) %>%
#   cleanup_droughtdat(statecode = '04', climdivcode = '02')

# alt: pull from local files
droughtdat = bind_rows(
  read_table('data/climdiv-pdsidv-v1.0.0-20210406',
           col_names = c('ID', 'JAN', 'FEB', 'MAR', 'APR', 
                         'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 
                         'OCT', 'NOV', 'DEC')) %>% 
    mutate(index = 'pdsi'),
  read_table('data/climdiv-zndxdv-v1.0.0-20210406',
             col_names = c('ID', 'JAN', 'FEB', 'MAR', 'APR', 
                           'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 
                           'OCT', 'NOV', 'DEC')) %>% 
    mutate(index = 'zndx')
  ) %>% 
  cleanup_droughtdat(statecode = '04', climdivcode = '02') %>% 
  filter(WY > 1995) %>% 
  write_csv('data/droughtdat.csv')

# generate breeding season summaries:
drought_sum = droughtdat %>% 
  # filter(season == 'breeding') %>% 
  filter(!is.na(season)) %>% 
  group_by(index, WY, season) %>%
  summarize(value = mean(value),
            .groups = 'drop') %>% 
  pivot_wider(names_from = c('index', 'season')) %>% 
  filter(WY >= 1999 & WY < 2021)

# correlations:
corrplot::corrplot.mixed(
  cor(drought_sum, method = 'pearson'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# zndx and pdsi strongly correlated within season (breeding = 0.92; rainy = 0.85);
# breeding and rainy season metrics not strongly correlated (zndx = 0.43; pdsi = 0.64)

corrplot::corrplot.mixed(
  cor(left_join(mflow_sum, drought_sum, by = 'WY') %>% 
        left_join(birddat, by = c('WY' = 'year')) %>% drop_na(), method = 'pearson'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt'
)
# flowtotal correlated with rainy season indices (pdsi = 0.82; zndx = 0.87)
# pgr correlated with zndx_rainy (0.82), flowtotal (0.72), and pdsi_rainy (0.67)

# Compiled predictors----------------
# check for correlations among predictors

modeldat = list(birddat %>% 
                  select(year, burrows, pgr) %>% 
                  mutate(burrows_tp1 = lead(burrows)),
                # keep flow14 data; line up with current water year's flow
                # data (starting the previous October)
                mflow_sum %>% 
                  select(year = WY, flow14 = flowtotal, flow14_3 = flowtotal_3) %>% 
                  # add leading effects (from year t+1, just before the next
                  # breeding season)
                  mutate(flow14_tp1 = lead(flow14),
                         flow14_3_tp1 = lead(flow14_3)),
                drought_sum %>% 
                  select(year = WY, pdsi = pdsi_breeding) %>% 
                  mutate(pdsi_tp1 = lead(pdsi))
                ) %>% 
  purrr::reduce(full_join) %>% 
  filter(year >= 1999) %>% 
  write_csv('data/modeldat.csv')
# missing only _tp1 data for 2020
  
# check correlations with birddat & across metrics (log-transform flow metrics)
corrplot::corrplot.mixed(
  cor(modeldat %>% mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
        drop_na(),
      method = 'pearson'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# strong positive correlations:
# - flow14_3 & burrows_tp1 (0.87)
# - pgr & flow14 (0.69)
# strongest negative correlation:
# - burrows & year (-0.63)
  
# Summary stats & figure------------
cov_stats = modeldat %>% select(year, flow14, flow14_3, pdsi) %>%
  mutate(t = year) %>% pivot_longer(-year) %>%
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

# ANALYSIS-----------

# input data: drop 2006 because can't having missing values in
# N_obs since it is a predictor
inputdat = list(
  # per-capita growth rate t to t+1
  y = modeldat %>% filter(year != 2006) %>% pull(pgr),
  # density-independent predictors of growth rate from t to t+1:
  Bpredictors = modeldat %>% 
    select(year, flow14, flow14_3, pdsi) %>%
    mutate(t = year) %>% 
    # log transform flow data
    mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
    # scale
    mutate_at(vars(flow14:t), ~scale(., center = TRUE, scale = TRUE)) %>% 
    filter(year != 2006) %>% select(-year),
  N = modeldat %>% filter(year != 2006) %>% 
    pull(burrows) %>% scale(center = TRUE, scale = TRUE) %>% as.vector(),
  # keep full list of unscaled N_obs for prediction purposes
  N_obs = modeldat %>% pull(burrows) %>% round(digits = 0),
  # range of values for each predictor
  predmatrix = tibble(flow14 = seq(-2.15, 1.45, length.out = 11),
                      flow14_3 = seq(-2, 1.4, length.out = 11),
                      pdsi = seq(-1.35, 2.05, length.out = 11),
                      t = (seq(1999, 2020, 2) - 2009.5)/6.493587,
                      N = (seq(-1.7, 1.4, length.out = 11)))
)

# fit logistic model: allows strength of density dependence to vary annually
# from a normal distribution with mean k0 and error var.k
mod = fit_BANSpopgrowth_model(
  inputdat, n.adapt = 5000, n.burnin = 15000, n.sample = 50000, n.chains = 3,
  vars = c('r0', 'beta', 'k0', 'k', 'var.o', 'var.p', 'var.k',
           'ysim', 'ysim.2006', 'Nsim', 'ypred', 'mu.ypred',
           'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  thin = 25)
save(mod, file = 'output/model.RData')

# MCMCvis::MCMCsummary(mod,
#                      c('r0', 'beta', 'k0', 'var.o', 'var.p', 'var.k',
#                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
# MCMCvis::MCMCtrace(mod,
#                    params = c('r0', 'beta', 'k0', 'var.o', 'var.p', 'var.k'),
#                    pdf = FALSE)
# fits well; n.eff sufficient; trace plots well-mixed;
  
# Summary stats & figures---------
effect_stats = get_sumstats(mod, param = c('beta', 'k0', 'var.k', 'r0'),
                            id = c(inputdat$Bpredictors %>% names(), 
                                   'k0', 'var.k', 'r0')) %>% 
  write_csv('output/effect_stats.csv')
# high probability of positive effect of flow14_3, negative effect of k0
  
# simulated vs. observed pgr & burrows:
fig_simulations = ggsave(filename = 'fig/simulations.png',
                         plot = plot_simulations(mod, obsdat = modeldat),
                         width = 6, height = 6, units = 'in')

# strength of effects
fig_effectsizes = ggsave(filename = 'fig/effect_size.png',
                         plot = plot_sumstats(mod, c('beta', 'k0', 'var.k'),
                                              id = c(inputdat$Bpredictors %>% 
                                                       names(), 'k0', 'var.k')),
                         width = 6, height = 3, units = 'in')
  
fig_covpred = ggsave(filename = 'fig/covariate_effects.png',
                     plot = plot_cov_relationships(mod, inputdat = inputdat,
                                                   obsdat = modeldat),
                     width = 6, height = 6, units = 'in')
  
  # # examine variability in k estimates over the years:
  # MCMCvis::MCMCplot(mod, params = 'k')
  # inputdat2$Bpredictors %>%
  #   mutate(k = MCMCvis::MCMCpstr(mod_growth_vark, 'k', func = median)$k) %>%
  #   pivot_longer(-k) %>%
  #   ggplot(aes(value, k)) + geom_point() + geom_smooth(method = 'gam') +
  #   facet_wrap(~name, scales = 'free')
  # # possibly stronger density dependence with highest flow14 values

