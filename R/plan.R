# Drake plan reference: https://books.ropensci.org/drake/plans.html

plan <- drake_plan(
  # update readme for github repo https://github.com/pointblue/bank_swallow
  readme_page = render_Rmd(file_in("Rmd/README.Rmd"),
                           file_out("README.md")),
  
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
           pgr = (lead(burrows) - burrows)/burrows),
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

  mflowdat = compile_waterdat(dir = file_in('data/CDEC'), type = 'MFLOW',
                              col_types = cols()) %>% 
    select(STATION_ID, `DATE TIME`, VALUE) %>% 
    mutate(VALUE = if_else(VALUE == '---' | VALUE < 0, NA_character_, VALUE),
           VALUE = as.numeric(VALUE),
           month = format(`DATE TIME`, '%m') %>% as.numeric(),
           year = format(`DATE TIME`, '%Y') %>% as.numeric(),
           WY = if_else(month < 10, year, year + 1)) %>% 
    # filter out incomplete annual data
    group_by(STATION_ID, WY) %>% 
    add_count() %>% 
    ungroup() %>% 
    filter(n >= 365) %>% 
    # subtract threshold values
    mutate(flow50 = if_else(VALUE > 50000, VALUE - 50000, 0), #daily flow above 50k cfs
           flow14 = if_else(VALUE > 14000, VALUE - 14000, 0)) %>%  #daily flow above 14k cfs
    # summarize annual flow totals above thresholds per station
    group_by(STATION_ID, WY) %>% 
    summarize(ndays = length(VALUE[!is.na(VALUE)]),
              flowtot = sum(VALUE, na.rm = TRUE),
              flow50 = sum(flow50, na.rm = TRUE),
              flow14 = sum(flow14, na.rm = TRUE),
              .groups = 'drop') %>% 
    # check those with data for <90% of days each year
    mutate_at(vars(flowtot:flow14),
              ~if_else(ndays < 0.9*365, NA_real_, .)),
  
  # Note: originally explored seasonal flow data, but annual totals highly
  # correlated with rainy season totals, and breeding season totals very low

  # # check correlations across stations & metrics:
  # # - each metric highly correlated across 4 stations (>0.95)
  # # - strong correlations across metrics also (>0.85)
  # corrplot::corrplot.mixed(
  #   cor(mflowdat %>% filter(WY >= 1999) %>% select(-ndays) %>% 
  #         pivot_wider(names_from = STATION_ID,
  #                     values_from = flowtot:flow14) %>%
  #         select(-WY) %>% drop_na(),
  #       method = 'pearson'),
  #   tl.col = 'black', tl.pos = 'lt')
  
  # ggplot(mflowdat, aes(WY, flowtot)) + geom_line() + geom_point() +
  #   geom_line(aes(WY, flow50), color = 'blue') +
  #   geom_line(aes(WY, flow14), color = 'red') +
  #   facet_wrap(~STATION_ID)

  # summarize average across stations
  mflow_sum = mflowdat %>% 
    group_by(WY) %>% 
    summarize(nstations = length(flowtot[!is.na(flowtot)]), #note sometimes only 3 stations available per time step
              flowtot = mean(flowtot, na.rm = TRUE),
              flow50 = mean(flow50, na.rm = TRUE),
              flow14 = mean(flow14, na.rm = TRUE),
              .groups = 'drop') %>% 
    select(-nstations) %>% #ranges 3-4 stations contributing data each water year
    # calculate running totals for current and 2 prior water years
    mutate(flowtot_3 = zoo::rollsum(flowtot, 3, fill = NA, align = 'right'),
           flow50_3 = zoo::rollsum(flow50, 3, fill = NA, align = 'right'),
           flow14_3 = zoo::rollsum(flow14, 3, fill = NA, align = 'right')) %>% 
    # filter after computing lag data
    filter(WY >= 1999),
  
  # check correlations:
  # corrplot::corrplot.mixed(
  #   cor(mflow_sum %>% drop_na(), method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # # - strong correlations among annual and cumulative metrics, but not across
  # # them; no strong trend/correlation with water year

  # 2. Breeding season conditions---------
  # loss of foraging habitat for BANS (and more broadly, aerial insectivores) is
  # likely an important factor in pop declines; drought/primary productivity may
  # play an important role in nestling survival

  # drought indices from NOAA for Sacramento Valley area:
  droughtdat = purrr::map_df(c('pdsidv', 'zndxdv'),
                          ~get_drought_indices(datname = .x, statecode = '04', 
                                               climdivcode = '02')) %>% 
    mutate(month = format(date, '%m') %>% as.numeric(),
           year = as.numeric(year),
           WY = if_else(month >= 10, year + 1, year),
           season = case_when(month >= 10 | month <= 3 ~ 'rainy',
                              month >= 4 & month <= 8 ~ 'breeding',
                              TRUE ~ NA_character_)),
  
  # generate breeding season summaries:
  drought_sum = droughtdat %>% 
    filter(season == 'breeding') %>% 
    group_by(index, WY) %>%
    summarize(value = mean(value),
              .groups = 'drop') %>% 
    pivot_wider(names_from = 'index') %>% 
    filter(WY >= 1999 & WY < 2021),

  # # correlations:
  # corrplot::corrplot.mixed(
  #   cor(drought_sum %>% filter(WY >= 1999), method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # # strong correlation across indices (0.92); no trend/water year correlation

  # Compiled predictors----------------
  # check for correlations among predictors

  modeldat = list(birddat %>% select(year, burrows, pgr) %>% 
                    mutate(burrows_tp1 = lead(burrows)),
                  # keep flow14 data; line up with current water year's flow
                  # data (starting the previous October)
                  mflow_sum %>% 
                    select(year = WY, flow14, flow14_3) %>% 
                    # add leading effects (from year t+1, just before the next
                    # breeding season)
                    mutate(flow14_tp1 = lead(flow14),
                           flow14_3_tp1 = lead(flow14_3)),
                  drought_sum %>% select(year = WY, pdsi = pdsidv) %>% 
                    mutate(pdsi_tp1 = lead(pdsi))
                  ) %>% 
    purrr::reduce(full_join) %>% 
    filter(year >= 1999),
  # missing only _tp1 data for 2020
  
  # # check correlations with birddat & across metrics (log-transform flow metrics)
  # corrplot::corrplot.mixed(
  #   cor(modeldat %>% mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
  #         drop_na(),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # strong positive correlations:
  # - flow14_3 & burrows_tp1 (0.87)
  # - pgr & flow14 (0.69)
  #
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
              .groups = 'drop'),
  
  fig_covariates = ggsave(filename = 'fig/covariates.png',
                          plot = plot_covariates(modeldat),
                          width = 6, height = 6, units = 'in'),

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
    predmatrix = tibble(flow14 = seq(-2.2, 1.5, length.out = 11),
                        flow14_3 = seq(-2, 1.4, length.out = 11),
                        pdsi = seq(-1.4, 2, length.out = 11),
                        t = (seq(1999, 2020, 2) - 2009.5)/6.493587,
                        N = (seq(-1.6, 1.4, length.out = 11)))
    ),

  # fit logistic model: allows strength of density dependence to vary annually
  # from a normal distribution with mean k0 and error var.k
  mod = fit_BANSpopgrowth_model(
    inputdat, n.adapt = 5000, n.burnin = 15000, n.sample = 50000, n.chains = 3,
    vars = c('r0', 'beta', 'k0', 'k', 'var.o', 'var.p', 'var.k',
             'ysim', 'ysim.2006', 'Nsim', 'ypred',
             'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
    thin = 25),
  
  # MCMCvis::MCMCsummary(mod,
  #                      c('r0', 'beta', 'k0', 'var.o', 'var.p', 'var.k',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod,
  #                    params = c('r0', 'beta', 'k0', 'var.o', 'var.p', 'var.k'),
  #                    pdf = FALSE)
  # fits well; n.eff sufficient; trace plots well-mixed;
  
  # Summary stats & figures---------
  effect_stats = get_sumstats(mod, c('beta', 'k0', 'var.k'),
                              id = c(inputdat$Bpredictors %>% names(), 
                                     'k0', 'var.k')),
  # strong positive effect of flow14_3, negative effect of k0
  
  # simulated vs. observed pgr & burrows:
  fig_simulations = ggsave(filename = 'fig/simulations.png',
                           plot = plot_simulations(mod, obsdat = modeldat),
                           width = 6, height = 6, units = 'in'),

  # strength of effects
  fig_effectsizes = ggsave(filename = 'fig/effect_size.png',
                           plot = plot_sumstats(mod, c('beta', 'k0', 'var.k'),
                                                id = c(inputdat$Bpredictors %>% 
                                                         names(), 'k0', 'var.k')),
                           width = 6, height = 3, units = 'in'),
  


  fig_covpred = ggsave(filename = 'fig/covariate_effects.png',
                       plot = plot_cov_relationships(mod, inputdat = inputdat,
                                                     obsdat = modeldat),
                       width = 6, height = 5, units = 'in')
  
  # # examine variability in k estimates over the years:
  # MCMCvis::MCMCplot(mod, params = 'k')
  # inputdat2$Bpredictors %>%
  #   mutate(k = MCMCvis::MCMCpstr(mod_growth_vark, 'k', func = median)$k) %>%
  #   pivot_longer(-k) %>%
  #   ggplot(aes(value, k)) + geom_point() + geom_smooth(method = 'gam') +
  #   facet_wrap(~name, scales = 'free')
  # # possibly stronger density dependence with highest flow14 values

)
