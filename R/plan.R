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
    filter(WY >= 1999) %>% 
    select(-nstations) %>% #ranges 3-4 stations contributing data each water year
    # calculate running totals for current and 2 prior water years
    mutate(flowtot_3 = zoo::rollsum(flowtot, 3, fill = NA, align = 'right'),
           flow50_3 = zoo::rollsum(flow50, 3, fill = NA, align = 'right'),
           flow14_3 = zoo::rollsum(flow14, 3, fill = NA, align = 'right')),
  
  # # check correlations:
  # corrplot::corrplot.mixed(
  #   cor(mflow_sum %>% drop_na(), method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # # - strong correlations among annual and cumulative metrics, but not across
  # # them; no trend/correlation with water year

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
                              TRUE ~ NA_character_)) %>% 
    filter(WY >= 1999),
  
  # generate breeding season summaries:
  drought_sum = droughtdat %>% 
    filter(season == 'breeding') %>% 
    group_by(index, WY) %>%
    summarize(value = mean(value),
              .groups = 'drop') %>% 
    pivot_wider(names_from = 'index') %>% 
    drop_na(),

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
  
  # # plot relationships with pgr:
  # modeldat %>%
  #   mutate_at(vars(starts_with('flow')), ~log(. + 1)) %>%
  #   pivot_longer(-starts_with(c('year', 'pgr'))) %>%
  #   ggplot(aes(value, pgr, color = year)) + geom_point() +
  #   scale_color_viridis_c() +
  #   geom_smooth(method = 'lm') +
  #   facet_wrap(~name, scales = 'free') +
  #   xlab(NULL)
  #
  # positive relationships with pgr:
  # - burrows_tp1 (of course)
  # - flow14
  # - flow14_3 (slight)
  # - flow14_3_tp1 (slight)
  # - pdsi
  #
  # negative relationship with pgr:
  # - burrows (density dependence)
  # - flow14_tp1 - but appears driven by 3 data points with unusually high
  #    growth
  # - pdsi_tp1 (slight)
  #

  
  # # plot relationships with #burrows: (carrying capacity?)
  # modeldat %>%
  #   mutate_at(vars(starts_with('flow')), ~log(. + 1)) %>%
  #   pivot_longer(-starts_with(c('year', 'burrows'))) %>%
  #   ggplot(aes(value, burrows, color = year)) + geom_point() +
  #   scale_color_viridis_c() +
  #   geom_smooth(method = 'lm') +
  #   facet_wrap(~name, scales = 'free') +
  #   xlab(NULL)
  # positive relationship with #burrows: 
  # - flow14_3
  # - flow14_3_tp1 (very slight)
  # - flow14_tp1 (very slight)
  #
  # negative relationship:
  # - pgr
  # - flow14 (slight)
  # - pdsi (slight)
  # - pdsi_tp1 (slight)

  
  # # check correlations with birddat & across metrics (log-transform flow metrics)
  # corrplot::corrplot.mixed(
  #   cor(modeldat %>% mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
  #         drop_na(),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # strong positive correlations:
  # - flow14_3 & burrows_tp1 (0.84)
  # - pgr & flow14 (0.7)
  #
  # strongest negative correlation:
  # - burrows & year (-0.6)

  # ANALYSIS-----------
  # 1. growth rate----------
  # may be easier than modeling N because it can possibly avoid autocorrelation
  # issues 
  # --> drop 1999 only if including EVI data (missing before 2000
  # --> keep 2020 even though missing growth data (can predict growth from covariates & number of burrows in 2021)

  inputdat = list(
    # per-capita growth rate t to t+1
    y = modeldat %>% pull(pgr),
    # density-independent predictors of growth rate from t to t+1:
    Bpredictors = modeldat %>% 
      select(year, flow14_annual, flow14_cum3, pdsi) %>%
      mutate(t = year) %>% 
      # log transform flow data
      mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
      # calculate anomaly:
      pivot_longer(-year) %>%
      group_by(name) %>%
      mutate(max = max(value, na.rm = TRUE),
             min = min(value, na.rm = TRUE),
             mean = mean(value, na.rm = TRUE),
             zvalue = (value - mean) / (max - min)) %>%
      select(year, name, zvalue) %>%
      pivot_wider(names_from = name, values_from = zvalue) %>%
      select(-year),
    N = modeldat %>% pull(burrows) %>% round(digits = 0)/10000
    ),

  scale_params = modeldat %>% 
    select(year, flow14_annual, flow14_cum3, pdsi) %>%
    mutate(t = year) %>% 
    # log transform flow data
    mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
    # calculate anomaly:
    pivot_longer(-year) %>%
    group_by(name) %>%
    summarize(max = max(value, na.rm = TRUE),
              min = min(value, na.rm = TRUE),
              center = mean(value, na.rm = TRUE),
              scale = max - min,
              .groups = 'drop') %>% 
    select(-max, -min),
  
  # corrplot::corrplot.mixed(
  #   cor(bind_cols(pgr = inputdat$y,
  #                 year = inputdat$year,
  #                 inputdat$Bpredictors) %>%
  #         drop_na(),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')


  # > simple growth model----
  # # # no density dependence
  # mod_growth = fit_BANSpopgrowth_model(
  #   inputdat, n.adapt = 5000, n.burnin = 10000, n.sample = 5000, n.chains = 3,
  #   holdout = 0, type = 'growth_simple',
  #   vars = c('rate', 'beta', 'sigma', 'ysim', 'Nsim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit')),
  # MCMCvis::MCMCsummary(mod_growth,
  #                      c('rate', 'beta', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_growth,
  #                    params = c('rate', 'beta', 'sigma'),
  #                    pdf = FALSE)
  # # fits well; n.eff sufficient; trace plots well-mixed
  # plot_model_predictions(mod_growth,
  #                        year = c(2000:2019, 2000:2020),
  #                        obsdat = modeldat, center = 0, scale = 10000)
  # # predictions ok, especially toward end of time series
  # get_sumstats(mod_growth, 'beta',
  #              id = c(inputdat$Bpredictors %>% names(), 'year'))
  # # strong effects of flow14_annual (+) and year (+)
  
  # > simple logistic model-----
  # including density dependent growth relative to a constant carrying capacity
  # (K) --> for this one, because N is now a predictor, can't have missing
  # values
  # inputdat2 = list(y = inputdat$y[c(1:6, 8:20)],
  #                  year = inputdat$year[c(1:6, 8:20)],
  #                  Bpredictors = inputdat$Bpredictors[c(1:6, 8:20), ],
  #                  # pop size estimate from t (center on 10000, in thousands)
  #                  N = inputdat$N[c(1:6, 8:20)]),
  # 
  # mod_growth_logistic = fit_BANSpopgrowth_model(
  #   inputdat2, n.adapt = 5000, n.burnin = 15000, n.sample = 50000, n.chains = 3,
  #   holdout = 0, type = 'growth_logistic',
  #   vars = c('rate', 'beta', 'K', 'sigma', 'ysim', 'Nsim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 10),
  # MCMCvis::MCMCsummary(mod_growth_logistic,
  #                      c('rate', 'beta', 'K', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_growth_logistic,
  #                    params = c('rate', 'beta', 'K', 'sigma'),
  #                    pdf = FALSE)
  # # fits well; n.eff low for K & trace plots a weird for K
  # plot_model_predictions(mod_growth_logistic,
  #                        year = c(2000:2005, 2007:2019, 2000:2006, 2008:2020),
  #                        obsdat = modeldat,
  #                        scale = 10000, center = 0)
  # # predictions reasonable
  # get_sumstats(mod_growth_logistic, 'beta',
  #              id = c(inputdat2$Bpredictors %>% names(), 'year'))
  # # strong effects of flow14_cum3 (+) only

  # > simple logistic 2-------
  # alt parameterization to directly estimate strength of density dependence
  # without having to estimate K; effect of density dependence still assumed to
  # be constant
  # mod_growth_logistic2 = fit_BANSpopgrowth_model(
  #   inputdat2, n.adapt = 5000, n.burnin = 15000, n.sample = 150000, n.chains = 3,
  #   holdout = 0, type = 'growth_logistic_alt',
  #   vars = c('rate', 'beta', 'dd', 'sigma', 'ysim', 'Nsim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 50),
  # MCMCvis::MCMCsummary(mod_growth_logistic2,
  #                      c('rate', 'beta', 'dd', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_growth_logistic2,
  #                    params = c('rate', 'beta', 'dd', 'sigma'),
  #                    pdf = FALSE)
  # # fits well; n.eff sufficient; trace plots well-mixed
  # plot_model_predictions(mod_growth_logistic2,
  #                       year = c(2000:2005, 2007:2019, 2000:2006, 2008:2020),
  #                       obsdat = modeldat,
  #                       center = 0, scale = 10000)
  # # predictions reasonable
  # get_sumstats(mod_growth_logistic2, c('beta', 'dd'),
  #              id = c(inputdat2$Bpredictors %>% names(), 'year', 'dd'))
  # # strong effects of flow14_cum3 (+) and dd (-); pdsi  and flow14_annual close-ish?
  
  # > logistic with variable density dependence------
  # allows strength of density dependence to vary annually from a normal
  # distribution with mean mu.dd and error sigma.d
  # mod_growth_vardd = fit_BANSpopgrowth_model(
  #   inputdat2, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
  #   holdout = 0, type = 'growth_logistic_variabledd',
  #   vars = c('rate', 'beta', 'sigma', 'dd', 'mu.d', 'sigma.d',
  #            'ysim', 'Nsim', 'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 100),
  # MCMCvis::MCMCsummary(mod_growth_vardd,
  #                      c('rate', 'beta', 'sigma', 'mu.d', 'sigma.d',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_growth_vardd,
  #                    params = c('rate', 'beta', 'mu.d', 'sigma', 'sigma.d'),
  #                    pdf = FALSE)
  # # fits well; n.eff low for mu.d, rate; trace plots a bit odd as well for these
  # plot_model_predictions(mod_growth_vardd,
  #                        year = c(2000:2005, 2007:2019, 2000:2006, 2008:2020),
  #                        obsdat = modeldat, scale = 10000, center = 0)
  # # predictions reasonable
  # get_sumstats(mod_growth_vardd, c('beta','mu.d'),
  #              id = c(inputdat2$Bpredictors %>% names(), 'year', 'dd'))
  # # strong effects of flow14_cum3 (+) and dd (-)
  
  # > logistic with trend in K-------
  # mod_growth_trendK = fit_BANSpopgrowth_model(
  #   inputdat2, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
  #   holdout = 0, type = 'growth_logistic_trenddd',
  #   vars = c('rate', 'beta', 'dd', 'kappa', 'dd_intercept', 'sigma',
  #            'ysim', 'Nsim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 100),
  # MCMCvis::MCMCsummary(mod_growth_trendK,
  #                      c('rate', 'beta', 'kappa', 'dd_intercept', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_growth_trendK,
  #                    params = c('rate', 'beta', 'kappa', 'dd_intercept', 'sigma'),
  #                    pdf = FALSE)
  # # fits well; n.eff low dd_intercept, kappa, beta[4], rate (confounded?)
  # plot_model_predictions(mod_growth_trendK,
  #                       year = c(2000:2005, 2007:2019, 2000:2006, 2008:2020),
  #                       obsdat = modeldat, scale = 10000, center = 0)
  # # predictions reasonable
  # get_sumstats(mod_growth_trendK, c('beta','kappa','dd_intercept'),
  #              id = c(inputdat2$Bpredictors %>% names(), 'year', 'dd_slope', 'dd'))
  # # strong effects of flow14_cum3 (+) and dd (-) [not dd_slope]
  
  # > logistic with variable K-------
  # let carrying capacity change with cumulative flow; now that N is a
  # predictor, have to drop missing values
  inputdat3 = list(y = inputdat$y[c(1:7,9:22)],
                   Bpredictors = inputdat$Bpredictors[c(1:7,9:22), c('flow14_annual', 'pdsi', 't')],
                   Kpredictors = inputdat$Bpredictors[c(1:7,9:22), 'flow14_cum3'],
                   # pop size estimate from t (center on 15000, scale in tens of
                   # thousands)
                   N = inputdat$N[c(1:7,9:22)],
                   predmatrix = tibble(flow14_annual = seq(-0.6, 0.4, length.out = 22),
                                       pdsi = seq(-0.4, 0.6, length.out = 22),
                                       t = inputdat$Bpredictors$t,
                                       flow14_cum3 = seq(-0.6, 0.4, length.out = 22))),

  mod_growth_varK = fit_BANSpopgrowth_model(
    inputdat3, n.adapt = 10000, n.burnin = 40000, n.sample = 600000, 
    n.chains = 3, holdout = 0, type = 'growth_logistic_varK',
    vars = c('r0', 'beta', 'kappa', 'k0', 'k', 
             'var.o', 'var.p', 'var.k', 'ysim', 'Nsim', 'ypred', 
             'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
    thin = 300),
  MCMCvis::MCMCsummary(mod_growth_varK,
                       c('r0', 'beta', 'kappa', 'k0', 
                         'var.o', 'var.p', 'var.k', 
                         'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  MCMCvis::MCMCtrace(mod_growth_varK,
                     params = c('r0', 'beta', 'kappa', 'k0',
                                'var.o', 'var.p', 'var.k'),
                     pdf = FALSE)
  # model fits well; n.eff sufficient; trace plots well mixed
  plot_model_predictions(mod_growth_varK,
                         year = c(1999:2005, 2007:2020, 1999:2021),
                         obsdat = modeldat, scale = 10000, center = 0)
  plot_cov_relationships(mod_growth_varK,
                         r0 = 0,
                         nm = names(inputdat3$predmatrix),
                         xvals = inputdat3$predmatrix %>% 
                           mutate(index = c(1:nrow(.))) %>% 
                           pivot_longer(-index) %>% 
                           left_join(scale_params) %>% 
                           mutate(value = value * scale + center,
                                  value = case_when(name %in% c('flow14_annual', 'flow14_cum3') ~ exp(value),
                                                    TRUE ~ value)) %>% 
                           select(index, name, value) %>% 
                           arrange(name, index),
                         obsdat = modeldat %>%
                           select(pgr, flow14_annual, pdsi, flow14_cum3, t = year) %>% 
                           pivot_longer(-pgr, names_to = 'cov'))
  # predictions reasonable
  get_sumstats(mod_growth_varK, c('beta', 'kappa', 'k0'),
               id = c(inputdat3$Bpredictors %>% names(), 
                      inputdat3$Kpredictors %>% names(), 'k0'))
  # strong effects of flow14_cum3 on strength of density dependence (+); dd
  # negative; other effects as expected - trending positive for flow14_annual &
  # pdsi, no trend in year
  plot_K_predictions(mod_growth_varK,
                     predictor = seq(400000, 10000000, 500000),
                     obsdat = modeldat,
                     scale = 10000, center = 0) + ylim(0, NA)
  
  
  
  
  
  # >>hold out version----------
  mod_growth_varK_holdout = fit_BANSpopgrowth_model(
    inputdat3, n.adapt = 10000, n.burnin = 40000, n.sample = 600000, 
    n.chains = 3, holdout = 3, type = 'growth_logistic_varK',
    vars = c('rate', 'beta', 'kappa', 'dd_intercept', 'dd', 
             'var.o', 'var.p', 'var.d', 'ysim', 'Nsim', 
             'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
    thin = 200),
  MCMCvis::MCMCsummary(mod_growth_varK_holdout,
                       c('rate', 'beta', 'kappa', 'dd_intercept', 'var.o', 'var.p',
                         'var.d', 'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  MCMCvis::MCMCtrace(mod_growth_varK_holdout,
                     params = c('rate', 'beta', 'kappa', 'dd_intercept',
                                'var.o', 'var.p', 'var.d'),
                     pdf = FALSE)
  # model fits well; n.eff sufficient; trace plots well mixed
  plot_model_predictions(mod_growth_varK_holdout,
                         year = c(1999:2005, 2007:2020, 1999:2006, 2008:2021),
                         obsdat = modeldat, scale = 10000, center = 0)
  # predictions reasonable
  get_sumstats(mod_growth_varK_holdout, c('beta', 'kappa', 'dd_intercept'),
               id = c(inputdat3$Bpredictors %>% names(), 
                      inputdat3$Kpredictors %>% names(), 'dd'))
  # strong effects of flow14_cum3 on strength of density dependence (+); dd
  # negative; other effects as expected - trending positive for flow14_annual &
  # pdsi, no trend in year
  plot_K_predictions(mod_growth_varK_holdout,
                     predictor = seq(400000, 10000000, 500000),
                     obsdat = modeldat,
                     scale = 10000, center = 0) + ylim(0, NA)
  
  # 2. N----------
  # these models are slower to run, and are parameterized differently to
  # directly take burrow counts in as the "y" observations, estimated with error
  # as a count of the true population N, and essentially calculate the growth
  # rate internally.
  
  # re-do input data: burrow counts need to be integers, so no scaling, but can
  # have missing values; include counts through 2020 to account for growth rate
  # estimated for 2019-20 interval
  inputdatN = list(y = modeldat %>% filter(year > 1999) %>%
                     pull(burrows) %>% round(digits = 0),
                   year = (modeldat %>% filter(year > 1999) %>%
                             pull(year) - 2010)/10,
                   Bpredictors = inputdat$Bpredictors),

  # > simple poisson-------
  # mod_N_poisson = fit_BANSpopgrowth_model(
  #   inputdatN, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
  #   holdout = 0, type = 'N_simple_poisson',
  #   vars = c('rate', 'beta', 'error', 'Nsim', 'ysim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 100),
  # MCMCvis::MCMCsummary(mod_N_poisson,
  #                      c('rate', 'beta',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_N_poisson,
  #                    params = c('rate', 'beta'),
  #                    pdf = FALSE)
  # # fit ok; n.eff very low; trace plots terrible; predictions suspiciously precise
  # plot_model_predictions(mod_N_poisson,
  #                        year = c(2000:2019, 2000:2020),
  #                        obsdat = modeldat, scale = 1, center = 0)
  # # predictions reasonable
  # get_sumstats(mod_N_poisson, c('beta'),
  #              id = c(inputdatN$Bpredictors %>% names(), 'year'))
  # # strong effects of flow14_annual (+) and year (+)
  
  # # > logistic poisson---------
  # mod_N_logistic_poisson = fit_BANSpopgrowth_model(
  #   inputdatN, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
  #   holdout = 0, type = 'N_logistic_poisson',
  #   vars = c('rate', 'beta', 'dd', 'sigma', 'Nsim', 'ysim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 100)
  # MCMCvis::MCMCsummary(mod_N_logistic_poisson,
  #                      c('rate', 'beta', 'dd', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_N_logistic_poisson,
  #                    params = c('rate', 'beta', 'dd', 'sigma'),
  #                    pdf = FALSE)
  # # fit ok; n.eff very low; dd tiny; trace plots terrible 
  # plot_model_predictions(mod_N_logistic_poisson,
  #                        year = c(2000:2019, 2000:2020),
  #                        obsdat = modeldat, scale = 1, center = 0)
  # # predictions suspiciously precise
  # get_sumstats(mod_N_logistic_poisson, c('beta', 'dd'),
  #              id = c(inputdatN$Bpredictors %>% names(), 'year', 'dd'))
  # # strong effects of flow14_annual (+) only

  # > simple negative binomial----------
  # allows more error in the count process than poisson

  # mod_N_negbin = fit_BANSpopgrowth_model(
  #   inputdatN, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
  #   holdout = 0, type = 'N_simple_negbin',
  #   vars = c('rate', 'beta', 'sigma', 'Nsim', 'ysim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 100),
  # MCMCvis::MCMCsummary(mod_N_negbin,
  #                      c('rate', 'beta', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_N_negbin,
  #                    params = c('rate', 'beta', 'sigma'),
  #                    pdf = FALSE)
  # # fit ok; n.eff a bit low for rate, with ok traceplots
  # plot_model_predictions(mod_N_negbin,
  #                        year = c(2000:2019, 2000:2020),
  #                        obsdat = modeldat, scale = 1, center = 0)
  # # predictions reasonable
  # get_sumstats(mod_N_negbin, c('beta'),
  #              id = c(inputdatN$Bpredictors %>% names(), 'year'))
  # # strong effects of flow14_annual (+) and year (+)

  # > logistic negative binomial----------
  # allows more error in the count process than poisson; add density dependence

  mod_N_logistic_negbin = fit_BANSpopgrowth_model(
    inputdatN, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
    holdout = 0, type = 'N_logistic_negbin',
    vars = c('rate', 'beta', 'dd', 'sigma', 'k', 'Nsim', 'ysim',
             'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
    thin = 100),
  # MCMCvis::MCMCsummary(mod_N_logistic_negbin,
  #                      c('rate', 'beta', 'dd', 'sigma', 'k',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_N_logistic_negbin,
  #                    params = c('rate', 'beta', 'dd', 'sigma', 'k'),
  #                    pdf = FALSE)
  # # fit very good; n.eff super low for rate, dd, beta[2] and beta[4]; trace plots for
  # # these also terrible -- looks like a longer burn-in needed?
  # plot_model_predictions(mod_N_logistic_negbin,
  #                        year = c(2000:2019, 2000:2020),
  #                        obsdat = modeldat, scale = 1, center = 0)
  # # predictions very good...
  # get_sumstats(mod_N_logistic_negbin, c('beta', 'dd'),
  #              id = c(inputdatN$Bpredictors %>% names(), 'year', 'dd'))
  # # strong effects of flow14_cum3 (+) and dd (-)

  # > logistic negbin variable K----------
  inputdatN2 = list(y = inputdatN$y,
                    year = inputdatN$year,
                    Bpredictors = inputdatN$Bpredictors[, c('flow14_annual', 'pdsi')],
                    Kpredictors = inputdatN$Bpredictors[, 'flow14_cum3']),

  mod_N_logistic_negbin_varK = fit_BANSpopgrowth_model(
    inputdatN2, n.adapt = 5000, n.burnin = 20000, n.sample = 100000, n.chains = 3,
    holdout = 0, type = 'N_logistic_negbin_varK',
    vars = c('rate', 'beta', 'kappa', 'dd_intercept', 'sigma', 'sigma.d', 'dd',
             'k', 'Nsim', 'ysim',
             'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
    thin = 50),
  MCMCvis::MCMCsummary(mod_N_logistic_negbin_varK,
                       c('rate', 'beta', 'kappa', 'dd_intercept',
                         'sigma', 'sigma.d', 'k',
                         'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  MCMCvis::MCMCtrace(mod_N_logistic_negbin_varK,
                     params = c('rate', 'beta', 'kappa', 'dd_intercept',
                                'sigma', 'sigma.d', 'k'),
                     pdf = FALSE)
  # fit ok; n.eff=0 for sigma.d? also low for rate, dd_intercept, kappa;
  # trace plots a bit weird too
  # --> note, the estimates for sigma.d and sigma are super tiny- need a different prior?
  # also really small for kappa
  plot_model_predictions(mod_N_logistic_negbin_varK,
                         year = c(2000:2019, 2000:2020),
                         obsdat = modeldat, scale = 1, center = 0)
  # predictions reasonable
  get_sumstats(mod_N_logistic_negbin_varK, c('beta', 'kappa', 'dd_intercept'),
               id = c(inputdatN2$Bpredictors %>% names(), 'year',
                      inputdatN2$Kpredictors %>% names(), 'dd'))
  # strong effects of flow14_cum3 on density dependence (+) and of dd (-);
  # expected suggestions of effects for flow14_annual (+ 82%), pdsi (+ 77%), and
  # year (- 76%)
  
  plot_K_predictions(mod_N_logistic_negbin_varK,
                     predictor = seq(400000, 10000000, 500000),
                     obsdat = modeldat,
                     scale = 1, center = 0)
)
