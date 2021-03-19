# Custom functions are an important part of a drake workflow.
# This is where you write them.
# Details: https://books.ropensci.org/drake/plans.html#functions

render_Rmd = function(pathin, pathout) {
  rmarkdown::render(pathin, 
                    output_file = here::here(pathout))
}

# FUNCTION FAILS: alert handshake failure
# get_cdec_waterdat = function(station = 'VIN',
#                              sensor = 41, #mean daily flow
#                             start.date = '2010-10-01',
#                             end.date = as.character(Sys.Date())) {
#   # if (length(stations) > 1) {
#   #   stations = paste(st, collapse = '%2C+')
#   # }
#   # URL = 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet'
#   URL = 'https://cdec.water.ca.gov/dynamicapp/req/JasonDataServlet'
#   query = paste0('?', station, '&SensorNums=', sensor, '&dur_code=D',
#                  '&Start=', start.date, '&End=', end.date)
#   RCurl::getURL(paste0(URL, query))
# }

compile_waterdat = function(dir, type, ...) {
  filelist = list.files(dir, pattern = type, full.names = TRUE)
  purrr::map_df(filelist, read_csv, ...)
}

get_drought_indices = function(datname, statecode = NULL, climdivcode = NULL) {
  # 04 = california
  
  url <- "ftp://ftp.ncdc.noaa.gov/pub/data/cirs/climdiv/"
  filelist = RCurl::getURL(url, dirlistonly = TRUE) %>% 
    strsplit("\\r\\n")
  df = tibble(name = filelist[[1]]) %>% 
    filter(grepl(datname, name))
  dat <- try(RCurl::getURL(paste0(url, df$name[1]))) %>% 
    strsplit("   \\n")
  res <- read_table(dat[[1]],
             col_names = c('ID', 'JAN', 'FEB', 'MAR', 'APR', 
                           'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 
                           'OCT', 'NOV', 'DEC')) %>% 
    mutate(index = datname,
           state = substr(ID, 1, 2),
           climdiv = paste0(substr(ID, 3, 4)),
           year = substr(ID, 7, 10)) %>%
    pivot_longer(JAN:DEC, names_to = 'month', values_to = 'value') %>%
    mutate(climname = recode(climdiv, 
                            `01` = "North Coast Drainage",
                            `02` = "Sacramento Drainage",
                            `03` = "Northeast Interior Basins",
                            `04` = "Central Coast",
                            `05` = "San Joaquin Drainage",
                            `06` = "South Coast Drainage",
                            `07` = "Southeast Desert Basin"),
           date = as.Date(paste0(year, '-', month, '-01'), format = '%Y-%b-%d'),
           value = case_when(value < -99.9 ~ NA_real_, 
                             TRUE ~ value)) %>%
    arrange(index, state, climdiv, date)
  
  if (!is.null(statecode)) {
    res <- res %>% filter(state == statecode)
  }
  if (!is.null(climdivcode)) {
    res <- res %>% filter(climdiv == climdivcode)
  }
  return(res)
  ## Climate divisions within California:
  ## 01: North Coast Drainage
  ## 02: Sacramento Drainage
  ## 03: Northeast Interior Basins
  ## 04: Central Coast **TomKat**
  ## 05: San Joaquin Drainage
  ## 06: South Coast Drainage
  ## 07: Southeast Desert Basin
}

get_sumstats = function(mod, param, id, prob = 0.95) {
  samples = MCMCvis::MCMCchains(mod, params = param, ISB = TRUE)
  
  bind_cols(
    id,
    MCMCvis::MCMCpstr(mod, params = param, func = median) %>% unlist(),
    do.call(rbind, 
            MCMCvis::MCMCpstr(mod, params = param, 
                              func = function(y) HDInterval::hdi(y, prob))) %>% 
      as_tibble(),
    apply(samples, 2, function(x) 1-ecdf(x)(0)) %>% set_names(NULL),
    apply(samples, 2, function(x) ecdf(x)(0)) %>% set_names(NULL)
  ) %>% 
    set_names('id', 'median', 'lci', 'uci', 'probg0', 'probl0')
}

plot_sumstats = function(mod, param, id, prob = 0.95) {
  get_sumstats(mod, param, id, prob) %>%
    filter(id %in% c('flow14', 'flow14_3', 'pdsi', 't', 'k0')) %>%
    mutate(id = factor(id, levels = c('flow14', 'flow14_3', 'pdsi', 't', 'k0')),
           id = recode(id, 
                       'pdsi' = 'mean.pdsi',
                       'flow14' = 'annual.flow',
                       'flow14_3' = 'cumulative.flow',
                       'k0' = 'mu.kt'),
           label = if_else(round(probl0, 2) >= 0.9 | round(probg0, 2) >= 0.9, 
                           '*', '')) %>%
    ggplot(aes(id, median, ymin = lci, ymax = uci)) +
    geom_hline(aes(yintercept = 0), linetype = 'dashed') +
    geom_point() + geom_errorbar(width = 0.25) +
    geom_text(aes(y = uci + 0.01, label = label), size = 4) +
    # geom_text(aes(y = uci + 0.03, label = paste0(round(probg0 * 100, 0), '%')),
    #           size = 2.5) +
    labs(x = NULL, y = 'effect size') + 
    theme_classic() +
    coord_flip()
}

plot_covariates = function(df) {
  a = ggplot(df, aes(year, flow14/1000000)) + geom_line() + geom_point() +
    labs(title = 'A', x = NULL, y = 'cfs (millions)') +
    scale_x_continuous(breaks = seq(1999, 2020, 1),
                       labels = c(1999, '', 2001, '', 2003, '', 2005, '', 
                                  2007, '', 2009, '', 2011, '', 2013, '',
                                  2015, '', 2017, '', 2019, '')) +
    theme_classic()
  b = ggplot(df, aes(year, flow14_3/1000000)) + geom_line() + geom_point() +
    labs(title = 'B', x = NULL, y = 'cfs (millions)') +
    scale_x_continuous(breaks = seq(1999, 2020, 1),
                       labels = c(1999, '', 2001, '', 2003, '', 2005, '', 
                                  2007, '', 2009, '', 2011, '', 2013, '',
                                  2015, '', 2017, '', 2019, '')) + 
    ylim(0, 11) + 
    theme_classic()
  c = ggplot(df, aes(year, pdsi)) + geom_line() + geom_point() +
    labs(title = 'C', x = NULL, y = 'PDSI') +
    scale_x_continuous(breaks = seq(1999, 2020, 1),
                       labels = c(1999, '', 2001, '', 2003, '', 2005, '', 
                                  2007, '', 2009, '', 2011, '', 2013, '',
                                  2015, '', 2017, '', 2019, '')) + 
    theme_classic()
  
  library(patchwork)
  a/b/c
}

plot_simulations = function(mod, obsdat) {
  
  # growth rate simulations:
  pgr = bind_cols(MCMCvis::MCMCpstr(mod, params = c('ysim', 'ysim.2006'), 
                                    func = function(x) median(x)) %>%
                    purrr::map_df(as_tibble),
                  MCMCvis::MCMCpstr(mod, params = c('ysim', 'ysim.2006'),
                                    func = function(x) HDInterval::hdi(x, 0.95)) %>%
                    purrr::map_df(as_tibble)) %>% 
    mutate(year = c(1999:2005, 2007:2020, 2006)) %>% 
    arrange(year) %>% 
    left_join(obsdat %>% select(year, pgr))
  
  cat('\nForecasted growth rate for', max(pgr$year), ': ', 
      pgr %>% filter(year == max(year)) %>% pull(value),
      '(95% HPDI: ',
      pgr %>% filter(year == max(year)) %>% pull(lower), '-', 
      pgr %>% filter(year == max(year)) %>% pull(upper), ')'
      )

  N = bind_cols(MCMCvis::MCMCpstr(mod, params = 'Nsim', 
                                    func = function(x) median(x)) %>%
                    purrr::map_df(as_tibble),
                  MCMCvis::MCMCpstr(mod, params = 'Nsim',
                                    func = function(x) HDInterval::hdi(x, 0.95)) %>%
                    purrr::map_df(as_tibble)) %>% 
    mutate(year = c(1999:2021)) %>% 
    arrange(year) %>% 
    left_join(obsdat %>% select(year, burrows))
  
  cat('\nForecasted burrow count for', max(N$year), ': ', 
      N %>% filter(year == max(year)) %>% pull(value),
      '(95% HPDI: ',
      N %>% filter(year == max(year)) %>% pull(lower), '-', 
      N %>% filter(year == max(year)) %>% pull(upper), ')'
  )
  
  a = ggplot(pgr, aes(year, value, ymin = lower, ymax = upper)) + 
    geom_ribbon(fill = 'gray80') + geom_line() + geom_point() +
    labs(title = 'A', x = NULL, y = 'per-capita growth rate') + 
    geom_line(aes(y = pgr), color = 'red') +
    geom_point(aes(y = pgr), color = 'red') +
    scale_x_continuous(breaks = seq(1999, 2020, 1),
                       labels = c(1999, '', 2001, '', 2003, '', 2005, '', 
                                  2007, '', 2009, '', 2011, '', 2013, '',
                                  2015, '', 2017, '', 2019, '')) +
    ylim(-1, 1) +
    theme_classic()
  
  b = ggplot(N, aes(year, value/1000, ymin = lower/1000, ymax = upper/1000)) + 
    geom_ribbon(fill = 'gray80') + geom_line() + geom_point() +
    labs(title = 'B', x = NULL, y = 'burrows (thousands)') + 
    geom_line(aes(y = burrows/1000), color = 'red') +
    geom_point(aes(y = burrows/1000), color = 'red') +
    scale_x_continuous(breaks = seq(1999, 2020, 1),
                       labels = c(1999, '', 2001, '', 2003, '', 2005, '', 
                                  2007, '', 2009, '', 2011, '', 2013, '',
                                  2015, '', 2017, '', 2019, '')) +
    ylim(0, 25) +
    theme_classic()
  
  library(patchwork)
  a/b
}

plot_cov_relationships = function(mod, inputdat, obsdat) {
  
  nm = names(inputdat$predmatrix)
  
  scale_params = purrr::map_df(nm[1:4] %>% set_names(), 
                               function(x) {
                                 tibble(center = inputdat$Bpredictors %>% pull(x) %>% 
                                          attr(., 'scaled:center'),
                                        scale = inputdat$Bpredictors %>% pull(x) %>% 
                                          attr(., 'scaled:scale'))
                               }, .id = 'cov') %>% 
    bind_rows(tibble(cov = 'N', center = mean(inputdat$N_obs, na.rm = TRUE),
                     scale = sd(inputdat$N_obs, na.rm = TRUE)))
    
  # pull other stats:
  predvalues = bind_cols(MCMCvis::MCMCpstr(mod, params = c('ypred'), 
                                           func = function(x) mean(x)) %>%
                           purrr::map_df(as_tibble),
                         MCMCvis::MCMCpstr(mod, params = c('ypred'),
                                           func = function(x) HDInterval::hdi(x, 0.95)) %>%
                           purrr::map_df(as_tibble)) %>%
    set_names(paste(rep(nm, 3), rep(c('median', 'lower', 'upper'), each = 5),
                    sep = '.')) %>% 
    mutate(index = c(1:nrow(.))) %>% 
    pivot_longer(-index) %>% 
    separate(name, into = c('cov', 'metric'), sep = '\\.') %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    arrange(cov, index) %>% 
    left_join(inputdat$predmatrix %>% mutate(index = c(1:nrow(.))) %>%
                pivot_longer(-index) %>%
                left_join(scale_params, by = c('name' = 'cov')) %>%
                mutate(value = value * scale + center,
                       value = case_when(name %in% c('flow14', 'flow14_3') ~ exp(value)/1000000,
                                         TRUE ~ value)
                ) %>%
                select(index, name, value) %>%
                arrange(name, index), 
              by = c('index', 'cov' = 'name')) %>% 
    mutate(cov = recode(cov, 
                        'pdsi' = 'mean.pdsi',
                        'flow14' = 'annual.flow',
                        'flow14_3' = 'cumulative.flow',
                        't' = 't'))
  
  a = ggplot(predvalues %>% filter(cov == 'annual.flow'), aes(value, median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
    geom_hline(aes(yintercept = 0), linetype = 'dashed') +
    geom_line() + geom_point() +
    geom_point(data = obsdat, aes(flow14/1000000, pgr), color = 'red') +
    scale_x_log10() + ylim(-1, 1) +
    labs(title = 'A', x = 'cfs (millions)', y = 'per-capita growth rate') +
    theme_classic()
  
  b = ggplot(predvalues %>% filter(cov == 'cumulative.flow'), aes(value, median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
    geom_hline(aes(yintercept = 0), linetype = 'dashed') +
    geom_line() + geom_point() +
    geom_point(data = obsdat, aes(flow14_3/1000000, pgr), color = 'red') +
    scale_x_log10() + ylim(-1, 1) +
    labs(title = 'B', x = 'cfs (millions)', y = 'per-capita growth rate') +
    theme_classic()
  
  c = ggplot(predvalues %>% filter(cov == 'mean.pdsi'), aes(value, median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
    geom_hline(aes(yintercept = 0), linetype = 'dashed') +
    geom_line() + geom_point() +
    geom_point(data = obsdat, aes(pdsi, pgr), color = 'red') +
    labs(title = 'C', x = 'PDSI', y = NULL) +
    ylim(-1, 1) +
    theme_classic()
  
  d = ggplot(predvalues %>% filter(cov == 't'), aes(value, median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
    geom_hline(aes(yintercept = 0), linetype = 'dashed') +
    geom_line() + geom_point() +
    geom_point(data = obsdat, aes(year, pgr), color = 'red') +
    labs(title = 'D', x = 'year', y = NULL) +
    ylim(-1, 1) +
    theme_classic()
  
  library(patchwork)
  (a/b)|(c/d)

}


plot_K_predictions = function(mod, predictor, obsdat, scale = 1, center = 0) {
  # # PREDICTED CARRYING CAPACITY (in same scale/units as N)
  # # for a range of the predictor values
  # for (j in 1:npred) {
  #   dsim[j] = dd_intercept + beta[3] * Kpred[j]
  #   K[j] = -1 * rate/dsim[j]
  # }
  dat = MCMCvis::MCMCchains(mod, params = c('dd_intercept', 'kappa', 'rate'))
  xvals = (log(predictor + 1) - 14.64198)/(15.96866-12.91454)
  
  predvalues = purrr::map(xvals,
                          function(x) -dat[,'rate'] / 
                            (dat[,'dd_intercept'] + x*dat[,'kappa']) * scale) %>% 
    purrr::map_df(function(x) {
      c('value' = median(x), HDInterval::hdi(x, 0.95))
    }) %>% 
    mutate(predictor = predictor)
  
  p = ggplot(predvalues, aes(predictor/1000, value)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') +
    geom_point() + geom_line() +
    geom_point(data = obsdat, aes(flow14_cum3/1000, burrows, 
                                  color = as.factor(year))) +
    scale_color_viridis_d() +
    xlab('Cumulative flow >14cfs (3 years) (thousands)') +
    ylab('carrying capacity')
  
  return(p)
}
