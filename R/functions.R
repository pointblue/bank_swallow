# Custom functions are an important part of a drake workflow.
# This is where you write them.
# Details: https://books.ropensci.org/drake/plans.html#functions

render_Rmd = function(pathin, pathout) {
  rmarkdown::render(pathin, 
                    output_file = here::here(pathout))
}

get_drought_indices = function(datname) {
  # 04 = california
  
  url <- "ftp://ftp.ncdc.noaa.gov/pub/data/cirs/climdiv/"
  filelist = RCurl::getURL(url, dirlistonly = TRUE) %>% 
    strsplit("\\r\\n")
  df = tibble(name = filelist[[1]]) %>% 
    filter(grepl(datname, name))
  dat <- try(RCurl::getURL(paste0(url, df$name[1]))) %>% 
    strsplit("   \\n")
  read_table(dat[[1]],
             col_names = c('ID', 'JAN', 'FEB', 'MAR', 'APR', 
                           'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 
                           'OCT', 'NOV', 'DEC')) %>% 
    mutate(index = datname)
}

cleanup_droughtdat = function(df, statecode = NULL, climdivcode = NULL) {
  res <- df %>% 
    mutate(state = substr(ID, 1, 2),
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
           month = format(date, '%m') %>% as.numeric(),
           year = as.numeric(year),
           WY = if_else(month >= 10, year + 1, year),
           season = case_when(month >= 10 | month <= 3 ~ 'rainy',
                              month >= 4 & month <= 8 ~ 'breeding',
                              TRUE ~ NA_character_),
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
    set_names(c('id', 'median', 'lci', 'uci', 'probg0', 'probl0'))
}

plot_sumstats = function(mod, param, id, prob = 0.95, 
                         var) {
  get_sumstats(mod, param, id, prob) %>%
    filter(id %in% id) %>%
    mutate(id = factor(id, levels = id),
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
  a = ggplot(df, aes(WY, flow_concurrent/1000000)) + geom_line() + geom_point() +
    labs(title = 'A', x = NULL, y = 'cfs (millions)') +
    scale_x_continuous(breaks = seq(1999, 2020, 1),
                       labels = c(1999, '', 2001, '', 2003, '', 2005, '', 
                                  2007, '', 2009, '', 2011, '', 2013, '',
                                  2015, '', 2017, '', 2019, '')) +
    theme_classic()
  b = ggplot(df, aes(WY, flow_prior/1000000)) + geom_line() + geom_point() +
    labs(title = 'B', x = NULL, y = 'cfs (millions)') +
    scale_x_continuous(breaks = seq(1999, 2020, 1),
                       labels = c(1999, '', 2001, '', 2003, '', 2005, '', 
                                  2007, '', 2009, '', 2011, '', 2013, '',
                                  2015, '', 2017, '', 2019, '')) + 
    #ylim(0, 8) + 
    theme_classic()
  c = ggplot(df, aes(WY, flow_prior3/1000000)) + geom_line() + geom_point() +
    labs(title = 'C', x = NULL, y = 'cfs (millions)') +
    scale_x_continuous(breaks = seq(1999, 2020, 1),
                       labels = c(1999, '', 2001, '', 2003, '', 2005, '', 
                                  2007, '', 2009, '', 2011, '', 2013, '',
                                  2015, '', 2017, '', 2019, '')) + 
    #ylim(0, 8) + 
    theme_classic()
  d = ggplot(df, aes(WY, pdsi_prior)) + geom_line() + geom_point() +
    labs(title = 'D', x = NULL, y = 'PDSI') +
    scale_x_continuous(breaks = seq(1999, 2020, 1),
                       labels = c(1999, '', 2001, '', 2003, '', 2005, '', 
                                  2007, '', 2009, '', 2011, '', 2013, '',
                                  2015, '', 2017, '', 2019, '')) + 
    theme_classic()
  
  library(patchwork)
  a/b/c
}

extract_simulations = function(mod, obsdat) {
  # growth rate simulations:
  pgr = bind_cols(MCMCvis::MCMCpstr(mod, params = c('ysim', 'ysim.2007'), 
                                    func = function(x) median(x)) %>%
                    purrr::map_df(as_tibble),
                  MCMCvis::MCMCpstr(mod, params = c('ysim', 'ysim.2007'),
                                    func = function(x) HDInterval::hdi(x, 0.95)) %>%
                    purrr::map_df(as_tibble)) %>%
    mutate(WY = c(1999:2005, 2007:2020, 2006)) %>%
    arrange(WY) %>%
    left_join(obsdat %>% select(WY, pgr))
  
  N = bind_cols(MCMCvis::MCMCpstr(mod, params = 'Nsim', 
                                  func = function(x) median(x)) %>%
                  purrr::map_df(as_tibble),
                MCMCvis::MCMCpstr(mod, params = 'Nsim',
                                  func = function(x) HDInterval::hdi(x, 0.95)) %>%
                  purrr::map_df(as_tibble)) %>% 
    mutate(WY = c(1998:2020)) %>% 
    arrange(WY) %>% 
    left_join(obsdat %>% select(WY, burrows))
  
  return(
    list(pgr = pgr, N = N)
  )
}

plot_simulations = function(sim, ylim1 = c(-1, 1), ylim2 = c(0, 25), maxyear = NULL) {
  
  if (!is.null(maxyear)) {
    pgr = sim$pgr %>% filter(WY <= maxyear)
    N = sim$N %>% filter(WY <= maxyear)
  } else {
    pgr = sim$pgr
    N = sim$N
  }

  a = ggplot(pgr, aes(WY, mean, ymin = `95%_HPDL`, ymax = `95%_HPDU`)) + 
    geom_ribbon(fill = 'gray85') + geom_line() + geom_point() +
    labs(title = 'A', x = NULL, y = 'per-capita growth rate') + 
    geom_line(aes(y = pgr), color = 'red') +
    geom_point(aes(y = pgr), color = 'red') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(limits = c(1999, max(N$WY), 2),
                       breaks = seq(1999, max(N$WY), 2),
                       minor_breaks = seq(2000, 2020, 2)) +
    scale_y_continuous(limits = ylim1) +
    theme_classic()
  
  b = ggplot(N, aes(WY, mean, ymin = `95%_HPDL`, ymax = `95%_HPDU`)) + 
    geom_ribbon(fill = 'gray85') + geom_line() + geom_point() +
    labs(title = 'B', x = NULL, y = 'burrows (thousands)') + 
    geom_line(aes(y = burrowst), color = 'red') +
    geom_point(aes(y = burrowst), color = 'red') +
    scale_x_continuous(limits = c(1999, max(N$WY), 2),
                       breaks = seq(1999, max(N$WY), 2),
                       minor_breaks = seq(2000, max(N$WY), 2)) +
    scale_y_continuous(limits = ylim2, breaks = seq(0, 25, 5)) +
    theme_classic()
  
  library(patchwork)
  a/b
}

generate_predictions = function(mod, inputdat, obsdat, logflow = TRUE, density = FALSE) {
  
  if (density) {
    scale_params = tibble(cov = 'burrows_prior', 
                          center = obsdat %>% filter(!WY %in% c(2007)) %>% 
                            pull(burrows_prior) %>% 
                            scale(center = TRUE, scale = TRUE) %>% 
                            attr(., 'scaled:center'),
                          scale = obsdat %>% filter(!WY %in% c(2007)) %>% 
                            pull(burrows_prior) %>% 
                            scale(center = TRUE, scale = TRUE) %>% 
                            attr(., 'scaled:scale'))
    
    predvalues_backtransform = inputdat$predmatrix2 %>%  
      mutate(index = c(1:nrow(.))) %>%
      pivot_longer(-index) %>%
      left_join(scale_params, by = c('name' = 'cov')) %>%
      mutate(value = value * scale + center)
    
    predvalues = MCMCvis::MCMCsummary(mod, params = 'mu.ypred2', HPD = TRUE) %>% 
      as_tibble(rownames = 'var') %>% 
      set_names(c('var', 'mean', 'sd', 'lower', 'upper', 'Rhat', 'n.eff')) %>% 
      mutate(index = c(1:11),
             predictor = rep('N', 11)) %>% 
      left_join(predvalues_backtransform, by = c('index', 'predictor' = 'name')) 
    
  } else {
    nm = inputdat$predmatrix %>% names()
    
    scale_params = purrr::map_df(nm %>% rlang::set_names(),
                                 function(x) {
                                   tibble(center = inputdat$Bpredictors %>% pull(x) %>%
                                            attr(., 'scaled:center'),
                                          scale = inputdat$Bpredictors %>% pull(x) %>%
                                            attr(., 'scaled:scale'))
                                 }, .id = 'cov') 
    
    predvalues_backtransform = inputdat$predmatrix %>%  
      mutate(index = c(1:nrow(.))) %>%
      pivot_longer(-index) %>%
      left_join(scale_params, by = c('name' = 'cov')) %>%
      mutate(value = value * scale + center)
  
    if (logflow) {
      predvalues_backtransform = predvalues_backtransform %>% 
        mutate(value = case_when(grepl('flow_', name) ~ exp(value)-1,
                                 TRUE ~ value)) %>%
        select(index, name, value) %>%
        arrange(name, index)
    } else {
      predvalues_backtransform = predvalues_backtransform %>% 
        select(index, name, value) %>%
        arrange(name, index)
    }
    
    predvalues = MCMCvis::MCMCsummary(mod, params = 'mu.ypred', HPD = TRUE) %>% 
      as_tibble(rownames = 'var') %>% 
      set_names(c('var', 'mean', 'sd', 'lower', 'upper', 'Rhat', 'n.eff')) %>% 
      mutate(index = rep(c(1:11), length(nm)),
             predictor = rep(nm, each = 11)) %>% 
      left_join(predvalues_backtransform, by = c('index', 'predictor' = 'name')) 
  }
  
  return(predvalues)
}

plot_cov_relationships = function(predvalues, obsdat, ylim = c(-.55, .55)) {
  obsdat = obsdat %>% 
    pivot_longer(-pgr, names_to = 'cov', values_to = 'observed')
  
  ggplot(data = predvalues, aes(value, mean)) + facet_wrap(~predictor, scales = 'free') +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
    geom_hline(aes(yintercept = 0), linetype = 'dashed') +
    geom_line() + 
    geom_point(data = obsdat, aes(observed, pgr), color = 'red') +
    # scale_x_log10() + 
    ylim(ylim) +
    #labs(title = 'A', x = 'cfs (millions)', y = 'per-capita growth rate') +
    theme_classic() + 
    labs(x = NULL, y = 'per-capita growth rate')
  
  # purrr::map2(
  #   predlist,
  #   obslist,
  #   .f = function(x, y) {
  #     ggplot(data = x, aes(value, median)) +
  #       geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
  #       geom_hline(aes(yintercept = 0), linetype = 'dashed') +
  #       geom_line() + 
  #       geom_point(data = y, aes(observed, pgr), color = 'red') +
  #       # scale_x_log10() + 
  #       ylim(ylim) +
  #       #labs(title = 'A', x = 'cfs (millions)', y = 'per-capita growth rate') +
  #       theme_classic()
  #   }
  #   )
  
  # a = ggplot(predvalues %>% filter(cov == 'flow_concurrent'), aes(value, median)) +
  #   geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
  #   geom_hline(aes(yintercept = 0), linetype = 'dashed') +
  #   geom_line() + 
  #   geom_point(data = obsdat, aes(currentflow/1000000, pgr), color = 'red') +
  #   # scale_x_log10() + 
  #   ylim(ylim) +
  #   labs(title = 'A', x = 'cfs (millions)', y = 'per-capita growth rate') +
  #   theme_classic()
  # 
  # b = ggplot(predvalues %>% filter(cov == 'flow_prior'), aes(value, median)) +
  #   geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
  #   geom_hline(aes(yintercept = 0), linetype = 'dashed') +
  #   geom_line() + 
  #   geom_point(data = obsdat, aes(laggedflow/1000000, pgr), color = 'red') +
  #   # scale_x_log10() +
  #   ylim(ylim) +
  #   labs(title = 'B', x = 'cfs (millions)', y = NULL) +
  #   theme_classic()
  # 
  # c = ggplot(predvalues %>% filter(cov == 'flow_prior3'), aes(value, median)) +
  #   geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
  #   geom_hline(aes(yintercept = 0), linetype = 'dashed') +
  #   geom_line() + 
  #   geom_point(data = obsdat, aes(laggedflow/1000000, pgr), color = 'red') +
  #   # scale_x_log10() +
  #   ylim(ylim) +
  #   labs(title = 'C', x = 'cfs (millions)', y = NULL) +
  #   theme_classic()
  # 
  # c = ggplot(predvalues %>% filter(cov == 'mean.pdsi'), aes(value, median)) +
  #   geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
  #   geom_hline(aes(yintercept = 0), linetype = 'dashed') +
  #   geom_line() + 
  #   geom_point(data = obsdat, aes(pdsi, pgr), color = 'red') +
  #   labs(title = 'C', x = 'PDSI', y = 'per-capita growth rate') +
  #   ylim(ylim) +
  #   theme_classic()
  # 
  # d = ggplot(predvalues %>% filter(cov == 'N'), aes(value, median)) +
  #   geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
  #   geom_hline(aes(yintercept = 0), linetype = 'dashed') +
  #   geom_line() + 
  #   geom_point(data = obsdat, aes(burrows, pgr), color = 'red') +
  #   labs(title = 'D', x = 'burrows', y = NULL) +
  #   ylim(ylim) +
  #   theme_classic()
  # 
  # e = ggplot(predvalues %>% filter(cov == 't'), aes(value, median)) +
  #   geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'gray80') + 
  #   geom_hline(aes(yintercept = 0), linetype = 'dashed') +
  #   geom_line() + 
  #   geom_point(data = obsdat, aes(WY, pgr), color = 'red') +
  #   labs(title = 'E', x = 'year', y = 'per-capita growth rate') +
  #   ylim(ylim) +
  #   theme_classic()
  # 
  # library(patchwork)
  # a + b + c + d + e + plot_layout(nrow = 3, byrow = TRUE)

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
