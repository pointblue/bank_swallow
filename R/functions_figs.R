plot_observed_predicted = function(mod, obsdat, params = c('N', 'R'), 
                                   years = c(1999:2023), exclude99 = TRUE, 
                                   colorval = c('observed' = 'red', 'model' = 'black'),
                                   fillval = c('model' = 'gray70'), 
                                   Rmethod = c('ratio', '%'), 
                                   ylim1 = c(0, 30), ylim2 = c(0, 2.5),
                                   ...) {
  
  preddat = MCMCvis::MCMCsummary(mod, params = params, ...) |> 
    as_tibble(rownames = 'var', .name_repair = 'universal') |> 
    rename(lci = `..95._HPDL`, uci = `..95._HPDU`) |> 
    mutate(var = gsub('\\[.*$', '', var),
           WY = rep(years, length(params)),
           source = 'model') |> 
    # add observed data
    bind_rows(obsdat |> select(WY, N = burrows_orig, R = agr) |> 
                mutate(source = 'observed') |> 
                pivot_longer(N:R, names_to = 'var', values_to = 'mean'))
  
  if (exclude99) {
    preddat = preddat |> filter(!(var == 'R' & WY == '1999'))
  }
  
  if (Rmethod == '%') {
    preddat = preddat |> 
      mutate(across(c(mean, lci, uci), 
                    ~if_else(var == 'R', (.x - 1) * 100, .x)))
  }
  a = ggplot(preddat |> filter(var == 'N'), aes(WY, mean/1000)) + 
    geom_ribbon(aes(ymin = lci/1000, ymax = uci/1000, fill = source), alpha = 0.5) +
    geom_line(aes(color = source)) + 
    geom_point(aes(color = source, shape = source)) + 
    scale_color_manual(values = colorval) +
    scale_fill_manual(values = fillval) +
    scale_x_continuous(limits = c(1999, 2023),
                       breaks = seq(1999, 2023, 2),
                       minor_breaks = seq(2000, 2022, 2),
                       labels = NULL,
                       #labels = seq(1999, 2023, 2),
                       guide = guide_axis(minor.ticks = TRUE)) +
    scale_y_continuous(limits = ylim1, expand = c(0, 0)) +
    theme_classic() + 
    theme(panel.grid = element_blank(), 
          legend.position = 'none',
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          axis.ticks.length = unit(5, 'points'),
          axis.minor.ticks.length = rel(0.5)) +
    labs(y = 'Burrow count (thousands)', x = NULL)
  
  b = ggplot(preddat |> filter(var == 'R'), aes(WY, mean)) + 
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = source), alpha = 0.5) +
    geom_line(aes(color = source)) + 
    geom_point(aes(color = source, shape = source)) + 
    scale_color_manual(values = colorval) +
    scale_fill_manual(values = fillval) +
    theme_classic() + 
    scale_x_continuous(limits = c(1999, 2023),
                       breaks = seq(1999, 2023, 2),
                       minor_breaks = seq(2000, 2022, 2),
                       labels = seq(1999, 2023, 2),
                       guide = guide_axis(minor.ticks = TRUE)) +
    scale_y_continuous(expand = c(0, 0), limits = ylim2) +
    theme(panel.grid = element_blank(), 
          legend.position = 'none',
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          axis.ticks.length = unit(5, 'points'),
          axis.minor.ticks.length = rel(0.5)) +
    labs(y = 'Annual growth rate (%)', x = NULL) +
    geom_hline(aes(yintercept = 1), linetype = 'dashed')
  
  library(patchwork)
  a/b + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A')
}

plot_covariates = function(df, minyear = 1999) {
  
  a = ggplot(df |> filter(WY >= minyear), aes(WY, flowt)) + 
    geom_line() + geom_point() +
    labs(x = NULL, y = 'cfs (millions)') + 
    scale_y_continuous(breaks = seq(0, 10, 2), 
                       limits = c(0, 10)) +
    scale_x_continuous(breaks = seq(1999, 2023, 2),
                       minor_breaks = seq(2000, 2022, 2),
                       labels = NULL,
                       guide = guide_axis(minor.ticks = TRUE)) +
    theme_classic() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          axis.ticks.length = unit(5, 'points'),
          axis.minor.ticks.length = rel(0.5))
  
  b = ggplot(df |> filter(WY >= minyear), aes(WY, drought)) + 
    geom_line() + geom_point() +
    labs(x = NULL, y = 'PDSI') + 
    geom_hline(aes(yintercept = -2), linetype = 'dashed') +
    geom_hline(aes(yintercept = 2), linetype = 'dashed') +
    scale_x_continuous(breaks = seq(1999, 2023, 2),
                       minor_breaks = seq(2000, 2022, 2),
                       labels = seq(1999, 2023, 2),
                       guide = guide_axis(minor.ticks = TRUE)) +
    scale_y_continuous(limits = c(-6,6), 
                       breaks = seq(-6, 6, 2)) +
    theme_classic() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          axis.ticks.length = unit(5, 'points'),
          axis.minor.ticks.length = rel(0.5))
  
  library(patchwork)
  a/b + plot_layout(axes = 'collect', guides = 'collect') +
    plot_annotation(tag_levels = 'A')
}

plot_effects = function(dat,
                        type = c('bar', 'distribution'),
                        mod = NULL,
                        varnames = c('flowt', 'flowt1', 'drought1', 'N', 'WY')) {
  
  dat = dat |>  
    select(varname, mean, lcl = `95%_HPDL`, ucl = `95%_HPDU`, pg0 = `p>0`) |> 
    mutate(problab = paste0('P{\u03b2 > 0} = ', format(pg0, nsmall = 2)),
           problab = if_else(round(pg0,2) >= 0.9 | round(pg0, 2) <= 0.1, 
                             paste0(problab, '*'), problab),
           varname = factor(varname, 
                            levels = c('flowt', 'flowt1', 'drought1', 'N', 'WY')),
           varlabels = recode(
             varname,
             flowt = 'Total stream flow\n(t, cfs, millions)', 
             flowt1 = 'Prior total stream flow\n(t-1, cfs, millions)', 
             drought1 = 'Breeding season drought\n(PDSI, t-1)', 
             N = 'Prior burrow count\n(t-1)', 
             WY = 'Water year')) |> 
    arrange(varname)
  
  if (type == 'bar') {
    ggplot(dat, aes(varlabels, mean)) +
      geom_hline(aes(yintercept = 0), linetype = 'dashed') +
      geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.25) +
      geom_point() + 
      scale_x_discrete(labels = parse(text = levels(df$varlabels))) +
      geom_text(aes(y = ucl + 0.01, label = label), size = 4) +
      labs(x = NULL, y = 'effect size') + 
      theme_classic() +
      theme(panel.grid = element_blank(),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12))
    
  } else if (type == 'distribution') {
    
    samples_df = MCMCvis::MCMCchains(mod, params = 'beta', ISB = TRUE) |> 
      as_tibble() |> 
      # original order of betas
      set_names(c('flowt', 'flowt1', 'drought1', 'WY', 'N')) |> 
      pivot_longer(everything(), names_to = 'varname') |> 
      mutate(
        varname = factor(varname, 
                         levels = rev(c('flowt', 'flowt1', 'drought1', 'N', 'WY'))),
        varlabels = recode(
          varname,
          flowt = 'Total stream flow\n(t, cfs, millions)', 
          flowt1 = 'Prior total stream flow\n(t-1, cfs, millions)', 
          drought1 = 'Breeding season drought\n(PDSI, t-1)', 
          N = 'Prior burrow count\n(t-1)', 
          WY = 'Water year'))
    
    ggplot(samples_df, aes(y = varlabels)) +
      ggridges::geom_density_ridges(
        aes(x = value), scale = 0.95, rel_min_height = 0.01,
        color = NA, fill = 'gray70') +
      scale_y_discrete(#labels = function(x) str_wrap(x, width = 25),
        expand = expansion(add = 0.25)) +
      geom_errorbar(data = dat, aes(xmin = lcl, xmax = ucl), width = 0.1) +
      geom_point(data = dat, aes(x = mean)) +
      geom_vline(data = dat, aes(xintercept = 0), linetype = 'dashed') +
      xlim(-0.4, 0.4) +
      labs(x = 'Effect size', y = NULL) +
      theme_classic() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 8),
            axis.title = element_text(size = 10)) +
      geom_text(data = dat |> filter(pg0 >= 0.5), 
                aes(x = 0.01, label = paste0(round(pg0 * 100, 0), '%')), 
                hjust = 0, vjust = 0, nudge_y = 0.08, size = 2.5) +
      geom_text(data = dat |> filter(pg0 < 0.5),
                aes(x = -0.01, label = paste0(round((1 - pg0) * 100, 0), '%')), 
                hjust = 1, vjust = 0, nudge_y = 0.08, size = 2.5)
    
  }
  
}

compile_partial_effects = function(mod, param = 'R.pred', 
                                   covariates, predicted,
                                   Rmethod = c('ratio', '%')) {
  
  scale_params = tibble(
    varname = covariates |> names(),
    center = covariates |> purrr::map_dbl(attr_getter("scaled:center")),
    scale = covariates |> purrr::map_dbl(attr_getter("scaled:scale"))
  )
  
  backtransform = predicted |> select(-N) |>  
    mutate(index = c(1:dim(predicted)[1])) |> 
    pivot_longer(-index, names_to = 'varname') |> 
    left_join(scale_params, by = 'varname') |> 
    mutate(
      value_orig = value * scale + center,
      value_orig = if_else(grepl('flow', varname), exp(value_orig), value_orig)) |> 
    arrange(varname, index) |> 
    # add values of N for density dependence
    bind_rows(
      predicted |> select(value_orig = N) |> 
        mutate(index = c(1:11), varname = 'N')
    )
  
  # prob > 0
  modsum = MCMCvis::MCMCsummary(mod, params = 'beta', HPD = TRUE, pg0 = TRUE) |> 
    as_tibble(rownames = 'var') |> 
    mutate(varname = names(predicted))
    
  res <- MCMCvis::MCMCsummary(mod, param = param, HPD = TRUE) |> 
    as_tibble(rownames = 'var', .name_repair = 'universal') |> 
    rename(lci = `..95._HPDL`, uci = `..95._HPDU`) |> 
    mutate(var = gsub('\\[.*$', '', var),
           varname = rep(c(names(predicted)), each = 11),
           index = rep(c(1:dim(predicted)[1]), 
                       times = dim(predicted)[2])) |> 
    left_join(backtransform, by = c('index', 'varname')) |> 
    left_join(modsum |> select(varname, pg0 = `p>0`), by = 'varname') |> 
    mutate(varname = gsub('log.', '', varname))
    
  if (Rmethod == '%') {
    res = res |> 
      mutate(across(c(mean, lci, uci), ~(.x - 1) * 100))
  }
  return(res)
}

plot_partial_effects = function(dat, obsdat, 
                                #var.order = c('flowt', 'flowt1', 'drought1', 'N', 'WY'),
                                var.label = c('Total streamflow (t, cfs, millions)', 
                                              'Prior total streamflow (t-1, cfs, millions)',
                                              'Prior breeding season drought (t-1, PDSI)',
                                              'Prior burrow count (t-1, thousands)'),
                                labels = c('A', 'B', 'C', 'D', 'E'),
                                ylim = c(0, 2),
                                pg0 = TRUE,
                                year_label = NULL) {
  
  problab = dat |> select(varname, pg0) |> distinct() |> 
    mutate(problab = paste0('P{\u03b2 > 0} = ', format(pg0, nsmall = 2)),
           problab = if_else(pg0 >= 0.9 | pg0 <= 0.1, paste0(problab, '*'), problab))
  
  obsdat1 = obsdat |> 
    select(WY, flowt, flowt1, drought1, N = burrows_orig, mean = agr) |> 
    mutate(N = lag(N)/1000) |> filter(WY > 1999) |> 
    pivot_longer(-c('mean', 'WY'), names_to = 'varname', values_to = 'value_orig') |> 
    drop_na()
  
  default.plot = list(
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = 'gray80'),
    geom_line(),
    geom_hline(aes(yintercept = 1), linetype = 'dashed'),
    scale_y_continuous(limits = ylim, expand = c(0, 0)),
    theme_classic(),
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 10)))
  
  a = dat |> filter(varname == 'flowt') |> ggplot(aes(value_orig, mean)) + 
    default.plot +
    geom_point(data = obsdat1 |> filter(varname == 'flowt'), color = 'black') + 
    scale_x_continuous(limits = c(1, 9), breaks = seq(0, 10, 1)) +
    labs(x = var.label[1], y = 'Annual growth rate (%)', title = 'A')
  
  b = dat |> filter(varname == 'flowt1') |> ggplot(aes(value_orig, mean)) + 
    default.plot +
    geom_point(data = obsdat1 |> filter(varname == 'flowt1'), color = 'black') +
    scale_x_continuous(limits = c(1, 9), breaks = seq(0, 10, 1)) +
    labs(x = var.label[2], y = 'Annual growth rate (%)', title = 'B')
  
  c = dat |> filter(varname == 'drought1') |> ggplot(aes(value_orig, mean)) + 
    default.plot +
    geom_point(data = obsdat1 |> filter(varname == 'drought1'), color = 'black') + 
    scale_x_continuous(limits = c(-6, 4), breaks = seq(-6, 4, 2)) +
    labs(x = var.label[3], y = 'Annual growth rate (%)', title = 'C')
  
  d = dat |> filter(varname == 'N') |> ggplot(aes(value_orig, mean)) + 
    default.plot +
    geom_point(data = obsdat1 |> filter(varname == 'N'), color = 'black') + 
    scale_x_continuous(limits = c(10, 21), breaks = seq(10, 25, 1)) +
    labs(x = var.label[4], y = 'Annual growth rate (%)',
         title = 'D')
  
  # >> duplicates observed_predicted plot
  # e = dat |> filter(varname == 'WY') |> ggplot(aes(value_orig, mean)) + 
  #   default.plot +
  #   geom_point(data = obsdat1 |> select(value_orig = WY, mean) |> distinct(), 
  #              color = 'black') + 
  #   geom_line(data = obsdat1 |> select(value_orig = WY, mean) |> distinct()) + 
  #   scale_x_continuous(limits = c(1999, 2023),
  #                      breaks = seq(1999, 2023, 4),
  #                      minor_breaks = seq(2000, 2023, 1),
  #                      labels = seq(1999, 2023, 4),
  #                      guide = guide_axis(minor.ticks = TRUE)) +
  #   labs(x = NULL, y = 'Annual growth rate (%)', title = 'E')
  
  if (pg0) {
    a = a + geom_text(data = problab |> filter(varname == 'flowt'), 
                      aes(x = -Inf, y = Inf, label = problab), 
                      size = 3, hjust = -.1, vjust = 1.5, fontface = 'plain')
    b = b + geom_text(data = problab |> filter(varname == 'flowt1'), 
                      aes(x = -Inf, y = Inf, label = problab), 
                      size = 3, hjust = -.1, vjust = 1.5, fontface = 'plain')
    c = c + geom_text(data = problab |> filter(varname == 'drought1'), 
                      aes(x = -Inf, y = Inf, label = problab), 
                      size = 3, hjust = -.1, vjust = 1.5, fontface = 'plain')
    d = d + geom_text(data = problab |> filter(varname == 'N'), 
                      aes(x = -Inf, y = Inf, label = problab), 
                      size = 3, hjust = -.1, vjust = 1.5, fontface = 'plain')
  }
  
  if (!is.null(year_label)) {
    library(ggrepel)
    a = a + 
      geom_point(data = obsdat1 |> 
                   filter(varname == 'flowt' & WY %in% year_label),
                 col = 'red') +
      geom_text_repel(data = obsdat1 |> 
                        filter(varname == 'flowt' & WY %in% year_label),
                      aes(label = WY), size = 3, box.padding = 0.15,
                      col = 'red')
    b = b + 
      geom_point(data = obsdat1 |> 
                   filter(varname == 'flowt1' & WY %in% year_label),
                 col = 'red') +
      geom_text_repel(data = obsdat1 |> 
                        filter(varname == 'flowt1' & WY %in% year_label),
                      aes(label = WY), size = 3, box.padding = 0.15,
                      col = 'red')
    c = c + 
      geom_point(data = obsdat1 |> 
                   filter(varname == 'drought1' & WY %in% year_label),
                 col = 'red') +
      geom_text_repel(data = obsdat1 |> 
                        filter(varname == 'drought1' & WY %in% year_label),
                      aes(label = WY), size = 3, box.padding = 0.15,
                      #point.size = NA, 
                      col = 'red')
    d = d + 
      geom_point(data = obsdat1 |> 
                   filter(varname == 'N' & WY %in% year_label),
                 col = 'red') +
      geom_text_repel(data = obsdat1 |> 
                        filter(varname == 'N' & WY %in% year_label),
                      aes(label = WY), size = 3, box.padding = 0.15,
                      col = 'red')
    # does not apply to plot E
  }
  
  # blanklabelplot = ggplot() + labs(y = "Annual growth rate (%)", x = NULL) + 
  #   theme_classic() + guides(x = "none", y = "none")
  
  library(patchwork)
  wrap_plots(a,b,c,d, ncol = 1, byrow = TRUE) 
}

