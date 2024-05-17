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
    theme_classic()
  
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
    theme_classic()
  
  library(patchwork)
  a/b + plot_layout(axes = 'collect', guides = 'collect') +
    plot_annotation(tag_levels = 'A')
}

plot_observed_predicted = function(mod, obsdat, params = c('N', 'R'), 
                                   years = c(1999:2023), exclude99 = TRUE, 
                                   colorval = c('observed' = 'red', 'model' = 'black'),
                                   fillval = c('model' = 'gray70'), ...) {
  
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
  
  a = ggplot(preddat |> filter(var == 'N'), aes(WY, mean/1000)) + 
    geom_ribbon(aes(ymin = lci/1000, ymax = uci/1000, fill = source), alpha = 0.5) +
    geom_line(aes(color = source)) + 
    geom_point(aes(color = source, shape = source)) + 
    scale_color_manual(values = colorval) +
    scale_fill_manual(values = fillval) +
    theme_bw() + 
    scale_x_continuous(breaks = seq(1999, 2023, 2),
                       minor_breaks = seq(2000, 2022, 2),
                       labels = NULL,
                       guide = guide_axis(minor.ticks = TRUE)) +
    scale_y_continuous(limits = c(0, 30), expand = c(0, 0)) +
    theme(panel.grid = element_blank(), 
          legend.position = 'none',
          axis.minor.ticks.length = rel(0.5)) +
    labs(y = 'Burrow count (thousands)', x = NULL)
  
  b = ggplot(preddat |> filter(var == 'R'), aes(WY, mean)) + 
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = source), alpha = 0.5) +
    geom_line(aes(color = source)) + 
    geom_point(aes(color = source, shape = source)) + 
    scale_color_manual(values = colorval) +
    scale_fill_manual(values = fillval) +
    theme_bw() + 
    scale_x_continuous(breaks = seq(1999, 2023, 2),
                       minor_breaks = seq(2000, 2022, 2),
                       labels = seq(1999, 2023, 2),
                       guide = guide_axis(minor.ticks = TRUE)) +
    scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0)) +
    theme(panel.grid = element_blank(), 
          legend.position = 'none',
          axis.minor.ticks.length = rel(0.5)) +
    labs(y = 'Annual growth rate (%)', x = NULL) +
    geom_hline(aes(yintercept = 1), linetype = 'dashed')
  
  library(patchwork)
  a/b + plot_layout(axes = 'collect', guides = 'collect') +
    plot_annotation(tag_levels = 'A')
}

compile_partial_effects = function(mod, param = 'R.pred', covariates, predicted) {
  
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
    
  MCMCvis::MCMCsummary(mod, param = param, HPD = TRUE) |> 
    as_tibble(rownames = 'var', .name_repair = 'universal') |> 
    rename(lci = `..95._HPDL`, uci = `..95._HPDU`) |> 
    mutate(var = gsub('\\[.*$', '', var),
           varname = rep(c(names(predicted)), each = 11),
           index = rep(c(1:dim(predicted)[1]), 
                       times = dim(predicted)[2])) |> 
    left_join(backtransform, by = c('index', 'varname')) |> 
    left_join(modsum |> select(varname, pg0 = `p>0`), by = 'varname') |> 
    mutate(varname = gsub('log.', '', varname))
    
  
}

plot_partial_effects = function(dat, obsdat, 
                                var.order = c('flowt', 'flowt1', 'drought1', 'N', 'WY'),
                                labels = c('A', 'B', 'C', 'D', 'E'),
                                ylim = c(0, 2),
                                year_label = NULL) {
  
  dat = dat |> 
    mutate(problab = paste0('P{\u03b2 > 0} = ', format(pg0, nsmall = 2)),
           problab = if_else(pg0 >= 0.9 | pg0 <= 0.1, paste0(problab, '*'), problab))
  
  obsdat1 = obsdat |> select(any_of(var.order), N = burrows_orig, mean = agr) |> 
    mutate(N = lag(N)) |> filter(WY > 1999) |> 
    pivot_longer(-c('mean', 'WY'), names_to = 'varname', values_to = 'value_orig') |> 
    drop_na()
  
  default.plot = list(
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = 'gray80'),
    geom_line(),
    geom_hline(aes(yintercept = 1), linetype = 'dashed'),
    scale_y_continuous(limits = ylim, expand = c(0, 0)),
    geom_text(aes(x = -Inf, y = Inf, label = problab), 
              size = 3, hjust = -.1, vjust = 1.2, fontface = 'plain'),
    theme_bw(),
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0)))
  
  a = dat |> filter(varname == 'flowt') |> ggplot(aes(value_orig, mean)) + 
    default.plot +
    geom_point(data = obsdat1 |> filter(varname == 'flowt'), color = 'black') + 
    scale_x_continuous(limits = c(1.5, 8.5), breaks = seq(0, 10, 2)) +
    labs(x = 'cfs (millions)', y = 'Annual growth rate (%)', title = 'A')
  
  b = dat |> filter(varname == 'flowt1') |> ggplot(aes(value_orig, mean)) + 
    default.plot +
    geom_point(data = obsdat1 |> filter(varname == 'flowt1'), color = 'black') +
    scale_x_continuous(limits = c(1.5, 8.5), breaks = seq(0, 10, 2)) +
    labs(x = 'cfs (t-1, millions)', y = 'Annual growth rate (%)', title = 'B')
  
  c = dat |> filter(varname == 'drought1') |> ggplot(aes(value_orig, mean)) + 
    default.plot +
    geom_point(data = obsdat1 |> filter(varname == 'drought1'), color = 'black') + 
    scale_x_continuous(limits = c(-6, 4), breaks = seq(-6, 4, 2)) +
    labs(x = 'PDSI (t-1)', y = 'Annual growth rate (%)', title = 'C')
  
  d = dat |> filter(varname == 'N') |> ggplot(aes(value_orig/1000, mean)) + 
    default.plot +
    geom_point(data = obsdat1 |> filter(varname == 'N'), color = 'black') + 
    scale_x_continuous(limits = c(10, 21), breaks = seq(10, 25, 5)) +
    labs(x = 'Burrow count (t-1, thousands)', y = 'Annual growth rate (%)',
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
  
  if (!is.null(year_label)) {
    library(ggrepel)
    a = a + geom_text_repel(data = obsdat1 |> 
                        filter(varname == 'flowt' & WY %in% year_label),
                      aes(label = WY), size = 3, box.padding = 0.1)
    b = b + geom_text_repel(data = obsdat1 |> 
                        filter(varname == 'flowt1' & WY %in% year_label),
                      aes(label = WY), size = 3, box.padding = 0.1)
    c = c + geom_text_repel(data = obsdat1 |> 
                        filter(varname == 'drought1' & WY %in% year_label),
                      aes(label = WY), size = 3, box.padding = 0.1)
    d = d + geom_text_repel(data = obsdat1 |> 
                        filter(varname == 'N' & WY %in% year_label),
                      aes(label = WY), size = 3, box.padding = 0.1)
    # does not apply to plot E
  }
  
  # blanklabelplot = ggplot() + labs(y = "Annual growth rate (%)", x = NULL) + 
  #   theme_classic() + guides(x = "none", y = "none")
  
  library(patchwork)
  wrap_plots(a,b,c,d, ncol = 2, byrow = TRUE) +
    plot_layout(axis_titles = 'collect', axes = 'collect')
}