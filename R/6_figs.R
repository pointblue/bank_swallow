# Manuscript figures

# load packages and functions
source('R/packages.R')
source('R/functions_figs.R')

modeldat = read_csv('data/modeldat.csv') |> filter(WY >= 1999)
load('output/model_results.RData')
model_stats = read_csv('output/model_stats.csv')

# Fig 3: observed vs. predicted-----------
# growth rates and burrow counts
plot_observed_predicted(mod_res_final, obsdat = modeldat, HPD = TRUE, pg0 = TRUE)
ggsave(filename = 'fig/observed_predicted.png', width = 6.5, height = 5, units = 'in')

# Fig 4: covariates-----------

plot_covariates(modeldat, minyear = 1999)
ggsave(filename = 'fig/covariates.png', width = 6.5, height = 5, units = 'in')


# Fig 5: effect size plot----------

# model_stats |> filter(grepl('beta', var)) |> 
#   mutate(varname = gsub('log.', '', varname),
#          across(c(mean, sd, `95%_HPDL`, `95%_HPDU`), 
#                 ~if_else(varname == 'N', .*1000, .))) |> 
#   plot_effects(type = 'bar')

model_stats |> filter(grepl('beta', var)) |> 
  mutate(varname = gsub('log.', '', varname)) |> 
  plot_effects(type = 'distribution', mod = mod_res_final)
ggsave(filename = 'fig/effect_size.png', width = 5.5, height = 4, units = 'in')

# partial effects plots---------
# original predictor values from inputdat1
partial_effects = compile_partial_effects(mod_res_final, param = 'R.pred',
                                          covariates = inputdat_final$predictors, 
                                          predicted = inputdat_final$predmatrix)

plot_partial_effects(partial_effects, obsdat = modeldat, ylim = c(0.5, 1.7),
                     year_label = c(2001, 2010, 2015, 2017, 2018, 2022, 2023), 
                     pg0 = FALSE)
ggsave(filename = 'fig/partial_effects.png', width = 5.5, height = 7.5, units = 'in')


# rip rap plot-----------
riprapdat_full = readxl::read_excel('data/riprap_reach2and3allwithoffchannel_w new yrs.xls',
                                    skip = 2) |> 
  drop_na() |> 
  group_by(INBYDATE) |> 
  summarize(LENGTH_M = sum(LENGTH_M), .groups = 'drop') |>
  mutate(YEAR = as.numeric(INBYDATE),
         CUMULATIVE_M = cumsum(LENGTH_M)) |> 
  filter(YEAR < 2002) |> #assume some missing at the tail lend of this data set
  # add more recent totals:
  bind_rows(
    tibble(YEAR = c(2002, 2012, 2018),
           CUMULATIVE_M = c(81134.76, 83429.76, 86761.23))
  )

ggplot(riprapdat_full, aes(YEAR, CUMULATIVE_M/1000)) + 
  #geom_smooth(method = 'glm') + 
  geom_line() + 
  geom_point() +
  labs(x = NULL, y = 'Length of riprap (km)') +
  scale_x_continuous(breaks = seq(1940, 2020, by = 10),
                     minor_breaks = seq(1935, 2023, by = 5),
                     guide = guide_axis(minor.ticks = TRUE)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.ticks.length = unit(5, 'points'),
        axis.minor.ticks.length = rel(0.5))
ggsave(filename = 'fig/riprap_timeseries.png', width = 5.5, height = 3, units = 'in')


