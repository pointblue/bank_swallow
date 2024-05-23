# Manuscript figures

# load packages and functions
source('R/packages.R')
source('R/functions_figs.R')

modeldat = read_csv('data/modeldat.csv') |> filter(WY >= 1999)
load('output/model_results.RData')
model_stats = read_csv('output/model_stats.csv')

# covariates-----------

plot_covariates(modeldat, minyear = 1999)
ggsave(filename = 'fig/covariates.png', width = 6.5, height = 5, units = 'in')

# observed vs. predicted-----------
# growth rates and burrow counts
plot_observed_predicted(mod_res_final, obsdat = modeldat, HPD = TRUE, pg0 = TRUE)
ggsave(filename = 'fig/observed_predicted.png', width = 6.5, height = 5, units = 'in')

# partial effects plots---------
# original predictor values from inputdat1
partial_effects = compile_partial_effects(mod_res_final, param = 'R.pred',
                                          covariates = inputdat_final$predictors, 
                                          predicted = inputdat_final$predmatrix)

plot_partial_effects(partial_effects, obsdat = modeldat, ylim = c(0.5, 1.7),
                     year_label = c(2001, 2010, 2015, 2017, 2018, 2022, 2023), 
                     pg0 = FALSE)
ggsave(filename = 'fig/partial_effects.png', width = 5.5, height = 8.5, units = 'in')

# effect size plot----------

model_stats |> filter(grepl('beta', var)) |> 
  mutate(varname = gsub('log.', '', varname),
         across(c(mean, sd, `95%_HPDL`, `95%_HPDU`), 
                ~if_else(varname == 'N', .*1000, .))) |> 
  plot_effects(type = 'bar')

model_stats |> filter(grepl('beta', var)) |> 
  mutate(varname = gsub('log.', '', varname)) |> 
  plot_effects(type = 'distribution', mod = mod_res_final)
ggsave(filename = 'fig/effect_size.png', width = 5.5, height = 4, units = 'in')

# rip rap plot-----------
riprap = tibble(WY = c(2002, 2012, 2018),
                length.m = c(81134.76, 83429.76, 86761.23)) |> 
  # incorporate older estimates from figure: (1935 = 0, but unsure when first riprap started)
  bind_rows(
    tibble(WY = c(1962, 1978, 1987),
           length.m = c(18, 29, 40) * 1609.344)) |> 
  arrange(WY)
ggplot(riprap, aes(WY, length.m/1000)) + geom_smooth(method = 'glm') + geom_point() 

