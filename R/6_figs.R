# Manuscript figures

# load packages and functions
source('R/packages.R')
source('R/functions_figs.R')

modeldat = read_csv('data/modeldat.csv') |> filter(WY >= 1999)
load('output/mod1_final_20240514.Rdata')

# covariates-----------

plot_covariates(modeldat, minyear = 1999)
ggsave(filename = 'fig/covariates.png', width = 6.5, height = 5, units = 'in')

# observed vs. predicted-----------
# growth rates and burrow counts
plot_observed_predicted(mod1_res, obsdat = modeldat, HPD = TRUE, pg0 = TRUE)
ggsave(filename = 'fig/observed_predicted.png', width = 6.5, height = 5, units = 'in')

# partial effects plots---------
# original predictor values from inputdat1
partial_effects = compile_partial_effects(mod1_res, param = 'R.pred',
                                          covariates = inputdat1$predictors, 
                                          predicted = inputdat1$predmatrix)

plot_partial_effects(partial_effects, obsdat = modeldat)
ggsave(filename = 'fig/partial_effects.png', width = 6.5, height = 7, units = 'in')

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
