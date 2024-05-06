# Compile covariates and check for correlations

# load packages and functions
source('R/packages.R')
source('R/functions.R')

birddat = read_csv('data/birddat.csv')

# Generate model predictors----------

## burrows and lag burrows-----
# > do not include burrow data prior to 1999, or growth rate from 1998 to 1999,
# since methods changed; so earliest legit pgr and lagged burrow count is in
# 2000

birddat = read_csv('data/birddat.csv') |> 
  select(WY = year, burrowst = burrows) |> 
  mutate(burrowst1 = lag(burrowst)) |> 
  # calculate annual per-capita growth rate as the change in burrow count 
  # >> *from the previous year* (this is a change from prior version!)
  mutate(agr = burrowst/burrowst1,
         pgr = (burrowst - burrowst1)/burrowst1) |> 
  filter(WY > 1999)

# check correlations:
corrplot::corrplot.mixed(
  cor(birddat |> drop_na(), method = 'spearman'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# burrowst moderately correlated with pgr (0.54)
# burrowst1 moderately negatively correlated with pgr (-0.48)

ggplot(birddat, aes(agr)) + geom_density() + 
  geom_vline(xintercept = 1, col = 'red')

mean(birddat$pgr, na.rm = TRUE) # mean = -0.017
median(birddat$pgr, na.rm = TRUE) # median = 0.005
var(birddat$pgr, na.rm = TRUE) # var = 0.055

## flow & lagged flow -------

# > decided to include both total flow and thresholded version even though
# highly correlated with total; will be helpful to report that model results
# were virtually identical, but focus on presenting total in figures (easier to
# understand)

# > decided to keep lagged flow years separate as well (rather than cumulative
# total); easier to plot/compare effect size and see if effects diminish

flowdat = read_csv('data/flowdat.csv') |> 
  filter(WY != 2024) |> # drop incomplete water year
  select(WY, flowt = flowtotal, flow14t = flowtotal_threshold) |> 
  # # extend to 2024--2026
  # bind_rows(tibble(WY = c(2024:2026))) |> 
  # find lag values
  mutate(flowt1 = lag(flowt),
         flowt2 = lag(flowt, n = 2),
         flowt3 = lag(flowt, n = 3),
         flow14t1 = lag(flow14t),
         flow14t2 = lag(flow14t, n = 2),
         flow14t3 = lag(flow14t, n = 3),
         # convert to cfs (millions)
         across(matches('^flow'), ~./1000000)) |> 
  filter(WY > 1999) # match to birddat

# check correlations:
corrplot::corrplot.mixed(
  cor(flowdat |> drop_na(), method = 'spearman'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# > flowt and flow14t highly correlated (>0.93)
# > all flow data moderately negatively correlated with WY
# > flow and lagged flow values not correlated

## drought--------
# Note: monthly PDSI values do not depict drought on timescales < ~12 months;
# it's designed to be highly correlated month to month to account for "land
# memory" (https://climatedataguide.ucar.edu/climate-data/palmer-drought-severity-index-pdsi)
# --> the Z index is an alternative that is more responsive to short term drought conditions

# drought_sum = bind_rows(
#   # seasonal means
#   read_csv('data/droughtdat.csv') |>
#     filter(WY >= 1995) |>
#     filter(!is.na(season)) |>
#     group_by(index, WY, season) |>
#     summarize(value = mean(value),
#               .groups = 'drop'),
#   # annual means
#   read_csv('data/droughtdat.csv') |>
#     filter(WY >= 1995) |>
#     group_by(index, WY) |>
#     summarize(value = mean(value),
#               .groups = 'drop') |>
#     mutate(season = 'WY')) |>
#   pivot_wider(names_from = c('index', 'season'))
# 
# # correlations:
# corrplot::corrplot.mixed(
#   cor(drought_sum %>% filter(WY >= 2000 & WY <= 2022) %>% drop_na(), method = 'pearson'),
#   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# # zndx and pdsi strongly correlated within season (breeding = 0.90; rainy =
# # 0.87; WY = 0.94)
# # zndx breeding and rainy season metrics not strongly correlated (0.51), but
# # each is strongly correlated with WY (0.87, 0.86)
# # pdsi breeding and rainy strongly correlated (0.73), and both strongly 
# # correlated with WY (0.90, 0.95)
# # no correlation with WY for any

# decision to focus on mean pdsi during the breeding season (Apr-Aug): strongly
# correlated pdsi during the rainy season, average pdsi for the year, and zndx
# within a season

droughtdat = read_csv('data/droughtdat.csv') |>
  filter(index == 'pdsidv' & season == 'breeding' & WY != 2024) |>
  group_by(WY) |>
  summarize(drought = mean(value), .groups = 'drop') |> 
  mutate(drought1 = lag(drought)) |> 
  filter(WY > 1999)

corrplot::corrplot.mixed(
  cor(droughtdat |> drop_na(), method = 'pearson'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# no strong correlations

## riprap---------

riprapdat = read_csv('data/riprapdat.csv') |> 
  select(WY, riprap = riprap_approx) |> 
  # convert m to km
  mutate(riprap = riprap/1000,
         riprap1 = lag(riprap)) |> 
  filter(WY > 1999 & WY < 2024)


# Compile predictors----------------

modeldat = list(birddat,
                flowdat,
                droughtdat,
                riprapdat) |> 
  purrr::reduce(full_join)
write_csv(modeldat, 'data/modeldat.csv')

# check for correlations among predictors
corrplot::corrplot.mixed(
  cor(modeldat %>% select(!starts_with('flow14t')) %>% drop_na(), method = 'spearman'),
  order = 'hclust', tl.col = 'black', tl.pos = 'lt')
# strongest: 
# >> WY and riprap perfectly correlated (of course)
# - flowt and drought (0.65)
# - flowt1 and: drought1 (0.65), burrowst (0.67), pgr (0.68)
# - pgr and burrowst (0.54), burrowst1 (-0.48)

# ggplot(modeldat, aes(flow14t1)) + geom_density()
# ggplot(modeldat, aes(log(flow14t1))) + geom_density()
# ggplot(modeldat, aes(flowt1)) + geom_density()
# ggplot(modeldat, aes(log(flowt1))) + geom_density()
# 
# mean(modeldat$pgr, na.rm = TRUE) # mean = -0.017
# var(modeldat$pgr, na.rm = TRUE) # var = 0.055

# simulate = function(sigma) {
#   draws = rnorm(10000, mean = 0, sd = sigma)
#   hist(draws, freq = F, main = paste('random draws for var =', sigma),
#        breaks = 20)
#   #hist(plogis(draws), freq = F, main = 'inverse logit of draws', breaks = 20)
# }

# Summary stats------------
# note: flow data is in millions of cfs; riprap is in km
cov_stats = modeldat |> select(WY, flowt:last_col()) |> 
  pivot_longer(-WY) |> 
  group_by(name) |> 
  summarize(mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            min = min(value, na.rm = TRUE),
            .groups = 'drop')
write_csv(cov_stats, 'output/cov_stats.csv')
