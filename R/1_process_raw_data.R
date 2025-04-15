# load packages and functions
source('R/packages.R')
source('R/functions_process_raw_data.R')

# RESPONSE DATA--------
# BANS survey data: annual burrow counts between river miles 144-243 (roughly
# Colusa to Red Bluff); counts are the average of two separate counters'
# totals, so sometimes not an integer; survey protocols changed ~1999, so
# analysis focuses on data from 1999 onward

birddat = read_csv('data/BANS_burrow_counts.csv', col_types = cols()) %>% 
  # check for / fill in any missing years (gap in 2006)
  complete(year = seq(min(year), max(year))) %>%
  # add new data for 2021 and 2022
  bind_rows(tibble(year = c(2021, 2022, 2023),
                   burrows = c(20299, 11738, 7836)))

write_csv(birddat, 'data/birddat.csv')
# Note: no count data from 2006, so no agr/pgr calculated for 2005-06 or 2006-07
# intervals

ggplot(birddat, aes(year, burrows)) + 
  geom_line() + geom_point() + ylim(0, NA) +
  geom_smooth()
ggplot(birddat |> filter(year >= 1999), aes(year, burrows)) + 
  geom_line() + geom_point() + ylim(0, NA) +
  geom_smooth()
ggplot(birddat |> filter(year >= 1999), aes(year, burrows)) + 
  geom_line() + geom_point() + ylim(0, NA) +
  geom_smooth(method = 'lm')
ggplot(birddat |> filter(year >= 1999), aes(year, log(burrows))) + 
  geom_line() + geom_point() + #ylim(0, NA) +
  geom_smooth(method = 'lm')

# long-term trend: (1999-2023)
# log-normal because burrow counts are not all integers:
m = lm(log(burrows) ~ year, birddat |> filter(year >= 1999) |> drop_na())
summary(m)
(exp(m$coefficients['year'])-1)*100 #-1.49% per year; p = 0.027


# COVARIATES---------
# Focus on effects of environmental conditions on reproductive success, and
# subsequent growth to the following year; consider also density-dependence

# 1. Nest habitat quantity & quality-----------
# total daily mean flow above a threshold is proportional to stream power and
# extent of erosion, increasing the quantity and quality of nest habitat;
# thresholds considered: 14000 cfs and 50000 cfs; 14-15,000 cited as a
# threshold above which erosion occurs; 50,000 cfs cited as causing widespread
# erosion throughout the system (Larsen et al. 2006; BANS-TAC 2013)

# > data available via: 
# https://cdec.water.ca.gov/dynamicapp/wsSensorData
# or from package CDECRetrieve

# relevant gauges (N to S): 
# - Vina Bridge (VIN): data available from 1993-01-01
# - Hamilton City (HMC): data available from 2001-01-01
# - Ord Ferry (ORD): data available from 1993-01-01
# - Butte City Gauge (BTC): data available from 1998-01-30

vin = cdec_query('VIN', '41', 'D', start_date = '1993-01-01') #broken
hmc = cdec_query('HMC', '41', 'D', start_date = '1993-01-01') #actual start: 2001-01-04
ord = cdec_query('ORD', '41', 'D', start_date = '1993-01-01')
btc = cdec_query('BTC', '41', 'D', start_date = '1993-01-01') #now only starting 1998-01-01

flowdat_raw = bind_rows(vin, hmc,ord,btc)
write_csv(flowdat_raw, 'data/flowdat_raw.csv') #archive raw input data
  
flowdat_format = flowdat_raw |> 
  mutate(month = format(datetime, '%m') %>% as.numeric(),
         year = format(datetime, '%Y') %>% as.numeric(),
         WY = if_else(month < 10, year, year + 1),
         date = format(datetime, '%Y-%m-%d') |> as.Date()) |> 
  group_by(location_id, WY) %>%
  mutate(n = length(parameter_value[!is.na(parameter_value)])) %>%
  ungroup()

flowdat_format %>% select(location_id, WY, n) %>% distinct() %>% filter(n < 360) %>% 
  print(n = Inf)
# >> lots of stations have less than 365 days of data in a WY

# calculate mean value across stations each day
flowdat_mean = flowdat_format %>% 
  group_by(WY, month, date) %>% 
  summarize(nstations = length(parameter_value[!is.na(parameter_value)]), #note sometimes only 3 stations available per time step
            mflow = mean(parameter_value, na.rm = TRUE),
            .groups = 'drop')

flowdat_mean %>% filter(nstations == 0) %>% arrange(date) %>% print(n = Inf)
# many days in fall 1995 through winter 1996, plus a few days fall 2013, otherwise sporadic

flowdat_format |> 
  filter(date >= '2022-10-01' & date <= '2023-09-30') |> 
  ggplot(aes(date, parameter_value)) + 
  geom_line(aes(color = location_id)) +
  geom_line(data = flowdat_mean |> 
              filter(date >= '2022-10-01' & date <= '2023-09-30'), 
            aes(date, mflow))

# calculate sum by water year
flowdat_sum = flowdat_mean %>% 
  group_by(WY) %>% 
  summarize(ndays = length(mflow[!is.na(mflow)]),
            flowtotal = sum(mflow, na.rm = TRUE),
            ndays_threshold = length(mflow[mflow > 15000 & !is.na(mflow)]),
            flowtotal_threshold = sum(mflow[mflow > 15000], na.rm = TRUE),
            .groups = 'drop')

flowdat_sum %>% filter(ndays < 365)
# as of 2024-04-11, after averaging across stations, 9 WY have less than 365
# days of data (excluding 2024); but all are >=336 (>90% of days), and since
# 1997, all are >=360 (>98%)
  
write_csv(flowdat_sum, 'data/flowdat.csv')

# >> Note: originally explored seasonal flow totals, but annual totals highly
# correlated with rainy season totals, and breeding season totals very low

# characterize typical timing of annual flows (for inclusion in methods/study
# area description)
flowdat_month = flowdat_mean |> 
  filter(WY >= 1999 & WY < 2024) |> 
  mutate(mflow14 = if_else(mflow > 14000, mflow-14000, 0),
         mflow15 = if_else(mflow > 15000, mflow-15000, 0),
         month_name = factor(month.abb[month],
                             levels = month.abb[c(10:12,1:9)])) |> 
  group_by(WY, month_name, month) |> 
  summarize(mflow = sum(mflow, na.rm = TRUE),
            mflow14 = sum(mflow14, na.rm = TRUE),
            mflow15 = sum(mflow15, na.rm = TRUE),
            .groups = 'drop') |> 
  group_by(WY) |> 
  mutate(mflow_cum = cumsum(mflow),
         mflow_tot = sum(mflow),
         month_prop = mflow_cum/mflow_tot,
         mflow14_cum = cumsum(mflow14),
         mflow14_tot = sum(mflow14),
         month_prop14 = if_else(mflow14_tot > 0, 
                                mflow14_cum/mflow14_tot,
                                NA),
         mflow15_cum = cumsum(mflow15),
         mflow15_tot = sum(mflow15),
         month_prop15 = if_else(mflow15_tot > 0, 
                                mflow15_cum/mflow15_tot,
                                NA))

ggplot(flowdat_month, aes(month_name, mflow/1000000, group = as.factor(WY))) + 
  geom_line(aes(color = as.factor(WY)))
ggplot(flowdat_month, aes(month_name, mflow_cum/1000000, group = as.factor(WY))) + 
  geom_line(aes(color = as.factor(WY)))
ggplot(flowdat_month, aes(month_name, month_prop, group = as.factor(WY))) + 
  geom_line(aes(color = as.factor(WY)))

ggplot(flowdat_month, aes(month_name, mflow14/1000000, group = as.factor(WY))) + 
  geom_line(aes(color = as.factor(WY)))
ggplot(flowdat_month, aes(month_name, mflow14_cum/1000000, group = as.factor(WY))) + 
  geom_line(aes(color = as.factor(WY)))
ggplot(flowdat_month, aes(month_name, month_prop14, group = as.factor(WY))) + 
  geom_line(aes(color = as.factor(WY)))


# what proportion of cumulative total prior to Apr 1?
flowdat_month |> 
  select(WY, month_name, month, month_prop, month_prop14, month_prop15) |> 
  filter(month == 3) |> 
  pivot_longer(month_prop:month_prop15) |> 
  group_by(name) |> 
  summarize(mean = mean(value, na.rm = TRUE),
            min = min(value, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE))

# means:
# month_prop   0.547
# month_prop14 0.880
# month_prop15 0.883

# 2. Breeding season conditions---------
# loss of foraging habitat for BANS (and more broadly, aerial insectivores) is
# likely an important factor in pop declines; drought/primary productivity may
# play an important role in nestling survival

# drought indices from NOAA for Sacramento Valley area:
# >> code to pull straight from NCDC FTP website:

droughtdat = purrr::map_df(c('pdsidv', 'zndxdv'),
                           ~get_drought_indices(datname = .x)) %>%
  cleanup_droughtdat(statecode = '04', climdivcode = '02')

write_csv(droughtdat, 'data/droughtdat.csv')

# 3. Riprap---------

# data points from Adam Henderson for total length, Red Bluff to Colusa
riprap = tibble(WY = c(2002, 2012, 2018),
                length.m = c(81134.76, 83429.76, 86761.23)) |> 
  # incorporate older estimates from figure: (1935 = 0, but unsure when first riprap started)
  bind_rows(
    tibble(WY = c(1962, 1978, 1987),
          length.m = c(18, 29, 40) * 1609.344)) |> 
  arrange(WY)
ggplot(riprap, aes(WY, length.m/1000)) + geom_smooth(method = 'glm') + geom_point() 

m1 = lm(length.m ~ WY, riprap) 
summary(m1) #avg increase of 1km per year, but maybe decelerating?

riprap_approx = tibble(
  WY = c(1997:2024),
  riprap_approx = approx(x = riprap$WY, y = riprap$length.m, 
                         xout = c(1997:2024), rule = 2)$y,
  riprap_spline = spline(x = riprap$WY, y = riprap$length.m, 
                         xout = c(1997:2024),
                         method = 'natural')$y) |> 
  left_join(riprap |> select(WY, observed = length.m), by = 'WY')

ggplot(riprap_approx, aes(WY, observed)) + geom_point() + 
  geom_line(aes(y = riprap_approx), color = 'blue') +
  geom_line(aes(y = riprap_spline), color = 'red') + ylim(0, NA)

write_csv(riprap_approx, 'data/riprapdat.csv')
