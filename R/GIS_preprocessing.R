## SCRIPT FOR PRE-PROCESSING GIS DATA (OUTSIDE OF MAIN WORKFLOW)
library(tidyverse)
library(sf)

# Extract BANS data from excel-----------
readxl::read_excel(
  file_in('data/BANS/SR R2 and R3 BANS ACTIVE BURROWS and OUTFLOWS 1986-2020.xlsx'),
  sheet = 2, range="A5:E38", 
  col_types = c('numeric', 'skip', 'skip', 'skip', 'numeric'),
  col_names = c('year', 'burrows')) %>% 
  write_csv('data/BANS_burrow_counts.csv')

# Primary study area (breeding season)-----------
# burrow locations:
shp = read_sf('GIS/ds6.shp') %>% filter(RIVMI >= 144 & RIVMI <= 243)
# ggplot(shp) + geom_sf()

# 10km buffer
shp_buff = shp %>% mutate(state = 'CA') %>%
  group_by(state) %>% summarize(.groups = 'drop') %>%
  st_buffer(dist = 10000) %>%
  write_sf('GIS/BANS_buffer_10km.shp')

shp_buff = shp %>% mutate(state = 'CA') %>%
  group_by(state) %>% summarize(.groups = 'drop') %>%
  st_buffer(dist = 2000) %>%
  write_sf('GIS/BANS_buffer_2km.shp')

# 1km buffer: "too complex" for AppEEARS
shp_buff_1k = shp %>% mutate(state = 'CA') %>%
  group_by(state) %>% summarize(.groups = 'drop') %>%
  st_buffer(dist = 1000) %>%
  write_sf('GIS/BANS_buffer_1km.shp')

# Possible wintering area-----------------
# from BOTW range maps
westhemi = c(xmax = -20, xmin = -140, ymin = -180, ymax = 180)
mx = c(xmax = -95, xmin = -110, ymin = 15, ymax = 30)

winter_range = read_sf('GIS/BOTW_BANS.shp') %>% 
  select(binomial, seasonal) %>% 
  filter(seasonal == 3) %>% 
  st_cast("POLYGON") %>% 
  st_make_valid() %>% 
  st_crop(westhemi)

mx_range = winter_range %>% st_crop(mx) %>% write_sf('GIS/BANS_winter_MX.shp')


# Vegetation & phenology data-----------
# Use AppEEARS to compile data and stats from within target areas.

# > breeding range-----------
read_csv('GIS/CA_2k/MOD13A3-006-Statistics.csv', col_types = cols()) %>% 
  mutate(metric = gsub('.*_', '', Dataset),
         year = as.numeric(format(Date, '%Y')),
         month = as.numeric(format(Date, '%m'))) %>% 
  select(metric, year, month, Count, Mean, Median, std = `Standard Deviation`) %>% 
  ggplot(aes(month, Median)) + 
  facet_wrap(~metric, scales = 'free_y') +
  geom_errorbar(aes(ymin = Median - std, ymax = Median + std)) +
  geom_point(aes(color = as.factor(year))) + 
  geom_line(aes(color = as.factor(year))) +
  geom_smooth(size = 3) +
  xlim(1, 7)

# variation by month
read_csv('GIS/CA_2k/MOD13A3-006-Statistics.csv', col_types = cols()) %>% 
  mutate(metric = gsub('.*_', '', Dataset),
         year = as.numeric(format(Date, '%Y')),
         month = as.numeric(format(Date, '%m'))) %>% 
  select(metric, year, month, Count, Mean, Median) %>% 
  filter(metric == 'EVI') %>% 
  ggplot(aes(Mean, as.factor(month))) + geom_violin()

# try calculating seasonal averages: "pre-breeding" (e.g. Mar-Apr) and "breeding" (May-Jun)
cadat = read_csv('GIS/CA_2k/MOD13A3-006-Statistics.csv', col_types = cols()) %>% 
  mutate(metric = gsub('.*_', '', Dataset),
         year = as.numeric(format(Date, '%Y')),
         month = as.numeric(format(Date, '%m'))) %>% 
         # season = case_when(month == 3 | month == 4 ~ 'MarApr',
         #                    month == 5 | month == 6 ~ 'MayJun',
         #                    TRUE ~ NA_character_)) %>% 
  # filter(!is.na(season)) %>% 
  filter(month >= 4 & month <= 6) %>% 
  select(metric, year, month, Count, Mean, Median, std = `Standard Deviation`) %>% 
  group_by(metric, year) %>% 
  summarize(n = length(Mean),
            Mean = mean(Mean),
            .groups = 'drop')

# compare EVI & NDVI
cadat %>% 
  pivot_wider(names_from = 'metric', values_from = 'Mean') %>% 
  ggplot(aes(EVI, NDVI)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = 'red') +
  geom_smooth() + ylim(0, .75) + xlim(0, .75) 
  # facet_wrap(~season)

# compare pre-breeding and breeding: no apparent relationship
cadat %>% select(-n) %>% 
  pivot_wider(names_from = 'season', values_from = 'Mean') %>% 
  ggplot(aes(MarApr, MayJun)) + geom_point() + 
  # geom_abline(aes(slope = 1, intercept = 0), color = 'red') +
  geom_smooth() + 
  # ylim(0, .8) + xlim(0, .8) + 
  facet_wrap(~metric)

# variation over years: increasing trend in MayJun EVI & NDVI
ggplot(cadat, aes(year, Mean)) + 
  facet_wrap(~metric, scales = 'free_y') +
  geom_point() + geom_smooth() + ylim(0, NA)



# > breeding range 10k---------
# try calculating seasonal averages: "pre-breeding" (e.g. Mar-Apr) and "breeding" (May-Jun)
cadat10k = read_csv('GIS/CA_10k/MOD13A3-006-Statistics.csv', col_types = cols()) %>% 
  mutate(metric = gsub('.*_', '', Dataset),
         year = as.numeric(format(Date, '%Y')),
         month = as.numeric(format(Date, '%m')),
         season = case_when(month == 3 | month == 4 ~ 'MarApr',
                            month == 5 | month == 6 ~ 'MayJun',
                            TRUE ~ NA_character_)) %>% 
  filter(!is.na(season)) %>% 
  select(metric, year, month, season, Count, Mean, Median, std = `Standard Deviation`) %>% 
  group_by(metric, year, season) %>% 
  summarize(n = length(Mean),
            Mean = mean(Mean),
            .groups = 'drop')

# compare EVI & NDVI
cadat10k %>% 
  pivot_wider(names_from = 'metric', values_from = 'Mean') %>% 
  ggplot(aes(EVI, NDVI)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = 'red') +
  geom_smooth() + ylim(0, .75) + xlim(0, .75) + 
  facet_wrap(~season)

# compare pre-breeding and breeding: no apparent relationship
cadat10k %>% select(-n) %>% 
  pivot_wider(names_from = 'season', values_from = 'Mean') %>% 
  ggplot(aes(MarApr, MayJun)) + geom_point() + 
  # geom_abline(aes(slope = 1, intercept = 0), color = 'red') +
  geom_smooth() + 
  # ylim(0, .8) + xlim(0, .8) + 
  facet_wrap(~metric)

# variation over time: slight increasing trend in MayJun EVI only; less extreme
# than 2k version
ggplot(cadat10k, aes(year, Mean, color = season)) + 
  facet_wrap(~metric, scales = 'free_y') +
  geom_point() + geom_smooth() + ylim(0, NA)



# > winter range-------------
# phenology indices for MX range
read_csv('GIS/MX/MCD12Q2-006-Statistics.csv', col_types = cols()) %>% 
  mutate(year = as.numeric(format(Date, '%Y'))) %>% 
  filter(Count > 20000 & Dataset != 'NumCycles') %>% #drop smaller second cycle of greenup
  select(year, Dataset, Mean, Median, std = `Standard Deviation`) %>% 
  filter(Dataset == 'Greenup') %>% 
  # conversion from days since 1970 to day of year:
  mutate(MedDate = if_else(Dataset %in% c('Greenup', 'Peak', 'Senescence'),
                           as.Date('1970-01-01') + Median,
                           as.Date('1970-01-01')),
         Mean = if_else(Dataset %in% c('Greenup', 'Peak', 'Senescence'),
                        as.numeric(format(as.Date('1970-01-01') + Mean, '%j')),
                        Mean),
         Median = if_else(Dataset %in% c('Greenup', 'Peak', 'Senescence'),
                          as.numeric(format(as.Date('1970-01-01') + Median, '%j')),
                          Median)) %>% 
  ggplot(aes(year, Median)) + facet_wrap(~Dataset, scales = 'free_y') + 
  geom_errorbar(aes(ymin = Median - std, ymax = Median + std)) +
  geom_point()

# Greenup starts late Apr through mid June
# Peak occurs late July-late Aug
# Senescence occurs late Aug through late Sept

# timing of greenup:
read_csv('GIS/MX/MCD12Q2-006-Statistics.csv', col_types = cols()) %>% 
  filter(Dataset == 'Greenup') %>% 
  mutate(year = as.numeric(format(Date, '%Y'))) %>% 
  filter(Count > 20000) %>% #drop smaller second cycle of greenup
  select(year, Date, Dataset, Median, UQ = `Upper Quartile`, 
         LQ = `Lower Quartile`) %>% 
  # conversion from days since 1970 to day of year:
  mutate_at(vars(Median:LQ),
            ~as.numeric(format(as.Date('1970-01-01') + ., '%j'))) %>% 
  pivot_longer(Median:LQ) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_point() + geom_line() + 
  geom_hline(aes(yintercept = 32)) + #Feb 1
  geom_hline(aes(yintercept = 120)) + # Apr 30
  # geom_hline(aes(yintercept = 151)) + #May 31
  ylim(0, NA)

# NDVI & EVI for MX range
read_csv('GIS/MX/MOD13A3-006-Statistics.csv', col_types = cols()) %>% 
  mutate(metric = gsub('.*_', '', Dataset),
         year = as.numeric(format(Date, '%Y')),
         month = as.numeric(format(Date, '%m'))) %>% 
  select(metric, year, month, Count, Mean, Median, std = `Standard Deviation`) %>% 
  ggplot(aes(month, Median, group = as.factor(year))) + 
  facet_wrap(~metric, scales = 'free_y') +
  geom_errorbar(aes(ymin = Median - std, ymax = Median + std)) +
  geom_point(aes(color = as.factor(year))) + 
  geom_line(aes(color = as.factor(year))) +
  xlim(1, 6)

# try calculating winter/spring averages
mxdat = read_csv('GIS/MX/MOD13A3-006-Statistics.csv', col_types = cols()) %>% 
  mutate(metric = gsub('.*_', '', Dataset),
         year = as.numeric(format(Date, '%Y')),
         month = as.numeric(format(Date, '%m'))) %>% 
  select(metric, year, month, Count, Mean, Median, std = `Standard Deviation`) %>% 
  filter(month >= 2 & month <= 4) %>% 
  group_by(metric, year) %>% 
  summarize(n = length(Mean),
            Mean = mean(Mean),
            .groups = 'drop')

# compare EVI & NDVI
mxdat %>% 
  pivot_wider(names_from = 'metric', values_from = 'Mean') %>% 
  ggplot(aes(EVI, NDVI)) + geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0), color = 'red') +
  geom_smooth() + ylim(0, .75) + xlim(0, .75)

# variation over time
ggplot(mxdat, aes(year, Mean)) + facet_wrap(~metric, scales = 'free_y') +
  geom_point() + geom_smooth() + ylim(0, 0.6)

# focus on breeding season, and early spring lead-up (Jan through July?)

# For test period, 2014-2020:
# - EVI MayJunJul strongly correlated with dry season PDSI and zndx (>0.8),
#     and annual PDSI (>0.7) and zndx (>.582)
# - EVI MarApr negatively correlated with all drought indices (?)
# - EVI JanApr negatively correlated with all drought indices, but relatively
#     weakly

# get rasters from MODIS, subset to area near BANS burrows:
# (do this step outside of drake plan, repeat for sets of years)

# 1. create a boundary shapefile


# 2. get EVI rasters from MODIS for buffer area, during breeding season (May-Jun)
MODIStsp::MODIStsp_get_prodnames()
MODIStsp::MODIStsp(
  gui = FALSE,
  out_folder = "GIS/MODIS",
  out_folder_mod = "C:/Users/kdybala/AppData/Local/Temp",
  selprod = "Vegetation Indexes_16Days_250m (M*D13Q1)",
  bandsel = "EVI", #don't need NDVI
  quality_bandsel = NULL,
  indexes_bandsel = NULL,
  sensor = 'both', #more data!
  user = Sys.getenv('EARTHDATA_ID'),
  password = Sys.getenv('EARTHDATA_PASSWORD'),
  start_date = "1999.05.01",
  end_date = "2020.06.30",
  download_range = 'Seasonal',
  spatmeth = 'file',
  spafile = 'GIS/BANS_buffer_1km.shp',
  delete_hdf = 'TRUE',
  ts_format = "R RasterStack",
  out_format = 'GTiff',
  compress = 'LZW',
  verbose = TRUE)