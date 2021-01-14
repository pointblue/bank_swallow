# Drake plan reference: https://books.ropensci.org/drake/plans.html

plan <- drake_plan(
  # update readme for github repo https://github.com/pointblue/bank_swallow
  readme_page = render_Rmd(file_in("Rmd/README.Rmd"),
                           file_out("README.md")),
  
  # Response variable--------
  # BANS survey data from Greg
  
  birddat = readxl::read_excel(
    file_in('data/from Greg/SR R2 and R3 BANS ACTIVE BURROWS and OUTFLOWS 1986-2020.xlsx'),
    sheet = 2, range="A5:E38", 
    col_types = c('numeric', 'skip', 'skip', 'skip', 'numeric'),
    col_names = c('year', 'burrows')) %>% 
    # calculate annual % change
    mutate(agr = (burrows - lag(burrows)) / lag(burrows) * 100),
  
  # # EXPLORE:
  # ggplot(birddat, aes(year, burrows)) + geom_line() + geom_point()
  # # long-term cycles evident? steep decline through 1995, increase through 2001,
  # # slower decline through 2015, then a rebound
  # 
  # ggplot(birddat, aes(year, agr)) + geom_line() + geom_point()
  # # dramatic annual % changes, alternating + and - through 1996, then switches
  # # to a pattern of steep increase, multi-year decrease through 2010, more
  # # variable since
  
  # Predictors---------
  # Consider long-term trends, as well as direct and indirect (lagged) effects
  # of (1) stream flow (erosion) on nest site quality; (2) primary productivity
  # and food availability during the breeding season; and (3) over-winter
  # habitat conditions, drought, and food availability
  
  # 1. Stream flow-----------
  # Larsen et al. (2006) used daily mean stream flow from USGS, sum total above
  # flow threshold of 425cms during a given water year (Oct 1-Sept 30)

  # no available API/token, so manually downloaded hourly and mean daily and
  # flow data for each gauge from:
  # https://cdec.water.ca.gov/dynamicapp/wsSensorData
  
  # relevant gauges (N to S): 
  # - Vina Bridge (VIN): mean daily flow data available from 1993-01-01; hourly flow from 1985-10-01
  # - Hamilton City (HMC): mean daily flow available from 2001-01-01; hourly from 1991-06-19
  # - Ord Ferry (ORD): mean daily flow available from 1993-01-01; hourly from 1985-10-01
  # - Butte City Gauge (BTC): available from 1993-01-01
  # - Moulton Weir: no mean daily flow data? hourly available from 1998-01-28

  mflowdat = compile_waterdat(dir = file_in('data/CDEC'), type = 'MFLOW') %>% 
    mutate(VALUE = if_else(VALUE == '---' | VALUE < 0, NA_character_, VALUE),
           VALUE = as.numeric(VALUE)) %>% 
    mutate_at(vars(STATION_ID, DURATION, SENSOR_TYPE, UNITS), as.factor) %>% 
    # add water year
    mutate(month = format(`DATE TIME`, '%m') %>% as.numeric(),
           year = format(`DATE TIME`, '%Y') %>% as.numeric(),
           WY = if_else(month < 10, year, year + 1)),
  
  # EXPLORE FLOW FOR ANY WATER YEAR:
  # mflowdat %>% filter(WY == 2006) %>% 
  #   ggplot(aes(`DATE TIME`, VALUE, color = STATION_ID)) + geom_line() +
  #   geom_hline(aes(yintercept = 15000), color = 'red')
  
  # calculate indicator of stream power as the sum of daily mean flow > 425cms
  # (~15000 cfs) by water year, for each gauge
  streampower = mflowdat %>% group_by(STATION_ID, WY) %>% add_count() %>% 
    ungroup() %>% filter(n >= 365) %>% 
    # add threshold value
    mutate(season = if_else(month >= 10 | month <= 3, 'rainy', 'dry'),
           tVALUE = VALUE - 15000,
           tVALUE = if_else(tVALUE < 0, 0, tVALUE)) %>%  
    group_by(STATION_ID, WY, season) %>% 
    summarize(totalflow = sum(VALUE, na.rm = TRUE),
              streampower = sum(tVALUE, na.rm = TRUE),
              .groups = 'drop'),
  # --> consider seasonal flows (rather than total annual)
  
  # EXPLORE STREAM POWER ACROSS YEARS & STATIONS:
  # ggplot(streampower, aes(WY, power/1000000, color = STATION_ID)) + 
  #   geom_line() + geom_point()
  # streampower %>% filter(power < 250000 & STATION_ID == 'VIN') %>% arrange(power)
  # # lowest years in 2020, 1994, 2014 (followed by 2018, 2012, 2007, 2008, 2001, 2009)
  # 
  # # check correlations across stations
  # corrplot::corrplot.mixed(
  #   cor(streampower %>%
  #         pivot_wider(names_from = 'STATION_ID', values_from = 'power') %>%
  #         filter(WY > 2001) %>% select(-WY, -n),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # # --> 4 stations are very strongly correlated; may be able to use just one to
  # # represent all (especially since HMC data not available until WY2002)
  #
  # RELATIONSHIP WITH TOTAL FLOW: very strong (+0.97)
  # ggplot(streampower, aes(streampower, totalflow)) + geom_point()
  # cor(streampower %>% select(totalflow:streampower))
  
  # EXPLORE DIRECT RELATIONSHIP TO BURROWS/AGR:
  # left_join(streampower %>% filter(STATION_ID == 'VIN'),
  #           birddat, by = c('WY' = 'year')) %>%
  #   # ggplot(aes(streampower, burrows)) + geom_point() + ylim(0, 30000) # weak negative
  #   # ggplot(aes(totalflow, burrows)) + geom_point() + ylim(0, 30000) # weak negative
  #   # ggplot(aes(streampower, agr)) + geom_point() + ylim(-60, 60) #weak negative
  #   ggplot(aes(totalflow, agr)) + geom_point() + ylim(-60, 60) #weak negative
  
  # EXPLORE LAGGED RELATIONSHIP TO BURROWS/AGR:
  # left_join(streampower %>% filter(STATION_ID == 'VIN'),
  #           birddat %>% mutate(WY = year - 1), by = 'WY') %>%
  #   # ggplot(aes(streampower, burrows)) + geom_point() + ylim(0, 30000) # pretty flat
  #   # ggplot(aes(totalflow, burrows)) + geom_point() + ylim(0, 30000) # pretty flat
  #   # ggplot(aes(streampower, agr)) + geom_point() + ylim(-60, 60) # positive
  #   # ggplot(aes(totalflow, agr)) + geom_point() + ylim(-60, 60) # positive
  # 
  # left_join(streampower %>% filter(STATION_ID == 'VIN'),
  #           birddat %>% mutate(WY = year - 2), by = 'WY') %>% 
  #   ggplot(aes(power, burrows)) + geom_point() + ylim(0, 25000) +
  #   xlab('streampower, t-2') + ylab('burrows, t')
  # # no obvious relationship
  
  
  # 2. Primary productivity---------
  
  # consider NDVI or EVI from MODIS
  
  #--> figure out how to add bounding box or shapefile, and what output
  # projection and resolution we need before testing this over a short timeframe!
  processMODIS = MODIStsp::MODIStsp(
    gui = FALSE, 
    out_folder = "data/MODIS", 
    selprod = "Vegetation_Indexes_16Days_250m (M*D13Q1)",
    bandsel = c("EVI", "NDVI"), 
    quality_bandsel = "QA_usef", 
    indexes_bandsel = NULL, 
    user = Sys.getenv('EARTHDATA_ID'),
    password = Sys.getenv('EARTHDATA_PASSWORD'),
    start_date = "2020.06.01", 
    end_date = "2020.06.15", 
    verbose = FALSE),
  
  # 3. Winter weather------
  
  # Lamanna et al. 2012 used El Nino as an indicator of winter survival,
  # specifically the ESPI index (El Nino Precipitation index), described here: 
  # https://psl.noaa.gov/enso/dashboard.html 
  
  espi = RCurl::getURL('http://eagle1.umd.edu/GPCP_ICDR/Data/ESPI.txt') %>% 
    read_table() %>% 
    select(year = YYYY, month = MM, ESPI),
  
  # # EXPLORE MONTHLY PATTERNS
  # espi %>% group_by(year) %>% filter(ESPI == max(ESPI)) %>% 
  #   ggplot(aes(month)) + geom_density(adjust = 0.25)
  # # annual peak can be in any month; most frequent is Jan, Feb, Apr
  # espi %>% group_by(year) %>% filter(ESPI == min(ESPI)) %>% 
  #   ggplot(aes(month)) + geom_density(adjust = 0.25)
  # # annual nadir can be in any month; most frequent is Feb, Apr (?)
  
  # consider also SOI: sign is opposite that of other indices, and noisier; no
  # API, so manually downloaded from:
  # https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/soi.long.data
  
  soi = read_table('data/soi.csv', skip = 1, n_max = 155,
                   col_names = c('year', '1', '2', '3', '4', '5', '6', '7', 
                                 '8', '9', '10', '11', '12'),
                   col_types = cols()) %>% 
    pivot_longer(-year, names_to = 'month', values_to = 'SOI') %>% 
    mutate(month = as.numeric(month),
           SOI = if_else(SOI < -99, NA_real_, SOI)),
  
  # # EXPLORE MONTHLY PATTERNS
  # soi %>% filter(year >= 1985) %>% 
  #   group_by(year) %>% filter(SOI == max(SOI)) %>%
  #   ggplot(aes(month)) + geom_density(adjust = 0.25)
  # # annual peak most frequently in Feb
  # soi %>% filter(year >= 1985) %>% 
  #   group_by(year) %>% filter(SOI == min(SOI)) %>%
  #   ggplot(aes(month)) + geom_density(adjust = 0.25)
  # # annual nadir most frequently in June
  
  # EXPLORE RELATIONSHIP BETWEEN INDICES:
  # inner_join(espi, soi) %>% ggplot(aes(ESPI, SOI)) + geom_point()
  # 
  # corrplot::corrplot.mixed(
  #   cor(inner_join(espi, soi) %>% select(ESPI, SOI) %>% filter(!is.na(SOI)), method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # # negatively correlated, but not super tightly (-0.7)
  # 
  # inner_join(espi, soi) %>% 
  #   mutate(date = as.Date(paste(year, month, '01', sep = '-'))) %>% 
  #   filter(year>=2000 & year <=2010) %>% 
  #   pivot_longer(ESPI:SOI) %>% 
  #   ggplot(aes(date, value, color = name)) + geom_line()
  
  
  # COMBINE THE TWO METRICS & SUMMARIZE
  # seasonal metrics assumed to be June-Oct "rainy" and otherwise "dry"; for now
  # explore "bioyear" as starting with start of rainy season in June, with the
  # same number as the following BANS breeding season
  climdat = inner_join(espi, soi, by = c('year', 'month')) %>% 
    mutate(season = if_else(month >= 6 & month <= 10, 'rainy', 'dry'),
           bioyear = if_else(month >= 6, year + 1, year)) %>% 
    group_by(bioyear, season) %>% 
    summarize(ESPI = mean(ESPI),
              SOI = mean(SOI),
              .groups = 'drop'),
  
  # ggplot(climdat, aes(bioyear, ESPI, color = season)) + geom_line()
  # ggplot(climdat, aes(bioyear, SOI, color = season)) + geom_line()
  # cor(climdat %>% select(bioyear:ESPI) %>% 
  #       pivot_wider(names_from = 'season', values_from = 'ESPI') %>% 
  #       select(-bioyear) %>% filter(complete.cases(.)))
  # cor(climdat %>% select(bioyear:season, SOI) %>% 
  #       pivot_wider(names_from = 'season', values_from = 'SOI') %>% 
  #       select(-bioyear) %>% filter(complete.cases(.)))
  # # ESPI and SOI during June-Oct rainy season strongly correlated with the
  # # subsequent dry season (Nov-May) (+0.78 for both metrics)
  
  # EXPLORE RELATIONSHIP WITH #BURROWS/AGR:
  # left_join(birddat, 
  #           climdat %>% filter(season == 'rainy'), 
  #           by = c('year' = 'bioyear')) %>% 
  # ggplot(aes(ESPI, burrows)) + geom_point() #weak negative?
  # ggplot(aes(SOI, burrows)) + geom_point() #very weak positive?
  # ggplot(aes(ESPI, agr)) + geom_point() #weak negative?
  # ggplot(aes(SOI, agr)) + geom_point() #weak positive?
  
  # LAGGED EFFECT?
  # left_join(birddat %>% mutate(bioyear = year - 1),
  #           climdat %>% filter(season == 'rainy'),
  #           by = 'bioyear') %>%
  # ggplot(aes(ESPI, burrows)) + geom_point() + ylim(0, 30000)
  # ggplot(aes(SOI, burrows)) + geom_point() + ylim(0, 30000)
  # ggplot(aes(ESPI, agr)) + geom_point() + ylim(-60, 60)
  # ggplot(aes(SOI, agr)) + geom_point() + ylim(-60, 60) 
  
  
  # Analysis-----------
  
  
  # data = generate_data(),
  # model = fit_model(data)
)
