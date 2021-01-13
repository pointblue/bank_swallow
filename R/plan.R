# Drake plan reference: https://books.ropensci.org/drake/plans.html

plan <- drake_plan(
  # update readme for github repo https://github.com/pointblue/bank_swallow
  readme_page = render_Rmd(file_in("Rmd/README.Rmd"),
                           file_out("README.md")),
  
  # Response variable--------
  # BANS survey data
  
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
  streampower = full_join(
    mflowdat %>% group_by(STATION_ID, WY) %>% count() %>% ungroup(),
    mflowdat %>% 
      mutate(VALUE = VALUE - 15000,
             VALUE = if_else(VALUE < 0, 0, VALUE)) %>%  
      group_by(STATION_ID, WY) %>% 
      summarize(power = sum(VALUE, na.rm = TRUE),
                .groups = 'drop'),
    by = c('WY', 'STATION_ID')) %>% 
    filter(n >= 365),
  
  # EXPLORE STREAM POWER ACROSS YEARS & STATIONS:
  # ggplot(streampower, aes(WY, power/1000000, color = STATION_ID)) + 
  #   geom_line() + geom_point()
  # streampower %>% filter(power < 250000 & STATION_ID == 'VIN') %>% arrange(power)
  # # lowest years in 2020, 1994, 2014 (followed by 2018, 2012, 2007, 2008, 2001, 2009)
  # 
  # # check correlations
  # corrplot::corrplot.mixed(
  #   cor(streampower %>%
  #         pivot_wider(names_from = 'STATION_ID', values_from = 'power') %>%
  #         filter(WY > 2001) %>% select(-WY, -n),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # # --> 4 stations are very strongly correlated; may be able to use just one to
  # # represent all (especially since HMC data not available until WY2002)
  
  # 2. Primary productivity---------
  
  # consider NDVI or EVI from MODIS
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
  
  # # EXPLORE RELATIONSHIP:
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
  

  # Analysis-----------
  
  
  # data = generate_data(),
  # model = fit_model(data)
)
