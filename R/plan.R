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
  
  # 3. Winter weather------
  
  # Analysis-----------
  
  
  # data = generate_data(),
  # model = fit_model(data)
)
