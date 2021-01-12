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

  # 4 stations are very well correlated; may be able to use just one
  
  # calculate indicator of stream power as the sum of daily mean flow > 425cms
  # (~15000 cfs) by water year, based on VIN
  VIN_cumpower = mflowdat %>% 
    filter(STATION_ID == 'VIN' & VALUE > 15000) %>% 
    group_by(WY) %>% 
    summarize(power = sum(VALUE),
              .groups = 'drop')
  
  # 2. Primary productivity---------
  
  # 3. Winter weather------
  
  # Analysis-----------
  
  
  # data = generate_data(),
  # model = fit_model(data)
)
