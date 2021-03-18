# Drake plan reference: https://books.ropensci.org/drake/plans.html

plan <- drake_plan(
  # update readme for github repo https://github.com/pointblue/bank_swallow
  readme_page = render_Rmd(file_in("Rmd/README.Rmd"),
                           file_out("README.md")),
  
  # RESPONSE--------
  # BANS survey data from Greg: annual burrow counts between river miles 144-243
  # (roughly Colusa to Red Bluff)
  
  # From Greg: Pre 1999 data have not be error proofed in the same way as latter
  # data, and yes some of the methods were different. For example surveys
  # weren’t completed in just a few days in one concerted effort like since
  # 1999. Instead they were done over the course of the summer, which would
  # certainly affect counts. I wasn’t involved in surveys back then.  I think
  # for the most part we should only analyze data from 1999 onwards.
  
  # > Burrow counts are the average of two separate counters' totals.
  
  birddat = readxl::read_excel(
    file_in('data/BANS/SR R2 and R3 BANS ACTIVE BURROWS and OUTFLOWS 1986-2020.xlsx'),
    sheet = 2, range="A5:E38", 
    col_types = c('numeric', 'skip', 'skip', 'skip', 'numeric'),
    col_names = c('year', 'burrows')) %>% 
    # check for / fill in any missing years
    complete(year = seq(min(year), max(year))) %>%
    # calculate annual change to the following year
    mutate(agr = lead(burrows)/burrows,
           pgr = (lead(burrows) - burrows)/burrows),
  # Note: no count data from 2006, so no agr/pgr calculated for 2005-06
  # interval, or 2006-07 interval

  # > any evidence for density dependence?
  # # when #burrows is higher (> ~14k), growth to the following year tends to be
  # # <1; relationship seems fairly linear, but not super strong
  # ggplot(birddat, aes(burrows, agr)) + geom_point() +
  #   geom_hline(aes(yintercept = 1), color = 'red') +
  #   geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'cs')) +
  #   scale_y_log10()
  # ggplot(birddat, aes(burrows, pgr)) + geom_point() +
  #   geom_hline(aes(yintercept = 0), color = 'red') +
  #   geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'cs'))
  
  # PREDICTORS---------
  # Consider long-term trends in population growth rate, as well as direct and
  # indirect (lagged) effects of environmental conditions on nest success,
  # survival, and immigration/recruitment
  
  # 1. Stream flow-----------
  # high flows cause bank erosion which creates new colony sites, but also
  # destroys existing colony sites; timing matters for impacts on the breeding
  # population: during winter = good for nest site availability (carrying
  # capacity, immigration/recruitment) and quality (nest success), during
  # breeding = bad (nest success)

  # Proposed metrics: total flow, flow above a certain cfs (related to "stream
  # power")
  
  # > river gauge data----------
  
  # - Larsen et al. (2006) used daily mean stream flow from USGS, sum total
  # above flow threshold of 425cms (~15,000 cfs) during a given water year (Oct
  # 1-Sept 30)
  # - BANS-TAC report indicates flows >50,000 cfs can cause widespread erosion
  # and loss of multiple colonies (but beneficial during nonbreeding season);
  # lower breeding season flows between 14,000-30,000 cfs can cause
  # local/partial colony collapses and elevated river stage can cause nest
  # failure; breeding season: Apr 1-Aug 31
  
  # no available API/token, so manually downloaded hourly and mean daily and
  # flow data for each gauge from:
  # https://cdec.water.ca.gov/dynamicapp/wsSensorData
  # -->(hourly flow data full of missing data points and times when flow was
  # above/below the rating table; assume mean daily flow data has accounted for
  # some of this?)

  # relevant gauges (N to S): 
  # - Vina Bridge (VIN): mean daily flow data available from 1993-01-01; hourly flow from 1985-10-01
  # - Hamilton City (HMC): mean daily flow available from 2001-01-01; hourly from 1991-06-19
  # - Ord Ferry (ORD): mean daily flow available from 1993-01-01; hourly from 1985-10-01
  # - Butte City Gauge (BTC): available from 1993-01-01
  # - Moulton Weir: no mean daily flow data? hourly available from 1998-01-28
  
  mflowdat = compile_waterdat(dir = file_in('data/CDEC'), type = 'MFLOW',
                              col_types = cols()) %>% 
    mutate(VALUE = if_else(VALUE == '---' | VALUE < 0, NA_character_, VALUE),
           VALUE = as.numeric(VALUE)) %>% 
    mutate_at(vars(STATION_ID, DURATION, SENSOR_TYPE, UNITS), as.factor) %>% 
    # add water year
    mutate(month = format(`DATE TIME`, '%m') %>% as.numeric(),
           year = format(`DATE TIME`, '%Y') %>% as.numeric(),
           WY = if_else(month < 10, year, year + 1)) %>% 
    # filter out incomplete annual data
    group_by(STATION_ID, WY) %>% 
    add_count() %>% 
    ungroup() %>% 
    filter(n >= 365) %>% 
    mutate(season = case_when(month >= 10 | month <= 3 ~ 'rainy',
                              month >= 4 & month <= 5 ~ 'breeding_early',
                              month >= 6 & month <= 8 ~ 'breeding_late',
                              TRUE ~ NA_character_),
           season = factor(season, 
                           levels = c('rainy', 'breeding_early', 
                                      'breeding_late')
                           ), #rainy comes first in a WY
           flow50 = if_else(VALUE > 50000, VALUE - 50000, 0), #daily flow above 50k cfs
           flow14 = if_else(VALUE > 14000, VALUE - 14000, 0) #daily flow above 14k cfs
    ),
  
  # calculate annual and seasonal flow totals per station
  mflow_tot = bind_rows(
    # annual totals
    mflowdat %>% 
      group_by(STATION_ID, WY) %>% 
      summarize(season = 'annual',
                n = length(VALUE[!is.na(VALUE)]),
                flowtot = sum(VALUE, na.rm = TRUE),
                flow50 = sum(flow50, na.rm = TRUE),
                flow14 = sum(flow14, na.rm = TRUE),
                .groups = 'drop') %>% 
      # check those with data for <90% of days each year
      mutate_at(vars(flowtot:flow14),
                ~if_else(n < 0.9*365, NA_real_, .)),
    # seasonal totals
    mflowdat %>% 
      group_by(STATION_ID, WY, season) %>% 
      summarize(n = length(VALUE[!is.na(VALUE)]),
                flowtot = sum(VALUE, na.rm = TRUE),
                flow50 = sum(flow50, na.rm = TRUE),
                flow14 = sum(flow14, na.rm = TRUE),
                .groups = 'drop') %>% 
      filter(!is.na(season)) %>% 
      # check those with data for <90% of days in each season
      mutate_at(vars(flowtot:flow14),
                ~case_when(season == 'breeding_early' & n < 0.9*61 ~ NA_real_,
                           season == 'breeding_late' & n < 0.9*61 ~ NA_real_,
                           season == 'rainy' & n < 0.9*182 ~ NA_real_,
                           TRUE ~ .)) %>% 
      select(-n)
  ),

  # ggplot(mflow_tot %>% filter(STATION_ID == 'VIN'),
  #        aes(WY, flowtot)) + geom_line() + geom_point() +
  #   geom_line(aes(WY, flow50), color = 'blue') +
  #   geom_line(aes(WY, flow14), color = 'red') +
  #   facet_wrap(~season, scales = 'free')
  # # decline in late breeding season flows over time (and high flows) over time
  # # --> less change in early breeding season flows
  
  # # check correlations across stations, seasons, and metrics
  # # > early breeding season: very strong correlations across stations & metrics
  # corrplot::corrplot.mixed(
  #   cor(mflow_tot %>% filter(WY >= 2000 & season == 'breeding_early') %>%
  #         pivot_wider(names_from = STATION_ID,
  #                     values_from = flowtot:flow14) %>%
  #         select(-WY, -season, -n) %>% drop_na(),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # 
  # # > late breeding season: zero flow50!  strong correlations across stations;
  # # little relationship between flow14 & flowtot
  # corrplot::corrplot.mixed(
  #   cor(mflow_tot %>% filter(WY >= 2000 & season == 'breeding_late') %>%
  #         pivot_wider(names_from = STATION_ID,
  #                     values_from = flowtot:flow14) %>%
  #         select(-starts_with('flow50')) %>% 
  #         select(-WY, -season, -n) %>% drop_na(),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # 
  # # > rainy season: very strong correlations across stations & all metrics
  # corrplot::corrplot.mixed(
  #   cor(mflow_tot %>% filter(WY >= 2000 & season == 'rainy') %>%
  #         pivot_wider(names_from = STATION_ID,
  #                     values_from = flowtot:flow14) %>%
  #         select(-WY, -season, -n) %>% drop_na(),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')

  
  # summarize average across stations
  mflow_sum = mflow_tot %>% 
    group_by(WY, season) %>% 
    summarize(n = length(flowtot[!is.na(flowtot)]), #note sometimes only 3 stations available per time step
              flowtot = mean(flowtot, na.rm = TRUE),
              flow50 = mean(flow50, na.rm = TRUE),
              flow14 = mean(flow14, na.rm = TRUE),
              .groups = 'drop') %>% 
    filter(WY >= 1997) %>% 
    select(-n) %>% #ranges 2-4 stations contributing data; mostly 4 from 2002 on
    pivot_wider(names_from = season, values_from = flowtot:flow14) %>% 
    # calculate running totals for current rainy season + 1-2 prior WYs
    mutate(flowtot_cum3 = flowtot_rainy + lag(flowtot_annual, 1) + 
             lag(flowtot_annual, 2),
             # zoo::rollsum(lag(flowtot_annual), 2, fill = NA, align = 'right'),
           flow50_cum3 = flow50_rainy + lag(flow50_annual, 1) + 
             lag(flow50_annual, 2),
             # zoo::rollsum(lag(flow50_annual), 2, fill = NA, align = 'right'),
           flow14_cum3 = flow14_rainy + lag(flow14_annual, 1) + 
             lag(flow14_annual, 2),
             # zoo::rollsum(lag(flow14_annual), 2, fill = NA, align = 'right')
           flowtot_cum2 = flowtot_rainy + lag(flowtot_annual, 1),
           flow50_cum2 = flow50_rainy + lag(flow50_annual, 1),
           flow14_cum2 = flow14_rainy + lag(flow14_annual, 1)
           ),
  
  # # check correlations:
  # corrplot::corrplot.mixed(
  #   cor(mflow_sum %>% select(-WY) %>%
  #         select(-(starts_with('flow50') & ends_with('late'))) %>%
  #         drop_na(), method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # # - strong correlations within cumulative metrics, but not with any
  # # seasonal/annual metrics
  # # - strong correlations among annual & rainy season metrics
  # # - strong correlations (>0.7) among breeding season metrics, but not really
  # # with rainy season or annual metrics

  
  # > unimpaired flow data ---------
  # Total "unimpaired" outflow and the "Sacramento River Index" or "4 River
  # Index" from: http://cdec.water.ca.gov/reportapp/javareports?name=WSIHIST
  
  # --> NOTE: this is "unimpaired flow", but it is strongly correlated with
  # river gauge data from 4 stations along reaches 2-3 (and these data go back
  # farther)
  
  sri = RCurl::getURL('http://cdec.water.ca.gov/reportapp/javareports?name=WSIHIST') %>% 
    strsplit("\\r\\n"),
  
  sri_sum = sri[[1]][416:538] %>%
    read_table(skip = 4,
               col_names = c('WY', 'Sac_Oct-Mar', 'Sac_Apr-Jul', 'Sac_WYsum',
                             'Sac_Index', 'Type_Sac',
                             'SJ_Oct-Mar', 'SJ_Apr-Jul', 'SJ_WYsum', 'SJ_Index',
                             'Type_SJ')) %>%
    select(WY, starts_with('Sac')) %>%
    rename(Runoff_OctMar = 'Sac_Oct-Mar',
           Runoff_AprJul = 'Sac_Apr-Jul',
           Runoff_WY = 'Sac_WYsum',
           SRI = 'Sac_Index'),
  
  
  # # check correlations:
  # corrplot::corrplot.mixed(
  #   cor(full_join(sri_sum %>% filter(WY >= 2000),
  #                 mflow_sum %>% select(WY, flow14_rainy, flow14_breeding, flow14_cum)) %>%
  #         select(-WY) %>% drop_na(),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # 
  # # very strong positive correlations between SRI, unimpaired runoff, and gauge
  # # data; cross-seasonal correlations not strong; cumulative flow data not
  # # correlated with annual/seasonal runoff or SRI

  # 2. Rip rap/revetment----------
  # total extent (miles) of rip rap installed on relevant stretches of the
  # Sacramento River directly related to total nest site availability (carrying
  # capacity), and thus limits growth, but annual data are lacking
  
  
  
  # 3. Foraging habitat (breeding season)---------
  # loss of foraging habitat for BANS (and more broadly, aerial insectivores) is
  # likely an important factor in pop declines; for BANS natural cover near
  # waterways/colonies is important; drought may also play an important role,
  # particularly for during winter just before migration with effects on
  # survival
  
  # > proposed metrics: indices of drought/primary productivity near the river
  
  # > Vegetation indices-----
  # from MODIS (see GIS pre-processing script)

  # NDVI/EVI within 1km of burrows, 2000-2020
  # --> mean value in Apr-Aug (breeding season) 
  cadat = read_csv('GIS/CA_1k/MOD13Q1-006-Statistics.csv', col_types = cols()) %>% 
    mutate(metric = gsub('.*_', '', Dataset),
           year = as.numeric(format(Date, '%Y')),
           month = as.numeric(format(Date, '%m')),
           WY = if_else(month >= 10, year + 1, year),
           # same seasons as for flow data above
           season = case_when(month >= 10 | month <= 3 ~ 'rainy',
                              month >= 4 & month <= 5 ~ 'breeding_early',
                              month >= 6 & month <= 8 ~ 'breeding_late',
                              TRUE ~ NA_character_),
           season = factor(season, 
                           levels = c('rainy', 'breeding_early', 
                                      'breeding_late'))) %>% 
    select(metric, Date, WY, month, season, 
           Mean, SD = `Standard Deviation`),
  
  # # increasing trend in EVI over time, with higher peak values in later years
  # ggplot(cadat %>% mutate(Date = as.numeric(format(Date, '%j'))) %>%
  #          filter(metric == 'EVI'),
  #        aes(Date, Mean, group = year, color = year)) +
  #   annotate("rect", xmin = 91, xmax = 243, ymin = 0.2, ymax = 0.5,
  #            fill = 'gray80', alpha = 0.5) +
  #   annotate("rect", xmin = 121, xmax = 181, ymin = 0.2, ymax = 0.5,
  #            fill = 'gray50', alpha = 0.5) +
  #   geom_point() + geom_line() + scale_color_viridis_c()
  # # May-June only:
  # ggplot(cadat %>% mutate(Date = as.numeric(format(Date, '%j'))) %>%
  #          filter(metric == 'EVI' & month >= 5 & month <= 6),
  #        aes(Date, Mean, group = year, color = year)) +
  #   geom_point() + geom_line() + scale_color_viridis_c()

  
  # summarize values by season/year
  cadat_sum = bind_rows(
    # by season
    cadat %>% 
      group_by(metric, WY, season) %>% 
      summarize(Mean = mean(Mean),
                SD = sqrt(sum(SD^2)), # error propagation
                .groups = 'drop') %>% 
      filter(!is.na(season)),
    # add combined breeding season
    cadat %>% 
      mutate(season = gsub('_early|_late', '', season)) %>% 
      filter(season == 'breeding') %>% 
      group_by(metric, WY, season) %>% 
      summarize(Mean = mean(Mean),
                SD = sqrt(sum(SD^2)), # error propagation
                .groups = 'drop')
    ),

  
  # # mean values increase over time, but not super strong relative to
  # # variance; breeding season variables very similar
  # ggplot(cadat_sum, aes(WY, Mean, color = season)) +
  #   geom_line() + geom_point() + geom_smooth() +
  #   geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), width = 0.1) +
  #   facet_wrap(~metric, nrow = 2)
  # 
  # # check correlations: strong correlation across seasonal EVI & NDVI
  # corrplot::corrplot.mixed(
  #   cor(cadat_sum %>% filter(season != 'rainy') %>% 
  #         select(metric, season, year, Mean) %>% 
  #         pivot_wider(names_from = c('metric', 'season'), values_from = Mean) %>%
  #         select(-year),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # 
  # # breeding early and late fairly correlated (0.59)
  # corrplot::corrplot.mixed(
  #   cor(cadat_sum %>% filter(metric == 'EVI') %>% 
  #         select(year, season, Mean) %>% 
  #         pivot_wider(names_from = season, values_from = Mean) %>%
  #         select(-year),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')

  # > Drought indices---------
  drought = purrr::map_df(c('pdsidv', 'zndxdv'),
                          ~get_drought_indices(datname = .x, statecode = '04', 
                                               climdivcode = '02')) %>% 
    mutate(month = format(date, '%m') %>% as.numeric(),
           year = as.numeric(year),
           WY = if_else(month >= 10, year + 1, year),
           season = case_when(month >= 10 | month <= 3 ~ 'rainy',
                              month >= 4 & month <= 5 ~ 'breeding_early',
                              month >= 6 & month <= 8 ~ 'breeding_late',
                              TRUE ~ NA_character_)) %>% 
    filter(WY >= 1985),
  
  # generate annual and seasonal summaries:
  drought_sum = bind_rows(
    drought %>%
      group_by(index, WY, season) %>%
      summarize(value = mean(value),
                .groups = 'drop') %>% 
      filter(!is.na(season)),
    drought %>% 
      mutate(season = gsub('_early|_late', '', season)) %>% 
      filter(season == 'breeding') %>% 
      group_by(index, WY, season) %>%
      summarize(value = mean(value),
                .groups = 'drop')) %>% 
    drop_na(),
  
  # # no clear trend; PDSI more variable?
  # ggplot(drought_sum %>% filter(WY > 1999 & season == 'breeding'), 
  #        aes(WY, value, color = season)) +
  #   geom_line() + geom_point() + geom_smooth() +
  #   facet_wrap(~index, nrow = 2)
  # 
  # # strong correlation
  # drought_sum %>% filter(WY>1999 & season == 'breeding') %>% 
  #   pivot_wider(names_from = index) %>% 
  #   ggplot(aes(pdsidv, zndxdv)) + geom_point() + geom_smooth()
  # 
  # # correlations:
  # corrplot::corrplot.mixed(
  #   cor(drought_sum %>% filter(WY >= 1999) %>%
  #         pivot_wider(names_from = c('index', 'season'),
  #                     values_from = value),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # # none with WY
  # # strong across indices within season, and across breeding (early/late)
  # # e.g. pdsi_breeding & zndx_breeding (0.89)
  
  
  # 4. Winter conditions------
  
  # > Vegetation indices from MODIS (see GIS pre-processing script)
  # NDVI/EVI within possible winter range in MX, 2000-2018
  # --> this is highly speculative!
  
  # take the mean value in Feb-Apr (ending just before migration and the arrival
  # in CA of the bulk of breeding birds in late Apr/early May )
  mxdat = read_csv('GIS/MX/MOD13A3-006-Statistics.csv', col_types = cols()) %>% 
    mutate(metric = gsub('.*_', '', Dataset),
           year = as.numeric(format(Date, '%Y')),
           month = as.numeric(format(Date, '%m')),
           WY = if_else(month >= 10, year + 1, year)) %>% 
    select(metric, WY, month, Mean, SD = `Standard Deviation`) %>% 
    filter(month >= 2 & month <= 4) %>% 
    group_by(metric, WY) %>% 
    summarize(Mean = mean(Mean),
              SD = sqrt(sum(SD^2)),
              .groups = 'drop'),
  
  # # slight increasing trend over time; 2010 stands out as particularly high
  # ggplot(mxdat, aes(year, Mean)) + geom_point() + geom_line() + 
  #   geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.1) +
  #   geom_smooth() +
  #   facet_wrap(~metric, nrow = 2)

  # # check correlations: very strong positive correlation (0.99) between NDVI & EVI
  # corrplot::corrplot.mixed(
  #   cor(mxdat %>% select(-year), method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')

  
  # Compiled predictors----------------
  # check for correlations among predictors

  modeldat = list(birddat %>% select(year, burrows, agr, pgr) %>% 
                    mutate(burrows1 = lead(burrows)),
                  # line up with current water year's flow data (starting the
                  # previous October)
                  mflow_sum %>% 
                    select(year = WY, flow14_cum3, flow14_cum2, flow14_rainy, 
                           flow14_breeding_early, flow14_breeding_late,
                           flow14_annual) %>% 
                    # add leading effects (from just before the next breeding season)
                    mutate(flow14_breeding = flow14_breeding_early + flow14_breeding_late,
                           flow14_cum31 = lead(flow14_cum3),
                           flow14_cum21 = lead(flow14_cum2),
                           flow14_rainy1 = lead(flow14_rainy),
                           flow14_breeding1 = lead(flow14_breeding),
                           flow14_breeding_early1 = lead(flow14_breeding_early),
                           flow14_breeding_late1 = lead(flow14_breeding_late),
                           flow14_annual1 = lead(flow14_annual)),
                  # line up with current breeding season EVI
                  cadat_sum %>% filter(metric == 'EVI') %>% 
                    select(-SD) %>% 
                    rename(year = WY) %>% 
                    pivot_wider(names_from = c('metric', 'season'),
                                values_from = Mean) %>% 
                    # consider effects early in the next breeding season
                    mutate(EVI_rainy1 = lead(EVI_rainy),
                           EVI_breeding_early1 = lead(EVI_breeding_early),
                           EVI_breeding_late1 = lead(EVI_breeding_late),
                           EVI_breeding1 = lead(EVI_breeding)),
                  drought_sum %>% filter(season == 'breeding') %>% 
                    pivot_wider(names_from = 'index') %>% 
                    select(year = WY, pdsi = pdsidv, zndx = zndxdv) %>% 
                    mutate(pdsi1 = lead(pdsi),
                           zndx1 = lead(zndx)),
                  mxdat %>% filter(metric == 'EVI') %>% 
                    select(year = WY, EVI_MX = Mean) %>% 
                    mutate(EVI_MX1 = lead(EVI_MX))
                  ) %>% 
    purrr::reduce(full_join) %>% 
    filter(year >= 1999),
  
  # # exploratory plots:
  # modeldat %>% select(pgr, year, burrows, burrows1:EVI_MX1) %>%
  #   mutate_at(vars(starts_with('flow')), ~log(. + 1)) %>%
  #   pivot_longer(-pgr) %>%
  #   ggplot(aes(value, pgr)) + geom_point() +
  #   geom_smooth(method = 'gam') +
  #   facet_wrap(~name, scales = 'free') +
  #   xlab(NULL)
  # #
  # # positive relationships with pgr:
  # # - burrows1 (of course)
  # # - EVI_breeding & EVI_breeding_late (not EVI_breeding_early); but few high values doing a lot of the work
  # # - zndx_breeding (&early/late); pdsidv_breeding (& early/late)
  # # - zndxdv_rainy (correlated with flow?)
  # # - flow14_rainy, flow14_annual, flow14_cum21
  # # - flow14_cum3 - slight
  # # - flow14_breeding_late1 - slight - could flows reduce the burrow counts in t-1 so growth appears bigger?
  # # - flow14_breeding (early & late) - nonlinear, increases at higher values
  # #
  # # negative relationship with pgr:
  # # - burrows (density dependence)
  # # - flow14_breeding_late1 - very slight -- reducing same-year burrow counts?
  # # - flow14_rainy1 & flow14_annual1 - nonlinear; flattens/starts increasing again at high values;
  # #    -->detrimental impacts of fresh cuts just before breeding season?
  # 
  # # nonlinear relationship with pgr:
  # # - zndxdv_breeding_early1 - intriguing
  # 
  # # no apparent effects on pgr:
  # # - EVI_breeding_early; EVI_breeding_early1; EVI_breeding_late1; EVI_breeding1
  # # - EVI_MX, EVI_MX1
  # # - EVI_rainy; EVI_rainy1
  # # - zndxdv_rainy1; zndxdv_breeding_late1; pdsidv_breeding1; pdsidv_rainy1
  # # - year
  # # - flow14_cum3; flow14_cum31
  # # - flow14_breeding1, flow14_breeding_early1
  
  
  # # check correlations with birddat & across metrics (log-transform flow metrics)
  # corrplot::corrplot.mixed(
  #   cor(modeldat %>% dplyr::select(-year) %>%
  #         mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
  #         mutate(agr = log(agr)) %>% drop_na(),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # # positive correlations:
  # # - flow14_annual/rainy & pgr/agr (0.68/0.68)
  # # - EVI_breeding_late & pgr (0.62)
  # # - flow14_rainy & flow14_annual (0.98)
  # # - ca_evi/ca_evi1 (0.86) (due to trend; correlation weaker when detrended)
  # # - flow14_breeding & flow14_breeding_early (0.88); with flow14_breeding_late (0.67)
  # # - flow14_cum3 & burrows1 (0.84) [but not agr/pgr] - reflects cumulative effect of flow14_rainy1
  # 
  # # negative correlations:
  # # - burrows & EVI_MX (-0.7); burrows1 & EVI_MX [but not agr/pgr]
  # # - flow14_cum3 & EVI_MX1 (-0.62) --> does this explain the apparent burrows/EVI_MX negative correlation?
  # #    --> negative relationship between wet winters in CA and MX? SOI?

  # corrplot::corrplot.mixed(
  #   cor(modeldat %>% dplyr::select(-year) %>%
  #         mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
  #         mutate(agr = log(agr)) %>% drop_na(),
  #       method = 'spearman'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # # positive correlations:
  # # - flow14_rainy & agr/pgr (0.65)
  # # - flow14_cum3 & burrows1 (0.83)
  # # - ca_evi & ca_evi1 (0.83) - autocorrelation (trend)
  # # negative correlations:
  # # - strongest is agr/pgr & burrows (-0.57) - some density dependence? or
  # # potential for growth is limited with #burrows already high

  
  # subset of predictors:
  modeldat %>%
    select(pgr, burrows, year,
           flow14_annual, flow14_annual1, flow14_cum3, flow14_breeding_early1,
           EVI_breeding, zndx, zndx1, pdsi, pdsi1,
           EVI_MX, EVI_MX1) %>%
    mutate_at(vars(starts_with('flow')), ~log(. + 1)) %>%
    pivot_longer(-pgr) %>%
    ggplot(aes(value, pgr)) + geom_point() +
    geom_smooth(method = 'gam') +
    facet_wrap(~name, scales = 'free_x') +
    xlab(NULL)
  # pgr (t to t+1) positively related to EVI_breeding, flow14_annual, pdsi, zndx (slight)
  # and negatively related to burrows, flow14_annual1 (nonlinear)

  # modeldat %>%
  #   select(pgr, burrows1, year,
  #          flow14_annual, flow14_annual1, flow14_cum3, flow14_breeding_early1,
  #          EVI_breeding, zndx, zndx1, pdsi, pdsi1,
  #          EVI_MX, EVI_MX1) %>%
  #   mutate_at(vars(starts_with('flow')), ~log(. + 1)) %>%
  #   pivot_longer(-burrows1) %>%
  #   ggplot(aes(value, burrows1)) + geom_point() +
  #   geom_smooth(method = 'gam') +
  #   facet_wrap(~name, scales = 'free_x') +
  #   xlab(NULL)
  # burrows(t+1) positively related to flow14_cum3, pgr, flow14_annual
  # and negatively related to EVI_MX/EVI_MX1, year (but nonlinear)
  # 
  # corrplot::corrplot.mixed(
  #   cor(modeldat %>%
  #         select(pgr, burrows, burrows1, year,
  #                flow14_annual, flow14_annual1, flow14_cum3, flow14_breeding_early1,
  #                EVI_breeding, zndx, zndx1, pdsi, pdsi1,
  #                EVI_MX, EVI_MX1) %>%
  #         mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
  #         drop_na(),
  #       method = 'spearman'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')
  # positive correlations:
  # - EVI_breeding & year (0.87)
  # - EVI_MX & year (0.59)
  # - flow14_annual & pgr (0.68)
  # - zndx & pdsi (0.91)
  #
  # # negative correlations:
  # - pgr & burrows (-0.57)
  # - EVI_MX & burrows (-0.59) (because of increasing trend?)
  # - year & burrows (-0.54)

  
  # ANALYSIS-----------
  # 1. growth rate----------
  # may be easier than modeling N because it can possibly avoid autocorrelation
  # issues 
  # --> drop 1999 only if including EVI data (missing before 2000
  # --> keep 2020 even though missing growth data (can predict growth from covariates & number of burrows in 2021)

  inputdat = list(
    # per-capita growth rate t to t+1
    y = modeldat %>% pull(pgr),
    # density-independent predictors of growth rate from t to t+1:
    Bpredictors = modeldat %>% 
      select(year, flow14_annual, flow14_cum3, pdsi) %>%
      mutate(t = year) %>% 
      # log transform flow data
      mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
      # calculate anomaly:
      pivot_longer(-year) %>%
      group_by(name) %>%
      mutate(max = max(value, na.rm = TRUE),
             min = min(value, na.rm = TRUE),
             mean = mean(value, na.rm = TRUE),
             zvalue = (value - mean) / (max - min)) %>%
      select(year, name, zvalue) %>%
      pivot_wider(names_from = name, values_from = zvalue) %>%
      select(-year),
    N = modeldat %>% pull(burrows) %>% round(digits = 0)/10000
    ),

  scale_params = modeldat %>% 
    select(year, flow14_annual, flow14_cum3, pdsi) %>%
    mutate(t = year) %>% 
    # log transform flow data
    mutate_at(vars(starts_with('flow')), ~log(.+1)) %>%
    # calculate anomaly:
    pivot_longer(-year) %>%
    group_by(name) %>%
    summarize(max = max(value, na.rm = TRUE),
              min = min(value, na.rm = TRUE),
              center = mean(value, na.rm = TRUE),
              scale = max - min,
              .groups = 'drop') %>% 
    select(-max, -min),
  
  # corrplot::corrplot.mixed(
  #   cor(bind_cols(pgr = inputdat$y,
  #                 year = inputdat$year,
  #                 inputdat$Bpredictors) %>%
  #         drop_na(),
  #       method = 'pearson'),
  #   order = 'hclust', tl.col = 'black', tl.pos = 'lt')


  # > simple growth model----
  # # # no density dependence
  # mod_growth = fit_BANSpopgrowth_model(
  #   inputdat, n.adapt = 5000, n.burnin = 10000, n.sample = 5000, n.chains = 3,
  #   holdout = 0, type = 'growth_simple',
  #   vars = c('rate', 'beta', 'sigma', 'ysim', 'Nsim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit')),
  # MCMCvis::MCMCsummary(mod_growth,
  #                      c('rate', 'beta', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_growth,
  #                    params = c('rate', 'beta', 'sigma'),
  #                    pdf = FALSE)
  # # fits well; n.eff sufficient; trace plots well-mixed
  # plot_model_predictions(mod_growth,
  #                        year = c(2000:2019, 2000:2020),
  #                        obsdat = modeldat, center = 0, scale = 10000)
  # # predictions ok, especially toward end of time series
  # get_sumstats(mod_growth, 'beta',
  #              id = c(inputdat$Bpredictors %>% names(), 'year'))
  # # strong effects of flow14_annual (+) and year (+)
  
  # > simple logistic model-----
  # including density dependent growth relative to a constant carrying capacity
  # (K) --> for this one, because N is now a predictor, can't have missing
  # values
  # inputdat2 = list(y = inputdat$y[c(1:6, 8:20)],
  #                  year = inputdat$year[c(1:6, 8:20)],
  #                  Bpredictors = inputdat$Bpredictors[c(1:6, 8:20), ],
  #                  # pop size estimate from t (center on 10000, in thousands)
  #                  N = inputdat$N[c(1:6, 8:20)]),
  # 
  # mod_growth_logistic = fit_BANSpopgrowth_model(
  #   inputdat2, n.adapt = 5000, n.burnin = 15000, n.sample = 50000, n.chains = 3,
  #   holdout = 0, type = 'growth_logistic',
  #   vars = c('rate', 'beta', 'K', 'sigma', 'ysim', 'Nsim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 10),
  # MCMCvis::MCMCsummary(mod_growth_logistic,
  #                      c('rate', 'beta', 'K', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_growth_logistic,
  #                    params = c('rate', 'beta', 'K', 'sigma'),
  #                    pdf = FALSE)
  # # fits well; n.eff low for K & trace plots a weird for K
  # plot_model_predictions(mod_growth_logistic,
  #                        year = c(2000:2005, 2007:2019, 2000:2006, 2008:2020),
  #                        obsdat = modeldat,
  #                        scale = 10000, center = 0)
  # # predictions reasonable
  # get_sumstats(mod_growth_logistic, 'beta',
  #              id = c(inputdat2$Bpredictors %>% names(), 'year'))
  # # strong effects of flow14_cum3 (+) only

  # > simple logistic 2-------
  # alt parameterization to directly estimate strength of density dependence
  # without having to estimate K; effect of density dependence still assumed to
  # be constant
  # mod_growth_logistic2 = fit_BANSpopgrowth_model(
  #   inputdat2, n.adapt = 5000, n.burnin = 15000, n.sample = 150000, n.chains = 3,
  #   holdout = 0, type = 'growth_logistic_alt',
  #   vars = c('rate', 'beta', 'dd', 'sigma', 'ysim', 'Nsim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 50),
  # MCMCvis::MCMCsummary(mod_growth_logistic2,
  #                      c('rate', 'beta', 'dd', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_growth_logistic2,
  #                    params = c('rate', 'beta', 'dd', 'sigma'),
  #                    pdf = FALSE)
  # # fits well; n.eff sufficient; trace plots well-mixed
  # plot_model_predictions(mod_growth_logistic2,
  #                       year = c(2000:2005, 2007:2019, 2000:2006, 2008:2020),
  #                       obsdat = modeldat,
  #                       center = 0, scale = 10000)
  # # predictions reasonable
  # get_sumstats(mod_growth_logistic2, c('beta', 'dd'),
  #              id = c(inputdat2$Bpredictors %>% names(), 'year', 'dd'))
  # # strong effects of flow14_cum3 (+) and dd (-); pdsi  and flow14_annual close-ish?
  
  # > logistic with variable density dependence------
  # allows strength of density dependence to vary annually from a normal
  # distribution with mean mu.dd and error sigma.d
  # mod_growth_vardd = fit_BANSpopgrowth_model(
  #   inputdat2, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
  #   holdout = 0, type = 'growth_logistic_variabledd',
  #   vars = c('rate', 'beta', 'sigma', 'dd', 'mu.d', 'sigma.d',
  #            'ysim', 'Nsim', 'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 100),
  # MCMCvis::MCMCsummary(mod_growth_vardd,
  #                      c('rate', 'beta', 'sigma', 'mu.d', 'sigma.d',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_growth_vardd,
  #                    params = c('rate', 'beta', 'mu.d', 'sigma', 'sigma.d'),
  #                    pdf = FALSE)
  # # fits well; n.eff low for mu.d, rate; trace plots a bit odd as well for these
  # plot_model_predictions(mod_growth_vardd,
  #                        year = c(2000:2005, 2007:2019, 2000:2006, 2008:2020),
  #                        obsdat = modeldat, scale = 10000, center = 0)
  # # predictions reasonable
  # get_sumstats(mod_growth_vardd, c('beta','mu.d'),
  #              id = c(inputdat2$Bpredictors %>% names(), 'year', 'dd'))
  # # strong effects of flow14_cum3 (+) and dd (-)
  
  # > logistic with trend in K-------
  # mod_growth_trendK = fit_BANSpopgrowth_model(
  #   inputdat2, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
  #   holdout = 0, type = 'growth_logistic_trenddd',
  #   vars = c('rate', 'beta', 'dd', 'kappa', 'dd_intercept', 'sigma',
  #            'ysim', 'Nsim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 100),
  # MCMCvis::MCMCsummary(mod_growth_trendK,
  #                      c('rate', 'beta', 'kappa', 'dd_intercept', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_growth_trendK,
  #                    params = c('rate', 'beta', 'kappa', 'dd_intercept', 'sigma'),
  #                    pdf = FALSE)
  # # fits well; n.eff low dd_intercept, kappa, beta[4], rate (confounded?)
  # plot_model_predictions(mod_growth_trendK,
  #                       year = c(2000:2005, 2007:2019, 2000:2006, 2008:2020),
  #                       obsdat = modeldat, scale = 10000, center = 0)
  # # predictions reasonable
  # get_sumstats(mod_growth_trendK, c('beta','kappa','dd_intercept'),
  #              id = c(inputdat2$Bpredictors %>% names(), 'year', 'dd_slope', 'dd'))
  # # strong effects of flow14_cum3 (+) and dd (-) [not dd_slope]
  
  # > logistic with variable K-------
  # let carrying capacity change with cumulative flow; now that N is a
  # predictor, have to drop missing values
  inputdat3 = list(y = inputdat$y[c(1:7,9:22)],
                   Bpredictors = inputdat$Bpredictors[c(1:7,9:22), c('flow14_annual', 'pdsi', 't')],
                   Kpredictors = inputdat$Bpredictors[c(1:7,9:22), 'flow14_cum3'],
                   # pop size estimate from t (center on 15000, scale in tens of
                   # thousands)
                   N = inputdat$N[c(1:7,9:22)],
                   predmatrix = tibble(flow14_annual = seq(-0.6, 0.4, length.out = 22),
                                       pdsi = seq(-0.4, 0.6, length.out = 22),
                                       t = inputdat$Bpredictors$t,
                                       flow14_cum3 = seq(-0.6, 0.4, length.out = 22))),

  mod_growth_varK = fit_BANSpopgrowth_model(
    inputdat3, n.adapt = 10000, n.burnin = 40000, n.sample = 600000, 
    n.chains = 3, holdout = 0, type = 'growth_logistic_varK',
    vars = c('r0', 'beta', 'kappa', 'k0', 'k', 
             'var.o', 'var.p', 'var.k', 'ysim', 'Nsim', 'ypred', 
             'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
    thin = 300),
  MCMCvis::MCMCsummary(mod_growth_varK,
                       c('r0', 'beta', 'kappa', 'k0', 
                         'var.o', 'var.p', 'var.k', 
                         'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  MCMCvis::MCMCtrace(mod_growth_varK,
                     params = c('r0', 'beta', 'kappa', 'k0',
                                'var.o', 'var.p', 'var.k'),
                     pdf = FALSE)
  # model fits well; n.eff sufficient; trace plots well mixed
  plot_model_predictions(mod_growth_varK,
                         year = c(1999:2005, 2007:2020, 1999:2021),
                         obsdat = modeldat, scale = 10000, center = 0)
  plot_cov_relationships(mod_growth_varK,
                         r0 = 0,
                         nm = names(inputdat3$predmatrix),
                         xvals = inputdat3$predmatrix %>% 
                           mutate(index = c(1:nrow(.))) %>% 
                           pivot_longer(-index) %>% 
                           left_join(scale_params) %>% 
                           mutate(value = value * scale + center,
                                  value = case_when(name %in% c('flow14_annual', 'flow14_cum3') ~ exp(value),
                                                    TRUE ~ value)) %>% 
                           select(index, name, value) %>% 
                           arrange(name, index),
                         obsdat = modeldat %>%
                           select(pgr, flow14_annual, pdsi, flow14_cum3, t = year) %>% 
                           pivot_longer(-pgr, names_to = 'cov'))
  # predictions reasonable
  get_sumstats(mod_growth_varK, c('beta', 'kappa', 'k0'),
               id = c(inputdat3$Bpredictors %>% names(), 
                      inputdat3$Kpredictors %>% names(), 'k0'))
  # strong effects of flow14_cum3 on strength of density dependence (+); dd
  # negative; other effects as expected - trending positive for flow14_annual &
  # pdsi, no trend in year
  plot_K_predictions(mod_growth_varK,
                     predictor = seq(400000, 10000000, 500000),
                     obsdat = modeldat,
                     scale = 10000, center = 0) + ylim(0, NA)
  
  
  
  
  
  # >>hold out version----------
  mod_growth_varK_holdout = fit_BANSpopgrowth_model(
    inputdat3, n.adapt = 10000, n.burnin = 40000, n.sample = 600000, 
    n.chains = 3, holdout = 3, type = 'growth_logistic_varK',
    vars = c('rate', 'beta', 'kappa', 'dd_intercept', 'dd', 
             'var.o', 'var.p', 'var.d', 'ysim', 'Nsim', 
             'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
    thin = 200),
  MCMCvis::MCMCsummary(mod_growth_varK_holdout,
                       c('rate', 'beta', 'kappa', 'dd_intercept', 'var.o', 'var.p',
                         'var.d', 'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  MCMCvis::MCMCtrace(mod_growth_varK_holdout,
                     params = c('rate', 'beta', 'kappa', 'dd_intercept',
                                'var.o', 'var.p', 'var.d'),
                     pdf = FALSE)
  # model fits well; n.eff sufficient; trace plots well mixed
  plot_model_predictions(mod_growth_varK_holdout,
                         year = c(1999:2005, 2007:2020, 1999:2006, 2008:2021),
                         obsdat = modeldat, scale = 10000, center = 0)
  # predictions reasonable
  get_sumstats(mod_growth_varK_holdout, c('beta', 'kappa', 'dd_intercept'),
               id = c(inputdat3$Bpredictors %>% names(), 
                      inputdat3$Kpredictors %>% names(), 'dd'))
  # strong effects of flow14_cum3 on strength of density dependence (+); dd
  # negative; other effects as expected - trending positive for flow14_annual &
  # pdsi, no trend in year
  plot_K_predictions(mod_growth_varK_holdout,
                     predictor = seq(400000, 10000000, 500000),
                     obsdat = modeldat,
                     scale = 10000, center = 0) + ylim(0, NA)
  
  # > 2. N----------
  # these models are slower to run, and are parameterized differently to
  # directly take burrow counts in as the "y" observations, estimated with error
  # as a count of the true population N, and essentially calculate the growth
  # rate internally.
  
  # re-do input data: burrow counts need to be integers, so no scaling, but can
  # have missing values; include counts through 2020 to account for growth rate
  # estimated for 2019-20 interval
  inputdatN = list(y = modeldat %>% filter(year > 1999) %>%
                     pull(burrows) %>% round(digits = 0),
                   year = (modeldat %>% filter(year > 1999) %>%
                             pull(year) - 2010)/10,
                   Bpredictors = inputdat$Bpredictors),

  # > simple poisson-------
  # mod_N_poisson = fit_BANSpopgrowth_model(
  #   inputdatN, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
  #   holdout = 0, type = 'N_simple_poisson',
  #   vars = c('rate', 'beta', 'error', 'Nsim', 'ysim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 100),
  # MCMCvis::MCMCsummary(mod_N_poisson,
  #                      c('rate', 'beta',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_N_poisson,
  #                    params = c('rate', 'beta'),
  #                    pdf = FALSE)
  # # fit ok; n.eff very low; trace plots terrible; predictions suspiciously precise
  # plot_model_predictions(mod_N_poisson,
  #                        year = c(2000:2019, 2000:2020),
  #                        obsdat = modeldat, scale = 1, center = 0)
  # # predictions reasonable
  # get_sumstats(mod_N_poisson, c('beta'),
  #              id = c(inputdatN$Bpredictors %>% names(), 'year'))
  # # strong effects of flow14_annual (+) and year (+)
  
  # # > logistic poisson---------
  # mod_N_logistic_poisson = fit_BANSpopgrowth_model(
  #   inputdatN, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
  #   holdout = 0, type = 'N_logistic_poisson',
  #   vars = c('rate', 'beta', 'dd', 'sigma', 'Nsim', 'ysim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 100)
  # MCMCvis::MCMCsummary(mod_N_logistic_poisson,
  #                      c('rate', 'beta', 'dd', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_N_logistic_poisson,
  #                    params = c('rate', 'beta', 'dd', 'sigma'),
  #                    pdf = FALSE)
  # # fit ok; n.eff very low; dd tiny; trace plots terrible 
  # plot_model_predictions(mod_N_logistic_poisson,
  #                        year = c(2000:2019, 2000:2020),
  #                        obsdat = modeldat, scale = 1, center = 0)
  # # predictions suspiciously precise
  # get_sumstats(mod_N_logistic_poisson, c('beta', 'dd'),
  #              id = c(inputdatN$Bpredictors %>% names(), 'year', 'dd'))
  # # strong effects of flow14_annual (+) only

  # > simple negative binomial----------
  # allows more error in the count process than poisson

  # mod_N_negbin = fit_BANSpopgrowth_model(
  #   inputdatN, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
  #   holdout = 0, type = 'N_simple_negbin',
  #   vars = c('rate', 'beta', 'sigma', 'Nsim', 'ysim',
  #            'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
  #   thin = 100),
  # MCMCvis::MCMCsummary(mod_N_negbin,
  #                      c('rate', 'beta', 'sigma',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_N_negbin,
  #                    params = c('rate', 'beta', 'sigma'),
  #                    pdf = FALSE)
  # # fit ok; n.eff a bit low for rate, with ok traceplots
  # plot_model_predictions(mod_N_negbin,
  #                        year = c(2000:2019, 2000:2020),
  #                        obsdat = modeldat, scale = 1, center = 0)
  # # predictions reasonable
  # get_sumstats(mod_N_negbin, c('beta'),
  #              id = c(inputdatN$Bpredictors %>% names(), 'year'))
  # # strong effects of flow14_annual (+) and year (+)

  # > logistic negative binomial----------
  # allows more error in the count process than poisson; add density dependence

  mod_N_logistic_negbin = fit_BANSpopgrowth_model(
    inputdatN, n.adapt = 5000, n.burnin = 15000, n.sample = 300000, n.chains = 3,
    holdout = 0, type = 'N_logistic_negbin',
    vars = c('rate', 'beta', 'dd', 'sigma', 'k', 'Nsim', 'ysim',
             'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
    thin = 100),
  # MCMCvis::MCMCsummary(mod_N_logistic_negbin,
  #                      c('rate', 'beta', 'dd', 'sigma', 'k',
  #                        'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  # MCMCvis::MCMCtrace(mod_N_logistic_negbin,
  #                    params = c('rate', 'beta', 'dd', 'sigma', 'k'),
  #                    pdf = FALSE)
  # # fit very good; n.eff super low for rate, dd, beta[2] and beta[4]; trace plots for
  # # these also terrible -- looks like a longer burn-in needed?
  # plot_model_predictions(mod_N_logistic_negbin,
  #                        year = c(2000:2019, 2000:2020),
  #                        obsdat = modeldat, scale = 1, center = 0)
  # # predictions very good...
  # get_sumstats(mod_N_logistic_negbin, c('beta', 'dd'),
  #              id = c(inputdatN$Bpredictors %>% names(), 'year', 'dd'))
  # # strong effects of flow14_cum3 (+) and dd (-)

  # > logistic negbin variable K----------
  inputdatN2 = list(y = inputdatN$y,
                    year = inputdatN$year,
                    Bpredictors = inputdatN$Bpredictors[, c('flow14_annual', 'pdsi')],
                    Kpredictors = inputdatN$Bpredictors[, 'flow14_cum3']),

  mod_N_logistic_negbin_varK = fit_BANSpopgrowth_model(
    inputdatN2, n.adapt = 5000, n.burnin = 20000, n.sample = 100000, n.chains = 3,
    holdout = 0, type = 'N_logistic_negbin_varK',
    vars = c('rate', 'beta', 'kappa', 'dd_intercept', 'sigma', 'sigma.d', 'dd',
             'k', 'Nsim', 'ysim',
             'pvalue.mean', 'pvalue.sd', 'pvalue.fit'),
    thin = 50),
  MCMCvis::MCMCsummary(mod_N_logistic_negbin_varK,
                       c('rate', 'beta', 'kappa', 'dd_intercept',
                         'sigma', 'sigma.d', 'k',
                         'pvalue.mean', 'pvalue.sd', 'pvalue.fit'))
  MCMCvis::MCMCtrace(mod_N_logistic_negbin_varK,
                     params = c('rate', 'beta', 'kappa', 'dd_intercept',
                                'sigma', 'sigma.d', 'k'),
                     pdf = FALSE)
  # fit ok; n.eff=0 for sigma.d? also low for rate, dd_intercept, kappa;
  # trace plots a bit weird too
  # --> note, the estimates for sigma.d and sigma are super tiny- need a different prior?
  # also really small for kappa
  plot_model_predictions(mod_N_logistic_negbin_varK,
                         year = c(2000:2019, 2000:2020),
                         obsdat = modeldat, scale = 1, center = 0)
  # predictions reasonable
  get_sumstats(mod_N_logistic_negbin_varK, c('beta', 'kappa', 'dd_intercept'),
               id = c(inputdatN2$Bpredictors %>% names(), 'year',
                      inputdatN2$Kpredictors %>% names(), 'dd'))
  # strong effects of flow14_cum3 on density dependence (+) and of dd (-);
  # expected suggestions of effects for flow14_annual (+ 82%), pdsi (+ 77%), and
  # year (- 76%)
  
  plot_K_predictions(mod_N_logistic_negbin_varK,
                     predictor = seq(400000, 10000000, 500000),
                     obsdat = modeldat,
                     scale = 1, center = 0)
)
