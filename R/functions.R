# Custom functions are an important part of a drake workflow.
# This is where you write them.
# Details: https://books.ropensci.org/drake/plans.html#functions

render_Rmd = function(pathin, pathout) {
  rmarkdown::render(pathin, 
                    output_file = here::here(pathout))
}

# FUNCTION FAILS: alert handshake failure
# get_cdec_waterdat = function(station = 'VIN',
#                              sensor = 41, #mean daily flow
#                             start.date = '2010-10-01',
#                             end.date = as.character(Sys.Date())) {
#   # if (length(stations) > 1) {
#   #   stations = paste(st, collapse = '%2C+')
#   # }
#   # URL = 'https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet'
#   URL = 'https://cdec.water.ca.gov/dynamicapp/req/JasonDataServlet'
#   query = paste0('?', station, '&SensorNums=', sensor, '&dur_code=D',
#                  '&Start=', start.date, '&End=', end.date)
#   RCurl::getURL(paste0(URL, query))
# }

compile_waterdat = function(dir, type, ...) {
  filelist = list.files(dir, pattern = type, full.names = TRUE)
  purrr::map_df(filelist, read_csv, ...)
}

get_drought_indices = function(datname, statecode = NULL, climdivcode = NULL) {
  # 04 = california
  
  url <- "ftp://ftp.ncdc.noaa.gov/pub/data/cirs/climdiv/"
  filelist = RCurl::getURL(url, dirlistonly = TRUE) %>% 
    strsplit("\\r\\n")
  df = tibble(name = filelist[[1]]) %>% 
    filter(grepl(datname, name))
  dat <- try(RCurl::getURL(paste0(url, df$name[1]))) %>% 
    strsplit("   \\n")
  res <- read_table(dat[[1]],
             col_names = c('ID', 'JAN', 'FEB', 'MAR', 'APR', 
                           'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 
                           'OCT', 'NOV', 'DEC')) %>% 
    mutate(index = datname,
           state = substr(ID, 1, 2),
           climdiv = paste0(substr(ID, 3, 4)),
           year = substr(ID, 7, 10)) %>%
    pivot_longer(JAN:DEC, names_to = 'month', values_to = 'value') %>%
    mutate(climname = recode(climdiv, 
                            `01` = "North Coast Drainage",
                            `02` = "Sacramento Drainage",
                            `03` = "Northeast Interior Basins",
                            `04` = "Central Coast",
                            `05` = "San Joaquin Drainage",
                            `06` = "South Coast Drainage",
                            `07` = "Southeast Desert Basin"),
           date = as.Date(paste0(year, '-', month, '-01'), format = '%Y-%b-%d'),
           value = case_when(value < -99.9 ~ NA_real_, 
                             TRUE ~ value)) %>%
    arrange(index, state, climdiv, date)
  
  if (!is.null(statecode)) {
    res <- res %>% filter(state == statecode)
  }
  if (!is.null(climdivcode)) {
    res <- res %>% filter(climdiv == climdivcode)
  }
  return(res)
  ## Climate divisions within California:
  ## 01: North Coast Drainage
  ## 02: Sacramento Drainage
  ## 03: Northeast Interior Basins
  ## 04: Central Coast **TomKat**
  ## 05: San Joaquin Drainage
  ## 06: South Coast Drainage
  ## 07: Southeast Desert Basin
}

get_sumstats = function(mod, param, id, prob = 0.95) {
  samples = MCMCvis::MCMCchains(mod, params = param, ISB = TRUE)
  
  bind_cols(
    id,
    MCMCvis::MCMCpstr(mod, params = param, func = median) %>% unlist(),
    do.call(rbind, 
            MCMCvis::MCMCpstr(mod, params = param, 
                              func = function(y) HDInterval::hdi(y, prob))) %>% 
      as_tibble(),
    apply(samples, 2, function(x) 1-ecdf(x)(0)) %>% set_names(NULL),
    apply(samples, 2, function(x) ecdf(x)(0)) %>% set_names(NULL)
  ) %>% 
    set_names('id', 'median', 'lci', 'uci', 'probg0', 'probl0')
}
