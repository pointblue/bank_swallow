get_drought_indices = function(datname) {
  # 04 = california
  
  url <- "ftp://ftp.ncdc.noaa.gov/pub/data/cirs/climdiv/"
  filelist = RCurl::getURL(url, dirlistonly = TRUE) %>% 
    strsplit("\\r\\n")
  df = tibble(name = filelist[[1]]) %>% 
    filter(grepl(datname, name))
  dat <- try(RCurl::getURL(paste0(url, df$name[1]))) %>% 
    strsplit("   \\n")
  read_table(dat[[1]],
             col_names = c('ID', 'JAN', 'FEB', 'MAR', 'APR', 
                           'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 
                           'OCT', 'NOV', 'DEC')) %>% 
    mutate(index = datname)
}

cleanup_droughtdat = function(df, statecode = NULL, climdivcode = NULL) {
  res <- df %>% 
    mutate(state = substr(ID, 1, 2),
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
           month = format(date, '%m') %>% as.numeric(),
           year = as.numeric(year),
           WY = if_else(month >= 10, year + 1, year),
           season = case_when(month >= 10 | month <= 3 ~ 'rainy',
                              month >= 4 & month <= 8 ~ 'breeding',
                              TRUE ~ NA_character_),
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
}
