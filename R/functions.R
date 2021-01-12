# Custom functions are an important part of a drake workflow.
# This is where you write them.
# Details: https://books.ropensci.org/drake/plans.html#functions

generate_data <- function() {
  tibble(x = rnorm(1e5), y = rnorm(1e5))
}

fit_model <- function(data) {
  summary(lm(y ~ x, data = data))
}

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

compile_waterdat = function(dir, type) {
  filelist = list.files(dir, pattern = type, full.names = TRUE)
  purrr::map_df(filelist, read_csv, col_types = cols())
}