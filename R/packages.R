# Load all your packages before calling make().

library(drake)
library(tidyverse)
library(MODIStsp) #for accessing MODIS data (e.g. NDVI)
# MODIStsp::MODIStsp_get_prodnames()
# MODIStsp::MODIStsp_get_prodlayers("M*D13Q1")