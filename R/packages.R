# Load packages 

library(tidyverse)
library(rjags)
library(MCMCvis)
library(HDInterval)
library(patchwork)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("lag", "dplyr")
