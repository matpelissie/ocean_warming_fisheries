###-###-###-###-###-###-###-###-###
#
# Load packages
#
# 00_Rpackages.R
#
###-###-###-###-###-###-###-###-###

library(tidyverse) # version 2.0.0
library(patchwork) # version 1.2.0
library(cowplot) # version 1.1.3
library(MuMIn) # version 1.47.5
library(chngpt) # version 2023.11-29
library(asdetect) # version 0.1.0 (pak::pak("caboulton/asdetect"))
library(ggtext) # version 0.1.2
library(ggrepel) # version 0.9.5
library(stringi) # version 1.8.3
library(plyr) # version 1.8.9
library(pracma) # version 2.4.4
library(lubridate) # version 1.9.3
library(scales) # verion 1.3.0
library(ncdf4) # version 1.22
library(sf) # version 1.0-15
library(lwgeom) # version 0.2-15
library(hier.part) # version 1.0-6
# to install the right version:
# devtools::install_version("hier.part", version = "1.0-6", repos = "http://cran.us.r-project.org")
library(rfishbase) # version 4.1.2 (previous versions miss functions and values)
# to install the right version:
# remotes::install_version("rfishbase", version="4.1.2", repos="https://cran.r-project.org", force=TRUE)
library(FishLife) # version 3.1.0 (remotes::install_github("James-Thorson-NOAA/FishLife"))
# Ignore the call to install earlier rfishbase version while loading FishLife
