###-###-###-###-###-###-###-###-###
#
# Download and prepare datasets
#
# 01_datasets.R
#
###-###-###-###-###-###-###-###-###

library(tidyverse)
dir.create("data/raw_data/", showWarnings = FALSE)

# Fisheries data -------

## Download RAM Legacy Stock Assessment Database (v4.61) ----------------------

download.file(
  "https://zenodo.org/records/7814638/files/RAMLDB%20v4.61.zip?download=1",
  destfile = "data/raw_data/RAMLDB v4.61.zip")
unzip("data/raw_data/RAMLDB v4.61.zip",
      exdir = "data/raw_data/RAMLDB v4.61")
unlink("data/raw_data/RAMLDB v4.61.zip")

dir.create("data/RAMLDB v4.61")
file.copy("data/raw_data/RAMLDB v4.61/R Data/", "data/RAMLDB v4.61", recursive=TRUE)


# Download FAO fish catch data --------------------------------------------

dir.create("data/FAO_data", showWarnings = FALSE)

# Catch data
download.file(
  "https://www.fao.org/fishery/static/Data/Capture_2024.1.0.zip",
  destfile = "data/raw_data/Capture_2024.1.0.zip")
unzip("data/raw_data/Capture_2024.1.0.zip",
      exdir = "data/FAO_data/Capture_2024.1.0")
unlink("data/raw_data/Capture_2024.1.0.zip")

# Taxonomy
download.file(
  "https://www.fao.org/fishery/static/ASFIS/ASFIS_sp.zip",
  destfile = "data/raw_data/ASFIS_sp.zip")
unzip("data/raw_data/ASFIS_sp.zip",
      exdir = "data/FAO_data/ASFIS_sp")
unlink("data/raw_data/ASFIS_sp.zip")



# Spatial data --------------------------------------------------

## Area polygons ------------------------------------------------------------

### FAO fishing areas -------------------------------------------------------

# The following download chunk is failing at the moment
# so the shape file is directly provided in the repository
# download.file("https://www.fao.org/fishery/geoserver/fifao/ows?service=WFS&request=GetFeature&version=1.0.0&typeName=FAO_AREAS_ERASE_LOWRES&outputFormat=SHAPE-ZIP",
#               destfile = "data/raw_data/FAO_AREAS_ERASE_LOWRES.zip")
# unzip("data/raw_data/FAO_AREAS_ERASE_LOWRES.zip",
#       exdir = "data/spatial_data/FAO_AREAS_ERASE_LOWRES")
# unlink("data/raw_data/FAO_AREAS_ERASE_LOWRES.zip")


### Download LME polygons -------------------------------

download.file("http://geonode.iwlearn.org/geoserver/wfs?format_options=charset%3AUTF-8&typename=geonode%3Almes66gcd&outputFormat=SHAPE-ZIP&version=1.0.0&service=WFS&request=GetFeature&access_token=689743d98d5611ef9567000d3ab6a624",
              destfile = "data/raw_data/LME66.zip")
unzip("data/raw_data/LME66.zip",
      exdir = "data/spatial_data/LME66")
unlink("data/raw_data/LME66.zip")


### Global oceans polygons already downloaded (after manual request) -------

# from https://www.marineregions.org/downloads.php#goas
# (Global Oceans and Seas v01 (2021-12-14, 88 MB) [Shapefile])
# and simplified using the following code to speed up computation and allow storage online

# oceans <- sf::read_sf(dsn = "data/raw_data/GOaS_v1_20211214/",
#                       layer = "goas_v01")
# oceans_simpl <- rmapshaper::ms_simplify(oceans, keep = 0.001,
#                                         keep_shapes = FALSE)
# # Save simplified ocean polygons
# dir.create("data/spatial_data/ocean_simpl/")
# sf::st_write(oceans_simpl, "data/spatial_data/ocean_simpl/ocean_simpl.shp")
# unlink("data/raw_data/GOaS_v1_20211214/")


### Assign stocks to LME ----------------------------------------------------

sf::sf_use_s2(FALSE) # To remove error for intersection

lme <- sf::read_sf(dsn = "data/spatial_data/LME66/", layer = "lmes66gcd")
oceans <- sf::read_sf(dsn = "data/spatial_data/ocean_simpl/",
                      layer = "ocean_simpl") %>%
  dplyr::filter(grepl("Ocean", name))

# Group all LMEs
lme_merge <- lme %>% dplyr::summarise()

# Substract LMEs from oceans
high_seas <- sf::st_difference(oceans, lme_merge)
high_seas <- high_seas %>%
  dplyr::rename(LME_NAME = name) %>%
  dplyr::mutate(LME_NUMBER = 67:74)

# Combine LMEs and high sea regions
lme_all <- dplyr::bind_rows(lme, high_seas) %>%
  dplyr::mutate(area_lme = sf::st_area(.)) %>%
  sf::st_make_valid()

# sf::write_sf(lme_all, dsn = "data/spatial_data/LME_highseas.shp")



## Stock polygons ----------------------------------------------------------

### Download stock polygons ----------------------

# Download shapefiles
download.file(
  "https://chrismfree.com/files/ramldb_boundaries.zip",
  destfile = "data/raw_data/ramldb_boundaries.zip")
unzip("data/raw_data/ramldb_boundaries.zip",
      exdir = "data/raw_data/ramldb_boundaries")
unlink("data/raw_data/ramldb_boundaries.zip")

file.copy("data/raw_data/ramldb_boundaries/ramldb_boundaries/",
          "data/spatial_data/stock_boundaries", recursive=TRUE)

# Download metadata
download.file(
  "https://chrismfree.com/files/ramldb_v3.8_stock_boundary_table_v2_formatted.xlsx",
  destfile = "data/spatial_data/stock_boundaries/ramldb_v3.8_stock_boundary_table_v2_formatted.xlsx")

# Additional stock polygons delineated were directly added
# to the repository "data/spatial_data/stock_boundaries/additional_boundaries"


### List polygons available -------------------------------------------------

# Read metadata
main <- readxl::read_xlsx("data/spatial_data/stock_boundaries/ramldb_v3.8_stock_boundary_table_v2_formatted.xlsx") %>%
  dplyr::mutate(batch="batch1") # Free's polygons
length(unique(main$stockid))  # 685 stocks

add <- readr::read_csv("data/spatial_data/stock_boundaries/additional_ramldb_stock_boundary_formatted.csv") %>%
  dplyr::mutate(batch="batch2") # New polygons
length(unique(add$stockid))  # 64 stocks

# Build merged dataset
fulldata <- main %>%
  dplyr::bind_rows(add) %>%
  dplyr::select(assessid, stockid, shp_source, shp_path, zone_col, zones, batch) %>%
  dplyr::rename(assessid_shapefile=assessid) %>%
  dplyr::arrange(stockid)
length(unique(fulldata$assessid_shapefile)) # 749 stocks


### Build spatial polygons --------------------------------------------------

shpdir1 <- "data/spatial_data/stock_boundaries/ramldb_boundaries/"
shpdir2 <- "data/spatial_data/stock_boundaries/additional_boundaries/"

# Read RAMLDB stock boundary shapefiles
assess1 <- unique(fulldata[fulldata$batch=="batch1",]$assessid_shapefile)
stocks1 <- unique(fulldata[fulldata$batch=="batch1",]$stockid)

boundaries1 <- mapply(function(x, y)
  sf::st_read(dsn=paste0(shpdir1, x,".shp"), quiet=TRUE) %>%
    dplyr::mutate(assessid = x, stockid = y),
  assess1, stocks1, SIMPLIFY=FALSE)


assess2 <- unique(fulldata[fulldata$batch=="batch2",]$assessid_shapefile)
stocks2 <- unique(fulldata[fulldata$batch=="batch2",]$stockid)

boundaries2 <- mapply(function(x, y)
  sf::st_read(dsn=paste0(shpdir2, x,".shp"), quiet=TRUE) %>%
    dplyr::mutate(assessid = x, stockid = y),
  assess2, stocks2, SIMPLIFY=FALSE)

# Pool stock boundaries from different origins
boundaries <- c(boundaries1, boundaries2)
stocks <- c(stocks1, stocks2)

boundaries_crs <-
  lapply(boundaries, function(x)
    sf::st_transform(x, crs="+proj=longlat +datum=WGS84") %>%
      dplyr::select(assessid, stockid))

stock_boundaries <-
  do.call(dplyr::bind_rows, boundaries_crs) %>%
  terra::vect()

rm(boundaries, boundaries1, boundaries2)
gc()



## Estimate LME overlap ----------------------------------------------------

frac_lme <- lapply(1:length(stocks), function(i){

  # Load polygon:
  pol <- boundaries_crs[[i]]

  # To avoid errors with latitude reaching 90:
  if(sf::st_bbox(pol)[[4]]>89.99){
    pol <- pol %>%
      sf::st_crop(xmin = -180, ymin = -90,
                  xmax = 180, ymax = 89.99)
  }

  # Estimate polygon area:
  pol <- pol %>%
    sf::st_make_valid() %>%
    dplyr::mutate(area_pol=sf::st_area(.))

  # Convert polygon CRS if different:
  if(sf::st_crs(pol)$input != "WGS 84"){
    pol <- pol %>%
      sf::st_transform(crs = "WGS84")
  }

  # Divide polygon into LMEs, estimate areas and proportions:
  frac <- sf::st_intersection(pol, lme_all) %>%
    dplyr::mutate(area_lme=sf::st_area(.),
                  prop_area = round(as.numeric(area_lme/area_pol), 5),
                  stockid = fulldata %>%
                    dplyr::filter(stockid==stocks[i]) %>%
                    dplyr::pull(stockid)) %>%
    suppressWarnings() %>%
    suppressMessages()

  # Simplify output:
  frac_simpl <- frac %>%
    sf::st_drop_geometry(frac) %>%
    dplyr::select(stockid, LME_NUMBER, LME_NAME, prop_area) %>%
    dplyr::arrange(desc(prop_area))

  print(i)

  return(frac_simpl)

}) %>% dplyr::bind_rows()

readr::write_csv(frac_lme, "data/spatial_data/frac_lme.csv")



# SST data ----------------------------------------------------------------

download.file(
  "https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_sst.nc.gz",
  destfile = "data/spatial_data/HadISST_sst.nc.gz")
R.utils::gunzip("data/spatial_data/HadISST_sst.nc.gz")


## Stock level -------------------------------------------------------------

sst <- ncdf4::nc_open("data/spatial_data/HadISST_sst.nc")

# Extract coordinates
lat <- ncdf4::ncvar_get(sst, varid="latitude")
lon <- ncdf4::ncvar_get(sst, varid="longitude")
range(lat); range(lon) # check -89.5 89.5 -179.5 179.5

# Extract dates
ncdf4::ncatt_get(sst, varid="time", attname="units") # check days since 1870-1-1
start_date <- as.Date("1870-01-01") # must inspect above
dates_orig <- ncdf4::ncvar_get(sst, varid="time")
dates <- start_date + dates_orig

# Extract SST into array
sst_array <- ncdf4::ncvar_get(sst, varid="sst")

# Replace missing values (NA=-1000.0)
sst_array[sst_array==-1000.0] <- NA

# Min and maximum temperatures
tmin <- min(sst_array, na.rm=T)
tmax <- max(sst_array, na.rm=T)

rm(sst)
gc()

# Create SST raster (broken down into four steps to avoid R session crashes)
sst_rast1 <- sst_array[,,1:500] %>%
  terra::rast(crs="+proj=longlat +datum=WGS84",
              extent=c(-90, 90, -180, 180)) %>%
  terra::t()

# Add dates as temp matrix names
names(sst_rast1) <- paste("a", dates[1:500], sep="")

sst_rast2 <- sst_array[,,501:1000] %>%
  terra::rast(crs="+proj=longlat +datum=WGS84",
              extent=c(-90, 90, -180, 180)) %>%
  terra::t()
names(sst_rast2) <- paste("a", dates[501:1000], sep="")

sst_rast3 <- sst_array[,,1001:1500] %>%
  terra::rast(crs="+proj=longlat +datum=WGS84",
              extent=c(-90, 90, -180, 180)) %>%
  terra::t()
names(sst_rast3) <- paste("a", dates[1001:1500], sep="")

sst_rast4 <- sst_array[,,1501:1855] %>%
  terra::rast(crs="+proj=longlat +datum=WGS84",
              extent=c(-90, 90, -180, 180)) %>%
  terra::t()
names(sst_rast4) <- paste("a", dates[1501:1855], sep="")

sst_rast <- c(sst_rast1, sst_rast2, sst_rast3, sst_rast4)

# Calculate mean monthly SST within stock boundaries (zonal stats)
sst.monthly <- terra::extract(sst_rast, stock_boundaries,
                              method="simple", fun=mean, na.rm=TRUE) %>%
  dplyr::select(-ID)

# Merge stock ids with temp matrix
sst.monthly.df <- data.frame(
  assessid=stock_boundaries$assessid,
  stockid=stock_boundaries$stockid,
  sst.monthly)

# Convert temp matrix from wide-to-long
mtemp <- sst.monthly.df %>%
  tidyr::pivot_longer(cols=-c(assessid, stockid),
                      names_to="date", values_to="sst_c") %>%
  dplyr::mutate(date = as.POSIXct(strptime(date, format="%Y.%m.%d")),
                year = as.numeric(format(date, "%Y")),
                month = as.numeric(format(date, "%m"))) %>%
  dplyr::select(assessid, stockid, year, month, date, season, sst_c) %>%
  dplyr::arrange(stockid, date)

# # Convert date to actual date
# mtemp$date <- gsub("a", "", mtemp$date)
# mtemp$date <- as.POSIXct(strptime(mtemp$date, format="%Y.%m.%d"))
# mtemp$year <- as.numeric(format(mtemp$date, "%Y"))
# mtemp$month <- as.numeric(format(mtemp$date, "%m"))
#
# # Rearrange columns
# mtemp <- mtemp %>%
#   dplyr::select(assessid, stockid, year, month, date, sst_c) %>%
#   dplyr::arrange(stockid, date)

# Calculate mean annual SST by stock boundary
atemp <- mtemp %>%
  dplyr::group_by(stockid, assessid, year) %>%
  dplyr::summarize(sst_c=mean(sst_c))

# Calculate mean SST by stock boundary
stemp <- atemp %>%
  dplyr::group_by(stockid, assessid) %>%
  dplyr::summarize(sst_c=mean(sst_c))

readr::write_csv(atemp,
                 "data/spatial_data/yearly_had_sst.csv")


## LME level ---------------------------------------------------------------

# Read LMEs:
lme <- sf::read_sf(dsn = "data/spatial_data/LME66/", layer = "lmes66gcd")
oceans <- sf::read_sf(dsn = "data/spatial_data/ocean_simpl/",
                      layer = "ocean_simpl") %>%
  dplyr::filter(grepl("Ocean", name)) %>%
  dplyr::rename(LME_NAME = name) %>%  # Not using high seas only, not working
  dplyr::mutate(LME_NUMBER = 67:73) # but difference expected to be marginal

lmes <- dplyr::bind_rows(oceans, lme)


# Calculate mean monthly SST within stock boundaries (zonal stats)
sst.monthly.lme <- terra::extract(sst_rast, lmes %>% terra::vect(),
                                  method="simple", fun=mean, na.rm=TRUE) %>%
  dplyr::select(-ID)


# Merge LME names with temp matrix
sst.monthly.lme.df <- data.frame(lme=lmes$LME_NAME, sst.monthly.lme)

# Convert temp matrix from wide-to-long
mtemp.lme <- sst.monthly.lme.df %>%
  tidyr::pivot_longer(cols=-lme, names_prefix="a",
                      names_to="date", values_to="sst_c") %>%
  dplyr::mutate(date = as.POSIXct(strptime(date, format="%Y.%m.%d")),
                year = as.numeric(format(date, "%Y")),
                month = as.numeric(format(date, "%m"))) %>%
  dplyr::select(lme, year, month, date, sst_c) %>%
  dplyr::arrange(lme, date)

# Calculate mean annual SST by LME
atemp.lme <- mtemp.lme %>%
  dplyr::group_by(lme, year) %>%
  dplyr::summarize(sst_c=mean(sst_c))

# Calculate mean SST by LME
stemp.lme <- atemp.lme %>%
  dplyr::group_by(lme) %>%
  dplyr::summarize(sst_c=mean(sst_c))

readr::write_csv(atemp.lme,
                 "data/spatial_data/lme_yearly_had_sst.csv")



# Remove unused raw data left
unlink("data/raw_data/", recursive = TRUE)
