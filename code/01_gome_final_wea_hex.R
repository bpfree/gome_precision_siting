#####################
### 01. Final WEA ###
#####################

# clear environment
rm(list = ls())

# calculate start time of code (determine how long it takes to complete all code)
start <- Sys.time()

#####################################
#####################################

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(docxtractr,
               dplyr,
               elsa,
               fasterize,
               fs,
               ggplot2,
               janitor,
               ncf,
               paletteer,
               pdftools,
               plyr,
               purrr,
               raster,
               RColorBrewer,
               reshape2,
               rgeoda,
               rmapshaper,
               rnaturalearth, # use devtools::install_github("ropenscilabs/rnaturalearth") if packages does not install properly
               sf,
               sp,
               stringr,
               terra, # is replacing the raster package
               tidyr)

#####################################
#####################################

# set directories
## define data directory (as this is an R Project, pathnames are simplified)
### input directories
#### study area grid
wea_dir <- "data/a_raw_data/gome_wea/Gulf of Maine WEA"

#### hex grid
hex_dir <- "data/a_raw_data/gome_hex_grid"

### output directories
gome_gpkg <- "data/b_model_scores/gome_model_scores.gpkg"

#####################################
#####################################

# set parameters
## designate region name
region <- "gome"

## coordinate reference system
### EPSG:26918 is NAD83 / UTM 18N (https://epsg.io/26918)
crs <- "EPSG:26919"

#####################################
#####################################

## study region
gome_wea <- sf::st_read(dsn = file.path(wea_dir, "Final_WEA_Poly.shp")) %>%
  # change projection to match AIS data coordinate reference system
  sf::st_transform(crs = crs)

### Inspect study region coordinate reference system
cat(crs(gome_wea))

## hex grid
gome_hex <- sf::st_read(dsn = file.path(hex_dir, "GoME_Grid_10ac.shp")) %>%
  dplyr::select(OBJECTID)

### Inspect study region coordinate reference system
cat(crs(gome_hex))

#####################################
#####################################

## get hexes that intersect with the final WEA
gome_hex_final <- gome_hex[gome_wea,]

#####################################
#####################################

# export data
sf::st_write(obj = gome_hex_final, dsn = gome_gpkg, layer = "gome_wea_final_hex")

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate