####################################
### 04. Levelized Cost of Energy ###
####################################

# Clear environment
rm(list = ls())

# Calculate start time of code (determine how long it takes to complete all code)
start <- Sys.time()

#####################################
#####################################

# Load packages
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
               rgdal,
               rgeoda,
               rgeos,
               rmapshaper,
               rnaturalearth, # use devtools::install_github("ropenscilabs/rnaturalearth") if packages does not install properly
               sf,
               sp,
               stringr,
               terra, # is replacing the raster package
               tidyr)

#####################################
#####################################

# Set directories
## Define data directory (as this is an R Project, pathnames are simplified)
### Input directories
lcoe_dir <- "data/a_raw_data/lcoe/GoME_Final_WEA_LCOE_data.gdb"

#### study area grid
wea_dir <- "data/a_raw_data/gome_wea/Gulf of Maine WEA"

#### hex grid
gome_gpkg <- "data/b_model_scores/gome_model_scores.gpkg"

#####################################
#####################################

# Set parameters
## designate region name
region <- "gome"

# set parameters
## coordinate reference system
### EPSG:26918 is NAD83 / UTM 18N (https://epsg.io/26918)
crs <- "EPSG:26919"

## layer names
layer <- "lcoe_2035"

#####################################
#####################################

# Load data
## hex grid
gome_hex <- sf::st_read(dsn = gome_gpkg, layer = sf::st_layers(dsn = gome_gpkg)[[1]][grep(pattern = "final_hex",
                                                                                         x = sf::st_layers(dsn = gome_gpkg)[[1]])])

## study region
### final WEA
gome_wea <- sf::st_read(dsn = file.path(wea_dir, "Final_WEA_Poly.shp")) %>%
  # change projection to match all other data coordinate reference system
  sf::st_transform(crs = crs)

## Levelized cost of energy (2035)
lcoe <- sf::st_read(dsn = lcoe_dir, layer = sf::st_layers(dsn = lcoe_dir)[[1]][grep(pattern = "lcoe",
                                                                                    x = sf::st_layers(dsn = lcoe_dir)[[1]])])
#####################################
#####################################

# Energy cost in Oregon call areas
gome_lcoe <- lcoe %>%
  rmapshaper::ms_clip(gome_wea) %>%
  dplyr::select(lcoe_2035_GW)

#####################################

# Normalize the levelized cost data
## The analysis seeks to have values get normalized with a linear function
## and place them between 0.8 and 1.0
### Linear function = (value - min) / (max - min)
### Linear function between values = ((value - min) / (max - min)) * (target_max - target_min) + target_min

## Targets
### target minimum
tmin <- 0.000858

### target maximum
tmax <- 1.0

## Values
### value minimum
min <- min(gome_lcoe$lcoe_2035_GW)

### value maximum
max <- max(gome_lcoe$lcoe_2035_GW)

gome_lcoe_norm <- gome_lcoe %>%
  # normalize data between 0.000858 and 1
  dplyr::mutate(lcoe_norm = 
                  # flip the normalization so high costs get low
                  # scores (0.000858) and low costs get higher scores (1.0)
                  (tmax + tmin) - 
                  # normalize the data (will become between 0.8 and 1.0)
                  (((lcoe_2035_GW - min) / (max - min)) *
                     # and then rescale to be between 0.8 and 1.0
                     (tmax - tmin) + tmin)) %>%
  # select fields of interest
  dplyr::select(lcoe_2035_GW,
                lcoe_norm)

#####################################
#####################################

# # Levelized cost of energy (2035) hex grid
gome_lcoe_hex <- gome_hex[gome_lcoe_norm, ] %>%
  # spatially join continental shelf values to Gulf of Maine final WEA hex cells
  sf::st_join(x = .,
              y = gome_lcoe_norm,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(OBJECTID,
                lcoe_norm) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(OBJECTID) %>%
  # summarise the normalized values for the levelized cost of energy in 2035
  ## take the minimum value of the index
  dplyr::summarise(lcoe_norm = min(lcoe_norm))

#####################################
#####################################

# Export data
## Wind submodel
sf::st_write(obj = gome_lcoe_hex, dsn = gome_gpkg, layer = paste(region, layer, "hex", sep = "_"), append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
