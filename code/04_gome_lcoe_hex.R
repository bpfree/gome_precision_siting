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

#### hex grid
gome_gpkg <- "data/b_model_scores/gome_model_scores.gpkg"


### Output directories

#####################################
#####################################

# Set parameters
## designate region name
region <- "gome"

## layer names
layer <- "lcoe_2025"

#####################################
#####################################

# Load data
## hex grid
gome_hex <- sf::st_read(dsn = gome_gpkg,
                        layer = sf::st_layers(dsn = gome_gpkg)[[1]][grep(pattern = "final_hex",
                                                                         x = sf::st_layers(dsn = gome_gpkg)[[1]])])

#####################################

## Levelized cost of energy (2025)
lcoe <- sf::st_read(dsn = lcoe_dir,
                        layer = sf::st_layers(dsn = lcoe_dir)[[1]][grep(pattern = "lcoe",
                                                                         x = sf::st_layers(dsn = lcoe_dir)[[1]])])

#####################################
#####################################

# Energy cost in Oregon call areas
gome_lcoe <- lcoe %>%
  sf::st_join(x = .,
              y = gome_hex,
              join = st_intersects) %>%
  dplyr::select(OBJECTID, lcoe_2035_GW) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(OBJECTID) %>%
  dplyr::summarise(lcoe_min = min(lcoe_2035_GW))

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
# min <- min(gome_lcoe$narw_min)
# 
# ### value maximum
# max <- max(gome_lcoe$narw_min)
# 
# ## Gulf of Maine levelized cost of energy normalized
# gome_lcoe_2035_norm <- gome_lcoe %>%
#   dplyr::mutate(lcoe_norm = (narw_min - min) / (max - min),
#                 layer = "levelized cost of energy (2035)") %>%
#   # select fields of interest
#   dplyr::select(OBJECTID,
#                 layer,
#                 narw_min,
#                 lcoe_norm)

gome_lcoe_2035_norm <- gome_lcoe %>%
  # normalize data between 0.000858 and 1
  dplyr::mutate(lcoe_norm = 
                  # flip the normalization so high costs get low
                  # scores (0.000858) and low costs get higher scores (1.0)
                  (tmax + tmin) - 
                  # normalize the data (will become between 0.8 and 1.0)
                  (((lcoe_mid - min) / (max - min)) *
                     # and then rescale to be between 0.8 and 1.0
                     (tmax - tmin) + tmin)) %>%
  # create field called "layer" and fill with "levelized cost of energy (2035)" for summary
  dplyr::mutate(layer = "levelized cost of energy (2035)") %>%
  # select fields of interest
  dplyr::select(OBJECTID,
                layer,
                narw_min,
                lcoe_norm)

#####################################
#####################################

# Levelized cost of energy (2025) hex grid
gome_hex_lcoe_2025 <- gome_hex[gome_lcoe_2025_norm, ] %>%
  # spatially join continental shelf values to Gulf of Maine final WEA hex cells
  sf::st_join(x = .,
              y = gome_lcoe_2025_norm,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(index, layer,
                lcoe_norm) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(index) %>%
  # summarise the normalized values for the levelized cost of energy in 2025
  ## take the minimum value of the index
  dplyr::summarise(lcoe_norm_index = min(lcoe_norm))

#####################################
#####################################

# Export data
## Wind submodel
sf::st_write(obj = gome_hex_lcoe_2025, dsn = wind_submodel, layer = paste0(region, "_hex_", layer), append = F)

## Wind geopackage
sf::st_write(obj = gome_lcoe_2025, dsn = wind_gpkg, layer = "gome_lcoe_2025", append = F)
sf::st_write(obj = gome_lcoe_2025_norm, dsn = wind_gpkg, layer = "gome_lcoe_2025_norm", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate