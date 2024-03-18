############################
### 03. VMS density data ###
############################

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
#### VMS density (2009 - 2021)
vms_dir <- "data/a_raw_data/gome_fisheries_submodel.gpkg"

#### hex grid
gome_gpkg <- "data/b_model_scores/gome_model_scores.gpkg"

#####################################
#####################################

# set parameters
## coordinate reference system
### EPSG:26918 is NAD83 / UTM 18N (https://epsg.io/26918)
crs <- "EPSG:26919"

#####################################
#####################################

# load data
## VMS data
vms_density <- sf::st_read(dsn = vms_dir, layer = sf::st_layers(dsn = vms_dir)[[1]][1])

testy <- sf::st_read(dsn = vms_dir, layer = sf::st_layers(dsn = vms_dir)[[1]][1]) %>%
  sf::st_drop_geometry()

## hex grid
gome_hex <- sf::st_read(dsn = gome_gpkg,
                        layer = sf::st_layers(dsn = gome_gpkg)[[1]][grep(pattern = "final_hex",
                                                                         x = sf::st_layers(dsn = gome_gpkg)[[1]])])

#####################################
#####################################

# VMS hex grids
test <- gome_hex %>%
  dplyr::left_join(x = .,
                   y = testy,
                   by = "OBJECTID")

gome_vms_hex <- vms_density %>%
  sf::st_join(x = .,
              y = gome_hex,
              join = st_intersects) %>%
  dplyr::select(OBJECTID.y, VMS_2009_2021_Z) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(OBJECTID.y) %>%
  dplyr::summarise(vms_min = min(VMS_2009_2021_Z)) %>%
  dplyr::rename("OBJECTID" = "OBJECTID.y") %>%
  dplyr::select(OBJECTID, vms_min)


# need to jon hex values back with overall grid

# replace any values that are NA/NULL with a value of 1
dplyr::mutate(across(2:6, ~replace(x = .,
                                   list = is.na(.),
                                   # replacement values
                                   values = 1)))

#####################################
#####################################

# export data
sf::st_write(obj = gome_vms_hex, dsn = gome_gpkg, layer = "gome_vms_hex")

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate