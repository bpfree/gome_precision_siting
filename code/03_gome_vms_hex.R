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
vms_dir <- "data/a_raw_data/GoME_Final_WEA_VMS_data.gdb"

#### study area grid
wea_dir <- "data/a_raw_data/gome_wea/Gulf of Maine WEA"

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

## vector datasets
zmf_function <- function(data){
  
  # calculate minimum value
  min <- min(data$VMS_2009_2021pub, na.rm = T)
  
  # calculate maximum value
  max <- max(data$VMS_2009_2021pub, na.rm = T)
  
  # calculate z-score minimum value
  ## this ensures that no value gets a value of 0
  z_max <- max + (max * 1 / 10000)
  
  # create a field and populate with the value determined by the z-shape membership scalar
  data <- data %>%
    # calculate the z-shape membership value (more desired values get a score of 1 and less desired values will decrease till 0.01)
    ## ***Note: in other words, habitats with higher richness values will be closer to 0
    dplyr::mutate(vms_z = ifelse(VMS_2009_2021pub == min, 1, # if value is equal to minimum, score as 1
                                   # if value is larger than minimum but lower than mid-value, calculate based on scalar equation
                                   ifelse(VMS_2009_2021pub > min & VMS_2009_2021pub < (min + z_max) / 2, 1 - 2 * ((VMS_2009_2021pub - min) / (z_max - min)) ** 2,
                                          # if value is lower than z_maximum but larger than than mid-value, calculate based on scalar equation
                                          ifelse(VMS_2009_2021pub >= (min + z_max) / 2 & VMS_2009_2021pub < z_max, 2 * ((VMS_2009_2021pub - z_max) / (z_max - min)) ** 2,
                                                 # if value is equal to maximum, value is equal to 0.01 [all other values should get an NA]
                                                 ifelse(VMS_2009_2021pub == z_max, 0.01, NA)))))
  
  # return the layer
  return(data)
}

#####################################
#####################################

# load data
## VMS data
vms_density <- sf::st_read(dsn = vms_dir, layer = sf::st_layers(dsn = vms_dir)[[1]][1])

## study region
### final WEA
gome_wea <- sf::st_read(dsn = file.path(wea_dir, "Final_WEA_Poly.shp")) %>%
  # change projection to match all other data coordinate reference system
  sf::st_transform(crs = crs)

## hex grid
gome_hex <- sf::st_read(dsn = gome_gpkg,layer = sf::st_layers(dsn = gome_gpkg)[[1]][grep(pattern = "final_hex",
                                                                                         x = sf::st_layers(dsn = gome_gpkg)[[1]])])

#####################################
#####################################

# VMS density
gome_vms <- vms_density %>%
  rmapshaper::ms_clip(target = .,
                      clip = gome_wea) %>%
  dplyr::select(VMS_2009_2021pub) %>%
  zmf_function(data = .) %>%
  # replace any values that are NA/NULL with a value of 1
  dplyr::mutate(across(vms_z,
                       ~replace(x = .,
                                list = is.na(.),
                                # replacement values
                                values = 1))) %>%
  dplyr::relocate(geometry, .after = vms_z)

gome_vms_hex <- gome_hex[gome_vms, ] %>%
  # spatially join North Atlantic right whale values to Gulf of Maine WEA hex cells
  sf::st_join(x = .,
              y = gome_vms,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(OBJECTID,
                VMS_2009_2021pub, vms_z) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(OBJECTID) %>%
  # summarise the North Atlantic right whale score values
  ## take the maximum value of the North Atlantic right whale score for any that overlap
  ## ***Note: this will provide the most conservation given that
  ##          high values are less desirable
  dplyr::summarise(vms_min = min(vms_z))

# gome_vms_hex <- vms_density %>%
#   sf::st_join(x = .,
#               y = gome_hex,
#               join = st_intersects) %>%
#   dplyr::select(OBJECTID, VMS_2009_2021pub) %>%
#   # group by the index values as there are duplicates
#   dplyr::group_by(OBJECTID) %>%
#   dplyr::summarise(vms_min = min(VMS_2009_2021pub))

# rescale data with z-shaped membership function
# gome_vms_z_hex <- gome_vms_hex %>%
#   sf::st_cast(x = .,
#               to = "POLYGON") %>%
#   zmf_function(data = .) %>%
#   # replace any values that are NA/NULL with a value of 1
#   dplyr::mutate(across(z_value,
#                        ~replace(x = .,
#                                 list = is.na(.),
#                                 # replacement values
#                                 values = 1)))

#####################################
#####################################

# export data
sf::st_write(obj = gome_vms_hex, dsn = gome_gpkg, layer = "gome_vms_hex", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate