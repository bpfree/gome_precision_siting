##############################################
### 01. North Atlantic right whale density ###
##############################################

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
#### North Atlantic right whale yearly density
narw_dir <- "data/a_raw_data/narw_yr_density"

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

# function
## z-membership function
### Adapted from https://www.mathworks.com/help/fuzzy/zmf.html
zmf_function <- function(raster){
  # calculate minimum value
  min <- terra::minmax(raster)[1,]
  
  # calculate maximum value
  max <- terra::minmax(raster)[2,]
  
  # calculate z-score minimum value
  ## this ensures that no value gets a value of 0
  z_max <- max + (max * 1 / 10000)
  
  # calculate z-scores (more desired values get score of 1 while less desired will decrease till 0)
  z_value <- ifelse(raster[] == min, 1, # if value is equal to minimum, score as 1
                    # if value is larger than minimum but lower than mid-value, calculate based on reduction equation
                    ifelse(raster[] > min & raster[] < (min + z_max) / 2, 1 - 2 * ((raster[] - min) / (z_max - min)) ** 2,
                           # if value is larger than mid-value but lower than maximum, calculate based on equation
                           ifelse(raster[] >= (min + z_max) / 2 & raster[] < z_max, 2*((raster[] - z_max) / (z_max - min)) ** 2,
                                  # if value is equal to maximum, score min - (min * 1 / 1000); otherwise give NA
                                  ifelse(raster[] == z_max, 0, NA))))
  
  # set values back to the original raster
  zvalues <- terra::setValues(raster, z_value)
  
  # return the raster
  return(zvalues)
}

#####################################
#####################################

# load data
## NARW density data
narw_density <- terra::rast(paste(narw_dir, "narw_yr_density.tif", sep = "/"))

### inspect data
#### plot data
plot(narw_density)

#### minimum and maximum values
terra::minmax(narw_density)[1]
terra::minmax(narw_density)[2]

### coordinate reference system
cat(crs(narw_density))

#####################################

## study region
### final WEA
gome_wea <- sf::st_read(dsn = file.path(wea_dir, "Final_WEA_Poly.shp")) %>%
  # change projection to match all other data coordinate reference system
  sf::st_transform(crs = crs) %>%
  ### buffered WEA (so all required hex grid cells get covered)
  sf::st_buffer(dist = 1600)

### Inspect study region coordinate reference system
cat(crs(gome_wea))

## hex grid
gome_hex <- sf::st_read(dsn = gome_gpkg,
                        layer = sf::st_layers(dsn = gome_gpkg)[[1]][grep(pattern = "final_hex",
                                                                         x = sf::st_layers(dsn = gome_gpkg)[[1]])])

### Inspect study region coordinate reference system
cat(crs(gome_hex))

#####################################
#####################################

# limit data to study region
gome_narw <- terra::crop(x = narw_density,
                         # crop using study region
                         y = gome_wea,
                         # mask using study region (T = True)
                         mask = T)

plot(gome_narw)

#####################################
#####################################

# rescale AIS values using a z-membership function
narw_z <- gome_narw %>%
  zmf_function()

## inspect rescaled data
plot(narw_z)

#####################################
#####################################

# convert raster to vector data (as polygons)
# convert to polygon
gome_narw_polygon <- terra::as.polygons(x = narw_z,
                                           # do not aggregate all similar values together as single feature
                                           aggregate = F,
                                           # use the values from original raster
                                           values = T) %>%
  # change to simple feature (sf)
  sf::st_as_sf() %>%
  # simplify column name to "narw" (this is the first column of the object, thus the colnames(.)[1] means take the first column name from the ais object)
  dplyr::rename(narw = colnames(.)[1]) %>%
  # add field "layer" and populate with "narw"
  dplyr::mutate(layer = "narw") %>%
  # limit to the study region
  rmapshaper::ms_clip(clip = gome_wea) %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 18N) that will convert units from degrees to meters
  sf::st_transform(crs = crs)

## inspect vectorized rescaled North Atlantic right whale data (***warning: lots of data, so will take a long time to load; comment out unless want to display data)
plot(gome_narw_polygon)

#####################################
#####################################

# North Atlantic right whale hex grids
gome_narw_hex <- gome_hex[gome_narw_polygon, ] %>%
  # spatially join North Atlantic right whale values to Gulf of Maine WEA hex cells
  sf::st_join(x = .,
              y = gome_narw_polygon,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(OBJECTID,
                layer, narw) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(OBJECTID) %>%
  # summarise the North Atlantic right whale score values
  ## take the maximum value of the North Atlantic right whale score for any that overlap
  ## ***Note: this will provide the most conservation given that
  ##          high values are less desirable
  dplyr::summarise(narw_min = min(narw))

#####################################
#####################################

# export data
sf::st_write(obj = gome_narw_hex, dsn = gome_gpkg, layer = "gome_narw_hex", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate