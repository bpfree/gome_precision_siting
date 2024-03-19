##################################################
### 05. Gulf of Maine model score calculations ###
##################################################

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

# region
region <- "gome"

# model geometric means
model2 <- 1/2
model3 <- 1/3

#####################################
#####################################

# set directories
## lcoe
## Gulf of Maine geopackage (with model scores)
gome_model_geopackage <- "data/b_model_scores/gome_model_scores.gpkg"


#####################################
#####################################

vms <- sf::st_read(dsn = gome_model_geopackage, layer = sf::st_layers(dsn = gome_model_geopackage)[[1]][grep(pattern = "vms",
                                                                                 x = sf::st_layers(dsn = gome_model_geopackage)[[1]])])
narw <- sf::st_read(dsn = gome_model_geopackage, layer = sf::st_layers(dsn = gome_model_geopackage)[[1]][grep(pattern = "narw",
                                                                           x = sf::st_layers(dsn = gome_model_geopackage)[[1]])])
lcoe <- sf::st_read(dsn = gome_model_geopackage, layer = sf::st_layers(dsn = gome_model_geopackage)[[1]][grep(pattern = "lcoe",
                                                                                    x = sf::st_layers(dsn = gome_model_geopackage)[[1]])])

gome_data <- vms %>%
  cbind(narw,
        lcoe) %>%
  dplyr::select(OBJECTID,
                z_value,
                narw_min,
                lcoe_norm) %>%
  dplyr::rename(vms_min = z_value)

# load data
gome_model <- gome_data %>%
  # calculate model scores
  ## model 1: VMS only data
  ## model 2: geometric mean of VMS and North Atlantic right whale
  ## model 3: geometric mean of VMS, North Atlantic right whale, and NREL LCOE data
  dplyr::mutate(model1 = vms_min,
                model2 = (vms_min ** model2) * (narw_min ** model2),
                model3 = (vms_min ** model3) * (narw_min ** model3) * (lcoe_norm ** model3)) %>%
  dplyr::relocate(geom, .after = model3)

#####################################
#####################################

# export data
sf::st_write(obj = gome_model, dsn = gome_model_geopackage, layer = paste(region, "model", sep = "_"), append = F)
