#########
### Gulf of Maine model score calculations ###

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
## submodel raw data directory
data_dir <- "data/a_raw_suitability/gome_precision_siting.gdb"

## Gulf of Maine geopackage (with model scores)
gome_model_geopackage <- "data/b_model_scores/gome_model_scores.gpkg"

#####################################
#####################################

# load data
gome_data <- sf::st_read(dsn = data_dir, layer = sf::st_layers(dsn = data_dir)[[1]][1]) %>%
  # remove any unneeded fields -- select all that do not start with "Shape" -- so will exclude "Shape_Length" and "Shape_Area"
  dplyr::select(!starts_with("Shape")) %>%
  # calculate model scores
  ## model 1: VMS only data
  ## model 2: geometric mean of VMS and North Atlantic right whale
  ## model 3: geometric mean of VMS, North Atlantic right whale, and NREL LCOE data
  dplyr::mutate(model1 = VMS_2009_2021_Z,
                model2 = (VMS_2009_2021_Z ** model2) * (NARWden_Z ** model2),
                model3 = (VMS_2009_2021_Z ** model3) * (NARWden_Z ** model3) * (lcoe_2035_GW_DL ** model3)) %>%
  dplyr::relocate(Shape, .after = model3)

#####################################
#####################################

# export data
gome_data <- sf::st_write(obj = gome_data, dsn = gome_model_geopackage, layer = paste(region, "model", sep = "_"), append = F)
