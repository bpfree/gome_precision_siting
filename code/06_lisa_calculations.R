##############################################
### 06. Local Index of Spatial Association ###
##############################################

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

# set directories
## lcoe
## Gulf of Maine geopackage (with model scores)
gome_gpkg <- "data/b_model_scores/gome_model_scores.gpkg"

#####################################
#####################################

# Inspect available layers and names within suitability model geopackage
sf::st_layers(dsn = gome_gpkg,
              do_count = F)

#####################################
#####################################

# Set parameters
## designate region name
region <- "gome"

## lisa
layer <- "lisa"
classification <- "highhigh"

## distance threshold
distance <- 250

## model
model <- "model3"

## permutations
permutations <- 9999

## confidence interval
ci_cutoff <- 0.05

## designate date
date <- format(Sys.time(), "%Y%m%d")

#####################################
#####################################

# Load data
## Gulf of Maine suitability areas
model_areas_sf <- sf::st_read(dsn = gome_gpkg, layer = "gome_model")

#####################################
#####################################

# Calculate local indicator of spatial autocorrelation (LISA)
## ***NOTE: Multiple packages exist for calculating in R. rgeoda was
##          selected as it was developed by the team of Luc Anselin,
##          who developed the LISA methodology and GeoDA software.

## Create distance weights (250 meters)
### Examine the minimum distance for the layer to have at least
### a single neighbor adjacent
dist_thres <- rgeoda::min_distthreshold(model_areas_sf)
dist_thres # every hex will have at least a single neighbor at a distance of 216.1691 meters

### In Gulf of Maine, it was previously determine that a distance of 250 meters
### returned the ideal results
weights <- rgeoda::distance_weights(sf_obj = model_areas_sf,
                                         # set distance weight to be 250 meters
                                         dist_thres = distance,
                                         # set false to calculate distance in a Euclidean space
                                         is_arc = F,
                                         # set false to keep in metric system (units are in miles)
                                         is_mile = F)

## Create the LISA product
start <- Sys.time() # set start time
lisa <- rgeoda::local_moran(w = weights, # weight is equal to distance weight
                            # analyze using the final model geometric mean values
                            df = model_areas_sf[model],
                            # run permutations to generate pseudo-p-values (default permutations is 999)
                            permutations = permutations,
                            # set cutoff to be significance limit
                            significance_cutoff = ci_cutoff,
                            # set seed so future runs can compare (default is 123456789)
                            seed = 561974)
print(Sys.time() - start) # print how long it takes to calculate the LISA results

#####################################

## Gather important fields
### p-values
maine_lisa_pvalues <- lisa$p_vals %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("p_vals" = ".")

### cluster values
maine_lisa_cvalues <- lisa$c_vals %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("c_vals" = ".")

### LISA values
maine_lisa_lisa_values <- lisa$lisa_vals %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("lisa_vals" = ".")

### number of neighbors
maine_lisa_neighbors <- lisa$nn_vals %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("num_neigh" = ".")

### labels
maine_lisa_labels <- lisa$labels %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("labels" = ".") %>%
  # add field to join with cluster values
  ## Predefined values
  ### 0 Not significant
  ### 1 High-High
  ### 2 Low-Low
  ### 3 Low-High
  ### 4 High-Low
  ### 5 Undefined
  ### 6 Isolated
  dplyr::mutate(c_vals = c(0, 1, 2, 3, 4, 5, 6))

### colors
maine_lisa_colors <- lisa$colors %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("colors" = ".") %>%
  # add field to join with cluster values
  ## Predefined values
  ### 0 #eeeeee (not significant)
  ### 1 #FF0000 (high-high)
  ### 2 #0000FF (low-low)
  ### 3 #a7adf9 (low-high)
  ### 4 #f4ada8 (high-low)
  ### 5 #464646 (undefined)
  ### 6 #999999 (isolated)
  dplyr::mutate(c_vals = c(0, 1, 2, 3, 4, 5, 6))

### joined cluster values with respective labels and colors
maine_lisa_cvalue_labels_colors <- maine_lisa_cvalues %>%
  # join cluster values with labels using cluster value
  dplyr::left_join(x = .,
                   # label dataset
                   y = maine_lisa_labels,
                   # join field is cluster values
                   by = "c_vals") %>%
  # join cluster values with colors using cluster value
  dplyr::left_join(x = .,
                   # colors dataset
                   y = maine_lisa_colors,
                   # join field is cluster values
                   by = "c_vals")

#####################################

# false discovery rate
maine_fdr <- rgeoda::lisa_fdr(gda_lisa = lisa, current_p = 0.05)
maine_bo <- rgeoda::lisa_bo(gda_lisa = lisa, current_p = 0.05)

#####################################

# Gulf of Maine hex grid joined with all new fields for LISA
maine_hex_lisa <- model_areas_sf %>%
  cbind(maine_lisa_pvalues,
        maine_lisa_cvalue_labels_colors,
        maine_lisa_lisa_values,
        maine_lisa_neighbors)

#####################################
#####################################

combinations <- maine_hex_lisa %>%
  dplyr::group_by(labels,
                  colors) %>%
  dplyr::summarise()
combinations

maine_lisa_highhigh <- maine_hex_lisa %>%
  dplyr::filter(labels == "High-High")

significance <- maine_hex_lisa %>%
  dplyr::mutate(fdr = maine_fdr) %>%
  dplyr::mutate(fdr_cluster = rgeoda::lisa_clusters(gda_lisa = lisa, cutoff = maine_fdr))

maine_new_labels <- lisa$labels %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("labels" = ".") %>%
  # add field to join with cluster values
  ## Predefined values
  ### 0 Not significant
  ### 1 High-High
  ### 2 Low-Low
  ### 3 Low-High
  ### 4 High-Low
  ### 5 Undefined
  ### 6 Isolated
  dplyr::mutate(fdr_cluster = c(0, 1, 2, 3, 4, 5, 6))

significance_new_label <- significance %>%
  dplyr::left_join(x = .,
                   y = maine_new_labels,
                   by = "fdr_cluster")

#####################################
#####################################

label_colors <- c("#eeeeee", # Not significant
                  "#f4ada8", # High-Low
                  "#0000FF", # Low-Low
                  "#FF0000", # High-High
                  "#a7adf9") # Low-High

g <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = maine_lisa_highhigh, aes(fill = "#FF0000"), color = NA)
g

#####################################
#################################f####

# Export data
## LISA
sf::st_write(obj = maine_lisa_highhigh, dsn = gome_gpkg, layer = paste(region, layer, model, ci_cutoff, classification, sep = "_"), append = F)
sf::st_write(obj = maine_hex_lisa, dsn = gome_gpkg, layer = paste(region, "hex", layer, model, ci_cutoff, sep = "_"), append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate