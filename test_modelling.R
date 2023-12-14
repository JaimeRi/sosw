library(hydrographr)
library(data.table)
library(dplyr)
library(terra)
library(tools)
library(stringr)
library(ranger)
library(leaflet)
library(leafem)

bbox=c(-4.1025352478027299,36.5484085083008026,1.9870934486389200,42.7649116516112997)

spdata = read.csv("/home/jaime/data/sosw/fish_roi_modelling.csv")

spdata = spdata[spdata[,2] == "Sander lucioperca",][,c(1:5)]

#tile_id <- get_tile_id(data = spdata, lon = "long", lat = "lat")
####tile_id = c("h16v04", "h18v04")
#tile_id = tile_id[c(2)]


# Variables in raster format
####vars_tif <- c("sub_catchment", "slope_curv_max_dw_cel", "spi", "basin", "segment")

# Download the .tif tiles of the desired variables
####download_tiles(variable = vars_tif, tile_id = tile_id, file_format = "tif",
####               download_dir = "/home/jaime/data/sosw/env_jucar")
####download_tiles(variable = "accumulation", tile_id=tile_id, file_format= "tif", download_dir = "/home/jaime/data/sosw/env_jucar")
# Extend timeout to 1000s to allow uninterrupted downloading
# options(timeout = 1000)
#
# # Download
# # Present, 1981-2010
# ####download.file("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio12_1981-2010_V.2.1.tif",
#               destfile = "/home/jaime/data/sosw/env_jucar/bio12_1981-2010.tif", mode = "wb")
# ####download.file("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio15_1981-2010_V.2.1.tif",
#               destfile = "/home/jaime/data/sosw/env_jucar/bio15_1981-2010.tif", mode = "wb")
# ####download.file("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio1_1981-2010_V.2.1.tif",
#               destfile = "/home/jaime/data/sosw/env_jucar/bio1_1981-2010.tif", mode = "wb")
#
# ####files_chelsa <- list.files("/home/jaime/data/sosw/env_jucar/", pattern = ".tif", full.names = TRUE)
#
# for(ifile in files_chelsa) {
#   crop_to_extent(
#     raster_layer = ifile,
#     bounding_box = bbox,
#     out_dir = "/home/jaime/data/sosw/env_jucar/final_layers",
#     file_name = basename(ifile),
#     read = FALSE,
#     quiet = TRUE)
# }
#
# dirs_h90m = list.dirs(paste0("/home/jaime", "/data/sosw/env_jucar"), recursive = TRUE, full.names = TRUE)
#
# dirs_h90m = dirs_h90m[grep("tiles20d", dirs_h90m)]
#
# for(idir in dirs_h90m) {
#   # only choose rasters
#   tiles <- list.files(idir, pattern = ".tif$", full.names = TRUE)
#   for(itile in tiles) {
#     crop_to_extent(
#       raster_layer = itile,
#       bounding_box = bbox,
#       out_dir = idir,
#       file_name = paste0(str_remove(basename(itile), ".tif"), "_crop.tif"),
#       read = FALSE,
#       quiet = TRUE)
#   }
# }
#
# for(idir in dirs_h90m) {
#   # Get input file extension
#   file_extension <- file_ext(list.files(idir, full.names = FALSE)[1])
#
#   # Assign file extension to output files
#   ivar_name <- paste0(
#     str_remove(basename(idir), "_tiles20d"), ".", file_extension
#   )
#
#   # Run the function
#   merge_tiles(tile_dir = idir,
#               tile_names = list.files(idir, full.names = FALSE,
#                                       pattern = "_crop.tif"),
#               out_dir = "/home/jaime/data/sosw/env_jucar/final_layers",
#               file_name = ivar_name,
#               read = FALSE,
#               bigtiff = TRUE)
# }

layer_dir="/home/jaime/data/sosw/env_jucar/final_layers"

spdata_ids <- extract_ids(data = spdata, lon = "long", lat = "lat",
                          id = "IGB_id", quiet = FALSE,
                          subc_layer = paste0(layer_dir, "/sub_catchment.tif"))

var_layers = list.files(layer_dir)[c(3,4,5,7,8)]
#var_layers = list.files(layer_dir)[c(8)]
report_no_data(data_dir = layer_dir, var_layer = var_layers)

# Run the function that returns the zonal statistics
# stats_table_zon <- extract_zonal_stat(
#   data_dir = layer_dir,
#   subc_layer = paste0(layer_dir, "/sub_catchment.tif"),
#   subc_id = "all",
#   var_layer = var_layers,
#   out_dir = "/home/jaime/data/sosw",
#   file_name = "zonal_stats8.csv",
#   n_cores = 1)

stats_table_zon = read.table("/home/jaime/data/sosw/zonal_stats.csv", sep = " " , header=TRUE)

slope_scale <- function(x, na.rm = F) (x * 0.000001)
clim_scale <- function(x, na.rm = F) (x * 0.1)
offset <- function(x, na.rm = F) (x - 273.15)

stats_table <- stats_table_zon  %>%
  mutate(across(contains("slope_curv_max_dw_cel"), slope_scale)) %>%
  mutate(across(starts_with("bio"), clim_scale))  %>%
  mutate(across(matches("bio1_.*_mean"), offset))

stats_table <- stats_table[!is.na(stats_table$slope_curv_max_dw_cel_mean),]

names(stats_table)[1] = "subcatchment_id"
pseudoabs <- stats_table %>%
  filter(!subcatchment_id %in% spdata_ids$subcatchment_id) %>%
  sample_n(10000)


presence <- left_join(spdata_ids, stats_table, by = "subcatchment_id")
data_model <- data.table::rbindlist(list(presence, pseudoabs), fill = TRUE)

#data_model$occurrence <- ifelse(!is.na(data_model$IGB_ID), 1, 0)

indx = which(!is.na(data_model$IGB_id))
data_model$occurrence = 0
data_model$occurrence[indx] = 1
data_model$occurrence <- as.factor(data_model$occurrence)

#data_model$occurrence <- ifelse(which(!is.na(data_model$IGB_ID)), 1, 0)
#data_model$occurrence <- as.factor(data_model$occurrence)

data_train <- data_model[data_model$subcatchment_id != 0,]
data_train = data_model

# number of presence records:
pres_num <- as.numeric(table(data_train$occurrence)["1"])
sample_size <- c("0" = pres_num / nrow(data_train),
                 "1" = pres_num / nrow(data_train))
sample_size

model <- ranger(data_train$occurrence ~ .,
                data = data_train[, 5:14],
                num.trees = 1000,
                mtry= 6,
                replace = T,
                sample.fraction = sample_size,
                oob.error = T,
                keep.inbag = T,
                num.threads = 4,
                importance ='impurity',
                probability = T)

model


#pred <- predict(model, data = stats_table[,-1])
n <- 100000

nr <- nrow(stats_table)
stats_table_split <- split(stats_table,
                                   rep(1:ceiling(nr/n),
                                       each=n, length.out=nr))

#rm(stats_table, stats_table_present)
gc()

# create the output list
pred_list <- list()

# run prediction in a loop.
# The rbindlist combines the i-item as a datatable for prediction
# [,!1] to remove the subcatchment column that is not needed
loop_length <- length(stats_table_split)

for(i in 1:loop_length) {
  cat("Now predicting on chunk", i, "\n")
  pred_list[[i]] <- predict(model, data = rbindlist(stats_table_split[i]) [,-1])
}

pred_list_combine <- list()
for(i in 1:loop_length) {
  cat("Now combining chunk", i, "\n")
  pred_sub <- setDT(as.data.frame(pred_list[[i]]))
  names(pred_sub) <- c("prob_0", "prob_1")
  pred_list_combine[[i]] <- pred_sub
  rm(pred_sub)
}

# combine
pred_all <- rbindlist(pred_list_combine)
names(pred_all) <- c("prob_0", "prob_1")


## Attach the subc_id
pred_all <- cbind(stats_table[,1], pred_all)

#--- multiply the probability values by 100, to convert them to integers---
pred_all$prob_1 <- as.integer(round(pred_all$prob_1, 2) * 100)

#---reclassify the sub-catchment raster---

reclass_raster(
  data = pred_all,
  rast_val = "V1",
  new_val = "prob_1",
  raster_layer = paste0(layer_dir, "/sub_catchment.tif"),
  recl_layer = paste0("/home/jaime/data/sosw/env_jucar", "/pred_all_Sander_lucioperca.tif"),
  read = FALSE)
