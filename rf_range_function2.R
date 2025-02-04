library(dplyr)
library(ranger)
library(hydrographr)
library(modEvA)
library(PresenceAbsence)
#library(terra)

# Function from https://rdrr.io/rforge/biomod2/src/R/SampleMat2.R
# Separates presence-absence data in two subsets for calibration and evaluation

SampleMat2 <- function(ref, ratio, as.logi=FALSE)
  {
    # set a new random seed to ensure that sampling is random (issue when CTA is involved and seed needs to be set to a fix number)
    set.seed(as.double(Sys.time()) + as.numeric(format(Sys.time(), "%OS6"))*1000000)
    
    ntot <- length(ref)
    npres<- sum(ref)    
    ncal <- ceiling(ntot*ratio)
    
    pres <- sample(which(ref==1), ceiling(npres*ratio))
    absc <- sample(which(ref==0), ncal-length(pres))
    
    if(as.logi){
      calib <- rep(FALSE, ntot)
      calib[c(pres,absc)] <- TRUE
      eval <- !calib
    } else{
      calib <- c(pres,absc)
      eval <- (1:ntot)[-c(pres,absc)]
    }
    
    return(list("calibration"=calib, "evaluation"=eval))
  }


# Make a folder to save the list if the folder doesn't exist 
#folder_path <- "/mnt/shared/danube/models/output/rf"
folder_path <- "/mnt/shared/sosw/bio_modelling/rhein/randomForest_output/"

if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
  cat(paste0("Folder created successfully.", "\n"))
} else {
  cat(paste0("Folder already exists.", "\n"))
}

# sub-catchment IDs Danube
load("/mnt/shared/danube/models/input/subcID_danube.RData")

# Predict table path
cat(paste0("Importing table with environmental predictors", "\n"))
predict_table_path <- "/mnt/shared/sosw/bio_modelling/danube/pc_scores.RData"
load(predict_table_path)

pred_table <- pc_signif |>
  filter(subcID %in% subcID_danube)

# environmental data for 
env_data <- pc_signif

# Import presences (coordinates)
occur <- read.csv("/mnt/shared/danube/models/input/danube_records_final_20-03-2024.csv") |>
  dplyr::select(subcatchment_id,species)

# Species list
all_species <- occur |>
  distinct(species) |>
  pull()

# Sub-catchment IDs in the study area
ids_study_area <- env_data |>
  pull(subcID)

for (sp in all_species) {
  
  cat(paste0("Start modeling ",sp,': ', which(all_species %in% sp) , 
              ' of ', length(all_species)), " species", "\n")
  
  # Record the start time
  start <- Sys.time() 

  # IDs for presences
  ids_p <- occur |>
    dplyr::filter(species == sp) |>
    pull(subcatchment_id)
  
  # Generate pseudo-absences. Randomly, use sub-catchment IDs. 10000 pseudo-absences
  ids_a <-  sample(ids_study_area[!ids_study_area %in% ids_p],
                   10000)
  
  # Bind presences and pseudo-absences
  pa <- data.frame("subcID" = c(ids_p, ids_a),
                   "occurrence" =  c(rep(1, length(ids_p)), rep(0, length(ids_a))))
    
  
  # Make data frame with environmental variables for presences and pseudo-absences
  env_data_pa <- env_data |>
    inner_join(pa, by = c("subcID" = "subcID"))
  
  # separate the original data in one subset  for training (80 % of records) and the other for 
  #  validation (20 % of records). 
  a <- SampleMat2(ref = env_data_pa$occurrence, ratio = 0.8) 
  data_train <- env_data_pa[a$calibration, ] 
  data_test <- env_data_pa[a$evaluation, ]
  
  # number of presence records:
  pres_num <- as.numeric(table(data_train$occurrence)["1"])
  sample_size <- c("0" = pres_num / nrow(data_train),
                   "1" = pres_num / nrow(data_train))
  
  # convert the response to factor for RF model to return probabilities
  data_train$occurrence <- as.factor(data_train$occurrence)
  
  # train model
  cat(paste0("Train model for ",sp, "\n") )
  model <- NULL
  model <- ranger(data_train$occurrence ~ .,
                  data = data_train |> dplyr::select(-subcID),
                  num.trees = 1000,
                  mtry= 6,
                  replace = T,
                  sample.fraction = sample_size,
                  oob.error = T,
                  keep.inbag = T,
                  num.threads = 4,
                  importance ='impurity',
                  probability = T)
  
  # test the model
    # predict on test data
    pred_test <- predict(model, data = data_test |> dplyr::select(-subcID, -occurrence))
    
    # Calculate AUC. Use AUC to asses model predictive ability
    auc <- AUC(obs = data_test$occurrence,
               pred = as.numeric(pred_test$predictions[,2]),
               plot = FALSE)$AUC
    
    # Calculate threshold
    th <- data.frame("ids" = data_test$subcID,
                     "obs" = data_test$occurrence,
                      "pred" = as.numeric(pred_test$predictions[,2])) |>
      optimal.thresholds(opt.methods = c("Sens=Spec"))
    
    # Calculate threshold dependent metrics
    th_metrics <- threshMeasures(obs = data_test$occurrence, 
                                 pred = as.numeric(pred_test$predictions[,2]),
                                 simplif = T, 
                                 thresh = th[,2], 
                                 standardize = T,
                                 measures = c("Sensitivity",
                                              "Specificity", 
                                              "Omission", 
                                              "Commission", 
                                              "TSS"))
  
  # predict
   cat(paste0("Prediction for ",sp, "\n"))
  pred <- predict(model, data = pred_table|> dplyr::select(-subcID))
  
  # extract predictions and make a data frame with sub-catchment IDs
  prediction <- data.frame("subcID" = pred_table$subcID,
                           "p" = as.numeric(pred$predictions[,2]))
  
  
  # Path to raster layer with sub-catchments
  path_subc <- "/mnt/shared/danube/models/input/sub_catchment_danube.tif"
  
  # Output path
  path_map_output <- paste0(
                            "/mnt/shared/danube/models/output/rf/",
                            sub(" ", "_", sp),
                            "_rf.tif")
  
  # Multiply the probability values by 100, to convert them to integers
  prediction4map <- prediction |>
    mutate(p = as.integer(round(prediction$p, 2) * 100))
  
  # Reclassify sub_catchments with predictions and save map
  # cat(paste0("Making map with predictions for ",sp, "\n"))
  # reclass_raster(
  #   data = prediction4map,
  #   rast_val = "subcID",
  #   new_val = "p",
  #   raster_layer = path_subc,
  #   recl_layer = path_map_output,
  #   read = FALSE,
  #   quiet = FALSE
  # )
  
  # list with evaluation metrics, threshold and predictions
  results <- list("auc" = auc,
                  "threshold" = th,
                  "threshold_metrics" = th_metrics,
                  "predictions" = prediction,
                  "predictions4map" = prediction4map,
                  "variable_importance" =  model$variable.importance)
  # export list 
  cat(paste0("Exporting results for ",sp, "\n"))
  save(results, file = paste0(folder_path,"/",
                       sub(" ", "_", sp),
                       "_rf.RData")
       )
  # free unused memory
  gc()
  
  # Record the end time
  end <- Sys.time() 
  
  cat(paste0("Execution time for ", sp, ": ", end - start, "\n"))
  
}


