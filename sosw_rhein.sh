#! /bin/bash

export DIR=/mnt/shared/sosw
export TMP=/mnt/shared/sosw/tmp


### STEP 1: download to a folder the environmental variables of interest

TILE=(h18v02)

for TL in ${TILE[@]}
do

	for file in $(find /mnt/shared/Environment90m_v.1.0_online/LandCover/ -name "*_2020_${TL}.zip")
	do 
	sudo unzip -j $file -d /mnt/shared/sosw/env_tl
    done

	for file in $(find /mnt/shared/Environment90m_v.1.0_online/Climate/present -name "*_${TL}.zip")
	do 
	sudo unzip -j $file -d /mnt/shared/sosw/env_tl
	done
	
	for file in $(find /mnt/shared/Environment90m_v.1.0_online/Hydrography90m -name "*_${TL}.zip")
	do 
	sudo unzip -j $file -d /mnt/shared/sosw/env_tl
	done
	
	for file in $(find /mnt/shared/Environment90m_v.1.0_online/Soil -name "*_${TL}.zip")
	do 
	sudo unzip -j $file -d /mnt/shared/sosw/env_tl
	done
	
done


### STEP 2: download subcatchment file and prepare for the ROI

# path to subcatchments
subc="/mnt/shared/hydrography90m_v.1.0_online/hydrography90m_v.1.0/r.watershed/sub_catchment_tiles20d/sub_catchment_h18v02.tif"

gdalwarp  -t_srs EPSG:4326  -te 4 46 12 53 \
    $subc  $DIR/rheinData/subcatchments_rhein.tif \
    -co COMPRESS=LZW -co ZLEVEL=9

#### preparation of the subcatchment IDs. The subcatchment tif file was already available
#### procedure run in R with terra package
#scp sv2:/mnt/shared/sosw/rheinData/subcatchments_rhein.tif ./
#library(terra)
#r = raster("subcatchments_rhein.tif")
#l = unique(r)
#write.table(l, file="~/data/sosw/subc_ids_rhein.txt", row.names = FALSE, col.names = FALSE)
#scp subc_ids_rhein.txt  sv2:/mnt/shared/sosw/rheinData



### STEP 2:  Create the master table for the area of interest

bio=( $(ls /mnt/shared/Environment90m_v.1.0_online/Climate/present) )

#hyd=( $(ls /mnt/shared/Environment90m_v.1.0_online/Hydrography90m) )

hyd=( chancurv chandistdwseg chandistupcel chandistupseg chanelvdwcel chanelvdwseg chanelvupcel chanelvupseg changraddwseg changradupcel changradupseg elev_drop flow flow_accum flowpos gradient length out_dist out_drop outdiffdwbasin outdiffdwscatch outdistdwbasin outdistdwscatch outlet_elev slopdiff slopgrad  source_elev strdiffdwnear strdiffupfarth strdiffupnear strdistdwnear strdistprox strdistupfarth strdistupnear stright )

soil=( soil_ACDWRB soil_AWCtS soil_BDRICM soil_BDRLOG soil_BLDFIE soil_CECSOL soil_CLYPPT soil_CRFVOL soil_HISTPR soil_ORCDRC soil_PHIHOX soil_SLGWRB soil_SLTPPT soil_SNDPPT soil_TEXMHT soil_WWP )

cob=( c10_2020 c20_2020 c30_2020 c40_2020 c50_2020 c60_2020 c70_2020 c80_2020 c90_2020 c100_2020 c110_2020 c120_2020 c130_2020 c140_2020 c150_2020 c160_2020 c170_2020 c180_2020 c190_2020 c200_2020 c210_2020 c220_2020 )

ENV=( ${bio[@]} ${hyd[@]} ${soil[@]} ${cob[@]} )
printf -v joined '%s,' "${ENV[@]}"

#### !!!!!!!!!!!!!!!!!!!!   check if the function works with the new names of the files!!!!!!!!

time  bash /mnt/shared/danube/dev_create_roi_table.sh \
  $(echo "${joined%,}") \
  mean,sd \
  h18v02 \
  /mnt/shared/sosw/env_tl  \
  /mnt/shared/sosw/rheinData/subc_ids_rhein.txt  \
  /mnt/shared/sosw/bio_modelling/rhein/rhein_predictTB.csv  \
  /mnt/shared/sosw/tmp \
  20

### Check resullts
awk -F, '{print NF}'  /mnt/shared/sosw/bio_modelling/rhein/rhein_predictTB.csv | sort | uniq


### STEP 3: merge previous table witht the water discharge data

########   DATA FROM UTRECHT   ################################################
# main directory
export DIS=/mnt/shared/sosw/proc_dan


# #  extrcat discharge last year 2019
gdal_translate -of GTiff -b 10 -a_srs EPSG:4326  \
     -co COMPRESS=DEFLATE -co ZLEVEL=9 \
     $DIS/2010_2019_yearly_discharge.nc $DIS/Utrecht_dis19_tmp.tif
#    vtr 0.0083333333 -0.008333333 \
#    -projwin 8 51 30 42  \


# RHEIN
gdalwarp  -t_srs EPSG:4326  $(pkinfo -i $DIR/rheinData/subcatchments_rhein.tif -te) \
    -tr 0.00833333333333 -0.00833333333333 \
    -tap  $DIS/Utrecht_dis19_tmp.tif \
    $DIR/rheinData/rhein_dis19.tif \
    -co COMPRESS=LZW -co ZLEVEL=9

# create a GRASS GIS Database
mkdir $DIR/grassdata
#sudo grass --text -c $DIR/rheinData/subcatchments_rhein.tif $DIR/grassdata/rhein

sudo grass --text $DIR/grassdata/rhein/PERMANENT

r.in.gdal --o input=/mnt/shared/sosw/rheinData/subcatchments_rhein.tif output=subcat
r.in.gdal --o input=/mnt/shared/sosw/rheinData/rhein_dis19.tif output=disch19

VARNAME=maf
VARNAME=disch19

    # Add header to the output .csv file
    echo "subc_id,${VARNAME}_data_cells,${VARNAME}_nodata_cells,${VARNAME}_min,${VARNAME}_max,${VARNAME}_range,${VARNAME}_mean,${VARNAME}_mean_abs,${VARNAME}_sd,${VARNAME}_var,${VARNAME}_cv,${VARNAME}_sum,${VARNAME}_sum_abs"  > $DIR/rheinData/stats_${VARNAME}_tmp.csv
    # Calculate zonal statistics and append to the output .csv
    r.univar --qq -t --o map=$VARNAME zones=subcat  separator=comma |  awk -F, 'BEGIN{OFS=",";} NR>1 {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}'  \
        >> $DIR/rheinData/stats_${VARNAME}_tmp.csv
##############################################################


# extracat the fields (columns) of interest
    awk -F, 'BEGIN{OFS=",";}  $2 > 0 {print $1, $7, $9}' \
        $DIR/rheinData/stats_${VARNAME}_tmp.csv > $DIR/rheinData/stats_${VARNAME}.csv

rm $DIR/rheinData/stats_${VARNAME}_tmp.csv

# based on subc iD de $DIR/stats_${VARNAME}.csv extrcat those subcatchment from
# $TB and the resultan should be the new predictionTB because for discharge there
# is only data on those subcatchment

# here is the table with all env. lavers
TB=/mnt/shared/sosw/bio_modelling/rhein/rhein_predictTB.csv


awk -F, 'NR > 1 {print $1}' $DIR/rheinData/stats_${VARNAME}.csv > $TMP/ids_subc_new.txt

awk -F, 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $TMP/ids_subc_new.txt $TB \
     > $TMP/subTB.csv

# here check if test has the same number of rows as $DIR/stats_${VARNAME}.csv.
# if not the run below

awk -F, 'NR > 1 {print $1}' $TMP/subTB.csv > $TMP/ids_subc_new2.txt

awk -F, 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $TMP/ids_subc_new2.txt $DIR/rheinData/stats_${VARNAME}.csv \
     > $TMP/sub_stats_${VARNAME}.csv

# and paste tables for final table with new variables
paste -d"," $TMP/subTB.csv \
    <(awk -F, 'BEGIN{OFS=",";} {print $2, $3}' $TMP/sub_stats_disch19.csv) \
    > $DIR/rheinData/predTBrhein.csv

# clean a bit
rm $TMP/sub_stats_${VARNAME}.csv $TMP/ids_subc_new2.txt $TMP/subTB.csv \
    $TMP/ids_subc_new.txt

#########################
#######   CHUNK TO CREATE PREDICTION/PROJECTION TABLE
#######   BASED ON 5/6 KM² FILE FOR RHINE BASIN ONLY

###  open grasss session
###  load subcatchment file
###  load projection(s) file
###  calculate r.univar statistics
###  extract file with only subcatchment IDs where infomration is available
###  close grass
###  subset big prediction (rhein_predictTB.csv) table based on subc IDs and
###      with environmental columns of interest
###  merge tables

export DIS=$DIR/utrecht/biodiversity_data 

for ssp in ssp3 ssp5
do
    for year in 25 86 
    do
    # #  extrcat discharge 2040 ssp3
    gdal_translate -of GTiff -b ${year} -a_srs EPSG:4326  \
         -co COMPRESS=DEFLATE -co ZLEVEL=9 \
         $DIS/${ssp}_rhine_yearly_discharge.nc $DIS/utrecht_dis${year}${ssp}.tif
     done
done


sudo grass --text $DIR/grassdata/rhein/PERMANENT

r.in.gdal --o input=/mnt/shared/sosw/utrecht/biodiversity_data/utrecht_dis25ssp3.tif output=dis25ssp3

VARNAME=dis25ssp3

    # Add header to the output .csv file
    echo "subc_id,${VARNAME}_data_cells,${VARNAME}_nodata_cells,${VARNAME}_min,${VARNAME}_max,${VARNAME}_range,${VARNAME}_mean,${VARNAME}_mean_abs,${VARNAME}_sd,${VARNAME}_var,${VARNAME}_cv,${VARNAME}_sum,${VARNAME}_sum_abs"  > /mnt/shared/sosw/utrecht/stats_${VARNAME}_tmp.csv
    # Calculate zonal statistics and append to the output .csv
    r.univar --qq -t --o map=$VARNAME zones=subcat  separator=comma |  awk -F, 'BEGIN{OFS=",";} NR>1 {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}'  \
        >> /mnt/shared/sosw/utrecht/stats_${VARNAME}_tmp.csv

# extract the fields (columns) of interest and only subcatchment with data
# for the time being extrcat only mean
    awk -F, 'BEGIN{OFS=",";}  $2 > 0 {print $1, $7}' \
        /mnt/shared/sosw/utrecht/stats_${VARNAME}_tmp.csv \
        > /mnt/shared/sosw/utrecht/stats_${VARNAME}.csv

    rm /mnt/shared/sosw/utrecht/stats_${VARNAME}_tmp.csv
#----------   exit GRASS

VARNAME=dis25ssp3
awk -F, 'NR > 1 {print $1}' $DIR/utrecht/stats_${VARNAME}.csv \
    > $TMP/ids_predArea.txt

awk -F, 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $TMP/ids_predArea.txt $DIR/rheinData/predTBrhein.csv \
     > $TMP/rhine_proj_tmp.csv

# if number of rows do not coincide then:
awk -F, 'NR > 1 {print $1}' $TMP/rhine_proj_tmp.csv > $TMP/ids_predArea2.txt

awk -F, 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $TMP/ids_predArea2.txt /mnt/shared/sosw/utrecht/stats_${VARNAME}.csv \
     > $TMP/sub_stats_${VARNAME}.csv

# and paste tables for final table with new variables
paste -d"," $TMP/rhine_proj_tmp.csv \
    <(awk -F, 'BEGIN{OFS=",";} {print $2}' $TMP/sub_stats_${VARNAME}.csv) \
    > $DIR/rheinData/projTBrhein.csv

# clean a bit
rm $TMP/sub_stats_${VARNAME}.csv $TMP/ids_predArea.txt $TMP/rhine_proj_tmp.csv  \
    $TMP/ids_predArea2.txt

#########################
#########################

#############    STEP 4:   create the tables per species

#scp ~/proyectos/sosw/species_list_rhein.txt sv2:/mnt/shared/sosw/rheinData
### file with species list:
splst=/mnt/shared/sosw/rheinData/species_list_rhein.txt

### database with species ocuurence data
occ=/mnt/shared/sosw/rheinData/fish_rhein.csv

for i in {1..31}
do
    spp=$(awk -v row="${i}" 'NR == row' $splst)
    awk -F, -v SP="${spp}"  '$2 ~ SP' $occ > $DIR/rheinData/spptmp_${i}.csv  
done

cat $(find $DIR/rheinData/spptmp*.csv) > $DIR/rheinData/test.csv 

awk -F, '{print $3, $4}' $DIR/rheinData/test.csv  | gdallocationinfo -valonly -wgs84 $DIR/rheinData/subcatchments_rhein.tif | awk '!NF{$0=0}1' > $DIR/rheinData/test2.csv

paste -d"," $DIR/rheinData/test.csv $DIR/rheinData/test2.csv > $DIR/rheinData/test3.csv

echo "Code,spp_name,Lon,Lat,Year,subcid" > $DIR/bio_modelling/rhein/species_db_rhine.csv
awk -F, '$6 > 0' $DIR/rheinData/test3.csv >> $DIR/bio_modelling/rhein/species_db_rhine.csv

rm $DIR/rheinData/spptmp_*.csv $DIR/rheinData/test.csv $DIR/rheinData/test2.csv \
    $DIR/rheinData/test3.csv

#### remove duplicated using helper function

bash $DIR/helpFunc_rmDuplicates.sh \
    $DIR/bio_modelling/rhein/species_db_rhine.csv \
    $DIR/rheinData/subcatchments_rhein.tif \
    $DIR \
    $TMP/spp_r \
    $DIR/bio_modelling/rhein/species_db_rhine_nodup.csv

### make a new list of species to model , that is, species that after removing duplicates still have more than 10 ocurrence points
awk -F, '{print $2}' $DIR/bio_modelling/rhein/species_db_rhine_nodup.csv
| sort | uniq -c | awk '$1 > 10 {print  $2, $3}' > $DIR/bio_modelling/rhein/spp_to_model.txt

###  prepare the tables for each species with ocurrence/absence and environmental
###  variables

# make sure folder exist
mkdir $DIR/bio_modelling/rhein/spp_input

for i in $(seq 2 $(wc -l <  $DIR/bio_modelling/rhein/spp_to_model.txt))
do
    spp=$(awk -v ROW="${i}" 'NR == ROW' $DIR/bio_modelling/rhein/spp_to_model.txt)
    spfnm=$(echo $spp | awk 'BEGIN{OFS="_";}{print $1, $2}')
    echo "species,X,Y" > $TMP/${spfnm}_rhein_tmp.csv
    awk -F, -v SPP="${spp}" 'BEGIN{OFS=",";}  $2 == SPP {print $2, $3, $4}' \
        $DIR/bio_modelling/rhein/species_db_rhine_nodup.csv \
        >> $TMP/${spfnm}_rhein_tmp.csv 

bash $DIR/dev_create_model_table.sh  \
   $TMP/${spfnm}_rhein_tmp.csv \
   $DIR/rheinData/predTBrhein.csv \
   /mnt/shared/sosw/rheinData/subcatchments_rhein.tif  \
   10000  \
   $TMP  \
   /mnt/shared/sosw/bio_modelling/rhein/spp_input/${spfnm}.csv

done


#### for some tables the first row has 0 as subcacthment id. Probably because
#### the point does not overlap with the subcatchments. delete that row

#for arch in $(find $DIR/bio_modelling/rhein/spp_input -name "*.csv")
#do
#
#    SUBID=$(awk -F, 'NR == 2 {print $1}' $arch)
#
#    if [ $SUBID -eq 0 ]
#    then
#        nombre=$(basename $arch .csv)
#        awk -F, 'NR != 2' $arch > $TMP/${nombre}_temp.csv
#        mv $TMP/${nombre}_temp.csv $arch
#    fi
#
#done

#################################
#################################
#######   R CODE TO RUN RANDOM FOREST


library(dplyr)
library(ranger)
library(hydrographr)
library(modEvA)
library(PresenceAbsence)

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

# Predict table path
cat(paste0("Importing table with environmental predictors", "\n"))
pred_tb_c = read.csv("/mnt/shared/sosw/rheinData/predTBrhein.csv")
pred_tb_f = read.csv("/mnt/shared/sosw/rheinData/projTBrhein.csv")

nm = names(pred_tb_c[,c(2,8,14,34,58,84,91,139,155)])
nmf = names(pred_tb_f[,c(2,8,14,34,58,84,91,139,157)])

#predict_table_path <- "/mnt/shared/sosw/bio_modelling/rhein/"
#load(predict_table_path)


# Import presences (coordinates)
#occur <- read.csv("/mnt/shared/danube/models/input/danube_records_final_20-03-2024.csv") |>
#  dplyr::select(subcatchment_id,species)

# Species list
#all_species <- occur |>
#  distinct(species) |>
#  pull()

all_species = read.table("/mnt/shared/sosw/bio_modelling/rhein/spp_to_model.txt")

#for (sp in all_species) {
for (i in 2:dim(all_species)[1])  {

  sp = paste0(all_species[i,1], "_", all_species[i,2])  
  
  #cat(paste0("Start modeling ",sp,': ', which(all_species %in% sp) , 
  #            ' of ', length(all_species)), " species", "\n")
  
  # Record the start time
  start <- Sys.time() 


  # IDs for presences
#  ids_p <- occur |>
#    dplyr::filter(species == sp) |>
#    pull(subcatchment_id)
#  
#  # Generate pseudo-absences. Randomly, use sub-catchment IDs. 10000 pseudo-absences
#  ids_a <-  sample(ids_study_area[!ids_study_area %in% ids_p],
#                   10000)
#  
#  # Bind presences and pseudo-absences
#  pa <- data.frame("subcID" = c(ids_p, ids_a),
#                   "PresAbs" =  c(rep(1, length(ids_p)), rep(0, length(ids_a))))
#    
#  
#  # Make data frame with environmental variables for presences and pseudo-absences
#  env_data_pa <- env_data |>
#    inner_join(pa, by = c("subcID" = "subcID"))
  
env_data_pa = read.csv(paste0("/mnt/shared/sosw/bio_modelling/rhein/spp_input/",
sp, ".csv"))
env_data_pa = env_data_pa[rowSums(is.na(env_data_pa)) == 0,]

  # separate the original data in one subset  for training (80 % of records) and the other for 
  #  validation (20 % of records). 
  a <- SampleMat2(ref = env_data_pa$PresAbs, ratio = 0.8) 
  data_train <- env_data_pa[a$calibration, ] 
  data_test <- env_data_pa[a$evaluation, ]
  
  # number of presence records:
  pres_num <- as.numeric(table(data_train$PresAbs)["1"])
  sample_size <- c("0" = pres_num / nrow(data_train),
                   "1" = pres_num / nrow(data_train))
  
  # convert the response to factor for RF model to return probabilities
  data_train$PresAbs <- as.factor(data_train$PresAbs)
  
  # train model
  cat(paste0("Train model for ",sp, "\n") )
  model <- NULL
  model <- ranger(data_train$PresAbs ~ .,
                  data = data_train[,nm],
#                  data = data_train |> dplyr::select(-subcID),
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
    #pred_test <- predict(model, data = data_test |> dplyr::select(-subcID, -PresAbs))
   # pred_test <- predict(model, data = data_test[,nm])
    
    # Calculate AUC. Use AUC to asses model predictive ability
    #auc <- AUC(obs = data_test$PresAbs,
    #           pred = as.numeric(pred_test$predictions[,2]),
    #           plot = FALSE)$AUC
    
    # Calculate threshold
    #th <- data.frame("ids" = data_test$subcID,
    #                 "obs" = data_test$PresAbs,
    #                  "pred" = as.numeric(pred_test$predictions[,2])) |>
    #  optimal.thresholds(opt.methods = c("Sens=Spec"))
    
    # Calculate threshold dependent metrics
    #th_metrics <- threshMeasures(obs = data_test$PresAbs, 
    #                             pred = as.numeric(pred_test$predictions[,2]),
    #                             simplif = T, 
    #                             thresh = th[,2], 
    #                             standardize = T,
    #                             measures = c("Sensitivity",
    #                                          "Specificity", 
    #                                          "Omission", 
    #                                          "Commission", 
    #                                          "TSS"))
  
  # predict
  cat(paste0("Prediction for ",sp, "\n"))
  #pred <- predict(model, data = pred_table|> dplyr::select(-subcID))
  predc <- predict(model, data = pred_tb_c[,nm])

  # create subset of prediction table and rename columns same as names of present table
subf = pred_tb_f[,nmf]
names(subf) = nm

  predf <- predict(model, data = subf)
#  predf <- predict(model, data = pred_tb_f[,nmf])
  
  # extract predictions and make a data frame with sub-catchment IDs
  predictionc <- data.frame("subcID" = pred_tb_c$subcID,
                           "p" = as.numeric(predc$predictions[,2]))
  predictionf <- data.frame("subcID" = pred_tb_f$subcID,
                           "p" = as.numeric(predf$predictions[,2]))
  
  
  # Path to raster layer with sub-catchments
  #path_subc <- "/mnt/shared/sosw/rheinData/subcatchments_rhein.tif"
  
  # Output path
  #path_map_output <- paste0(
  #                          "/mnt/shared/danube/models/output/rf/",
  #                          sub(" ", "_", sp),
  #                          "_rf.tif")
  
  # Multiply the probability values by 100, to convert them to integers
  prediction4mapc <- predictionc |>
    mutate(p = as.integer(round(predictionc$p, 2) * 100))

  outc = paste0("/mnt/shared/sosw/bio_modelling/rhein/reclassTB/",sp,"_c.txt")  
  write.table(prediction4mapc, outc, row.names=FALSE, col.names=FALSE)
 
  prediction4mapf <- predictionf |>
    mutate(p = as.integer(round(predictionf$p, 2) * 100))
  
  outf = paste0("/mnt/shared/sosw/bio_modelling/rhein/reclassTB/",sp,"_f.txt")  
  write.table(prediction4mapf, outf, row.names=FALSE, col.names=FALSE)
}


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

#########################################################
##############################################################################
##############################################################################
##############################################################################



##############
#######   reclassification of predictions

sudo grass --text /mnt/shared/sosw/grassdata/rhein/PERMANENT

dir=/mnt/shared/sosw
rcl=$dir/bio_modelling/rhein/reclassTB
out=$dir/bio_modelling/rhein/randomForest_output

#r.in.gdal --o input=/mnt/shared/sosw/mekongData/subcatchments_mekong.tif output=subcat

paths=( $(ls $rcl/*.txt) )
#tiff=${paths[0]}

for tiff in ${paths[@]:2}
do
spp=$(basename $tiff .txt)
outn=$(echo ${spp}_prob.tif)

awk '{print $1, "=", $2}' $rcl/${spp}.txt > $rcl/${spp}_rcl.txt

r.reclass --o input=subcat output=$outn rules=$rcl/${spp}_rcl.txt

r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=Byte \
    format=GTiff nodata=255 input=$outn output=$out/${outn}

done


for sp in $(ls $out | awk -F_ 'BEGIN{OFS="_";}{print $1, $2}' | sort | uniq)
do
    r.mapcalc --o "${sp}_diff = ${sp}_f_prob.tif - ${sp}_c_prob.tif"

    r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=Int16 \
        format=GTiff nodata=-9999 input=${sp}_diff output=${sp}_diff.tif

done



