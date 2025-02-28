#! /bin/bash

### STEP 1: download to a folder the environmental variables of interest
ENV=/mnt/shared/Environment90m_v.1.0_online/

TILE=(h18v02 h18v04 h20v02 h20v04)

for TL in ${TILE[@]}
do
    for var in c60 c80 
	do 
    file=$(find $ENV/esa_cci_landcover_v2_1_1/${var} -name "${var}_2020_${TL}.zip")
	sudo unzip -j $file -d /mnt/shared/sosw/danube/env
	done

    for var in bio01 bio12 bio04 bio15 
	do 
    file=$(find $ENV/chelsa_bioclim_v2_1/1981-2010_observed/${var} -name "${var}_*_${TL}.zip")    
	sudo unzip -j $file -d /mnt/shared/sosw/danube/env
	done
    
    for var in bio1 bio12 bio4 bio15
    do 
    file=$(find $ENV/chelsa_bioclim_v2_1/2071-2100_projected/${var} -name "${var}_*_ukesm1-0-ll_ssp585_*_${TL}.zip")    
	sudo unzip -j $file -d /mnt/shared/sosw/danube/env
	done
	
    for var in flow length slope_grad_dw_cel
	do 
    file=$(find $ENV/hydrography90m_v1_0/${var} -name "${var}_${TL}.zip")    
	sudo unzip -j $file -d /mnt/shared/sosw/danube/env
	done
	
done


### STEP 2:  Create the master table for the area of interest


# present variables
ENV=( c60_2020 c80_2020 bio01_1981-2010_observed bio12_1981-2010_observed bio04_1981-2010_observed bio15_1981-2010_observed flow length slope_grad_dw_cel)
printf -v joined '%s,' "${ENV[@]}"

time  bash /mnt/shared/sosw/dev_create_roi_table.sh \
  $(echo "${joined%,}") \
  mean \
  h18v02,h18v04,h20v02,h20v04 \
  /mnt/shared/sosw/danube/env  \
  /mnt/marquez/marquez/vignette/data/subc_IDs.txt  \
  /mnt/shared/sosw/danube/danube_predictTB.csv  \
  /mnt/shared/tmp/danube \
  25

# future variables
FUT=(  bio12_2071-2100_ukesm1-0-ll_ssp585_V.2.1 bio15_2071-2100_ukesm1-0-ll_ssp585_V.2.1 bio1_2071-2100_ukesm1-0-ll_ssp585_V.2.1 bio4_2071-2100_ukesm1-0-ll_ssp585_V.2.1 )
printf -v joined '%s,' "${FUT[@]}"


time  bash /mnt/shared/sosw/dev_create_roi_table.sh \
  $(echo "${joined%,}") \
  mean \
  h18v02,h18v04,h20v02,h20v04 \
  /mnt/shared/sosw/danube/env  \
  /mnt/marquez/marquez/vignette/data/subc_IDs.txt  \
  /mnt/shared/sosw/danube/danube_predictFT.csv  \
  /mnt/shared/tmp/danube \
  16

### Check resullts
#awk -F, '{print NF}' /mnt/shared/sosw/danube/danube_predictTB.csv | sort | uniq -c

####  change the name of the bio climate variables for the present and future

sed 's/bio01_1981-2010_observed_mean/bio1/g' $DIR/sdm/danube_predictTB.csv > $TMP/predTB_tmp2.csv
sed 's/bio04_1981-2010_observed_mean/bio4/g' $TMP/predTB_tmp2.csv > $TMP/predTB_tmp3.csv
sed 's/bio12_1981-2010_observed_mean/bio12/g' $TMP/predTB_tmp3.csv > $TMP/predTB_tmp4.csv
sed 's/bio15_1981-2010_observed_mean/bio15/g' $TMP/predTB_tmp4.csv > $DIR/sdm/danube_predictTB.csv

rm $TMP/predTB_tmp*.csv 

sed 's/bio1_2071-2100_ukesm1-0-ll_ssp585_V.2.1_mean/bio1/g' $DIR/sdm/danube_predictFT.csv > $TMP/predTB_tmp2.csv
sed 's/bio4_2071-2100_ukesm1-0-ll_ssp585_V.2.1_mean/bio4/g' $TMP/predTB_tmp2.csv > $TMP/predTB_tmp3.csv
sed 's/bio12_2071-2100_ukesm1-0-ll_ssp585_V.2.1_mean/bio12/g' $TMP/predTB_tmp3.csv > $TMP/predTB_tmp4.csv
sed 's/bio15_2071-2100_ukesm1-0-ll_ssp585_V.2.1_mean/bio15/g' $TMP/predTB_tmp4.csv > $DIR/sdm/danube_predictFT.csv

rm $TMP/predTB_tmp*.csv 

# prepare the predction table for the future by adding the missing columns

paste -d"," \
    $DIR/sdm/danube_predictFT.csv \
    <(awk -F, 'BEGIN{OFS=","} {print $6, $7, $8, $9, $10}' $DIR/sdm/danube_predictTB.csv) \
    > $DIR/sdm/danube_predictFT2.csv

mv $DIR/sdm/danube_predictFT2.csv $DIR/sdm/danube_predictFT.csv 



## STEP 3: merge previous table witht the water discharge data

########   DATA FROM UTRECHT   ################################################
# main directory
export DIR=/mnt/shared/sosw/proc_dan


# #  extrcat discharge last year 2019
gdal_translate -of GTiff -b 10 -a_srs EPSG:4326  \
     -co COMPRESS=DEFLATE -co ZLEVEL=9 \
     $DIR/2010_2019_yearly_discharge.nc $DIR/Utrecht_dis19_tmp.tif
#    vtr 0.0083333333 -0.008333333 \
#    -projwin 8 51 30 42  \


gdalwarp  -t_srs EPSG:4326  -te 8 42 30 51 \
    -tr 0.00833333333333 -0.00833333333333 \
    -tap  $DIR/Utrecht_dis19_tmp.tif \
    $DIR/Utrecht_dis19.tif \
    -co COMPRESS=LZW -co ZLEVEL=9



# create a GRASS GIS Database
mkdir $DIR/grassdata
#sudo grass --text -c /mnt/shared/sosw/tmp/basins_sub.tif $DIR/grassdata/danube 

sudo grass --text $DIR/grassdata/danube/PERMANENT

r.in.gdal --o input=/mnt/shared/sosw/tmp/danube_subcatchments.tif output=subcat
r.in.gdal --o input=/mnt/shared/sosw/proc_dan/danube_MAF_oce_iasa.tif output=maf
r.in.gdal --o input=/mnt/shared/sosw/proc_dan/Utrecht_dis19.tif output=disch19

VARNAME=maf
VARNAME=disch19

    # Add header to the output .csv file
    echo "subc_id,${VARNAME}_data_cells,${VARNAME}_nodata_cells,${VARNAME}_min,${VARNAME}_max,${VARNAME}_range,${VARNAME}_mean,${VARNAME}_mean_abs,${VARNAME}_sd,${VARNAME}_var,${VARNAME}_cv,${VARNAME}_sum,${VARNAME}_sum_abs"\
        > $DIR/stats_${VARNAME}_tmp.csv
    # Calculate zonal statistics and append to the output .csv
    r.univar --qq -t --o map=$VARNAME zones=subcat  separator=comma |  awk -F, 'BEGIN{OFS=",";} NR>1 {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}'  \
        >> $DIR/stats_${VARNAME}_tmp.csv
##############################################################

# extrcat the fields (columns) of interest
    awk -F, 'BEGIN{OFS=",";}  $2 > 0 {print $1, $7, $9}' \
        $DIR/stats_${VARNAME}_tmp.csv > $DIR/stats_${VARNAME}.csv

rm $DIR/stats_${VARNAME}_tmp.csv

# based on subc iD de $DIR/stats_${VARNAME}.csv extrcat those subcatchment from
# $TB and the resultan should be the new predictionTB because for discharge there
# is only data on those subcatchment

# here is the table with all env. lavers
#TB=/mnt/shared/danube/out/danube_predictTB.csv
TB=/mnt/shared/danube/out/danube_predictTB.csv

# location of all env. layers / all tiles
# /mnt/shared/danube/env_tiles
# all land cover from 110 and bigger are missing
# current table is here /mnt/shared/danube/out/danube_predictTB.csv

awk -F, 'NR > 1 {print $1}' $DIR/stats_${VARNAME}.csv > $DIR/ids_subc_new.txt

awk -F, 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $DIR/ids_subc_new.txt $TB \
     > $DIR/subTB.csv

# here check if test has the same number of rows as $DIR/stats_${VARNAME}.csv.
# if not the run below

awk -F, 'NR > 1 {print $1}' $DIR/subTB.csv > $DIR/ids_subc_new2.txt

awk -F, 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $DIR/ids_subc_new2.txt $DIR/stats_${VARNAME}.csv \
     > $DIR/sub_stats_${VARNAME}.csv

# and paste tables for final table with new variables
paste -d"," $DIR/subTB.csv \
    <(awk -F, 'BEGIN{OFS=",";} {print $2, $3}' $DIR/sub_stats_disch19.csv) \
    > $DIR/predTBdanube.csv

################################################################################
#############    STEP 4:   create the tables per species

export DIR=/mnt/shared/sosw/danube

### file with species list:
splst=$DIR/fish/species_list.csv
awk -F, 'NR > 1 {print $1}' $splst > $DIR/fish/spp_list.txt

### database with species ocuurence data
#occ1=$DIR/fish/fish_danube.csv
sed 's/\"//g' $DIR/fish/danube_records_final_04-11-2024.csv > $DIR/fish/danube4all.csv
occ=$DIR/fish/danube4all.csv


for i in $(seq 1 $(wc -l <  $DIR/fish/spp_list.txt))
do
    spp=$(awk -v ROW="${i}" 'NR == ROW' $DIR/fish/spp_list.txt)
    spfnm=$(echo $spp | awk 'BEGIN{OFS="_";}{print $1, $2}')

    cat \
        <(echo "species,lon,lat") \
        <(awk -F, -v SPP="${spp}" 'BEGIN{OFS=",";}  $1 == SPP && $2 > 0 {print $1, $3, $4}' $occ) \
        > $DIR/fish/${spfnm}_tmp.csv
    
bash $DIR/../dev_create_model_table.sh  \
   $DIR/fish/${spfnm}_tmp.csv \
   $DIR/sdm/danube_predictTB.csv \
   /mnt/marquez/marquez/vignette/data/subcatchments.tif  \
   10000  \
   /mnt/shared/tmp/danube  \
   $TMP/${spfnm}_tmp1.csv

#sed 's/bio01_1981-2010_observed_mean/bio1/g' $TMP/${spfnm}_tmp1.csv > $TMP/${spfnm}_tmp2.csv
#sed 's/bio04_1981-2010_observed_mean/bio4/g' $TMP/${spfnm}_tmp2.csv > $TMP/${spfnm}_tmp3.csv
#sed 's/bio12_1981-2010_observed_mean/bio12/g' $TMP/${spfnm}_tmp3.csv > $TMP/${spfnm}_tmp4.csv
#sed 's/bio15_1981-2010_observed_mean/bio15/g' $TMP/${spfnm}_tmp4.csv > $DIR/sdm/spp_input/${spfnm}.csv

rm $TMP/${spfnm}_*.csv
rm $DIR/fish/${spfnm}_tmp.csv

done


#mkdir $DIR/cog_layers
#
## runs also on server3
#OUT=$DIR/cog_layers
#
#for basin in jucar danube rhein mekong
#do
#    mkdir $OUT/${basin}
#
#    for FILE in $DIR/bio_modelling/${basin}/randomForest_output/*.tif ; do
#        filename="${FILE##*/}"
#
#        gdaladdo -r average $FILE 2 4 6 8 10 12 14 16 18 20 22 24 26 28 32 64 128
#
#        ### create cloud-optimised tif
#        gdal_translate $FILE   $DIR/cog_layers/${basin}/COG_${filename}  \
#            -co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=LZW 
#    done
#done




#################################
#################################
#######   R CODE TO RUN RANDOM FOREST


library(dplyr)
library(ranger)
library(hydrographr)
library(modEvA)
library(PresenceAbsence)


folder_path <- "/mnt/shared/sosw/danube/sdm/randomForest_out/"

# Predict table path
pred_tb_c = read.csv("/mnt/shared/sosw/danube/sdm/danube_predictTB.csv")
pred_tb_f = read.csv("/mnt/shared/sosw/danube/sdm/danube_predictFT.csv")

nmc = names(pred_tb_c[,c(2:10)])
nmf = names(pred_tb_f[,c(2:10)])


all_species = read.table("/mnt/shared/sosw/danube/fish/species_to_model.txt")


#for (i in 2:dim(all_species)[1])  {
#for (i in c(2:9,24,25)){
for (i in c(10:23,26:39)){

    sp = strsplit(all_species[i,1], split=c("\\."))[[1]][1]  
  # Record the start time
  start <- Sys.time() 

    env_data_pa = read.csv(paste0("/mnt/shared/sosw/danube/sdm/spp_input/",
sp, ".csv"))
    env_data_pa = env_data_pa[rowSums(is.na(env_data_pa)) == 0,]

    data_train = env_data_pa

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
                  data = data_train[,nmc],
                  num.trees = 1000,
                  mtry= 6,
                  replace = T,
                  sample.fraction = sample_size,
                  oob.error = T,
                  keep.inbag = T,
                  num.threads = 4,
                  importance ='impurity',
                  probability = T)
  
  # predict
  cat(paste0("Prediction for ",sp, "\n"))
  predc <- predict(model, data = pred_tb_c[,nmc])
  predf <- predict(model, data = pred_tb_f[,nmf])
  
  # extract predictions and make a data frame with sub-catchment IDs
  predictionc <- data.frame("subcID" = pred_tb_c$subcID,
                           "p" = as.numeric(predc$predictions[,2]))
  predictionf <- data.frame("subcID" = pred_tb_f$subcID,
                           "p" = as.numeric(predf$predictions[,2]))
  
  
  # Multiply the probability values by 100, to convert them to integers
  prediction4mapc <- predictionc |>
    mutate(p = as.integer(round(predictionc$p, 2) * 100))

  outc = paste0("/mnt/shared/sosw/danube/sdm/reclassTB/",sp,"_c.txt")  
  write.table(prediction4mapc, outc, row.names=FALSE, col.names=FALSE)
 
  prediction4mapf <- predictionf |>
    mutate(p = as.integer(round(predictionf$p, 2) * 100))
  
  outf = paste0("/mnt/shared/sosw/danube/sdm/reclassTB/",sp,"_f.txt")  
  write.table(prediction4mapf, outf, row.names=FALSE, col.names=FALSE)
}

################################################################################


##############
#######   reclassification of predictions

sudo grass --text /mnt/shared/sosw/grassdata/danube/PERMANENT

dir=/mnt/shared/sosw/danube
rcl=$dir/sdm/reclassTB
out=$dir/sdm/randomForest_out

#r.in.gdal --o input=/mnt/shared/sosw/danube/gis/msk_danube.tif output=maskt

paths=( $(ls $rcl/*_[cf].txt) )
#tiff=${paths[0]}

for tiff in ${paths[@]}
do

    spp=$(basename $tiff .txt)
    [[ -f $out/${spp}_prob.tif ]] && continue
    outn=$(echo ${spp}_prob.tif)

    awk '{print $1, "=", $2}' $rcl/${spp}.txt > $rcl/${spp}_rcl.txt

    r.reclass --o input=subcat output=$outn rules=$rcl/${spp}_rcl.txt

    r.mapcalc --o "${outn}_msk = if(maskt == 1, $outn, null())"

    r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=Byte \
        format=GTiff nodata=255 input=${outn}_msk output=$out/${outn}

done


for sp in $(ls $out | awk -F_ 'BEGIN{OFS="_";}{print $1, $2}' | sort | uniq)
do

    [[ -f $out/${sp}_diff.tif ]] && continue
    r.mapcalc --o "${sp}_diff = ${sp}_f_prob.tif_msk - ${sp}_c_prob.tif_msk"

    r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=Int16 \
        format=GTiff nodata=-9999 input=${sp}_diff output=$out/${sp}_diff.tif
done

## sum up all the files with the difference between future and present
test=( $(g.list rast | grep "_diff") )
printf -v testn '%s + ' "${test[@]}" 
r.mapcalc --o "all_diff = $(echo "${testn% + }")" 

    r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=Int16 \
        format=GTiff nodata=-9999 input=all_diff output=$out/allspp_diff.tif


