#! /bin/bash

### script to create the sdm workflow for the Danube using the danube mask and
### discharge variables provided by IIASA

# main directory
export DIR=/mnt/shared/sosw/proc_dan


#######   DATA FROM IIASA #####################################################

# danube mask procided by IIASA
export MASK=$DIR/basin.tif

# variable of interest
# MAF_once: Mean of discharge for all days : 1/1/1990-31/12/2022
# Q90_once: 10% of the lowest discharge for all days : 1/1/1990-31/12/2022
# EF_VMF_12month.nc: Environmental flow requirement with Variable Monthly Flow 
# discharge_yearly.nc
export VAR=$DIR/MAF_once.nc

# this data need to be crop to the following extnetion:
gdalwarp  -t_srs EPSG:4326  -te 8 42 30 51 \
    $VAR \
    $DIR/danube_MAF_oce_iasa.tif \
    -co COMPRESS=LZW -co ZLEVEL=9

#############################################################
########   DATA FROM UTRECHT   ################################################


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
TB=/mnt/shared/env_test/projectionTB.csv

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


