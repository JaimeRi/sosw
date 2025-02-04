#! /bin/bash

### STEP 1: download to a folder the environmental variables of interest

TILE=(h18v02 h18v04 h20v02 h20v04)

for TL in ${TILE[@]}
do

	for file in $(find /mnt/shared/Environment90m_v.1.0_online/LandCover/ -name "*_2020_${TL}.zip")
	do 
	sudo unzip -j $file -d /mnt/shared/danube/env_tl
	done

	for file in $(find /mnt/shared/Environment90m_v.1.0_online/Climate/present -name "*_${TL}.zip")
	do 
	sudo unzip -j $file -d /mnt/shared/danube/env_tl
	done
	
	for file in $(find /mnt/shared/Environment90m_v.1.0_online/Hydrography90m -name "*_${TL}.zip")
	do 
	sudo unzip -j $file -d /mnt/shared/danube/env_tl
	done
	
	for file in $(find /mnt/shared/Environment90m_v.1.0_online/Soil -name "*_${TL}.zip")
	do 
	sudo unzip -j $file -d /mnt/shared/danube/env_tl
	done
	
done

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
  h18v02,h18v04,h20v02,h20v04 \
  /mnt/shared/danube/env_tl  \
  /mnt/shared/danube/subc_IDs.txt  \
  /mnt/shared/danube/out/danube_predictTB.csv  \
  /mnt/shared/tmp/danube \
  20

### Check resullts
awk -F, '{print NF}' /mnt/shared/danube/out/danube_predictTB.csv | sort | uniq




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

#############    STEP 4:   create the tables per species


### file with species list:
splst=/mnt/shared/danube/species_list.csv
awk -F, 'NR > 1 {print $1}' $splst > $DIR/spp_list.txt

### database with species ocuurence data
occ=/mnt/shared/danube/fish_danube.csv

spp="Sander lucioperca"

awk -F, -v SPP="${spp}" 'BEGIN{OFS=",";}  $2 == SPP {print $2, $3, $4}' $occ > $DIR/spp1.csv 


 time bash dev_create_model_table.sh  \
   $DIR/spp1.csv \
   $DIR/predTBdanube.csv \
   /mnt/shared/sosw/tmp/danube_subcatchments.tif  \
   10000  \
   /mnt/shared/danube/tmp  \
   /mnt/shared/sosw/bio_modelling/danube/spp_input/spp_table.csv

for i in $(seq 1 $(wc -l <  $DIR/spp_list.txt))
do
    spp=$(awk -v ROW="${i}" 'NR == ROW' $DIR/spp_list.txt)
    spfnm=$(echo $spp | awk 'BEGIN{OFS="_";}{print $1, $2}')
    
    echo "species,X,Y" > /mnt/shared/danube/tmp/${spfnm}_tmp.csv
    awk -F, -v SPP="${spp}" 'BEGIN{OFS=",";}  $2 == SPP {print $2, $3, $4}' $occ \
    >> /mnt/shared/danube/tmp/${spfnm}_tmp.csv 

bash dev_create_model_table.sh  \
   /mnt/shared/danube/tmp/${spfnm}_tmp.csv \
   $DIR/predTBdanube.csv \
   /mnt/shared/sosw/tmp/danube_subcatchments.tif  \
   10000  \
   /mnt/shared/danube/tmp  \
   /mnt/shared/sosw/bio_modelling/danube/spp_input/${spfnm}.csv

done




mkdir $DIR/cog_layers

# runs also on server3
OUT=$DIR/cog_layers

for basin in jucar danube rhein mekong
do
    mkdir $OUT/${basin}

    for FILE in $DIR/bio_modelling/${basin}/randomForest_output/*.tif ; do
        filename="${FILE##*/}"

        gdaladdo -r average $FILE 2 4 6 8 10 12 14 16 18 20 22 24 26 28 32 64 128

        ### create cloud-optimised tif
        gdal_translate $FILE   $DIR/cog_layers/${basin}/COG_${filename}  \
            -co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=LZW 
    done
done
