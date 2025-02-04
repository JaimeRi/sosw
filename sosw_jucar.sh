#! /bin/bash

export DIR=/mnt/shared/sosw
export TMP=/mnt/shared/sosw/tmp


### STEP 1: download to a folder the environmental variables of interest

##   TILE=(h16v04 h18v04) tiles for Jucar but h18v04 has been downloaded previously
TILE=(h16v04)

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

### STEP 2:  Create the master table for the area of interest

bio=( $(ls /mnt/shared/Environment90m_v.1.0_online/Climate/present) )

#hyd=( $(ls /mnt/shared/Environment90m_v.1.0_online/Hydrography90m) )

hyd=( chancurv chandistdwseg chandistupcel chandistupseg chanelvdwcel chanelvdwseg chanelvupcel chanelvupseg changraddwseg changradupcel changradupseg elev_drop flow flow_accum flowpos gradient length out_dist out_drop outdiffdwbasin outdiffdwscatch outdistdwbasin outdistdwscatch outlet_elev slopdiff slopgrad  source_elev strdiffdwnear strdiffupfarth strdiffupnear strdistdwnear strdistprox strdistupfarth strdistupnear stright )

soil=( soil_ACDWRB soil_AWCtS soil_BDRICM soil_BDRLOG soil_BLDFIE soil_CECSOL soil_CLYPPT soil_CRFVOL soil_HISTPR soil_ORCDRC soil_PHIHOX soil_SLGWRB soil_SLTPPT soil_SNDPPT soil_TEXMHT soil_WWP )

cob=( c10_2020 c20_2020 c30_2020 c40_2020 c50_2020 c60_2020 c70_2020 c80_2020 c90_2020 c100_2020 c110_2020 c120_2020 c130_2020 c140_2020 c150_2020 c160_2020 c170_2020 c180_2020 c190_2020 c200_2020 c210_2020 c220_2020 )

ENV=( ${bio[@]} ${hyd[@]} ${soil[@]} ${cob[@]} )
printf -v joined '%s,' "${ENV[@]}"

#### !!!!!!!!!!!!!!!!!!!!   check if the function works with the new names of the files!!!!!!!!


#### preparation of the subcatchment IDs. The subcatchment tif file was already available
#### procedure run in R with terra package
#library(terra)
#r = raster("~/data/sosw/env_jucar/final_layers/sub_catchment.tif")
#l = unique(r)
#write.table(l, file="~/data/sosw/env_jucar/final_layers/subc_ids.txt", row.names = FALSE, col.names = FALSE)
#scp subc_ids.txt  sv2:/mnt/shared/sosw/jucarData

time  bash /mnt/shared/danube/dev_create_roi_table.sh \
  $(echo "${joined%,}") \
  mean,sd \
  h16v04,h18v04 \
  /mnt/shared/sosw/env_tl  \
  /mnt/shared/sosw/jucarData/subc_ids.txt  \
  /mnt/shared/sosw/bio_modelling/jucar/jucar_predictTB.csv  \
  /mnt/shared/sosw/tmp \
  20

### Check resullts
awk -F, '{print NF}'  /mnt/shared/sosw/bio_modelling/jucar/jucar_predictTB.csv | sort | uniq


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


# Júcar  POLYGON((-9.86 46.4, 6.29 46.4, 6.29 34.74, -9.86 34.74, -9.86 46.4))
gdalwarp  -t_srs EPSG:4326  $(pkinfo -i $DIR/jucarData/sub_catchment.tif -te) \
    -tr 0.00833333333333 -0.00833333333333 \
    -tap  $DIS/Utrecht_dis19_tmp.tif \
    $DIR/jucarData/jucar_dis19.tif \
    -co COMPRESS=LZW -co ZLEVEL=9

# create a GRASS GIS Database
mkdir $DIR/grassdata
#sudo grass --text -c $DIR/jucarData/sub_catchment.tif $DIR/grassdata/jucar 

sudo grass --text $DIR/grassdata/jucar/PERMANENT

r.in.gdal --o input=/mnt/shared/sosw/jucarData/sub_catchment.tif output=subcat
r.in.gdal --o input=/mnt/shared/sosw/jucarData/jucar_dis19.tif output=disch19

VARNAME=maf
VARNAME=disch19

    # Add header to the output .csv file
    echo "subc_id,${VARNAME}_data_cells,${VARNAME}_nodata_cells,${VARNAME}_min,${VARNAME}_max,${VARNAME}_range,${VARNAME}_mean,${VARNAME}_mean_abs,${VARNAME}_sd,${VARNAME}_var,${VARNAME}_cv,${VARNAME}_sum,${VARNAME}_sum_abs"  > $DIR/jucarData/stats_${VARNAME}_tmp.csv
    # Calculate zonal statistics and append to the output .csv
    r.univar --qq -t --o map=$VARNAME zones=subcat  separator=comma |  awk -F, 'BEGIN{OFS=",";} NR>1 {print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}'  \
        >> $DIR/jucarData/stats_${VARNAME}_tmp.csv
##############################################################

# extrcat the fields (columns) of interest
    awk -F, 'BEGIN{OFS=",";}  $2 > 0 {print $1, $7, $9}' \
        $DIR/jucarData/stats_${VARNAME}_tmp.csv > $DIR/jucarData/stats_${VARNAME}.csv

rm $DIR/jucarData/stats_${VARNAME}_tmp.csv

# based on subc iD de $DIR/stats_${VARNAME}.csv extrcat those subcatchment from
# $TB and the resultan should be the new predictionTB because for discharge there
# is only data on those subcatchment

# here is the table with all env. lavers
TB=/mnt/shared/sosw/bio_modelling/jucar/jucar_predictTB.csv


awk -F, 'NR > 1 {print $1}' $DIR/jucarData/stats_${VARNAME}.csv > $TMP/ids_subc_new.txt

awk -F, 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $TMP/ids_subc_new.txt $TB \
     > $TMP/subTB.csv

# here check if test has the same number of rows as $DIR/stats_${VARNAME}.csv.
# if not the run below

awk -F, 'NR > 1 {print $1}' $TMP/subTB.csv > $TMP/ids_subc_new2.txt

awk -F, 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $TMP/ids_subc_new2.txt $DIR/jucarData/stats_${VARNAME}.csv \
     > $TMP/sub_stats_${VARNAME}.csv

# and paste tables for final table with new variables
paste -d"," $TMP/subTB.csv \
    <(awk -F, 'BEGIN{OFS=",";} {print $2, $3}' $TMP/sub_stats_disch19.csv) \
    > $DIR/jucarData/predTBjucar.csv

# clean a bit
rm $TMP/sub_stats_${VARNAME}.csv $TMP/ids_subc_new2.txt $TMP/subTB.csv \
    $TMP/ids_subc_new.txt


#############    STEP 4:   create the tables per species


### file with species list:
#splst=/mnt/shared/danube/species_list.csv
#awk -F, 'NR > 1 {print $1}' $splst > $DIR/spp_list.txt
printf  %"s %s\n"  Anguilla anguilla Salaria ﬂuviatilis Cobitis paludica \
    Achondrostoma arcasii Barbus haasi Luciobarbus guiraonis Parachondrostoma arrigonis \
    Parachondrostoma turiense Squalius valentinus Salmo trutta Valencia hispanica \
    > $DIR/jucarData/species_jucar_list.txt

### database with species ocuurence data
#occ=/mnt/shared/sosw/bio_tb/fish_jucar.csv
occ=/mnt/shared/sosw/jucarData/fish_jucar_mart_fn.csv

awk -F, 'NR > 1 {print $2}' $occ | sort | uniq > $DIR/jucarData/species_jucar_list.txt 

#spp="Sander lucioperca"

#awk -F, -v SPP="${spp}" 'BEGIN{OFS=",";}  $2 == SPP {print $2, $3, $4}' $occ > $DIR/spp1.csv 
#
#
# time bash dev_create_model_table.sh  \
#   $DIR/spp1.csv \
#   $DIR/predTBdanube.csv \
#   /mnt/shared/sosw/tmp/danube_subcatchments.tif  \
#   10000  \
#   /mnt/shared/danube/tmp  \
#   /mnt/shared/sosw/bio_modelling/danube/spp_input/spp_table.csv

for i in $(seq 1 $(wc -l <  $DIR/jucarData/species_jucar_list.txt))
do
    spp=$(awk -v ROW="${i}" 'NR == ROW' $DIR/jucarData/species_jucar_list.txt)
    spfnm=$(echo $spp | awk 'BEGIN{OFS="_";}{print $1, $2}')
    echo "species,X,Y" > $TMP/${spfnm}_jucar_tmp.csv
    awk -F, -v SPP="${spp}" 'BEGIN{OFS=",";}  $2 == SPP {print $2, $3, $4}' $occ \
    >> $TMP/${spfnm}_jucar_tmp.csv 

bash /mnt/shared/danube/dev_create_model_table.sh  \
   $TMP/${spfnm}_jucar_tmp.csv \
   $DIR/jucarData/predTBjucar.csv \
   /mnt/shared/sosw/jucarData/sub_catchment.tif  \
   10000  \
   $TMP  \
   /mnt/shared/sosw/bio_modelling/jucar/spp_input/${spfnm}.csv

done


#### for some tables the first row has 0 as subcacthment id. Probably because
#### the point does not overlap with the subcatchments. delete that row

for arch in $(find $DIR/bio_modelling/jucar/spp_input -name "*.csv")
do

    SUBID=$(awk -F, 'NR == 2 {print $1}' $arch)

    if [ $SUBID -eq 0 ]
    then
        nombre=$(basename $arch .csv)
        awk -F, 'NR != 2' $arch > $TMP/${nombre}_temp.csv
        mv $TMP/${nombre}_temp.csv $arch
    fi

done


# create richness map


sudo grass --text $DIR/grassdata/jucar/PERMANENT

files=( $(find randomForest_output -name "*.tif") )
for indx in {0..34}
do
    r.in.gdal --o input=${files[${indx}]} output=spp_${indx}
done

r.mapcalc "suma = spp_0 + spp_1 + spp_2 + spp_3 + spp_4 + spp_5 + spp_6 + spp_7 + spp_8 + spp_9 + spp_10 + spp_11 + spp_12 + spp_13 + spp_14 + spp_15 + spp_16 + spp_17 + spp_18 + spp_19 + spp_20 + spp_21 + spp_22 + spp_23 + spp_24 + spp_25 + spp_26 + spp_27 + spp_28 + spp_29 + spp_30 + spp_31 + spp_32 + spp_33 + spp_34"

r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=Int32 \
    format=GTiff nodata=-9999 input=suma output=richness.tif
