#! /bin/bash

export DIR=/mnt/shared/sosw
export TMP=/mnt/shared/sosw/tmp


### STEP 1: download to a folder the environmental variables of interest

TILE=(h28v06)

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
subc="/mnt/shared/hydrography90m_v.1.0_online/hydrography90m_v.1.0/r.watershed/sub_catchment_tiles20d/sub_catchment_h28v06.tif"

gdalwarp  -t_srs EPSG:4326  -te 101 8 110 24 \
    $subc  $DIR/mekongData/subcatchments_mekong.tif \
    -co COMPRESS=LZW -co ZLEVEL=9

#### preparation of the subcatchment IDs. The subcatchment tif file was already available
#### procedure run in R with terra package
#scp sv2:/mnt/shared/sosw/mekongData/subcatchments_mekong.tif ./
#library(terra)
#r = raster("subcatchments_mekong.tif")
#l = unique(r)
#write.table(l, file="~/data/sosw/subc_ids_mekong.txt", row.names = FALSE, col.names = FALSE)
#scp subc_ids_mekong.txt  sv2:/mnt/shared/sosw/mekongData



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
  h28v06 \
  /mnt/shared/sosw/env_tl  \
  /mnt/shared/sosw/mekongData/subc_ids_mekong.txt  \
  /mnt/shared/sosw/bio_modelling/mekong/mekong_predictTB.csv  \
  /mnt/shared/sosw/tmp \
  20

### Check resullts
awk -F, '{print NF}'  /mnt/shared/sosw/bio_modelling/mekong/mekong_predictTB.csv | sort | uniq

#######   SPECIES DATA  

#############    STEP 4:   create the tables per species

#scp ~/proyectos/sosw/species_list_mekong.txt sv2:/mnt/shared/sosw/mekongData
### file with species list:
splst=/mnt/shared/sosw/mekongData/species_list_mekong.txt

### database with species ocuurence data
occ=/mnt/shared/sosw/bio_tb/fish_mekong_gbif.csv

for i in {1..34}
do
    spp=$(awk -v row="${i}" 'NR == row' $splst)
    awk -F, -v SP="${spp}"  '$6 ~ SP' $occ > $DIR/mekongData/spptmp_${i}.csv  
done

cat $(find $DIR/mekongData/spptmp*.csv) > $DIR/mekongData/test.csv 

awk -F, '{print $8, $7}' $DIR/mekongData/test.csv  | gdallocationinfo -valonly -wgs84 $DIR/mekongData/subcatchments_mekong.tif | awk '!NF{$0=0}1' > $DIR/mekongData/test2.csv

paste -d"," $DIR/mekongData/test.csv $DIR/mekongData/test2.csv > $DIR/mekongData/test3.csv

echo "Code,spp_name,Lon,Lat,Year,subcid" > $DIR/bio_modelling/mekong/species_db_mekong.csv
awk -F, 'BEGIN{OFS=",";} $10 > 0 {print $1, $6, $8, $7, $9,$10}' $DIR/mekongData/test3.csv >> $DIR/bio_modelling/mekong/species_db_mekong.csv

rm $DIR/mekongData/spptmp_*.csv $DIR/mekongData/test.csv $DIR/mekongData/test2.csv \
    $DIR/mekongData/test3.csv

####### remove duplicates
export table="$DIR/bio_modelling/mekong/species_db_mekong.csv"

rmDuplic(){

export spp=$1
export nm=$(echo "$spp" | awk 'BEGIN{OFS="_"}{print $1,$2}')

awk -F, -v SPP="$spp" '$2 == SPP {print $1, $3, $4}' $table \
    > $DIR/tmp/spp/fd_${nm}.txt

paste -d" " \
    $DIR/tmp/spp/fd_${nm}.txt \
    <(awk '{print $2, $3}' $DIR/tmp/spp/fd_${nm}.txt | gdallocationinfo -valonly \
    -wgs84 $DIR/mekongData/subcatchments_mekong.tif) \
    > $DIR/tmp/spp/fd_${nm}_sid.txt

awk '{print $4}' $DIR/tmp/spp/fd_${nm}_sid.txt | sort | uniq -c \
    > $DIR/tmp/spp/fd_${nm}_sid_summ.txt

for i in $(seq 1 $(wc -l < $DIR/tmp/spp/fd_${nm}_sid_summ.txt))
do
    num=$(awk -v row="$i" 'NR == row {print $1}' $DIR/tmp/spp/fd_${nm}_sid_summ.txt)
    cid=$(awk -v row="$i" 'NR == row {print $2}' $DIR/tmp/spp/fd_${nm}_sid_summ.txt)

    if [[ $num -eq 1 ]]
    then
    awk -v cid="$cid" '$4 == cid {print $0}' $DIR/tmp/spp/fd_${nm}_sid.txt \
        > $DIR/tmp/spp/tmp_${nm}_${i}.txt
    else
    #  here add code to select the record with the coordinates with higher
    # accuracy, i.e., more decimal points
    # num=0.000001
    # decimals=${num#*.}              #Removes the integer part and the dot (=000001)
    # decimalPointsCount=${#decimals} #Counts the length of resulting string (=6)
    awk -v cid="$cid" '$4 == cid {print $0}' $DIR/tmp/spp/fd_${nm}_sid.txt \
        | awk 'NR == 1' > $DIR/tmp/spp/tmp_${nm}_${i}.txt

    fi
done

cat $(find $DIR/tmp/spp -name "tmp_${nm}_*.txt") | awk '{print $1}' \
    > $DIR/tmp/spp/retain_indx_${nm}.txt 

rm $DIR/tmp/spp/tmp_${nm}_*.txt $DIR/tmp/spp/fd_${nm}*.txt

}

export -f rmDuplic
awk -F, 'NR > 1 {print $2}' $table | sort | uniq  \
    |  parallel -j 50 rmDuplic 

# join all tables together
cat $(find $DIR/tmp/spp -name "retain_indx*.txt") > $DIR/tmp/retain_indx.txt

# filter database to remove duplicates 
awk -F, 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
    $DIR/tmp/retain_indx.txt  $table \
    > $DIR/bio_modelling/mekong/species_db_mekong_nodup.csv

#### delete species with less than 5 occurrences

# check number of ocuurrences per species
awk -F, 'NR > 1 {print $2}' $DIR/bio_modelling/mekong/species_db_mekong_nodup.csv | sort | uniq -c

# check if they are eliminated
awk -F, 'NR == 1 || !/Anguilla marmorata/ && !/Helostoma temminckii/  && !/Tor sinensis/' $DIR/bio_modelling/mekong/species_db_mekong_nodup.csv | awk -F, '{print $2}' | sort | uniq -c \

# save new file    
awk -F, 'NR == 1 || !/Anguilla marmorata/ && !/Helostoma temminckii/  && !/Tor sinensis/' $DIR/bio_modelling/mekong/species_db_mekong_nodup.csv  \
    > $DIR/bio_modelling/mekong/species_db_mekong_final.csv

rm $DIR/bio_modelling/mekong/species_db_mekong.csv $DIR/bio_modelling/mekong/species_db_mekong_nodup.csv

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

#for i in $(seq 1 $(wc -l <  $DIR/mekongData/species_mekong_list.txt))
#do
#    spp=$(awk -v ROW="${i}" 'NR == ROW' $DIR/mekongData/species_mekong_list.txt)
#    spfnm=$(echo $spp | awk 'BEGIN{OFS="_";}{print $1, $2}')
#    echo "species,X,Y" > $TMP/${spfnm}_mekong_tmp.csv
#    awk -F, -v SPP="${spp}" 'BEGIN{OFS=",";}  $2 == SPP {print $2, $3, $4}' $occ \
#    >> $TMP/${spfnm}_mekong_tmp.csv 
#
#bash /mnt/shared/danube/dev_create_model_table.sh  \
#   $TMP/${spfnm}_mekong_tmp.csv \
#   $DIR/mekongData/predTBmekong.csv \
#   /mnt/shared/sosw/mekongData/sub_catchment.tif  \
#   10000  \
#   $TMP  \
#   /mnt/shared/sosw/bio_modelling/mekong/spp_input/${spfnm}.csv
#
#done
#
#
##### for some tables the first row has 0 as subcacthment id. Probably because
##### the point does not overlap with the subcatchments. delete that row
#
#for arch in $(find $DIR/bio_modelling/mekong/spp_input -name "*.csv")
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



#### procedure to reclasy the maps

## in r 

obj=list.files("randomForest_output/")

for (arch in obj){
    name=strsplit(arch, split="[.]")[[1]][1]
    load(paste0("randomForest_output/", arch))
    out=paste0("randomForest_output/", name, "_reclas.txt")
    write.table(results$predictions4map, out, row.names=FALSE, col.names=FALSE)
}

#sudo grass --text -c $DIR/mekongData/subcatchments_mekong.tif $DIR/grassdata/mekong

sudo grass --text $DIR/grassdata/jucar/PERMANENT

dir=/mnt/shared/sosw/bio_modelling/mekong/randomForest_output

r.in.gdal --o input=/mnt/shared/sosw/mekongData/subcatchments_mekong.tif output=subcat

paths=( $(ls $dir/*.RData) )
#tiff=${paths[0]}

for tiff in ${paths[@]:1}
do
spp=$(basename $tiff .RData)
outn=$(echo ${spp}_prob.tif)

awk '{print $1, "=", $2}' $dir/${spp}_reclas.txt > $dir/${spp}_reclas2.txt

r.reclass --o input=subcat output=$outn rules=$dir/${spp}_reclas2.txt

r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=Byte \
    format=GTiff nodata=255 input=$outn output=$dir/${outn}

done







