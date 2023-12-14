#! /bin/bash

### Given a selection of variables available in the computational Units zip 
### folders, the script unzip the wished files to a given destination folder 


export DIR=/mnt/shared/sosw/danube
export SOS=/mnt/shared/sosw

export BIO=/mnt/shared/regional_unit_bio
export BIOF=/mnt/shared/regional_unit_tables_bio_fut
export VAR=/mnt/shared/regional_unit_tables


# global basin of computational units
#export GLCOMPUNITS=$2
export GLCOMPUNITS=$SOS/tmp/danube_compUnit.tif

#export sppTB=$1
export sppTB=$SOS/sppTB/fish_danube.csv

### check which RU
RUS=( $(awk -F, 'FNR > 1 {print $3, $4}' $sppTB | gdallocationinfo -valonly -geoloc $GLCOMPUNITS | sort | uniq) )

### check which variables for which the unzip has to be done
ENV=( LCprop spi slopcmax )


for i in ${RUS[@]}; do for j in ${ENV[@]}; do echo $i $j; done; done | xargs -n 2 -P 6 bash -c '
CU=$1
envv=$2
sudo unzip -j $VAR/CU_${CU}.zip "*/stats_${CU}_${envv}.txt" -d $DIR/env
' _ 
