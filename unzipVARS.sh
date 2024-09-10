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



##########################################
##########################################
DIR=/mnt/shared/danube

bio=( $(ls /mnt/shared/EnvTablesTiles/Climate/present) )
hyd=( $(ls /mnt/shared/EnvTablesTiles/Hydrography90m) )
cob=( c10 c20 c30 c40 c50 c60 c70 c80 c90 c100 )

ENV=( ${bio[@]} ${hyd[@]} ${cob[@]} )


for t in h18v02 h18v04 h20v02 h20v04
do
    for v in ${bio[@]}
    do
        sudo unzip  /mnt/shared/EnvTablesTiles/Climate/present/${v}/${t}_${v}.zip \
           -d $DIR/env_tiles
    done

    for v in ${hyd[@]}
    do
        sudo unzip /mnt/shared/EnvTablesTiles/Hydrography90m/${v}/${t}_${v}.zip \
           -d  $DIR/env_tiles
    done

    for v in ${cob[@]}
    do
        sudo unzip /mnt/shared/EnvTablesTiles/LandCover/${v}/${t}_${v}_2020.zip \
           -d  $DIR/env_tiles
    done
done
