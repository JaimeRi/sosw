#!/bin/bash

export DAN=/mnt/shared/sosw/danube/
export DIR=/mnt/shared/sosw/danube/subc
export PRED=/mnt/shared/sosw/danube/env/pred

export BIO=/mnt/shared/regional_unit_bio
export BIOF=/mnt/shared/regional_unit_tables_bio_fut



#### aquí empezar un for loop para correr por RU
for CU in $(find $DIR -name "subcatchment*.txt" | awk -F[_.] '{print $2}' )
#for CU in 66 159
do
    #export CU=159

# make a list of subcatchment ids 
export IDS=$DIR/subcatchmentIDs_${CU}.txt

# CHUNK FOR THE BIOCLIMATE VARIABLES

time for k in bio1 bio12 bio15
do
# extract the data for the variable of interest for the list of subcatchments
awk 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
    $IDS $BIO/CU_$CU/stats_${CU}_${k}.txt \
    | awk '{print $1, $5, $6}' >  $PRED/ENV_${CU}_${k}.txt
done

# CHUNK FOR LAND COVER DATA
#awk 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
#    $IDS $DIR/env/stats_${CU}_LCprop.txt \
#    | awk 'NR == 1 { for (i=1; i<=NF; i++) {f[$i] = i}  } \
#        {print $(f["c90_y2020"]),$(f["c110_y2020"]), $(f["c80_y2020"]) }' \
#    > $PRED/ENV_${CU}_LCprop.txt

# CHUNK FOR ALL OTHER VARIABLES

for z in spi slopcmax
do
# extract the data for the variable of interest for the list of subcatchments
awk 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
    $IDS $DAN/env/stats_${CU}_${z}.txt \
    | awk '{print $1, $5, $6}' >  $PRED/ENV_${CU}_${z}.txt
done

done


################
### join tables of different RU for same variable
variables=( $(find $PRED -name "ENV_*.txt" | awk -F[_.] '{print $3}' \
    | sort | uniq) )

# the chunk below assumes the mean and sd were extracted
# add if statement for Land Cover data
echo ${variables[@]} | xargs -n 1 -P 5 bash -c $'
X=$1
echo "subCid ${X}_mean ${X}_sd" > $PRED/ENVaggreg_${X}.txt
cat $(find $PRED -name "ENV_*_${X}.txt") >> $PRED/ENVaggreg_${X}.txt
sed \'/subcID/d\' $PRED/ENVaggreg_${X}.txt > $PRED/ENVaggreg_${X}f.txt
sort -g $PRED/ENVaggreg_${X}f.txt > $PRED/ENVaggreg_${X}.txt
rm $PRED/ENVaggreg_${X}f.txt 
' _

paste -d" " $(find $PRED/ENVaggreg_*.txt) > $PRED/ENVpred.txt

################
###  join tables of all variables together and match with species ocuurrence
###  table

## Chunk to delete the subCid column for each table before joining

read -a header < $PRED/ENVpred.txt # read first line into array "header"
declare -a arr=() # array to store the position in which the subCid name is

for i in ${!header[@]}               # iterate through array indexes
do
    if [ "${header[i]}" = "subCid" ]    # find column equal the pattern
    then
        arr+=( "$[++i]"  )
    fi
done

printf -v joined '%s,' "${arr[@]}"

paste -d"," \
    <(printf "%s\n" subcID $(awk 'NR > 1 {print $1}' $PRED/ENVpred.txt)) \
    <(tr -s ' ' ',' < $PRED/ENVpred.txt | \
    cut -d"," --complement -f $(echo "${joined%,}")) \
    > $PRED/danube_predict.csv

rm $PRED/ENV*.txt  


