#! /bin/bash

#  path to env tables
export ENVTB=/data/marquez/vignette/env_tables

# path, including the file, of the subcatchment for the roi
export SUBC=/mnt/shared/sosw/tmp/danube_subcatchments.tif

# output folder
export OUT=/data/marquez/vignette/out


grass -f --text --tmp-location $SUBC

r.in.gdal --o input=$SUBC output=subcat
r.describe -1 -n subcat > $OUT/subc_IDs.txt



#### run per tile

# identify tiles
tiles=$(ls $ENVTB | awk -F_ '{print $1}' | sort | uniq )


time for TL in $tiles
do
    for k in bio1 spi
    do
    # extract the data for the variable of interest for the list of subcatchments
    awk 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $OUT/subc_IDs.txt $ENVTB/${TL}_${k}.txt \
     >  $OUT/ENV_${TL}_${k}.txt
    done

done


################
### join tables of different itiles for same variable
variables=( $(find $OUT -name "ENV_*.txt" | awk -F[_.] '{print $3}' \
    | sort | uniq) )

# the chunk below assumes the mean and sd were extracted
# add if statement for Land Cover data
echo ${variables[@]} | xargs -n 1 -P 2 bash -c $'
X=$1
echo "subCid ${X}_min ${X}_max ${X}_range ${X}_mean ${X}_sd" > $OUT/ENVaggreg_${X}.txt
cat $(find $OUT -name "ENV_*_${X}.txt") >> $OUT/ENVaggreg_${X}.txt
sed \'/subcID/d\' $OUT/ENVaggreg_${X}.txt > $OUT/ENVaggreg_${X}f.txt
sort -g $OUT/ENVaggreg_${X}f.txt > $OUT/ENVaggreg_${X}.txt
rm $OUT/ENVaggreg_${X}f.txt 
' _

paste -d" " $(find $OUT/ENVaggreg_*.txt) > $OUT/ENVpred.txt

################
###  join tables of all variables together


## Chunk to delete the subCid column for each table before joining

read -a header < $OUT/ENVpred.txt # read first line into array "header"
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
    <(printf "%s\n" subcID $(awk 'NR > 1 {print $1}' $OUT/ENVpred.txt)) \
    <(tr -s ' ' ',' < $OUT/ENVpred.txt | \
    cut -d"," --complement -f $(echo "${joined%,}")) \
    > $OUT/danube_predict_md.csv

#### remove duplicates
awk -F, '!a[$0]++'  $OUT/danube_predict_md.csv > $OUT/danube_predict.csv

# remove temporal files
rm $OUT/ENV*.txt  

