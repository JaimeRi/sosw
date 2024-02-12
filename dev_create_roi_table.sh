#! /bin/bash


#####  INPUTS

#  path to env tables
export ENVTB=/data/marquez/vignette/env_tables

# path, including the file, of the subcatchment for the roi
export SUBC=/mnt/shared/sosw/tmp/danube_subcatchments.tif

# output file
export OUTFILE=/data/marquez/vignette/out/projectionTB.csv

# file with the list of subcatchments IDs
export SUBCIDS=/data/marquez/vignette/out/subc_IDs.txt

# folder for temporal files
export TMP=/data/marquez/vignette/out

# path to file with list of variables (1st row) and tiles (2nd row)
export VT=$TMP/var_tiles.txt

# variables of interest
#export var=( bio1 c10_2020 spi stright )
export var=( awk 'NR == 1' $VT )

# tiles of interest
#export tiles=( h18v02 h18v04 h20v02 h20v04 )
export tiles=( awk 'NR == 2' $VT )


##### ANALYSIS

# extract the data for the variables of interest for the list of subcatchments
# of interest for all tiles
for TL in ${tiles[@]}
do
    for k in ${var[@]}
    do
    awk 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $SUBCIDS $ENVTB/${TL}_${k}.txt \
     | awk 'NR > 1 {for(i=1; i<=NF; i++) $i+=0}1' CONVFMT="%.3f" \
     >  $TMP/ENV_${TL}_${k}.txt
    done
done

###  Additional step to trim the number of decimal places to 3
awk -F, 'BEGIN{OFS=",";} NR > 1 {for(i=1; i<=NF; i++) $i+=0}1' CONVFMT="%.3f" \
    $TMP/pres_abs_tmp.csv > $OUTF

### join tables of different tiles for same variable
echo ${var[@]} | xargs -n 1 -P 1 bash -c $'

X=$1

listf=( $(find $TMP -name "ENV_*_${X}.txt")  )

cols=$( awk \'NR==1 {print NF}\' ${listf[0]} )

if [[ "$cols" -eq 6 ]]
then
    echo "subCid ${X}_min ${X}_max ${X}_range ${X}_mean ${X}_sd" > $TMP/ENVaggreg_${X}.txt
    cat ${listf[@]} >> $TMP/ENVaggreg_${X}.txt
    sed \'/subcID/d\' $TMP/ENVaggreg_${X}.txt > $TMP/ENVaggreg_${X}f.txt
    sort -g $TMP/ENVaggreg_${X}f.txt > $TMP/ENVaggreg_${X}.txt
    rm $TMP/ENVaggreg_${X}f.txt
else
    echo "subCid ${X}" > $TMP/ENVaggreg_${X}.txt
    cat ${listf[@]} >> $TMP/ENVaggreg_${X}.txt
    sed \'/subcID/d\' $TMP/ENVaggreg_${X}.txt > $TMP/ENVaggreg_${X}f.txt
    sort -g $TMP/ENVaggreg_${X}f.txt > $TMP/ENVaggreg_${X}.txt
    rm $TMP/ENVaggreg_${X}f.txt
fi

' _

## join all the environmental variables together
paste -d" " $(find $TMP/ENVaggreg_*.txt) > $TMP/ENVpred.txt


###  join tables of all variables together

## Chunk to delete the subCid column for each table before joining
read -a header < $TMP/ENVpred.txt # read first line into array "header"
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
    <(printf "%s\n" subcID $(awk 'NR > 1 {print $1}' $TMP/ENVpred.txt)) \
    <(tr -s ' ' ',' < $TMP/ENVpred.txt | \
    cut -d"," --complement -f $(echo "${joined%,}")) \
    > $TMP/projectionTBtmp.csv

#### remove duplicates
awk -F, '!a[$0]++'  $TMP/projectionTBtmp.csv > $OUTFILE

# remove temporal files
rm $TMP/projectionTBtmp.csv
rm $TMP/ENV*.txt  
rm $VT

exit
