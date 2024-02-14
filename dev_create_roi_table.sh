#! /bin/bash


#### R parameters

###  1. list of variables of interest
###  2. list of statistics of interest :  default "ALL"
###  3. list of neccesary tiles IDs
###  The R function needs to get these previous parameters and create a file
###  with three rows to store each of the three parameters consecutively
###  4. path to environmental tables
###  5. path, including the file, of the subcatchment for the region of interest
###  6. file with the list of subcatchments IDs
###  7. path to output file
###  8. path to temporal file
###  9. number of cores (if possible) to run internal process in parallel: 
###     Highest number needed = (n.tiles * n.variables)

#####  PARAMETERS

# path to file with list of variables (1st row) and tiles (2nd row) and 
# statistics of interest (e.g. mean, sd)
export VT=$1
#export VT=$TMP/var_stats_tiles.txt

# variables of interest
#export var=( bio1 c10_2020 spi stright )
export var=( $(awk 'NR == 1' $VT) )

# select summary statistics
# ALL
# c(mean, sd)
export SS=( $(awk 'NR == 2' $VT) )

# tiles of interest
#export tiles=( h18v02 h18v04 h20v02 h20v04 )
export tiles=( $(awk 'NR == 3' $VT) )

#  path to environmental tables for each tile
export ENVTB=$2
#export ENVTB=/data/marquez/vignette/env_tables

#        cp /mnt/shared/EnvTablesTiles/LandCover/c10/h18v02_c10_2020.zip \
#            /data/marquez/vignette/env_tables

# path, including the file, of the subcatchment for the roi
export SUBC=$3
#export SUBC=/mnt/shared/sosw/tmp/danube_subcatchments.tif

# file with the list of subcatchments IDs
export SUBCIDS=$4
#export SUBCIDS=/data/marquez/vignette/out/subc_IDs.txt

#library(terra)
#library(raster)
#r = raster("/mnt/shared/sosw/tmp/danube_subcatchments.tif")
#u = terra::unique(r)
#write.table(u, file="/data/marquez/vignette/out/subc_IDs.txt", col.names=FALSE, row.names=FALSE)

# output file
export OUTFILE=$5
#export OUTFILE=/data/marquez/vignette/out/projectionTB2.csv

# folder to store temporal files: this folder is defined within the R function
export TMP=$6
#export TMP=/data/marquez/vignette/out

# number of cores to run the extraction of information (rows) from tile tables
# (n.tiles * n.variables)
export NCORES=$7
#export NCORES=16


##### ANALYSIS

##################
# Move through each tiles and extract the subcatchment of interests for
# all variables 

subsetTB(){
    TL=$1  # tile
    k=$2   # variable
    awk 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
     $SUBCIDS $ENVTB/${TL}_${k}.txt \
     | awk 'NR > 1 {for(i=1; i<=NF; i++) $i+=0}1' CONVFMT="%.3f" \
     >  $TMP/ENV_${TL}_${k}.txt
}

export -f subsetTB
time parallel -j $NCORES subsetTB ::: ${tiles[@]} ::: ${var[@]}


##################
#   Subset the tables if user is only interested in few statistics (e.g., mean)

if [[ "${SS[@]}" != 'ALL' ]]    # run only if user do not select ALL
then
    for TB in $(find $TMP -name "ENV_*.txt")
    do
        NR=$(awk 'NR == 1 {print NF}' $TB)
        [[ "$NR" != 6 ]] && continue
        
        RAND_STRING=$(xxd -l 8 -c 32 -p < /dev/random)
       
        read -a header <  $TB   # read first line into array "header"
        declare -a arr=() # array to store the position in which the subCid name is

        for s in ${SS[@]}
        do
            for i in ${!header[@]}               # iterate through array indexes
            do
                if [ "${header[i]}" = "$s" ]    # find column equal the pattern
                then
                    arr+=( "$[++i]"  )
                fi
            done
        done

        printf -v joined '%s,' "${arr[@]}"

        cut -d" " -f $(echo "1,${joined%,}") $TB \
        > $TMP/ENV_${RAND_STRING}.txt

        mv $TMP/ENV_${RAND_STRING}.txt $TB
    done
fi


##################
### join all tables for the same variable available in the different tiles
time echo ${var[@]} | xargs -n 1 -P $NCORES bash -c $'

X=$1

listf=( $(find $TMP -name "ENV_*_${X}.txt")  )

cat ${listf[@]} > $TMP/aggreg_${X}.txt
awk \'FNR == 1; FNR > 1 && /^[0-9]/\' $TMP/aggreg_${X}.txt \
    > $TMP/aggreg_${X}_tmp1.txt
sort -g $TMP/aggreg_${X}_tmp1.txt > $TMP/aggreg_${X}.txt

read -a header < $TMP/aggreg_${X}.txt

# set the right name for the header (e.g. bio1_mean)
if [[ ${#header[@]} -gt 2 ]]; then

    nof=( ${header[@]:1} )
    allh=( ${header[0]} $(echo "${nof[@]/#/${X}_}") )
    echo ${allh[@]} > $TMP/aggreg_${X}_tmp2.txt
    awk \'NR>1\' $TMP/aggreg_${X}.txt  >> $TMP/aggreg_${X}_tmp2.txt
    mv $TMP/aggreg_${X}_tmp2.txt $TMP/aggreg_${X}.txt
fi

rm $TMP/aggreg_${X}_tmp*.txt

' _


#################
## join all the environmental variables together
paste -d" " $(find $TMP/aggreg_*.txt) > $TMP/aggreg_all.txt

## the previous line creates repetition of the subcID column
## Chunk to delete the subCid column 
read -a header < $TMP/aggreg_all.txt  # read first line into array "header"
declare -a arr=() # array to store the position in which the subCid columns are

for i in ${!header[@]}               # iterate through array indexes
do
    if [ "${header[i]}" = "subcID" ]    # find column equal the pattern
    then
        arr+=( "$[++i]"  )
    fi
done

printf -v joined '%s,' "${arr[@]:1}"  # remove from array the first column

##  remove subcID columns (except column 1) and make the table comma separated
cut -d" " --complement -f $(echo "${joined%,}") $TMP/aggreg_all.txt \
    | tr -s ' ' ',' > $TMP/aggreg_all_trim.csv

##  remove duplicates and create final output table
awk -F, '!a[$0]++'  $TMP/aggreg_all_trim.csv > $OUTFILE



#########################
# remove temporal files
rm $TMP/aggreg*
rm $TMP/ENV*  
rm $VT 

exit
