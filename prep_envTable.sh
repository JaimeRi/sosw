#! /bin/bash

####  the script prepares a table of environmental variables for a particular species

export DIR=/mnt/shared/sosw/danube
export TMP=$DIR/tmp
export SOS=/mnt/shared/sosw

export BIO=/mnt/shared/regional_unit_bio
export BIOF=/mnt/shared/regional_unit_tables_bio_fut
export VAR=/mnt/shared/regional_unit_tables


# global basin of computational units
#export GLCOMPUNITS=$2
export GLCOMPUNITS=$SOS/tmp/danube_compUnit.tif

## path to microbasins
#export MICROB=$4
export MICROB=$SOS/tmp/danube_subcatchments.tif


#export sppTB=$1
export sppTB=$SOS/sppTB/fish_danube.csv

### check which RU
#RUS=( $(awk -F, 'FNR > 1 {print $3, $4}' $sppTB | gdallocationinfo -valonly -geoloc $GLCOMPUNITS | sort | uniq) )

####  aqui empieza la función para correr por especie

export sp=$3
export sp="Huso huso"

# nombre pegado por under score
export nm=$(echo $sp | awk 'BEGIN{OFS="_"}{print $1,$2}')

mkdir $TMP/${nm}
export SPPDIR=$TMP/${nm}

# crear archivo solo con datos de la especie de interes
awk -F, -v SPP="$sp" 'NR == 1 || $2 == SPP' $sppTB > $SPPDIR/tmp1.csv 

# line to calculate the number of rows in the presence table (no counting header)
C=$( awk -F, '{print NR}' $SPPDIR/tmp1.csv | tail -n2 | head -n1 )

# adicionar al archivo anterioir las columnas de compUnit y subcatchment id
paste -d "," $SPPDIR/tmp1.csv \
    <(printf "%s\n" CompUnit $(awk -F, 'FNR > 1 {print $3, $4}' $SPPDIR/tmp1.csv | gdallocationinfo -valonly -geoloc $GLCOMPUNITS)) \
    <(printf "%s\n" SubcatchID $(awk -F, 'FNR > 1 {print $3, $4}' $SPPDIR/tmp1.csv | gdallocationinfo -valonly -geoloc $MICROB)) \
    <(printf "%s\n" PresAbs $(printf '1%.0s\n' $(eval "echo {1.."$(($C))"}") )) \
    > $SPPDIR/tmp2.csv

#######  Generar pseudo-Ausencias

# access the grass sesion
sudo grass --text --tmp-location $GLCOMPUNITS <<'EOF' 

export sp="Huso huso"
export nm=$(echo $sp | awk 'BEGIN{OFS="_"}{print $1,$2}')

r.external input=/mnt/shared/sosw/tmp/danube_compUnit.tif out=comp
r.external input=/mnt/shared/sosw/tmp/danube_subcatchments.tif out=sub
r.random -s --o input=sub cover=comp npoints=10000 vector=abs
v.out.ascii -c input=abs columns=* separator=comma precision=0 \
    > /mnt/shared/sosw/danube/tmp/${nm}/abs1.csv
EOF

# remove subcid replicates
awk -F, 'NR > 1 {print  $4}' $SPPDIR/abs1.csv | sort | uniq -c | \
    awk '$1 == 1 {print $2}' > $SPPDIR/rmabs.txt 

awk -F, 'NR==FNR {a[$1]; next} FNR==1 ||  $4 in a' \
    $SPPDIR/rmabs.txt $SPPDIR/abs1.csv \
    >  $SPPDIR/abs2.csv

rm $SPPDIR/rmabs.txt


# prepare format
# add IGB id if the previous is the final version
Aseq=$(echo "$(wc -l < $SPPDIR/abs2.csv) - 1 " | bc)

paste -d"," \
    <(printf "%s\n" IGB_id $(printf "A%1g\n" $(seq 1 $Aseq))) \
    <(printf "%s\n" PresAbs $(printf "0%.0s\n" $(eval "echo {1.."$(("$Aseq"))"}"))) \
    <(awk -F, 'BEGIN{OFS=",";}{print $5, $4}' $SPPDIR/abs2.csv) \
    > $SPPDIR/abs3.csv

#### Join presence and absence

# columns presence $1, $8, $6, $7
# columns absence $3, $5, $4
cat <(awk -F, 'BEGIN{OFS=",";}{print $1, $8, $6, $7}' tmp2.csv) \
    <(awk -F, 'NR>1' abs3.csv) > $SPPDIR/pa1.csv

# check again for no duplicated in subCid
awk -F, '{print $4}' $SPPDIR/pa1.csv | sort | uniq -c | awk '$1 > 1'

cp $SPPDIR/pa1.csv $SPPDIR/pa.csv

rm $SPPDIR/pa1.csv $SPPDIR/tmp*.csv $SPPDIR/abs*.csv


#### aquí empezar un for loop para correr por RU

for CU in $(awk -F, 'FNR > 1 {print $3}' $SPPDIR/pa.csv | sort | uniq)
do
    #export CU=159

# subset of occurrences for each RU and extrcat only subcatchment id
# make a list of subcatchment ids but without duplicating ids
awk -F, -v CU="$CU"  '$3 == CU {print $4}' $SPPDIR/pa.csv | sort | uniq \
    > $SPPDIR/subcid_${CU}.txt


# CHUNK FOR THE BIOCLIMATE VARIABLES

for k in bio1 bio12 bio15
do
# extract the data for the variable of interest for the list of subcatchments
awk 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
    $SPPDIR/subcid_${CU}.txt $BIO/CU_$CU/stats_${CU}_${k}.txt \
    | awk '{print $1, $5, $6}' >  $SPPDIR/ENV_${CU}_${k}.txt
done

# CHUNK FOR LAND COVER DATA
#awk 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
#    $SPPDIR/subcid_${CU}.txt $DIR/env/stats_${CU}_LCprop.txt \
#    | awk 'NR == 1 { for (i=1; i<=NF; i++) {f[$i] = i}  } \
#        {print $(f["c90_y2020"]),$(f["c110_y2020"]), $(f["c80_y2020"]) }' \
#    > $SPPDIR/ENV_${CU}_LCprop.txt

# CHUNK FOR ALL OTHER VARIABLES

for z in spi slopcmax
do
# extract the data for the variable of interest for the list of subcatchments
awk 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
    $SPPDIR/subcid_${CU}.txt $DIR/env/stats_${CU}_${z}.txt \
    | awk '{print $1, $5, $6}' >  $SPPDIR/ENV_${CU}_${z}.txt
done

done

################
### join tables of different RU for same variable
variables=( $(find $SPPDIR -name "ENV_*.txt" | awk -F[_.] '{print $4}' \
    | sort | uniq) )

# the chunk below assumes the mean and sd were extracted
# add if statement for Land Cover data
echo ${variables[@]} | xargs -n 1 -P 5 bash -c $'
X=$1
echo "subCid ${X}_mean ${X}_sd" > $SPPDIR/ENVaggreg_${X}.txt
cat $(find $SPPDIR -name "ENV_*_${X}.txt") >> $SPPDIR/ENVaggreg_${X}.txt
sed \'/subcID/d\' $SPPDIR/ENVaggreg_${X}.txt > $SPPDIR/ENVaggreg_${X}f.txt
sort -g $SPPDIR/ENVaggreg_${X}f.txt > $SPPDIR/ENVaggreg_${X}.txt
rm $SPPDIR/ENVaggreg_${X}f.txt 
' _

paste -d" " $(find $SPPDIR/ENVaggreg_*.txt) > $SPPDIR/ENVall.txt

################
###  join tables of all variables together and match with species ocuurrence
###  table

## Chunk to delete the subCid column for each table before joining

read -a header < $SPPDIR/ENVall.txt # read first line into array "header"
declare -a arr=() # array to store the position in which the subCid name is

for i in ${!header[@]}               # iterate through array indexes
do
    if [ "${header[i]}" = "subCid" ]    # find column equal the pattern
    then
        arr+=( "$[++i]"  )
    fi
done

printf -v joined '%s,' "${arr[@]}"

# JOIN
paste -d"," \
    <(sort -t, -g -k4 $SPPDIR/pa.csv) \
    <(tr -s ' ' ',' < $SPPDIR/ENVall.txt | cut -d"," --complement -f $(echo "${joined%,}")) \
    > $SPPDIR/pa_env_${nm}.csv

rm $SPPDIR/ENV_*.txt $SPPDIR/ENV*.txt  $SPPDIR/pa.csv $SPPDIR/subcid_*.txt

