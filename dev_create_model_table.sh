#! /bin/bash


### INPUT

## file with table with species name and coordinates
## be sure no subcatchment duplicates!!!
export SPP=$1
export SPP=/data/marquez/vignette/out/spp.csv

# crear archivo solo con datos de la especie de interes
#awk -F, -v SPP="Huso huso" 'BEGIN{OFS=",";} NR == 1 || $2 == SPP {print $2, $3, $4}' \
#    /mnt/shared/sosw/sppTB/fish_danube.csv > /data/marquez/vignette/out/spp.csv 

## projection table created with function create_roi_table.sh
export PTB=$2
export PTB=/data/marquez/vignette/out/danube_predict.csv

## path to subcatchment file
export SUBC=$3
export SUBC=/mnt/shared/sosw/tmp/danube_subcatchments.tif

### number of pseudoabsences
export ABS=$4
export ABS=10000

## temporal folder
export TMP=$5
export TMP=/data/marquez/vignette/out

## output file, the model table
export OUTF=/data/marquez/vignette/out/model_table.csv


### ANALYSIS

## Presence Data

# species name
export SPPN=$(awk -F, 'NR == 2 {print $1}' $SPP)

# line to calculate the number of rows in the presence table (no counting header)
C=$( awk -F, '{print NR}' $SPP | tail -n2 | head -n1 )

# adicionar al archivo anterioir las columnas de subcatchment id
# add new columns: subcatchment id and presence/absence column (add 1 presences)
paste -d "," $SPP \
    <(printf "%s\n" SubcatchID $(awk -F, 'FNR > 1 {print $2, $3}' $SPP | gdallocationinfo -valonly -geoloc $SUBC)) \
    <(printf "%s\n" PresAbs $(printf '1%.0s\n' $(eval "echo {1.."$(($C))"}") )) \
    > $TMP/tmp1.csv

# extrcat from projection table the environmental data for the presences 
awk -F, 'NR==FNR {a[$4]; next} FNR==1 ||  $1 in a' \
    $TMP/tmp1.csv $PTB  >  $TMP/tmp2.csv

# join tables
paste -d"," \
    <(sort -t, -g -k4 $TMP/tmp1.csv) \
    <(sort -t, -g -k1 $TMP/tmp2.csv) \
    | cut -d"," --complement -f 6 \
    > $TMP/pa_env_tmp.csv

## (Pseudo)absences

grass --text --tmp-location $SUBC <<'EOF' 
r.external input=$SUBC out=sub
r.random -s input=sub vector=abs npoints=$ABS
v.out.ascii -c input=abs column=value precision=3 separator=space \
    | awk 'BEGIN{OFS=","} NR > 1  {print int($4), $1, $2}' \
    > $TMP/abs_tmp1.csv
EOF

# remove duplicates
awk -F, '!a[$1]++'  $TMP/abs_tmp1.csv > $TMP/abs_tmp2.csv

# remove duplicates in absences if subcatchment id also exist in presence
awk -F, 'NR==FNR {a[$4]; next} $1 in a' \
    $TMP/tmp1.csv $TMP/abs_tmp2.csv > $TMP/abs_rm.txt

if [[ $(wc -l < $TMP/abs_rm.txt) -gt 0  ]]
then
    awk -F, '{print $1}' $TMP/abs_rm.txt > $TMP/abs_indx_rm.txt
    awk -F, 'NR==FNR {a[$1]; next} ! ($1 in a)' \
        $TMP/abs_indx_rm.txt $TMP/abs_tmp2.csv > $TMP/abs_tmp3.csv
else
    cp $TMP/abs_tmp2.csv $TMP/abs_tmp3.csv
fi

#  number of rows in the absence table 
A=$( wc -l < $TMP/abs_tmp3.csv )

# create preliminary table with same format as tmp1
paste -d"," \
    <(printf "${SPPN}%.0s\n" $(eval "echo {1.."$(($A))"}")) \
    <(awk -F, 'BEGIN{OFS=","} {print $2, $3, $1}' $TMP/abs_tmp3.csv) \
    <(printf '0%.0s\n' $(eval "echo {1.."$(($A))"}")) \
    > $TMP/abs_tmp4.csv

# extract from projection table the environmental data for the presences 
awk -F, 'NR==FNR {a[$4]; next} $1 in a' \
    $TMP/abs_tmp4.csv $PTB  >  $TMP/abs_tmp5.csv

# join tables
paste -d"," \
    <(sort -t, -g -k4 $TMP/abs_tmp4.csv) \
    <(sort -t, -g -k1 $TMP/abs_tmp5.csv) \
    | cut -d"," --complement -f 6 \
    > $TMP/abs_env_tmp.csv

### Join presences and absences
cat $TMP/pa_env_tmp.csv $TMP/abs_env_tmp.csv > $OUTF 

### remove temporal files
rm $TMP/*tmp* 


exit