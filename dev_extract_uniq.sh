#! /bin/bash


# path, including the file, of the subcatchment for the roi
#export SUBC=/mnt/shared/sosw/tmp/danube_subcatchments.tif
#export SUBC=/mnt/shared/sosw/tmp/sub_test.tif
export  SUBC=$1

# path to output folder
#export OUT=/data/marquez/vignette/out
export OUT=$2

grass -f --text --tmp-location $SUBC  <<'EOF'

r.in.gdal --o input=$SUBC output=subcat
r.stats -n input=subcat  > $OUT/subc_IDs.txt

EOF

