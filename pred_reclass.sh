#! /bin/bash

export DIR=/mnt/shared/sosw/danube

grass -f --text /data/marquez/grassdata/danube/PERMANENT <<'EOF'

#export DIR=/mnt/shared/sosw/danube
export TB=/mnt/shared/sosw/danube/env/pred/prediction.txt

cat \
    <(awk 'NR > 1 {print $1, "=" , $3}' $TB) \
    <(echo "* = NULL") > $DIR/env/pred/rules.txt

r.reclass input=subcat output=recl_raster rules=$DIR/env/pred/rules.txt --overwrite

r.out.gdal input=recl_raster output=$DIR/env/pred/pred_Huso_huso_test.tif type=Int32  format=GTiff nodata=-9999  --o -f -m -c createopt="COMPRESS=DEFLATE,ZLEVEL=9"

EOF
