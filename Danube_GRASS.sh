#! /bin/bash

export STREAM=/mnt/shared/hydrography90m_v1_2022_all_data_OLD/stream_tiles_final20d_1p/all_streams_uniqID.vrt
gdalwarp  -t_srs EPSG:4326  -te 8 42 30 51 $STREAM $DIR/../tmp/danube_stream.tif \
    -co COMPRESS=LZW -co ZLEVEL=9

# global vrt file of flow accumulation
export FLOW=/mnt/shared/hydrography90m_v1_2022_all_data_OLD/CompUnit_flow_pos_noenlarge/flow_global.vrt

# global basin of computational units
export GLCOMPUNITS=/mnt/shared/hydrography90m_v1_2022_all_data_OLD/lbasin_compUnit_overview/lbasin_compUnit.tif
gdalwarp  -t_srs EPSG:4326  -te 8 42 30 51 $GLCOMPUNITS $DIR/tmp/danube_compUnit.tif \
    -co COMPRESS=LZW -co ZLEVEL=9

## path to microbasins
export MICROB=/mnt/shared/cog_layers_for_GeoNode/layers/sub_catchment_cog.tif
gdalwarp  -t_srs EPSG:4326  -te 8 42 30 51 $MICROB $DIR/tmp/danube_subcatchments.tif \
    -co COMPRESS=LZW -co ZLEVEL=9

###############################################################################
###############################################################################
###############################################################################

# create a GRASS GIS Database
sudo grass78 -c -text /mnt/shared/sosw/tmp/basins_sub.tif /mnt/shared/grassdata/danube 


###############################################################################
###############################################################################
###############################################################################


# acces the GRASS GIS Database after creation
sudo grass --text /mnt/shared/grassdata/danube/PERMANENT
 
r.in.gdal --o input=/mnt/shared/sosw/tmp/danube_subcatchments.tif output=subcat
r.in.gdal --o input=/mnt/shared/sosw/tmp/danube_compUnit.tif output=compunit

for cu in $(r.describe -n compunit) 
do
    r.mapcalc --overwrite "indx = if(compUnit==${cu}, subcat, null())"
    r.describe -1 -n indx > /mnt/shared/sosw/tmp/subcatchmentIDs_${cu}.txt
done



