export MNT=/mnt/shared
export CONN=/mnt/shared/danube/connect

gdalwarp  -co COMPRESS=LZW -co ZLEVEL=9 $(pkinfo  -te -i /mnt/shared/sosw/msk_danube.tif) \
    /mnt/shared/hydrography90m_v1_2022_all_data_OLD/dir_tiles_final20d_1p/all_dir_dis.vrt \
    /mnt/shared/sosw/dir_danube.tif -overwrite


grass --text --tmp-location /mnt/shared/sosw/dir_danube.tif <<'EOF'

# read binary lake layer with same extent as direction
r.external --o input=/mnt/shared/sosw/dir_danube.tif  output=dir_danube

	## calculate the sub-basin
r.water.outlet --overwrite input=dir_danube output=upstr \
    coordinates=19.35622,44.89216


# zoom to the region of interest (only upstream basin extent)
g.region -a --o zoom=upstr

    ##  Export the basin as tif file
    r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" \
        type=Byte  format=GTiff nodata=0 \
        input=upstr output=/mnt/shared/sosw/upstream.tif

EOF

# crreate ppoligon with border of area of interest (does not work in servers)
gdal_polygonize.py \
    $MNT/danube/connect/upstream.tif \
    -b 1 -f "GPKG" \
    $MNT/danube/connect/upstream.gpkg \
    upstream DN
#https://gis.stackexchange.com/questions/147820/st-intersect-with-ogr2ogr-and-spatialite


###  crop and mask data to Drina basin

export EXT=$(pkinfo -te -i /mnt/shared/danube/connect/upstream.tif)

####    SUBCATCHMENTS

gdalwarp  -co COMPRESS=LZW -co ZLEVEL=9 $(echo $EXT) \
$MNT/sosw/tmp/danube_subcatchments.tif $MNT/danube/connect/subcatchments_tmp.tif 

gdal_calc.py --co=COMPRESS=DEFLATE --co=ZLEVEL=9  \
    --NoDataValue=0 --overwrite --type=UInt32 \
    -A $MNT/danube/connect/subcatchments_tmp.tif \
    -B $MNT/danube/connect/upstream.tif \
    --outfile=$MNT/danube/connect/subcatchments.tif --calc="A*B"

rm $MNT/danube/connect/subcatchments_tmp.tif

####    STREAM NETWORK

gdalwarp  -co COMPRESS=LZW -co ZLEVEL=9 $(echo $EXT) \
$MNT/sosw/tmp/danube_stream.tif $MNT/danube/connect/stream1.tif 

gdal_calc.py --co=COMPRESS=DEFLATE --co=ZLEVEL=9  \
    --NoDataValue=0 --overwrite --type=UInt32 \
    -A $MNT/danube/connect/stream1.tif \
    -B $MNT/danube/connect/upstream.tif \
    --outfile=$MNT/danube/connect/stream.tif --calc="A*B"

rm $MNT/danube/connect/stream1.tif

####  SPECIES PREDICTION
gdalwarp  -co COMPRESS=LZW -co ZLEVEL=9 $(echo $EXT) \
    /mnt/shared/danube/models/output/rf/Barbus_barbus_rf_prob.tif \
    $MNT/danube/connect/Bb_tmp.tif

gdal_calc.py --co=COMPRESS=DEFLATE --co=ZLEVEL=9  \
    --NoDataValue=0 --overwrite --type=UInt32 \
    -A $MNT/danube/connect/Bb_tmp.tif\
    -B $MNT/danube/connect/upstream.tif \
    --outfile=$MNT/danube/connect/Barbus_barbus_rf_prob.tif --calc="A*B"

rm $MNT/danube/connect/Bb_tmp.tif


## copy files created in qgis locally
#scp /home/jaime/data/sosw/dams_drina.gpkg sv2:/mnt/shared/danube/connect
#scp /home/jaime/data/sosw/Barbus_barbus_drina.gpkg  sv2:/mnt/shared/danube/connect

# create temporal grass session
grass --text --tmp-location $CONN/upstream.tif

# read stream raster file
r.in.gdal input=$CONN/stream.tif out=stream --o
r.thin stream out=stream_thin
# make raster file a vector
r.to.vect -v input=stream_thin out=streamlines type=line
#v.out.ogr -s input=streamlines type=line out=$CONN/streamlines.gpkg format=GPKG

# import spp points
v.in.ogr --o input=$CONN/Barbus_barbus_drina.gpkg layer=Barbus_barbus_drina \
    out=bb type=point key=fid 

# import dam points
v.in.ogr --o input=$CONN/dams_drina.gpkg layer=dams_drina \
    out=dams type=point key=OBJECTID 

# joint stream vector network and species points
v.net --o input=streamlines points=bb output=stream_net operation=connect thresh=1
# calculate distance bewteen all points
v.net.allpairs -g --o input=stream_net out=all

## procedure to break to stream network and recalculate distances

v.to.rast --o input=dams type=point out=damr use=attr attribute_column=OBJECTID
r.buffer --o input=damr out=damrb distances=150

#r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" \
#        type=UInt32  format=GTiff nodata=0 \
#        input=damrb output=$CONN/damrb.tif

r.mapcalc --o "streamcut = if(isnull(damrb), stream_thin, null())"

#    r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" \
#        type=UInt32  format=GTiff nodata=0 \
#        input=streamcut output=$CONN/streamcut.tif
r.to.vect -v input=streamcut out=streamlinescut type=line

v.net --o input=streamlinescut points=bb output=stream_netcut operation=connect thresh=1
v.net.allpairs -g --o input=stream_netcut out=allcut


ids=( $(awk 'NR > 1 {print $2}' $CONN/all_pairs_bb.txt | sort | uniq))

for i in ${ids[@]}
do
    #awk -v n="$i" '$2 == n {print $3}' $CONN/all_pairs_bb_cut.txt | sort | uniq | wc -l
    n1=$(awk -v n="$i" '$2 == n {print $3}' $CONN/all_pairs_bb_cut.txt | sort | uniq | wc -l)
    echo "$i > $n1"
done | sort -n -k1


#####################################



setwd("~/data/connect")

t1 = read.table("all_pairs_bb.txt", header=TRUE)
t1$km = t1$cost / 1000

t2 = read.table("all_pairs_bb_cut.txt", header=TRUE)
t2$km = t2$cost / 1000

a = data.frame('connected' = t1$km)
b = data.frame('disconnected' = t2$km)
boxplot(dplyr::bind_rows(a, b))
boxplot(dplyr::bind_rows(a, b), ylab="Distance (km)", cex.lab=1.5, cex.axis=2)

# procedure in R to create the patches out of the broken network
library(terra)
r = rast("streamcut.tif")
p = patches(r, 8)
writeRaster(p, "patches8_01.tif")
unique(p[]) #   a total of 23 patches

# 1 13 4603.9 km
# 2872 9  5050.8 km
# 7044 3   318.9 km
# 8248 0 68.4 km
# 6147 0 2549.6 km
# 6329 1 14253.3 km
# 13084 2 7054.3 km
# 12800 9 611.7 km
# 15623 7 2318.4 km

#####################################################


r.in.gdal input=$CONN/subcatchments.tif output=subcatch --o
r.in.gdal input=$CONN/patches8_01.tif output=patches --o
r.in.gdal input=$CONN/Barbus_barbus_rf_prob.tif output=bbprob --o

### transform probability values to 1 and 0
r.mapcalc "bbbin = if(bbprob < 50, 0, 1)"

#    r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" \
#        type=Byte  format=GTiff nodata=255 \
#        input=bbbin output=$CONN/bbbin.tif


## extrcat stream segments of interest
r.mapcalc --o "allstr = if(bbbin == 1, subcatch, null())"

#    r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" \
#        type=UInt32  format=GTiff nodata=0 \
#        input=allstr output=$CONN/allstr1.tif

r.stats -n in=allstr > $CONN/ids_basin.txt

# extract from projection table the environmental data for the presences 
awk -F, 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
    $CONN/ids_basin.txt $CONN/../out/danube_predictTB.csv  >  $CONN/patch_all.csv


# patches to analize

for P in 1 2872 7044 8248 6157 6329 13084 12800 15623
do

    ###  identify segemt ids for each patch

    # si los parches = 1 
    r.mapcalc  --o "stridp_$P = if(patches == $P && bbbin == 1, subcatch, null())"

#    r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" \
#        type=UInt32  format=GTiff nodata=0 \
#        input=stridp output=$CONN/str_patch_1.tif

    r.stats -n in=stridp_$P > $CONN/ids_patch_$P.txt
    
    # extract from projection table the environmental data for the presences 
    awk -F, 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
        $CONN/ids_patch_$P.txt $CONN/../out/danube_predictTB.csv  >  $CONN/patch_$P.csv

done

####################################################

pa = read.csv("patch_all.csv")
p1 = read.csv("patch_1.csv") 
p2 = read.csv("patch_12800.csv") 
p3 = read.csv("patch_13084.csv") 
p4 = read.csv("patch_15623.csv") 
p5 = read.csv("patch_2872.csv") 
p6 = read.csv("patch_6157.csv") 
p7 = read.csv("patch_6329.csv") 
p8 = read.csv("patch_7044.csv") 
p9 = read.csv("patch_8248.csv") 

plot(pa[,2], pa[,8], 
        xlab="mean annual temperature", 
        ylab="annual precipitation",
        xlim=c(2800,2855),
        ylim=c(7500,20000))
points(p1[,2], p1[,8], col="blue")
points(p2[,2], p2[,8], col="blue")
points(p3[,2], p3[,8], col="orange")  ### upstream 2 connections
points(p4[,2], p4[,8], col="red")
points(p5[,2], p5[,8], col="red")
points(p6[,2], p6[,8], col="red")
points(p7[,2], p7[,8], col="green")  ### larger patch
points(p8[,2], p8[,8], col="red")
points(p9[,2], p9[,8], col="red")  ### extreme case  with 0 connection 68 km

