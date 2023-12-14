export DIR=/mnt/shared/sosw

###   make masks of ROI for the study cases

### RHINE
printf "gid,WKT\n1,\"POLYGON((-0.55 55.01, -0.55 43.92, 15.35 43.92, 15.35 55.01, -0.55 55.01))\"\n" \
    > $DIR/roi_rhine.csv

ogr2ogr -f "GPKG" $DIR/roi_rhine.gpkg -dialect sqlite -sql "SELECT gid, ST_GeomFromText(WKT,4326) FROM roi_rhine" -nln roi -a_srs EPSG:4326 $DIR/roi_rhine.csv

gdal_rasterize -a gid -l roi -tr 0.083 0.083 $DIR/roi_rhine.gpkg $DIR/msk_rhine.tif -ot Byte 

### JUCAR
# Júcar  POLYGON((-9.86 46.4, 6.29 46.4, 6.29 34.74, -9.86 34.74, -9.86 46.4))

printf "gid,WKT\n1,\"POLYGON((-9.86 46.4, 6.29 46.4, 6.29 34.74, -9.86 34.74, -9.86 46.4))\"\n" \
    > $DIR/roi_jucar.csv

ogr2ogr -f "GPKG" $DIR/roi_jucar.gpkg -dialect sqlite -sql "SELECT gid, ST_GeomFromText(WKT,4326) FROM roi_jucar" -nln roi -a_srs EPSG:4326 $DIR/roi_jucar.csv

gdal_rasterize -a gid -l roi -tr 0.083 0.083 $DIR/roi_jucar.gpkg $DIR/msk_jucar.tif -ot Byte 

###############################################################################
#######   Fish  #######################
###############################################################################

######################  Fishnet2

# CHeck on http://www.fishnet2.net/search.aspx?c=OS using the polygon below
# POLYGON((-9.86 46.4, -9.86 34.74, 6.29 34.74, 6.29 46.4, -9.86 46.4))

export fishnet=$DIR/bio_down/SearchResults--8585111637345333768_fishnet2_jucar.txt

# check number of columns in table and if they are consistent given the field separator
#awk -F"\t" '{print NF}' $fishnet | sort | uniq -c

# check name of headers and their position
#printf "%s\n" $(awk -F"\t" 'NR == 1'  $fishnet | sed 's/,/ /g') | cat -n  

paste -d"," \
    <( awk -F"\t"  'BEGIN{OFS=","}{print $1, $2, $3}' ${fishnet} ) \
    <(awk -F"\t" '{print $5}' $fishnet | awk -F[,\(\)] '{print $1, $2}' | awk '{print $1, $2}') \
    <(awk -F"\t"  'BEGIN{OFS=","}{print $6, $10, $9, $21}' $fishnet) \
    > $DIR/tmp/fisnnet_sub1.csv

# remove records with empty species names
# remove records where species names are given as sp or sp.
# remove records with species names as Unkown
# remove records where species names are given as x
# remove strange names

### RHINE
awk -F, '$4 != "" && $4 !~ / sp/ && $4 !~ /?/ && $4 !~ /Unknown/ && $4 !~ / x/ && $4 != "LeuciscusÊleuciscus LeuciscusÊleuciscus" && $4 != "RutilusÊrutilus RutilusÊrutilus"' \
    $DIR/tmp/fisnnet_sub1.csv > $DIR/tmp/fisnnet_sub2.csv

### JUCAR
awk -F, '$4 != " " && $4 !~ / sp/ && $4 !~ /?/ && $4 !~ /Unknown/ && $4 !~ / x/  && $4  !~ /Unidentified/' \
    $DIR/tmp/fisnnet_sub1.csv > $DIR/tmp/fisnnet22_sub2.csv


# removes records where species names only have the genus
paste -d"," $DIR/tmp/fisnnet22_sub2.csv \
    <(awk -F, '{print $4}' $DIR/tmp/fisnnet22_sub2.csv | awk '{print NF}') \
    | awk -F, 'BEGIN{OFS=","} $9 > 1  {print $1, $2, $3, $4, $5, $6, $7, $8}' \
    > $DIR/tmp/fisnnet_sub3.csv 

cat <(awk -F, 'NR == 1' $DIR/tmp/fisnnet22_sub2.csv) \
    $DIR/tmp/fisnnet_sub3.csv > $DIR/tmp/fisnnet_sub4.csv

# add IGB id if the revious is the final version
Fseq=$(echo "$(wc -l < $DIR/tmp/fisnnet_sub4.csv) - 1 " | bc)

paste -d"," \
    <(printf "%s\n" IGB_id $(printf "F%1g\n" $(seq 1 $Fseq))) \
    $DIR/tmp/fisnnet_sub4.csv \
    > $DIR/bio_tb/fish_jucar_fishnet2.csv

# remove temporal files
rm $DIR/tmp/fisnnet_sub3.csv $DIR/tmp/fisnnet_sub1.csv $DIR/tmp/fisnnet_sub2.csv $DIR/tmp/fisnnet_sub4.csv




###################### GBIF

# jucar
unzip $DIR/bio_down/0116072-230530130749713.zip -d $DIR/bio_down
gbif=$DIR/bio_down/0116072-230530130749713.csv

#awk -F'\t' '{print NF}' $gbif | sort |  uniq -c  
#awk -F'\t' 'NR == 1 {printf $0}'  $gbif  > borrar.txt  
#printf "%s\n" $(awk -F'\t' 'NR == 1 {printf $0}'  $gbif) | cat -n
awk 'BEGIN{FS="\t"; OFS=","}  {print $1, $6, $7, $8, $10, $22, $23, $33}' \
    $gbif > $DIR/tmp/fish_gbif1.csv

# remove records with empty fields for species names
awk -F, '$5 != ""' $DIR/tmp/fish_gbif1.csv > $DIR/tmp/fish_gbif2.csv

# add IGB id if the revious is the final version
Gseq=$(echo "$(wc -l < $DIR/tmp/fish_gbif2.csv) - 1 " | bc)

paste -d"," \
    <(printf "%s\n" IGB_id $(printf "G%1d\n" $(seq 1 $Gseq))) \
    $DIR/tmp/fish_gbif2.csv \
    > $DIR/bio_tb/fish_jucar_gbif.csv


# remove temporal files
rm $DIR/tmp/fish_gbif1.csv $DIR/tmp/fish_gbif2.csv


##################### idigbio

# initial table created in gbif.R
#scp -r /home/jaime/data/sosw/bio_down/fish_idigbio.txt      sv3:/mnt/shared/sosw/bio_down/

idigbio=$DIR/bio_down/fish_jucar_idigbio.txt

# paste together the useful and fixed columns
paste -d"," \
    <(awk -F"|" 'BEGIN{OFS=","}{print $1, $2}' $idigbio) \
    <(awk -F[,\|] '{print $4}' $idigbio | awk '{print $1, $2}') \
    <(awk -F"|" 'BEGIN{OFS=","}{print $5, $6, $7}' $idigbio) \
    <(awk -F"|" '{print $8}' $idigbio | awk -F"-" '{print $1}') \
    > $DIR/tmp/fish_idigbio.csv


# delete numbers
# delete sp x ? cf.
awk -F, '$3 != "" && $3 !~ / sp/ && $3 !~ /?/ && $3 !~ /cf./ && $3 !~ / x/ && $3 !~ /[+0-9]/'  \
    $DIR/tmp/fish_idigbio.csv > $DIR/tmp/fish_idigbio2.csv

# delete only genus
paste -d"," $DIR/tmp/fish_idigbio2.csv \
    <(awk -F, '{print $3}' $DIR/tmp/fish_idigbio2.csv | awk '{print NF}') \
    | awk -F, 'BEGIN{OFS=","} $8 > 1  {print $1, $2, $3, $4, $5, $6, $7}' \
    > $DIR/tmp/fish_idigbio3.csv

# add header again
cat <(awk -F, 'NR == 1' $DIR/tmp/fish_idigbio2.csv) \
    $DIR/tmp/fish_idigbio3.csv > $DIR/tmp/fish_idigbio4.csv

# primer letra mayuscula
paste -d"," \
<(awk -F, 'BEGIN{OFS=","} {print $1, $2}' $DIR/tmp/fish_idigbio4.csv) \
<(awk -F, '{print $3}' $DIR/tmp/fish_idigbio4.csv \
    | awk '{print toupper(substr($0,1,1)) tolower(substr($0,2)) }') \
<(awk -F, 'BEGIN{OFS=","} {print $4, $5, $6, $7}' $DIR/tmp/fish_idigbio4.csv) \
> $DIR/tmp/fish_idigbio5.csv

# add IGB id if the revious is the final version
Iseq=$(echo "$(wc -l < $DIR/tmp/fish_idigbio5.csv) - 1 " | bc)

paste -d"," \
    <(printf "%s\n" IGB_id $(printf "I%1g\n" $(seq 1 $Iseq))) \
    $DIR/tmp/fish_idigbio5.csv \
    > $DIR/bio_tb/fish_jucar_idigbio.csv

rm $DIR/tmp/fish_idigbio.csv $DIR/tmp/fish_idigbio2.csv \
    $DIR/tmp/fish_idigbio3.csv $DIR/tmp/fish_idigbio4.csv \
    $DIR/tmp/fish_idigbio5.csv

##################### Vertnet

#  check names and order of column names
#awk -F, 'NR == 1 {print $0}'  $VER/vertnet_latest_fishes.csv  > $DIR/borrar.txt
#printf "%s\n" $(sed 's/,/ /g' $DIR/borrar.txt) | cat -n  

cp $DIR/bio_down/vertnet_latest_fishes.csv \
    $DIR/tmp/vernet_fish_filter.csv

# remove fossils
#sed -i '/FossilSpecimen/d' $DIR/tmp/vernet_fish_filter.csv 

# subset the big table
awk 'BEGIN{FPAT = "([^,]*)|(\"[^\"]*\")"; OFS=","} NF == 194 {print $8, $36, $37, $130, $146, $152, $160, $161, $162, $163, $164, $165, $177 }' $DIR/bio_down/vertnet_latest_fishes.csv \
    > $DIR/tmp/vertnet_fish_sub.csv
#awk 'BEGIN{FPAT = "([^,]*)|(\"[^\"]*\")"; OFS=","} NF == 194 {print $8, $36, $37, $130, $146, $152, $160, $161, $162, $163, $164, $165, $177 }' $DIR/tmp/vernet_fish_filter.csv \
#    > $DIR/tmp/vertnet_fish_sub.csv

# filter by leaving only mappable == 1
awk -F, 'NR == 1 || $13 == 1' $DIR/tmp/vertnet_fish_sub.csv \
    > $DIR/tmp/fish_vernet.csv

# identify which point overlap with region of interest and add column (0 and 1)
# to main table
paste -d"," \
    $DIR/tmp/fish_vernet.csv \
    <(printf "%s\n" "inbbox" $(awk -F, 'NR > 1 {print $3, $2}' $DIR/tmp/fish_vernet.csv  | gdallocationinfo -valonly -wgs84 $DIR/msk_jucar.tif | awk '!NF{$0=0}1')) \
    > $DIR/tmp/fish_vernet2.csv

# select only points in the region of interest
awk -F, 'NR == 1 || $14 == 1' $DIR/tmp/fish_vernet2.csv \
    > $DIR/tmp/fish_vernet3.csv

# delete space, sp, ?
awk -F, '$6 != "" && $6 !~ / sp/ && $6 !~ /?/'  \
    $DIR/tmp/fish_vernet3.csv > $DIR/tmp/fish_vernet4.csv

# delete only genus
paste -d"," $DIR/tmp/fish_vernet4.csv \
    <(awk -F, '{print $6}' $DIR/tmp/fish_vernet4.csv | awk '{print NF}') \
    | awk -F, 'BEGIN{OFS=","} $15 > 1  \
    {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' \
    > $DIR/tmp/fish_vernet5.csv

# add header again
cat <(awk -F, 'NR == 1' $DIR/tmp/fish_vernet4.csv) \
    $DIR/tmp/fish_vernet5.csv > $DIR/tmp/fish_vernet6.csv

# add IGB id if the revious is the final version
Vseq=$(echo "$(wc -l < $DIR/tmp/fish_vernet6.csv) - 1 " | bc)

paste -d"," \
    <(printf "%s\n" IGB_id $(printf "V%1g\n" $(seq 1 $Vseq))) \
    $DIR/tmp/fish_vernet6.csv \
    > $DIR/bio_tb/fish_jucar_vernet.csv

# remove temp files
rm $DIR/tmp/vernet_fish_filter.csv \ #$DIR/tmp/vertnet_fish_sub.csv \
    $DIR/tmp/fish_vernet.csv  $DIR/tmp/fish_vernet2.csv \
    $DIR/tmp/fish_vernet3.csv $DIR/tmp/fish_vernet4.csv \
    $DIR/tmp/fish_vernet5.csv $DIR/tmp/fish_vernet6.csv

###############
############### MERGE

cat \
    <(echo "IGB_id,species,long,lat,year") \
    <(awk -F, 'BEGIN{OFS=","} NR > 1 {print $1, $5, $7, $8, $9}' $DIR/bio_tb/fish_jucar_fishnet2.csv) \
    <(awk -F, 'BEGIN{OFS=","} NR > 1 {print $1, $6, $8, $7, $9}' $DIR/bio_tb/fish_jucar_gbif.csv) \
    <(awk -F, 'BEGIN{OFS=","} NR > 1 {print $1, $4, $6, $7, $8}' $DIR/bio_tb/fish_jucar_idigbio.csv) \
    <(awk -F, 'BEGIN{OFS=","} NR > 1 {print $1, $7, $4, $3, $5}' $DIR/bio_tb/fish_jucar_vernet.csv) \
    > $DIR/bio_tb/fish_jucar.csv

################
################

### Create MASK of Júcar

gdal_rasterize \
    -l F162C175_Demarcacion \
    -a OBJECTID \
    -ot Byte \
    -tr 50 50 \
    -a_nodata 0 \
    $DIR/gis/F162_Ambito_de_la_Demarcacion_Hidrografica_del_Jucar/F162C175_Demarcacion.shp  \
    $DIR/gis/msk_jucar_tmp.tif

gdalwarp -t_srs EPSG:4326 -srcnodata 0 -dstnodata 0 \
    $DIR/gis/msk_jucar_tmp.tif $DIR/gis/msk_jucar.tif


paste -d"," \
    $DIR/bio_tb/fish_jucar.csv \
    <(printf "%s\n" "inbbox" $(awk -F, 'NR > 1 {print $3, $4}' $DIR/bio_tb/fish_jucar.csv  | gdallocationinfo -valonly -wgs84 $DIR/gis/msk_jucar.tif | awk '!NF{$0=0}1')) \
    > $DIR/tmp/fish_jucar_bbox.csv



# select only points in the region of interest
awk -F, 'NR == 1 || $6 == 1' $DIR/tmp/fish_jucar_bbox.csv \
    > $DIR/tmp/fish_jucar_inbbox.csv

#awk -F, 'NR > 1 {print $2}' $DIR/tmp/fish_jucar_inbbox.csv | sort | uniq \
#    >  $DIR/tmp/especiesJucar.txt

# select from main table only the ocurrences of the species present in Júcar
# basin, otherwise, delete
awk -F, 'NR==FNR {a[$1]; next} FNR==1 || $2 in a' \
    $DIR/tmp/especiesJucar.txt  $DIR/bio_tb/fish_jucar.csv \
    > $DIR/tmp/allinjucar.csv

#awk -F, 'NR > 1 {print $2}' $DIR/tmp/allinjucar.csv | sort | uniq -c \
#    > $DIR/tmp/corregirNombres.txt

###   Some species are marine --  procedure to delete them

# rasterize a vector hand-made poligon of marine areas
gdal_rasterize \
    -l marine -a id \
    -ot Byte -tr 0.005 -0.005 \
    -a_nodata 0 \
    $DIR/tmp/marine.gpkg \
    $DIR/gis/mks_marine.tif

# create new table with new column with 1 for marine and 0 for freshwater
paste -d"," \
    $DIR/tmp/allinjucar.csv \
    <(printf "%s\n" "marine" $(awk -F, 'NR > 1 {print $3, $4}' $DIR/tmp/allinjucar.csv | gdallocationinfo -valonly -wgs84 $DIR/gis/mks_marine.tif | awk '!NF{$0=0}1')) \
    > $DIR/tmp/fish_marine_tmp.csv


#awk -F, '$6 ==1' $DIR/tmp/fish_marine_tmp.csv | sort | uniq | wc -l 
#awk -F, '$6 ==1 {print $2}' $DIR/tmp/fish_marine_tmp.csv | sort | uniq | wc -l

# create file with names of marine species
awk -F, '$6 ==1 {print $2}' $DIR/tmp/fish_marine_tmp.csv > $DIR/tmp/marine_sppnames.txt

# remove marine species from all species table
awk -F, 'NR==FNR {a[$1]; next} FNR==1 || ! ($2 in a)' \
    $DIR/tmp/marine_sppnames.txt $DIR/tmp/fish_marine_tmp.csv \
    > $DIR/tmp/fish_jucar_freshw.csv


paste -d"," \
    $DIR/tmp/fish_jucar_freshw.csv \
    <(printf "%s\n" "inbbox" $(awk -F, 'NR > 1 {print $3, $4}' $DIR/tmp/fish_jucar_freshw.csv  | gdallocationinfo -valonly -wgs84 $DIR/gis/msk_jucar.tif | awk '!NF{$0=0}1')) \
    > $DIR/tmp/fish_jucar_freshw_bbox.csv

# select only points in the region of interest
awk -F, 'NR == 1 || $7 == 1' $DIR/tmp/fish_jucar_freshw_bbox.csv \
    > $DIR/tmp/fish_jucar_freshw_inbbox.csv

#awk -F, 'NR > 1 {print $2}' $DIR/tmp/fish_jucar_freshw_inbbox.csv | sort | uniq -c \
#    >  $DIR/tmp/especiesJucar.txt
#awk -F, 'NR > 1 {print $2}' $DIR/tmp/fish_jucar_freshw.csv | sort | uniq -c \
#    > $DIR/tmp/corregirNombres.txt


###  add species not included in the previous analysis but in the Jucar papers

#awk -F, 'BEGIN{OFS=","}{print $1, $2, $3, $4, $5}' $DIR/tmp/fish_jucar_freshw_inbbox.csv \
#    > $DIR/tmp/fish_jucar_plus.csv 
#cp $DIR/tmp/fish_jucar_freshw_inbbox.csv $DIR/tmp/fish_jucar_plus.csv 

#awk -F, 'BEGIN{OFS=","} /Luciobarbus guiraonis/  {print $1, $2, $3, $4, $5}' \
#    $DIR/tmp/fish_marine_tmp.csv >> $DIR/tmp/fish_jucar_plus.csv 


#add=(Luciobarbus guiraonis, Cobitis paludica, Squalius pyrenaicus, Anguilla anguilla, Squalius valentinus )

# get from $DIR/tmp/fish_marine_tmp.csv

awk -F, '/Luciobarbus guiraonis/ || /Cobitis paludica/ || \
        /Squalius pyrenaicus/ || /Anguilla anguilla/ || /Squalius valentinus/' \
        $DIR/tmp/fish_marine_tmp.csv > $DIR/tmp/fish_add.csv

# then run mask to jucar

paste -d"," \
    $DIR/tmp/fish_add.csv \
    <(awk -F, 'NR > 1 {print $3, $4}' $DIR/tmp/fish_add.csv  | gdallocationinfo -valonly -wgs84 $DIR/gis/msk_jucar.tif | awk '!NF{$0=0}1') \
    > $DIR/tmp/fish_add_bbox.csv

# select only points in the region of interest
awk -F, '$7 == 1' $DIR/tmp/fish_add_bbox.csv \
    > $DIR/tmp/fish_add_inbbox.csv

# merge last table with previous table
cat $DIR/tmp/fish_jucar_freshw_inbbox.csv $DIR/tmp/fish_add_inbbox.csv \
    > $DIR/tmp/test.csv

awk -F, 'NR > 1 {print $2}' $DIR/tmp/test.csv | sort | uniq -c


####   based on the new roi crop the biological data only with the species of interests

## make the manual made polygon a raster
gdal_rasterize -l roi_modelling -a roi -tr 0.0083 -0.0083 -a_nodata 0 -ot Byte \
    -of GTiff gis/roi_modelling.gpkg gis/roi_modelling.tif


# create new table with new column with 1 for marine and 0 for freshwater
paste -d"," \
    $DIR/bio_tb/fish_jucar.csv \
    <(printf "%s\n" "marine" $(awk -F, 'NR > 1 {print $3, $4}' $DIR/bio_tb/fish_jucar.csv | gdallocationinfo -valonly -wgs84 $DIR/gis/mks_marine.tif | awk '!NF{$0=0}1')) \
    > $DIR/tmp/fish_jucar_marine_tmp.csv

# select only terrestrial points 
awk -F, 'NR == 1 || $6 == 0' $DIR/tmp/fish_jucar_marine_tmp.csv \
    > $DIR/tmp/fish_jucar_terrestrial.csv

# select only points for the roi_modelling area
paste -d"," \
    $DIR/tmp/fish_jucar_terrestrial.csv \
    <(printf "%s\n" "roi" $(awk -F, 'NR > 1 {print $3, $4}' \
    $DIR/tmp/fish_jucar_terrestrial.csv | gdallocationinfo -valonly -wgs84 \
    $DIR/gis/roi_modelling.tif | awk '!NF{$0=0}1')) \
    | awk -F, 'NR == 1 || $7 == 1' > $DIR/tmp/fish_roi_tmp.csv

# subset data based on species name

awk -F, 'NR == 1 || /Scomberomorus maculatus/ || /Sciaena umbrina/ || /Pylodictis olivaris/ || /Pterygoplichthys pardalis/ || /Poecilia reticulata/ || /Luciobarbus sclateri/ || /Luciobarbus bocagei/ || /Lepisosteus oculatus/ || /Ictiobus bubalus/ || /Exoglossum laurae/ || /Chondrostoma toxostoma/ || /Carpiodes cyprinus/ || /Campogramma glaycos/ || /Barbus barbus/ || /Luciobarbus guiraonis/ || /Cobitis paludica/ || /Squalius pyrenaicus/ || /Anguilla anguilla/ || /Squalius valentinus/ || /Achondrostoma arcasii/ || /Barbus haasi/ || /Gobio gobio/ || /Liparis atlanticus/ || /Luciobarbus graellsii/ || /Misgurnus anguillicaud/ || /Parachondrostoma arrigonis/ || /Parachondrostoma miegii/ || /Parachondrostoma turiensis/ || /Pseudochondrostoma polylepis/ || /Salmo trutta/ || /Sander lucioperca/ || /Squalius alburnoides/ || /Tinca tinca/' $DIR/tmp/fish_roi_tmp.csv > $DIR/tmp/fish_roi_modelling.csv

###  crop roi data for only within Jucar basin 

paste -d"," \
    $DIR/tmp/fish_roi_modelling.csv \
    <(printf "%s\n" "jucar" $(awk -F, 'NR > 1 {print $3, $4}' $DIR/tmp/fish_roi_modelling.csv  | gdallocationinfo -valonly -wgs84 $DIR/gis/msk_jucar.tif | awk '!NF{$0=0}1')) \
    > $DIR/tmp/fish_jucar_tmp.csv

# select only points in the region of interest
awk -F, 'NR == 1 || $8 == 1' $DIR/tmp/fish_jucar_tmp.csv \
    > $DIR/tmp/fish_jucar_freshw.csv

paste -d" " \
    <(awk -F, 'NR > 1 {print $2}' $DIR/tmp/fish_jucar_freshw.csv | sort | uniq -c) \
    <(awk -F, 'NR > 1 {print $2}' $DIR/tmp/fish_roi_modelling.csv | sort | uniq -c) \
    > $DIR/tmp/table.csv




###############################################################################

##  How many species
##  How many occurrences per species  (inbbox and outside)
##  How many genera, families, order?

##   Anaylsis per species:
# 1 remove duplicates
