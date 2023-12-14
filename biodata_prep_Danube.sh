export DIR=/mnt/shared/sosw

###   make masks of ROI for the study cases

### DANUBE
printf "gid,WKT\n1,\"POLYGON((8 51, 30 51, 30 42, 8 42, 8 51))\"\n" \
    > $DIR/roi_danube.csv

ogr2ogr -f "GPKG" $DIR/roi_danube.gpkg -dialect sqlite -sql "SELECT gid, ST_GeomFromText(WKT,4326) FROM roi_danube" -nln roi -a_srs EPSG:4326 $DIR/roi_danube.csv

gdal_rasterize -a gid -l roi -tr 0.083 0.083 $DIR/roi_danube.gpkg $DIR/msk_danube.tif -ot Byte 


###############################################################################
#######   Fish  #######################
###############################################################################

######################  Fishnet2

# CHeck on http://www.fishnet2.net/search.aspx?c=OS using the polygon below
# POLYGON((8 51, 30 51, 30 42, 8 42, 8 51))

export fishnet=$DIR/bio_down/fishnet2_danube.txt

# check number of columns in table and if they are consistent given the field separator
#awk -F"\t" '{print NF}' $fishnet | sort | uniq -c

# check name of headers and their position
#printf "%s\n" $(awk -F"\t" 'NR == 1'  $fishnet | sed 's/,/ /g') | cat -n  

paste -d"," \
    <( awk -F"\t"  'BEGIN{OFS=","}{print $1, $2, $3}' ${fishnet} ) \
    <(awk -F"\t" '{print $5}' $fishnet | awk -F[,\(\)] '{print $1, $2}' | awk '{print $1, $2}') \
    <(awk -F"\t"  'BEGIN{OFS=","}{print $6, $10, $9, $21}' $fishnet) \
    > $DIR/tmp/fisnnet_sub1.csv

# Check errors in names
# awk -F, '{print $4}' $DIR/tmp/fisnnet_sub1.csv  | sort | uniq -c | more

# remove records with empty species names
# remove records where species names are given as sp or sp.
# remove records with species names as Unkown
# remove records where species names are given as x
# remove strange names

### DANUBE 
awk -F, '$4 != " " && $4 !~ / sp/ && $4 !~ /?/ && $4 !~ /Unknown/ && $4 !~ / x/ && $4 != "LeuciscusÊleuciscus LeuciscusÊleuciscus" && $4 != "RutilusÊrutilus RutilusÊrutilus" && $4 != "GobioÊgobio GobioÊgobio" && $4  !~ /Unidentified/' \
    $DIR/tmp/fisnnet_sub1.csv > $DIR/tmp/fisnnet_sub2.csv


# removes records where species names only have the genus
paste -d"," $DIR/tmp/fisnnet_sub2.csv \
    <(awk -F, '{print $4}' $DIR/tmp/fisnnet_sub2.csv | awk '{print NF}') \
    | awk -F, 'BEGIN{OFS=","} $9 > 1  {print $1, $2, $3, $4, $5, $6, $7, $8}' \
    > $DIR/tmp/fisnnet_sub3.csv 

cat <(awk -F, 'NR == 1' $DIR/tmp/fisnnet_sub2.csv) \
    $DIR/tmp/fisnnet_sub3.csv > $DIR/tmp/fisnnet_sub4.csv

# add IGB id if the previous is the final version
Fseq=$(echo "$(wc -l < $DIR/tmp/fisnnet_sub4.csv) - 1 " | bc)

paste -d"," \
    <(printf "%s\n" IGB_id $(printf "F%1g\n" $(seq 1 $Fseq))) \
    $DIR/tmp/fisnnet_sub4.csv \
    > $DIR/bio_tb/fish_danube_fishnet2.csv

# remove temporal files
rm $DIR/tmp/fisnnet_sub3.csv $DIR/tmp/fisnnet_sub1.csv $DIR/tmp/fisnnet_sub2.csv $DIR/tmp/fisnnet_sub4.csv




###################### GBIF

#Danube 
unzip $DIR/bio_down/0028469-231002084531237.zip -d $DIR/bio_down
gbif=$DIR/bio_down/0028469-231002084531237.csv

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
    > $DIR/bio_tb/fish_danube_gbif.csv


# remove temporal files
rm $DIR/tmp/fish_gbif1.csv $DIR/tmp/fish_gbif2.csv
rm $DIR/bio_down/0028469-231002084531237.csv

##################### idigbio

# initial table created in gbif.R
#scp -r /home/jaime/data/sosw/bio_down/fish_idigbio.txt      sv3:/mnt/shared/sosw/bio_down/

# paste the individual tables and remove the header (except the first line) from the global table
cat $(find $DIR/tmp -name "idigbio_*.txt") > $DIR/tmp/idigbio_danube_tmp.txt
awk 'NR == 1 || !/uuid/' $DIR/tmp/idigbio_danube_tmp.txt > $DIR/bio_down/fish_danube_idigbio.txt
rm $DIR/tmp/idigbio_danube_tmp.txt

idigbio=$DIR/bio_down/fish_danube_idigbio.txt

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
    > $DIR/bio_tb/fish_danube_idigbio.csv

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
    <(printf "%s\n" "inbbox" $(awk -F, 'NR > 1 {print $3, $2}' $DIR/tmp/fish_vernet.csv  | gdallocationinfo -valonly -wgs84 $DIR/msk_danube.tif | awk '!NF{$0=0}1')) \
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
    > $DIR/bio_tb/fish_danube_vernet.csv

# remove temp files
rm $DIR/tmp/vernet_fish_filter.csv \ #$DIR/tmp/vertnet_fish_sub.csv \
    $DIR/tmp/fish_vernet.csv  $DIR/tmp/fish_vernet2.csv \
    $DIR/tmp/fish_vernet3.csv $DIR/tmp/fish_vernet4.csv \
    $DIR/tmp/fish_vernet5.csv $DIR/tmp/fish_vernet6.csv

###############
############### MERGE

cat \
    <(echo "IGB_id,species,long,lat,year") \
    <(awk -F, 'BEGIN{OFS=","} NR > 1 {print $1, $5, $7, $8, $9}' $DIR/bio_tb/fish_danube_fishnet2.csv) \
    <(awk -F, 'BEGIN{OFS=","} NR > 1 {print $1, $6, $8, $7, $9}' $DIR/bio_tb/fish_danube_gbif.csv) \
    <(awk -F, 'BEGIN{OFS=","} NR > 1 {print $1, $4, $6, $7, $8}' $DIR/bio_tb/fish_danube_idigbio.csv) \
    <(awk -F, 'BEGIN{OFS=","} NR > 1 {print $1, $7, $4, $3, $5}' $DIR/bio_tb/fish_danube_vernet.csv) \
    > $DIR/tmp/fish_danube_all.csv

################
################

### FILTER:

## get all the names of the species located within the Danube basin only AND
## filter the table by extracting the records of those species only. Afterwards
## the table could be filtered by using the mask of the continental area


### Create MASK of Danube
export GLBASINS=/mnt/shared/hydrography90m_v1_2022_all_data_OLD/lbasin_tiles_final20d_ovr/all_lbasin.tif

gdalwarp  -t_srs EPSG:4326  -te 8 42 30 51 $GLBASINS $DIR/tmp/basins_sub.tif \
    -co COMPRESS=LZW -co ZLEVEL=9

# danube 1291835   # small island 1200647
grass78 -f -text --tmp-location -c $DIR/tmp/basins_sub.tif <<'EOF'
r.external input=$DIR/tmp/basins_sub.tif output=basins
r.mapcalc "danube = if(basins == 1291835 || basins == 1200647, 1, null())"
g.region zoom=danube
r.out.gdal --o -f -c -m  createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=Byte \
    format=GTiff nodata=0 input=danube output=$DIR/msk_danube.tif
EOF

# identify names of species within danube basin only
paste -d"," \
    $DIR/tmp/fish_danube_all.csv \
    <(printf "%s\n" "indanube" $(awk -F, 'NR > 1 {print $3, $4}' $DIR/tmp/fish_danube_all.csv | gdallocationinfo -valonly -wgs84 $DIR/msk_danube.tif | awk '!NF{$0=0}1')) \
    |  awk -F, '$6 == 1 {print $2}' | sort | uniq  > $DIR/tmp/names_danube.txt


# select from main table only the ocurrences of the species present in Danube
# basin, otherwise, delete
awk -F, 'NR==FNR {a[$1]; next} FNR==1 || $2 in a' \
    $DIR/tmp/names_danube.txt  $DIR/tmp/fish_danube_all.csv \
    > $DIR/tmp/spp_in_danube.csv

# Now mask for points in continental area

paste -d"," \
    $DIR/tmp/spp_in_danube.csv \
    <(printf "%s\n" "inbbox" $(awk -F, 'NR > 1 {print $3, $4}' $DIR/tmp/spp_in_danube.csv | gdallocationinfo -valonly -wgs84 $DIR/tmp/basins_sub.tif | awk '!NF{$0=0}1')) \
    | awk -F, 'BEGIN{OFS=","} NR == 1 || $6 > 1 {print $1,$2,$3,$4,$5}' \
    > $DIR/tmp/fish_danube_inroi.csv

## identify species with less than 5 ocurrences and delete



## Check for errors in names

# create list with names of species with more ten or more records
awk -F, 'NR > 1 {print $2}' $DIR/tmp/fish_danube_inroi.csv | sort | uniq -c \
 | awk '$1 > 9 {print substr($0, index($0,$2)) }' > $DIR/tmp/spp_remove.txt

# and remove those names from main table
awk -F, 'NR==FNR {a[$1]; next} FNR==1 ||  $2 in a' \
    $DIR/tmp/spp_remove.txt  $DIR/tmp/fish_danube_inroi.csv \
    > $DIR/tmp/fish_danube_fil1.csv

# Visually check species names and correct
awk '{gsub("truttafario","trutta"); \
    gsub("cernuus","cernua"); \
    gsub("encrasicholus","encrasicolus"); \
    gsub("valdykovi","vladykovi");print}' \
    $DIR/tmp/fish_danube_fil1.csv > $DIR/tmp/fish_danube_fil2.csv 

##   DELETE DUPLICATES
# run function to delete duplicates. ocurrences are duplicates if the are in the
# same subcatchment
time bash $DIR/helpFunc_rmDuplicates.sh $DIR/tmp/fish_danube_fil2.csv 


# create list with names of species with ten or more records
awk -F, 'NR > 1 {print $2}' $DIR/tmp/fish_danube_fil3.csv | sort | uniq -c \
 | awk '$1 > 9 {print substr($0, index($0,$2)) }' > $DIR/tmp/spp_keep.txt

# and keep those species from main table
awk -F, 'NR==FNR {a[$1]; next} FNR==1 ||  $2 in a' \
    $DIR/tmp/spp_keep.txt  $DIR/tmp/fish_danube_fil3.csv \
    > $DIR/tmp/fish_danube_fil4.csv

rm $DIR/tmp/spp_keep.txt

cp $DIR/tmp/fish_danube_fil4.csv $DIR/sppTB
mv $DIR/sppTB/fish_danube_fil4.csv $DIR/sppTB/fish_danube.csv



### remove species not in the danube

# plant
# Woodsia ilvensis

# Marine
#Trachinus draco
#Scorpaena porcus
#Salvelinus malma





###############################################################################
###############################################################################
###############################################################################
###############################################################################


