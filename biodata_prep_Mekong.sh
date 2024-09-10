export DIR=/mnt/shared/sosw
gdal_rasterize -l vietnam -a fid -tr 0.0833 0.083333 -a_nodata 0 -ot Int32 -of GTiff /home/jaime/data/sosw/mekongGIS/vietnam.gpkg /home/jaime/data/sosw/mekongGIS/vietnam.tif
###   make masks of ROI for the study cases


### VIETNAM
printf "gid,WKT\n1,\"POLYGON((101.921002148438 23.4160971210938,102.505014916992 12.2723738056641,102.732922338867 12.1536720234375,102.87061640625 11.8569175678711,103.167370861816 11.7073533222657,103.150752612305 11.3749883320313,103.173305950928 11.1910005695801,103.176867004395 10.9702152546387,103.311000018311 10.953597005127,103.410709515381 11.1043482685547,103.535346386719 11.1921875874024,103.652861151123 11.082981947754,103.752570648193 10.8681317219239,103.720521166992 10.7636741535645,103.705757632828 10.7323665585023,103.57043760109 10.6397791683655,103.668366571426 10.5471917782288,103.887371359634 10.6308765346986,104.127742468643 10.6255349544984,104.307575668716 10.5934854732972,104.412626745987 10.4759707088929,104.569313098526 10.3851638454896,104.626289953995 10.3032596157533,104.752707352066 10.2783322414857,104.852416849136 10.2409411800843,104.936101605606 10.1625980038148,105.091007431412 10.0806937740785,105.153325867081 9.98098427700816,105.14798428688 9.84210319180309,105.073202164078 9.81539529080211,105.04293320961 9.84388371853649,104.907613177872 9.71212474026501,104.905832651138 9.62665945706189,104.852416849136 8.92157087063615,104.969931613541 8.78625083889787,104.966370560074 8.6473697536928,105.308231692886 8.91088771023576,105.717752841568 9.28479832424941,106.095224509048 9.45216783718885,105.689264413834 10.1074016750794,106.800313095474 10.7519523525696,107.281055313492 10.605949160431,107.398570077896 10.5204838772279,108.384981888199 11.1579124477845,108.851479892349 11.410747243927,109.097192581558 12.0837863491516,109.15060838356 13.0773202663879,109.043776779556 14.3592995144347,108.691232486343 15.2353186672667,108.103658664322 15.9617735744932,107.334471115494 16.5920800381162,107.152857388687 16.8765191837766,107.014644001007 17.0988624596097,106.589988375091 17.4514067528226,106.413716228485 17.7298366207578,106.435750246811 17.9381582485654,106.453778079986 17.956186081741,106.391682210159 18.0403159698941,106.188702162552 18.1818678451993,105.87087814064 18.4462760651089,105.683922833633 18.7186966553189,105.550383328628 18.9857756653286,105.670568883133 19.1994388733364,105.764046536636 19.5806941601253,105.921623152542 20.0420731499172,106.14062794075 20.1221968529201,106.378328259659 20.2637487282253,106.525221715164 20.4373500847316,106.597333047867 20.8533256428218,106.594662257767 21.0830135914302,106.832362576675 21.0963675419307,107.23565188179 21.0956998444057,107.323787955093 21.3173754227138,107.459998250198 21.3974991257167,107.679003038406 21.4555888103938,107.881983086014 21.5731035747981,108.05291365242 21.679935178802,108.261235280227 21.7734128323054,108.426824266434 22.0591873730158,107.353166646194 23.0006408833002,106.044479497147 23.4386504597162,101.921002148438 23.4160971210938))\"\n" \
    > $DIR/mekongGIS/roi_vietnam.csv

ogr2ogr -f "GPKG" $DIR/mekongGIS/roi_vietnam.gpkg -dialect sqlite -sql "SELECT gid, ST_GeomFromText(WKT,4326) FROM roi_vietnam" -nln roi -a_srs EPSG:4326 $DIR/mekongGIS/roi_vietnam.csv

gdal_rasterize -a gid -l roi -tr 0.083 0.083 $DIR/mekongGIS/roi_vietnam.gpkg $DIR/mekongGIS/roi_vietnam.tif -ot Byte 
gdal_rasterize -a fid -l roi_vietnam -tr 0.083 0.083 /home/jaime/data/sosw/mekongGIS/roi_vietnam.gpkg /home/jaime/data/sosw/mekongGIS/roi_vietnam.tif -ot Byte 


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
gbif=$DIR/mekongGIS/gbif_mekong.csv

#awk -F'\t' '{print NF}' $gbif | sort |  uniq -c  
#awk -F'\t' 'NR == 1 {printf $0}'  $gbif  > borrar.txt  
#printf "%s\n" $(awk -F'\t' 'NR == 1 {printf $0}'  $gbif) | cat -n
awk 'BEGIN{FS="\t"; OFS=","}  {print $1, $6, $7, $8, $10, $22, $23, $33}' \
    $gbif > $DIR/tmp/fishm_gbif1.csv

# remove records with empty fields for species names
awk -F, '$5 != ""' $DIR/tmp/fishm_gbif1.csv > $DIR/tmp/fishm_gbif2.csv

# add IGB id if the revious is the final version
Gseq=$(echo "$(wc -l < $DIR/tmp/fishm_gbif2.csv) - 1 " | bc)

paste -d"," \
    <(printf "%s\n" IGB_id $(printf "G%1d\n" $(seq 1 $Gseq))) \
    $DIR/tmp/fishm_gbif2.csv \
    > $DIR/bio_tb/fish_mekong_gbif.csv


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


