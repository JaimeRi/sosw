cd ~/proyectos/sosw/
 head DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01.txt 
 awk -F, '{print NF}' DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01.csv | sort | uniq -c
 awk  '{print NF}' DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01.txt | sort | uniq -c
 awk  'NF ==  7 {print $0}' DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01.txt | sort | uniq -c
 awk -F, 'NF == 6 {print $1, $2}' DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01.csv | sort | uniq -c
 awk -F, 'NF == 6 {print $1, $2, $3}' DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01.csv | sort | uniq -c
 awk -F, 'NF == 6 {print $0}' DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01.csv | sort | uniq -c

 awk -F"|" '{print NF}' DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01_jaime.csv | sort |
uniq -c
 
 awk -F"|" 'NF == 6 {print $0, NR}' DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01_jaime.csv | sort | uniq -c
 
 awk -F"|" 'NR == 822 {print $0, NR}' DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01_jaime.csv

 wc -l < DB_FISH_RICHNESS_JRBD_Martinez-Capel_ver01_jaime.csv
 wc -l < db_fish_jucar.csv

 awk -F, '{print NF}' db_fish_jucar.csv | sort | uniq -c
 awk -F, 'NF == 6 {print $0, NR}' db_fish_jucar.csv | sort | uniq -c

#### select records with coordinates

 awk -F, '$11 != "NA"' db_fish_jucar.csv > db_fish_jucar_coord.csv

 ogr2ogr -f "gpkg" -nln db_jucar -a_srs EPSG:32630 db_jucar_all.gpkg db_fish_jucar_coord.csv -oo x_possible_names=X -oo y_possible_names=Y -oo autodetect_type=yes

 ogr2ogr -f "gpkg" -nln db_jucar -a_srs EPSG:32630 db_jucar_all.gpkg db_fish_jucar_coord.csv -oo x_possible_names=X -oo y_possible_names=Y -oo autodetect_type=yes

 ogr2ogr -s_srs EPSG:32630 -t_srs EPSG:4326 db_jucar_ll.gpkg db_jucar_all.gpkg

### Use the OGR SQL http://www.gdal.org/ogr_sql.html and add X and Y columns with

ogrinfo db_jucar_ll.gpkg -sql "alter table db_jucar add column Xll double"
ogrinfo db_jucar_ll.gpkg -sql "alter table db_jucar add column Yll double"

####Then switch to SQLite dialect http://www.gdal.org/ogr_sql_sqlite.html

ogrinfo db_jucar_ll.gpkg -dialect SQLite -sql "update db_jucar set Xll=ST_X(geom)"
ogrinfo db_jucar_ll.gpkg -dialect SQLite -sql "update db_jucar set Yll=ST_Y(geom)"

ogr2ogr -f "csv" db_jucar_ll.csv db_jucar_ll.gpkg

#awk -F, -v occname="Anguilla.anguilla" -v lon="Xll" -v lat="Yll" \
#    'NR == 1 { for (i=1; i<=NF; i++) {f[$i] = i} } \
#    BEGIN{OFS=",";} $(f[occname])==1 {print $(f[lon]),$(f[lat])}' db_jucar_ll.csv


for spp in $(cat spp_list_jucar_martinez.txt)
do

awk -F, -v occname="${spp}"  -v lon="Xll" -v lat="Yll" \
    'NR == 1 { for (i=1; i<=NF; i++) {f[$i] = i} } \
    BEGIN{OFS=",";} gsub(/"/, "", $(f[occname])) \
    {print occname, $(f[occname]), $(f[lon]), $(f[lat])}' db_jucar_ll.csv \
        | awk -F, '$2 > 0' \
        | awk -F, 'BEGIN{OFS=",";} gsub("\\.", " ", $1) {print($1,$3,$4)}' \
        > spp_tmp/${spp}_jucar.csv
done

cat spp_tmp/* > all_spp_jucar_martinez.csv

# add IGB id if the previous is the final version
Fseq=$(echo "$(wc -l < all_spp_jucar_martinez.csv)")

paste -d"," \
    <(printf "MAR%1g\n" $(seq 1 $Fseq)) \
    all_spp_jucar_martinez.csv  \
    > jdbid.csv

######### remove duplicates

export table="jdbid.csv"

rmDuplic(){

export spp=$1
export nm=$(echo "$spp" | awk 'BEGIN{OFS="_"}{print $1,$2}')

awk -F, -v SPP="$spp" '$2 == SPP {print $1, $3, $4}' $table \
    > tmp/fd_${nm}.txt

paste -d" " \
    tmp/fd_${nm}.txt \
    <(awk '{print $2, $3}' tmp/fd_${nm}.txt | gdallocationinfo -valonly \
    -wgs84 ~/data/sosw/env_jucar/final_layers/sub_catchment.tif) \
    > tmp/fd_${nm}_sid.txt

awk '{print $4}' tmp/fd_${nm}_sid.txt | sort | uniq -c \
    > tmp/fd_${nm}_sid_summ.txt

for i in $(seq 1 $(wc -l < tmp/fd_${nm}_sid_summ.txt))
do
    num=$(awk -v row="$i" 'NR == row {print $1}' tmp/fd_${nm}_sid_summ.txt)
    cid=$(awk -v row="$i" 'NR == row {print $2}' tmp/fd_${nm}_sid_summ.txt)

    if [[ $num -eq 1 ]]
    then
    awk -v cid="$cid" '$4 == cid {print $0}' tmp/fd_${nm}_sid.txt \
        > tmp/tmp_${nm}_${i}.txt
    else
    #  here add code to select the record with the coordinates with higher
    # accuracy, i.e., more decimal points
    # num=0.000001
    # decimals=${num#*.}              #Removes the integer part and the dot (=000001)
    # decimalPointsCount=${#decimals} #Counts the length of resulting string (=6)
    awk -v cid="$cid" '$4 == cid {print $0}' tmp/fd_${nm}_sid.txt \
        | awk 'NR == 1' > tmp/tmp_${nm}_${i}.txt

    fi
done

cat $(find tmp/ -name "tmp_${nm}_*.txt") | awk '{print $1}' \
    > tmp/retain_indx_${nm}.txt 

rm tmp/tmp_${nm}_*.txt tmp/fd_${nm}*.txt

}

export -f rmDuplic
awk -F, '{print $2}' $table | sort | uniq  \
    |  parallel -j 2 rmDuplic 

# join all tables together
cat $(find tmp/ -name "retain_indx*.txt") > tmp/retain_indx.txt

# filter database to remove duplicates 
awk -F, 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
    tmp/retain_indx.txt  $table \
    > fish_jucar_mart_fn.csv


