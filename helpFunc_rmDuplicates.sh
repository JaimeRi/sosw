#! /bin/bash


# chunk to delete records located in the same subcatchment

#export spp=$1
#export spp="Huso huso"
#export nm=$(echo "$spp" | awk 'BEGIN{OFS="_"}{print $1,$2}')
#export nm=$1
#export spp=$(echo "$nm" | awk -F_ '{print $1, $2}')

#export table=$2
#export table=$DIR/tmp/fish_danube_fil2.csv

export table=$1

rmDuplic(){

export spp=$1
export nm=$(echo "$spp" | awk 'BEGIN{OFS="_"}{print $1,$2}')

awk -F, -v SPP="$spp" '$2 == SPP {print $1, $3, $4}' $table \
    > $DIR/tmp/spp/fd_${nm}.txt

paste -d" " \
    $DIR/tmp/spp/fd_${nm}.txt \
    <(awk '{print $2, $3}' $DIR/tmp/spp/fd_${nm}.txt | gdallocationinfo -valonly \
    -wgs84 $DIR/tmp/danube_subcatchments.tif) \
    > $DIR/tmp/spp/fd_${nm}_sid.txt

awk '{print $4}' $DIR/tmp/spp/fd_${nm}_sid.txt | sort | uniq -c \
    > $DIR/tmp/spp/fd_${nm}_sid_summ.txt

for i in $(seq 1 $(wc -l < $DIR/tmp/spp/fd_${nm}_sid_summ.txt))
do
    num=$(awk -v row="$i" 'NR == row {print $1}' $DIR/tmp/spp/fd_${nm}_sid_summ.txt)
    cid=$(awk -v row="$i" 'NR == row {print $2}' $DIR/tmp/spp/fd_${nm}_sid_summ.txt)

    if [[ $num -eq 1 ]]
    then
    awk -v cid="$cid" '$4 == cid {print $0}' $DIR/tmp/spp/fd_${nm}_sid.txt \
        > $DIR/tmp/spp/tmp_${nm}_${i}.txt
    else
    #  here add code to select the record with the coordinates with higher
    # accuracy, i.e., more decimal points
    # num=0.000001
    # decimals=${num#*.}              #Removes the integer part and the dot (=000001)
    # decimalPointsCount=${#decimals} #Counts the length of resulting string (=6)
    awk -v cid="$cid" '$4 == cid {print $0}' $DIR/tmp/spp/fd_${nm}_sid.txt \
        | awk 'NR == 1' > $DIR/tmp/spp/tmp_${nm}_${i}.txt

    fi
done

cat $(find $DIR/tmp/spp -name "tmp_${nm}_*.txt") | awk '{print $1}' \
    > $DIR/tmp/spp/retain_indx_${nm}.txt 

rm $DIR/tmp/spp/tmp_${nm}_*.txt $DIR/tmp/spp/fd_${nm}*.txt

}

export -f rmDuplic
awk -F, 'NR > 1 {print $2}' $table | sort | uniq  \
    |  parallel -j 50 rmDuplic 

# join all tables together
cat $(find $DIR/tmp/spp -name "retain_indx*.txt") > $DIR/tmp/retain_indx.txt

# filter database to remove duplicates 
awk -F, 'NR==FNR {a[$1]; next} FNR==1 ||  $1 in a' \
    $DIR/tmp/retain_indx.txt  $table \
    > $DIR/tmp/fish_danube_fil3new.csv


