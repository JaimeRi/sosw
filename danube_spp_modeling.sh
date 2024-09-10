export DIR=/mnt/shared/danube


awk -F, 'NR > 1 {print $1}' $DIR/species_list.csv > $DIR/spp_danube_prio.txt

SPP=$(awk 'NR == 1' $DIR/spp_danube_prio.txt)
SPPN=$(echo $SPP | tr " " "_")

awk -F, -v SP="$SPP" 'BEGIN{OFS=",";} NR == 1 ||  $2 == SP  {print $2, $3, $4}' $DIR/fish_danube.csv \
    > $DIR/spp_$SPPN.txt

SPPN=Zingel_streber
awk -F, 'BEGIN{OFS=",";}   /"Zingel streber"/{print $2, $3, $4}' $DIR/danube_records_final_20-03-2024.csv | tr -d '"' > $DIR/spp_$SPPN.txt


bash $DIR/dev_create_model_table.sh \
    $DIR/spp_$SPPN.txt \
    $DIR/out/danube_predictTB.csv  \
    /mnt/shared/sosw/tmp/danube_subcatchments.tif \
    10000  \
    $DIR/tmp \
    $DIR/model_table.csv 


### to split the data given a number

awk 'NR==1{
        header=$0; 
        count=1; 
        print header > "predTB_" count; 
        next 
     } 

     !( (NR-1) % 500000){
        count++; 
        print header > "predTB_" count;
     } 
     {
        print $0 > "predTB_" count
     }' $DIR/out/danube_predictTB.csv
