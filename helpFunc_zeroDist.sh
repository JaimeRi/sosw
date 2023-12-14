#! /bin/bash

# function to delete duplicates when records have the same coordinates

# Parameters:
# spp= name of species, e.g., Anguilla anguilla

zeroDist(){

## Check for duplicated (per species)

export spp=$1
export nm=$(echo $spp | awk 'BEGIN{OFS="_"}{print $1,$2}')

awk -F, -v SPP="$spp" '$2 == SPP {print $1, $3, $4}' $DIR/tmp/fish_danube_fil2.csv \
    > $DIR/tmp/spp/fd_${nm}.txt

export tabla=$DIR/tmp/spp/fd_${nm}.txt

R --vanilla --no-readline -q  << "EOF"

library(sp)

DIR = Sys.getenv(c("DIR"))
nm = Sys.getenv(c("nm"))
tabla = Sys.getenv(c("tabla"))

tb = read.table(tabla, sep =" ")

coordinates(tb) = c("V2", "V3")
zd = zerodist(tb)
tb2 <- tb[-zd[,2], ]
out = paste0("/mnt/shared/sosw/tmp/spp/indx_", nm, ".txt")
write.table(as.data.frame(tb2)[,1], out, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)

EOF

rm $tabla

}

export -f zeroDist
#awk -F, 'NR >1 {print $2}' $DIR/tmp/fish_danube_fil2.csv | sort | uniq \
#    | parallel -j 30 zeroDist
