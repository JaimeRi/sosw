library(rgbif)
setwd("/home/jaime/data/sosw")
setwd("/mnt/shared/sosw")
# Rhine  POLYGON((-0.55 55.01, 15.35 55.01, 15.35 43.92, -0.55 43.92, -0.55 55.01))

#1 -0.55 55.01
#2 -0.55 43.92
#3 15.35 43.92
#4 15.35 55.01

# Júcar  POLYGON((-9.86 46.4, 6.29 46.4, 6.29 34.74, -9.86 34.74, -9.86 46.4))
# Danube POLYGON((3.79 54.31, 40.14 54.31, 40.14 40.14, 3.79 40.14, 3.79 54.31))
# Danube POLYGON((8 51, 30 51, 30 42, 8 42, 8 51))
# Mekong POLYGON((96.89 25.04, 110.76 25.04, 110.76 0.84, 96.89 0.84, 96.89 25.04))


####  Download data for fish using the list of families

Fam = read.table("~/Nextcloud/ERC/FamilyFishes.txt")
Fam = read.table("FamilyFishes.txt")

## Find out the key of each family and populate the table
Fam$Key = NA

for (i in 1:dim(Fam)[1]){
  Fam$Key[i] = name_suggest(q = Fam[i,1], rank = "family")[[1]]$key[1]
}
# save(Fam, file="FishFamilyKeys.RData")
# load(file="FishFamilyKeys.RData")

# 545 Families of fishes
KEY = Fam[,2]

# prepare query to download
occ_download(
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"),
  pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
  pred_within("POLYGON((8 51, 8 42, 30 42, 30 51, 8 51))"),
  pred_in("taxonKey", KEY),
  format = "SIMPLE_CSV",
  user="jaime_ric", pwd ="2ruZNXfWaP6L.58", email="marquez@igb-berlin.de" 
)

# DANUBE = GBIF Occurrence Download https://doi.org/10.15468/dl.qc7e3s Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-10-26
# RHINE = GBIF Occurrence Download https://doi.org/10.15468/dl.r695p8 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-06-06
# JUCAR = GBIF Occurrence Download https://doi.org/10.15468/dl.q53db3 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-07-26

# Check status
occ_download_wait('0028469-231002084531237') # danube

# Once the request is finished, run to download
d <- occ_download_get('0028469-231002084531237') %>%
  occ_download_import()

system("mv 0028469-231002084531237.zip bio_down")
#Warning message:
#  In data.table::fread(targetpath, data.table = FALSE, fill = fill,  :
#                         Found and resolved improper quoting out-of-sample. First healed line 1526448: <<1265784867	2f8a1bdb-1ef9-418b-abc7-a5535746f316	http://id.snsb.info/zsm/collection_zsm/488684/527893/379481	Animalia	Chordata		Perciformes	Blenniidae	Salaria	Salaria pavo		SPECIES	Salaria pavo (Risso, 1810)	Salaria pavo (Risso, 1810)		FR	"Salsas" [probably erroneous for Étang de Salses, situated between Leucate / Le Barcarès]; Meditarrenean Sea		PRESENT		0674aea0-a7e1-11d8-9534-b8a03c50a862	42.85	3.0							1951-10-21T00:00:00	21	10	1951	5211367	5211367	PRESERVED_SPECIMEN	SNSB-ZSM	ZSMpiscescoll>>. If the fields are not quoted (e.g. field separator does not appear within any field), try quote="" to avoid this warning.
####  Download data for amphibians

# prepare query to download
#occ_download(
#  pred("hasGeospatialIssue", FALSE),
#  pred("hasCoordinate", TRUE),
#  pred("occurrenceStatus","PRESENT"),
#  pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
#  pred_within("POLYGON((-0.55 55.01, -0.55 43.92, 15.35 43.92, 15.35 55.01, -0.55 55.01))"),
#  pred("taxonKey", 131),
#  format = "SIMPLE_CSV"
#)
## GBIF Occurrence Download https://doi.org/10.15468/dl.4hs495 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-06-07
#occ_download_wait('0014202-230530130749713')
#
#d <- occ_download_get('0014202-230530130749713') %>%
#  occ_download_import()

# Warning message:
#   In data.table::fread(targetpath, data.table = FALSE, fill = fill,  :
#                          Found and resolved improper quoting out-of-sample. First healed line 103923: <<3751555398	b00e720c-6a61-40ee-b52b-f195e310c11a	Pysanets.amph.664	Animalia	Chordata	Amphibia		Urodela				FAMILY	Urodela	Lissotriton montandoni (Boulenger, 1880)		UA		Lvivs'ka Oblast'	PRESENT		ca2fd897-6108-4361-91f8-b39dc8d12d13	49.016667	23.566667							2001-05-01T00:00:00	1	5	2001	8439630		MATERIAL_CITATION							CC_BY_4_0		"Усне повідомлення: Гринчишин Т.			2023-04-21T06:00:38.464Z		COORDINATE_ROUNDED;COUNTRY_INVALID;CONTINENT_DERIVED_FROM_COORDINATES;TAXON_MATCH_HIG>>. If the fields are not quoted (e.g. field separator does not appear within any field), try quote="" to avoid this warning.

library("ridigbio")

# Júcar  POLYGON((-9.86 46.4, 6.29 46.4, 6.29 34.74, -9.86 34.74, -9.86 46.4))
# Danube POLYGON((8 51, 30 51, 30 42, 8 42, 8 51))


### FISHES

for(i in 1:dim(Fam)[1]){

rq <- list(family=Fam[i,1], geopoint=list(
  type="geo_bounding_box",
  top_left=list(lat=51, lon=8),
  bottom_right=list(lat=42, lon=30)
))

test =  idig_search_records(rq)
test = test[,c(1,4,5,6,7,9,10,11)]
nombre = paste0("/mnt/shared/sosw/tmp/tmp_idigbio/idigbio_", i, ".txt")
write.table(test, file=nombre, quote=FALSE, row.names = FALSE, sep="|")
}

#rq <- list(family=Fam[1,1], geopoint=list(
#  type="geo_bounding_box",
#  top_left=list(lat=8, lon=51),
#  bottom_right=list(lat=30, lon=42)
#))
#
##idig_count_records(rq)   # 4833
#test = idig_search_records(rq)
#test = test[,c(1,4,5,6,7,9,10,11)]
#write.table(test, file="~/data/sosw/bio_down/fish_jucar_idigbio.txt", quote=FALSE, row.names = FALSE, sep="|")
#
####  AMPHIBIANS
#rq <- list(class="amphibia", geopoint=list(
#  type="geo_bounding_box",
#  top_left=list(lat=55.01, lon=-0.55),
#  bottom_right=list(lat=43.92, lon=15.35)
#))
#
#idig_count_records(rq)   # 3747
#test = idig_search_records(rq)
#write.table(test, file="idigbio_amphibia.txt", quote=FALSE, row.names = FALSE)

