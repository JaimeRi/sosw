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
  pred_within("POLYGON((101.921002148438 23.4160971210938,102.505014916992 12.2723738056641,102.732922338867 12.1536720234375,102.87061640625 11.8569175678711,103.167370861816 11.7073533222657,103.150752612305 11.3749883320313,103.173305950928 11.1910005695801,103.176867004395 10.9702152546387,103.311000018311 10.953597005127,103.410709515381 11.1043482685547,103.535346386719 11.1921875874024,103.652861151123 11.082981947754,103.752570648193 10.8681317219239,103.720521166992 10.7636741535645,103.705757632828 10.7323665585023,103.57043760109 10.6397791683655,103.668366571426 10.5471917782288,103.887371359634 10.6308765346986,104.127742468643 10.6255349544984,104.307575668716 10.5934854732972,104.412626745987 10.4759707088929,104.569313098526 10.3851638454896,104.626289953995 10.3032596157533,104.752707352066 10.2783322414857,104.852416849136 10.2409411800843,104.936101605606 10.1625980038148,105.091007431412 10.0806937740785,105.153325867081 9.98098427700816,105.14798428688 9.84210319180309,105.073202164078 9.81539529080211,105.04293320961 9.84388371853649,104.907613177872 9.71212474026501,104.905832651138 9.62665945706189,104.852416849136 8.92157087063615,104.969931613541 8.78625083889787,104.966370560074 8.6473697536928,105.308231692886 8.91088771023576,105.717752841568 9.28479832424941,106.095224509048 9.45216783718885,105.689264413834 10.1074016750794,106.800313095474 10.7519523525696,107.281055313492 10.605949160431,107.398570077896 10.5204838772279,108.384981888199 11.1579124477845,108.851479892349 11.410747243927,109.097192581558 12.0837863491516,109.15060838356 13.0773202663879,109.043776779556 14.3592995144347,108.691232486343 15.2353186672667,108.103658664322 15.9617735744932,107.334471115494 16.5920800381162,107.152857388687 16.8765191837766,107.014644001007 17.0988624596097,106.589988375091 17.4514067528226,106.413716228485 17.7298366207578,106.435750246811 17.9381582485654,106.453778079986 17.956186081741,106.391682210159 18.0403159698941,106.188702162552 18.1818678451993,105.87087814064 18.4462760651089,105.683922833633 18.7186966553189,105.550383328628 18.9857756653286,105.670568883133 19.1994388733364,105.764046536636 19.5806941601253,105.921623152542 20.0420731499172,106.14062794075 20.1221968529201,106.378328259659 20.2637487282253,106.525221715164 20.4373500847316,106.597333047867 20.8533256428218,106.594662257767 21.0830135914302,106.832362576675 21.0963675419307,107.23565188179 21.0956998444057,107.323787955093 21.3173754227138,107.459998250198 21.3974991257167,107.679003038406 21.4555888103938,107.881983086014 21.5731035747981,108.05291365242 21.679935178802,108.261235280227 21.7734128323054,108.426824266434 22.0591873730158,107.353166646194 23.0006408833002,106.044479497147 23.4386504597162,101.921002148438 23.4160971210938))"),
  #pred_within("POLYGON((8 51, 8 42, 30 42, 30 51, 8 51))"),
  pred_in("taxonKey", KEY),
  format = "SIMPLE_CSV",
  user="jaime_ric", pwd ="2ruZNXfWaP6L.58", email="marquez@igb-berlin.de" 
)

# DANUBE = GBIF Occurrence Download https://doi.org/10.15468/dl.qc7e3s Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-10-26
# RHINE = GBIF Occurrence Download https://doi.org/10.15468/dl.r695p8 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-06-06
# JUCAR = GBIF Occurrence Download https://doi.org/10.15468/dl.q53db3 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-07-26
# MEKONG = GBIF Occurrence Download https://doi.org/10.15468/dl.b785v8 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-03-18
# Check status
occ_download_wait('0028469-231002084531237') # danube
occ_download_wait('0023072-240314170635999') # mekong

# Once the request is finished, run to download
d <- occ_download_get('0023072-240314170635999') %>%
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

