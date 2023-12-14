library(sp)
library(raster)

setwd("/mnt/shared/sosw")
fish = read.csv("tmp/fish_danube_fil4.csv")
coordinates(fish) = ~long+lat
fish$Occu = 1

r <- raster( xmn=8, xmx=30, ymn=42, ymx=51)
res(r) = 0.25
x <- rasterize(fish, r, fun='count')
xl = log(x)
writeRaster(x[["Occu"]], "tmp/Counts2D.tif", overwrite=TRUE)
writeRaster(xl[["Occu"]], "tmp/Counts2Dl.tif", overwrite=TRUE)


####   How many species ?
length(levels(fish$species))   ## 164 species

# plant
Woodsia ilvensis

# Marine
Trachinus draco
Scorpaena porcus
Salvelinus malma

