library(ranger)
library(data.table)




spp = read.csv("/mnt/shared/danube/model_table.csv")
spp$PresAbs = as.factor(spp$PresAbs)

# number of presence records:
pres_num <- as.numeric(table(spp$PresAbs)["1"])
sample_size <- c("0" = pres_num / nrow(spp),
                 "1" = pres_num / nrow(spp))

nm = names(spp[,c(3,4,9,10,15,16,41,42,73,77,78,80,81,91,94,114,143)])



model <- ranger(spp$PresAbs ~ .,
                data = spp[,nm[c(3,4,5,6,7,8,11,13,14,15,16,17)]],
                #data = spp[,c(3,4,9,10,15,16,41,42,73,77,78,80,81,91,94,114,143)],
                #data = spp[,5:14],
                num.trees = 1000,
                mtry= 6,
                replace = T,
                sample.fraction = sample_size,
                oob.error = T,
                keep.inbag = T,
                num.threads = 4,
                importance ='impurity',
                probability = T)



# create the output list
pred_list <- list()
subc_list <- list()

for(i in 1:length(nm)) {
    cat("Now predicting on chunk", i, "\n")
    arc = paste0("/mnt/shared/danube/out/predTB_", i)
    predtb = read.csv(arc)
    predtb = na.omit(predtb)

    #pred_list[[i]] <- predict(model, data = rbindlist(stats_table_split[i]) [,-1])
    pred_list[[i]] <- predict(model, data = predtb[,nm])
    subc_list[[i]] = predtb[,1]
}

pred_list_combine <- list()
subc_list_combine = list()

for(i in 1:length(nm)) {
  cat("Now combining chunk", i, "\n")
  pred_sub <- setDT(as.data.frame(pred_list[[i]]))
  subc_sub <- setDT(as.data.frame(subc_list[[i]]))
  names(pred_sub) <- c("prob_0", "prob_1")
  names(subc_sub) <- "subcID"
  pred_list_combine[[i]] <- pred_sub
  subc_list_combine[[i]] <- subc_sub
  rm(pred_sub)
}

# combine
pred_all <- rbindlist(pred_list_combine)
names(pred_all) <- c("prob_0", "prob_1")

subcALL = rbindlist(subc_list_combine)

## Attach the subc_id
#pred_all <- cbind(predtb[,1], pred_all)

#--- multiply the probability values by 100, to convert them to integers---
pred_all$prob_1 <- as.integer(round(pred_all$prob_1, 2) * 100)

reclassTB = subcALL
reclassTB[,"prob_1"] = pred_all[,2]


write.table(reclassTB, file="/mnt/shared/danube/prediction_rule.txt", 
            sep=" ", row.names=FALSE, col.names=FALSE)

awk '{print $1, "=", $2}' prediction_rule.txt > rules.txt

export RASTER=/mnt/shared/sosw/tmp/danube_subcatchments.tif
export RULES=/mnt/shared/danube/rules.txt
export OUTPUT=/mnt/shared/danube/prediction_Zingel_streber_rfy.tif
export NODATA=-9999
export TYPE=Int32
export COMPRESSION=DEFLATE
export LEVEL=2
export BTIFF=NO


# Start GRASS GIS session
sudo grass -f --gtext --tmp-location  $RASTER   <<'EOF'

    # Load raster input file
    r.in.gdal --o input=$RASTER  output=raster    --overwrite

    # Reclassify the raster according to the rules
    r.reclass input=raster output=recl_raster rules=$RULES --overwrite

    # Export reclassified raster map
    r.out.gdal input=recl_raster output=$OUTPUT type=$TYPE  format=GTiff nodata=$NODATA  --o -f -m -c createopt="COMPRESS=$COMPRESSION,ZLEVEL=$LEVEL,BIGTIFF=$BTIFF"

EOF

exit
