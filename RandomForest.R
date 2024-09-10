library(ranger)
library(data.table)
library(hydrographr)

spp = read.csv("/mnt/shared/sosw/danube/tmp/Huso_huso/pa_env_Huso_huso.csv")
spp = read.csv("../model_table.csv")
spp = read.csv("~/data/sosw/model_table.csv")
spp$PresAbs = as.factor(spp$PresAbs)

# number of presence records:
pres_num <- as.numeric(table(spp$PresAbs)["1"])
sample_size <- c("0" = pres_num / nrow(spp),
                 "1" = pres_num / nrow(spp))

model <- ranger(spp$PresAbs ~ .,
                data = spp[,c(3,4,9,10,15,16,41,42,73,77,78,80,81,91,94,114,143)],
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

#save(model, file="/mnt/shared/sosw/danube/tmp/Huso_huso/rf.RData")


##########
######    PREDICTION

### load the prediction table

predtb = read.csv("/mnt/shared/sosw/danube/env/pred/danube_predict.csv")
predtb = read.csv("~/data/sosw/danube_predictTB.csv")
predtb = na.omit(predtb)

predtb = read.csv("predTB_1")

## prediction has to be run in tiles or subtables

# define the size of the tiles
n <- 500000

nr <- nrow(predtb)
stats_table_split <- split(predtb,
                                   rep(1:ceiling(nr/n),
                                       each=n, length.out=nr))
gc()

# create the output list
pred_list <- list()

# run prediction in a loop.
# The rbindlist combines the i-item as a datatable for prediction
# [,!1] to remove the subcatchment column that is not needed
loop_length <- length(stats_table_split)

for(i in 1:loop_length) {
  cat("Now predicting on chunk", i, "\n")
  pred_list[[i]] <- predict(model, data = rbindlist(stats_table_split[i]) [,-1])
}

pred = predict(model, data = predtb[,c(2,3,8,9,14,15,40,41,72,76,77,79,80,90,93,113,142)])


pred_list_combine <- list()
for(i in 1:loop_length) {
  cat("Now combining chunk", i, "\n")
  pred_sub <- setDT(as.data.frame(pred_list[[i]]))
  names(pred_sub) <- c("prob_0", "prob_1")
  pred_list_combine[[i]] <- pred_sub
  rm(pred_sub)
}

# combine
pred_all <- rbindlist(pred_list_combine)
names(pred_all) <- c("prob_0", "prob_1")


## Attach the subc_id
pred_all <- cbind(predtb[,1], pred_all)

#--- multiply the probability values by 100, to convert them to integers---
pred_all$prob_1 <- as.integer(round(pred_all$prob_1, 2) * 100)
#write.table(pred_all , file="/mnt/shared/sosw/danube/env/pred/prediction.txt", 
#            sep=" ", row.names=FALSE)
#pred_all = read.table("/mnt/shared/sosw/danube/env/pred/prediction.txt", header=TRUE)

# Add Row using rbind()
new_row = c(V1 = "*", prob_0 = 0, prob_1 = NULL)
pred_all2 = rbind(pred_all,new_row)

#---reclassify the sub-catchment raster---

reclass_raster(
  data = pred_all,
  rast_val = "V1",
  new_val = "prob_1",
  raster_layer = "/mnt/shared/sosw/tmp/danube_subcatchments.tif",
  recl_layer = "/mnt/shared/sosw/danube/env/pred/pred_Huso_huso_sub.tif",
  read = FALSE,
  bigtiff = TRUE)

