library(plyr)
library(dplyr)
library(ggplot2)

#setwd("/Users/eric.crandall/Desktop/tng/fluidigm/data/git/sturgeon_fluidigm/")
source("R/chip_scoring_functions.R")

chips <- c("1381905043",  "1381935339", "1381992037", "1381992302")
gN <<- 0

four_plates <- lapply(chips, function(chip) {
  gN <<- gN + 1  # increment chip counter  
  read.csv(paste0("data/four_scored_chips/recall1_",chip,".csv"), stringsAsFactors = FALSE) %>% 
    cbind(plate = gN, long_plate_name = chip, .)  
  }) %>%
  rbind.fill()

plates.reorg <- ReOrganizeFile(four_plates)

# verify consistency in regional labels

meta.file <- read.csv("data/meta/acipenser.mer",
                      stringsAsFactors=FALSE)

#for i in *_raw.txt; do awk -v f=`echo $i|cut -d "_" -f 1` 'd[$2]++==0 && /AM/{print f"\t" $2}' $i; done > raw_sample.txt

crandall.rawlabel <- read.table("data/crandall_files/raw_sample.txt")
crandall.tbl <- tbl_df(crandall.rawlabel)
colnames(crandall.tbl) <- c("region", "name")

ggplot(crandall.tbl, aes(x=name, y=region))+
  geom_point()

meta.slim <- meta.file %>%
  rename(name = NMFS_DNA_ID,
         water = WATER_NAME,
         stage = REPORTED_LIFE_STAGE,
         batch.id = BATCH_ID) %>% 
  select(name, water, stage, batch.id)

crandall.reg.call <- read.csv("/Users/eric.crandall/Desktop/tng/fluidigm/data/preprocess/regional.tbl")
colnames(crandall.reg.call) <- c("name", "define.region")
