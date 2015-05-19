library(plyr)
library(dplyr)
library(ggplot2)

#setwd("/Users/eric.crandall/Desktop/tng/fluidigm/data/git/sturgeon_fluidigm/")
source("R/chip_scoring_functions.R")


#### Step 1: Read in the data for the four "training" plates ####
chips <- c("1381905043",  "1381935339", "1381992037", "1381992302")
gN <<- 0

four_plates <- lapply(chips, function(chip) {
  ff <- paste(chip, "csv", sep = ".")   # filename
  path <- file.path("data/four_training_chips", ff)  # file path
  gN <<- gN + 1  # increment chip counter
  
  read.csv(path, skip = 15, stringsAsFactors = FALSE) %>% # read it in 
    cbind(plate = gN, long_plate_name = chip, .)  #%>% # add the plate number and name columns
    #dplyr::rename(., rel.dye1 = Allele.X.1, rel.dye2 = Allele.Y.1)  # correctly name the relative dye intensity columns
}) %>%
  do.call(what = rbind, args = .)  # in the end, rbind them all together.

# at this point "four_plates" is the long format data that we need


#### Step 2:  Plot the four training plates together at each locus, connecting re-genotyped individuals ####
if(!file.exists("outputs")) dir.create("outputs")
MultiChipRelativeIntensityPlot(four_plates, prefix = "outputs/plate_x_y")
