library(plyr)
library(dplyr)
library(ggplot2)

source("R/chip_scoring_functions.R")


#### Step 1: Read in the data for the four "training" plates ####
chips <- c("1381905043",  "1381935339", "1381992037", "1381992302")

four_plates <- lapply(chips, function(chip) {
  ff <- paste(chip, "_raw.csv", sep = "")   # filename
  path <- file.path("data/training_chips", ff)  # file path
  tmp <- read.csv(path, skip = 15, stringsAsFactors = FALSE) 
  bottom_line <- which(tmp$ID == "Dose Meter Reading Data") - 1  # remove the bottom, irrelevant part of each file.
  tmp <- tmp[1:bottom_line,]
  tmp$long_plate_name = chip
  tmp
  }) %>%
  bind_rows(.id = "plate") %>%     # in the end, bind them all together.
  tbl_df

# at this point "four_plates" is the long format data that we need


#### Step 2:  Plot the four training plates together at each locus, connecting re-genotyped individuals ####
if(!file.exists("outputs")) dir.create("outputs")
MultiChipRelativeIntensityPlot(four_plates, prefix = "outputs/plate_x_y")




#### Step 3:  Read in how those four chips were scored using the fluidigm software ####
# Fluidigm is expecting only 3 genotypes, so, if we are scoring five genotypes, we call
# them one of three genotypes, but if we call genotype XX way up in the upper corner,
# we post-process that later and realize it is a 4th cluster.  That is taken care of
# by the spanningK function.

scored_plates <- read.csv("data/training_chips/training_chips_scored.csv", stringsAsFactors = FALSE, row.names = 1) %>%
  tbl_df


data.reorg <- ReOrganizeFile(scored_plates)

# reassign K value: the number of cluster due to limitations of fluidigm' scoring capability
data.K <- data.reorg %>%
  group_by(assay, plate) %>%
  mutate(new.k = SpanningK(k, rel.dye1, rel.dye2)) %>%
  ungroup() %>%
  group_by(assay) %>%
  mutate(total.k = max(new.k))


# and now we can get a data frame of the number of clusters that we
# score at each locus:
num_geno_clusts_df <- data.K %>%
  group_by(assay.name) %>%
  summarize(nClusters = max(total.k)) 



MultiChipRelativeIntensityPlot(data.K, 
                               color.by = "genotype", 
                               dont_reorganize = TRUE, 
                               prefix = "outputs/plate_x_y_by_final_genotype")
