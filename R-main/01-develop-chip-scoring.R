library(plyr)
library(dplyr)
library(ggplot2)

source("R/chip_scoring_functions.R")


#### Step 1: Read in the data for the four "training" plates ####
chips <- c("1381905043",  "1381935339", "1381992037", "1381992302")

four_plates <- lapply(chips, function(chip) {
  ff <- paste(chip, "_raw.csv", sep = "")   # filename
  path <- file.path("data/training_chips", ff)  # file path
  read_fluidigm_detailed_csv(path) %>%
    mutate(long_plate_name = chip)
  }) %>%
  bind_rows(.id = "plate")      # in the end, bind them all together.
  

# at this point "four_plates" is the long format data that we need


#### Step 2:  Plot the four training plates together at each locus, connecting re-genotyped individuals ####
# and we also color according to plate
if(!file.exists("outputs")) dir.create("outputs")
MultiChipRelativeIntensityPlot(four_plates, prefix = "outputs/plate_x_y")




#### Step 3:  Read in how those four chips were scored by eye upon the fluidigm software and translate that into up to 5 clusters ####
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


#### Step 4. I make four pages of plots like the plate_x_y plots, but color things by genotype  ####
MultiChipRelativeIntensityPlot(data.K, 
                               prefix = "outputs/plate_x_y_by_final_genotype_plates_combined",
                               alreadyOrganized = TRUE,
                               color.by = "new.k",
                               exclude.seg=TRUE)




# and then make them for each plate separately
for(i in unique(data.K2$plate.name)) {
  tmp <- data.K2 %>% filter(plate.name == i)
  MultiChipRelativeIntensityPlot(tmp, 
                                 color.by = "new.k", 
                                 alreadyOrganized = TRUE, 
                                 exclude.seg=TRUE,
                                 prefix = paste("outputs/plate_x_y_by_final_genotype_plate", i, "only", sep = "_")
  )
}



#### Report Number of Clusters And store in Variable for a Table ####
# and now we can get a data frame of the number of clusters that we
# score at each locus:
num_geno_clusts_df <- data.K %>%
  group_by(assay.name) %>%
  summarize(nClusters = max(total.k)) 

# and we write that to the output
num_geno_clusts_df %>%
  write.csv(file = "outputs/number-of-genotype-clusters.csv", row.names = FALSE)

# and we store a variable for using in the paper
REPORT_num_geno_clusts_table <- num_geno_clusts_df %>%
  group_by(nClusters) %>% tally()
