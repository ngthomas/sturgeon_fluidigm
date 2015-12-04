library(ggplot2)
library(plyr)
library(dplyr)
library(stringr)
source("R/chip_scoring_functions.R")



# Read in the scored raw intensity data for the four "training" plates
scored_training_plates <- read.csv("data/training_chips/training_chips_scored.csv", stringsAsFactors = FALSE, row.names = 1) %>%
  tbl_df

# read in raw intensity data for plate 5
plate5 <- read_fluidigm_detailed_csv("data/more_chips/1382136064_raw.csv") %>%
  mutate(plate = 5, long_plate_name = 1382136064)

# combine those
Dat <- bind_rows(scored_training_plates, plate5)



# call the first four plates 1234 and the fifth plate 5:
Dat2 <- Dat %>% 
  mutate( plate = ifelse(plate < 5, 1234, 5),
          #plate = as.character(1*(plate<5)+2*(plate>=5)),
         long_plate_name = as.character(long_plate_name)) %>%
  select(plate:Allele.Y.1, long_plate_name)



# make the plots 
MultiChipRelativeIntensityPlot(Dat2, 
                  prefix="outputs/plate_x_y_plate_five_and_training_plates_", 
                  n.pages=4, 
                  self.exclude = T)



#### Now, call the actual genotypes on all these, and write it to a file for future use ####
# For this we need to read in the scored version of plate 5.  This has the NTC normalized intensities,
# but that is OK.  
# read in raw intensity data for plate 5
plate5_scored <- read_fluidigm_detailed_csv("data/more_chips/1382136064_scored.csv") %>%
  mutate(plate = 5, long_plate_name = 1382136064)

# combine those
scory_dat <- bind_rows(scored_training_plates, plate5_scored)

Dat.reorg <- ReOrganizeFile(scory_dat)

# reassign K value: the genotype cluster. This allows us to use
# Fluidigm's 3-genotype scoring system to score loci with 5 genotypes.
Dat.K <- Dat.reorg %>%
  group_by(assay, plate) %>%
  mutate(new.k = SpanningK(k, rel.dye1, rel.dye2)) %>%
  ungroup() %>%
  group_by(assay) %>%
  mutate(total.k = max(new.k))


saveRDS(ungroup(Dat.K), file = "outputs/genotype_from_five_chips.rds", compress = "xz")
