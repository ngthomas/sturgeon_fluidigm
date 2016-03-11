source("R/load-packages.R")
source("R/chip_scoring_functions.R")
source("R/analysis_funcs.R")



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
  mutate(total.k = max(new.k)) %>%
  ungroup


saveRDS(Dat.K, file = "outputs/genotype_from_five_chips.rds", compress = "xz")


#### Count up and plot the number of loci typed per sample ####
# we want a histogram of the number of samples typed at 0,...,74 SNPs
# faceted by plate.  And we also want to color the bars by short sample prefix

# so, first, get the sample prefixes
sams <- read.csv("data/meta/sample_sheet.csv", stringsAsFactors = FALSE) %>% 
  tbl_df %>%
  rename(full.name = NMFS_DNA_ID) %>%
  select(full.name, group_name) 

sams$group_name <- factor(sams$group_name, levels = c("reSac", "rjSac", "rnSac",  "rjKla", "nnEel",  "nnByc"))

tmp <- Dat.K %>%
  filter(total.k > 0) %>%
  left_join(sams) %>%
  rename(Chip = plate) %>%
  group_by(Chip, full.name, group_name) %>%
  summarise(`Number of assays with called genotype categories` = sum(new.k != 0))

# and now we have to get the group of the duplicate samples.  Hassle...
byctmp <- readRDS("data/meta/bycatch_IDS.rds") %>%
  filter(!is.na(Duplicate_Tissue)) %>%
  rename(full.name = NMFS_DNA_ID) %>%
  left_join(sams)
bychash <- byctmp$group_name %>% 
  setNames(byctmp$Duplicate_Tissue)

trytochange <- which(is.na(tmp$group_name))
tmp$group_name[trytochange] <- bychash[tmp$full.name[trytochange]]

# now there are only three sample with unknown group.  Let's just toss em...
tmp <- tmp %>%
  filter(!is.na(group_name))

tmp$Chip <- paste("Chip", tmp$Chip)
  
g <- ggplot(tmp, aes(x = `Number of assays with called genotype categories`, fill = group_name, order = -as.numeric(group_name))) + 
  geom_histogram(breaks = c(0:75) - 0.5) +
  facet_wrap(~ Chip, ncol = 1) +
  ylab("Number of samples") +
  scale_fill_discrete(name="Group\nShort\nName") +
  theme_bw()

ggsave(g, filename = "outputs/successful-assay-histogram.pdf", width = 6, height = 6)

system("pdfcrop outputs/successful-assay-histogram.pdf")


#### Count up  discordant pairs:  I did this but don't report on it in the paper. ##### 
boing <- count_discordant_genotype_cats(Dat.K)

# now, I could see how things change when we toss individuals that are missing more than 10 genotype calls, but I haven't done that yet.



