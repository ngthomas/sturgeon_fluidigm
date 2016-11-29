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


#### Intermediate Step:  Explore missing data amongst bycatch ####
# one of the referees seemed to be concerned that bycatch individuals might
# have variation that was not part of the original chip-scoring scheme and
# thought that it was possible that some of the missing data amongst bycatch samples
# (which i think is pretty small anyway) could be the result of unscored 
# true variation.  I doubt it.  But I wanted to look here at the distribution of 
# distances from the origin for the genotypes that are not called as a fraction of that
# same distance for called genotypes, to see how many of the uncalled genotypes
# just have low illuminance.  

bycids_here <- read.csv("data/meta/sample_sheet.csv", stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  filter(category == "non-reference") %>%
  select(NMFS_DNA_ID) %>%
  unlist() %>%
  unname()
  

# get just the bycatch genos and add a z-score for the radius on it
byc_genos <- Dat.K %>%
  filter(full.name %in% bycids_here) %>%
  filter(total.k > 0)  %>% # only loci that we try to call
  mutate(Geno_called = new.k > 0,  # whether or not it is called
         radius = sqrt(rel.dye1^2 + rel.dye2))  %>%
  group_by(assay.name) %>%     # group by assay to compute a radius Z-score for each point
  mutate(zmean = mean(radius[Geno_called == TRUE]),
         zsd = sd(radius[Geno_called == TRUE]),
         z_radius = (radius - zmean) / zsd)  # this line computes the z-score of the radius


# now we can plot histograms of all these
ggplot(byc_genos, aes(x = z_radius, fill = Geno_called)) +
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity")
  
# this shows that some of the no-calls clearly are not just cases of low illuminance, which
# is good to know.



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



