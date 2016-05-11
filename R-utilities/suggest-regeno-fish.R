# having scored 4 training chips, and then a 5th chip (that was not a training chip)
# we are now going to be doing more chips.  When we run these new chips we are going to
# want to run some individuals that we have previously run, so that we can identify
# clusters correctly.  At the same time, although it is probably not totally necessary,
# we should regenotype some individuals that we have only typed once before (and which
# were not on the training chips) so we can put them on without "knowing" they are
# duplicates, so that we can verify we can re-genotype them accurately (and cosistently.)

# I propose we do 12 individuals of the first kind and 12 of the second.  Obviously we
# want the individuals to be ones with high quality DNA.  So, here we do that.

library(dplyr)

# read in the data for all the plates we have done to this point:
dat <- readRDS("outputs/genotype_from_five_chips.rds")


# get the sample sheet:
sams <- read.csv("data/meta/sample_sheet.csv", stringsAsFactors = FALSE) %>% tbl_df


# it would be good to know which bycatch individuals are duplicates so we can toss the dupie ones out
byc <- readRDS("data/meta/bycatch_IDS.rds")


# so, just keep the fish that we want to keep. 
dat2 <- dat %>%
  filter(full.name %in% sams$NMFS_DNA_ID & !(full.name %in% byc$Duplicate_Tissue)) 



# now, we want to count up how many times each of these fish has been genotyped. We want the number
# of training plates and the number of times on plate 5
gcounts <- dat2 %>%
  group_by(full.name) %>%
  summarize(training_genos = sum(unique(plate) %in% 1:4),
            plate5genos = sum(unique(plate) %in% 5))


# and we also want individuals with a high average number of 
# markers that were not missing.  So, let's count that up to
not_miss <- dat2 %>%
  group_by(assay.name) %>%
  filter(sum(new.k) > 0) %>% # this tosses out loci that are not scored for anyone
  group_by(full.name) %>%
  summarise(ave_non_miss = mean(new.k != 0))

# check it out:
hist(not_miss$ave_non_miss, breaks = 20)

# so, I think we should shoot for ave_non_miss > 0.95

# so, we whittle it down to some candidates and then we get some DPS assignments on there from structure, too.
candidates <- left_join(gcounts, not_miss) %>%
  filter(ave_non_miss > 0.95) %>%
  left_join(., readRDS("outputs/structure_dps_assigns.rds")) %>%
  rename(NMFS_DNA_ID = full.name)



# and now, let us order these so that if we go down the list we get
# 12 genotyped only training chips, 12 genotyped only on plate 5 and
# we will alternate them north and south DPSs.  So, in fact, we should
# go north-north, south-south in order and so on.
set.seed(78)
final <- candidates %>%
  filter((training_genos == 1 & plate5genos == 0) |  (training_genos == 0 & plate5genos == 1)) %>%
  mutate(group = paste(training_genos == 1, DPS == "north")) %>%
  group_by(group) %>%
  mutate(sort_val = sample(1:n())) %>%
  ungroup %>%
  arrange(sort_val, group)


# print them out to send to carlos and ellen
write.csv(final, file = "outputs/regeno-suggestions.csv", row.names = FALSE)

