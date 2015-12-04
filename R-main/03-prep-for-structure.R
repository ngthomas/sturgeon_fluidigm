####  Prepare scored sturgeon genotypes for a structure run or two ####

library(dplyr)
library(reshape2)
library(stringr)


#### Data reading and munging  ####

# first, read the data in
#g <- tbl_df(read.csv(gzfile("data/new_score_allplates.csv.gz"), stringsAsFactors = FALSE))
g<- readRDS("outputs/genotype_from_five_chips.rds")

# now, we have an issue to deal with. Some individuals were scored on multiple
# plates.  Here is what I propose we do:
# 1. Toss out all genotypes that are 0's  
gc <- g %>% filter(new.k > 0)

# 2. Now count up the number of individuals typed on more than one plate at the same locus
#    that got genotypes that are discordant
discord <- gc %>% 
  group_by(full.name, assay.name) %>%
  summarise(num = n(), nd = n_distinct(new.k)) %>%
  filter(num >1, nd > 1) %>%
  ungroup

# count how many per individaul that is:
discord %>%
  group_by(full.name) %>%
  tally()

# note that the inconsistencies are few, apart from some that are clustered
# into some clearly lousy samples. We are going to toss out the 
# discordant ones.  

# 3. Toss out the genotype calls (at 17 indiv x locus combos) that were discordant
gc2 <- gc %>% anti_join(discord)


# now, we can summarise those multiple calls of the remaining ones by just the 
# first genotypes (since they are all the same)
gc3 <- gc2 %>%
  group_by(full.name, assay.name) %>%
  summarise(single_call = first(new.k)) %>%
  ungroup


# see how much missing data different individuals have
how_many_loci_typed <- gc3 %>%
  group_by(full.name) %>%
  tally()

# and toss out indivs with fewer than 60 loci typed
gc4 <- how_many_loci_typed %>%
  filter(n < 60) %>%
  anti_join(gc3, .)


# so, we have tossed about 44 individuals because they had too much missing data.
# We could be less stringent in the future
length(unique(gc3$full.name))
length(unique(gc4$full.name))

length(unique(gc3$full.name)) - length(unique(gc4$full.name))


# now that we have the genotypes we want to use, we will need to put them into a wide-format
# data frame.  A job for dcast...
wide <- dcast(data = gc4, full.name ~ assay.name, value.var = "single_call")


# set the NAs to -9s
wide[is.na(wide)] <- -9


#### Attaching some meta data ####

# we are going to want to put a PopID on here, so we need the meta data
meta <- tbl_df(read.table("data/meta/AM001_AM006.tab", sep = "\t", header = T, stringsAsFactors = FALSE))

meta2 <- meta %>% 
  select(NMFS_DNA_ID, WATERSHED, LOCATION_COMMENTS_M, SAMPLE_COMMENTS)

meta3 <- meta2 %>% mutate(origin = paste(WATERSHED, LOCATION_COMMENTS_M, sep = ""))
meta3$origin[meta3$origin == ""] <- "Unknown Origin"
meta3$origin = factor(meta3$origin,
                      levels = c("Sacramento River",
                                 "Klamath River",
                                 "Groundfish fishery",
                                 "Unknown Origin",
                                 "Mystery"))

meta4 <- meta3 %>% select(NMFS_DNA_ID, origin)

merged <- left_join(wide, meta4, by = c("full.name" = "NMFS_DNA_ID")) %>%
  select(full.name, origin, starts_with("ame")) %>%
  tbl_df() %>%
  arrange(origin)

# now, some of the fish in wide weren't in the meta data file...call them Mystery
merged$origin[is.na(merged$origin)] <- "Mystery"

#### Then write that out into a structure file ####
to_write <- merged
to_write$origin <- as.integer(to_write$origin)

outf <- "struct_input.txt"
cat(names(to_write)[-(1:2)], "\n", file = outf)
write.table(to_write, file = outf, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")


#### Now, as long as I am at this, I should try to write to gsi_sim files ####

# Since it is pretty obvious that all the fish (adults and juveniles) from the 
# sacramento and the Klamath are from the right spot, we will include all of them 
# into the baseline (juvies and adults)

# get the locus stuff to write:
locs <- cbind("PLOIDY 1", names(merged)[-(1:2)])



glist <- split(merged, merged$origin)

# write the baseline file
outf <- "gs_baseline.txt"
cat(nrow(glist$`Sacramento River`) + nrow(glist$`Klamath River`), nrow(locs), "\n", file = outf)
write.table(locs, quote = FALSE, row.names = FALSE, col.names = FALSE, file = outf, append = TRUE)
cat("POP  Sacto\n", file = outf, append = TRUE)
write.table(glist$`Sacramento River`[, -2], quote = FALSE, row.names = FALSE, col.names = FALSE, file = outf, append = TRUE)
cat("POP  Klamath\n", file = outf, append = TRUE)
write.table(glist$`Klamath River`[, -2], quote = FALSE, row.names = FALSE, col.names = FALSE, file = outf, append = TRUE)


# write the mixture file
outf <- "gs_mixure.txt"
cat(nrow(glist$`Groundfish fishery`) + nrow(glist$`Unknown Origin`) + nrow(glist$`Mystery`), nrow(locs), "\n", file = outf)
write.table(locs, quote = FALSE, row.names = FALSE, col.names = FALSE, file = outf, append = TRUE)
cat("POP  Mixture\n", file = outf, append = TRUE)
lapply(levels(merged$origin)[-(1:2)], function(z) 
  write.table(glist[[z]][, -2], quote = FALSE, row.names = FALSE, col.names = FALSE, file = outf, append = TRUE)
)





