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

# now that we have the genotypes we want to use, we will need to put them into a wide-format
# data frame.  A job for dcast...
wide <- dcast(data = gc3, full.name ~ assay.name, value.var = "single_call")


# set the NAs to 0s
wide[is.na(wide)] <- 0


#### Attaching some meta data ####

# we are going to want to put short labels on there, and some pop columns, sorted the way
# we want. So lets do that.
needed_ids <- read.csv("data/meta/sample_sheet.csv", stringsAsFactors = FALSE) %>%
  tbl_df %>%
  mutate(group_name = factor(group_name, levels = c("reSac", "rjSac", "rnSac", "rjKla", "nnEel", "nnByc"))) %>%
  arrange(group_name, pipe_name) %>%
  mutate(gn1 = as.integer(group_name),
         gn2 = as.integer(group_name),
         full.name = NMFS_DNA_ID) %>%
  select(full.name, pipe_name, group_name, gn1, gn2) 


# and left join them on there
AllOfEm <- left_join(needed_ids, wide)


# then toss the ones with too much missing data
NoHiMissers <- how_many_loci_typed %>%
  filter(n < 60) %>%
  anti_join(AllOfEm, .)

# then toss one more that had some NAs and sort them correctly and tweeze off uneeded columns
NoHiMissers2 <- NoHiMissers[complete.cases(NoHiMissers), ] %>%
  arrange(gn1, pipe_name) %>%
  select(-full.name, -group_name)



if(!file.exists("StructureArea/data")) create.dir("StructureArea/data")
if(!file.exists("StructureArea/input")) create.dir("StructureArea/input")



# now write to a file and prepare all the commands to run structure
outf <- "StructureArea/data/genos_slg_pipe.txt_dat001"
cat(names(NoHiMissers2)[-(1:3)], "\n", file = outf)
write.table(NoHiMissers2, file = outf, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

# now write a few other funny things we need to do:
TheN <- nrow(NoHiMissers2)
cat(TheN, file = "StructureArea/input/N.txt", eol = "\n")

TheL <- ncol(NoHiMissers2) - 3
cat(TheL, file = "StructureArea/input/L.txt", eol = "\n")

cat(outf, file = "StructureArea/input/InputFileNames.txt", eol = "\n")

write.table(cbind(1:length(levels(NoHiMissers$group_name)), levels(NoHiMissers$group_name)),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", 
            file = "StructureArea/data/pop_idxs.txt")

write.table(cbind(1:length(levels(NoHiMissers$group_name)), levels(NoHiMissers$group_name)),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", 
            file = "StructureArea/clump_and_distruct/pop_names.txt")

# Now the hard part, the commands
Kvals <- 2:4
Reps <- 3
KvalR <- rep(Kvals, each = Reps)
reppy <- rep(1:length(Kvals), Reps)
set.seed(5)
struct_seeds <- round(runif(n = length(KvalR), min = 1, max = 10^6))

Comms <- paste(
  "echo \"Starting Rep ",
  reppy,
  " of K = ", 
  KvalR, 
  " for data set genos_slg_pipe.txt_dat001 at $(date)\"; ../bin/structure  -K ",
  KvalR, 
  " -i ../data/genos_slg_pipe.txt_dat001  -N ",
  TheN, 
  " -L ",
  TheL,
  " -D ", 
  struct_seeds,
  " -o StructOuput_genos_slg_pipe.txt_dat001_k00",
  KvalR,
  "_Rep00",
  reppy,
  ".txt > StdoutStruct_genos_slg_pipe.txt_dat001_k00",
  KvalR,
  "_Rep00",
  reppy,
  ".txt;  echo \"Done with Rep ", 
  reppy, 
  " with K = ",
  KvalR,
  " for data set genos_slg_pipe.txt_dat001 at $(date)\"",
  sep = ""
  )

cat(Comms, sep = "\n", file = "StructureArea/input/Commands.txt")


## Now, fire that thing off, here it calls for 6 processors
system("cd StructureArea/arena; nohup ../script/ExecuteStructureRuns.sh  6  > BIG_LOG.txt  2>&1 &")

#### Now, as long as I am at this, I should try to write to gsi_sim files ####
if(FALSE) {
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

}



