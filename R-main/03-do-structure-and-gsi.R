

####  Prepare scored sturgeon genotypes for a structure run or two ####
source("R/load-packages.R")

#### Data reading and munging  ####

# first, read the data in
#g <- tbl_df(read.csv(gzfile("data/new_score_allplates.csv.gz"), stringsAsFactors = FALSE))
g<- readRDS("outputs/genotype_from_five_chips.rds")

# right here we do a funny thing---tally up how many loci there are of different
# numbers of genotype categories (we use this in the results...)
tmp <- g %>%
  filter(total.k > 0) %>%
  group_by(assay, total.k) %>%
  tally

num_loc_by_num_cat <- table(tmp$total.k)

# and let's also count stuff up as far as number of typed loci
tmp <- g %>%
  filter(total.k > 0, new.k > 0) %>%
  group_by(full.name, plate) %>%
  summarise(numL = n())

num_60_locs_or_more <- tmp %>% filter(numL >= 60) %>% dim %>% "["(1)

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


# then toss the ones with too much missing data.  Currently that means >= 40 markers
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
N_for_structure <- TheN

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
if(ReDoStructureRuns == TRUE) {
  system("cd StructureArea/arena;  ../script/ExecuteStructureRuns.sh  6  > BIG_LOG.txt  2>&1")
  system("cd StructureArea/clump_and_distruct; ./script/ClumpAndDistructAll.sh 6 ")
}

#### Now, as long as I am at this, I should try to write to gsi_sim files and run the analyses ####
# Since it is pretty obvious that all the fish (adults and juveniles) from the 
# sacramento and the Klamath are from the right spot, we will include all of them 
# into the baseline (juvies and adults).


# take all the individuals:
ForGSI <- AllOfEm[complete.cases(AllOfEm), ]  # this just removed one individual that is all NAs---apparently was not typed at all
  
  
NorthBase <- ForGSI %>%
  filter(group_name == "rjKla") %>%
  select(-pipe_name, -group_name, -gn1, -gn2) %>%
  mutate(full.name = paste("Klamath", 1:n(), sep = "_")) %>%  # name the baseline fish something informative for analyzing the self-assignment later
  select(full.name, everything())
  
SouthBase <- ForGSI %>%
  filter(group_name %in%  c("reSac", "rjSac", "rnSac")) %>%
  select(-pipe_name, -group_name, -gn1, -gn2) %>%
  mutate(full.name = paste("Sacramento", 1:n(), sep = "_")) %>%  # name em for self-assignment analysis
  select(full.name, everything())

GSI_mix <- ForGSI %>%
  filter(!group_name %in%  c("rjKla", "reSac", "rjSac", "rnSac")) %>%
  select(-pipe_name, -group_name, -gn1, -gn2)

# get the locus preamble stuff to write:
locs <- cbind("PLOIDY 1", names(NorthBase[-1]))


# write the baseline file
outf <- "outputs/gs_baseline.txt"
cat(nrow(SouthBase) + nrow(NorthBase), nrow(locs), "\n", file = outf)
write.table(locs, quote = FALSE, row.names = FALSE, col.names = FALSE, file = outf, append = TRUE)
cat("POP  SouthernDPS\n", file = outf, append = TRUE)
write.table(SouthBase, quote = FALSE, row.names = FALSE, col.names = FALSE, file = outf, append = TRUE)
cat("POP  NorthernDPS\n", file = outf, append = TRUE)
write.table(NorthBase, quote = FALSE, row.names = FALSE, col.names = FALSE, file = outf, append = TRUE)


# write the mixture file
outf <- "outputs/gs_mixure.txt"
cat(nrow(GSI_mix), nrow(locs), "\n", file = outf)
write.table(locs, quote = FALSE, row.names = FALSE, col.names = FALSE, file = outf, append = TRUE)
cat("POP  Mixture\n", file = outf, append = TRUE)
write.table(GSI_mix, quote = FALSE, row.names = FALSE, col.names = FALSE, file = outf, append = TRUE)


# now, fire off some gsi_sim analyses
# here, analyze the mixture sample
system("cd outputs; ../bin/gsi_sim-Darwin -b gs_baseline.txt -t gs_mixure.txt > gsi_dumpola.txt")

# here, do self-assignment, and process the output
system("cd outputs; ../bin/gsi_sim-Darwin -b gs_baseline.txt --self-assign | awk -F\";\" 'BEGIN {print \"FromPop  SouthScore NorthScore NumL\"} /^UNSORTED_SELF_ASS_LIKE_GC_CSV/ {print $1, $3, $6, $(10)}' | sed 's/UNSORTED_SELF_ASS_LIKE_GC_CSV:\\///g; s/_[0-9]*//g;' > self-assigment-scores.txt")


# read in the self-assignment
selfies <- read.table("outputs/self-assigment-scores.txt", header = TRUE) %>%
  tbl_df

selfies %>% 
  filter(FromPop == "Sacramento", NorthScore > 50)

selfies %>%
  filter(FromPop == "Klamath")

# make a figure of self assignment scores vs Num Loci in Sacto:
selfies_for_plot <- selfies %>%
  filter(FromPop == "Sacramento") %>%
  mutate(Assignment = ifelse(SouthScore > 50, "Correct", "Incorrect"))

selfy_plot <- ggplot(selfies_for_plot, aes(x = NumL, y = SouthScore, colour = Assignment)) +
  geom_point(alpha = 0.8) +
  xlab("Number of loci scored") +
  ylab("Posterior probability of Southern DPS") +
  theme_bw()


ggsave(selfy_plot, filename = "outputs/self-ass-plot.pdf", width = 6, height = 4)

system("pdfcrop outputs/self-ass-plot.pdf")

## Now let us also read in the mixture results and save them as an rds.
mixture_ass <- read.table("outputs/pop_pofz_full_em_mle.txt", header = T, stringsAsFactors = FALSE) %>%
  tbl_df

# then make a column for which DPS they are assigned to and write it
mixture_ass %>%
  mutate(DPS_gsi_sim = ifelse(SouthernDPS > 0.5, "south", "north")) %>%
  saveRDS(file = "outputs/gsi-sim-DPS-assignments-of-bycatch.rds")
  
