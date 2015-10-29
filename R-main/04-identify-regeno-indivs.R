# Our goal here is to identify some individuals that would be good to 
# regenotype.  Our goals might be manifold here.  I think it would be
# good to have individuals that gave good genotypes, but have only been
# genotyped once before (and so were not already used as re-scores
# to develop the scoring method).  But, we might want to regenotype
# individuals that didn't score well, too.  In that case, it will be best
# to focus on individuals from populations from which we have
# few samples.  And we should keep in mind that Plate 3 was just sort
# of crappy, with really low success.



# get all the multiple scores:
g <- tbl_df(read.csv(gzfile("data/new_score_allplates.csv.gz"), stringsAsFactors = FALSE))

# figure out which loci are ones that we don't score at all because they are always 0's
tossers <- g %>% 
  group_by(assay, new.k) %>% 
  tally() %>%
  group_by(assay) %>%
  mutate(num_genos = n()) %>%
  filter(num_genos == 1) %>%
  select(assay)


# and then remove those loci from our consideration
gclean <- anti_join(g, tossers)


# now, we summarize this to show, for each individual (on each plate) the proportion
# of loci that were scored:
prop_scored <- gclean %>%
  group_by(name, plate) %>%
  summarise(prop_success = sum(new.k != 0) / n()) %>%
  group_by(name) %>%
  mutate(num_plates = n()) %>%
  ungroup
  
# and now we would also like to associate some meta data with those
meta <- read.csv("data/meta/acipenser.mer", stringsAsFactors = FALSE)

ppsm <- left_join(prop_scored, meta, by = c("name" = "NMFS_DNA_ID"))

# now, just out of curiosity, look at gtyping success as a function of sample comments:
success_by_plate_and_comment <- ppsm %>%
  group_by(plate, SAMPLE_COMMENTS) %>%
  summarize(ave_success = mean(prop_success),
            num_samps = n(),
            min_succ = min(prop_success),
            max_succ = max(prop_success)
            ) %>%
  ungroup


# Look at this:
nasties <- success_by_plate_and_comment %>%
  arrange(ave_success) %>% as.data.frame

# it shows us that we should just plain avoid ones that have sample comments like: "GST Egg", "died in trap",
# etc.  So, I am going to just avoid the ones that have low average values.
avoid_comms <- nasties %>%
  filter(ave_success < .81) %>%
  select(SAMPLE_COMMENTS)


# now, we also are going to want to try to pump up our baselines from the northern DPS,
# so let us see how many of these things we have:
ppsm %>% 
  group_by(WATER_NAME, REACH_SITE) %>% tally()

# so, we are really low on Klamath,   let's see how many of those Klamath fish yielded good genotypes
ppsm %>%
  filter(WATER_NAME == "Klamath River") %>%
  select(name:num_plates) %>% as.data.frame

# that merely tells us that our Klamath fish have had great success and many
# have been genotyped multiple times.  I don't think there is much that we can do
# there, except possibly regenotype the ones that have been genotyped
# only once (to make sure that Klamath fish are represented in our "newly regenotyped samples")


# now, I think we should try to associate a cluster (i.e. DPS from Structure) for all the fish here
# so that we can be sure to represent them both in the freshly-regenotyped fish.
gsi <- read.table("outputs/hasty_gsi_sim_results.txt", header = TRUE, stringsAsFactors = FALSE) %>%
  tbl_df

ppsm_w_gsi <- left_join(ppsm, gsi) %>%
  select(name:num_plates, assigned_to, ass_score, WATER_NAME, everything())


# now, let's see how many individuals have been genotyped only once and
# sort them out by proportion of success, and then remove individuals that
# are from the "dodgy SAMPLE_COMMENTS" categories and see what it looks like:
Candidates <- ppsm_w_gsi %>% 
  filter(num_plates == 1) %>%
  arrange(prop_success) %>%
  filter(!(SAMPLE_COMMENTS %in% avoid_comms$SAMPLE_COMMENTS ))

# holy smokes! When we do that we see that there are only two indivdiduals with
# lower than 85% genotyping success. So, let's include those individuals and then 
# just randomly order them with so that they alternate by Klamath and Sacramento
# so that we have even representation of those two groups in our
# "newly regenotyped" fish.  We make "origin" to represent the origin of the fish...
set.seed(345)
Final <- Candidates %>%
  mutate(origin = ifelse(is.na(assigned_to), 
                         ifelse(WATER_NAME == "Klamath River", "Klamath", 
                                ifelse(WATER_NAME == "Sacramento River", "Sacto", NA)),
                         assigned_to)
  ) %>%
  group_by(origin) %>%
  mutate(sort_val = sample(1:n())) %>%
  ungroup %>%
  arrange(sort_val, origin) %>%
  select(name:WATER_NAME, origin, sort_val, everything()) %>%
  rename(NMFS_DNA_ID = name)
  


write.csv(Final, "outputs/suggested_fish_to_regenotype.csv", row.names = FALSE)
