
# The goal here is to compute power for individual identification
# given genotype frequency categories by simulation.  So, we need
# to get those frequencies, then do the simulations, and it would be
# nice to do simulations for all the loci and for a decreasing number
# of successfully genotyped ones, and also maybe do a few different
# genotyping eror rate values

source("R/load-packages.R")
source("R/self-id.R")

#### Get allele frequencies ####
# For this we will just count up alleles from fish assigned to northern or southern
# DPS by the structure run

so <- Rrunstruct:::slurp_results("StructureArea/arena/StructOuput_genos_slg_pipe.txt_dat001_k002_Rep001.txt_f", K = 2, Rep = 1)

# now, we need to get the full names on there with a few joins...
sam <- read.csv("data/meta/sample_sheet.csv", stringsAsFactors = FALSE) %>% tbl_df

# here is the join, and some renaming and adding whether it is north or south
dps_df <- so %>%
  left_join(., sam %>% select(pipe_name, NMFS_DNA_ID), by = c("Label" = "pipe_name")) %>%
  rename(full.name = NMFS_DNA_ID) %>%
  mutate(DPS = ifelse(`2` > 0.5, "south", "north")) %>%
  select(full.name, DPS)


# now, we stick DPS designations onto the genotype categories and count the up
geno <- readRDS("outputs/genotype_from_five_chips.rds")

# first filter to just the individuals 
tmp_counts <- geno %>%
  left_join(dps_df) %>%
  filter(!is.na(DPS), total.k > 0)


# these are the counts of the observed genotypes
gc_counts <- tmp_counts %>%
  filter(new.k > 0) %>%  # chuck out ones that we don't have a DPS for (too much missing data, basically)
  group_by(DPS, assay.name, new.k) %>%
  summarise(gc_count = n(), tot.k = first(total.k))


gcc_north <- gc_counts %>% filter(DPS == "north") %>% droplevels
gcc_south <- gc_counts %>% filter(DPS == "south") %>% droplevels


# now, make a list with the gc freqs in order (1, 2, 3,...)
# and having all the necessary entries.  If none were observed in one of the 
# DPSs, we still add a bit to avoid zeroes
north_gcf <- split(gcc_north, gcc_north$assay.name) %>%
  lapply(., function(x) {
    y<-rep(1/x$tot.k[1], x$tot.k[1]);
    y[x$new.k] <- y[x$new.k] + x$gc_count
    y/sum(y)
    })

south_gcf <- split(gcc_south, gcc_south$assay.name) %>%
  lapply(., function(x) {
    y<-rep(1/x$tot.k[1], x$tot.k[1]);
    y[x$new.k] <- y[x$new.k] + x$gc_count
    y/sum(y)
  })




# here we are going to get the frequencies of unobserved genotypes
# at least in all of these fish...
tmp_miss_rate <- tmp_counts %>%
  group_by(DPS, assay.name) %>%
  summarise(miss_rate = (0.5 + sum(new.k == 0)) / (n() + 0.5))


mr_north <- tmp_miss_rate %>% filter(DPS == "north")
north_miss <- mr_north$miss_rate %>% setNames(mr_north$assay.name)
north_miss <- north_miss[names(north_gcf)]  # just to make sure they are in the same, consistent order


mr_south <- tmp_miss_rate %>% filter(DPS == "south")
south_miss <- mr_south$miss_rate %>% setNames(mr_south$assay.name)
south_miss <- south_miss[names(south_gcf)]  # just to make sure they are in the same, consistent order




full_north <- full_logl_sims(GP = north_gcf, miss_rate = north_miss, epsilon = 0.05, reps = 10^5, max_num_miss_loc = 34)
full_south <- full_logl_sims(GP = south_gcf, miss_rate = south_miss, epsilon = 0.05, reps = 10^5, max_num_miss_loc = 34)




# then look at this:
lapply(full_north, function(x) {mm <- max(x$U); nn <- mean(x$S < mm); list(mm, nn)})
lapply(full_south, function(x) {mm <- max(x$U); nn <- mean(x$S < mm); list(mm, nn)})


sims_df <- bind_rows(
  lapply(names(full_north), function(x) {
    y <- as.data.frame(full_north[[x]])
    y$L <- as.numeric(x)
    y$DPS <- "North"
    y}) %>%
    bind_rows %>%
    tbl_df,
  lapply(names(full_south), function(x) {
    y <- as.data.frame(full_south[[x]])
    y$L <- as.numeric(x)
    y$DPS <- "South"
    y}) %>%
    bind_rows %>%
    tbl_df
) %>%
  tidyr::gather(key = "Relat", value = logl, S, U)


g <- ggplot(sims_df %>% filter(L %in% c(74, 60, 50, 40)), aes(x = logl, fill = factor(L))) + 
  geom_density(alpha = 0.2) +
  facet_wrap(~ DPS, nrow = 2) +
  scale_fill_discrete(name = "Number\nof loci\ntyped in\nboth\nmembers\nof the\npair") +
  xlab("Log likelihood ratio statistic") +
  theme_bw() +
  xlim(-100, 70)


ggsave(g, filename = "outputs/lambda_densities.pdf", width = 10, height = 3)

system("pdfcrop outputs/lambda_densities.pdf")

# So, it seems clear to me that there should be a number for the logl threshold which is like
# 10.0, which will toss out all non-related indivs and have a super low false negative rate
# even for as few as 40 loci.  So, just set the cutoff at that and go for it, it seems to me.

# So, now all I need to do are those pairwise comparisons...


#### Doing the pairwise comparisons  ####

# now, we need to augment out dps_df data frame to include the North or South origin
# of the bycatch samples that were duplicately sampled because we are going to define
# pairs as follows:
# 1. Same NMFS_DNA_ID --- i.e. multiply-genotyped DNA (MGD)
# 2. Same individual but different tissue   (SIDT)
# 3. Not known to be the same  (NKS)
# once again this is going to involve dealing with the dupie-tissues
dupies <- readRDS("data/meta/bycatch_IDS.rds") %>%
  tbl_df %>%
  filter(!is.na(Duplicate_Tissue))

dps_df_augmented <- dupies %>%
  left_join(., dps_df, by = c("NMFS_DNA_ID" = "full.name")) %>%
  rename(full.name = Duplicate_Tissue) %>%
  select(-NMFS_DNA_ID) %>%
  bind_rows(., dps_df)
  
  

# here are the genotypes that include all the different times they were genotyped
mgeno <- readRDS("outputs/genotype_from_five_chips.rds") %>%
  filter(total.k > 0) %>%
  left_join(dps_df_augmented) %>%
  filter(!is.na(DPS)) %>%
  droplevels


# look -- no one multiply-genotyped on the same plate
mgeno %>%
  group_by(plate, full.name, assay.name) %>%
  tally %>%
  filter(n>1)


# but how about multiply genotyped?  Plenty
mgeno %>%
  group_by(full.name, assay.name) %>%
  tally %>%
  filter(n>1)


# so, name the indivs with full.name-plate, and prepare to pull them into a list of pairs of
# genotypes for indexing out of the matrices of lr's
bungo <- mgeno %>%
  mutate(bung_name = paste(full.name, plate, sep = "-")) %>%
  mutate(new.k = ifelse(new.k == 0, NA, new.k)) %>%
  arrange(assay.name, bung_name)


# now split that into a bunch of little data frames on locus
ntmp <- bungo %>%
  filter(DPS == "north") %>%
  select(assay.name, bung_name, new.k) 

stmp <- bungo %>%
  filter(DPS == "south") %>%
  select(assay.name, bung_name, new.k) 

n_split <- split(ntmp, ntmp$assay.name)
s_split <- split(stmp, stmp$assay.name)

# now, we want to get the loglratio matrices for these
# so we can do a massive, hideous subscripting of those values.
lr_north <- lapply(north_gcf, function(x) gcat_pair_probs(x, e = 0.05)$lr)
lr_south <- lapply(south_gcf, function(x) gcat_pair_probs(x, e = 0.05)$lr)



# now, this is an interesting operation...we are going to grab the logl-ratio
# values out by making a matrix of all pairwise comparisons of individuals for
# each locus, then storing those in a big matrix and summing things over loci.
# this would probably start to bog down with more than 10^5 pairs, but should
# work fine for the size of problem we have here.
mighty_pairs_subscript <- function(i, lr, glist) {
  lr[[i]][
    cbind(
      rep(glist[[i]]$new.k, each = length(glist[[i]]$new.k)),
      rep(glist[[i]]$new.k, times = length(glist[[i]]$new.k))
    )]
}


# get the northern DPS pairs
nsplat <- sapply(names(lr_north), function(i) mighty_pairs_subscript(i, lr_north, n_split))
np_df <- data.frame(
  id1 = rep(n_split[[1]]$bung_name, each = length(n_split[[1]]$bung_name)),
  id2 = rep(n_split[[1]]$bung_name, times = length(n_split[[1]]$bung_name)),
  LLR = rowSums(nsplat, na.rm = TRUE),
  NumL = rowSums(!is.na(nsplat)),
  stringsAsFactors = FALSE
) %>%
  tbl_df %>%
  filter(id1 < id2)  %>% # remove instances of the indiv genotype to itself and also to ensure no double-counting of pairs, report them only in sorted their order version
  tidyr::separate(id1, into = c("name1", "chip1"), sep = "-") %>%
  tidyr::separate(id2, into = c("name2", "chip2"), sep = "-") %>%
  select(name1:chip2, LLR, NumL)
  

# and here get the southern DPS pairs
ssplat <- sapply(names(lr_south), function(i) mighty_pairs_subscript(i, lr_south, s_split))
sp_df <- data.frame(
  id1 = rep(s_split[[1]]$bung_name, each = length(s_split[[1]]$bung_name)),
  id2 = rep(s_split[[1]]$bung_name, times = length(s_split[[1]]$bung_name)),
  LLR = rowSums(ssplat, na.rm = TRUE),
  NumL = rowSums(!is.na(ssplat)),
  stringsAsFactors = FALSE
) %>%
  tbl_df %>%
  filter(id1 < id2)  %>% # remove instances of the indiv genotype to itself and also to ensure no double-counting of pairs, report them only in sorted their order version
  tidyr::separate(id1, into = c("name1", "chip1"), sep = "-") %>%
  tidyr::separate(id2, into = c("name2", "chip2"), sep = "-") %>%
  select(name1:chip2, LLR, NumL)




# now, put those into one big tidy data frame
# and only look at those pairs with >=40 loci in common
allpairs <- bind_rows(
  sp_df %>% mutate(DPS = "South"),
  np_df %>% mutate(DPS = "North")
) %>%
  filter(NumL >= 40)



apairs2 <- allpairs %>% 
  left_join(., dupies, by = c("name1" = "NMFS_DNA_ID")) %>%
  left_join(., dupies, by = c("name2" = "NMFS_DNA_ID")) %>%
  mutate(pairtype = ifelse(
    name1 == name2, "MGD", ifelse(
      (!is.na(Duplicate_Tissue.x) & name2 == Duplicate_Tissue.x) | (!is.na(Duplicate_Tissue.y) & name1 == Duplicate_Tissue.y)    ,  "SIDT", "NKS"))
  ) %>%
  filter(!(name1 == "AM000129" & name2 == "AM000147" & chip1 == 5))  # manually remove a duplicate row (dupie genotyped twice...)


# here, note that AM000120      AM000149  AM000131           AM000137   are bolluxed up.  
apairs2 %>% filter(LLR > 0, name1 != name2) %>% as.data.frame()



greensmear <- apairs2 %>% 
  filter(pairtype == "NKS", LLR < 0)
the_rest <- apairs2 %>%
  filter(!(pairtype == "NKS" & LLR < 0))
the_num_pairs <- apairs2 %>%
  group_by(DPS, pairtype) %>%
  tally %>%
  ungroup %>%
  mutate(pretty_nums = formatC(n, big.mark = ","))

set.seed(5)
g <- ggplot(greensmear, aes(x = LLR, y = pairtype, colour = pairtype)) + 
  geom_point(position = position_jitter(width = 0, height = 0.3), alpha = 0.1) +
  geom_point(data = the_rest, position = position_jitter(width = 0, height = 0.3), alpha = 1.0) +
  theme_bw() +
  xlab("Log likelihood ratio statistic") +
  ylab("Pairtype") +
  scale_color_discrete(name = "Pairtype") +
  geom_text(data = the_num_pairs, mapping = aes(label = pretty_nums), 
            x = -95, colour = "black", hjust = 1, size = 3.3) +
  facet_wrap(~ DPS, ncol = 1) +
  xlim(-100, 70) 

ggsave(g, filename = "outputs/obs-self-id-logls.pdf", width = 10, height = 3)

system("pdfcrop outputs/obs-self-id-logls.pdf")

# write out the_num_pairs for future reference
saveRDS(the_num_pairs, filename = "outputs/the_num_pairs.rds")
