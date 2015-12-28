source("R/load-packages.R")
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
# here we are going to summarize it into a table so that we know how many 
# individuals from each of the categories is on each plate and how many 
# individuals intersect across plates
# before we do that we have to add the bycatch duplicate individuals so they are
# counted up there.
dupies <- readRDS("data/meta/bycatch_IDS.rds") %>%
  filter(!is.na(Duplicate_Tissue)) %>%
  mutate(NMFS_DNA_ID = Duplicate_Tissue,
         collection_location = "Bycatch") %>%
  select(-Duplicate_Tissue)

samsheet <- read.csv("data/meta/sample_sheet.csv", stringsAsFactors = FALSE)

sswd <- plyr::rbind.fill(dupies, samsheet) %>%
  tbl_df

plate_tmp <- four_plates %>%
  group_by(Name, plate) %>%
  filter(!(Name == "NTC")) %>%
  tally() %>% 
  ungroup %>%
  mutate(NMFS_DNA_ID = Name) %>%
  select(-n, -Name) %>%
  left_join(., sswd) %>%
  mutate(collection_location = ifelse(is.na(collection_location), "Bycatch", collection_location))  # this is a kluge for those 3 bycatch fish that had screwed up IDs 


plate_loc_counts <- plate_tmp %>%
  group_by(plate, collection_location) %>%
  tally %>%
  tidyr::spread(collection_location, n) %>%
  setNames(c("Chip", "Bycatch", "Klamath", "Sacramento")) %>%
  select(Chip, Sacramento, Klamath, Bycatch)

# this is a table of how many fish from each location were on each plate
plate_loc_counts[is.na(plate_loc_counts)] <- 0


# now we want to see how many individuals are shared between each plate:
# get a vector of the IDs of the fish on each plate and then intersect
ids_on_plates <- lapply(as.character(1:4), function(x) plate_tmp$NMFS_DNA_ID[plate_tmp$plate == x] )
intsctf <- function(plate) {sapply(ids_on_plates, function(x) length(intersect(ids_on_plates[[plate]], x)))}
tmp <- lapply(1:4, function(x) intsctf(x)) %>%
  do.call(rbind, .)
tmp[lower.tri(tmp)] <- ""
colnames(tmp) <- 1:4

# so these are the number of shared samples between chips.
shared_samps <- cbind(Chip = 1:4, tmp)
colnames(shared_samps)[5] <- "4"

# now write those out to latex tables
# location counts
outf <- "outputs/chip_location_counts.tex"
cat(names(plate_loc_counts), sep = " & ", file = outf)
cat("\\\\\\hline\n", file = outf, append = TRUE)
write.table(plate_loc_counts, sep = " & ", eol = "\\\\\n", 
            row.names = F, col.names = F, file = outf, 
            quote = FALSE, append = TRUE)


outf <- "outputs/cross_chip_shared_samples.tex"
cat(colnames(shared_samps), sep = " & ", file = outf)
cat("\\\\\\hline\n", file = outf, append = TRUE)
write.table(shared_samps, sep = " & ", eol = "\\\\\n", 
            row.names = F, col.names = F, file = outf, 
            quote = FALSE, append = TRUE)






#### Step 2:  Plot the four training plates together at each locus, connecting re-genotyped individuals ####
# and we also color according to plate.  I attach the total.k to each locus (to print that out
# on the plots) before plotting.
if(!file.exists("outputs")) dir.create("outputs")
MultiChipRelativeIntensityPlot(four_plates, prefix = "outputs/plate_x_y")



#### Step 3: Read in how those four chips were scored by eye upon the fluidigm software and translate that into up to 5 clusters ####
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



# once we have done that we can make plots just like before, but we are going to
# want to put the number of clusters on there.
### These plots will go into the supplement!
MultiChipRelativeIntensityPlot(data.K, 
                               prefix = "outputs/plate_x_y_with_num_clusts_",
                               alreadyOrganized = TRUE,
                               label_with_num_clusts = TRUE,
                               lineSegAlpha = 0.45)


#### and now I want to make a figure of four example loci for the actual paper. ####
# We want the top line colored by plate and with the line segs.
# the bottom row will be colored by genotype with no line segs.
# I will just do a separate ggplot for each row and put them together in latex

# here are the loci that I want to include.  They have 2, 3, 4, and 5 clusters, respectively
foci_loci <- c("ame_21745", "ame_20992", "ame_40574", "ame_3058")

slim_dat <- data.K %>%
  ungroup %>%
  filter(as.character(assay.name) %in% foci_loci) %>%
  droplevels %>%
  mutate(assay.name = factor(assay.name, levels = foci_loci)) %>%
  mutate(assay = as.numeric(assay.name))

toprow <- MultiChipRelativeIntensityPlot(slim_dat, 
                               prefix = "outputs/junker_for_paper_by_plate",
                               alreadyOrganized = TRUE,
                               label_with_num_clusts = TRUE,
                               lineSegAlpha = 0.45,
                               n.pages = 1,
                               color.by = "plate",
                               num.columns = 4,
                               returnPlot = TRUE)

toprow <- toprow + 
  xlab("Raw intensity of dye 1") +
  ylab("Raw intensity of dye 2") +
  guides(color = guide_legend(
    title = "Chip\nor\nChip_Pair",
    override.aes = list(
      linetype = c(0, 5, 5, 5, 0, 5, 5, 0, 5, 0),
      shape = c(16, NA, NA, NA, 16, NA, NA, 16, NA, 16)))) +
  theme(axis.title = element_text(size=17),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))

ggsave(toprow, filename = "outputs/four_loci_by_plate.pdf", width = 20, height = 5)



bottomrow <- MultiChipRelativeIntensityPlot(slim_dat, 
                                            prefix = "outputs/junker_for_paper_by_geno",
                                            alreadyOrganized = TRUE,
                                            label_with_num_clusts = TRUE,
                                            lineSegAlpha = 0.45,
                                            n.pages = 1,
                                            color.by = "new.k",
                                            num.columns = 4,
                                            returnPlot = TRUE,
                                            exclude.seg=TRUE)

bottomrow <- bottomrow + 
  xlab("Raw intensity of dye 1") +
  ylab("Raw intensity of dye 2") +
  guides(color = guide_legend(
    title = "Genotype\nCategory"
    )) +
  theme(axis.title = element_text(size=17),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))

ggsave(bottomrow, filename = "outputs/four_loci_by_geno.pdf", width = 20, height = 5)

# crop the white space around those
system("pdfcrop outputs/four_loci_by_geno.pdf")
system("pdfcrop outputs/four_loci_by_plate.pdf")

#### Step 4. I make four pages of plots like the plate_x_y plots, but color things by genotype  ####
MultiChipRelativeIntensityPlot(data.K, 
                               prefix = "outputs/plate_x_y_by_final_genotype_plates_combined",
                               alreadyOrganized = TRUE,
                               color.by = "new.k",
                               exclude.seg=TRUE)




# and then make them for each plate separately
for(i in unique(data.K$plate.name)) {
  tmp <- data.K %>% filter(plate.name == i)
  MultiChipRelativeIntensityPlot(tmp, 
                                 color.by = "new.k", 
                                 alreadyOrganized = TRUE, 
                                 exclude.seg=TRUE,
                                 prefix = paste("outputs/plate_x_y_by_final_genotype_plate", i, "only", sep = "_")
  )
}


#### Make two figures for the 



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
