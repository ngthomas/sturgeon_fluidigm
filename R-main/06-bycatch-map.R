
# quick little ggplot map to see where fish from different DPSs are caught
library(maps)
library(mapdata)
library(maptools)
library(dplyr)
library(ggplot2)

# get the DPS assignments
dps_df <- readRDS("outputs/structure_dps_assigns.rds") %>%
  rename(NMFS_DNA_ID = full.name)

# get the gsi_sim assignments too
gsi_dps <- readRDS("outputs/gsi-sim-DPS-assignments-of-bycatch.rds") %>%
  rename(NMFS_DNA_ID = IndivName)

# get the sample sheet
sam <- read.csv("data/meta/sample_sheet.csv", stringsAsFactors = FALSE) %>% tbl_df


# only evaluate these if you have the private meta-data from the Region
if(PRIVATE_ACCESS == TRUE) {
  # get the meta data (not on repo)
  tmp <- readxl::read_excel("data/meta/private-green-sturgeon-reconciled.xlsx", sheet = 1, skip = 0) 
  region_meta <- tmp[-1, names(tmp) != ""] # messy excel file---have to rip of extra unnamed columns, and the first row has nada in it
  
  
  #### Some Reporting and Output Stuff First  ####
  # now, before we do anything else we are going to join to each of these assignents the reconciled meta
  # data file from the region, so carlos can send it back to them.
  Results4Reg <- gsi_dps %>%
    left_join(dps_df) %>%
    rename(DPS_structure = DPS,
           SouthernDPS_prob_gsi_sim = SouthernDPS,
           NorthernDPS_prob_gsi_sim = NorthernDPS) %>%
    full_join(region_meta, .)
  
  # and before we send that back we will want to make a note of the duplicated
  # samples, because the duplicate tisssues don't get DPS assignments (because
  # their other half has it.)
  bycid_rev <- readRDS("data/meta/bycatch_IDS.rds") %>%
    filter(!is.na(Duplicate_Tissue)) %>%
    rename(Duplicated_Tissue_Of = NMFS_DNA_ID,
           NMFS_DNA_ID = Duplicate_Tissue)
  
  Res4Reg_final <- left_join(Results4Reg, bycid_rev)
  
  write.csv(Res4Reg_final, file = "outputs/private-dps-assignments-with-meta-data-for-region.csv", row.names = FALSE)
  
  
  ## In our first pass through these data we found a handful of matching samples
  ## amongst the bycatch that we were not told by the Region were duplicately sampled
  ## tissues.  I used the following lines to grab the meta-data for those individuals
  ## and discover that those pairs were all sampled at the exact same time and were the
  ## same size, so clearly they were just duplicate tissues from the same individual and
  ## I updated that accordingly.
  if(FALSE) {
    ## Now, while we are at it, we might as well get a summary of the 
    ## matching samples that were not known as matching beforehand.
    nks_matchers <- readRDS("outputs/nks_matching_pairs") 
    
    ## make a data frame that has the matching pair followed by three NAs so that we 
    ## can left_join other stuff onto them
    nksv <- nks_matchers %>% 
      select(name1, name2) %>%
      as.matrix %>%
      t() %>% 
      rbind(NA, NA, NA) %>%
      as.vector
    
    nksdf <- data.frame(NMFS_DNA_ID = nksv, stringsAsFactors = FALSE) %>%
      left_join(region_meta %>% filter(!is.na(NMFS_DNA_ID)))
    
    write.csv(nksdf, file = "outputs/previously-unannounced-matching-genotypes-for-region.csv", 
              na = "", row.names = FALSE)
  }
  
  
  #### Now, start making maps ####
  
  # now put gsi_dps together with the region's lat-long data
  DF <- region_meta %>%
    left_join(gsi_dps, .) %>%
    filter(!is.na(RETRIEVE_LAT), !is.na(RETRIEVE_LONG))  %>%
    select(NMFS_DNA_ID:DPS_gsi_sim, DISSECTION_BARCODE_NMFS, RETRIEVE_LAT, RETRIEVE_LONG) %>%
    rename(DPS = DPS_gsi_sim)
  
  # here I am going to make a data frame of permuted lats and longs (within DPS)
  df2 <- DF %>% 
    group_by(DPS) %>%
    mutate(perm_idx = sample(1:n()),
           permuted_lat = RETRIEVE_LAT[perm_idx],
           permuted_long = RETRIEVE_LONG[perm_idx])
  
  # and I will write that out
  df2 %>%
    ungroup() %>%
    select(NMFS_DNA_ID, permuted_lat, permuted_long) %>%
    saveRDS(file = "data/meta/permuted_lat_longs.rds")
} else {  # If you don't have PRIVATE_ACCESS get the permuted lat-longs from here
  permy <- readRDS("data/meta/permuted_lat_longs.rds")
  DF <- left_join(gsi_dps, sam) %>%
    left_join(permy) %>%
    rename(RETRIEVE_LAT = permuted_lat,
           RETRIEVE_LONG = permuted_long,
           DPS = DPS_gsi_sim) %>%
    filter(!is.na(RETRIEVE_LAT) & !is.na(RETRIEVE_LONG))
}
# now we just have to map those
# look at our limits:
range(DF$RETRIEVE_LAT)
range(DF$RETRIEVE_LONG)



# let us get some hires map data. see http://eriqande.github.io/2014/12/17/compare-map-resolutions.html
# This is no longer executed, because I have saved sf1 into an RDS with the repository
# but I leave the code here so someone can see what I did in this step to make that RDS.
if(FALSE) {
  if (!rgeosStatus()) gpclibPermit()
  gshhs.f.b <- "/Users/eriq/Maps/gshhg-bin-2.3.3/gshhs_f.b"
  sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-128, -120), ylim = c(36, 48)) %>%
    fortify()
  saveRDS(sf1, file = "data/maps/RgshhMap_sf1.rds", compress = "xz")
}
sf1 <- readRDS("data/maps/RgshhMap_sf1.rds")

# plot up a map at full large scale  
# make a column of jittered lats and longs.  We do this instead of geom_jitter because
# we don't want things to jitter onto land and we want the jitter to be consistent
# between the different "zoom" levels of each plot.
set.seed(5)
DF2 <- DF %>%
  mutate(jitlat = RETRIEVE_LAT + runif(n = nrow(DF), min = -0.01, max = 0.01),
         jitlong = RETRIEVE_LONG - runif(n = nrow(DF), min = 0, max = 0.02))


g <- ggplot() + geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey80") + 
  coord_fixed(1.3) + 
  geom_point(data = DF2, mapping = aes(x = jitlong, 
                                       y = jitlat, 
                                       colour = DPS), 
              alpha = 0.6) +
  scale_colour_manual(values = c(north = "red", south = "blue")) +
  xlab("Latitude") +
  ylab("Longitude") +
  theme_bw()


# now, we want to make two sub-figures that are blow-ups of the two
# obvious regions of interest.  We will make those without the color guide
# then print them out on the same pages via grid
long_skinny <- g + theme(legend.position="top")

wa_coast <- g + 
  coord_fixed(ratio = 1.3, xlim = c(-122.5, -124.5), ylim = c(45.75, 47)) +
  theme(legend.position="none")

ca_coast <- g + 
  coord_fixed(ratio = 1.3, xlim = c(-121.5, -123.5), ylim = c(37, 38.25)) +
  theme(legend.position="none")


#### Define a function multiplot ####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
# ... are extra arguments passed to grid.layout
multiplot <- function(plotlist=NULL, file, cols=1, layout=NULL, ...) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- plotlist
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), ...)))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#### Now plot stuff ####
pdf(file = "outputs/bycatch_map.pdf", width = 15, height = 12)
multiplot(plotlist = list(long_skinny, wa_coast, ca_coast), 
          layout = matrix(c(1,1,2,3), nrow=2),
          widths = unit(c(1,1), "null"))
dev.off()
system("cd outputs; pdfcrop bycatch_map.pdf")


