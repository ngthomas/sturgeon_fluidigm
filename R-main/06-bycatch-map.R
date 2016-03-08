
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

# get the meta data (not on repo)
tmp <- readxl::read_excel("data/meta/private-green-sturgeon-reconciled.xlsx", sheet = 1, skip = 0) 
region_meta <- tmp[-1, names(tmp) != ""] # messy excel file---have to rip of extra unnamed columns, and the first row has nada in it



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



# now put gsi_dps together with the region's lat-long data
DF <- region_meta %>%
  left_join(gsi_dps, .) %>%
  filter(!is.na(RETRIEVE_LAT), !is.na(RETRIEVE_LONG))  %>%
  select(NMFS_DNA_ID:DPS_gsi_sim, DISSECTION_BARCODE_NMFS, RETRIEVE_LAT, RETRIEVE_LONG) %>%
  rename(DPS = DPS_gsi_sim)


# now we just have to map those
# look at our limits:
range(DF$RETRIEVE_LAT)
range(DF$RETRIEVE_LONG)



# let us get some hires map data. see http://eriqande.github.io/2014/12/17/compare-map-resolutions.html
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "/Users/eriq/Maps/gshhg-bin-2.3.3/gshhs_f.b"
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-128, -120), ylim = c(36, 48)) %>%
  fortify()


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


