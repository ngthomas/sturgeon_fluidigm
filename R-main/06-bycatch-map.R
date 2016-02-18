
# quick little ggplot map to see where fish from different DPSs are caught
library(maps)
library(mapdata)
library(maptools)
library(dplyr)
library(ggplot2)

# get the DPS assignments
dps_df <- readRDS("outputs/structure_dps_assigns.rds") %>%
  rename(NMFS_DNA_ID = full.name)

# get the sample sheet
sam <- read.csv("data/meta/sample_sheet.csv", stringsAsFactors = FALSE) %>% tbl_df

# get the lat long data (not on repo)
lat_long <- readxl::read_excel("data/meta/private.xls", sheet = 2, skip = 4) %>%
  select(TISSUE_SAMPLE_ID, NMFS_DNA_ID, RETRIEVE_LAT, RETRIEVE_LONG)


# while we are at it, get all the meta data for the bycatch:
all_meta <- readxl::read_excel("data/meta/private.xls", sheet = 2, skip = 4)
AssigAll <- sam %>%
  filter(collection_location == "Bycatch") %>%
  select(NMFS_DNA_ID) %>%
  left_join(dps_df) %>%
  left_join(all_meta)
write.csv(AssigAll, file = "outputs/bycaught_sturgeon_with_DPS.csv")

# now put these together to get just the bycatch samples
DF <- sam %>%
  filter(collection_location == "Bycatch") %>%
  select(NMFS_DNA_ID) %>%
  left_join(dps_df) %>%
  left_join(lat_long) %>%
  filter(!is.na(RETRIEVE_LAT), !is.na(RETRIEVE_LONG))  

# note, we don't have lat-longs for the recently received samples

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
set.seed(7)
g <- ggplot() + geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey80") + 
  coord_fixed(1.3) + 
  geom_point(data = DF, mapping = aes(x = RETRIEVE_LONG, 
                                       y = RETRIEVE_LAT, 
                                       colour = DPS), 
              alpha = 0.7,
             position=position_jitter(w = 0.05, h = 0.05)) +
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
