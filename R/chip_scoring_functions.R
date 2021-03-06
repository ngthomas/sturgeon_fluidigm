

#### HELPER FUNCTIONS ####
#' read a fluidigm detailed table results and lop off the garbage at the bottom
#' 
#' this skips the preamble info and chucks the output
#' past the lines with "Dose Meter Reading Data" if that line is present.
#' If it isn't, then it just returns everything. Wraps output in tbl_df
#' @param path  path to the file
read_fluidigm_detailed_csv <- function(path) {
  tmp <- tbl_df(read.csv(path, skip = 15, stringsAsFactors = FALSE))
  
  DoseMeterLine <- which(tmp$ID == "Dose Meter Reading Data")
  
  if(length(DoseMeterLine) == 0) {
    return(tmp)
  }
  bottom_line <- DoseMeterLine - 1  # remove the bottom, irrelevant part of each file.
  
  tmp[1:bottom_line,]
}




#' Organize the combined data.  Put on appropriate factor levels, etc
#' @param x  the combined_file (like that created in Step 1 of 01-develop-chip-scoring.R)
#' 
#' 
ReOrganizeFile <- function (x) {
  x %>%
    filter(Type=="Unknown") %>% ## remove NTC 
    mutate(
      Name = Name,
      Final = factor(Final, levels=c("No Call","XX","XY", "YY", "Invalid")),
      Assay = factor(Assay, levels = unique(Assay)),
      plate = as.integer(plate),
      assay = as.numeric(Assay),
      name = as.numeric(factor(Name)),
      k = as.numeric(Final) - 1
    ) %>%
    rename(
      rel.dye1 = Allele.X.1,
      rel.dye2 = Allele.Y.1,
      full.name = Name,
      assay.name = Assay,
      plate.name = long_plate_name
    ) %>%
    arrange(assay, name) %>%
    select(plate, assay, name, rel.dye1, rel.dye2, k, full.name, assay.name, plate.name)
}


#### MULTI-CHIP RELATIVE INTENSITY PLOTTING FUNCTION ####
#' Plot relative intensities from multiple chips together on same axes
#' 
#' This creates n.pages pages of plots.  Each one has a grid of facets
#' on it, each one representing a locus.  Each point shows the relative
#' intensity of a single individual, colored according to the plate it was
#' on.  If the same individual was genotyped multiple times, then a line 
#' segment is drawn between its different scores.  The color of the segment
#' shows which plates the individual was scored on.   
#' @param DF the long format data frame that has all the intensities.  It
#' should have the columns "plate", "Assay", "Name", "rel.dye1", "rel.dye2",
#' and "Final".
#' @param prefix  The filename prefix to be given the pdf output. Can 
#' include a directory prefix like "outputs/plate_x_y"
#' @param n.pages  The number of pages to break the plots over.  Default
#' is 4 which gives 24 loci per page
#' @param num.columns  How many columns per page.  Default is 6.
#' @param self.exclude A boolean option for omitting lines drawn between identical
#' samples within the same plate. Default is False.
#' @param color.by 
#' @param exclude.seg if it is true, segments are drawn to show how the same sample 
#' labels are traveled within plates and plates 
MultiChipRelativeIntensityPlot <- function (DF, 
                                            prefix="plate_x_y", 
                                            n.pages = 4,
                                            num.columns = 6,
                                            self.exclude = FALSE,
                                            color.by = "plate",
                                            exclude.seg = FALSE,
                                            alreadyOrganized = FALSE,
                                            label_with_num_clusts = FALSE,
                                            lineSegAlpha = 1.0,
                                            returnPlot = FALSE,
                                            LegendForCalledGenotypes = FALSE) {
  
  if (!alreadyOrganized) {
    DF <- ReOrganizeFile(DF) }
  
  # make a data frame of all those individuals that appear on more than one plate
  # it turns out that there are 56 individuals that were regenotyped on more
  # than one plate.
  repeat.data <- DF %>%
    group_by(assay,name) %>%
    summarise(count = n()) %>%
    filter(count > 1) %>%
    select(assay, name) %>%
    inner_join(DF) 
  
  
  # here, Thomas wrote some code to get a data 
  # frame of the segments that should be connected between different points
  # that represent re-genotyped individuals' values on the final plot
  list.seg<-ddply(repeat.data, .(assay.name, name), function(x) {
    num.rep <- dim(x)[1]
    all.combn <- utils::combn(1:num.rep, 2)  
    seg.pt <- t(apply(all.combn, 2, function(i) {
      c(as.numeric(x$rel.dye1[i[1]]), 
        x$rel.dye2[i[1]],
        x$rel.dye1[i[2]],
        x$rel.dye2[i[2]],
        paste(sort(c(x$plate[i[1]], 
                     x$plate[i[2]])), collapse="_")
      )}))
    
    colnames(seg.pt)<- c("x.start", "y.start", "x.end", "y.end", "plate.pair")
    seg.pt
  }) %>%
    tbl_df
  
  
  
  # if self exclude is true, identical sample segments are not drawn if those samples
  # are from the same plate
  if (self.exclude) {
    # add on columns of plate.1 and plate.2 (I don't know where these were supposed to come from...Thomas had a call out to them.)
    list.seg <- list.seg %>%
      mutate(
        plate.1 = str_split_fixed(plate.pair, "_", 2)[,1],
        plate.2 = str_split_fixed(plate.pair, "_", 2)[,2]
      )
    
    list.seg <- list.seg %>% filter(as.integer(plate.1) != as.integer(plate.2)) %>%
      select(-plate.1, -plate.2)
  }
  

  # now, because the apply function returns a matrix, we have to coerce all these
  # things back to numeric or factors with correct levels.  We might be able to 
  # avoid this if we adply-ed rather than apply-ed.  But it's done and working now!
  list.seg$x.start <- as.numeric(as.character(list.seg$x.start))
  list.seg$y.start <- as.numeric(as.character(list.seg$y.start))
  list.seg$x.end <- as.numeric(as.character(list.seg$x.end))
  list.seg$y.end <- as.numeric(as.character(list.seg$y.end))
  

  # Do the ggplotting, breaking over n.pages pages
  n.loci <- length(unique(DF$assay))
  for (i in 1:(n.pages)) {
    
    min.intv <- round(n.loci/n.pages)*(i-1)
    max.intv <- min.intv + round(n.loci/n.pages)
    
    gdattt <- DF %>% 
      filter(as.integer(assay)>min.intv, as.integer(assay) <= max.intv) %>% 
      droplevels()
    g <- ggplot(data = gdattt, 
                aes_string(x="rel.dye1", y="rel.dye2", color=paste0("factor(", color.by,")")))

    g <- g + geom_point(alpha=0.7)
    
    if (color.by == "new.k" || label_with_num_clusts == TRUE) {
      g <- g + geom_text(aes(label=total.k, x=1, y=1), color="black")
    }
    

    seg.data <- list.seg %>% 
      filter(as.integer(assay.name)>min.intv, as.integer(assay.name) <= max.intv) %>%
      droplevels()

    if (!exclude.seg) {
      if(nrow(seg.data)>0) {
        g<- g + geom_segment(data=seg.data,
                             aes(x=x.start, 
                                 y=y.start, 
                                 xend=x.end, 
                                 yend=y.end, 
                                 color=plate.pair),
                             linetype=5,
                             alpha = lineSegAlpha)
      }}
    
    # color-blind-friendly palette: http://www.cookbook-r.com/Graphs/Colors_%28ggplot2%29/
    # though eric added some more colors to it and put black near the end cuz the black was just fugly.
    cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999", "#9933FF", "#33FF99")
    
    if (color.by == "new.k") {
      cbPalette[1] <- "light grey"
    }
    
    g<- g + facet_wrap(~assay.name, ncol = num.columns )+
      theme_bw() + 
      scale_color_manual(values=cbPalette)
    
    g <- g + 
      xlab("Raw intensity of dye 1") +
      ylab("Raw intensity of dye 2") +
      theme(axis.title = element_text(size=17),
            axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 13))
    
    
    if(LegendForCalledGenotypes == FALSE) {
      g <- g + 
        guides(color = guide_legend(
          title = "Chip\nor\nChip_Pair",
          override.aes = list(
            linetype = c(0, 5, 5, 5, 0, 5, 5, 0, 5, 0),
            shape = c(16, NA, NA, NA, 16, NA, NA, 16, NA, 16))))
    } else {
      g <- g + 
        guides(color = guide_legend(
          title = "Genotype\nCategory"
        )) 
    }
    
    
    ggsave(paste0(prefix,i,".pdf"),g, width = 34, height = 22)
  }
  
  if(returnPlot == TRUE) return(g)
}


#### ERIC HAS WORKED THINGS OVER AND CLEANED THEM UP DOWN TO THIS POINT ####

#### After manually scoring, these functions are used to name the genotype clusters at each locus ####
# Since the fluidigm only allow things to be labeled into three groups 
# (four if using "invalid" label), for loci with more than three or four clusters, I reuse the same
# label in full rotation. I need to correct the k value to match up the actual cluster ID. For now
# the max k:8

#### CORRECTING the K - CLUSTER VALUE ####
#' The fluidigm software calls a maximum number of three clusters. This function
#' increases the ceiling of k if there are more than three hand-label clusters 
#' @param orig.k the k or cluster number that the software assigned
#' @param rel.dye1 the relative intensity level of the first dye  
#' @param rel.dye2 the relative intensity level of the second dye

SpanningK <- function(orig.k ,rel.dye1, rel.dye2) 
{
  k <- orig.k
  
  if(max(k) > 3) {
    sort.k.ind <- order(desc(atan(rel.dye2/rel.dye1)))
    sort.k <- orig.k[sort.k.ind]
    
    sort.k.gt0.ind <- sort.k.ind[which(sort.k>0)]
    sort.k.gt0 <- orig.k[sort.k.gt0.ind]
    
    track.counter <- 0
    inc.switch <- 0
    for ( i in 2:length(sort.k.gt0)) {
      if(sort.k.gt0[i-1] - sort.k.gt0[i] >= 2) {
        inc.switch <- 1
        track.counter <- sort.k.gt0[i]
      }
      
      if (inc.switch == 1 && abs(sort.k.gt0[i] - track.counter) %in% c(0,1)) {
        k[sort.k.gt0.ind[i]] <- k[sort.k.gt0.ind[i]] + 4
        track.counter <- max(sort.k.gt0[i], track.counter)
      }
    }}
  k
}

#### Calculating Predictive Poster of the Leave-one-out sample ####
#' This function removes one of the dataset.
#' @param orig.k the k or cluster number that the software assigned
#' @param data
#' @param data.rm 


cal.pred.one.y <- function(data, data.rm){
  
  y.tbl <- as.matrix(data[,3:7])-data.rm
  okay.loci <- data$max.k > 0
  occupy.matrix <- okay.loci*t(sapply(1:nrow(data), function(i){
    blank<-rep(0,5)
    blank[1:data$max.k[i]] <- 1
    blank
  }))
  
  #hyperparam value : for now, assume uniform
  hyper.alpha <- ifelse(data$max.k >0, 1/data$max.k, 0)
  
  #data$max.
  new.alpha<- (hyper.alpha + y.tbl)*occupy.matrix
  
  ## This is the longer way to express predictive Polya-Eggenberg urn scheme
  
  #lg.new.alpha <- lgamma(new.alpha)
  #lg.new.alpha[occupy.matrix==0]<-0
  #p1<- rowSums(lg.new.alpha)
  #p2 <- lgamma(rowSums(new.alpha))
  #p2[okay.loci==0]<-0
  #p3 <- lgamma(rowSums(new.alpha)+1)
  #p4 <- (lgamma(new.alpha+1)+p1-lg.new.alpha)*occupy.matrix  
  #exp((p4+(p2-p1-p3))*occupy.matrix)
  
  # the quicker way is this
  prob.new <- new.alpha/rowSums(new.alpha)
  prob.new[occupy.matrix==0]<-0
  prob.new
}

get.regional.prob <- function(assay, new.k, region, for.region="North", take.one.out=F, n.loci=96, n.k=5) {
  
  #n.loci<- dim(prob.matrix.north)[1]
  #n.k <- dim(prob.matrix.north)[2]
  region <- region[1]
  
  rm.matrix<-matrix(0, ncol=n.k, nrow=n.loci)
  
  if(take.one.out && for.region==region) {
    rm.matrix[cbind(assay,new.k)] <- 1
  }
  
  if(for.region=="North") {
    prob.matrix.north <- cal.pred.one.y(north.obs, data.rm<-rm.matrix)
    exp(sum(log(prob.matrix.north[cbind(assay, new.k)])))
  }
  else {
    prob.matrix.south <- cal.pred.one.y(south.obs, data.rm<-rm.matrix)
    exp(sum(log(prob.matrix.south[cbind(assay, new.k)])))
  }
}

