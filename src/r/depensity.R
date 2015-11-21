# This function makes relative intensity plot from a combined CSV plate file

MakeIntensityPlot <- function (combined.FILE, prefix="plate_x_y", n.pages = 4, self.exclude=FALSE) {
  

ReOrganizeFile <- function (test.FILE) {
  clean.file <- test.FILE[test.FILE$Type=="Unknown",]
  clean.file$Name <- factor(clean.file$Name)
  clean.file$Final <-factor(clean.file$Final, levels=c("YY","XY", "XX"))
  clean.file$Assay <- factor(clean.file$Assay, levels = unique(clean.file$Assay))
  
  core.data <- data.frame(
    plate = clean.file$plate,
    assay= as.numeric(clean.file$Assay),
    name = as.numeric(clean.file$Name),
    rel.dye1 = clean.file$Allele.X.1,
    rel.dye2 = clean.file$Allele.Y.1,
    k=as.numeric(clean.file$Final),
    full.name=clean.file$Name,
    assay.name=clean.file$Assay
  )
  
  core.data <- arrange(core.data, assay, name)
  core.data
}

core.data <- ReOrganizeFile(combined.FILE)

core.pdata<-tbl_df(core.data)

repeat.data <- core.pdata %>%
  group_by(assay,name) %>%
  summarise(count = n()) %>%
  filter(count > 1) %>%
  select(assay, name) %>%
  inner_join(core.pdata) 

list.seg<-ddply(repeat.data, .(assay, name), function(x) {
  num.rep <- dim(x)[1]
  all.combn <- combn(1:num.rep, 2)  
  seg.pt <- t(apply(all.combn, 2, function(i) {
    c(x$rel.dye1[i[1]], 
      x$rel.dye2[i[1]],
      x$rel.dye1[i[2]],
      x$rel.dye2[i[2]],
      x$plate[i[1]],
      x$plate[i[2]],
      paste(sort(c(x$plate[i[1]], 
                   x$plate[i[2]])), collapse="_"),
      as.character(x$assay.name[1])
    )}))
  
  colnames(seg.pt)<- c("x.start", "y.start", "x.end", "y.end", "plate.1","plate.2","plate.pair", "assay.name")
  seg.pt
})

if (self.exclude) list.seg <- list.seg %>% filter(as.integer(plate.1)!=as.integer(plate.2))


list.seg$x.start <- as.numeric(as.character(list.seg$x.start))
list.seg$y.start <- as.numeric(as.character(list.seg$y.start))
list.seg$x.end <- as.numeric(as.character(list.seg$x.end))
list.seg$y.end <- as.numeric(as.character(list.seg$y.end))
list.seg$assay.name <- factor(list.seg$assay.name,levels(core.data$assay.name))

n.loci <- length(unique(core.data$assay.name))
for (i in 1:(n.pages)) {
  
  min.intv <- round(n.loci/n.pages)*(i-1)
  max.intv <- min.intv + round(n.loci/n.pages)
  g <- ggplot(data=core.data %>% 
                filter(as.integer(assay.name)>min.intv, as.integer(assay.name) <= max.intv) %>% 
                droplevels(), 
              aes(x=rel.dye1, y=rel.dye2, color=factor(plate))) +
    geom_point(alpha=0.7)
  
  seg.data <- list.seg %>% 
                   filter(as.integer(assay.name)>min.intv, as.integer(assay.name) <= max.intv) %>%
                   droplevels()
  if(nrow(seg.data)>0){
    g<- g + geom_segment(data=seg.data,
                 aes(x=x.start, y=y.start, xend=x.end, yend=y.end, color=plate.pair), linetype=5)
    }
    g<- g+ facet_wrap(~assay.name, ncol = 6 )+
    theme_bw()
  
  ggsave(paste0(prefix,i,".pdf"),g, width = 34, height = 22)
}
}



# Since the fluidigm software allows users to labeled points into three groups (four if we also use the "invalid" label), 
# So for any loci contains more than three clusters, I reuse some of the labels in rotations. 
# This function corrects the current k (max: 3) value to the actual cluster value. 
# For now, the largest number of clusters this able to correct up to is 8.

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

ReOrganizeFile <- function (test.FILE) {
  clean.file <- test.FILE[test.FILE$Type=="Unknown",] ## remove NTC 
  clean.file$Name <- factor(clean.file$Name)
  clean.file$Final <-factor(clean.file$Final, levels=c("No Call","XX","XY", "YY", "Invalid"))
  clean.file$Assay <- factor(clean.file$Assay, levels = unique(clean.file$Assay))
  
  core.data <- data.frame(
    plate = clean.file$plate,
    assay= as.numeric(clean.file$Assay),
    name = as.numeric(clean.file$Name),
    rel.dye1 = clean.file$Allele.X.1,
    rel.dye2 = clean.file$Allele.Y.1,
    k=as.numeric(clean.file$Final)-1,
    full.name=clean.file$Name,
    assay.name=clean.file$Assay,
    plate.name=clean.file$long_plate_name
  )
  
  core.data <- arrange(core.data, assay, name)
  core.data %>% tbl_df
}


ReOrganizeFile_DP <- function (test.FILE) {
  test.FILE %>%
    filter(Type=="Unknown") %>% ## remove NTC 
    mutate(
      Name = factor(Name),
      Final = factor(Final, levels=c("No Call","XX","XY", "YY", "Invalid")),
      Assay = factor(Assay, levels = unique(Assay)),
      plate = as.integer(plate),
      assay = as.numeric(Assay),
      name = as.numeric(Name),
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

