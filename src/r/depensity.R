MakeIntensityPlot <- function (combined.FILE, prefix="plate_x_y", n.groups = 4) {
  

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
      paste(sort(c(x$plate[i[1]], 
                   x$plate[i[2]])), collapse="_"),
      as.character(x$assay.name[1])
    )}))
  
  colnames(seg.pt)<- c("x.start", "y.start", "x.end", "y.end", "plate.pair", "assay.name")
  seg.pt
})

list.seg$x.start <- as.numeric(as.character(list.seg$x.start))
list.seg$y.start <- as.numeric(as.character(list.seg$y.start))
list.seg$x.end <- as.numeric(as.character(list.seg$x.end))
list.seg$y.end <- as.numeric(as.character(list.seg$y.end))
list.seg$assay.name <- factor(list.seg$assay.name,levels(core.data$assay.name))

n.loci <- length(unique(core.data$assay.name))
for (i in 1:(n.groups)) {
  
  min.intv <- round(n.loci/n.groups)*(i-1)
  max.intv <- min.intv + round(n.loci/n.groups)
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
