library(ggplot2)
library(plyr)
library(dplyr)
library(gridExtra)
library(scales)
library(grid)


source("R/chip_scoring_functions.R")

# Meta-value: regional labels
# regional_5.tbl is a combination of regional.tbl and fifth plates sample info drawn from xls file
# (AM006_ST001.xld)
region.FILE=read.csv("data/assess_5thChip_rr/annotate/regional_5.tbl", stringsAsFactors = FALSE)
colnames(region.FILE) <- c("full.name", "region")


#### Step 1: Read in the data for the four "training" plates ####
# chip names
chip <- c("1381905043",  "1381935339", "1381992037", "1381992302", "1382136064")


sturgeon.file <- lapply(chip, function(x) {
    
  read.csv(paste0("data/assess_5thChip_rr/annotate/", x, ".csv", sep = ""), 
           skip=15,
           stringsAsFactors = FALSE) %>%
    tbl_df %>%
    mutate(long_plate_name = x)
}) %>%
bind_rows(.id = "plate")



#### Step 2: Regional Assignment

data.reorg <- ReOrganizeFile(sturgeon.file)

# reassign K value: the number of cluster due to limitations of fluidigm' scoring capability
data.K <- data.reorg %>%
group_by(assay, plate) %>%
mutate(new.k = SpanningK(k, rel.dye1, rel.dye2)) %>%
ungroup() %>%
group_by(assay) %>%
mutate(total.k = max(new.k))

# merging the scoring file with the regional files
make_region_unk <- function(region,plate) {ifelse(is.na(region),paste0("unk_", plate),region)}
data.w.region <- left_join(data.K, region.FILE) %>%
  mutate(region = make_region_unk(as.character(region), plate))

# Since there are duplicated samples across four plates, I need to
# reduce any duplicated samples and pick the one that with k>0 and the highest sum rel-xy signal
data.reduce <- data.w.region %>%
group_by(assay, name) %>%
arrange(.,desc(k>0), desc(rel.dye1+rel.dye2)) %>%
filter(row_number()==1)


# creating a known-regional training dataset (North and South)
north.obs <- data.reduce %>%
filter(region %in% "North") %>%
group_by(assay.name) %>%
summarise(ct.0=sum(new.k==0),
ct.1=sum(new.k==1),
ct.2=sum(new.k==2),
ct.3=sum(new.k==3),
ct.4=sum(new.k==4),
ct.5=sum(new.k==5),
max.k =max(total.k) )  # shouldn't need max

south.obs <- data.reduce %>%
filter(region %in% "South") %>%
group_by(assay.name) %>%
summarise(ct.0=sum(new.k==0),
ct.1=sum(new.k==1),
ct.2=sum(new.k==2),
ct.3=sum(new.k==3),
ct.4=sum(new.k==4),
ct.5=sum(new.k==5),
max.k =max(total.k))

# Leave-one-out
# Take one of the "North" sample out and re-calculate its predictive label for the "North" training dataset
n.region.lou <- data.reduce %>%
  filter(region %in% "North",total.k > 0) %>%
  arrange(name, assay) %>%
  group_by(name) %>%
  summarise(prob.north =get.regional.prob(assay, new.k, region, for.region="North", take.one.out=T), 
            prob.south =get.regional.prob(assay, new.k, region, for.region="South", take.one.out=T),
            region.lab = region[1],
            n.unassign.cluster=sum(new.k==0)) %>%
  mutate(prob.being.north=prob.north/(prob.north+prob.south))

# LOU
#Take one South sample out and re-calculate its predictive label
s.region.lou <- data.reduce %>%
  filter(region %in% "South",total.k > 0) %>%
  arrange(name, assay) %>%
  group_by(name) %>%
  summarise(prob.north =get.regional.prob(assay, new.k, region, for.region="North", take.one.out=T), 
            prob.south =get.regional.prob(assay, new.k, region, for.region="South", take.one.out=T),
            region.lab = region[1],
            n.unassign.cluster=sum(new.k==0)) %>%
  mutate(prob.being.north=prob.north/(prob.north+prob.south))

# Calculate the predictive probabilities for by-catch label
region.unk.result <- data.reduce %>%
  filter(!(region %in% c("North","South")),total.k > 0) %>%
  arrange(name, assay) %>%
  group_by(name) %>%
  summarise(prob.north =get.regional.prob(assay, new.k, region, for.region="North"), 
            prob.south =get.regional.prob(assay, new.k, region, for.region="South"),
            region.lab = region[1],
            n.unassign.cluster=sum(new.k==0)) %>%
  mutate(prob.being.north=prob.north/(prob.north+prob.south))

s.region.lou$region.lab <- factor(s.region.lou$region.lab,
                                  levels=sort(c(" ",unique(s.region.lou$region.lab))))


# Function: creates a new log log-likelihood scale transformation object for the probability axis
log_split <- function(x) {
  val<-ifelse(x<0.5, log(x,100)-log(0.5,100),-log(1-x,100)+log(0.5,100) )
  val[is.infinite(val)] <- 7
  val
}
inv_log_split <- function(x ) ifelse(x<0,(100^x)*0.5, 1-2*(100^(-x))) 
log_split_trans <- function () trans_new("log_split",log_split,inv_log_split)


p<- ggplot()+
  geom_rect(data=rbind(n.region.lou, s.region.lou, region.unk.result), aes(xmax=Inf, xmin=-Inf, ymin=0.01, ymax=0.99,alpha=0.002), fill="light pink")+
  geom_point(data=n.region.lou, aes(x=as.numeric(as.factor(name)), y=prob.being.north, color=n.unassign.cluster))+
  geom_point(data=s.region.lou, aes(x=as.numeric(as.factor(name)), y=prob.being.north, color=n.unassign.cluster))+
  geom_point(data=region.unk.result, aes(x=as.numeric(as.factor(name)), y=prob.being.north, color=n.unassign.cluster))+
  xlab("Individual Sample")+
  ylab("Prob of being in North DPS")+
  theme_bw()+
  #annotation_logticks(base = 10, sides="l")+
  scale_y_continuous(trans="log_split", breaks = c(1e-16, 1e-6, 0.01, 0.1, 0.5, 0.9, 0.99, 1-(1e-6), 1-(1e-16)),
                     labels=c(1e-16, 1e-6, 0.01, 0.1, 0.5, 0.9, 0.99, 1-(1e-6), 1-(1e-16)))+
  scale_colour_gradientn(colours = rev(rainbow(7)), breaks = seq(0, 60, 10)) +
  facet_wrap(~region.lab,drop=F, scale="free_x",ncol=3)
p

g<-ggplotGrob(p)
## remove empty panels
g$grobs[names(g$grobs) %in% c("panel2", "strip_t2")] = NULL
## remove them from the layout
g$layout = g$layout[!(g$layout$name %in% c("panel-2","strip_t-2")),]
## move axis closer to panel
#g$layout[g$layout$name == "axis_l-2", c("l", "r")] = c(3,3)
grid.newpage()
grid.draw(g)

ggsave2 <- ggplot2::ggsave; body(ggsave2) <- body(ggplot2::ggsave)[-2]
ggsave2("outputs/DPS-prediction.pdf", height=17, width=20)


## pass genotype df to step 03
saveRDS(data.reduce,
        file="outputs/genotype_from_five_chips.rds",
        compress = "xz")
