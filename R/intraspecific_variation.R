## Accounting for intraspecific variation 

source("~/Desktop/Phylomethods/PhyloPlasticitySals/R/packages.R")  
setwd("/Users/morganfleming/Desktop/Phylomethods/PhyloPlasticitySals/data")
sal.tree<-read.nexus("~/Desktop/Phylomethods/PhyloPlasticitySals/data/amphibian_tree.nex")
sal.data<-read.csv("~/Desktop/Phylomethods/PhyloPlasticitySals/data/salamanders.csv")
dim(sal.data)
unique(sal.data$Genus)

#data validation and cleaning 
sals.dat$Genus[which(sal.data$Genus=='Leurognathus')]<-'Desmognathus marmoratus' 
sals.dat <- sal.data
sals.dat$gen2<-gsub(' ', '_',sals.dat$Genus)
sals.dat$gen2[which(sals.dat$gen2=='Desmognathus_monticola_')] <- 'Desmognathus_monticola'
sals.dat$gen2[which(sals.dat$gen2=='Plethodon_cinereus_')] <- 'Plethodon_cinereus'
sals.dat$gen2[which(sals.dat$gen2=='Eurycea_bislineata_')] <- 'Eurycea_bislineata'
sals.dat$gen2[which(sals.dat$gen2=='Plethodon_glutinosis')] <- 'Plethodon_glutinosus'
sals.dat<-sals.dat[!(sals.dat$gen2=="Desmognathus_ochrophaeus"),]

unique(sals.dat$Genus)


PruneTree <- function( phy, data){
  tk <- which(phy$tip.label%in%data[,22]=='TRUE')
  ready.tree <- keep.tip(phy, tk) ## 14 tips
  plot(ready.tree) ## okay looks good.
  return(ready.tree)
}

tree.go <- PruneTree(sal.tree,sals.dat); print(tree.go)


#function for avg body size and variance to use for different subsets of species samples
View(sal.dat)
body.size.species.sub <- function(data){
  calc.data<-matrix(nrow=14, ncol=4)
  rownames(calc.data) <- unique(data[,6])
  colnames(calc.data) <- c('avg.size', 'variance')
  for (i in 1:length(unique(data[,6]))){
    species <-unique(data[,6])[i]
    data.per.species <- data[which(data[,6]==species),]
    female.only<-data.per.species[which(data.per.species$Sex=='female'), ]
    female.avg <-round(mean(female.only$SVL))
    female.size.variance <- sd(female.only$SVL)
    male.only<-specdat[which(specdat$Sex=='male'), ]
    male.avg <-round(mean(male.only$SVL))
    male.size.variance <- sd(female.only$SVL)
    calc.data[i,1]<-female.avg
    calc.data[i,2]<-female.size.variance 
    calc.data[i,3]<-male.avg 
    calc.data[i,4]<-male.size.variance 
  }
  return(calc.data)
}
body.size.species.sub(sal.dat)


#Function for pulling random samples 
Get.subset.sal <- function (data){
  if (sals.dat[which(sals.dat$Sex=='female')]){
    female <- sample(1:nrow(sals.dat.size.female$SVL), 5, replace = FALSE)
  }
  sal.sub <- subset(sals.dat,Sex=='female', select = SVL)
  summary(sals.sub[,1:7]) 
sub <- sample(1:nrow(sals.dat.size.female$SVL), 5, replace = FALSE)
sub
#repeat for species. change size for each BM 



###
#data= sals.dat
row.names(sals.dat)
SVL<-sals.dat$SVL
SVL<-sals.dat$Elevation
names(sals.dat$Genus)<-row.names(sals.dat)
names(SVL)<-NULL
names(SVL)<-row.names(dat)


# Brownian motion
BM1 <- geiger::fitContinuous(tree.go, Elev.class, model="BM") 
print(BM1)

OU1 <- fitContinuous(tree.go, Elev.class, model="OU")
print(OU1)
par(mfcol=(c(1,2)))
plot(tree.go, show.tip.label=F)

ou.tree <- rescale(tree.go, model="OU", OU1$opt$alpha)
plot(ou.tree, show.tip.label = F)

AIC.BM1 <- BM1$opt$aic
AIC.OU1 <- OU1$opt$aic

#Elev interval class -- high, med, low
P.cinereus <- Adults[Adults$Species=="P.c.cinereus",] #just subsetting the species of interest
range(sals.dat$Elevation, na.rm = TRUE)
plot(sals.dat$Genus,sals.dat$Elevation)


newcol<-matrix(nrow=13, ncol=1)
rownames(newcol) <- unique(disdatsal$gen2)
colnames(newcol) <- c('size_class')
SizeClass <- function(x) {
  for (i in 1:length(cleaned.data)){
    spec<-cleaned.data[i]
    newcol[i]<-ifelse(spec$SVLF>50,"big","small")
    return(newcol)
    }
}




















##########################################################


library(tidyverse)
library(tibble)
ORGAN <- read_csv("Organ_march_21.csv", quoted_na = TRUE, col_names = TRUE)

#Selecting desired columns and subsetting for only the adults (Male/Female)
sals<- ORGAN %>%
  select("Species","Genus","Species.Name","Sex","SVL","Total","Habitat","Stream.Distance","Microhabitat","Date.Collected","Month","Elevation","Year Collected","Study.Areas")%>%
  filter(Sex == "female" | Sex == "male")%>%
  filter(Species %in% c("D.fuscus","D.monticola","D.o.carolinensis","D.quadramaculatus","D.wrighti","E.b.wilderae","P.c.cinereus","P.glutinosus","P.j.metcalfi","P.richmondi","P.welleri","P.yonahlossee","Pseudotriton.r.nitidus"))%>%
  mutate(Species = as.factor(Species),
         Sex = as.factor(Sex),
         Genus = as.factor(Genus),
         Species.Name = as.factor(Species.Name),
         Date.Collected = as.Date(Date.Collected, format="%m/%d/%Y"), 
         Month = as.factor(Month), 
         Habitat = as.factor(Habitat),
         Study.Areas = as.factor(Study.Areas)
         
  ) 

#--------------------------- Transformations and Prep -------------------------------#
hist(sals$Elevation)
ELEVATION_ <- log(sals$Elevation)
hist(ELEVATION_)

hist(sals$SVL)
SVL_ <- log(sals$SVL)
hist(SVL_)

qqplot(ELEVATION_,SVL_)
plot(ELEVATION_,SVL_)

summary(sals)
str(sals)
levels(sals$Species)




#--------------------------------- Species Ranges ----------------------------------#

ggplot(data = sals, mapping = aes(x = Species, y = Elevation)) + 
  geom_boxplot()+
  ggtitle("Species Ranges")+
  theme(plot.title = element_text(hjust = 0.5))

    #classify elev ranges with ~ equal sample size/species rich

#----------------------------- SVL across all Species ---------------------------------# 

summary((lm(SVL_ ~ Species, data = sals)))

ggplot(sals, aes(Species,SVL_))+
  geom_boxplot(alpha=0.5)+
  ggtitle("SVL across Species")+
  ylab("SVL")+ 
  xlab("Species")+
  geom_smooth(method='lm', se=FALSE)

Species <-aov(SVL_ ~ Species, data = sals)
TUKEY<-TukeyHSD(Species) #multiple comparison

library(multcompView)
AOV<- aov(SVL_ ~ Species, data = sals)
multcompLetters4(AOV, TUKEY)
# obvious species differences 


#-------------------------------SVL across Genus -----------------------------------#

summary((lm(SVL_ ~ Genus, data = sals)))

ggplot(sals, aes(Genus,SVL_))+
  geom_boxplot(alpha=0.5)+
  ggtitle("SVL across Genus")+
  ylab("SVL")+ 
  xlab("Genus")+
  geom_smooth(method='lm', se=FALSE)

mean.diff(x=SVL_,y=sals$Genus)
# obvious size distinctions across Genus


#----------------------------- SVL across ELEVATIONS --------------------------------#

summary(lm(SVL_~ELEVATION_, data=sals))

ggplot(sals, aes(ELEVATION_,SVL_))+
  geom_point(alpha=0.5)+
  ggtitle("SVL across Elevation_")+
  ylab("SVL")+ 
  xlab("Elevation_")+
  geom_smooth(method='lm', se=FALSE)
#Negative correlation before accounting for species and Sex 


#------------- SVL across ELEVATION accounting for SPECIES and SEX--------------------# 

ggplot(sals, aes(ELEVATION_,SVL_))+
  geom_point(alpha=0.5)+
  ggtitle("SVL across Elevationals")+
  ylab("SVL")+ 
  xlab("Elevation")+
  facet_grid(Sex ~ Species)+
  geom_smooth(method='lm', se=FALSE)

library(car)
leveneTest(SVL_~Species,sals) #homogeneity of variance: checks that the variances in different 
#groups created by the categorical independent variable are equal 
#uneven variance in SVL_ across species 


