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