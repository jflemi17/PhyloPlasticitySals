#setwd("/Users/morganfleming/Desktop/Phylomethods/Phylo_Final_Project/data")
tre<-read.nexus("~/Desktop/Phylomethods/Phylo_Final_Project/data/amphibian_tree.nex")

sals.data<-read.csv("~/Desktop/Phylomethods/Phylo_Final_Project/data/salamanders.csv")

sals.data$Genus[which(sals.data$Genus=='Leurognathus')]<-'Desmognathus marmoratus' 
salamander <- sals.data
salamander$gen2<-gsub(' ', '_',salamander$Genus)
salamander$gen2[which(salamander$gen2=='Desmognathus_monticola_')] <- 'Desmognathus_monticola'
salamander$gen2[which(salamander$gen2=='Plethodon_cinereus_')] <- 'Plethodon_cinereus'
salamander$gen2[which(salamander$gen2=='Eurycea_bislineata_')] <- 'Eurycea_bislineata'
salamander$gen2[which(salamander$gen2=='Plethodon_glutinosis')] <- 'Plethodon_glutinosus'
salamander<-salamander[!(salamander$gen2=="Desmognathus_ochrophaeus"),]

unique(salamander$Genus)

salamander<- salamander[!is.na(salamander$Elevation),]

salamander <- salamander[ ,c("Genus", "Sex", "SVL", "Elevation")]

# #Adults
# salamander <-salamander[salamander$Sex == "female" | salamander$Sex == "male",]
# 
# str(salamander)
# salamander$Genus <- as.factor(salamander$Genus)
# salamander$Sex <- as.factor(salamander$Sex)
# str(salamander)
