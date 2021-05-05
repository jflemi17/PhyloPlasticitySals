setwd("/Users/morganfleming/Desktop/Phylomethods/PhyloPlasticitySals/data")
tre<-read.nexus("amphibian_tree.nex")
sal.data<-read.csv("salamanders.csv")
dim(sal.data)
unique(sal.data$Genus)

sal.data$Genus[which(sal.data$Genus=='Leurognathus')]<-'Desmognathus marmoratus' 
sal.data2 <- sal.data
sal.data2$gen2<-gsub(' ', '_',sal.data2$Genus)
sal.data2$gen2[which(sal.data2$gen2=='Desmognathus_monticola_')] <- 'Desmognathus_monticola'
sal.data2$gen2[which(sal.data2$gen2=='Plethodon_cinereus_')] <- 'Plethodon_cinereus'
sal.data2$gen2[which(sal.data2$gen2=='Eurycea_bislineata_')] <- 'Eurycea_bislineata'
sal.data2$gen2[which(sal.data2$gen2=='Plethodon_glutinosis')] <- 'Plethodon_glutinosus'
sal.data2<-sal.data2[!(sal.data2$gen2=="Desmognathus_ochrophaeus"),]

unique(sal.data2$Genus)
