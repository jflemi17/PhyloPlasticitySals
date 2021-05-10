#Prune Tree
PruneTree <- function( phy, data){
  tk <- which(phy$tip.label%in%data[,1]=='TRUE')
  tr2 <- keep.tip(phy, tk) ## 14 tips
  plot(tr2)
  return(tr2)
}

#Visualize Data 
VisualizeData <- function(phy, data) {
   tree_print <- plot(phy)
   return(data)
}

#Plot Tree
plot_tree <- function(phy, file) {
  pdf(file=file)
  plot(phy)
  dev.off()
}


#Subset samples by three elevational intervals -- High, Mid, Low. Elevations were selected so as to evenly distribute samples across the groups 
elev.class <- function(data){
  data$elev<- NA
  elevation <- data$Elevation
  for(i in 1:length(data)){
    if(elevation[i] <= 3500){data$elev[i] <- "Low"} else{
      if(elevation[i] > 3500 & data$Elevation[i] <= 4500){ data$elev[i] <-"Mid"} else{
        if(elevation[i] > 4500){data$elev[i] <- "High"} 
      }
    }  
  }
  return(data$elev)
  # sals.and.elev.class<- cbind(data$elev, data)
  # return(sals.and.elev.class)
}



#function for avg body size and variance to use for different subsets of species samples
Mean.SD.Species <- function(data){
  pheno.plast<- c()
  avg.svl<- c()
  all.species<- unique(data[,1])
  for(i in 1:length(unique(data$Genus))){
    species<- unique(data[,1])[i]
    spec <- data[which(data[,1]==species),]
    avg.SVL <-round(mean(spec$SVL))
    plasticity <- round(sd(spec$SVL),3)
    avg.svl[i]<-avg.SVL
    pheno.plast[i]<-plasticity 
  }
  all.dat <- cbind(all.species,pheno.plast, avg.svl)
  return(all.dat)
}
  

#Function for pulling random samples 
Get.subset.sal <- function (data){
  if (sals.dat[which(sals.dat$Sex=='female')]){
    female <- sample(1:nrow(sals.dat.size.female$SVL), 5, replace = FALSE)
  }
  sal.sub <- subset(sals.dat,Sex=='female', select = SVL)
  summary(sals.sub[,1:7]) 
  sub <- sample(1:nrow(sals.dat.size.female$SVL), 5, replace = FALSE)
  sub





#Graveyard 
#   
#   calc.data<-matrix(nrow=14, ncol=4)
#   rownames(calc.data) <- unique(data[,6])
#   colnames(calc.data) <- c('avg.size', 'variance')
#   for (i in 1:length(unique(data[,6]))){
#     species <-unique(data[,6])[i]
#     data.per.species <- data[which(data[,6]==species),]
#     female.only<-data.per.species[which(data.per.species$Sex=='female'), ]
#     female.avg <-round(mean(female.only$SVL))
#     female.size.variance <- sd(female.only$SVL)
#     male.only<-specdat[which(specdat$Sex=='male'), ]
#     male.avg <-round(mean(male.only$SVL))
#     male.size.variance <- sd(female.only$SVL)
#     calc.data[i,1]<-female.avg
#     calc.data[i,2]<-female.size.variance 
#     calc.data[i,3]<-male.avg 
#     calc.data[i,4]<-male.size.variance 
#   }
#   return(calc.data)
# }
# body.size.species.sub(sal.dat)