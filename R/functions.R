PruneTree <- function( phy, data){
  tk <- which(phy$tip.label%in%data[,22]=='TRUE')
  tr2 <- keep.tip(phy, tk) ## 14 tips
  plot(tr2) ## okay looks good.
  return(tr2)
}


male.female.avg.size <- function(data){
  cleaned.data<-matrix(nrow=14, ncol=2)
  rownames(cleaned.data) <- unique(data[,22])
  colnames(cleaned.data) <- c('female', 'male')
  for (i in 1:length(unique(data[,22]))){
    spec<-unique(data[,22])[i]
    specdat<-sal.data2[which(data[,22]==spec),]
    SVLF<-specdat[which(specdat$Sex=='female'), ]
    SVLF<-round(mean(SVLF$SVL))
    SVLM<-specdat[which(specdat$Sex=='male'), ]
    SVLM<-round(mean(SVLM$SVL))
    cleaned.data[i,1]<-SVLF
    cleaned.data[i,2]<-SVLM
  }
  return(cleaned.data)
}

Size.Class.Female <- function(data){
  # female.size.class<-matrix(nrow=14, ncol=1)
  # rownames(female.size.class) <- unique(data)
  # colnames(female.size.class) <- c('size class')
  for (i in 1:length(data[,1])){
   
    female.size.class<- ifelse (data[,1] > 50, "Large","Small")
  }
  return(female.size.class)
  # all.discrete <- cbind(data, female.size.class)
  # return(all.discrete)
}



VisualizeData <- function(phy, data) {
   tree_print <- plot(phy)
   return(data)
 }




)







# plot_tree <- function(phy, file) {
#   pdf(file=file)
#   plot(phy)
#   dev.off()
# }
# 
# 
# plot_tree <- function(phy, file) {
#   pdf(file=file, width=50, height=50)
#   plot(phy)
#   plot(phy, type="fan")
#   plot(phy, type="fan", show.tip.label=FALSE, edge.width=0.1)
#   dev.off()
# }
# 
# CleanData <- function(phy, data) {
#   pruned<- geiger::treedata(phy=phy, data=data, sort=TRUE, warning=FALSE)
#   pruned_phy <- pruned$phy
#   pruned_data <- pruned$data
# }
# 
# 
# 
# GetSizeForSpecies <- function(species){
#   print(species)
#   saldata <- read.csv("data/salamanders.csv")
#   size <- subset(saldata, Genus == c(species) & Sex == c("male","female"), select = c(SVL))
#   mean_size <- mean(size$SVL, na.rm=TRUE)
#   return(mean_size)
# }
# 
# 
# GetSizeForAllSpecies <- function(tree) {
#   result <- data.frame(species=tree$tip.label, SVL=NA)
#   for (tip_index in seq_along(tree$tip.label)) {
#     SVL <- GetSizeForSpecies(tree$tip.label[tip_index])
#     result$SVL[tip_index] <- SVL[1]
#   }
#   return(result)
# }
# 
# 
# 
# GetSingleGenusSpecies <- function(x) {
#   return(paste(strsplit(x, " |_")[[1]][1:2], collapse=" "))
# }
# GetAllGenusSpecies <- function(x) {
#   sapply(x, GetSingleGenusSpecies)
# }
# 
# 
# VisualizeData <- function(phy, data) {
#   tree_print <- plot(phy)
#   return(data)
# }
# 
# newcol<-matrix(nrow=15, ncol=1)
# rownames(newcol) <- unique(disdatsal$gen2)
# colnames(newcol) <- c('size_class')
# SizeClass <- function(x) {
#   for (i in 1:length(cleaned.data)){
#     spec<-cleaned.data[i]
#     newcol[i]<-ifelse(spec$SVLF>50,"big","small")
#     return(newcol)
#     }
#   }
# 
# 
# 
# GetLatDat <- function(x){
#   for (i in 1:length(unique(x))) {
#     
#     mycoinfo <- mycoportal(taxon = i)
#     latlong <- mycoinfo@records$Locality[1:1]
#     Species_lat <- strsplit(latlong, " ")[[1]][1]
#     
#   }
# }
