# #plan <- drake_plan(
#   #taxon_names = read.csv(file=file_in("data/taxa.csv")),
#   #tree.id = rotl::tnrs_match_names(taxon_names$Taxon[1])$ott_id,
#   #tree = rotl::tol_subtree(ott_id=tree.id),
#   #tree_print = plot_tree(tree, file=file_out("results/plethodontid_tree.pdf")),
#   ##phy = ape::rcoal(20+round(stats::runif(1,1,20))),
#   ##even = is_even(phy),
#   ##plotted = plot_tree(phy, file_out("results/tree.pdf")),
#   ##save_even = save(phy, even, file=file_out("Results/even.rda"))
# #)
# 
# 
# #plan <- drake_plan(
#   #phy = ape::read.tree("data/pyron_2011.tre"), # Tree from https://tree.opentreeoflife.org/curator/study/view/pg_423?, Pyron, R.A., & Wiens J.J. 2011. A large-scale phylogeny of Amphibia including over 2800 species, and a revised classification of extant frogs, salamanders, and caecilians. Molecular Phylogenetics Evolution 61 (2): 543-583. Accessed on Feb 26 2021
#   #discrete_ORGANdata = read.csv(file="data/salamanders.csv", stringsAsFactors=FALSE),
#   #pruned_objects = geiger::treedata(phy, data= discrete_ORGANdata, sort=TRUE, warning=FALSE),
#   #pruned_phy = pruned_objects$phy,
#   #pruned_data = pruned_objects$data, 
#   #tree_print = plot_tree(pruned_phy, file=file_out("Results/plethodontid_tree.pdf"))
#   #even = is_even(phy),
#   #plotted = plot_tree(phy, file_out("Results/plethodontidae_tree.pdf"))
#   #save_even = save(phy, even, file=file_out("results/even.rda"))
#   #)
# 
# 
# #plan <- drake_plan(
#   #taxon_names = read.csv(file=file_in("data/taxa.csv")),
#   #tree.id = rotl::tnrs_match_names(taxon_names$Taxon[1])$ott_id,
#   #phy = rotl::tol_subtree(ott_id=tree.id),
#   #discrete_ORGANdata = read.csv(file="data/salamanders.csv", stringsAsFactors=FALSE),
#   #pruned_objects = geiger::treedata(phy, data= discrete_ORGANdata, sort=TRUE, warning=FALSE),
#   #pruned_phy = pruned_objects$phy,
#   #pruned_data = pruned_objects$data, 
#   #tree_print = plot_tree(pruned_phy, file=file_out("Results/plethodontid_tree.pdf"))
#   #even = is_even(phy),
#   #plotted = plot_tree(phy, file_out("Results/plethodontidae_tree.pdf"))
#   #save_even = save(phy, even, file=file_out("results/even.rda"))
# #)
# 
# 
# 
# plan <- drake_plan(
#   #phy <- ape::read.tree("data/pyron_2011.tre") # Tree from https://tree.opentreeoflife.org/curator/study/view/pg_423?, Pyron, R.A., & Wiens J.J. 2011. A large-scale phylogeny of Amphibia including over 2800 species, and a revised classification of extant frogs, salamanders, and caecilians. Molecular Phylogenetics Evolution 61 (2): 543-583. Accessed on Feb 26 2021
#   phy = read.nexus("data/amphibian_tree.nex")
#   phy$tip.label = unname(GetAllGenusSpecies(phy$tip.label))
#   species_discretedata = read.csv("data/salamanders.csv")
#   pruned = geiger::treedata(phy = phy, data = species_discretedata, sort=TRUE)
#   pruned_phy = pruned$phy
#   pruned_data = pruned$discrete_ORGANdata
#   all_SVL = GetSizeForAllSpecies(phy)
#   #even = is_even(phy),
#   #plotted = plot_tree(phy, file_out("Results/plethodontidae_tree.pdf"))
#   #save_even = save(phy, even, file=file_out("results/even.rda"))
# )
# 
# 

plan <- drake_plan(
  PruneTree(tre,sal.data2),
  avg.size = male.female.avg.size(sal.data2),
  size.class = Size.Class.Female(avg.size),
)


# # plan <- drake_plan(
# #   tk = which(tre$tip.label%in%sal.data2$gen2=='TRUE'),
# #   tr2 = keep.tip(tre, tk), ## 14 tips
# #   plot(tr2) ## okay looks good.
# # )
# 
# CleanData(tre,sal.data2)
# 
# cleaned.data<-matrix(nrow=14, ncol=2)
# rownames(cleaned.data) <- unique(sal.datasal$gen2)
# colnames(cleaned.data) <- c('female', 'male')
# for (i in 1:length(unique(sal.datasal$gen2))){
#   spec<-unique(sal.datasal$gen2)[i]
#   specdat<-sal.datasal[which(sal.datasal$gen2==spec),]
#   SVLF<-specdat[which(specdat$Sex=='female'), ]
#   SVLF<-round(mean(SVLF$SVL))
#   SVLM<-specdat[which(specdat$Sex=='male'), ]
#   SVLM<-round(mean(SVLM$SVL))
#   cleaned.data[i,1]<-SVLF
#   cleaned.data[i,2]<-SVLM
# }
# 
# 
# cleaned.data.size.class<-matrix(nrow=14, ncol=1)
# rownames(cleaned.data.size.class) <- unique(disdatsal$gen2)
# colnames(cleaned.data.size.class) <- c('size class')
# for (i in 1:length(cleaned.data[,1])){
#   
#   cleaned.data.size.class<- ifelse (cleaned.data[,1] > 50, "Large","Small")
#   print(cleaned.data.size.class)
# }
# 
# all.discrete<- cbind(cleaned.data, cleaned.data.size.class)

