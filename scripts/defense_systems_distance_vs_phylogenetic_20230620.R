library(stringr)
library(dplyr)
library(tidyr)
library(reshape2)
library(phytools)
library(ggplot2)
library(ggpubr)
library(writexl)


path<-getwd()
setwd(path)
setwd("../20230522_defense_distance_vs_phylogenetic/")


#basic parameters
subsetsize<-1000
folderForResults<-"results20230803"

if (!dir.exists(folderForResults)){
  dir.create(folderForResults)
}else{
  print("dir exists")
}
####
#read Mash distance matrix
getCorrelationPlot<-function(mashmatrix,tree,defensesystemsfile,subset=0,rx=0.01,ry=0.97){
  # mashmatrix<-"burktree_matrix.phylip"
  # tree<-"burktree_resavediTOL_newick.txt"
  # defensesystemsfile<-"merged_burk3.csv"
  # subset<-subsetsize
  
  DistanceMatrix<-read.csv(mashmatrix, sep=" ", skip=1, header=F,
                           stringsAsFactors = F)
  DistanceMatrix$V1<-str_replace(DistanceMatrix$V1,'.fna','')
  colnames(DistanceMatrix)<-c("Genome", DistanceMatrix$V1)
  rownames(DistanceMatrix)<-DistanceMatrix$Genome
  # ###The code below is to check that indeed I removed the most weird cases from the tree
  # DistmatRownamesColnames<-DistanceMatrix[,c(2:ncol(DistanceMatrix))]
  # 
  # SumDistancesDf<-data.frame(sumdist=colSums(DistmatRownamesColnames), 
  #                       genome = colnames(DistmatRownamesColnames))
  #read phylogenetic tree
  Tree<-read.tree(tree)
  newtiplabels<-gsub(".fna","",Tree$tip.label)
  
  # ##I do that to check that I remove the most weird distances if I use tree tips
  # SumDistancesDfFiltered<-subset(SumDistancesDf,
  #                                SumDistancesDf$genome %in% newtiplabels)
  # ##It is indeed the case, I don't have to think about it again
  
  ###subset genomes
  GenomesSubset<-newtiplabels
  if (subset > 0)
  {
    GenomesSubset<-sample(newtiplabels,subsetsize)
  }
  
  
  #subset mash distances
  DistmatSubsetPre<-subset(DistanceMatrix, DistanceMatrix$Genome %in%
                             GenomesSubset)
  DistmatSubset<-subset(DistmatSubsetPre, select=GenomesSubset)
  DistmatSubsetOrdered<-DistmatSubset[, sort(colnames(DistmatSubset))]
  
  #Convert subset to long
  PhylDistDf<-melt(as.matrix(DistmatSubset),varnames = c("row", "col"))
  # PhylDistWithoutRepeats<-PhylDistDf[as.numeric(PhylDistDf$row) > as.numeric(PhylDistDf$col), ]
  
  #############
  #read defense systems info
  EcoliDefense<-read.csv(defensesystemsfile, header = T)
  
  #rename column for other dataset
  if ("defense_system2" %in% colnames(EcoliDefense))
  {
    colnames(EcoliDefense)[colnames(EcoliDefense) == 'defense_system'] <- 'incorrect_names'
    colnames(EcoliDefense)[colnames(EcoliDefense) == 'defense_system2'] <- 'defense_system'
  }
  
  #get subset
  EcoliDefenseSubset<-subset(EcoliDefense,
                             EcoliDefense$genome %in% GenomesSubset)
  
  #transform
  DefenseBySystem<-EcoliDefenseSubset %>% group_by(genome,defense_system) %>%
    count(defense_system)
  
  DefenseBySystemWide<-as.data.frame(DefenseBySystem %>% 
                                       pivot_wider(names_from = defense_system, values_from = n))
  #binary
  DefenseBySystemWideBinary<-DefenseBySystemWide[,c(2:ncol(DefenseBySystemWide))]
  rownames(DefenseBySystemWideBinary)<-DefenseBySystemWide$genome
  DefenseBySystemWideBinary[is.na(DefenseBySystemWideBinary)]<-0
  DefenseBySystemWideBinary[DefenseBySystemWideBinary>0]<-1
  #filter out systems that are rarely found
  #I set 5% filter
  Systemsum<-colSums(DefenseBySystemWideBinary)
  mincount<-0.005 #System was found in 0.5% of all genomes
  SystemsToKeep<-names(Systemsum[Systemsum > nrow(DefenseBySystemWideBinary)*mincount])
  
  DefenseBySystemWideBinaryFiltered<-DefenseBySystemWideBinary[,colnames(DefenseBySystemWideBinary) %in%
                                                                 SystemsToKeep]
  
  #########
  #Calculate distances
  GenDist<-dist(data.matrix(DefenseBySystemWideBinaryFiltered,rownames.force = T),
                method ="binary")
  #melt distances
  GenDistDf<-melt(as.matrix(GenDist), varnames = c("row", "col"))
  
  ##########Plot
  ForPlotWithRepeats<-merge(PhylDistDf,GenDistDf, by=c('row','col'))
  
  #is it correct???
  ForPlot<-ForPlotWithRepeats[as.numeric(ForPlotWithRepeats$row) > as.numeric(ForPlotWithRepeats$col), ]
  ##remove outliers
  ForPlotNoOutliers<-subset(ForPlot,ForPlot$value.x < 0.95)
  
  ToSave<-ForPlotNoOutliers
  colnames(ToSave)[3:4]<-c("MashDistance","DefenseSystemsContentDistance")
  
  write_xlsx(ToSave,path =paste(folderForResults,"/",tree,"_data_for_plot.xlsx",sep=""))
  
  ###
  CorCoef<-cor.test(ForPlotNoOutliers$value.x, ForPlotNoOutliers$value.y, 
                    method="spearman",
                    exact = F)
  pvalue<-ifelse(CorCoef$p.value == 0, "p-value < 2.2e-16",
                 paste("p-value =",format(CorCoef$p.value, 2, digits=3)))
  
  #2d density plot
  CorrelationPlotDensity<-ggplot(ForPlotNoOutliers,
                          aes(x=value.x,
                              y=value.y))+
    geom_hex(aes(fill=log10(after_stat(count))),
             bins=80)+
    geom_smooth(method = "lm",
                formula = y~x,
                color="#993404",
                fill="#fec44f",
                se = T,
                level=0.99)+
    scale_fill_gradientn(colors = c("#f7fcf0","#e0f3db","#ccebc5",
                                 "#a8ddb5","#7bccc4","#4eb3d3",
                                 "#2b8cbe","#0868ac","#084081"),
                          name="Count(log10)")+
    ylim(c(-0.01,1.03))+
    annotate("text", label=paste("r = ",format(round(CorCoef$estimate, 2), nsmall = 2),
                                 "\n",pvalue,
                                 sep=""),
             x=rx,y=ry,
             hjust=0,
             size=3.5)+
    xlab("Mash distance")+
    ylab("Distance between defense system composition (Jaccard)")+
    theme_classic()+
    theme(axis.text = element_text(size=10))
  # #plot
  # CorrelationPlot<-ggplot(ForPlotNoOutliers,
  #        aes(x=value.x,
  #            y=value.y))+
  #   geom_point(alpha=.2,
  #              color="#02818a",
  #              size=0.5)+
  #   geom_smooth(method = "lm",
  #               formula = y~x,
  #               color="#993404",
  #               fill="#fec44f")+
  #   ylim(c(-0.01,1.03))+
  #   annotate("text", label=paste("r = ",format(round(CorCoef$estimate, 2), nsmall = 2),
  #                                "\n",pvalue,
  #                                sep=""),
  #             x=rx,y=ry,
  #           hjust=0,
  #           size=3.5)+
  #   xlab("Mash distance")+
  #   ylab("Distance between defense system composition (Jaccard)")+
  #   theme_classic()+
  #   theme(axis.text = element_text(size=10))
  return(CorrelationPlotDensity)
}
####################
#Doing it for all datasets
#26k Ecoli
EcoliPlot<-getCorrelationPlot("matrix.phylip","Ecoli_tree_rapidnj.rM2.treeshrink_corrected.nwk",
                              "complete_prophage_platon_26k_all7.csv",subset=1000)
ggsave("Ecoli26k_sub1000_correlation_phylogenetic_distance_vs_defense_systems_densities.png",
       plot=EcoliPlot, path=folderForResults,
       width=20, height=20, units="cm", dpi=300)
ggsave("Ecoli26k_sub1000_correlation_phylogenetic_distance_vs_defense_systems_densities.svg",
       plot=EcoliPlot, path=folderForResults,
       width=20, height=20, units="cm", dpi=300)
#Enetrobacteriales dataset
EnterPlot<-getCorrelationPlot("entertree_matrix.phylip","entertree_resavediTOL_newick.txt",
                              "merged_enter3.csv",subset=1000)
ggsave("Enter_sub1000_correlation_phylogenetic_distance_vs_defense_systems_densities.png",
       plot=EnterPlot, path=folderForResults,
       width=20, height=20, units="cm", dpi=300)
ggsave("Enter_sub1000_correlation_phylogenetic_distance_vs_defense_systems_densities.svg",
       plot=EnterPlot, path=folderForResults,
       width=20, height=20, units="cm", dpi=300)
#Pseumonadales dataset
PseuPlot<-getCorrelationPlot("pseutree_matrix.phylip","pseutree_resavediTOL_newick.txt",
                              "merged_pseu3.csv",subset=1000)
ggsave("Pseu_sub1000_correlation_phylogenetic_distance_vs_defense_systems_densities.png",
       plot=PseuPlot, path=folderForResults,
       width=20, height=20, units="cm", dpi=300)
ggsave("Pseu_sub1000_correlation_phylogenetic_distance_vs_defense_systems_densities.svg",
       plot=PseuPlot, path=folderForResults,
       width=20, height=20, units="cm", dpi=300)
#Bacillales dataset
BaciPlot<-getCorrelationPlot("bacitree_matrix.phylip","bacitree_resavediTOL_newick.txt",
                             "merged_baci3.csv",subset=1000)
ggsave("Baci_sub1000_correlation_phylogenetic_distance_vs_defense_systems_densities.png",
       plot=BaciPlot, path=folderForResults,
       width=20, height=20, units="cm", dpi=300)
ggsave("Baci_sub1000_correlation_phylogenetic_distance_vs_defense_systems_densities.svg",
       plot=BaciPlot, path=folderForResults,
       width=20, height=20, units="cm", dpi=300)

#Burkholdriales dataset
BurkPlot<-getCorrelationPlot("burktree_matrix.phylip","burktree_resavediTOL_newick.txt",
                             "merged_burk3.csv",subset=1000)
ggsave("Burk_sub1000_correlation_phylogenetic_distance_vs_defense_systems_densities.png",
       plot=BurkPlot, path=folderForResults,
       width=20, height=20, units="cm", dpi=300)
ggsave("Burk_sub1000_correlation_phylogenetic_distance_vs_defense_systems_densities.svg",
       plot=BurkPlot, path=folderForResults,
       width=20, height=20, units="cm", dpi=300)









