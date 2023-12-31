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
setwd("../")


#basic parameters
subsetsize<-1000
folderForFigures<-"./figures/defence_systems_vs_phylogenetic_distance/"
folderForData<-"./data/defence_systems_vs_phylogenetic_distance/"

if (!dir.exists(folderForFigures)){
  dir.create(folderForFigures)
}else{
  print("dir exists")
}

if (!dir.exists(folderForData)){
  dir.create(folderForData)
}else{
  print("dir exists")
}
####
#read Mash distance matrix
getCorrelationPlot<-function(mashmatrix,tree,defensesystemsfile,subset=0,prefix = "",rx=0.01,ry=0.97){
  # mashmatrix<-"./data/mash_distances/Ecoli_matrix.phylip"
  # tree<-"./data/Ecoli_phylogroups/Ecoli_E1.nwk"
  # defensesystemsfile<-"./data/26k_Ecoli_with_prophages.csv"
  # subset<-0
  # prefix<-"Ecoli_E1"
  # rx<-0.005
  
  DistanceMatrix<-read.csv(mashmatrix, sep=" ", skip=1, header=F,
                           stringsAsFactors = F)
  DistanceMatrix$V1<-str_replace(DistanceMatrix$V1,'.fna','')
  colnames(DistanceMatrix)<-c("Genome", DistanceMatrix$V1)
  rownames(DistanceMatrix)<-DistanceMatrix$Genome
  
  #read phylogenetic tree in order to remove contaminated genomes
  Tree<-read.tree(tree)
  newtiplabels<-gsub(".fna","",Tree$tip.label)

  ###subset genomes
  GenomesSubset<-newtiplabels
  if (subset > 0)
  {
    GenomesSubset<-sample(newtiplabels,subsetsize)
  }
  
  
  #subset mash distances
  DistmatSubsetPre<-subset(DistanceMatrix, DistanceMatrix$Genome %in%
                             GenomesSubset)
  DistmatSubset<-DistmatSubsetPre[,colnames(DistmatSubsetPre) %in% GenomesSubset]
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
  
  #####Save data in excel file
  ToSave<-ForPlotNoOutliers
  colnames(ToSave)[3:4]<-c("MashDistance","DefenseSystemsContentDistance")
  
  
  if(prefix =="")
  {
    preprefix<-strsplit(mashmatrix,"/")[[1]]
    prefix<-str_replace(preprefix[length(preprefix)],"_matrix.phylip","")
  }
 
  filename<-paste0(prefix,"_data_for_plot.xlsx")
  write_xlsx(ToSave,path =paste0(folderForData,"/",filename))
  #####
  
  
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
                alpha=0.5,
                se = TRUE,
                level=0.95,
                linewidth=0.5)+
    scale_fill_gradientn(colors = c("#f7fcf0","#e0f3db","#ccebc5",
                                 "#a8ddb5","#7bccc4","#4eb3d3",
                                 "#2b8cbe","#0868ac","#084081"),
                          name="Count (log10)")+
    ylim(c(-0.01,1.03))+
    annotate("text", label=paste0("r = ",format(round(CorCoef$estimate, 2), nsmall = 2),
                                 "\n",pvalue),
             x=rx,y=ry,
             hjust=0,
             size=1.5,
             family="ArialMT")+
    xlab("Mash distance")+
    ylab("Distance between defense system composition (Jaccard)")+
    theme_classic()+
    theme(#axis.text = element_text(size=5),
          text = element_text(family="ArialMT",size=5),
          #axis.title=element_text(size=5),
          legend.title = element_text(size=5),
          legend.key.size = unit(0.2,"cm"))
  CorrelationPlotDensity
  
  #######Save figures
  plotname<-paste0(folderForFigures,"/",prefix,"_sub",subset,
                   "_correlation_phylogenetic_distance_vs_defense_systems_densities")
  ggsave(paste0(plotname,".png"),
         plot=CorrelationPlotDensity,
         width=10, height=8, units="cm", dpi=300)
  ggsave(paste0(plotname,".svg"),
         plot=CorrelationPlotDensity,
         width=10, height=8, units="cm", dpi=300)
}
####################
#Doing it for all datasets
#Whole 26k Ecoli
getCorrelationPlot("./data/mash_distances/Ecoli_matrix.phylip",
                   "./data/Ecoli_tree_rapidnj.rM2.treeshrink_corrected.nwk",
                   "./data/26k_Ecoli_with_prophages.csv",subset=subsetsize)

#Enterobacteriales dataset
getCorrelationPlot("./data/mash_distances/entertree_matrix.phylip",
                   "./data/entertree_resavediTOL_newick.txt",
                   "./data/merged_enter3.csv",subset=subsetsize)
#Pseumonadales dataset
getCorrelationPlot("./data/mash_distances/pseutree_matrix.phylip",
                   "./data/pseutree_resavediTOL_newick.txt",
                   "./data/merged_pseu3.csv",subset=subsetsize)
#Bacillales dataset
getCorrelationPlot("./data/mash_distances/bacitree_matrix.phylip",
                   "./data/bacitree_resavediTOL_newick.txt",
                   "./data/merged_baci3.csv",subset=subsetsize)
#Burkholdriales dataset
getCorrelationPlot("./data/mash_distances/burktree_matrix.phylip",
                   "./data/burktree_resavediTOL_newick.txt",
                   "./data/merged_burk3.csv",subset=subsetsize)

########
#Do the same for the individual E.coli phylogroups

#the total number of leaves in those groups is  ~1000, that is why no subsetting here
getCorrelationPlot("./data/mash_distances/Ecoli_matrix.phylip",
                   "./data/Ecoli_phylogroups/Ecoli_E2.nwk",
                   "./data/26k_Ecoli_with_prophages.csv",prefix = "Ecoli_E2",rx=0.005)
getCorrelationPlot("./data/mash_distances/Ecoli_matrix.phylip",
                   "./data/Ecoli_phylogroups/Ecoli_E1.nwk",
                   "./data/26k_Ecoli_with_prophages.csv",prefix = "Ecoli_E1",rx=0.005)
getCorrelationPlot("./data/mash_distances/Ecoli_matrix.phylip",
                   "./data/Ecoli_phylogroups/Ecoli_C.nwk",
                   "./data/26k_Ecoli_with_prophages.csv",prefix = "Ecoli_C",rx=0.005)

#sub-setting for larger groups
getCorrelationPlot("./data/mash_distances/Ecoli_matrix.phylip",
                   "./data/Ecoli_phylogroups/Ecoli_A.nwk",
                   "./data/26k_Ecoli_with_prophages.csv",subset=subsetsize,prefix = "Ecoli_A",
                   rx=0.005)
getCorrelationPlot("./data/mash_distances/Ecoli_matrix.phylip",
                   "./data/Ecoli_phylogroups/Ecoli_B1.nwk",
                   "./data/26k_Ecoli_with_prophages.csv",subset=subsetsize,prefix = "Ecoli_B1",
                   rx=0.005)
getCorrelationPlot("./data/mash_distances/Ecoli_matrix.phylip",
                   "./data/Ecoli_phylogroups/Ecoli_B21.nwk",
                   "./data/26k_Ecoli_with_prophages.csv",subset=subsetsize,prefix = "Ecoli_B21",
                   rx=0.005)
getCorrelationPlot("./data/mash_distances/Ecoli_matrix.phylip",
                   "./data/Ecoli_phylogroups/Ecoli_B22.nwk",
                   "./data/26k_Ecoli_with_prophages.csv",subset=subsetsize,prefix = "Ecoli_B22",
                   rx=0.005)










