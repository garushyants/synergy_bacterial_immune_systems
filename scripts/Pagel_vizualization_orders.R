library(phytools)
library(stringr)
library(dplyr)
library(tidyr)
library(castor)
library(ggplot2)
library(ggpubr)
library(writexl)


path<-getwd()
setwd(path)
setwd("../")

datasets<-c("pseu","baci","burk","enter")
folderWithRawData<-"../20230314_pagel_other_datasets/"
for (dt in datasets) {
  #dt<-datasets[1]
  folderwithmatrices<-paste0(folderWithRawData,"pagel_",dt)
  defensedata<-paste(folderWithRawData,"merged_",dt,"3.csv",sep="")
  treefile<-paste("./data/",dt,"tree_resavediTOL_newick.txt",sep="")
  #merge p-values together and write to the local folder
  system(paste("cat ",folderwithmatrices,"/*.txt > ","./data/",dt,"pagel_fitDiscrete_all.tsv", sep=""))
  
  #######################
  #Open file with pagel p-values
  pvaluecutoff<-0.01
  
  #Read in the data
  PagelResults<-read.csv(paste("./data/",dt,"pagel_fitDiscrete_all.tsv",sep=""), header =F, sep="\t", fill =T,
                         col.names = paste0("V",seq_len(5)))
  
  
  
  
  #####read counts data in order to check results stability
  DefenseDt<-read.csv(defensedata, header = T)
  
  DefenceBySystem<-DefenseDt %>% group_by(genome,defense_system2) %>%
    count(defense_system2)
  DefenceBySystemWide<-as.data.frame(DefenceBySystem %>% 
                                       pivot_wider(names_from = genome, values_from = n))
  #filter
  tree <- read.tree(treefile)
  ###do proper names on the tree that match the data
  newtipnames<-tree$tip.label
  rep_str = c("'"='','.fna'='')
  newtipnames<-str_replace_all(newtipnames,rep_str)
  
  DefenceBySystemWideFiltered<-DefenceBySystemWide[, colnames(DefenceBySystemWide) %in% newtipnames]
  DefenceBySystemWideBinary<-DefenceBySystemWideFiltered
  DefenceBySystemWideBinary[is.na(DefenceBySystemWideBinary)]<-0
  DefenceBySystemWideBinary[DefenceBySystemWideBinary>0]<-1
  
  DefenceBySystemWideBinary$sum<-rowSums(DefenceBySystemWideBinary)
  
  DefenceBySystemWideBinary$System<-DefenceBySystemWide$defense_system2
  
  DefenseCounts<-DefenceBySystemWideBinary[,c("System","sum")]
  
  #Get a list of very rare systems to remove
  mincount<-0.01 #Filtering most rare and most common systems
  SystemsToKeep<-DefenseCounts[(DefenseCounts$sum > length(newtipnames)*mincount) &
                                 (DefenseCounts$sum < length(newtipnames)*(1-mincount)),]$System
  
  # filter<-sqrt(length(newtipnames))
  # SystemsToKeep<-DefenseCounts[(DefenseCounts$sum > filter) &
  #                                (DefenseCounts$sum < (length(newtipnames)-filter)),]$System
  
  ###
  #Get vector of systems ordered by abundance
  #It is use later on to sort the data on heatmap
  SystemOrderDf<-subset(DefenseCounts, DefenseCounts$System %in% SystemsToKeep)
  SystemOrder<-SystemOrderDf[order(-SystemOrderDf$sum),]$System
  
  ######################################
  #Filter and do correction for multiple testing
  PagelResultsFiltered<-subset(PagelResults, (PagelResults$V1 %in% SystemsToKeep) & (PagelResults$V2 %in% SystemsToKeep))
  PagelResultsSorted<-PagelResultsFiltered[order(PagelResultsFiltered$V3),]
  PagelResultsSorted$rank<-c(1:length(PagelResultsSorted$V1))
  PagelResultsSorted$Bonferroni<-ifelse(PagelResultsSorted$V3<pvaluecutoff/length(PagelResultsSorted$V1), "Y","")
  PagelResultsSorted$Benjamini.Hochberg<-ifelse(PagelResultsSorted$V3<pvaluecutoff*PagelResultsSorted$rank/length(PagelResultsSorted$V1), "Y","")
  colnames(PagelResultsSorted)[1:5]<-c("System.I","System.II","Pagel.p.value","P.value.System.I.dep", "P.value.System.II.dep")
  
  
  ########
  #Merge with significant values
  PreCount<-merge(PagelResultsSorted,DefenseCounts, by.x = "System.I", by.y = "System")
  colnames(PreCount)[9]<-"Sum.System.I"
  PagelResultsSortedCounts<-merge(PreCount, DefenseCounts, by.x ="System.II", by.y ="System")
  colnames(PagelResultsSortedCounts)[10]<-"Sum.System.II"
  #
  
  #select only once that are somehow significant
  PagelResultsSignifCounts<-subset(PagelResultsSortedCounts,PagelResultsSortedCounts$Pagel.p.value<pvaluecutoff)
  
  DirectionDf<-data.frame(rank=NULL,direction=NULL)
  #add transitions
  for (j in PagelResultsSignifCounts$rank)
  {
    #j<-1
    task<-PagelResultsSignifCounts[PagelResultsSignifCounts$rank==j,]
    #these two strings below are added to deal with AbiO/Nhi
    task$System.I<-ifelse(grepl("/",task$System.I), str_replace(task$System.I,"/","_"),task$System.I)
    task$System.II<-ifelse(grepl("/",task$System.II), str_replace(task$System.II,"/","_"),task$System.II)
    #most probable model
    bestmodel<-which.min(task[3:5])
    prefix<-ifelse(bestmodel == 1, paste(task$System.I,"_",task$System.II, sep=""),
                      ifelse(bestmodel == 2, paste(task$System.I,"_",task$System.II,"_",task$System.I, sep=""),
                             paste(task$System.I,"_",task$System.II,"_",task$System.II, sep="")))
    #prefix<-paste(task$System.I,"_",task$System.II, sep="")
    depmatrixfile<-paste(folderwithmatrices,"/",paste(task$System.I,"_",task$System.II, sep=""),"/",
                      prefix,"_dependent.matrix", sep="")
    
    if (file.exists(depmatrixfile))
    {
      DependentMatrixFile<-read.csv(depmatrixfile)
      DependentMatrix<-(DependentMatrixFile[,c(2:5)]) %>% mutate_if(is.numeric, round, digits=4)
      
      #calculate normalized ratios
      # ratios<-c((DependentMatrix[2,4]-DependentMatrix[1,3])/DependentMatrix[2,4],
      #           (DependentMatrix[4,2]-DependentMatrix[3,1])/DependentMatrix[4,2],
      #           (DependentMatrix[3,4]-DependentMatrix[1,2])/DependentMatrix[3,4],
      #           (DependentMatrix[4,3]-DependentMatrix[2,1])/DependentMatrix[4,3])
      
      #DependentMatrix<-DependentMatrix+0.0001
      
      ratio<-function(v1,v2){
        smallvalue<-0.00000001
        v1*(v2+smallvalue)/(v2*v2+smallvalue)
      }
      ratios<-c(ratio(DependentMatrix[2,4],DependentMatrix[4,2]),
      ratio(DependentMatrix[1,3],DependentMatrix[3,1]),
      ratio(DependentMatrix[3,4],DependentMatrix[4,3]),
      ratio(DependentMatrix[1,2],DependentMatrix[2,1]))
      
      direction<-ifelse((ratios[2] + ratios[4]) > (ratios[1]+ratios[3]),-1,1)
      
      # #I do that to remove very small values that sometimes occur
      # directionnum<-ifelse(ratios[1]==0 | ratios[2]==0, ratios[3]-ratios[4],
      #                   ifelse(ratios[3]==0 | ratios[4]==0, ratios[1]-ratios[2],
      #                          ratios[1]-ratios[2]+ ratios[3]-ratios[4]))
      # direction<-ifelse(directionnum >0 ,1,-1)
      
      taskoutdf<-data.frame(rank=j,ratios[1], ratios[2], ratios[3], ratios[4],direction)
      DirectionDf<-rbind(DirectionDf,taskoutdf)
  
  
    }
    else
    {
      print(j)
    }
  }
  
  #write results to Excel table
  ExcelDf<-merge(PagelResultsSortedCounts,DirectionDf, by="rank",all.x=T)
  
  write_xlsx(ExcelDf,paste("./data/",dt,"_pagel_all_results_with_direction_",pvaluecutoff,".xlsx",sep=""))
  
  
  ####
  PagelResultsDirectionF<-as.data.frame(merge(PagelResultsSignifCounts,DirectionDf, by="rank"))
  
  #do symmetric
  PagelResultsSymHalf1<-PagelResultsDirectionF
  PagelResultsSymHalf2<-PagelResultsDirectionF
  colnames(PagelResultsSymHalf2)[2:3]<-c("System.I","System.II")
  PagelResultsSym<-rbind(PagelResultsSymHalf1,PagelResultsSymHalf2)
  
  #do full
  PagelWide<-as.data.frame(pivot_wider(PagelResultsSym[,c(2,3,15)],
                                       names_from = System.II, 
                                       values_from = direction))
  
  PagelWideSortedRow<-PagelWide[order(match(PagelWide$System.I,SystemOrder)),]
  PrePagelWideSorted<-PagelWideSortedRow[,c(2:ncol(PagelWideSortedRow))]
  PagelWideSorted<-PagelWideSortedRow[,c(2:ncol(PagelWideSortedRow))]%>%
    select(order(match(colnames(PagelWideSortedRow[,c(2:(ncol(PagelWideSortedRow)))]),SystemOrder)))
  rownames(PagelWideSorted)<-PagelWideSortedRow$System.I
  
  #do triangles
  Positive<-PagelWideSorted %>% mutate(across(everything(), function(x){replace(x, which(x==-1), NA)}))
  Negative<-PagelWideSorted %>% mutate(across(everything(), function(x){replace(x, which(x==1), NA)}))
  
  Posmat<-as.matrix(Positive)
  Posmat<-replace(Posmat, is.na(Posmat),0)
  Negmat<-as.matrix(Negative)
  Negmat<-replace(Negmat, is.na(Negmat),0)
  
  CombinedMat <- matrix(NA, nrow = nrow(PagelWideSorted), ncol = ncol(PagelWideSorted))
  CombinedMat[upper.tri(CombinedMat)] <- Posmat[upper.tri(Posmat)]
  CombinedMat[lower.tri(CombinedMat)] <- Negmat[lower.tri(Negmat)]
  
  ##back to df
  PagelWideTri<-as.data.frame(CombinedMat)
  colnames(PagelWideTri)<-colnames(PagelWideSorted)
  PagelWideTri$System.I<-rownames(PagelWideSorted)
  
  #back to long
  PreResultsToPlot<-PagelWideTri %>% 
    pivot_longer(!System.I, names_to="System.II", values_to = "direction")
  ResultsToPlot<-merge(PreResultsToPlot,PagelResultsSym, by=c("System.I","System.II","direction"),
                       all.x=T)
  ResultsToPlot$Signif<-ifelse(ResultsToPlot$Bonferroni == "Y","**",
                               ifelse(ResultsToPlot$Benjamini.Hochberg == "Y","*",""))
  
  #Order labels by abundance
  ResultsToPlot$System.I<-factor(ResultsToPlot$System.I,levels=SystemOrder)
  ResultsToPlot$System.II<-factor(ResultsToPlot$System.II,levels=SystemOrder)
  
  #Plot results
  heatmap<-ggplot(data = ResultsToPlot, aes(System.I,System.II))+
    geom_tile(color = "#d9d9d9",
              aes(fill = direction, 
              alpha = 1-Pagel.p.value))+
    geom_text(aes(label=Signif),
              size=4,
              color="white",
              vjust=0.8,
              hjust=0.5,
              inherit.aes = TRUE)+
    guides(alpha="none")+
    theme_classic()+
    scale_fill_gradient2(high= "#bf812d", mid="white", low="#35978f", name ="",
                         labels= c("Mutually exclusive","","","","Co-occuring"))+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=12),
          axis.text.y = element_text(size=12),
          axis.title = element_blank(),
          legend.key.size = unit(0.5, 'cm'))
  heatmap
  
  ggsave(paste("./figures/",dt,"_heatmap_pagel_all_models_001.png",sep=""), plot= heatmap,
         height = 31, width =35, units ="cm", dpi=200)
  ggsave(paste("./figures/",dt,"_heatmap_pagel_all_models_001.svg",sep=""), plot= heatmap,
         height = 31, width =35, units ="cm", dpi=200)
}
