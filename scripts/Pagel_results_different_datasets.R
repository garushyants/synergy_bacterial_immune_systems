library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)

path<-getwd()
setwd(path)
setwd("../")
###Read xls files for all five datasets

datasets<-c("pseu","baci","burk","enter")
postfix<-"_pagel_all_results_20230329_all_models_001.xlsx"

files<-paste(datasets,postfix,sep="")
files<-c(files,"Ecoli_pagel_rescaled_all_results_20230329.xlsx")

CorrelData=lapply(files,read_excel)
names(CorrelData)<-c(datasets,"ecoli")

CorelDataDf<-data_frame(id = names(CorrelData), CorrelData) %>%
  unnest(cols = c(CorrelData))


CorelDataDfWithDir<-subset(CorelDataDf, CorelDataDf$Benjamini.Hochberg == "Y")
##
#Doing this procedure below to be sure that pairs are in the same order
CorelDataDfWithDir$sysA<-ifelse(CorelDataDfWithDir$System.I == order(CorelDataDfWithDir$System.I,CorelDataDfWithDir$System.II)[1],
                               CorelDataDfWithDir$System.I, CorelDataDfWithDir$System.II)
CorelDataDfWithDir$sysB<-ifelse(CorelDataDfWithDir$System.II == order(CorelDataDfWithDir$System.I,CorelDataDfWithDir$System.II)[2],
                               CorelDataDfWithDir$System.II, CorelDataDfWithDir$System.I)

getsubset<-function(datasets,selectOpDir=F, selectSaDir=F){
  #datasets<-c("ecoli","enter")
  #selectOpDir<-T
  SubsetDataset<-subset(CorelDataDfWithDir,CorelDataDfWithDir$id %in% datasets)
  dupe = SubsetDataset[,c("sysA","sysB")] # select columns to check duplicates
  Duplicated<-SubsetDataset[duplicated(dupe) | duplicated(dupe, fromLast=TRUE),]
  
  #Draw the heatmaps
  ForHeatmap<-Duplicated
  ForHeatmap[is.na(ForHeatmap)]<-""
  ForHeatmap$signif<-ifelse(ForHeatmap$Bonferroni=="Y","**",ifelse(ForHeatmap$Benjamini.Hochberg == "Y","*",""))
  if (selectOpDir)
  {
    Df<-ForHeatmap
    ForHeatmap<-Df %>% group_by(sysA,sysB)%>% filter(length(unique(direction))>1)
  }
  if (selectSaDir)
  {
    Df<-ForHeatmap
    ForHeatmap<-Df %>% group_by(sysA,sysB)%>% filter(length(unique(direction))==1)
  }
  return (ForHeatmap)
}


draw_heatmap<-function(df){
  heatmap<-ggplot(data = df, aes(id,paste(sysA,sysB)))+
    geom_tile(color = "#d9d9d9",
              aes(fill = as.factor(direction)))+
    geom_text(aes(label=signif),
              size=3.5,
              color="white",
              vjust=0.8,
              hjust=0.5,
              inherit.aes = TRUE)+
    scale_fill_manual(values=c("#80b1d3","#b3de69"))+
    scale_y_discrete(limits=rev, name="")+
    guides(alpha="none")+
    theme_classic()
  return(heatmap)
}

################
#Ecoli and Enter
ForHeatmapEcoliEnter<-getsubset(c("ecoli","enter"))
EcoliEnter<-draw_heatmap(ForHeatmapEcoliEnter)

ggsave("Ecoli_vs_Enter_20230329.png",plot=EcoliEnter,
       height = 28, width = 13, dpi =300, units = "cm")
#only opposite directions
ForHeatmapEcoliEnterOpDir<-getsubset(c("ecoli","enter"),selectOpDir = T)
EcoliEnterOpDir<-draw_heatmap(ForHeatmapEcoliEnterOpDir)
EcoliEnterOpDir
ggsave("Ecoli_vs_Enter_OpDir_20230329.png",plot=EcoliEnterOpDir,
       height = 18, width = 13, dpi =300, units = "cm")
#only same direction
ForHeatmapEcoliEnterSaDir<-getsubset(c("ecoli","enter"),selectSaDir = T)
EcoliEnterSaDir<-draw_heatmap(ForHeatmapEcoliEnterSaDir)
EcoliEnterSaDir
ggsave("Ecoli_vs_Enter_SameDir_20230329.png",plot=EcoliEnterSaDir,
       height = 25, width = 17, dpi =300, units = "cm")
################
#Four datasets
ForHeatmapFour<-getsubset(datasets)
FourDatasets<-draw_heatmap(ForHeatmapFour)

ggsave("Fourdatasets_20230329.png",plot=FourDatasets,
       height = 18, width = 13, dpi =300, units = "cm")
#only opposite directions
ForHeatmapFourOpDir<-getsubset(datasets,selectOpDir = T)
FourDatasetsOpDir<-draw_heatmap(ForHeatmapFourOpDir)
FourDatasetsOpDir
ggsave("Fourdatasets_OpDir_20230329.png",plot=FourDatasetsOpDir,
       height = 15, width = 13, dpi =300, units = "cm")
#same dir
ForHeatmapFourSaDir<-getsubset(datasets,selectSaDir = T)
FourDatasetsSaDir<-draw_heatmap(ForHeatmapFourSaDir)
FourDatasetsSaDir
ggsave("Fourdatasets_SameDir_20230329.png",plot=FourDatasetsSaDir,
       height = 15, width = 13, dpi =300, units = "cm")
################
#All datasets together
ForHeatmapFive<-getsubset(unique(CorelDataDfWithDir$id))

write_xlsx(ForHeatmapFive, path="FiveDatasets_compare_directions.xlsx")

FiveDatasets<-draw_heatmap(ForHeatmapFive)
FiveDatasets
ggsave("Fivedatasets_20230329.png",plot=FiveDatasets,
       height = 45, width = 25, dpi =300, units = "cm")
ggsave("Fivedatasets_20230616.svg",plot=FiveDatasets,
       height = 60, width = 20, dpi =300, units = "cm")
#only opposite directions
ForHeatmapFiveOpDir<-getsubset(unique(CorelDataDfWithDir$id),selectOpDir = T)
FiveDatasetsOpDir<-draw_heatmap(ForHeatmapFiveOpDir)
FiveDatasetsOpDir
ggsave("Fivedatasets_OpDir_20230329.png",plot=FiveDatasetsOpDir,
       height = 25, width = 20, dpi =300, units = "cm")

#same direction
ForHeatmapFiveSaDir<-getsubset(unique(CorelDataDfWithDir$id),selectSaDir = T)
FiveDatasetsSaDir<-draw_heatmap(ForHeatmapFiveSaDir)
FiveDatasetsSaDir
ggsave("Fivedatasets_SameDir_20230329.png",plot=FiveDatasetsSaDir,
       height = 25, width = 20, dpi =300, units = "cm")



# DupAgregated<-Duplicated %>% group_by(sysA,sysB) %>% summarise(across(c(id,direction),toString))
# 
# DupAggregatedDiffDirections<-subset(DupAgregated, grepl(" 1, -1",DupAgregated$direction) |
#          grepl("-1, 1",DupAgregated$direction) |
#            grepl("^1, -1",DupAgregated$direction))
# write.csv(DupAggregatedDiffDirections, "Four_datasets_all_models_001_opposite_directions_20230322.csv", row.names = F)
# 
# DupAggregatedSameDirection<-DupAgregated %>% anti_join(DupAggregatedDiffDirections)
# write.csv(DupAggregatedSameDirection, "Four_datasets_all_models_001_same_direction_20230322.csv", row.names = F)
# 
# ###
# #Add the tests outcome to figure out 
# PreDuplicatedStrict<-subset(CorelDataDfSignif,CorelDataDfSignif$Bonferroni == "Y")
# dupestrict<-PreDuplicatedStrict[,c("sysA","sysB")]
# DuplicatedStrict<-PreDuplicatedStrict[duplicated(dupestrict) | duplicated(dupestrict, fromLast=TRUE),]
# 
# DupAgregatedStrict<-DuplicatedStrict %>% group_by(sysA,sysB) %>% summarise(across(c(id,direction),toString))
# 
# DupAggregatedDiffDirections<-subset(DupAgregated, grepl(" 1, -1",DupAgregated$direction) |
#                                       grepl("-1, 1",DupAgregated$direction) |
#                                       grepl("^1, -1",DupAgregated$direction))



