library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(ggpubr)
library(ggh4x)

path<-getwd()
setwd(path)
setwd("../")
###Read xls files for all five datasets

datasets<-c("pseu","baci","burk","enter")
postfix<-"_pagel_all_results_with_direction_0.01.xlsx"

files<-paste("./data/",datasets,postfix,sep="")
files<-c(files,"./data/Ecoli_pagel_all_results_with_direction.xlsx")

CorrelData=lapply(files,read_excel)
names(CorrelData)<-c(datasets,"ecoli")

CorelDataDf<-data_frame(id = names(CorrelData), CorrelData) %>%
  unnest(cols = c(CorrelData))
CorelDataDf$sysA<-ifelse(CorelDataDf$System.I == order(CorelDataDf$System.I,CorelDataDf$System.II)[1],
                                CorelDataDf$System.I, CorelDataDf$System.II)
CorelDataDf$sysB<-ifelse(CorelDataDf$System.II == order(CorelDataDf$System.I,CorelDataDf$System.II)[2],
                                CorelDataDf$System.II, CorelDataDf$System.I)

#####Get in wide format everything that has no pagel calculations
CorelDataWide<-spread(CorelDataDf[,c(1,5,17,18)],key="id",value="Pagel.p.value")
PairsNotPresent<-CorelDataWide %>% pivot_longer(cols = 'baci':'pseu') %>%
  filter(is.na(value))


CorelDataDfWithDir<-subset(CorelDataDf, CorelDataDf$Benjamini.Hochberg == "Y")

getsubset<-function(datasets,selectOpDir=F, selectSaDir=F){
  #datasets<-c("ecoli","enter")
  #selectOpDir<-T
  SubsetDataset<-subset(CorelDataDfWithDir,CorelDataDfWithDir$id %in% datasets)
  dupe = SubsetDataset[,c("sysA","sysB")] # select columns to check duplicates
  Duplicated<-SubsetDataset[duplicated(dupe) | duplicated(dupe, fromLast=TRUE),]
  
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
  #df<-Fivep1
  MissingValues<-subset(PairsNotPresent,
                        paste(PairsNotPresent$sysA,PairsNotPresent$sysB)%in% df$pair)
  MissingValues$direction<-rep("No pair",length(MissingValues$sysA))
  MissingValues$signif<-rep("",length(MissingValues$sysA))
  colnames(MissingValues)[3]<-c("id")
  ForHeatmapPlot<-rbind(MissingValues[,c(1:3,5,6)],
                    df[,c(1,16:19)])
  #In order to keep the order of the columns uniform
  ForHeatmapPlot$id<-factor(ForHeatmapPlot$id,
                               levels=c("ecoli",
                                        "enter",
                                        "baci",
                                        "burk",
                                        "pseu"))
  
  heatmap<-ggplot(data = ForHeatmapPlot, aes(id,interaction(sysB,sysA, sep="&")))+
    geom_tile(color = "#d9d9d9",
              aes(fill = as.factor(direction),
                  width=0.98, height=0.95))+
    # coord_equal(expand=T)+
    geom_text(aes(label=signif),
              size=2.1,
              color="white",
              vjust=0.8,
              hjust=0.5,
              inherit.aes = TRUE)+
    guides(
      y = guide_axis_nested(delim = "&"))+
    geom_vline(xintercept=1.5)+
    scale_fill_manual(values=c("-1" = "#35978f","1" = "#bf812d","No pair" = "#d9d9d9"), name="",
                      labels = c("-1" = "Mutually exclusive","1" = "Co-occuring",
                                 "No pair" = "No pair"),
                      na.value = "#d9d9d9")+
    scale_y_discrete(limits=rev, name="")+
    scale_x_discrete(labels = c("E. coli",
                                "Enterobacteriales",
                                "Bacillales",
                                "Burkholderiales",
                                "Pseudomonadales"), name ="")+
    guides(alpha="none")+
    theme_classic()+
    theme(axis.text = element_text(family="ArialMT", size=4),
          legend.text = element_text(family="ArialMT", size=4),
          legend.key.size = unit(0.15, 'cm'),
          axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,
                                     face="italic"),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          ggh4x.axis.nestline = element_line(linewidth = 0.15))
    
  heatmap
  return(heatmap)
}


################
#All datasets together
ForHeatmapFive<-getsubset(unique(CorelDataDfWithDir$id))

#write_xlsx(ForHeatmapFive, path="FiveDatasets_compare_directions.xlsx")

ForHeatmapFive$pair<-paste(ForHeatmapFive$sysA,ForHeatmapFive$sysB)

#do multiple heatmaps and merge
FivePairs<-sort(unique(ForHeatmapFive$pair))
parts<-4
parlen<-length(FivePairs)%/%parts
Fivep1<-subset(ForHeatmapFive,ForHeatmapFive$pair %in% FivePairs[1:parlen])
Fivep2<-subset(ForHeatmapFive,ForHeatmapFive$pair %in% FivePairs[(parlen+1):(parlen*2)])
Fivep3<-subset(ForHeatmapFive,ForHeatmapFive$pair %in% FivePairs[(parlen*2+1):(parlen*3)])
Fivep4<-subset(ForHeatmapFive,ForHeatmapFive$pair %in% FivePairs[(parlen*3+1):length(FivePairs)])

FiveHeat1<-draw_heatmap(Fivep1)
FiveHeat2<-draw_heatmap(Fivep2)
FiveHeat3<-draw_heatmap(Fivep3)
FiveHeat4<-draw_heatmap(Fivep4)

FiveDatasets<-ggarrange(FiveHeat1,
                        FiveHeat2,
                        FiveHeat3,
                        FiveHeat4,
                        ncol=4,
                        common.legend = T)
FiveDatasets

###Save
ggsave("./figures/Pagel_compare_all_datasets.png",plot=FiveDatasets,
       height = 12, width = 19, dpi =300, units = "cm")
ggsave("./figures/Pagel_compare_all_datasets.svg",plot=FiveDatasets,
       height = 12, width = 19, dpi =300, units = "cm")




########
#Verifying that directions are not randomly assigned
#for E.coli and Enter datasets

EcoliEnter<-getsubset(c("ecoli","enter"))
Ecoli<-subset(EcoliEnter[EcoliEnter$id=="ecoli",])
Enter<-subset(EcoliEnter[EcoliEnter$id=="enter",])
#counting number of real mismatches
MergedReal<-merge(Ecoli,Enter, by=c("sysA","sysB"))
RealCount<-xtabs(~MergedReal$direction.y+MergedReal$direction.x)
RealSum<-sum(RealCount[2,1]+RealCount[1,2])

#shuffling directions in datasets
NotMatching<-vector()

for(j in c(1:1000)){
  Ecoli<-subset(EcoliEnter[EcoliEnter$id=="ecoli",])
  Enter<-subset(EcoliEnter[EcoliEnter$id=="enter",])
  Ecoli$direction<-sample(Ecoli$direction)
  Enter$direction<-sample(Enter$direction)
  Merged<-merge(Ecoli,Enter, by=c("sysA","sysB"))
  Res<-xtabs(~Merged$direction.y+Merged$direction.x)
  Sum<-sum(Res[1,2]+Res[2,1])
  NotMatching<-c(NotMatching,Sum)
}

#plotting the outcome
ggplot()+
  geom_bar(aes(x=NotMatching),
           width=2)+
  geom_vline(xintercept=RealSum,
             color="red")+
  theme_classic()





