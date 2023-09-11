library(readxl)
library(dplyr)
library(tidyr)
# library(ggplot2)
# library(ggnetwork)
# library(ggpubr)
library(igraph)

path<-getwd()
setwd(path)
setwd("../")

######
#create folder to output figures

folderForResults<-"./figures/all_datasets_interaction_graph"

if (!dir.exists(folderForResults)){
  dir.create(folderForResults)
}else{
  print("dir exists")
}

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


######
#read raw data on system co-occurence
rawfiles<-paste0("./data/merged_",datasets,"3.csv")
rawfiles<-c(rawfiles,"./data/26k_Ecoli_with_prophages.csv")

RawData<-lapply(rawfiles,read.csv)
names(RawData)<-names(CorrelData)
RawDataDf<-data_frame(id = names(RawData), RawData) %>%
  unnest(cols = c(RawData))
RawDataDf$defense_system2<-ifelse(is.na(RawDataDf$defense_system2),
                                  RawDataDf$defense_system,
                                  RawDataDf$defense_system2)
RawDataDfRed<-RawDataDf[,c("id","start","end","genome","defense_system2")]

##generate counts for nodes
SystemCountsDf<-RawDataDf[,c("id","genome","defense_system2")] %>%
  group_by(id,defense_system2) %>% count()

######
#generate graph for each dataset
generate_graph<-function(dataset)
{
  # dataset<-"enter"
  # layparam<-0.0000001
  DatasetCorelDf<-subset(CorelDataDf,CorelDataDf$id == dataset &
                           CorelDataDf$Benjamini.Hochberg =="Y")
  DatasetCorelDf$signif<-ifelse(is.na(DatasetCorelDf$Bonferroni),"BH","B")
  count_coocurrences<-function(systs)
  {
    sysA<-systs[1]
    sysB<-systs[2]
    # sysA<-unlist(DatasetCorelDf[1,"sysA"])
    # sysB<-unlist(DatasetCorelDf[1,"sysB"])
    
    InitialCounts<-RawDataDf[RawDataDf$id==dataset & RawDataDf$defense_system2 %in% c(sysA,sysB),] %>%
      group_by(genome) %>% summarise(sysc=length(unique(defense_system2)))
    count<-nrow(InitialCounts[InitialCounts$sysc==2,])
    return(count)
  }
  DatasetCorelDf$CooccurCounts<-apply(DatasetCorelDf[,c("sysA","sysB")],MARGIN=1,count_coocurrences)
  
  ########
  edges<-DatasetCorelDf[,c("sysA","sysB","direction","signif","CooccurCounts")]
  edges$direction<-as.character(edges$direction)
  names(edges)<-c("from","to","direction","type","weight")
  nodes<-SystemCountsDf[SystemCountsDf$id ==dataset & 
                          SystemCountsDf$defense_system2 %in% c(DatasetCorelDf$sysA,
                                                                DatasetCorelDf$sysB),c(2:3)]
  names(nodes)<-c("id","count")
  net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
  
  #get direction color
  dircolors <- c("#35978f","#bf812d")
  names(dircolors)<-c("-1","1")
  E(net)$color <- dircolors[E(net)$direction]
  V(net)$size <-log(V(net)$count)/log(max(V(net)$count))*10+1
  E(net)$width<-E(net)$weight/max(E(net)$weight)*5+0.5
  E(net)$lty<-ifelse(E(net)$type == "BH",3,1)
  
  # lay<-layout_with_fr(net,
  #                     weights = log2(E(net)$width),
  #                     niter=1000)#layparam)
  #lay<-layout_with_graphopt(net, charge=layparam)
  png(file=paste0(folderForResults,"/",dataset,"_pagel_graph.png"),
      width=19, height=19, res=300,units = "cm")
  plot(net,edge.arrow.size=.4,
       vertex.label.color="gray20",
       vertex.label.family="ArielMT",
       vertex.frame.color = "#c7e9b4",
       vertex.color ="#ffffd9",
       vertex.label.cex=log(V(net)$count)/log(max(V(net)$count))/3,
       layout = layout_in_circle)
  legend(x="bottomleft", c("Mutually exclusive","Co-occuring"), pch=21,
         col="#777777", pt.bg=dircolors, pt.cex=.8, cex=.8, bty="n", ncol=1)
  dev.off()
  #######
  lay<-layout_with_fr(net,
                      weights = E(net)$width)
  png(file=paste0(folderForResults,"/",dataset,"_pagel_graph_fr.png"),
      width=19, height=19, res=300,units = "cm")
  plot(net,edge.arrow.size=.4,
       vertex.label.color="gray20",
       vertex.label.family="ArielMT",
       vertex.frame.color = "#c7e9b4",
       vertex.color ="#ffffd9",
       vertex.label.cex=log(V(net)$count)/log(max(V(net)$count))/3,
       layout = lay)
  legend(x="bottomleft", c("Mutually exclusive","Co-occuring"), pch=21,
         col="#777777", pt.bg=dircolors, pt.cex=.8, cex=.8, bty="n", ncol=1)
  dev.off()
}

generate_graph("pseu")
generate_graph("enter")
generate_graph("ecoli")
generate_graph("burk")
generate_graph("baci")

# 
# 
# PlNet<-ggnetwork(net, layout = in_circle())
# PlNet$linewidth<-PlNet$weight /max(PlNet$weight,na.rm=T)
# 
# ggplot(PlNet, aes(x=x,y=y,yend=yend,xend=xend)) +
#   geom_edges(aes(linetype = type,
#                  color = as.factor(direction),
#                  linewidth=linewidth))+
#   geom_nodes(aes(size=log(count)/log(max(count))),
#              color = "#c7e9b4",pch=21,
#              fill="#ffffd9")+
#   geom_nodetext(aes(label = name),
#                 family = "ArialMT")+
#   theme_blank()











