library(phytools)
library(geiger)
library(stringr)
library(dplyr)
library(tidyr)
library(castor)
library(ggtree)
library(ggplot2)
library(ggpubr)


path<-getwd()
setwd(path)
setwd("../")

args = commandArgs(trailingOnly = T)
#args=c(6067,"./pagel_enter","entertree_resavediTOL_newick.txt","merged_enter3.csv",scalingFactor=0.1)
# second argument is folder
# test if there is at least one argument: if not, return an error
if (length(args) < 5) {
  stop("Five arguments must be supplied (task number; folder name; tree filename, defense info filename and scaling factor (0.1 is suggested value)).n", call.=FALSE)
} 

print(args)

####
#Read phylogenetic tree
tree <- read.tree(args[3])
###do proper names on the tree that match the data
newtipnames<-tree$tip.label
rep_str = c("'"='','.fna'='')
newtipnames<-str_replace_all(newtipnames,rep_str)
tree$tip.label<-newtipnames

#adjust edge lengths
newedgelengths<-ifelse(tree$edge.length>0,tree$edge.length, 0.0)
tree$edge.length<-newedgelengths

#####
#Read data
DefenseFile<-read.csv(args[4], header = T)

#column defense_system2 contains uniform system names that allows to work with padloc and defensefinder outputs
DefenceBySystem<-DefenseFile %>% group_by(genome,defense_system2) %>%
  count(defense_system2)
DefenceBySystemWide<-as.data.frame(DefenceBySystem %>% 
                                     pivot_wider(names_from = defense_system2, values_from = n))

#Filter data
DefenceBySystemWideFiltered<-subset(DefenceBySystemWide, DefenceBySystemWide$genome %in% newtipnames)

#subset tree to include only tips for which at least one defence system was found
tipnameswithdefense<-newtipnames[newtipnames %in% DefenceBySystemWideFiltered$genome]
treedef<-get_subtree_with_tips(tree, only_tips = tipnameswithdefense,
                               collapse_monofurcations = TRUE,
                               force_keep_root = F)$subtree

####add small value to change zero length branches
treedef$edge.length<-treedef$edge.length+0.0000001
####scale the tree branches as recommended in order to make fiting parameters stable
# default value is to scale to the mean of 1
branchlengthmeans<-mean(treedef$edge.length)
scalingfactor<-as.numeric(args[5])/branchlengthmeans
treedefscaled<-treedef
#also get rid of zero length branches
treedefscaled$edge.length<-treedefscaled$edge.length*scalingfactor

#check mean and median
mean(treedefscaled$edge.length)
median(treedefscaled$edge.length)

#change only to presence/absence without counts
DefenceBySystemWideBinary<-DefenceBySystemWideFiltered[,c(2:ncol(DefenceBySystemWideFiltered))]
DefenceBySystemWideBinary[is.na(DefenceBySystemWideBinary)]<-0
DefenceBySystemWideBinary[DefenceBySystemWideBinary>0]<-1

rownames(DefenceBySystemWideBinary)<-DefenceBySystemWideFiltered$genome

#this is important step that was missing before 20230131
DefenceBySystemWideBinary <- DefenceBySystemWideBinary[match(treedef$tip.label,rownames(DefenceBySystemWideBinary)),]

name.check(treedef, DefenceBySystemWideBinary)


# #Pagel correction
# ##################
rows<-ncol(DefenceBySystemWideBinary)

#number to task at sge farm
tasks<-c(1:((rows*rows-rows)/2))
ivec<-vector()
jvec<-vector()
for(i in colnames(DefenceBySystemWideBinary)[1:rows]) {
  for(j in colnames(DefenceBySystemWideBinary)[1:rows]) {
    if (i>j) {
      ivec<-append(ivec,i)
      jvec<-append(jvec,j)
    }
  }
}
tasksdf<-data.frame(tasks=tasks,i=ivec,j=jvec)

#write.csv(tasksdf, file="tasksdf.csv")
##################
taskToPerform <-tasksdf[tasksdf$tasks==args[1],]

i <- taskToPerform$i
print(i)
j <- taskToPerform$j
print(j)
var1 <- DefenceBySystemWideBinary %>% pull(i)
names(var1) <- rownames(DefenceBySystemWideBinary)
var2 <- DefenceBySystemWideBinary %>% pull(j)
names(var2) <- rownames(DefenceBySystemWideBinary)

res <- fitPagel(treedefscaled, x = var1, y = var2, method="fitDiscrete")

ResToWrite<-data.frame(v1=i,v2=j,v3=res$P[1])



#write model for prominent cases
if (res$P[1] < 0.01) {
  imod<-str_replace(i,'/','_')
  jmod<-str_replace(j,'/','_')
  dirfordetails<-paste(args[2],"/",imod,"_",jmod,"/",sep="")
  
  dir.create(dirfordetails)
  
  ###
  resi <- fitPagel(treedefscaled, x = var1, y = var2, dep.var = "x", method="fitDiscrete")
  resj <- fitPagel(treedefscaled, x = var1, y = var2, dep.var = "y", method="fitDiscrete")
  ResToWrite$v4<-resi$P[1]
  ResToWrite$v5<-resj$P[1]

  png(file=paste(dirfordetails,"/",imod,"_",jmod,"_pagel_probabilities.png", sep=""),
      width=2600,height=2600, units ="px",res=300)

  plot(res,lwd.by.rate=TRUE,
       main=c("Independent model",
              "Full dependent model"))
  dev.off()
  
  modeloutput<-paste(dirfordetails,"/",imod,"_",jmod,"_model.output",sep="")
  sink(modeloutput)
  print(res)
  sink()
  
  if (resi$P[1]<0.01){
    write.csv(resi$dependent.Q, file = paste(dirfordetails,"/",imod,"_",jmod,"_",imod,"_dependent.matrix",sep=""),
              quote = F)
    png(file=paste(dirfordetails,"/",imod,"_",jmod,"_",imod,"_pagel_probabilities.png", sep=""),
        width=2600,height=2600, units ="px",res=300)
    
    plot(resi,lwd.by.rate=TRUE,
         main=c("Independent model",
                "Full dependent model"))
    dev.off()
  }
  if (resj$P[1]<0.01){
    write.csv(resj$dependent.Q, file = paste(dirfordetails,"/",imod,"_",jmod,"_",jmod,"_dependent.matrix",sep=""),
              quote = F)
    png(file=paste(dirfordetails,"/",imod,"_",jmod,"_",jmod,"_pagel_probabilities.png", sep=""),
        width=2600,height=2600, units ="px",res=300)
    
    plot(resj,lwd.by.rate=TRUE,
         main=c("Independent model",
                "Full dependent model"))
    dev.off()
  }

  write.csv(res$dependent.Q, file = paste(dirfordetails,"/",imod,"_",jmod,"_dependent.matrix",sep=""),
            quote = F)

  #draw plots

  DFForVizualization<-DefenceBySystemWideBinary
  DFForVizualization$genome<-rownames(DefenceBySystemWideBinary)
  
  vizualizetree<-function(dsys1,dsys2)
  {
    #dsys1<-"Zorya II"
    #dsys2<-"ietAS"
    DsysDf<-as.data.frame(subset(DFForVizualization, DFForVizualization[[dsys1]]== 1 |
                                   DFForVizualization[[dsys2]]== 1)[,c("genome",dsys1,dsys2)])
    #subset tree
    subtree<-get_subtree_with_tips(treedef, only_tips = DsysDf$genome)$subtree
    #orderDf
    DsysDfforplot<-subset(DsysDf, DsysDf$genome %in% subtree$tip.label)
    #Trees are drawn with ggtree
    Sys1Plot<-ggtree(subtree, size =0.1, color = "#525252") %<+% DsysDfforplot +
      geom_tippoint(aes(color = as.factor(get(dsys1))),
                    size =0.7,
                    show.legend = F) +
      scale_color_manual(values = c("black","#31a354"), breaks=c(0,1))+
      ggtitle(dsys1)+
      theme(plot.title = element_text(hjust = 0.5),
            plot.margin = unit(c(0.5,0,0.5,0),"mm"))
    #
    Sys2Plot<-ggtree(subtree, size =0.1, color = "#525252") %<+% DsysDfforplot +
      geom_tippoint(aes(color = as.factor(get(dsys2))),
                    size =0.7,
                    show.legend = F) +
      scale_color_manual(values = c("black","#fd8d3c"), breaks=c(0,1)) +
      scale_x_reverse() +
      ggtitle(dsys2)+
      theme(plot.title = element_text(hjust = 0.5),
            plot.margin = unit(c(0.5,0,0.5,0),"mm"))
    #combine
    PlotToSave<-ggarrange(Sys1Plot, Sys2Plot,
                          align = "h",
                          ncol = 2)
    PlotToSave
    
    ggsave(paste(str_replace(dsys1,"/","_"),"_",str_replace(dsys2,"/","_"),"_tree.png", sep=""), path=dirfordetails,
           PlotToSave,
           width=3600,height=2400, units ="px",dpi=300)
    
  }
  
  ##
  vizualizetreefull<-function(dsys1,dsys2)
  {
    #dsys1<-"Gabija"
    #dsys2<-"tmn"
    DsysDf<-as.data.frame(DFForVizualization[,c("genome",dsys1,dsys2)])
    
    Sys1Plot<-ggtree(treedef, size =0.1, color = "#525252") %<+% DsysDf +
      geom_tippoint(aes(color = as.factor(get(dsys1))),
                    size =0.7,
                    show.legend = F) +
      scale_color_manual(values = c("black","#31a354"), breaks=c(0,1))+
      xlab(dsys1)+
      layout_dendrogram()
    #
    Sys2Plot<-ggtree(treedef, size =0.1, color = "#525252") %<+% DsysDf +
      geom_tippoint(aes(color = as.factor(get(dsys2))),
                    size =0.7,
                    show.legend = F) +
      scale_color_manual(values = c("black","#fd8d3c"), breaks=c(0,1)) +
      coord_flip() +
      xlab(dsys2)
    #combine
    PlotToSavefull<-ggarrange(Sys1Plot, Sys2Plot,
                              align = "h",
                              ncol = 1)
    PlotToSavefull
    ggsave(paste(str_replace(dsys1,"/","_"),"_",str_replace(dsys2,"/","_"),"_treefull.png", sep=""), path=dirfordetails,
           PlotToSavefull,
           width=6000,height=2400, units ="px",dpi=300)
  }
  vizualizetree(i,j)
  vizualizetreefull(i,j)
}

###Write the main results down
filename<-paste(args[2],"/",args[1],".txt",sep="")

write.table(ResToWrite, file=filename, sep="\t",row.names=F,
            quote=F,col.names = F)


