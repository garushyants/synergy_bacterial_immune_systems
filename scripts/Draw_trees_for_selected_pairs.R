  library(phytools)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(castor)
  library(ggtree)
  library(ggplot2)
  library(ggpubr)
  
  
  path<-getwd()
  setwd(path)
  
  ####
  #Read phylogenetic tree
  tree <- read.tree("../data/Ecoli_tree_rapidnj.rM2.treeshrink_corrected.nwk")
  ###do proper names on the tree that match the data
  newtipnames<-tree$tip.label
  rep_str = c("'"='','.fna'='')
  #rename tips
  newtipnames<-str_replace_all(newtipnames,rep_str)
  tree$tip.label<-newtipnames
  
  #####
  #Read occurence data
  EcoliDefense<-read.csv("../data/26k_Ecoli_with_prophages.csv", header = T)
  
  DefenseBySystem<-EcoliDefense %>% group_by(genome,defense_system) %>%
    count(defense_system)
  DefenseBySystemWide<-as.data.frame(DefenseBySystem %>% 
                                       pivot_wider(names_from = defense_system, values_from = n))
  #Filter data to the set only present on the tree
  DefenseBySystemWideFiltered<-subset(DefenseBySystemWide, DefenseBySystemWide$genome %in% newtipnames)
  
  #subset tree to include only tips for which at least one defence system was found
  tipnameswithdefense<-newtipnames[newtipnames %in% DefenseBySystemWideFiltered$genome]
  treedef<-get_subtree_with_tips(tree, only_tips = tipnameswithdefense)$subtree
  
  #change only to presence/absence without counts
  #this table is required for the phylogenetic tree presence/absence
  DefenseBySystemWideBinary<-DefenseBySystemWideFiltered[,c(2:ncol(DefenseBySystemWideFiltered))]
  DefenseBySystemWideBinary[is.na(DefenseBySystemWideBinary)]<-0
  DefenseBySystemWideBinary[DefenseBySystemWideBinary>0]<-1
  
  DefenseBySystemWideBinary$genome<-DefenseBySystemWideFiltered$genome
  
  #sort
  DefenseBySystemWideBinary <- DefenseBySystemWideBinary[match(treedef$tip.label,DefenseBySystemWideBinary$genome),]
  
  ###Draw phylogenetic trees
  locationpallete<-c("#1b9e77","#d95f02","#7570b3","#999999")
  locationpallete_n<-setNames(object=locationpallete,
                              nm<-c("chromosome","plasmid","prophage","unknown"))
  
  vizualizetree<-function(dsys1,dsys2)
  {
    #dsys1<-"Zorya II"
    #dsys2<-"Druantia III"
    DsysDf<-as.data.frame(subset(DefenseBySystemWideBinary, DefenseBySystemWideBinary[[dsys1]]== 1 |
                                   DefenseBySystemWideBinary[[dsys2]]== 1)[,c("genome",dsys1,dsys2)])
    #subset tree
    subtree<-get_subtree_with_tips(treedef, only_tips = DsysDf$genome)$subtree
    #orderDf
    DsysDfforplot<-subset(DsysDf, DsysDf$genome %in% subtree$tip.label)
    
    ##get info about location in the genome
    get_location<-function(sys)
    {
      #sys<-dsys2
      Dsys1Location<-subset(EcoliDefense,EcoliDefense$defense_system %in% sys & 
                              EcoliDefense$genome %in% DsysDfforplot$genome)
      Dsys1Location$Location<-ifelse(Dsys1Location$prophage_within == "full",
                                     "prophage",
                                     Dsys1Location$seqid_type)
      Dsys1LocationLong<-Dsys1Location[,c("genome","Location","defense_system")]
      Dsys1LocationWide<-dcast(Dsys1LocationLong,genome~Location)
      values<-colnames(Dsys1LocationWide)[2:ncol(Dsys1LocationWide)]
      for (v in values)
      {
        Dsys1LocationWide[[v]]<-ifelse(Dsys1LocationWide[[v]]== 0,NA,v)
      }
      Dsys1LocationForPlot<-as.data.frame(Dsys1LocationWide[,c(2:ncol(Dsys1LocationWide))])
      rownames(Dsys1LocationForPlot)<-Dsys1LocationWide$genome
      return(Dsys1LocationForPlot)
    }
    
    halfplot<-function(sys,direction=1,color)
    {
      Sys1Plot<-ggtree(subtree, size =0.1, color = "#525252")%<+% DsysDfforplot +
        geom_tippoint(aes(color = as.factor(get(sys))),
                      size =0.7,
                      show.legend = F) +
        scale_color_manual(values = c("black",color), breaks=c(0,1))+
        ggtitle(sys) +
        theme(plot.title = element_text(hjust = 0.5),
              plot.margin = unit(c(0.5,1,0.5,0),"mm")) +
              geom_treescale(fontsize=3, x=direction*0.01,y=-120,
                       offset=20)
        
      if(direction == -1)
      {
        Sys1Plot<-Sys1Plot + scale_x_reverse()
          
      }
      loc1df<-get_location(sys)
      Sys1PlotHeat<-gheatmap(Sys1Plot,loc1df, width=0.02*ncol(loc1df),offset=.0005,
                             color=NULL,
                             colnames_position="bottom", 
                             colnames_angle=90,
                             hjust=1.1,
                             font.size = 3)+
        scale_fill_manual(name="",values=locationpallete_n,na.value="white")+
        theme(legend.position = "none")+
        ggtree::vexpand(.1, -1)+
        ggtree::vexpand(.02, 1)
      return(Sys1PlotHeat)
    }
    Sys1PlotHeat<-halfplot(dsys1,1,"#ff7f00")
    Sys2PlotHeat<-halfplot(dsys2,-1,"#377eb8")
    
    #combine
    PlotToSave<-ggarrange(Sys1PlotHeat, Sys2PlotHeat,
                          align = "h",
                          ncol = 2)
    PlotToSave
    
    ggsave(paste("../figures/Ecoli_phylogenetic_trees/",dsys1,"_",dsys2,"_tree.png", sep=""),
           PlotToSave,
           width=3600,height=2700, units ="px",dpi=300)
    
    ggsave(paste("../figures/Ecoli_phylogenetic_trees/",dsys1,"_",dsys2,"_tree.svg", sep=""),
           PlotToSave,
           width=3600,height=2700, units ="px",dpi=300)
    
  }
  
  vizualizetree("Zorya II","ietAS")
  vizualizetree("Zorya II","Druantia III")
  vizualizetree("tmn","Gabija")
  vizualizetree("Kiwa","ietAS")
  vizualizetree("BREX I","RM IV")
  
  
  
