library(phytools)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ape)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(castor)

path<-getwd()
setwd(path)


#Read phylogenetic tree
tree <- read.tree("../data/Ecoli_tree_rapidnj.rM2.treeshrink_corrected.nwk")
###do proper names on the tree that match the data
newtipnames<-tree$tip.label
rep_str = c("'"='','.fna'='')
newtipnames<-str_replace_all(newtipnames,rep_str)
tree$tip.label<-newtipnames

#Read assembly_summary data
AssemblySummary<-read.csv("../data/Ecoli_numContigs_vs_numDefence/assembly_summary_refseq_ecoli.csv", header = T)
AssemblySummaryFiltered<-subset(AssemblySummary, AssemblySummary$assembly_accession %in% newtipnames)
###Read info about phylogroups
TenThouEcoliData<-read.csv("../data/Ecoli_phylogroups/Abram_metadata.csv", header=T)

PhyloGenomesInDataset<-subset(TenThouEcoliData,
                              TenThouEcoliData$id %in% AssemblySummaryFiltered$paired_asm_comp)[,c("id","Phylogroup")]

PreAlltipsCorrectNames<-merge(PhyloGenomesInDataset, AssemblySummaryFiltered,
                           by.x="id", by.y="paired_asm_comp", all.y = T)[,c(2,3)]
#add all other names
NTdf<-as.data.frame(newtipnames)
AlltipsCorrectNames<-merge(PreAlltipsCorrectNames,NTdf, by.x="assembly_accession", 
                           by.y ="newtipnames", all.y=T)

AlltipsCorrectNames$alpha<-ifelse(is.na(AlltipsCorrectNames$Phylogroup),0,1)

AlltipsCorrectNames<- AlltipsCorrectNames %>% replace(is.na(.),"none")

#########################
#get labels of interest
#read leaves names to add
# LeavesOfInterest<-read.csv("Ecoli_interesting strain.csv",
#                            stringsAsFactors = F)
# LeavesOfInterestShort<-subset(LeavesOfInterest[,c(3,5)],
#                               LeavesOfInterest$selected_strains !="")


#Load TreeCluster clusters for preliminary branch grouping
Clusters003<-read.csv("../data/Ecoli_phylogroups/TreeCluster_0.03_clusters", header=T, sep="\t")
Clusters003$SequenceName<-str_replace(Clusters003$SequenceName,".fna","")
ClustersOnly<-subset(Clusters003,Clusters003$ClusterNumber > 0)

ClusterNodes<-data.frame(node=character(), cluster=character())
clustersun<-unique(ClustersOnly$ClusterNumber)
clusnum<-length(clustersun)
getPalette=colorRampPalette(brewer.pal(9, "Set1"))
cluscol<-getPalette(clusnum)

for (i in c(1:clusnum))
{
  #cl = unique(Clusters003$ClusterNumber)[5]
  leaves=ClustersOnly[ClustersOnly$ClusterNumber==clustersun[i],1]
  
  node<-findMRCA(tree,leaves,type="node")
  cln<-data.frame(node=node,cluster=clustersun[i])
  ClusterNodes<-rbind(ClusterNodes,cln)
}


#########################################
#########################################
####Get trees for individual phylogroups from TreeCluster data

plot_clade<-function(clade)
{
  #vizualize to check that things were selected correctly
  cladeplot<-ggtree(clade, layout = 'circular') %<+% AlltipsCorrectNames +
    geom_tippoint(aes(color = Phylogroup, alpha =alpha),
                  size =0.8, show.legend = T) +
    scale_alpha(guide = 'none')+
    theme(legend.position = "right",
          plot.margin = unit(c(0,0,0,0), "mm"))
  return(cladeplot)
}
#############################
#Separate phylogroups E2 and E1
#E1 and E2 are cluster 11
Ephylo<-extract.clade(tree,ClusterNodes[ClusterNodes$cluster==11,"node"])
#plot_clade(Ephylo)

#E1
E2leaves<-AlltipsCorrectNames[AlltipsCorrectNames$Phylogroup=='E2',1]
E1refined_clade<-drop.clade(Ephylo,E2leaves)
#plot_clade(E1refined_clade)
#E2
E2clade<-extract.clade(Ephylo,findMRCA(Ephylo,E2leaves,type="node"))
#plot_clade(E2clade)

##############################
#Separate clades for phylogroups B1 and C
#B1 and C cluster 17
B1Cphylo<-extract.clade(tree,ClusterNodes[ClusterNodes$cluster==17,"node"])
#plot_clade(B1Cphylo)

#B1
Cleaves<-AlltipsCorrectNames[AlltipsCorrectNames$Phylogroup=='C',1]
B1refined_clade<-drop.clade(B1Cphylo,Cleaves)
#plot_clade(B1refined_clade)
#C
Cclade<-extract.clade(B1Cphylo,findMRCA(B1Cphylo,Cleaves,type="node"))
#plot_clade(Cclade)

# #Other important clades
# B21phylo<-extract.clade(tree,ClusterNodes[ClusterNodes$cluster==10,"node"])
# #plot_clade(B21phylo)
# 
# B22phylo<-extract.clade(tree,ClusterNodes[ClusterNodes$cluster==9,"node"])
# #plot_clade(B22phylo)
# 
# Aphylo<-extract.clade(tree,ClusterNodes[ClusterNodes$cluster==16,"node"])
# #plot_clade(Aphylo)

####
#Get updated node ids for phylogroups of interest
Ccommonnode<-findMRCA(tree,Cleaves,type="node")
B1commonnode<-findMRCA(tree,subset(tree$tip.label,tree$tip.label %in% B1refined_clade$tip.label), type="node")
E2commonnode<-findMRCA(tree,E2leaves,type="node")
E1commonnode<-findMRCA(tree,subset(tree$tip.label,tree$tip.label %in% E1refined_clade$tip.label), type="node")

ClustersRefined<-data.frame(node=c(Ccommonnode,B1commonnode,
                                   E2commonnode,E1commonnode,
                                   ClusterNodes[ClusterNodes$cluster==10,"node"],
                                   ClusterNodes[ClusterNodes$cluster==9,"node"],
                                   ClusterNodes[ClusterNodes$cluster==16,"node"]),
                            cluster = c("C","B1","E2","E1","B2-1","B2-2","A"))
#####get sets of tips for major phylogroups
MajorPhyloTips<-data.frame(phylo=c(rep("A",length(Aphylo$tip.label)),
                                   rep("B1",length(B1refined_clade$tip.label)),
                                   rep("B21",length(B21phylo$tip.label)),
                                   rep("B22",length(B22phylo$tip.label)),
                                   rep("E1",length(E1refined_clade$tip.label)),
                                   rep("E2",length(E2clade$tip.label)),
                                   rep("C",length(Cclade$tip.label))),
                           genome=c(Aphylo$tip.label,
                                    B1refined_clade$tip.label,
                                    B21phylo$tip.label,
                                    B22phylo$tip.label,
                                    E1refined_clade$tip.label,
                                    E2clade$tip.label,
                                    Cclade$tip.label))

###############################################
#Vizualize the whole E.coli phylogenetic tree

pupd<-ggtree(tree, layout = 'circular', open.angle=10, size=0.15) +
  geom_hilight(data = ClustersRefined, mapping = aes(node=node,fill = as.factor(cluster)),
               alpha=0.3)+
  #geom_tiplab(align=T)+
  geom_treescale(y=300, x=0.02, fontsize=3, linesize=0.7, offset=1.5)+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a",
                             "#33a02c","#fb9a99","#cab2d6",
                             "#6a3d9a"),
                    guide="none") 


#combine with phylogroups
p2upd<- pupd %<+% AlltipsCorrectNames + 
  geom_tippoint(aes(color = Phylogroup, alpha =alpha),
                size =0.8, show.legend = T) +
  scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a",
                              "#33a02c","#fb9a99","#e31a1c",
                              "#fdbf6f","#ff7f00","#cab2d6",
                              "#6a3d9a","#ffff99","#b15928",
                              "#000000")) +
  scale_alpha(guide = 'none')+
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), "mm"))+
  guides(color = guide_legend(override.aes = list(size=8)))

# pTLab<-p2upd %<+% LeavesOfInterestShort +
#   geom_tiplab2(aes(label=selected_strains))
# pTLab

ggsave("../figures/phylogenetic_trees/Ecoli_tree_with_phylogroups.png",
       plot=p2upd, width = 35, height=35, units = "cm", dpi=300)


###################################
###################################
#Draw tree subsample with
#information about defense systems
###################################

Ecoli26kDefensedata<-read.csv("../data/26k_Ecoli_with_prophages.csv", header = T)
Ecoli26kDefensedataFiltered<-subset(Ecoli26kDefensedata,
                                    Ecoli26kDefensedata$genome %in% newtipnames)
Ecoli26kDefensedataFiltered$location<-ifelse(Ecoli26kDefensedataFiltered$prophage_within=="full",
                                             paste("prophage",Ecoli26kDefensedataFiltered$seqid_type), 
                                             Ecoli26kDefensedataFiltered$seqid_type)
EcoliDefenseFiltRestructured<-Ecoli26kDefensedataFiltered %>% group_by(genome,location)%>%
  summarise(SysCount=n())
table(Ecoli26kDefensedataFiltered$location)

######################
###subsampling the whole tree
DrawSubsample<-function(genomes)
{
  EcoliDefenseSubsample<-subset(Ecoli26kDefensedataFiltered, Ecoli26kDefensedataFiltered$genome %in% genomes)
  EcoliDefenseSubsampleRestr<-EcoliDefenseSubsample %>% group_by(genome,location)%>% 
    summarise(syscount=n())
  EcolidefenseSubPhylo<-subset(MajorPhyloTips,
                               MajorPhyloTips$genome %in% EcoliDefenseSubsampleRestr$genome)[,c(2,1)]
  subtree<-get_subtree_with_tips(tree, only_tips = genomes)$subtree
  #Plot
  TreeSubsample<-ggtree(subtree, layout = 'circular', open.angle=5,
                        size=0.15) %<+% EcolidefenseSubPhylo + 
    geom_tippoint(aes(color = phylo),
                  size =1, show.legend = T) +
    scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a",
                                "#33a02c","#fb9a99","#cab2d6",
                                "#6a3d9a","#ffff99","#b15928",
                                "#000000"))+
    geom_treescale(y=2, x=0.02, fontsize=2.5, linesize=0.7, offset=1.3)+
    guides(color = guide_legend(override.aes = list(size=8)))
  TreeSubsample
  
  ##
  TreeSubsampleDefSys<-TreeSubsample +
    geom_fruit(data=EcoliDefenseSubsampleRestr,
               geom=geom_bar,
               mapping=aes(y=genome,
                           x=syscount,
                           fill=location),
               orientation="y",
               stat="identity",
               axis.params=list(
                 axis       = "x",
                 text.size  = 1.8,
                 hjust      = 1,
                 vjust      = 0.5,
                 nbreak     = 3,
               ),
               offset = 0.04,
               pwidth=0.5) +
    scale_fill_manual(values=c("#742c24",
                               "#ebddd3",
                               "#ba4535",
                               "#88727b",
                               "#feb24c",
                               "#ffeda0")) +
    theme(legend.position = "right",
          plot.margin = unit(c(0,0,0,0), "mm"))
  return(TreeSubsampleDefSys)
}
#Do 750 random genomes
GenomeSubsample<-sample(unique(EcoliDefenseFiltRestructured$genome),750)
Treesub750<-DrawSubsample(GenomeSubsample)
Treesub750
ggsave("../figures/phylogenetic_trees/Ecoli_tree_with_phylogroups_and_DefSystems_20230614_subsample750.png",
       plot=Treesub750, width = 35, height=35, units = "cm", dpi=300)
ggsave("../figures/phylogenetic_trees/Ecoli_tree_with_phylogroups_and_DefSystems_20230614_subsample750.svg",
       plot=Treesub750, width = 35, height=35, units = "cm", dpi=300)

######################
##Plot only for complete genomes
AssemblySummaryFilteredComplete<-subset(AssemblySummaryFiltered,
                                        AssemblySummaryFiltered$release_type == "Complete Genome")
CompleteSubsample<-sample(AssemblySummaryFilteredComplete$assembly_accession, 750)
TreeCompl750<-DrawSubsample(CompleteSubsample)
TreeCompl750

ggsave("../figures/phylogenetic_trees/Ecoli_tree_with_phylogroups_and_DefSystems_20230614_Complete750.png",
       plot=TreeCompl750, width = 35, height=35, units = "cm", dpi=300)
ggsave("../figures/phylogenetic_trees/Ecoli_tree_with_phylogroups_and_DefSystems_20230614_Complete750.svg",
       plot=TreeCompl750, width = 35, height=35, units = "cm", dpi=300)
