library(phytools)
#library(treeio)
library(castor)
library(ggtree)
library(stringr)
library(ggplot2)
library(ape)

path<-getwd()
setwd(path)


###Basic tree function
drawTreeBasic<-function(filename,rerootmiddle = F)
{
  Tree<-read.tree(filename)
  if (rerootmiddle) {
    Tree<-midpoint.root(Tree)
  }
  newtipnames<-Tree$tip.label
  rep_str = c("'"='')
  newtipnames<-str_replace_all(newtipnames,rep_str)
  Tree$tip.label<-newtipnames
  p<-ggtree(Tree, layout="fan", open.angle=10)
  myreturnlist<-list("tips"=newtipnames,
                     "tree"=Tree,
                     "plot"=p)
  return(myreturnlist)
}
###add taxonomy info on the plot
AddTaxInfoToTips<-function(AssemblyInfoFile,plot,legend=F)
{
  AssemblyInfo<-read.csv(AssemblyInfoFile, sep="\t", header=F)
  AssemblyInfo$genus<-word(AssemblyInfo$V8,1)
  unique(AssemblyInfo$genus)
  TipsDf<-AssemblyInfo[,c("V1","genus")]
  
  Plot<-plot %<+% TipsDf +
    geom_tippoint(aes(color=genus), show.legend = legend)+
    theme(legend.text=element_text(size=5))
  #geom_tiplab(size=1, as_ylab = T)
  
  return(Plot)
}

PseuTree<-drawTreeBasic("../data/pseutree_resavediTOL_newick.txt")

PseuTreeFigure<-AddTaxInfoToTips("../data/phylogenetic_trees_order_level/assembly_summary_genbank_pseu.txt",
                                 PseuTree$plot, legend =T)
PseuTreeFigure
BaciTreeFigure<-AddTaxInfoToTips("../data/phylogenetic_trees_order_level/assembly_summary_genbank_baci.txt",
                                 drawTreeBasic("../data/bacitree_resavediTOL_newick.txt")$plot, legend =T)
BaciTreeFigure
BurkTreeFigure<-AddTaxInfoToTips("../data/phylogenetic_trees_order_level/assembly_summary_genbank_burk.txt",
                                 drawTreeBasic("../data/burktree_resavediTOL_newick.txt")$plot, legend =T)
BurkTreeFigure
EnterTreeFigure<-AddTaxInfoToTips("../data/phylogenetic_trees_order_level/assembly_summary_genbank_enter.txt",
                                 drawTreeBasic("../data/entertree_resavediTOL_newick.txt")$plot, legend =F)
EnterTreeFigure


























