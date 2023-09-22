library(stringr)
library(dplyr)
library(phytools)
library(ggplot2)

path<-getwd()
setwd(path)


###read tree
tree<-read.tree("../data/Ecoli_tree_rapidnj.rM2.treeshrink_corrected.nwk")
GenomesToKeep<-tree$tip.label
rep_str = c("'"='','.fna'='')
GenomesToKeep<-str_replace_all(GenomesToKeep,rep_str)

###read contig count data
EcoliContigs<-read.csv("../data/Ecoli_numContigs_vs_numDefence/Ecoli_contigs.txt", sep="\t", header=F)
EcoliContigs$genome<-str_replace(EcoliContigs$V1,".fna","")
EcoliContigsFilt<-subset(EcoliContigs,
                         EcoliContigs$genome %in% GenomesToKeep)

AssemblySummary<-read.csv("../data/Ecoli_numContigs_vs_numDefence/assembly_summary_refseq_ecoli.csv", header=T)

######Read defence data
EcoliDefense<-read.csv("../data/26k_Ecoli_with_prophages.csv")

EcoliDefenseWithContigsCounts<-subset(EcoliDefense,
                                      EcoliDefense$genome %in% EcoliContigsFilt$genome &
                                        EcoliDefense$genome %in% AssemblySummary$assembly_accession)

EcoliDefCounts<-EcoliDefenseWithContigsCounts %>% group_by(genome) %>%
  summarize(NumDefPred=n())

#####
#Plot
ForPlot<-merge(EcoliDefCounts,EcoliContigs[,c("genome","V2")], by = "genome")


###
CorCoef<-cor.test(ForPlot$V2, ForPlot$NumDefPred,
                  method="spearman",
                  alternative = "two.sided",
                  conf.level=0.95,
                  exact = F)
pvalue<-ifelse(CorCoef$p.value == 0, "p-value < 2.2e-16",
               paste("p-value =",format(CorCoef$p.value, 2, digits=3)))


Plot<-ggplot(data=ForPlot,
             aes(x=V2,
                 y=NumDefPred))+
  geom_point(alpha=.1,
             shape=16,
             size=1.2)+
  geom_smooth(method = "lm",
              formula = y~x,
              color="#993404",
              fill="#fec44f",
              se = T,
              level=0.95)+
  annotate("text", label=paste("r = ",format(round(CorCoef$estimate, 2), nsmall = 2),
                               "\n",pvalue,
                               sep=""),
           x=1000,y=17,
           hjust=0,
           size=3.5,
           family="ArielMT")+
  xlab("Number of contigs")+
  ylab("Number of found defence systems")+
  theme_classic()+
  theme(axis.text = element_text(size=10, family="ArielMT"))

Plot


ggsave("../figures/Ecoli_numContigs_vs_numDefence/Ecoli_numPredictedDefence_vs_numContigs.svg",
       plot=Plot,
       width=20, height=18, units="cm", dpi=300)
ggsave("../figures/Ecoli_numContigs_vs_numDefence/Ecoli_numPredictedDefence_vs_numContigs.png",
       plot=Plot,
       width=20, height=18, units="cm", dpi=300)














