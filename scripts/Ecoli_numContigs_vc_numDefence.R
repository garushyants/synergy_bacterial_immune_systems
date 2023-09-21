library(stringr)
library(dplyr)
library(phytools)
library(ggplot2)

path<-getwd()
setwd(path)
setwd("../../20230920_Ecoli_number_of_contigs_vs_defence_systems/")

###read tree
tree<-read.tree("Ecoli_tree_rapidnj.rM2.treeshrink_corrected.nwk")
GenomesToKeep<-tree$tip.label
rep_str = c("'"='','.fna'='')
GenomesToKeep<-str_replace_all(GenomesToKeep,rep_str)

###read BV-BRC data
EcoliBVBRC<-read.csv("Ecoli_BVBRC_genome_20230920.csv", header=T)
AssemblySummary<-read.csv("assembly_summary_refseq_ecoli.csv", header=T)
AssemblySummary$IDtoMatch<-ifelse(AssemblySummary$paired_asm_comp !="", AssemblySummary$paired_asm_comp,
                                  AssemblySummary$assembly_accession)
DataToProcess<-merge(AssemblySummary[,c("assembly_accession","IDtoMatch")],
                     EcoliBVBRC[,c("Assembly.Accession","Contigs","Size","Plasmids")],
                     by.y="Assembly.Accession", by.x="IDtoMatch")
DataToProcessFinal<-subset(DataToProcess,
                           DataToProcess$assembly_accession %in% GenomesToKeep)

######Read defence data
EcoliDefense<-read.csv("complete_prophage_platon_26k_all7.csv")

EcoliDefenseWithContigsCounts<-subset(EcoliDefense,
                                      EcoliDefense$genome %in% DataToProcess$assembly_accession)

EcoliDefCounts<-EcoliDefenseWithContigsCounts %>% group_by(genome) %>%
  summarize(NumDefPred=n())

#####
#Plot
ForPlot<-merge(EcoliDefCounts,DataToProcess, by.x = "genome", by.y = "assembly_accession")
ForPlot$Plasmids[is.na(ForPlot$Plasmids)]<-0
ForPlot$allcontigs<-ForPlot$Contigs+ForPlot$Plasmids

###
CorCoef<-cor.test(ForPlot$Contigs, ForPlot$NumDefPred,
                  method="spearman",
                  alternative = "two.sided",
                  conf.level=0.99,
                  exact = F)
pvalue<-ifelse(CorCoef$p.value == 0, "p-value < 2.2e-16",
               paste("p-value =",format(CorCoef$p.value, 2, digits=3)))


Plot<-ggplot(data=ForPlot,
       aes(x=Contigs,
           y=NumDefPred))+
  geom_point(alpha=.2,
             shape=16,
             size=0.8)+
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
  ylab("Number of found defense systems")+
  theme_classic()+
  theme(axis.text = element_text(size=10, family="ArielMT"))
  
Plot


ggsave("Ecoli_numPredictedDefence_vs_numContigs.svg",
       plot=Plot,
       width=20, height=20, units="cm", dpi=300)
ggsave("Ecoli_numPredictedDefence_vs_numContigs.png",
       plot=Plot,
       width=20, height=20, units="cm", dpi=300)









