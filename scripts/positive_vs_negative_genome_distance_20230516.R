library(readxl)
library(dplyr)
library(phytools)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(castor)
library(reshape2)
library(writexl)


path<-getwd()
setwd(path)
setwd("../20230427_distances_on_genome/")

###Create folder to store results
folderForOutput<-"results20230610"

if (!dir.exists(folderForOutput)){
  dir.create(folderForOutput)
}else{
  print("dir exists")
}
######################################
######################################
###Compile data

###read Pagel results for E. coli dataset
Pagel_results<-read_xlsx("Ecoli_pagel_rescaled_all_results_20230329.xlsx")

Pagel_results$pair<-paste(Pagel_results$System.I,
                                    Pagel_results$System.II)
Pagel_results$pairrev<-paste(Pagel_results$System.II,
                                       Pagel_results$System.I)

#select only the ones about which we are absolutely sure (Bonferroni correction)
Pagel_results_signifPos<-subset(Pagel_results, (Pagel_results$Bonferroni=="Y" &
                                  Pagel_results$direction == 1))
#select negative pairs because they can be a good dataset to compare with
Pagel_results_signifNeg<-subset(Pagel_results, (Pagel_results$Bonferroni=="Y" &
                                                  Pagel_results$direction == -1))

################################
#We work later on with Enterobacteriales dataset, in particular
#with subset of full E.coli genomes, Shigellas are considered E.coli
###Read main file with data
SystemsLocation<-read.csv("merged_enter3.csv", header=T)
###Get final tree for Enterobacteriales to get read of contaminated and strange genomes
tree <- read.tree("entertree_resavediTOL_newick.txt")
tipnames<-tree$tip.label
###read Genbank assembly summary to filter E.coli from whole Enterobacteriales dataset 
AssemblySummary<-read.csv("Escherichia_assembly_summary_genbank.txt", header = F,
                          sep="\t")
#remove contaminated
SystemsLocationNoCont<-subset(SystemsLocation, SystemsLocation$genome %in% tipnames)


#################################
###Because chromosomes are circular I need to know all contigs lengths to be able to calculate distances correctly
###This lengths are precalculated from fasta files with awk oneliner 
ChromosomeLengths<-read.csv("../20230215_download_genomes/enter_genomes_chromosome_lengths.tsv",
                            header = F, sep="\t")
ChromosomeLengthsSep<-separate(ChromosomeLengths,V1,into = c("seqid","description"),sep = " ",extra = "merge")
#In addition get info about location on plasmid
ChromosomeLengthsSep$plasmid<-as.integer(grepl("plasmid",ChromosomeLengthsSep$description))


###Merge main dataset with the length info
DefSystemsWithChrInfo<-merge(SystemsLocationNoCont, ChromosomeLengthsSep,
                                  by="seqid", all.x=T)

#Subset only E.coli and Shigella
SystemsLocationEcoli<-subset(DefSystemsWithChrInfo, DefSystemsWithChrInfo$genome %in%
                               AssemblySummary$V1 &
                               #remove to genomes with incorrect Genbank IDs,
                               #because they are not E.coli
                               !(DefSystemsWithChrInfo$genome %in% c("GCA_000210475.1",
                                                                     "GCA_000009565.2")))
#################################
#Read data from E.coli dataset to compare results
###systems
Ecoli26KSystemsLocation<-read.csv("complete_prophage_platon_26k_all7.csv", header=T)
#tree
Ecoli26ktree <- read.tree("Ecoli_tree_rapidnj.rM2.treeshrink_corrected.nwk")
Finaltips<-gsub(".fna","",Ecoli26ktree$tip.label)
Ecoli26KSystemsLocationFiltered<-subset(Ecoli26KSystemsLocation, 
                                        Ecoli26KSystemsLocation$genome %in% Finaltips)

#############################################
#############################################
###Start analysis

#############
#Check how often systems are present on chromosomes
ChrVsPlsmPlot<-ggplot(SystemsLocationEcoli)+
  geom_bar(aes(x=defense_system2,
               y = ..count..,
               group=as.factor(plasmid),
               fill=as.factor(plasmid)),
           position="dodge")+
  scale_fill_manual(name="Location",
                    labels=c("chromosome","plasmid"),
                    values=c("#4daf4a","#984ea3"))+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ChrVsPlsmPlot
#Interesting systems are CBASS I, Gabihja,tmn, PifA and ppl
#Also there are some RMs on plasmids, but rarely

ggsave("EcoliComplete_ChromosmeVsPlasmid.png",
       plot=ChrVsPlsmPlot, path=folderForOutput,
       width=40, height=13, dpi=300, units="cm")

###Draw the same for 26k Ecoli

Ecoli26kChrVsPlsmPlot<-ggplot(Ecoli26KSystemsLocationFiltered)+
  geom_bar(aes(x=defense_system,
               y = ..count..,
               group=as.factor(seqid_type),
               fill=as.factor(seqid_type)),
           position="dodge")+
  scale_fill_manual(name="Location",
                    values=c("#4daf4a","#984ea3","#969696"))+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
Ecoli26kChrVsPlsmPlot
#Results are similar for the whole dataset
#(if unknown location is included)

ggsave("Ecoli26k_ChromosmeVsPlasmid.png",
       plot=Ecoli26kChrVsPlsmPlot, path=folderForOutput,
       width=40, height=13, dpi=300, units="cm")

###check the same plasmid thing for the whole Enterobacteriales dataset
EnterChrVsPlsmPlot<-ggplot(DefSystemsWithChrInfo)+
  geom_bar(aes(x=defense_system2,
               y = ..count..,
               group=as.factor(plasmid),
               fill=as.factor(plasmid)),
           position="dodge")+
  scale_fill_manual(name="Location",
                    labels=c("chromosome","plasmid"),
                    values=c("#4daf4a","#984ea3"))+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
EnterChrVsPlsmPlot
#behavior is similar but systems on chromosomes are more rare except for RMs


ggsave("Enterobacteriales_ChromosmeVsPlasmid.png",
       plot=EnterChrVsPlsmPlot, path=folderForOutput,
       width=40, height=13, dpi=300, units="cm")

####Check how often systems are within prophages for 26k E.coli dataset 
####for which we have such data
Ecoli26kChrVsPrPhPlot<-ggplot(Ecoli26KSystemsLocationFiltered)+
  geom_bar(aes(x=defense_system,
               y = ..count..,
               group=as.factor(prophage_within),
               fill=as.factor(prophage_within)),
           position="dodge")+
  scale_fill_manual(name="prophage",
                    values=c("#e41a1c","#377eb8","#969696"))+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
Ecoli26kChrVsPrPhPlot
#pure prophages are rare
#BstA is interesting in this sense
#and also ietAS sometimes close
ggsave("Ecoli26k_ChromosmeVsProphage.png",
       plot=Ecoli26kChrVsPrPhPlot, path=folderForOutput,
       width=40, height=13, dpi=300, units="cm")


###################################################
####We know now that most systems are on chromosomes and start asking questions about distances

####Let's calculate some background distances as it was done previously
#Calculate distances between all systems
ContigsToConsiderAll<-SystemsLocationEcoli %>% group_by(seqid) %>% summarize(count=length(defense_system2))
#remove contigs with one system
ContigsToConsider<-subset(ContigsToConsiderAll, ContigsToConsiderAll$count >1)
#get data only for contigs with more than one system
SetToCalculateDistances<-subset(SystemsLocationEcoli, SystemsLocationEcoli$seqid %in% ContigsToConsider$seqid)

#calculate all distances to the next system
NearestDistances<-SetToCalculateDistances %>% 
  group_by(seqid) %>%
  arrange(seqid,start) %>%
  mutate(distanced = start - lag(end, n=1, default =NA), 
         distancel = V2 +lag(end, n=1, default =NA)-start) %>%
  mutate(distance = pmin(distancel,distanced)) %>%#to clean the lengths
  #with code below I select the distance to closest system,
  #but this approach produces duplicates
  mutate(mindist=ifelse(is.na(lead(distance)), distance,
                        ifelse(distance>lead(distance,n=1, default=NA),
                        lead(distance,n=1, default=NA),distance)),
         paired=ifelse(is.na(lead(distance)), lag(defense_system2),
                       ifelse(distance>lead(distance,n=1, default=NA),
                        lead(defense_system2),
                       lag(defense_system2)))) %>%
  mutate(pairOrdered=paste(pmin(defense_system2,paired),pmax(defense_system2,paired)))


NearestDistancesShort<-NearestDistances[,c("genome","seqid","mindist","pairOrdered")]
NearestDistancesNoDupl<-NearestDistancesShort[!duplicated(NearestDistancesShort),]

NearestDistancesClean<-subset(NearestDistancesNoDupl, !is.na(NearestDistancesNoDupl$mindist) &
                                !(NearestDistancesNoDupl$pairOrdered %in% c(Pagel_results_signifPos$pair,
                                                                      Pagel_results_signifPos$pairrev)))
#############!!!!!Check again it seems that I am missing distances for first element in each first element



#Distance to nearest characterized system in dataset
Allp1<-ggplot(NearestDistancesClean, aes(x=mindist,
                                         after_stat(count*100/sum(count))))+
  geom_histogram(bins=100,
                 alpha = 0.7, fill = "#023858",color = "#023858")+
  geom_vline(xintercept=median(NearestDistancesClean$mindist),
             color="dark green", alpha=.6)+
  ylab("% of all distances")+
  xlab("distance")+
  scale_x_continuous(limits= c(-30000,2000000))+
  theme_bw()
Allp1
##############################################################
################Get distances for pairs
#############################
#Get distances for pairs of interest

getPairDistance<-function(sysI, sysII)
{
  #sysI<-"RM I"
  #sysII<-"RM II"
  SystemPairsDfAll<-subset(SetToCalculateDistances, SetToCalculateDistances$defense_system2 %in% c(sysI,sysII))
  #filter contigs with only one system
  PairContigsAll<-SystemPairsDfAll %>% group_by(seqid,genome) %>%summarize(count=length(unique(defense_system2)))
  PairContigsForCount<-subset(PairContigsAll, PairContigsAll$count>1)
  
  SystemPairsDf<-subset(SystemPairsDfAll, SystemPairsDfAll$seqid %in% PairContigsForCount$seqid)
  PairNearestDistances<-SystemPairsDf %>% 
    group_by(seqid) %>%
    arrange(seqid,start) %>%
    mutate(distanced = start - lag(end, n=1, default =NA),
           distancel = ChromosomeLengthsSep[ChromosomeLengthsSep$seqid == seqid,]$V2 +lag(end, n=1, default =NA)-start) %>%
    mutate(distance = pmin(distancel,distanced))
  PairNearestDistances$System.I<-rep(sysI, length(PairNearestDistances$seqid))
  PairNearestDistances$System.II<-rep(sysII, length(PairNearestDistances$seqid))
  #remove distances between systems of the same type
  PairNearestDistances$prev<-lag(PairNearestDistances$defense_system2)
  PairNearestDistancesClear<-subset(PairNearestDistances,
                                    PairNearestDistances$defense_system2 != PairNearestDistances$prev &
                                      !is.na(PairNearestDistances$distance))
  return(PairNearestDistancesClear)
}

###execute
DistancesListOfDFs<-apply(Pagel_results_signifPos, 1, function(a) getPairDistance(a[3], a[2]))
PairDistancesDf<-do.call("rbind",DistancesListOfDFs)
PairDistancesDf$pair<-paste(PairDistancesDf$System.II,
                            PairDistancesDf$System.I)
write_xlsx(PairDistancesDf, path=paste(folderForOutput,"/Distances_between_positive_significant_pairs.xlsx",sep=""))

#Draw all systems together
AllPairsPlotp1<-ggplot(PairDistancesDf, aes(x=distance,color=pair,fill=pair,
                                         after_stat(count*100/sum(count))))+
  geom_histogram(bins=100,
                 alpha = 0.4)+
  geom_vline(xintercept=median(PairDistancesDf$distance),
             color="dark green", alpha=.4)+
  ylab("% of all distances")+
  scale_x_continuous(limits= c(-30000,2000000))+
  theme_bw()+
  theme(legend.position = "none")
AllPairsPlotp1

AllPairsPlotCombined<-ggarrange(Allp1,
          AllPairsPlotp1,
          nrow=2,
          align="v")
AllPairsPlotCombined

ggsave("EcoliEnterAllPositive_vs_allclosest.png",
       path=folderForOutput,
       AllPairsPlotCombined, width=20, height=15, units="cm",
       dpi=300)


NearestDistancesClean$pair<-rep(" All", length(NearestDistancesClean$seqid))

write_xlsx(NearestDistancesClean, path=paste(folderForOutput,"/Distances_to_closest_system_control_All.xlsx", sep=""))

###Ploting distribution of distances
LongPlot<-ggplot()+
  geom_boxplot(data=NearestDistancesClean, aes(x=mindist,
                                               y=pair,
                                                 fill=pair),
              fill = "#023858")+
  geom_vline(xintercept=median(NearestDistancesClean$mindist),
             color="dark green", alpha=.4)+
  geom_boxplot(data=PairDistancesDf, 
              aes(x=distance,
                  y=pair,
                  fill= pair))+
  xlim(c(-10000,500000))+
  theme_bw()+
  theme(legend.position = "none")
LongPlot

ggsave("EcoliEnter_AllPositive_long_barplot.png",
       path=folderForOutput,
       LongPlot, width=10, height=40, units="cm",
       dpi=300)
#There can be any distances, but there are systems that are located close to each other

#How to look into that?

#Check whick pairs are close..
#some indeed are
SystemsMedian<-PairDistancesDf %>% group_by(System.I,System.II) %>%
  summarize(median=median(distance),
            count = length(distance),
            test=wilcox.test(distance,NearestDistancesClean$mindist, alternative="less")$p.value)
SystemsSignif<-subset(SystemsMedian, SystemsMedian$test*length(SystemsMedian$System.I) < 0.001 &
                        SystemsMedian$count > 10)

##save to file
write_xlsx(SystemsSignif, path=paste(folderForOutput,"PositiveSystemsSignifCloser.xlsx", sep="/"))

ForPlotSubset<-subset(PairDistancesDf, paste(PairDistancesDf$System.I,PairDistancesDf$System.II) %in% 
                        paste(SystemsSignif$System.I,SystemsSignif$System.II))


p2<-ggplot(ForPlotSubset, 
           aes(x=distance, fill= paste(System.I, System.II)))+
  geom_histogram(position="identity", alpha=.5, binwidth=10000)+
  geom_vline(xintercept=median(NearestDistancesClean$mindist))+
  facet_wrap(~pair, scales = "free_y")+
  xlim(c(-10000,500000))+
  theme_bw()+
  theme(legend.position = "none")
p2

ggsave("EcoliEnter_ByPositiveSystem_signifless.png",
       path=folderForOutput,
       p2, width=20, height=20, units="cm",
       dpi=300)

#Do I have to plot on the tree for tmn Gabija and Druantia III Zorya II
#where those are found
#DndABCDE and DndFGH is an obvious case

DrawPairTree<-function(sysI,sysII){
  #sysI<-"Zorya II"
  #sysII<-"Druantia III"
  #subset initial system locations
  PairSystemLocationLong<-subset(SystemsLocationEcoli,
                                 SystemsLocationEcoli$defense_system2 %in% c(sysI,sysII))[,c("genome","seqid","defense_system2")]
  
  PairSystemLocationLongForTree<-PairSystemLocationLong %>% 
    group_by(genome) %>%
    dplyr::summarise(system2 = ifelse(length(unique(defense_system2))>1,"both",defense_system2),
                     system = ifelse(length(unique(defense_system2))>1,0,
                                     ifelse(defense_system2[1] == sysI,1,2)),
                     chromosomes2= ifelse(length(unique(seqid))>1,"different","same"),
                     chromosomes = ifelse(length(unique(seqid))>1,0,1))
  names(PairSystemLocationLongForTree)[1]<-"label"

  #subset distances
  SelectedPairDistances<-subset(PairDistancesDf, PairDistancesDf$System.I %in% c(sysI,sysII) &
                                  PairDistancesDf$System.II %in% c(sysI,sysII))[,c("genome","seqid","plasmid","distance")]
  SelectedPairDistances$logdist<-log10(SelectedPairDistances$distance)
  SelectedPairDistances$plasmid2<-ifelse(SelectedPairDistances$plasmid == 1, "plasmid", "chromosome")
  names(SelectedPairDistances)[1]<-"label"
  #subset tree
  LeavesToKeep<-unique(SystemsLocationEcoli$genome)
  subtree<-get_subtree_with_tips(tree, only_tips = LeavesToKeep)$subtree
  
  #draw tree
  EcolFromEnterDataSubtree<-ggtree(subtree, size =0.1, color = "#525252")


  PairPlot<-EcolFromEnterDataSubtree +
    geom_facet(data=PairSystemLocationLongForTree,
               panel="system",
               geom=geom_tile,
               mapping = aes(x=system, fill=system2),
               width=0.7) +
    scale_fill_manual(values = c("#6a3d9a","#e31a1c","#1f78b4"),na.value="white")+
    new_scale_fill()+
    ##
    geom_facet(data=PairSystemLocationLongForTree,
               panel = "chromosome",
               geom=geom_tile,
               mapping=aes(x=chromosomes, fill=chromosomes2),
               width=0.7)+
    scale_fill_manual(values = c("#8c510a","#01665e"), name="chromosomes")+
    ##
    new_scale_fill()+
    geom_facet(data=SelectedPairDistances,
               panel="location",
               geom=geom_tile,
               mapping=aes(x=plasmid, fill=plasmid2),
               width=0.7)+
    scale_fill_manual(values = c("#fdb462","#b3de69"), name="location")+
    geom_facet(
      panel="Distance",
      data = SelectedPairDistances,
      geom = geom_col,
      aes(x = distance),
      orientation = 'y',
      color = "#6a3d9a",
      width=1)+
    theme_tree2()
  PairPlot2<-facet_widths(PairPlot, widths = c(1, 0.3,0.2,0.2,0.7))
  return(PairPlot2)
}

#draw and save
for (i in 1:nrow(SystemsSignif))
{
  treeplot<-DrawPairTree(SystemsSignif[i,"System.I"],SystemsSignif[i,"System.II"])
  ggsave(paste(SystemsSignif[i,"System.I"],SystemsSignif[i,"System.II"],"tree_with_distances.png",sep="_"),
         plot=treeplot,path=folderForOutput,
         width=30,height=40,dpi =400,units="cm")
}

#################################################
#################################################

























