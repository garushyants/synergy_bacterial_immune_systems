library(phytools)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(dplyr)
library(tidyr)
library(writexl)

path<-getwd()
setwd(path)


#####
#Read clades
setwd("./data/Ecoli_phylogroups/")
treesfiles<-c("Ecoli_A.nwk","Ecoli_B1.nwk","Ecoli_B21.nwk",
              "Ecoli_B22.nwk","Ecoli_C.nwk","Ecoli_E1.nwk",
              "Ecoli_E2.nwk")
treenames<-c("A","B1","B21","B22","C","E1","E2")

PhyloTrees=lapply(treesfiles,read.tree)
names(PhyloTrees)<-treenames

#Read tips for clades
PhyloTips=lapply(PhyloTrees, function(x) x$tip.label)

PhyloTipsDf<-data_frame(id = names(PhyloTips), PhyloTips) %>%
  unnest(cols = c(PhyloTips))

#save df with tips


PhyloTipsVec<-PhyloTipsDf$PhyloTips

#Read weights
weightsfiles<-c("Ecoli_A_leaf_weigts.tab","Ecoli_B1_leaf_weigts.tab",
                "Ecoli_B21_leaf_weigts.tab","Ecoli_B22_leaf_weigts.tab",
                "Ecoli_C_leaf_weigts.tab","Ecoli_E1_leaf_weigts.tab",
                "Ecoli_E2_leaf_weigts.tab")
TipWeights=lapply(weightsfiles, function(x) read.csv(x,sep="\t",header=F))
names(TipWeights)<-treenames

TipsWeightsDf<-data_frame(id = names(TipWeights), TipWeights) %>%
  unnest(cols = c(TipWeights))

#Read defense systems data
setwd("../")

EcoliDefense<-read.csv("ecoli_filtered.csv", header = T)

EcoliDefense$location<-ifelse(EcoliDefense$prophage_within=="full",
                              paste("prophage",EcoliDefense$seqid_type),
                              EcoliDefense$seqid_type)

EcoliDefense$length<-abs(EcoliDefense$end - EcoliDefense$start) +1

DefenseBySystem<-EcoliDefense %>% group_by(genome,immune_system) %>%
  count(immune_system)
DefenseBySystemWide<-as.data.frame(DefenseBySystem %>% 
                                     pivot_wider(names_from = immune_system, values_from = n))
DefenseBySystemWide[is.na(DefenseBySystemWide)]<-0

#filter
DefenseBySystemWideOnly7phylo<-subset(DefenseBySystemWide,
                                      DefenseBySystemWide$genome %in% PhyloTipsVec)
#Get to binary
DefenseBinary<-DefenseBySystemWideOnly7phylo[,c(2:ncol(DefenseBySystemWideOnly7phylo))]
DefenseBinary[is.na(DefenseBinary)]<-0
DefenseBinary[DefenseBinary>0]<-1



#drop columns with cutoff
cutoff<-0.005 #the system has to be present in 0.5% of genomes

DefenseBySystem7phyloFilt<-DefenseBinary[,colSums(DefenseBinary) > cutoff*nrow(DefenseBinary)]

#get system order for vizualization
SysOrderDf<-data.frame(System=colnames(DefenseBySystem7phyloFilt[,c(1:71)]),
                       count=colSums(DefenseBySystem7phyloFilt[,c(1:71)]))
SystemOrder<-SysOrderDf[order(SysOrderDf$count),]$System
SysOrderDf$System<-factor(SysOrderDf$System,levels=SystemOrder)

########
#Finally can work with dataset
DefenseBySystem7phyloFilt$genome<-DefenseBySystemWideOnly7phylo$genome
DefenseBySystemWithPhylo<-merge(DefenseBySystem7phyloFilt,PhyloTipsDf, by.x="genome",
                                by.y = "PhyloTips")
####
#the same but without cutoff
DefenseBySystemWithPhyloWithoutCutoff<-merge(DefenseBySystemWideOnly7phylo,PhyloTipsDf, by.x="genome",
                                by.y = "PhyloTips")
#save 
# write_xlsx(DefenseBySystemWithPhyloWithoutCutoff, path="Ecoli_7phylogroups_DefenseSystems_counts_raw.xlsx")
# 

#####Transform back to long
DefWithoutCutoffLong<-gather(DefenseBySystemWithPhyloWithoutCutoff,immune_system,count,
       `CRISPR I-E`:`Druantia II`)

DefFiltLong<-gather(DefenseBySystemWithPhylo,immune_system,count,
                             `CRISPR I-E`:`Lamassu I`)

CountPerPhylogroup<-DefWithoutCutoffLong %>% group_by(id)%>%
  summarise(total=length(unique(genome)))
CountSystemsPerPhylo<-DefFiltLong %>% group_by(id,immune_system)%>%
  summarise(systemtotal=sum(count))


ForPlotRaw<-merge(CountSystemsPerPhylo,CountPerPhylogroup,by="id")
ForPlotRaw$perc<-ForPlotRaw$systemtotal*100/ForPlotRaw$total
ForPlotRaw$immune_system<-factor(ForPlotRaw$immune_system,
                                     levels=SystemOrder)
###Systems counts by genome per phylogroup
SysPerGenomePlot<-DefWithoutCutoffLong %>% group_by(id,genome)%>%
  summarise(syspergenome=sum(count)) %>%
  group_by(id)%>%
  mutate(median=median(syspergenome))

#do stat
SysPerGenomePlot %>% group_by(id)%>%
  summarise(median = median(syspergenome),
            sd=sd(syspergenome))

doWilcoxTest<-function(group1, group2)
{
  
  Res<-wilcox.test(pull(SysPerGenomePlot[SysPerGenomePlot$id==group1,][3]),
              pull(SysPerGenomePlot[SysPerGenomePlot$id==group2,][3]),
              alternative = "two.sided")
  return(Res)
}

doWilcoxTest("A","C")
doWilcoxTest("A","E2")
doWilcoxTest("A","B21")
doWilcoxTest("B1","B21")
doWilcoxTest("B22","B21")
doWilcoxTest("C","B22")
doWilcoxTest("E1","E2")

#Plot
SysPerGenomePlotObj<-ggplot(data=SysPerGenomePlot,
       aes(x=syspergenome, group=id, fill=id))+
  geom_histogram(bins=15,
           color="white")+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a",
                             "#33a02c","#fb9a99","#cab2d6",
                             "#6a3d9a"),
                    name = "")+
  facet_wrap(~id, ncol=1, scales = "free_y")+
  geom_vline(aes(xintercept = median),
             linewidth=1.5,
             color="#67000d")+
  theme_classic()+
  scale_x_continuous(breaks=seq(1,20,1), limits=c(0,18))+
  xlab("Defence systems per genome")
SysPerGenomePlotObj  
####Calculate sum length of systems per genome on average
EcoliDefLengthForPlot<-subset(EcoliDefense,
                         EcoliDefense$genome %in% unique(SysPerGenomePlot$genome))
PreDefLengthForPlot<-merge(EcoliDefLengthForPlot, SysPerGenomePlot[,c(1,2)], by="genome")

#the problem is that some systems are very short, because the genomes are not full
#in some cases I have systems for which I need to know the genome length
#the easiest way is to remove outliers
minmaxlength<-quantile(PreDefLengthForPlot$length, probs=c(0.01, 0.99))
PreDefLengthForPlot<-subset(PreDefLengthForPlot,
                            PreDefLengthForPlot$length > minmaxlength[1] &
                              PreDefLengthForPlot$length < minmaxlength[2])

DefLengthForPlot<- PreDefLengthForPlot %>% group_by(id,genome)%>%
  summarise(lengthpergenome=sum(length))%>%
  group_by(id)%>%
  mutate(median=median(lengthpergenome))

#Plot
LengthPerGenomePlot<-ggplot(data=DefLengthForPlot,
       aes(x=lengthpergenome, group=id, fill=id))+
  geom_histogram(bins=30,
                 color="white")+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a",
                             "#33a02c","#fb9a99","#cab2d6",
                             "#6a3d9a"),
                    name = "")+
  facet_wrap(~id, ncol=1, scales = "free_y")+
  geom_vline(aes(xintercept = median),
             linewidth=1.5,
             color="#67000d")+
  theme_classic()+
  scale_x_continuous(breaks=seq(0,150000,10000))+
  xlab("Defence systems length per genome")
LengthPerGenomePlot

SysPerGenomeNumAndLength<-ggarrange(SysPerGenomePlotObj, LengthPerGenomePlot,
          ncol=2,
          common.legend = T,
          legend = "right")

#save
ggsave("../figures/Ecoli_phylogroups/SystemsNumAndLengthPerGenomeByPhylogroup_raw.png", 
       plot=SysPerGenomeNumAndLength,
       width=25, height = 28, units="cm",dpi=300)

#correlation between number and mean length of system
MedianDefLength<- PreDefLengthForPlot %>% group_by(id,genome)%>%
  summarise(medianlengthpergenome=median(length))
MedianLenDefVsNum<-merge(MedianDefLength, SysPerGenomePlot, by=c("genome","id"))

MedianLenDefVsNumPlot<-ggplot(data=MedianLenDefVsNum,
       aes(x=syspergenome,
           y=medianlengthpergenome,
           fill=id))+
  geom_point(shape=21,
             size=3,
             stroke=0.4,
             color="#bdbdbd")+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a",
                             "#33a02c","#fb9a99","#cab2d6",
                             "#6a3d9a"),
                    name = "")+
  scale_x_continuous(name="# of defence systems per genome",
                     limits=c(0,18),
                     breaks =seq(0,20,1))+
  scale_y_continuous(breaks=seq(0,20000,1000), 
                     name="median defence system length per genome")+
  stat_poly_eq(formula = y ~ x, data=MedianLenDefVsNum, inherit.aes = F,
               aes(x = syspergenome, y = medianlengthpergenome, label = paste(..rr.label..,
                                                                        ..p.value.label..,
                                                                        sep = "~~~~")), 
               parse = TRUE)+
  geom_smooth(data=MedianLenDefVsNum,
              aes(x=syspergenome,
                  y=medianlengthpergenome),
              inherit.aes=F,
              method='lm',
              color = "#49006a")+
  theme_classic2(base_family = "ArielMT")
MedianLenDefVsNumPlot

ggsave("../figures/Ecoli_phylogroups/SystemsNumVsMedianLengthPerGenome_raw.png", 
       plot=MedianLenDefVsNumPlot,
       width=25, height = 28, units="cm",dpi=300)



####Plot Raw counts
LongPlotFacet<-ggplot(ForPlotRaw,aes(y=immune_system,x=perc, fill=id))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a",
                             "#33a02c","#fb9a99","#cab2d6",
                             "#6a3d9a"))+
  scale_x_continuous(breaks=c(0,50,100))+
  facet_wrap(~id,nrow=1)+
  xlab("%")+
  ylab("")+
  theme_classic()+
  theme(panel.grid.minor = element_line(linewidth = 0.2),
        panel.grid.major = element_line(linewidth = 0.2),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        legend.position = "none")
LongPlotFacet

SysOrderDf$title<-rep("# of genomes",length(SysOrderDf$count))
SystemCountPlot<-ggplot(SysOrderDf,
       aes(y=System,x=-count))+
  geom_bar(stat="identity", position=position_dodge())+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks=seq(-25000,0,10000),
                     labels=seq(25000,0,-10000),
                     limits=c(-23500,5000))+
  geom_text(aes(y=System, label=count), x=500, hjust=0,
            size=3)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.2),
        panel.grid.major.y = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), 
        legend.position = "none") +
  facet_wrap(~title)
SystemCountPlot

LongRawToSave<-ggarrange(SystemCountPlot,
                         LongPlotFacet,
                         ncol = 2,
                         widths = c(1,4))
LongRawToSave

ggsave("../figures/Ecoli_phylogroups/SystemsByPhylogroup_LongPlot_raw_cutoff005.png", plot=LongRawToSave,
       width=30, height = 28, units="cm",dpi=300)
ggsave("../figures/Ecoli_phylogroups/SystemsByPhylogroup_LongPlot_raw_cutoff005.svg", plot=LongRawToSave,
       width=30, height = 28, units="cm",dpi=300)

###########
#Do statistical tests for selected systems
#The important value to report is effect size
do_chisquare<-function(system)
{
  #system<-"Druantia III"
  SysLong<-subset(ForPlotRaw, ForPlotRaw$immune_system == system)
  SysLong$absent<-SysLong$total-SysLong$systemtotal
  
  #Do Chi-Squared test for homogeneity
  data.table<-rbind(system=SysLong$systemtotal,absent=SysLong$absent)
  #data.table2<-rbind(sysperc=SysLong$perc, absperc=(100-SysLong$perc))
  TestRes<-chisq.test(data.table,
                      simulate.p.value = T,
                      B=2000)
  OddsRatio<-TestRes$observed[1,]/TestRes$expected[1,]
  names(OddsRatio)<-SysLong$id
  return(c(def_sys=system,TestRes$statistic,p.value=TestRes$p.value,OddsRatio))

}

write.csv(rbind(do_chisquare("Zorya II"),
do_chisquare("Druantia III"),
do_chisquare("Retron I-C"),
do_chisquare("AbiE"),
do_chisquare("CRISPR I-F"),
do_chisquare("CRISPR I-E"),
do_chisquare("Thoeris I"),
do_chisquare("Septu I"),
do_chisquare("qatABCD"),
do_chisquare("BREX I"),
do_chisquare("PsyrTA"),
do_chisquare("RM IV"),
do_chisquare("RM IIG")),"../data/Ecoli_phylogroups/chisquared_systems_difference.csv",
quote = F, row.names = F)

#########
##Add info about location

DefFiltLongWithLoc<-merge(DefFiltLong,EcoliDefense[,c("genome","immune_system","location")],
      by=c("genome","immune_system"))
SysSumPhylo<-DefFiltLongWithLoc%>% group_by(id,immune_system)%>%
  summarise(sum=sum(count))
PreForPlotRawWithLoc<-DefFiltLongWithLoc%>%
  group_by(id,immune_system,location) %>%
  summarise(loccount=n())
ForPlotRawWithLoc<-merge(PreForPlotRawWithLoc,SysSumPhylo, by=c("id","immune_system"))

ForPlotRawWithLoc$perc<-ForPlotRawWithLoc$loccount*100/ForPlotRawWithLoc$sum
ForPlotRawWithLoc$immune_system<-factor(ForPlotRawWithLoc$immune_system,
                                         levels=SystemOrder)

##plot
SystemsByPhylogroupByLocation<-ggplot(data=ForPlotRawWithLoc)+
  geom_col(aes(y=immune_system,x=perc,fill=location,
               color=location))+
  facet_wrap(~id, nrow=1)+
  scale_fill_manual(values=c("#742c24",
                             "#ebddd3",
                             "#ba4535",
                             "#88727b",
                             "#feb24c",
                             "#ffeda0"))+
  scale_color_manual(values=c("#742c24",
                              "#ebddd3",
                              "#ba4535",
                              "#88727b",
                              "#feb24c",
                              "#ffeda0"))+
  xlab("%")+
  ylab("")+
  scale_x_continuous(breaks=c(0,50,100))+
  theme_classic()+
  theme(panel.grid.minor = element_line(linewidth = 0.2),
        panel.grid.major = element_line(linewidth = 0.2),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"))
SystemsByPhylogroupByLocation

ggsave("../figures/Ecoli_phylogroups/SystemsByPhylogroupByLocation_LongPlotFacet_raw_cutoff005.png", plot=SystemsByPhylogroupByLocation,
       width=25, height = 28, units="cm",dpi=300)
ggsave("../figures/Ecoli_phylogroups/SystemsByPhylogroupByLocation_LongPlotFacet_raw_cutoff005.svg", plot=SystemsByPhylogroupByLocation,
       width=25, height = 28, units="cm",dpi=300)


###Do tests

###########################################################
###########################################################
###Do the same but with weights

DfForPlotPercOfGenomes<-data.frame(Counts=NULL,
                      Systems=NULL,
                      Phylogroup=NULL,
                      Perc=NULL)

for (phylogroup in treenames)
{
  #phylogroup<-"C"
  GroupWeights<-subset(TipsWeightsDf, TipsWeightsDf$id == phylogroup)
  sumWeights<-sum(GroupWeights$V2)
  GroupDefense<-subset(DefenseBySystemWithPhylo, DefenseBySystemWithPhylo$id == phylogroup)
  Combined<-merge(GroupDefense,GroupWeights, by.x ="genome",
                  by.y = "V1")
  #multiply by weights
  syscols<-colnames(Combined)[2:(ncol(Combined)-3)]
  MultipliedCombined<-Combined %>% mutate(across(syscols, ~.*V2))
  GroupValues<-colSums(MultipliedCombined[,c(2:(ncol(Combined)-3))])
  OutDf<-data.frame(Counts=GroupValues,Systems=syscols,Phylogroup=rep(phylogroup,length(syscols)))
  OutDf$Perc<-OutDf$Counts*100/sumWeights

  DfForPlotPercOfGenomes<-rbind(DfForPlotPercOfGenomes,OutDf)
}
# ###
# write_xlsx(DfForPlotPercOfGenomes, path = "26k_Ecoli_SystemsByPhylogroup_data.xlsx")
DfForPlotPercOfGenomes$Systems<-factor(DfForPlotPercOfGenomes$Systems,
                                       levels=SystemOrder)

LongPlotFacetWeights<-ggplot(DfForPlotPercOfGenomes,aes(y=Systems,x=Perc, fill=Phylogroup))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a",
                             "#33a02c","#fb9a99","#cab2d6",
                             "#6a3d9a"))+
  scale_x_continuous(breaks=c(0,50,100))+
  facet_wrap(~Phylogroup,nrow=1)+
  xlab("%")+
  ylab("")+
  theme_classic()+
  theme(panel.grid.minor = element_line(linewidth = 0.2),
        panel.grid.major = element_line(linewidth = 0.2),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        legend.position = "none")
LongPlotFacetWeights
# 
ggsave("../figures/Ecoli_phylogroups/SystemsByPhylogroup_LongPlotFacet_weight_cutoff005.png", plot=LongPlotFacetWeights,
       width=15, height = 28, units="cm",dpi=300)
ggsave("../figures/Ecoli_phylogroups/SystemsByPhylogroup_LongPlotFacet_weight_cutoff005.svg", plot=LongPlotFacetWeights,
       width=15, height = 28, units="cm",dpi=300)

#######################################
#Long Facet but based on systems location


DefSysLocWithPhylo<-merge(EcoliDefense[,c("genome","immune_system","location")],TipsWeightsDf, by.x="genome",
                                by.y = "V1")
DefSysLocWithPhyloFilt<-subset(DefSysLocWithPhylo,DefSysLocWithPhylo$immune_system %in% colnames(DefenseBySystem7phyloFilt))
SysSumPhyloWithW<-DefSysLocWithPhyloFilt%>% group_by(id,immune_system)%>%
  summarise(sum=sum(V2))

DefSysLocWithPhyloPerc<-merge(DefSysLocWithPhyloFilt,SysSumPhyloWithW, by=c("id","immune_system"))
DefSysLocWithPhyloPerc$perc<-DefSysLocWithPhyloPerc$V2*100/DefSysLocWithPhyloPerc$sum

SystemsByPhylogroupByLocationWeights<-ggplot(data=DefSysLocWithPhyloPerc)+
  geom_col(aes(y=immune_system,x=perc,fill=location,
               color=location))+
  facet_wrap(~id, nrow=1)+
  scale_fill_manual(values=c("#742c24",
                             "#ebddd3",
                             "#ba4535",
                             "#88727b",
                             "#feb24c",
                             "#ffeda0"))+
  scale_color_manual(values=c("#742c24",
                             "#ebddd3",
                             "#ba4535",
                             "#88727b",
                             "#feb24c",
                             "#ffeda0"))+
  xlab("%")+
  ylab("")+
  scale_x_continuous(breaks=c(0,50,100))+
  theme_classic()+
  theme(panel.grid.minor = element_line(linewidth = 0.2),
        panel.grid.major = element_line(linewidth = 0.2),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"))
SystemsByPhylogroupByLocationWeights

# ggsave("SystemsByPhylogroupByLocation_LongPlotFacetWeights_cutoff005.png", plot=SystemsByPhylogroupByLocationWeights,
#        width=22, height = 30, units="cm",dpi=300)
# ggsave("SystemsByPhylogroupByLocation_LongPlotFacetWeights_cutoff005.svg", plot=SystemsByPhylogroupByLocationWeights,
#        width=22, height = 30, units="cm",dpi=300)

#######
#Trying other plot
#I have to turn raw counts into weighted counts

PieChartsPerPhylogroupWide<-merge(DefenseBySystemWithPhyloWithoutCutoff, TipsWeightsDf, by.x="genome",by.y="V1", all.x=T)
PieChartsPerPhylogroupWideWeights<-PieChartsPerPhylogroupWide[,c(2:116)]*PieChartsPerPhylogroupWide[,119]
PieChartsPerPhylogroupWideWeights$genome<-PieChartsPerPhylogroupWide$genome
PieChartsPerPhylogroupWideWeights$id<-PieChartsPerPhylogroupWide$id.x

#write_xlsx(PieChartsPerPhylogroupWideWeights,"Ecoli_7phylogroups_DefenseSystems_counts_weights.xlsx")

PieChartsPerPhylogroupSummary<-PieChartsPerPhylogroupWideWeights[,c(1:115,117)] %>% group_by(id) %>%
  summarise(across(everything(), sum))
PieChartsPerPhylogroupSummary$sum<-rowSums(PieChartsPerPhylogroupSummary[,c(2:116)])
PieChartsPerPhylogroupSummaryLong<-gather(PieChartsPerPhylogroupSummary, 'System','Weight',2:116)
PieChartsPerPhylogroupSummaryLongFilt<-subset(PieChartsPerPhylogroupSummaryLong,
                                              PieChartsPerPhylogroupSummaryLong$Weight/PieChartsPerPhylogroupSummaryLong$sum > 0.03)

Accounted<-PieChartsPerPhylogroupSummaryLongFilt %>% group_by(id) %>% summarize(AccountedSum=sum(Weight))
NotAccounted<-merge(Accounted,PieChartsPerPhylogroupSummary[,c(1,117)], by="id")
colnames(NotAccounted)[3]<-"AllSum"
NotAccounted$System<-rep("other",length(NotAccounted$id))
NotAccounted$sum<-NotAccounted$AllSum - NotAccounted$AccountedSum

PieChartForPlot<-rbind(PieChartsPerPhylogroupSummaryLongFilt[,c("id","System","sum")],
                       NotAccounted[,c("id","System","sum")])

#write_xlsx(PieChartForPlot,"Ecoli_7phylogroups_DefenseSystems_weights_filtered_forplot0.03.xlsx")

palette<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#02818a",
           "#e31a1c","#fdbf6f","#045a8d","#969696","#cab2d6",
           "#6a3d9a",
           "#fa9fb5","#f768a1","#dd3497","#ae017e","#7a0177",
           "#f6e8c3","#4393c3","#66c2a5")


#Plot itself
PieChartsPerPhylogroup<-ggplot(PieChartForPlot,aes(x=factor(1),y=sum, fill=System))+
  geom_bar(stat="identity", position=position_fill(vjust=1))+
  coord_polar("y")+
  scale_fill_manual(values=palette)+
  facet_wrap(~id,nrow=1)+
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(legend.position = "top",
        axis.text = element_blank())
PieChartsPerPhylogroup

ggsave("../figures/Ecoli_phylogroups/Ecoli_7phylogroups_DefenseSystems_weights_filtered_0.03.png",
       plot = PieChartsPerPhylogroup,
       height =10, width = 25, dpi=300,
       units="cm")
ggsave("../figures/Ecoli_phylogroups/Ecoli_7phylogroups_DefenseSystems_weights_filtered_0.03.svg",
       plot = PieChartsPerPhylogroup,
       height =10, width = 25, dpi=300,
       units="cm")
