library(readxl)
library(ggplot2)
#remotes::install_github("coolbutuseless/ggpattern")
#library(ggpattern)
library(ggpubr)
library(reshape2)
library(dplyr)
library(tidyr)
library(stringr)

setwd("../")

folderForResults<-"./figures/EOP_and_AUC"

if (!dir.exists(folderForResults)){
  dir.create(folderForResults)
}else{
  print("dir exists")
}

rangelimits<-c(-5,5)

###EOP data for 4 pairs with three replicates
EOPdata<-read_xlsx("./data/EOP_20230531.xlsx")

EOPPlaqueFormation<-read_xlsx("./data/EOP_NoPlaqueFormation_20230531.xlsx")
EOPPlaqueFormationFull<-EOPPlaqueFormation[rep(seq_len(nrow(EOPPlaqueFormation)), each = 3), ]
EOPPlaqueFormationFull$Replicate<-rep(c(1:3),length(EOPPlaqueFormation$`Gabija+YFP`))

#Plot raw data for each pair
ForRawDataPlot<- gather(EOPdata, System, measurement, 3:12, factor_key=TRUE)

EOPRawPlot<-ggplot(data=ForRawDataPlot,aes(x=Phage, y=measurement))+
  geom_bar(stat="summary") +
  geom_jitter(shape=21)+
  facet_wrap(~System,ncol=2)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
EOPRawPlot

##Calculate the effect
getEffects<-function(dfin)
{
  df<-as.data.frame(dfin)
  #add small numbers so I won't deal with zeros
  df$adjs1<-ifelse(df[,3] == 0, df[,3]+0.000000001, df[,3])
  df$adjs2<-ifelse(df[,4] == 0, df[,4]+0.000000001, df[,4])
  df$AdditivExp<-as.vector(df$adjs1*df$adjs2)
  df$Effect<-ifelse(df$AdditivExp > df[,5],"synergistic",
                       ifelse(df$AdditivExp < df[,5], "antagonistic", "additive"))
  #I normalize effect size here by expected additive effect size
  df$EffectSizePre<-log10(df[,5]/df$AdditivExp)
  #dealing with infinite cases
  df$EffectSize<-ifelse(is.infinite(df$EffectSizePre),
                        ifelse(df$EffectSizePre > 0, rangelimits[2],
                               rangelimits[1]),df$EffectSizePre)
  return(df)
}

###This is to make it similar to what Yi is plotting normally
mergeReplicates<-function(dfeffects,systems){
  #dfeffects<-GabijaTmnEffects
  SDcutoff<-3
  tmp<-dfeffects %>% group_by(Phage)%>%
    summarise(AdditivExpMean= mean(AdditivExp),
                  AdditivExpSD=sd(AdditivExp),
                  EffectSizeMean=mean(EffectSize),
                  EffectSizeSD=sd(EffectSize))
  effectsizeout<-tmp %>% group_by(Phage) %>%
    #I select the largest SD because I expect SD to be the same for experiment and control
    #I use 3 SD cutoff, because I expect error to be normally distributed
    #and 3 SD cover 99.7% of distribution
    mutate(AdditivExpUp=AdditivExpMean+SDcutoff*max(AdditivExpSD,EffectSizeSD),
              AdditivExpDo= AdditivExpMean-SDcutoff*max(AdditivExpSD,EffectSizeSD)) %>%
    mutate(EffectSign=ifelse(EffectSizeMean< AdditivExpDo, "synergistic",
                             ifelse(EffectSizeMean > AdditivExpUp, "antagonistic",
                                    "additive")))
  effectsizeout$pair<-systems
  return(effectsizeout)
}

#Gabija +tmn
EOPGT<-EOPdata[,c(1:5)]
GabijaTmnEffects<-getEffects(EOPGT)
GabijaTmnEffectsMean<-mergeReplicates(GabijaTmnEffects,"Gabija + tmn")


#Druantia III + Zorya II
EOPDZ<-EOPdata[,c(1,2,6:8)]
DruantiaZoryaEffects<-getEffects(EOPDZ)
DruantiaZoryaEffectsMean<-mergeReplicates(DruantiaZoryaEffects,"Druantia III + Zorya II")

#ietAS + ZoryaII
EOPIZ<-EOPdata[,c(1,2,9,7,12)]
ietASZoryaEffects<-getEffects(EOPIZ)
ietASZoryaEffectsMean<-mergeReplicates(ietASZoryaEffects,"ietAS + Zorya II")

#ietAS + Kiwa

EOPIK<-EOPdata[,c(1,2,9,10,11)]
ietASKiwaEffects<-getEffects(EOPIK)
ietASKiwaEffectsMean<-mergeReplicates(ietASKiwaEffects,"ietAS + Kiwa")

#####Plot one value per phage
#####I am going here with means
MeansForPlot<-do.call("rbind", list(GabijaTmnEffectsMean,DruantiaZoryaEffectsMean,
                                    ietASZoryaEffectsMean,ietASKiwaEffectsMean))
#plot
PlotAllPairsMean<-ggplot()+
  geom_tile(data =MeansForPlot,
            aes(y = pair, x = Phage, fill = EffectSizeMean,
                colour = EffectSign),
            width =0.9,
            height=0.9,
            linewidth=1.1)+
  scale_fill_gradient2(low="#1a9850", mid="white",
                       high="#d73027",
                       limits = rangelimits, name = "log10(Effect size)")+
  scale_color_manual(values = c("additive"="#d1e5f0",
                                "antagonistic"="#d73027",
                                "synergistic"="#1a9850"))+
  ylab("")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust =0.5),
        legend.position = "bottom")
PlotAllPairsMean
#save
ggsave("GT_DZ_IZ_IK_EOPEffects_ReplicatesMerged_3SD.png",
       plot=PlotAllPairsMean,
       path=folderForResults,
       width = 40,height=10, dpi=300,
       units = "cm")
ggsave("GT_DZ_IZ_IK_EOPEffects_ReplicatesMerged_3SD.svg",
       plot=PlotAllPairsMean,
       path=folderForResults,
       width = 40,height=10, dpi=300,
       units = "cm")

##########################################################
#####Plot individual plots for systems with all replicates
#####It is the initial plot that I was doing
##Plot effect size
DrawEffectSizePlot<-function(df,colplaque,title)
{
  # df<-GabijaTmnEffects
  # colplaque<-"Gabija+tmn"
  # title<-"Gabija + tmn"
  dfforplot<-merge(df,EOPPlaqueFormationFull[,c("Phage", "Replicate",colplaque)],
                   by=c("Phage","Replicate"))
  colnames(dfforplot)[12]<-"PlaqueFormation"
  
  Plot<-ggplot()+
    geom_tile(data =dfforplot,
              aes(y = Replicate, x = Phage, fill = EffectSize,
                  colour = Effect),
              width =0.9,
              height=0.9,
              size=1.1)+
    scale_fill_gradient2(low="#1a9850", mid="white",
                         high="#d73027",
                         limits = rangelimits, name = "log10(Effect size)")+
    scale_color_manual(values = c("additive"="#d1e5f0",
                                  "antagonistic"="#d73027",
                                  "synergistic"="#1a9850"))+
    ggtitle(title)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust =0.5),
          legend.position = "bottom")
  
  ######Patterns on top of tiles to show 
  ######the ones that are not forming plaques
  #   Plot<-ggplot()+
  #   geom_tile_pattern(data =dfforplot,
  #             aes(y = Replicate, x = Phage, fill = EffectSize, 
  #                 colour = Effect,
  #                 pattern = as.factor(PlaqueFormation)),
  #             width =0.9,
  #             height=0.9,
  #             size=1.1,
  #             pattern_fill = "#bdbdbd",
  #             pattern_color = NA,
  #             pattern_angle = 45,
  #             pattern_density = 0.1,
  #             pattern_spacing = 0.025,
  #             pattern_key_scale_factor = 0.4)+
  #   scale_pattern_manual(values = c("none","stripe"),
  #                        name="No plaques")+
  #   scale_fill_gradient2(low="#1a9850", mid="white",
  #                        high="#d73027",
  #                        limits = rangelimits, name = "log10(Effect size)")+
  #   scale_color_manual(values = c("additive"="#d1e5f0",
  #                                 "antagonistic"="#d73027",
  #                                 "synergistic"="#1a9850"))+
  #   ggtitle(title)+
  #   theme_classic()+
  #   theme(axis.text.x = element_text(angle = 90, hjust=1, vjust =0.5),
  #         legend.position = "bottom")
  
  return(Plot)
}

DrawEffectDensity<-function(df,title)
{
  denplot<-ggplot(data=df)+
    geom_density(aes(x=EffectSize,
                     fill=as.factor(Replicate)),
                 alpha =0.2)+
    labs(fill = "Replicate", x = "log10(Effect Size)")+
    xlim(c(-5,5))+
    ggtitle(title)+
    theme_classic()
  return(denplot)
}

GTdensity<-DrawEffectDensity(GabijaTmnEffects,"Gabija + tmn")
DZdensity<-DrawEffectDensity(DruantiaZoryaEffects,"Druantia III + Zorya II")
IZdensity<-DrawEffectDensity(ietASZoryaEffects,"ietAS + Zorya II")
IKdensity<-DrawEffectDensity(ietASKiwaEffects,"ietAS + Kiwa")

GabijaTmn<-DrawEffectSizePlot(GabijaTmnEffects,"Gabija+tmn","Gabija + tmn")
DruantiaZorya<-DrawEffectSizePlot(DruantiaZoryaEffects,"DruantiaIII+ZoryaII","Druantia III + Zorya II")
ietASZorya<-DrawEffectSizePlot(ietASZoryaEffects,"ZoryaII+ietAS","ietAS + Zorya II")
ietASKiwa<-DrawEffectSizePlot(ietASKiwaEffects,"ietAS+Kiwa","ietAS + Kiwa")

ThreeEffectPlotsCombined<-ggarrange(GabijaTmn,
                                    DruantiaZorya,
                                    ietASZorya,
                                    ietASKiwa,
                                    ncol = 1,
                                    common.legend = T)+
  bgcolor("white")
#ThreeEffectPlotsCombined

ggsave("GT_DZ_IZ_IK_EOPEffects.png",
       plot=ThreeEffectPlotsCombined,
       path=folderForResults,
       width = 40,height=25, dpi=300,
       units = "cm")

ThreeDensitiesCombined<-ggarrange(GTdensity,
                                  DZdensity,
                                  IZdensity,
                                  IKdensity,
                                  ncol=1,
                                  common.legend = T)+
  bgcolor("white")

ggsave("GT_DZ_IZ_IK_EOPEffectDensity.png",
       plot=ThreeDensitiesCombined,
       path=folderForResults,
       width=15, height=25, dpi=300,
       units="cm")

#Combine in one plot
BothCombined<-ggarrange(ThreeEffectPlotsCombined,
          ThreeDensitiesCombined,
          ncol=2,
          widths = c(1,0.3))

ggsave("GT_DZ_IZ_IK_EOPEffect_Density_combined.png",
       plot=BothCombined,
       path=folderForResults,
       width = 55,height=25, dpi=300,
       units = "cm")
ggsave("GT_DZ_IZ_IK_EOPEffect_Density_combined.svg",
       plot=BothCombined,
       path=folderForResults,
       width = 55,height=25, dpi=300,
       units = "cm")


#################################################################################
#################################################################################
##Calculate AUC
##deal with raw data

LiquidEssayRawfile<-"./data/liquid_assay_4pairs.xlsx"
LiquidEssayrawSheets <- readxl::excel_sheets(LiquidEssayRawfile)
LiquidEssayRawtibble <- lapply(LiquidEssayrawSheets, function(x) readxl::read_excel(LiquidEssayRawfile, sheet = x,
                                                                    col_types = "numeric"))
names(LiquidEssayRawtibble) <- LiquidEssayrawSheets

maxmin<-200

GetDfFromRaw<-function(System,lim=maxmin){
  #System<-"Z"
  Gdf<-as.data.frame(LiquidEssayRawtibble[[System]])
  Gdf<-subset(Gdf,!is.na(Gdf$`Time/min`))
  GdfTrans<-as.data.frame(t(Gdf))[c(2:ncol(Gdf)),]
  ###remove NA rows
  
  ExpName<-rownames(GdfTrans)[c(TRUE,rep(FALSE, 2)) ]
  Time<-t(Gdf)[1,]
  colnames(GdfTrans)<-Time
  GdfTrans<-GdfTrans[,c(1:(lim/5+1))]#this is to fix difference in experiment length
  GdfTrans$replicate<-rep(c(1,2,3),length(GdfTrans$`5`)/3)
  GdfTrans$experimentLN<-rep(ExpName, each=3)
  GdfTransFull<-GdfTrans %>% separate(experimentLN,
                       into = c("Phage","MOI"),
                       sep=paste("-",System,"-MOI",sep=""))
  GdfTransFull$system<-rep(System,length(GdfTransFull$`5`))
  return(GdfTransFull)
}

###DZ
DRaw<-GetDfFromRaw("D")
ZRaw<-GetDfFromRaw("Z")
DZRaw<-GetDfFromRaw("DZ")

DZfull<-do.call("rbind", list(DRaw,ZRaw,DZRaw))

###IZ
IRaw<-GetDfFromRaw("I")
ZRaw<-GetDfFromRaw("Z")
IZRaw<-GetDfFromRaw("ZI")
IZRaw$system<-rep("IZ",length(IZRaw$`0`))

IZfull<-do.call("rbind", list(IRaw,ZRaw,IZRaw))

###GT
GRaw<-GetDfFromRaw("G")
TRaw<-GetDfFromRaw("T")
GTRaw<-GetDfFromRaw("GT")

GTfull<-do.call("rbind", list(GRaw,TRaw,GTRaw))

###IK
KRaw<-GetDfFromRaw("K")
IKRaw<-GetDfFromRaw("IK")

IKfull<-do.call("rbind", list(IRaw,KRaw,IKRaw))

######Read tab with controls
##
ControlsDf<-as.data.frame(LiquidEssayRawtibble[["Controls"]])
ControlsDf<-subset(ControlsDf,!is.na(ControlsDf$`Time/min`))
#Fix experiment length
ControlsDfTrans<-as.data.frame(t(ControlsDf))[c(2:ncol(ControlsDf)),]
ControlsName<-rownames(ControlsDfTrans)[c(TRUE,rep(FALSE, 2)) ]
ControlsNameSh<-str_replace(ControlsName,"-Control","")
ControlsNameSh<-str_replace(ControlsNameSh,"ZI","IZ")
###Fix to use I control for IK
ControlsNameSh<-recode(ControlsNameSh,I="IK")
ControlsTime<-t(ControlsDf)[1,]
colnames(ControlsDfTrans)<-ControlsTime
ControlsDfTrans$replicate<-rep(c(1,2,3),length(ControlsDfTrans$`5`)/3)
ControlsDfTrans$system<-rep(ControlsNameSh, each=3)
ControlsLong<-gather(ControlsDfTrans,Time,OD,`0`:`400`, factor_key = T)
ControlsLong$minutes<-strtoi(ControlsLong$Time)


################################
######Draw raw data on the plots

####filter

filterFullDf<-function(Df)
{
  #Retain only phages for which we have observations in all three experiments
  CheckPhages<-unique(Df[,c("Phage","system")]) %>% group_by(Phage)%>% summarise(count=n())
  PhagesToRetain<-CheckPhages[CheckPhages$count == 3,]$Phage
  DfadjPhage<-subset(Df, Df$Phage %in% PhagesToRetain)
  #Only keep MOI that is measured for all
  CheckMOI<-unique(DfadjPhage[,c("Phage","MOI")]) %>% group_by(MOI) %>%
    summarise(count=n())
  MOIToRetain<-CheckMOI[CheckMOI$count == length(PhagesToRetain),]$MOI
  DfadjPhageMOI<-subset(DfadjPhage,
                        DfadjPhage$MOI %in% MOIToRetain)
  return(DfadjPhageMOI)
}
####plot
getPhagePlot<-function(Df,System,mins,adj=T)
{
  # Df<-DZfull
  # System<-"DZ"
  # mins<-maxmin
  
  DfadjPhageMOI<-Df
    if(adj)
  {
    DfadjPhageMOI<-filterFullDf(Df)
  }
  ControlsSys<-subset(ControlsLong, ControlsLong$system == System)
  ControlsSys<-subset(ControlsSys, as.vector(ControlsSys$Time) <= mins)
  DfphageLong<-gather(DfadjPhageMOI,Time,OD,`0`:`200`, factor_key = T)
  DfphageLong$minutes<-strtoi(DfphageLong$Time)
  
  Phageplot<-ggplot()+
    geom_line(data=DfphageLong,
              aes(x=minutes,y=OD,color=system, linetype=as.factor(replicate)),
              linewidth=1.1)+
    scale_color_manual(values=c("#d95f02","#1b9e77","#7570b3"))+
    geom_line(data=ControlsSys,
              aes(x=minutes,y=OD,linetype=as.factor(replicate)),
              color = "grey",
              linewidth=0.9)+
    ylim(c(0,1.5))+
    xlim(c(0,mins))+
    facet_grid(Phage~MOI)+
    #ggtitle(PhageName)+
    theme_classic()+
    guides(linetype=guide_legend(title="replicate"))+
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
  return(Phageplot)
}

####data

DZPlot<-getPhagePlot(DZfull,"DZ",maxmin)
IZPlot<-getPhagePlot(IZfull,"IZ",maxmin)
GTPlot<-getPhagePlot(GTfull,"GT",maxmin)
IKPlot<-getPhagePlot(IKfull,"IK",maxmin)
#save
ggsave(paste("LiquidAssay_DZ_raw_",maxmin,".png", sep=""),
       plot=DZPlot,
       path=folderForResults,
       height=30, width=30, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_DZ_raw_",maxmin,".svg", sep=""),
       plot=DZPlot,
       path=folderForResults,
       height=30, width=30, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_IZ_raw_",maxmin,".png", sep=""),
       plot=IZPlot,
       path=folderForResults,
       height=30, width=30, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_IZ_raw_",maxmin,".svg", sep=""),
       plot=IZPlot,
       path=folderForResults,
       height=30, width=30, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_GT_raw_",maxmin,".png", sep=""),
       plot=GTPlot,
       path=folderForResults,
       height=30, width=30, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_GT_raw_",maxmin,".svg", sep=""),
       plot=GTPlot,
       path=folderForResults,
       height=30, width=30, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_IK_raw_",maxmin,".png", sep=""),
       plot=IKPlot,
       path=folderForResults,
       height=30, width=30, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_IK_raw_",maxmin,".svg", sep=""),
       plot=IKPlot,
       path=folderForResults,
       height=30, width=30, dpi=300, units ="cm")

###################
####Transform to AUC
####This is only necessary to lock at the raw data before switching to means
# getAUCplot<-function(Df,System,maxobs=maxmin)
# {
#   # Df<-IZfull
#   # System<-"IZ"
#   DZadj<-filterFullDf(Df)
#   DZadjGrowth<-DZadj[,c(2:(ncol(DZadj)-4))]-DZadj[,c(1)]
#   DZadjGrowthMin<-DZadjGrowth[,c(1:which(colnames(DZadjGrowth)==toString(maxobs)))]
#   DZadjNoNeg<-DZadjGrowthMin
#   DZadjNoNeg[DZadjNoNeg<0]<-0
#   DZadjNoNegMult<-5*DZadjNoNeg
#   DZadjNoNegMult$AUC<-rowSums(DZadjNoNegMult)
#   
#   DZAUCDf<-DZadj[,c((ncol(DZadj)-3):(ncol(DZadj)))]
#   DZAUCDf$AUC<-DZadjNoNegMult$AUC
#   
#   DZAUCwide<-spread(DZAUCDf, system, AUC)
#   DZAUCwide$AdditivExp<-DZAUCwide[,4]+DZAUCwide[,6]
#   
#   DZAUCwide$logMOI<-log10(as.double(DZAUCwide$MOI))
#   ForPlot<-gather(DZAUCwide,Line,AUC,4:7, factor_key = T)
#   
#   ControlsSys<-subset(ControlsLong, ControlsLong$system == System)
#   ControlSysLim<-subset(ControlsSys,as.numeric(as.character(ControlsSys$Time)) <= maxmin)
#   ControlSysLimWide<-spread(ControlSysLim[,c(1:4)],Time,OD)
#   ControlSysLimWideGrowth<-ControlSysLimWide[,c(4:ncol(ControlSysLimWide))]-
#     ControlSysLimWide[,c(3)]
#   ControlSysLimWideGrowth[ControlSysLimWideGrowth<0]<-0
#   ControlMult<-5*ControlSysLimWideGrowth
#   ControlMult$AUC<-rowSums(ControlMult)
#   
#   
#   legendcolors<-c("#d95f02","#1b9e77","#7570b3","black")
#   names(legendcolors)<-c(colnames(DZAUCwide)[4:6],"AdditivExp")
#   
#   AUCPlot<-ggplot()+
#     geom_line(data=ForPlot,aes(x=logMOI, y=AUC,
#                   color=Line,
#                   linetype=as.factor(replicate)))+
#     facet_wrap(~Phage, nrow = 1)+
#     geom_hline(yintercept = ControlMult$AUC,
#                color="grey")+#,linetype=ControlSysLimWide$replicate)+
#     ylab("AUC")+
#     xlab("MOI(log10)")+
#     theme_classic()+
#     theme(strip.background = element_blank(),
#           strip.text = element_text(face = "bold"))+
#     scale_color_manual(values=legendcolors)
#   
#   return(AUCPlot)
# }
# 
# DZAUCPlot<-getAUCplot(DZfull,"DZ")
# IZAUCPlot<-getAUCplot(IZfull,"IZ")
# GTAUCPlot<-getAUCplot(GTfull,"GT")
# 
# ggsave(paste("LiquidAssay_DZ_AUC_",minfrst,".png", sep=""),
#        plot=DZAUCPlot,
#        path=folderForResults,
#        height=8, width=35, dpi=300, units ="cm")
# ggsave(paste("LiquidAssay_IZ_AUC_",minfrst,".png", sep=""),
#        plot=IZAUCPlot,
#        path=folderForResults,
#        height=8, width=35, dpi=300, units ="cm")
# ggsave(paste("LiquidAssay_GT_AUC_",minfrst,".png", sep=""),
#        plot=GTAUCPlot,
#        path=folderForResults,
#        height=8, width=35, dpi=300, units ="cm")

###get AUC plot mean
getAUCMeanplot<-function(Df,System,maxobs=maxmin,ci)
{
  DZadj<-filterFullDf(Df)
  DZadjGrowth<-DZadj[,c(2:(ncol(DZadj)-4))]-DZadj[,c(1)]
  DZadjGrowthMin<-DZadjGrowth[,c(1:which(colnames(DZadjGrowth)==toString(maxobs)))]
  DZadjNoNeg<-DZadjGrowthMin
  DZadjNoNeg[DZadjNoNeg<0]<-0
  DZadjNoNegMult<-5*DZadjNoNeg
  DZadjNoNegMult$AUC<-rowSums(DZadjNoNegMult)
  
  DZAUCDf<-DZadj[,c((ncol(DZadj)-3):(ncol(DZadj)))]
  DZAUCDf$AUC<-DZadjNoNegMult$AUC
  
  #Calculate Expected sum
  DZAUCDfwide<-spread(DZAUCDf,system,AUC)
  DZAUCDfwide$AdditivExp<-DZAUCDfwide[,4]+DZAUCDfwide[,6]
  AdditivExpMeanSD<-DZAUCDfwide %>% group_by(Phage,MOI) %>%
    summarise(Mean=mean(AdditivExp),
              SD=sd(AdditivExp))
  AdditivExpMeanSD$system<-rep("Expected", length(AdditivExpMeanSD$Phage))
  
  ###function to calculate confidence intervals
  GetConfidenceInterval<-function(Df,r=3,confint=ci)
  {
    Df$error<-qt(confint,df=(r-1))*Df$SD/sqrt(r)
    Df$upper<-Df$Mean+Df$error
    Df$lower<-Df$Mean-Df$error
    return(Df)
  }
  ###
  AdditivExpForPlot<-GetConfidenceInterval(AdditivExpMeanSD)
  
  ##Calculate mean over replicates
  DZAUCDfMean<-DZAUCDf %>% group_by(Phage,MOI,system)%>%
    summarise(Mean=mean(AUC),
              SD=sd(AUC))
  DZAUCForPlot<-GetConfidenceInterval(DZAUCDfMean)
  
  #Actually now I have all the controls for induvidual systems and pairs of systems
  #But I am only using the ones for pairs because they all are similar
  #I need supplementary figure comparing controls
  ControlsSys<-subset(ControlsLong, ControlsLong$system == System)
  ControlSysLim<-subset(ControlsSys,as.numeric(as.character(ControlsSys$Time)) <= maxobs)
  ControlSysLimWide<-spread(ControlSysLim[,c(1:4)],Time,OD)
  ControlSysLimWideGrowth<-ControlSysLimWide[,c(4:ncol(ControlSysLimWide))]-
    ControlSysLimWide[,c(3)]
  ControlSysLimWideGrowth[ControlSysLimWideGrowth<0]<-0
  ControlMult<-5*ControlSysLimWideGrowth
  ControlMult$AUC<-rowSums(ControlMult)
  ControlMean<-mean(ControlMult$AUC)
  ControlSD<-sd(ControlMult$AUC)
  PreplotControlDf<-AdditivExpForPlot[,c("Phage","MOI")]
  PreplotControlDf$Mean<-rep(ControlMean,length(AdditivExpForPlot$Phage))
  PreplotControlDf$SD=rep(ControlSD,length(AdditivExpForPlot$Phage))
  PreplotControlDf$system<-rep("Control",length(AdditivExpForPlot$Phage))
  ControlDfForPlot<-GetConfidenceInterval(PreplotControlDf)
  
  
  ForPlot<-rbind(AdditivExpForPlot,ControlDfForPlot,DZAUCForPlot)
  ForPlot$logMOI<-log10(as.double(ForPlot$MOI))
  
  ##Plot results
  #Select colors
  legendcolors<-c("#d95f02","#1b9e77","#7570b3","black","grey")
  names(legendcolors)<-c(colnames(DZAUCDfwide)[4:6],"Expected","Control")
  #Plot
  AUCPlotWithConf<-ggplot(data=ForPlot)+
    geom_ribbon(aes(x=logMOI,ymin=lower,ymax=upper, fill=system,
                    color=system),
                linewidth=0.1,
                alpha=0.2)+
    geom_line(aes(x=logMOI, y=Mean,
                               color=system),
              linewidth=0.3,
              linetype="dotted")+
    geom_point(aes(x=logMOI, y=Mean,
                  color=system),
              shape=16)+
    facet_wrap(~Phage, nrow = 1)+
    ylab("AUC")+
    xlab("MOI(log10)")+
    theme_classic()+
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"))+
    scale_color_manual(values=legendcolors)+
    scale_fill_manual(values=legendcolors)
  
  return(AUCPlotWithConf)
}

confidenceinterval<-.95
DZAUCMeanPlot<-getAUCMeanplot(DZfull,"DZ",ci=confidenceinterval)
IZAUCMeanPlot<-getAUCMeanplot(IZfull,"IZ",ci=confidenceinterval)
GTAUCMeanPlot<-getAUCMeanplot(GTfull,"GT",ci=confidenceinterval)
IKAUCMeanPlot<-getAUCMeanplot(IKfull,"IK",ci=confidenceinterval)

ggsave(paste("LiquidAssay_DZ_AUCMean_",maxmin,"_ci_",confidenceinterval,".svg", sep=""),
       plot=DZAUCMeanPlot,
       path=folderForResults,
       height=8, width=35, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_IZ_AUCMean_",maxmin,"_ci_",confidenceinterval,".svg", sep=""),
       plot=IZAUCMeanPlot,
       path=folderForResults,
       height=8, width=35, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_GT_AUCMean_",maxmin,"_ci_",confidenceinterval,".svg", sep=""),
       plot=GTAUCMeanPlot,
       path=folderForResults,
       height=8, width=35, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_IK_AUCMean_",maxmin,"_ci_",confidenceinterval,".svg", sep=""),
       plot=IKAUCMeanPlot,
       path=folderForResults,
       height=8, width=35, dpi=300, units ="cm")

#
ggsave(paste("LiquidAssay_DZ_AUCMean_",maxmin,"_ci_",confidenceinterval,".png", sep=""),
       plot=DZAUCMeanPlot,
       path=folderForResults,
       height=8, width=35, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_IZ_AUCMean_",maxmin,"_ci_",confidenceinterval,".png", sep=""),
       plot=IZAUCMeanPlot,
       path=folderForResults,
       height=8, width=35, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_GT_AUCMean_",maxmin,"_ci_",confidenceinterval,".png", sep=""),
       plot=GTAUCMeanPlot,
       path=folderForResults,
       height=8, width=35, dpi=300, units ="cm")
ggsave(paste("LiquidAssay_IK_AUCMean_",maxmin,"_ci_",confidenceinterval,".png", sep=""),
       plot=IKAUCMeanPlot,
       path=folderForResults,
       height=8, width=35, dpi=300, units ="cm")