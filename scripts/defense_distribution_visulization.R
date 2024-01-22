library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
library(readr)
library(stringr)

path<-getwd()
setwd(path)

##Figure 1B-the frequency of the amount of defense system in strains by phylogroup--------------------------
data = read_csv("/Users/wubaobei/Desktop/test/defence/results/26k/latest/phylogroup_vs_number_of_defense.csv")

result <- data |>
  group_by(phylogroup) |>
  dplyr::summarize(weighted_median = weighted.median(def_num, frequency)+0.5)


p <- ggplot(data, aes(x = def_num, y = frequency)) +
  geom_bar(stat = 'identity') +
  facet_wrap(phylogroup ~ ., scales = "free_x",ncol=2) +
  geom_vline(data = result, aes(xintercept = weighted_median), color = "red",  size = 0.5) +
  scale_color_manual(values = "red")+
  theme(axis.text.x = element_text(family = "ArialMT",size=6),
        axis.text.y = element_text(family = "ArialMT",size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #legend.text = element_text(family = "ArialMT",size=6),
        #title = element_text(family = "ArialMT",size=8),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color="black",linewidth=0.4),
        panel.border = element_rect(fill=NA,colour = "black",linewidth=0.5,linetype = "solid"),
        panel.background = element_blank(),
        legend.position="none",
        strip.text.y = element_text(family = "ArialMT",size=6))+
  coord_cartesian(xlim = c(0, 20))
p
ggsave(plot=p,filename = "../figures/Ecoli_defense_distribution/F1B_ecoli_defense_amount_frequency.pdf",
       width=7,height = 11,units = "cm",device = "pdf")



##Figure 1C-defense system percentage by phylogroup--------------------------
data = read_csv("../data/phylogroup_defense_raw_percentage.csv")
data_long = data |> pivot_longer(cols=c("unassigned","E2","E1","C","B22","B21","B1","A"),
                                 names_to="phylogroup",
                                 values_to="percentage")
data1 = read_excel("../data/Table_S1_defence_occurrence_in_datasets.xlsx",sheet = "E. coli")
colnames(data1)[colnames(data1) == "immune_system"] = "defense_system"
defense_class = distinct(data1[c("defense_system","class")])
#count the frequency of defense system
frequency_table = table(data1$defense_system)
sorted_frequency = sort(frequency_table,decreasing = TRUE)
sorted_df = data.frame(Value=names(sorted_frequency),Frequency=sorted_frequency)
sorted_df$Value = factor(sorted_df$Value,levels=sorted_df$Value)
colnames(sorted_df) = c("defense_system","defense2","frequency")
sorted_df = merge(sorted_df,defense_class,by="defense_system",all.x=TRUE)
sorted_df = sorted_df |> arrange(desc(frequency))
#draw
#--barplot for the defense system frequency
#--also used for plotting the defense system frequency for 4 orders
p = ggplot(sorted_df, aes(x=defense_system,y=frequency,fill=class))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,family = "ArialMT",size=6),
        axis.text.y = element_text(family = "ArialMT",size=6),
        legend.position="none")+
  scale_fill_manual(values=c("Abi"="#9112b4","RM"="#94b447","RM-like"="#94b447",
                             "Phage nucleic acid cleavage"="#0072b5",
                             "Unknown"="#333333","TA"="#ef87be",
                             "Inhibition of phage transcription"="#6f99ad"))+
  geom_text(aes(x=defense_system,y=4e4,label=frequency),
            angle=90,vjust = -0.5,
            family = "ArialMT",size=6)
p
ggsave(plot=p,filename = "../figures/Ecoli_defense_distribution/F1C_ecoli_defense_frequency.pdf",
       width=21,height = 4,units = "cm",device = "pdf")
#--facet percentage plot (based on phylogroup)
data_long = merge(data_long,defense_class,by="defense_system",all.x=TRUE)
p = ggplot(data_long, aes(x=defense_system,y=percentage,fill=class))+
  geom_bar(stat="identity")+
  facet_grid(phylogroup~.)+
  scale_fill_manual(values=c("Abi"="#9112b4","RM"="#94b447","RM-like"="#94b447",
                             "Phage nucleic acid cleavage"="#0072b5",
                             "Unknown"="#333333","TA"="#ef87be",
                             "Inhibition of phage transcription"="#6f99ad"))+
  scale_x_discrete(limits=sorted_df$defense_system)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,family = "ArialMT",size=6),
        axis.text.y = element_text(family = "ArialMT",size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #legend.text = element_text(family = "ArialMT",size=6),
        #title = element_text(family = "ArialMT",size=8),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color="black",linewidth=0.4),
        panel.border = element_rect(fill=NA,colour = "black",linewidth=0.5,linetype = "solid"),
        panel.background = element_blank(),
        legend.position="none",
        strip.text.y = element_text(family = "ArialMT",size=6))
p
ggsave(plot=p,filename = "../figures/Ecoli_defense_distribution/F1C_phylogroup_defense_percentage_raw.pdf",
       width=21,height = 14.5,units = "cm",device = "pdf")



##Figure 2C and 2D-histogram of defense system pairs--------------------------
##positive pairs-----
data = read_excel("./data/Distances_between_positive_significant_pairs.xlsx")
p = ggplot(data,aes(x=distance))+
  geom_histogram(bins = 500)+
  scale_x_continuous(breaks=seq(0,3000000,by=500000))+
  geom_vline(xintercept = median(data$distance),lwd=0.5,color="red",alpha=0.5)+
  theme(panel.border = element_rect(fill=NA,colour = "black",size=1,linetype = "solid"),
        panel.background = element_blank(),
        axis.text.x = element_text(family = "ArialMT",size=8),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,family = "ArialMT",size=6),
        axis.text.y = element_text(family = "ArialMT",size=8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color="black",linewidth=0.4))
p
ggsave(plot=p,filename = "./figures/Ecoli_defense_distribution/F2D_positive_distance_histogram.png",
       width=10,height = 3,units = "cm",device = "png")
##all------
data = read_excel("./data/Distances_to_closest_system_control_All.xlsx")
p = ggplot(data,aes(x=mindist))+
  geom_histogram(bins = 500)+
  scale_x_continuous(breaks=seq(0,4000000,by=500000))+
  geom_vline(xintercept = median(data$mindist),lwd=0.5,color="red")+
  theme(panel.border = element_rect(fill=NA,colour = "black",size=1,linetype = "solid"),
        panel.background = element_blank(),
        axis.text.x = element_text(family = "ArialMT",size=8),
        axis.text.y = element_text(family = "ArialMT",size=8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color="black",size=0.4))
p
ggsave(plot=p,filename = "./figures/Ecoli_defense_distribution/F2C_all_distance_histogram.png",
       width=10,height = 3,units = "cm",device = "png")



##FS1B - defense system location and percentage--------------------------
data = read_csv("./data/FS1B_defense_location_percentage.csv")
data1 = read_csv("./data/ecoli_filtered.csv")
#count the frequency of defense system
frequency_table = table(data1$system)
sorted_frequency = sort(frequency_table,decreasing = TRUE)
sorted_df = data.frame(Value=names(sorted_frequency),Frequency=sorted_frequency)
sorted_df$Value = factor(sorted_df$Value,levels=sorted_df$Value)
p = ggplot(data, aes(x=system,y=percentage,fill=location))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("chromosome"="#772E25","prophage in chromosome"="#C44536",
                             "plasmid"="#EDDDD4","prophage in plasmid"="#7f646f"))+
  scale_x_discrete(limits=sorted_df$Value)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,family = "ArialMT",size=6),
        axis.text.y = element_text(family = "ArialMT",size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #legend.text = element_text(family = "ArialMT",size=6),
        #title = element_text(family = "ArialMT",size=8),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color="black",linewidth=0.4),
        panel.border = element_rect(fill=NA,colour = "black",linewidth=0.5,linetype = "solid"),
        panel.background = element_blank(),
        legend.position="none",
        strip.text.y = element_text(family = "ArialMT",size=6))
#legend.text = element_text(family = "ArialMT",size=6))
p
ggsave(plot=p,filename = "./figures/Ecoli_defense_distribution/FS1B_defense_location_percentage_raw.pdf",
       width=19.5,height = 5,units = "cm",device = "pdf")



##Figure S2C-the distance between cooccurring systems--------------------------
data = read_excel("../data/F2D_Distances_between_positive_significant_pairs.xlsx")
p = ggplot(data,aes(x=pair2,y=distance))+
  geom_violin(aes(x=pair2,y=distance),
              trim=TRUE,scale="width",draw_quantiles = c(0.5))+
  #geom_point(size=0.5)+
  theme(panel.border = element_rect(fill=NA,colour = "black",size=0.5,linetype = "solid"),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,family = "ArialMT",size=6),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1,family = "ArialMT",size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color="black",size=0.4))
p
p1 = p+annotate("text",1:length(table(data$pair2)),
                y=3e6,label=table(data$pair2),
                angle=90,size=6/.pt,family="ArialMT")

p1
ggsave(plot=p1,filename = "../figures/FS2B_positive_distance.pdf",
       width=21,height = 7,units = "cm",device = "pdf")



##Figure S2D-the location of cooccurring defense system pairs--------------------------
data = read.csv("../data/FS2C_Distances_between_positive_significant_pairs_location2.csv")
p = ggplot(data, aes(x = combination, y = percentage, fill = location)) +
  geom_bar(stat = "identity") +
  theme(panel.border = element_rect(fill=NA,colour = "black",linewidth=0.5,linetype = "solid"),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,family = "ArialMT",size=8),
        axis.text.y = element_text(angle = 90,family = "ArialMT",size=8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color="black",linewidth=0.4),
        legend.position = "None")+
  scale_fill_manual(values=c("chromosome"="#772E25","prophage-chromosome"="#eab05f",
                             "plasmid"="#EDDDD4","prophage-plasmid"="#FFEBAD",
                             "prophage"="#A4514F","plasmid-chromosome"="#7f646f"))
p
ggsave(plot=p,filename = "../figures/Ecoli_defense_distribution/FS2C_positive_location.pdf",
       width=20.5,height = 7,units = "cm",device = "pdf")



##Figure S4B-venn plot to show 5 datasets defense composition difference-----------------------------------------
library(VennDiagram)

data = read_excel("../data/FS4B_5_dataset_defense_frequency1.xlsx")
sets_list <- list(Ba = na.omit(data$ba), Bu = na.omit(data$bu), ec = na.omit(data$ec), 
                  en = na.omit(data$en), ps = na.omit(data$ps))
mycol = c("#90A5A7","#F0DFA7","#C2B099","#C18265","#AB6355")
venn.plot <- venn.diagram(
  x = sets_list,
  category.names = c("Bacillales", "Burkholderiales", "E. coli", "Enterobacterales", "Pseudomonadales"),
  filename = NULL,  # Use NULL to display the diagram instead of saving to a file
  output = TRUE,
  # Output features
  imagetype="pdf" ,
  height = 210 , 
  width = 210 , 
  #resolution = 300,
  
  # Circles
  #lwd = 2,
  lty = 'blank',
  fill = mycol,
  
  # Numbers
  cex = .6,
  #fontface = "regular",
  fontfamily = "ArialMT",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "ArialMT",
  #rotation = 1
)

ggsave("../figures/FS4B_Ecoli_defense_distribution/5_dataset_def_composition_venn.pdf",
       plot = venn.plot,height = 4, width =4, units ="cm", dpi=300,device="pdf")
grid.draw(venn.plot)