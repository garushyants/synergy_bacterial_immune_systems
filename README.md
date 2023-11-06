### Original code for Wu, Garushyants et al. "Synergistic anti-phage activity of bacterial immune systems"

#### **scripts**
Contain scripts to reproduce the analysis and visualize data

*Ecoli_numContigs_vc_numDefence.R* builds correlation between number of contigs in genome and number of predicted immune systems

Scripts *Ecoli26k_tree_phylogroups_defense_systems.R* and *Visualize_phylogenetic_trees_orders.R* visualize phylognetic trees and add data to them

*Ecoli_defsystems_presence_by_phylogroup.R* allows to analyze immune systems content in E. coli phylogroups and also performs enrichment analysis

*defense_systems_content_vs_phylogenetic_distance.R* builds correlations between immune systems content and phylogenetic distance

*Pagel_calculate_sge_farm.R* is the main script to run the Pagel correlation test for the binary traits

*Pagel_visualization_Ecoli.R* and *Pagel_visualization_orders.R* calculate directionality of interactions for the significant results and visualize the data for 26k E. coli dataset and order level datasets accordingly.

*Pagel_all_datasets_igraph.R* builds graphs for Pagel results

*Pagel_compare_results_all_datasets.R* script that compares and vizualizes results of Pagel test for all five datasets

*Draw_trees_for_selected_pairs.R* visualizes phylogenetic trees for selected pairs of immune systems.

*positive_vs_negative_genome_distance.R* analyses genomic distances for the co-occuring immune systems pairs and compares it with background distributions

*104_plasmids_cactus_visualization.R* visualizes [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) output for 104 plasmids carrying tmn and Gabija systems

*SynergisticVsAntagonistic_4pairs.R* Analyses EOP and liquid assay experiments data for four pairs of immune systems: Zorya II + Druantia III, tmn + Gabija, Kiwa + ietAS, and ietAS + Zorya II

*defense_distribution_visulization.R* Visuslize the data for defense system distribution, including Figure 1B, 1C, S2B, S2C, and S4.

#### **data**
Contains data required to run the scripts above and some of the outputs

*Table_S1_defence_occurrence_in_datasets.xlsx* contains the defense system distribution in *E. coli*, Enterobacteriales, Burkholderiales, Bacillales and Pseudomonadales, which is used for the defense system co-occurrence calculation. This is the combination of *26k_Ecoli_with_prophages.csv*, *merged_baci3.csv*, *merged_burk3.csv*, *merged_enter3.csv* and *merged_pseu3.csv*. [Figure 1A and S4]

*F1B_phylogroup_vs_number_of_defense.csv* provides information on the frequency of the amount of defense system within different *E. coli* strains, categorized by their phylogroups. [Figure 1B]

*F1C_phylogroup_defense_raw_percentage.csv* shows the percentage of the location of defense systems, categorized by their phylogroups. [Figure 1C]

*Table_S3_defense_cooccurrence_in_datasets.xlsx* contains the results of defense system co-occurrence calculation for *E. coli*, Enterobacteriales, Burkholderiales, Bacillales and Pseudomonadales. This is the combination of *Ecoli_pagel_all_results_with_direction.xlsx*, *enter_pagel_all_results_with_direction_0.01.xlsx*, *pseu_pagel_all_results_with_direction_0.01.xlsx*, *baci_pagel_all_results_with_direction_0.01.xlsx* and *burk_pagel_all_results_with_direction_0.01.xlsx*. [Figure 2B and 4]

*F2C_Distances_to_closest_system_control_All.xlsx* calculates the distance between a defense system and its closest systems in the *E.coli* dataset with complete genomes only. [Figure 2C]

*F2D_Distances_between_positive_significant_pairs.xlsx* calculates the distance between a defense system and its positive co-occurreing systems in the *E.coli* dataset with complete genomes only. [Figure 2D and S2B]

*EOP_20230531.xlsx* contains the results of Efficiency of Plating assay. [Figure 3B and S3A]

*F3C_TPI.xlsx* contains the results of time post infection assay. [Figure 3C and S3B]

*liquid_assay_4pairs.xlsx* contains the results of liquid assay. [Figure 3D, S3C and S3D]

*F5_color_plasmid.xlsx* contains the results of color plasmid assay. [Figure 5]

*F6_EOP.xlsx* contains the results of Efficiency of Plating assay investigating the cooperation of tmn and other ATPase domain containing systems. [Figure 6]

*FS1B_26k_location_system.csv* contains the raw count and percentage of defense systems' location. [Figure S1B]

*FS2A_ecoli_significance_pair_frequency.csv* quantifies the co-occurring pattern between different systems, highlighting their co-occurrence and significance. [Figure S2A] 

*FS2C_Distances_between_positive_significant_pairs_location2.csv* quantifies the location pattern of each co-occuring pair. [Figure S2C]

*FS4B_5_dataset_defense_frequency1.xlsx* shows the difference of defense system composition across four orders and *E. coli*. It's used for Venn plot. [Figure S4B]

*F5S_plasmid_loss.xlsx* contains the results of plasmid loss assay. [Figure S5]

#### **figures**
This folder contains figures produced by scripts described above

*Ecoli_pagel_heatmap_pagel_0.005.png* main figure with Pagel correlation analysis results for 26k *E. coli* dataset

*[dataset]_heatmap_pagel_all_models_001.png* Pagel correlation analysis results for four order level datasets, where [dataset] can be either enter (Enterobacteriales) or burk (Burkholderiales) or baci (Bacillales) or pseu (Pseudomonadales)

*Pagel_compare_all_datasets.png* Comparison of correlations for all datasets

*EOP_and_AUC* contains figures generated by SynergisticVsAntagonistic_4pairs.R

*Ecoli_numContigs_vs_numDefence* contains output of Ecoli_numContigs_vc_numDefence.R

*Ecoli_phylogenetic_trees* contains outputs of Ecoli26k_tree_phylogroups_defense_systems.R and Draw_trees_for_selected_pairs.R

*Ecoli_phylogroups* contains output of Ecoli_defsystems_presence_by_phylogroup.R

*all_datasets_interaction_graph* contains output of Pagel_all_datasets_igraph.R

*cooccurence_vs_genome_distance* contains output of positive_vs_negative_genome_distance.R

*defence_systems_vs_phylogenetic_distance* contains output of defense_systems_content_vs_phylogenetic_distance.R

*tmn_Gabija_plasmids* contains output of 104_plasmids_cactus_visualization.R









