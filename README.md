### Original code for Wu, Garushyants et al. "Synergistic anti-phage activity of bacterial immune systems"

#### **scripts**
Contain scripts to reproduce the analysis and vizualize data

*Ecoli_numContigs_vc_numDefence.R* builds correlation between number of contigs in genome and number of predicted immune systems

Scripts *Ecoli26k_tree_phylogroups_defense_systems.R* and *Visualize_phylogenetic_trees_orders.R* visualize phylognetic trees and add data to them.

*Ecoli_defsystems_presence_by_phylogroup.R* allows to analyze immune systems content in E. coli phylogroups and also performs enrichment analysis.

*defense_systems_content_vs_phylogenetic_distance.R* builds correlations between immune systems content and phylogenetic distance

*Pagel_calculate_sge_farm.R* is the main script to run the Pagel correlation test for the binary traits

*Pagel_visualization_Ecoli.R* and *Pagel_visualization_orders.R* calculate directionality of interactions for the significant results and visualize the data for 26k E. coli dataset and order level datasets accordingly.

*Pagel_all_datasets_igraph.R* builds graphs for Pagel results

*Pagel_compare_results_all_datasets.R* script that compares and vizualizes results of Pagel test for all five datasets

*Draw_trees_for_selected_pairs.R* visualizes phylogenetic trees for selected pairs of immune systems.

*104_plasmids_cactus_visualization.R* visualizes [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) output for 104 plasmids carrying tmn and Gabija systems








