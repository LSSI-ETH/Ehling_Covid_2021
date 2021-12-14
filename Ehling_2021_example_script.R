####
# Here we present a collection of code for the various analyses presented in Ehling 2021.
# For readability the visualizations are grouped by type (barplots, boxplots, network graphs)
# Additonally, a summary of the sequence extraction is presented.
#

library(Ehling2021)
library(ggplot2)
library(dplyr)
library(igraph)

FIGURE_PATH<-"./figures"


###
# Fig 1 Dataset overview 
###

#Print two-column plot including barplot and table 
plot_patients(patients=patients_metadata,
  figure_path=FIGURE_PATH)

###
# Dotplots: OD ratio
###





###
# For simplicity's sake the rest of the analysis are presented with subsets of the full data (50 seqs from each dataset)
###

###
# Comparative Barplots: CDR3 length, V-Usage, J usage, V identity
###


###
# Barplots CDRH3 length
###

cdr_length_across_reps(list_repertoires=example_data,
      rep_category=sample_overview,
      length_range=c(6:28),
      chain="H",
      unique_seqs=TRUE, 
      size_x_axis=10,
      width_height=c(x=15,y=15), 
      color_scale=c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6',"blue4"), 
      sample_name="",
      figure_path=FIGURE_PATH)

###
# Barplots VJ usage
###

vj_across_reps(list_repertoires=example_data,
      rep_category=sample_overview,
      chain="H",
      freq_threshold=0.01,
      unique_rows=TRUE,
      size_x_axis=12,
            size_y_axis=12,
      figure_path=FIGURE_PATH,
      plot_name="covid_vs_healthy",
      color_scale=c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6',"blue4"))


###
# Boxplot V_identity
###

boxplot_v_identity(list_repertoires=example_data[1:16], 
  rep_category=sample_overview,
  chain = "H", 
  color_scale = c("#d7191c","#fdae61", "#2c7bb6", "#ffffbf", "#abd9e9",'black'), 
  plot_name = "healthy_v_cov",
  figure_path = FIGURE_PATH)


###
# Barplot Isotype Usage
###

barplot_isotype_usage(list_repertoires=example_data[1:16],
  order_ids=c("F5922470","F5922432","F5921513","F5921620","F5921915","F5922132","F5921597","F5921393","F5921788","F5921795","F5921407","F5921819","F5922956","F5921486","F5921495","F5922145"),
  chain="H",
  sample_name="scCovid",
  text_size=3,
  width_height=c(x=10,y=10),
  figure_path=FIGURE_PATH,
  color_scale_category=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'))


###
# Piecharts: Figure Isotype Usage
###

pie_chart_isotype_usage(list_repertoires=example_data[1:16], 
                              chain = "H", 
                              order_ids=c("scCovid_1", "scCovid_2", "scCovid_3", "scCovid_4", "scCovid_5", "scCovid_6", "scCovid_7", "scCovid_8", "scCovid_9", "scCovid_10", "scCovid_11", "scCovid_12", "scCovid_13", "scCovid_14", "scCovid_15", "scCovid_16"),
                              sample_name='test', 
                              text_size = 4,
                              width_height=c(x=10,y=10), 
                              figure_path = FIGURE_PATH,
                              color_scale_category=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'))





###
# Heavy - Light chain comparison
###

heavy_light_comparisons(paired_df=paired_df_example,
  subset="all",
  figure_path=FIGURE_PATH)





##
# Clonal Expansion Rank
##
clonal_expansion_rank(dataset_summary,figure_path=FIGURE_PATH)


###
# Heatmaps: mAb enrichment
###

enrichment_heatmap_geom_tile(compare_rep_names=c("nCoV36_S1","nCoV36_S2"),
  reference_rep_name="nCoV36_Ab(Background)",
  ref_seqs_name="ref_n36",
  ordered_by="S2",
  names_cohort="",
  figure_path=FIGURE_PATH)



###
# Binder data analysis
###

supp_pool_analysis_132(chosen_132=chosen_132,
  figure_path=FIGURE_PATH)


###
# Network: Sequence similarity of binders and closest neighbors 
###
network_figure(binders_all=binders_all,
  binder_matches=binder_matches,
  figure_path=FIGURE_PATH)









