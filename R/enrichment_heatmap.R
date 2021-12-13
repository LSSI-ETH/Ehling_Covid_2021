#' Heatmap of mAb enrichment
#' @param compare_rep_names 2d vector containing names of samples
#' @param reference_rep_name Name of dataset to compare against
#' @param ordered_by Entry by which heatmap should be ordered
#' @param ref_seqs_name Reference file (background detection)
#' @param names_cohort ID to add to filename
#' @param mode Way in which analysis is performed (pairwise or multi)
#' @param figure_path Path to which figures are saved
#'
#' @return TRUE
#' #' @examples
#' enrichment_heatmap_standard(compare_rep_names=c("nCoV36_S1","nCoV36_S2"),reference_rep_name="nCoV36_Ab(Background)",ordered_by="S2",ref_seqs_name="ref_n36",names_cohort="default",mode="multi",figure_path)


enrichment_heatmap_standard<-function(compare_rep_names=c("nCoV36_S1","nCoV36_S2"),reference_rep_name="nCoV36_Ab(Background)",ordered_by="S2",ref_seqs_name="ref_n36",names_cohort="default",mode="multi",figure_path){

   #calculate enrichment for S1 and S2 (36)
   enrichment_df_S1 <- .calculate_enrichment(rep_1_name=compare_rep_names[[1]],rep_2_name=reference_rep_name)
   enrichment_df_S2 <- .calculate_enrichment(rep_1_name=compare_rep_names[[2]],rep_2_name=reference_rep_name)

   ref_seqs<-enrichment_data[[ref_seqs_name]]
   enrichment_df_S1<-enrichment_df_S1[enrichment_df_S1$aaSeqCDR3 %in% ref_seqs,]
   enrichment_df_S2<-enrichment_df_S2[enrichment_df_S2$aaSeqCDR3 %in% ref_seqs,]

   enrichment_df_S1 <- enrichment_df_S1[,c("aaSeqCDR3", "enrichment")]
   enrichment_df_S2 <- enrichment_df_S2[,c("aaSeqCDR3", "enrichment")]

   names(enrichment_df_S1)<-c("aaSeqCDR3", "enrichment_S1_Ab")
   names(enrichment_df_S2)<-c("aaSeqCDR3", "enrichment_S2_Ab")

   if(ordered_by=="S2"){
     enrichment_df<-merge(enrichment_df_S1,enrichment_df_S2,by="aaSeqCDR3")[,c("aaSeqCDR3", "enrichment_S2_Ab", "enrichment_S1_Ab")]
     names(enrichment_df)<-c("aaSeqCDR3","S2_v_Ab","S1_v_Ab")
   }else{
     enrichment_df<-merge(enrichment_df_S1,enrichment_df_S2,by="aaSeqCDR3")[,c("aaSeqCDR3", "enrichment_S1_Ab", "enrichment_S2_Ab")]
     names(enrichment_df)<-c("aaSeqCDR3","S1_v_Ab","S2_v_Ab")
  }

  if(mode=="pairwise"){
    sequences_rep_1_ordered<-as.character(enrichment_df[order(-enrichment_df$cloneCount_rep1),"aaSeqCDR3"])
    sequences_rep_2_ordered<-as.character(enrichment_df[order(-enrichment_df$cloneCount_rep2),"aaSeqCDR3"])

    seqs_rep1_rep2 = matrix(0, length(sequences_rep_1_ordered), length(sequences_rep_2_ordered)) # gives you a 2 x 3 matrix filled with the numbers 1 to 6
    rownames(seqs_rep1_rep2)<-sequences_rep_1_ordered
    colnames(seqs_rep1_rep2)<-sequences_rep_2_ordered

    for(i in 1:nrow(enrichment_df)){

      seqs_rep1_rep2[as.character(enrichment_df$aaSeqCDR3[i]),as.character(enrichment_df$aaSeqCDR3[i])]<-enrichment_df$enrichment[i]
    }


    ht_map<-ComplexHeatmap::Heatmap(seqs_rep1_rep2,name = "Enrichment",
                                    cluster_rows = FALSE,
                                    cluster_row_slices = FALSE,
                                    cluster_columns = FALSE,
                                    cluster_column_slices = FALSE,
                                    col = circlize::colorRamp2(c(min(seqs_rep1_rep2),mean(seqs_rep1_rep2),max(seqs_rep1_rep2)), c("blue","white","red")),
                                    row_names_gp = grid::gpar(fontsize = 5),
                                    column_names_gp = grid::gpar(fontsize = 5),
                                    rect_gp = grid::gpar(col = "white", lwd = 2))

    pdf(paste("figures/heatmap_enrichment_",names_cohort,"_flip.pdf",sep=""), width = 7,  height = 6.5,onefile=FALSE)
    ComplexHeatmap::draw(ht_map,heatmap_legend_side = "right")
    dev.off()
  }else{

    enrichment_df<-enrichment_df[order(enrichment_df[,2]),]

    enrichment_df<-enrichment_df[,c(1,3,2)]
    enrichment_heat = matrix(0, length(enrichment_df$aaSeqCDR3), ncol(enrichment_df)-1) # gives you a 2 x 3 matrix filled with the numbers 1 to 6
    rownames(enrichment_heat)<-enrichment_df$aaSeqCDR3
    enrichment_cols<-names(enrichment_df)[2:ncol(enrichment_df)]
    colnames(enrichment_heat)<-enrichment_cols



    for(i in 1:nrow(enrichment_df)){
      for(j in 2:ncol(enrichment_df)){
        enrichment_heat[as.character(enrichment_df$aaSeqCDR3[i]),enrichment_cols[j-1]]<-enrichment_df[i,j]
      }
    }

    enrichment_heat<-t(enrichment_heat)


    ht_map<-ComplexHeatmap::Heatmap(enrichment_heat,name = "Enrichment (log2)",
                                    cluster_rows = FALSE,
                                    cluster_row_slices = FALSE,
                                    cluster_columns = FALSE,
                                    cluster_column_slices = FALSE,
                                    col = circlize::colorRamp2(c(min(enrichment_heat),mean(enrichment_heat),max(enrichment_heat)), c("blue","white","red")),
                                    row_names_gp = grid::gpar(fontsize = 12),
                                    column_names_gp = grid::gpar(fontsize = 8),
                                    rect_gp = grid::gpar(col = "white", lwd = 2))

    pdf(file.path(figure_path,paste("enrichment_multi_",names_cohort,"_flip.pdf",sep="")), width = 11,  height = 5,onefile=FALSE)
    ComplexHeatmap::draw(ht_map,heatmap_legend_side = "right")
    dev.off()

  }

  return(TRUE)
}
