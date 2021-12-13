#' Heatmap of mAb enrichment using geom_tile
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
#' enrichment_heatmap_geom_tile(compare_rep_names=c("nCoV36_S1","nCoV36_S2"),reference_rep_name="nCoV36_Ab(Background)",ordered_by="S2",ref_seqs_name="ref_n36",names_cohort="default",mode="multi",figure_path)


enrichment_heatmap_geom_tile<-function(compare_rep_names=c("nCoV36_S1","nCoV36_S2"),reference_rep_name="nCoV36_Ab(Background)",ref_seqs_name="ref_n36",names_cohort="default",ordered_by="S2",text_size=c(5,5),figure_path){
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

  s1_enrich<-enrichment_df_S1[,c("enrichment_S1_Ab","aaSeqCDR3")]
  s2_enrich<-enrichment_df_S2[,c("enrichment_S2_Ab","aaSeqCDR3")]
  s1_enrich$cat<-"S1"
  s2_enrich$cat<-"S2"
  names(s1_enrich)<-c("enrichment","aaSeqCDR3","cat")
  names(s2_enrich)<-c("enrichment","aaSeqCDR3","cat")

  if(ordered_by=="S2"){
    order_col<-s2_enrich[order(s2_enrich$enrichment),"aaSeqCDR3"]
  }else{
    order_col<-s1_enrich[order(s1_enrich$enrichment),"aaSeqCDR3"]
  }

  htmap_enrich<-bind_rows(s1_enrich,s2_enrich)

  htmap_enrich$aaSeqCDR3<-factor(htmap_enrich$aaSeqCDR3,levels=order_col)

  htmap_enrich$cat<-factor(htmap_enrich$cat,levels=c("S2","S1"))

  ht_map<-ggplot2::ggplot(htmap_enrich, ggplot2::aes(aaSeqCDR3, cat, fill= enrichment)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient2(low="blue", mid="white", high="red",
                                  midpoint = 0,
                                  limits=c(-11,7),
                                  breaks=c(-11,0,7)
    )+
    ggplot2::labs(fill = "Enrichment") +
    ggplot2::theme_minimal()+
    ggplot2::coord_equal()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=text_size[1],angle = 45, vjust = 1.1, hjust=1),  axis.text.y = ggplot2::element_text(size=text_size[2],angle = 0, vjust = 0.5),axis.title = ggplot2::element_blank())


  grDevices::pdf(file.path(figure_path, paste("enrichment_test_",names_cohort,".pdf",sep="")), family="Helvetica", width = 10,  height = 7)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(ht_map, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()



  return(TRUE)
}
