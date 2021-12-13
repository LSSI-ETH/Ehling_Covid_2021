#' PLots comparing heavy and light chain usage
#' @param paired_df Dataframe with paired seqdata
#' @param subset ID for naming
#' @param figure_path Path to which figures should be saved.
#' @return TRUE (plots saved as pdfs into path_figures)
#' @examples
#' heavy_light_comparisons(paired_df=paired_df_example,subset="all",figure_path)
#'

heavy_light_comparisons<-function(paired_df=paired_df_example,subset="all",figure_path){

  paired_df$v_call<-as.character(paired_df$v_call)
  paired_df$v_call.light<-as.character(paired_df$v_call.light)

  unique_vh_calls<-sort(unique(as.character(paired_df$v_call)))
  unique_vl_calls<-sort(unique(as.character(paired_df$v_call.light)))

  vh_vl_mat = matrix(0, length(unique_vh_calls), length(unique_vl_calls)) # gives you a 2 x 3 matrix filled with the numbers 1 to 6
  rownames(vh_vl_mat)<-unique_vh_calls
  colnames(vh_vl_mat)<-unique_vl_calls

  for(i in 1:nrow(paired_df)){
    curr_vh<-paired_df$v_call[i]
    curr_vl<-paired_df$v_call.light[i]

    vh_vl_mat[curr_vh,curr_vl]<-vh_vl_mat[curr_vh,curr_vl]+1
  }

  vh_vl_mat<-t(vh_vl_mat)


  ht_map<-ComplexHeatmap::Heatmap(vh_vl_mat,name = "Counts",
                                  cluster_rows = TRUE,
                                  cluster_row_slices = TRUE,
                                  cluster_columns = TRUE,
                                  cluster_column_slices = TRUE,
                                  col = circlize::colorRamp2(c(0,max(vh_vl_mat)/2,max(vh_vl_mat)), c("blue","yellow","red")),
                                  row_names_gp = grid::gpar(fontsize = 5),
                                  column_names_gp = grid::gpar(fontsize = 5))

  pdf(file.path(figure_path,paste("heatmap_v_hl_clustered_",subset,".pdf",sep="")), width = 7,  height = 6.5,onefile=FALSE)
  ComplexHeatmap::draw(ht_map,heatmap_legend_side = "right")
  dev.off()


  col_freqs<-function(col){
    col/sum(col)
  }

  vh_vl_mat_Freq<-apply(vh_vl_mat, 2, col_freqs)
  col_counts<-apply(vh_vl_mat, 2, sum)
  row_counts<-apply(vh_vl_mat, 1, sum)

  row.names(vh_vl_mat_Freq)<-paste(row.names(vh_vl_mat_Freq)," (",row_counts," seqs)",sep="")
  colnames(vh_vl_mat_Freq)<-paste(colnames(vh_vl_mat_Freq)," (",col_counts," seqs)",sep="")


  ht_map<-ComplexHeatmap::Heatmap(vh_vl_mat_Freq,name = "Frequency\nper VH-gene",
                                  cluster_rows = FALSE,
                                  cluster_row_slices = FALSE,
                                  cluster_columns = FALSE,
                                  cluster_column_slices = FALSE,
                                  col = circlize::colorRamp2(c(0,0.5,1), c("blue","yellow","red")),
                                  row_names_gp = grid::gpar(fontsize = 5),
                                  column_names_gp = grid::gpar(fontsize = 5))

  pdf(file.path(figure_path,paste("heatmap_v_hl_colfreqs_",subset,".pdf",sep="")), width = 8,  height = 6.5,onefile=FALSE)
  ComplexHeatmap::draw(ht_map,heatmap_legend_side = "right")
  dev.off()


  len_cdr3s<-data.frame(CDRH3_len=nchar(as.character(paired_df$junction_aa.x)),
                        CDRL3_len=nchar(as.character(paired_df$junction_aa.y)),
                        c_call=paired_df$c_call,
                        dataset_id="example", stringsAsFactors=FALSE)

  #duplicate dataframe to add category all that combines all constant regions
  len_cdr3s_all<-len_cdr3s
  len_cdr3s_igg_iga<-len_cdr3s[len_cdr3s_all$c_call %in% c("IGHM","IGHA"),]
  len_cdr3s_all$c_call<-"All"
  len_cdr3s_full_plot<-rbind(len_cdr3s_igg_iga,len_cdr3s_all)

  scale_xy <- max(1.2*max(len_cdr3s_full_plot$CDRH3_len),1.2*max(len_cdr3s_full_plot$CDRL3_len))
  my_spectral <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,'Spectral'))(14)

  HL_length_scatter<-ggplot2::ggplot(len_cdr3s_full_plot, ggplot2::aes(x=CDRH3_len,y=CDRL3_len)) +
    ggplot2::geom_point(position=ggplot2::position_jitter(h=0.1, w=0.1),shape = 21,alpha = 0.4, size = 2, color=my_spectral[3])+
    ggplot2::theme_bw()+
    ggplot2::facet_wrap(.~c_call)+
    ggplot2::scale_x_continuous(breaks=seq(0,max(len_cdr3s_full_plot$CDRH3_len),1))+
    ggplot2::scale_y_continuous(breaks=seq(0,max(len_cdr3s_full_plot$CDRL3_len),1))#+


  grDevices::pdf(file.path(figure_path,"length_cdrhl_example.pdf"), family="Helvetica", width = 13,  height = 11)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(HL_length_scatter, vp=vplayout(1,1))   # plot for row 1, column 1
  grDevices::dev.off()


  return(TRUE)

}
