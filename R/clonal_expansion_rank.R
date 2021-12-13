#' Faceted barplot showing isotype and frequency of top 10 expanded clones per sample
#' @param dataset_summary Data of sequneces and clone count/freq
#' @param color_scale Color scale for isotypes
#' @param path_figures Path to which figures should be saved.
#' @return TRUE (plots saved as pdfs into path_figures)
#' @examples
#' clonal_expansion_rank(dataset_summary,,color_scale,path_figures="./figures")
#'

clonal_expansion_rank <- function(dataset_summary,color_scale=c("#edf8fb", "#b3cde3", "#8c96c6", "#88419d", "black"),figure_path){

  order_ids<-c("F5922470","F5922432","F5921513","F5921620","F5921915","F5922132","F5921597","F5921393","F5921788","F5921795","F5921407","F5921819","F5922956","F5921486","F5921495","F5922145")

  #pick top 10 from each
  dataset_summary_top_10 <- dataset_summary %>% group_by(dataset_id) %>% top_n(-10)

  dataset_summary_top_10 <- reshape2::melt(data.table::as.data.table(dataset_summary_top_10),
                                           id.vars = c('dataset_id', 'clonotype', 'Expansion', 'Rank'), measure.vars = c('IGHA', 'IGHD', 'IGHG', 'IGHM', 'IGHE', 'Missing'),
                                           value.name = 'Count', variable.name = 'Isotype')

  dataset_summary_top_10$dataset_id <- as.factor(dataset_summary_top_10$dataset_id)

  dataset_summary_top_10<-dataset_summary_top_10[dataset_summary_top_10$Isotype!="IGHE",]

  dataset_summary_top_10$c_call<-factor(dataset_summary_top_10$Isotype,levels=c("IGHA","IGHD","IGHG","IGHM","Missing"))

  dataset_summary_top_10$dataset_id<-factor(dataset_summary_top_10$dataset_id,levels=order_ids)

  #rename from fid to other datasetid for consistency in package plots
  dict_id <- dataset_summary[,c("dataset_id","dataset_id_old")]
  dict_id <- dict_id[!duplicated(dict_id),]
  dataset_summary_top_10_renamed <- merge(dataset_summary_top_10,dict_id)


  plt_rank <- ggplot2::ggplot(dataset_summary_top_10_renamed, ggplot2::aes(x = Rank, y = Count, fill=Isotype))+
    ggplot2::facet_wrap(~dataset_id_old, scales='free') +
    ggplot2::geom_bar(stat = 'identity', position = 'stack',color="black")+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=color_scale)+
    ggplot2::scale_x_continuous(breaks = c(1:10))+
    ggplot2::scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    ggplot2::geom_hline(yintercept=1, linetype="dashed", color = "black")+
    ggplot2::labs(y = 'Cell Count')+
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12),axis.text.x = ggplot2::element_text(size=15, angle = 0),axis.title.x = ggplot2::element_text(size=20),axis.title.y = ggplot2::element_text(size=20),axis.text.y = ggplot2::element_text(size=15, angle = 0))


  grDevices::pdf(file.path(figure_path,'ClonalExpansionRank.pdf'), family="Helvetica", width = 12,  height = 8)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(plt_rank, vp=vplayout(1,1))   # plot for row 1, column 1
  grDevices::dev.off()

  return(TRUE)

}


