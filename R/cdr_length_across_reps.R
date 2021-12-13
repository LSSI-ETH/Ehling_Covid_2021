#' Plot comparing CDR3 lengths across repertoies
#' @param list_repertoires List of repertoires
#' @param rep_category Metadata
#' @param length_range CDR lengths to be considered
#' @param chain Heavy (H) or Light (L)
#' @param unique_seqs If true sequences are uniqued before length dist is calculated
#' @param size_x_axis Size of text on x-axis
#' @param width_height Width and height of plot
#' @param color_scale Color scale to be applied
#' @param sample_name ID of analysis for plot naming
#' @param figure_path Path to which figure is saved
#'
#' @return TRUE (plots saved as pdfs into path_figures)
#' @examples
#' cdr_length_across_reps(list_repertoires, rep_category, length_range=c(6:28),chain="H",unique_seqs=TRUE, size_x_axis=10,width_height=c(x=15,y=15), color_scale, sample_name="",figure_path="")
#'

cdr_length_across_reps <- function(list_repertoires, rep_category, length_range=c(6:28),chain="H",unique_seqs=TRUE, size_x_axis=10,width_height=c(x=15,y=15), color_scale, sample_name="",figure_path=""){

  list_cdr3_length_df<-list()
  for(i in 1:length(list_repertoires)){
    repertoire_name<-names(list_repertoires)[i]
    if(chain=="H"){
      sequences_curr <- as.character(list_repertoires[[i]]$junction_aa)
    }else if(chain=="L"){
      sequences_curr <- as.character(list_repertoires[[i]]$junction_aa.light)
    }
    if(unique_seqs==TRUE){
      sequences_curr <- unique(sequences_curr)
    }
    data<-nchar(sequences_curr)
    lengths_junction_aa<-reshape2::melt(setNames(tabulate(data),1:max(data)))
    lengths_junction_aa$length<-row.names(lengths_junction_aa)
    lengths_junction_aa$Freq<-100*lengths_junction_aa$value/sum(lengths_junction_aa$value)
    lengths_junction_aa$name_repertoire<-repertoire_name
    list_cdr3_length_df[[i]]<-lengths_junction_aa
  }
  cdr3_length_df<-do.call(rbind,list_cdr3_length_df)

  cdr3_length_df<-cdr3_length_df[as.integer(cdr3_length_df$length) %in% length_range,]

  cdr3_length_df$length<-factor(cdr3_length_df$length,levels=unique(cdr3_length_df$length))

  names(cdr3_length_df)<-c("value","length","Freq","repertoire")

  cdr3_length_df<-merge(cdr3_length_df,rep_category,by="repertoire")

  cdr3_length_df_sum<-as.data.frame(cdr3_length_df %>%
                                      dplyr::group_by(category,length) %>%
                                      dplyr::summarize(mean_freq=mean(Freq),sd=sd(Freq),se=plotrix::std.error(Freq)))

  names(cdr3_length_df_sum)<-c("Category","length","Freq","sd","se")

  dodge <- position_dodge(width=0.9)
  limits <- aes(ymax = Freq + se, ymin=Freq - se) #shell for calculation of limit bars.

  cdr3_length <- ggplot2::ggplot(data=cdr3_length_df_sum, ggplot2::aes(length, Freq,fill=Category))+
    ggplot2::geom_bar(stat="identity",position="dodge",colour="black") +
    ggplot2::theme_bw()+
    ggplot2::geom_errorbar(limits,position=dodge,width=0.25)+
    ggplot2::scale_fill_manual(values=color_scale)+#3rd pos))+
    ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=size_x_axis, angle = 0,vjust=0.5))

  grDevices::pdf(file.path(figure_path, paste("CDR3_length_",sample_name,".pdf",sep="")),  width = width_height[["x"]],  height = width_height[["y"]])
  grid::pushViewport( grid::viewport(layout=grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(cdr3_length, vp=vplayout(1,1))
  grDevices::dev.off()


}
