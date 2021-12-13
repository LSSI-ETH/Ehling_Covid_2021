#' Plot comparing V,J usage across repertoires
#' @param list_repertoires List of repertoires
#' @param rep_category Metadata
#' @param chain Heavy (H) or Light (L)
#' @param freq_threshold threshold above which V,J genes are plotted
#' @param unique_rows If true all duplicated rows are excluded
#' @param size_x_axis Size of text on x-axis
#' @param size_y_axis Size of text on y-axis
#' @param figure_path Path to which figure is saved
#' @param plot_name ID of analysis for plot naming
#' @param color_scale Color scale to be applied
#'
#' @return TRUE (plots saved as pdfs into path_figures)
#' @examples
#' vj_across_reps(list_repertoires,rep_category,chain="H",freq_threshold=0.01,unique_rows=TRUE,size_x_axis=5,size_y_axis=5,figure_path,plot_name="",color_scale)
#'

vj_across_reps <- function(list_repertoires,rep_category,chain="H",freq_threshold=0.01,unique_rows=TRUE,size_x_axis=5,size_y_axis=5,figure_path,plot_name="",color_scale){

  germline_genes_present<-list()
  germline_genes_present[["V"]]<-c()
  germline_genes_present[["J"]]<-c()
  for(j in 1:length(list_repertoires)){
    germline_genes_present_curr<-list()
    if(chain=="H"){
      germline_genes_present_curr[["V"]]<-unique(as.character(list_repertoires[[j]]$v_call))
      germline_genes_present_curr[["J"]]<-unique(as.character(list_repertoires[[j]]$j_call))
    }else if(chain =="L"){
      germline_genes_present_curr[["V"]]<-unique(as.character(list_repertoires[[j]]$v_call.light))
      germline_genes_present_curr[["J"]]<-unique(as.character(list_repertoires[[j]]$j_call.light))
    }
    germline_genes_present[["V"]]<-c(germline_genes_present[["V"]],germline_genes_present_curr[["V"]][!(germline_genes_present_curr[["V"]] %in% germline_genes_present[["V"]])])
    germline_genes_present[["J"]]<-c(germline_genes_present[["J"]],germline_genes_present_curr[["J"]][!(germline_genes_present_curr[["J"]] %in% germline_genes_present[["J"]])])
  }

  germline_genes_present[["V"]]<-sort(germline_genes_present[["V"]])
  germline_genes_present[["J"]]<-sort(germline_genes_present[["J"]])


  list_germline_genes_per_sample <- list()
  for(z in 1:length(list_repertoires)){

    repertoire_name<-names(list_repertoires)[z]

    curr_repertoire<-list_repertoires[[z]]

    if(unique_rows==TRUE){
      curr_repertoire<- curr_repertoire[!duplicated(curr_repertoire),]
    }

    curr_repertoire$name_repertoire<-repertoire_name

    light<-FALSE
    list_order_genes<-list()
    #evaluate VDJ occurrence
    if(chain=="H"){
      v_genes<-factor(curr_repertoire$v_call,levels=germline_genes_present[["V"]])
      j_genes<-factor(curr_repertoire$j_call,levels=germline_genes_present[["J"]])
    }else if(chain=="L"){
      v_genes<-factor(curr_repertoire$v_call.light,levels=germline_genes_present[["V"]])
      j_genes<-factor(curr_repertoire$j_call.light,levels=germline_genes_present[["J"]])
    }
    occurrence_v <- as.data.frame(table(v_genes))
    occurrence_v$Freq <- occurrence_v$Freq/sum(occurrence_v$Freq)
    occurrence_v$gene <- "V"
    list_order_genes[["V"]]<-germline_genes_present[["V"]]

    names(occurrence_v)<-c("germline_gene","Freq","gene")

    occurrence_j <- as.data.frame(table(j_genes))
    occurrence_j$Freq <- occurrence_j$Freq/sum(occurrence_j$Freq)
    occurrence_j$gene <- "J"
    list_order_genes[["J"]]<-germline_genes_present[["J"]]

    names(occurrence_j)<-c("germline_gene","Freq","gene")

    vdj_occurrence <- base::rbind(occurrence_v,occurrence_j)
    vdj_occurrence$repertoire <- as.character(unique(curr_repertoire$name_repertoire))[!is.na(as.character(unique(curr_repertoire$name_repertoire)))]

    name_rep_vdj<-vdj_occurrence$repertoire[[1]]

    #make ready for plot by ordering V,D,J
    vdj_occurrence_plot_df <- vdj_occurrence
    vdj_occurrence_plot_df$gene <- factor(vdj_occurrence_plot_df$gene,levels=c("V","J"))

    #for each of V,D,J create plot
    vdj_occurrence_list<-list()
    for(i in 1:length(unique(vdj_occurrence_plot_df$gene))){
      curr_vdj_occurrence_plot_df <- vdj_occurrence_plot_df[vdj_occurrence_plot_df$gene==unique(vdj_occurrence_plot_df$gene)[i],]

      curr_vdj_occurrence_plot_df$germline_gene<-factor(curr_vdj_occurrence_plot_df$germline_gene,levels=list_order_genes[[unique(vdj_occurrence_plot_df$gene)[i]]])

      #add variable definitions pre plotting.
      germline_gene <- Freq <- NULL
      vdj_occurrence_list[[i]] <- curr_vdj_occurrence_plot_df


    }
    vdj_occurrence_list_df <- do.call(rbind,vdj_occurrence_list)

    list_germline_genes_per_sample[[name_rep_vdj]] <- vdj_occurrence_list_df


  }
  single_df <- do.call(rbind,list_germline_genes_per_sample)

  single_df_merged<-merge(single_df,rep_category,by="repertoire")

  single_sum<-single_df_merged %>%
    group_by(category,germline_gene,gene) %>%
    summarize(mean_freq=mean(Freq),sd=sd(Freq),se=plotrix::std.error(Freq))


  top_gene_df<-as.data.frame(single_df_merged %>%
                               dplyr::group_by(germline_gene,gene) %>%
                               dplyr::summarize(mean_freq=mean(Freq),sd=sd(Freq),se=plotrix::std.error(Freq)))

  top_genes<-as.character(top_gene_df[top_gene_df$mean_freq>freq_threshold,"germline_gene"])


  single_sum_df<-as.data.frame(single_sum)
  names(single_sum_df)<-c("Category","Gene","VJ","Freq","sd","se")


  single_sum_df$sd[is.na(single_sum_df$sd)]<-0
  single_sum_df$se[is.na(single_sum_df$se)]<-0

  dodge <- ggplot2::position_dodge(width=0.9)

  limits <- ggplot2::aes(ymax = Freq + se, ymin=Freq - se) #shell for calculation of limit bars.

  single_sum_df<-single_sum_df[single_sum_df$Gene %in% top_genes,]

  single_sum_df_V<-single_sum_df[single_sum_df$VJ=="V",]
  single_sum_df_J<-single_sum_df[single_sum_df$VJ=="J",]

  vdj_occurrence_plot_comp_V <- ggplot2::ggplot(data=single_sum_df_V, ggplot2::aes(Gene, Freq,fill=Category))+
    ggplot2::geom_bar(stat="identity",position="dodge",colour='black') +
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=color_scale)+#3rd pos))+
    ggplot2::geom_errorbar(limits,position=dodge,width=0.25)+
    ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=size_x_axis, angle = 90,vjust=0.5), axis.text.y = ggplot2::element_text(size=size_y_axis))

  vdj_occurrence_plot_comp_J <- ggplot2::ggplot(data=single_sum_df_J, ggplot2::aes(Gene, Freq,fill=Category))+
    ggplot2::geom_bar(stat="identity",position="dodge",colour='black') +
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=color_scale)+#3rd pos))+
    ggplot2::geom_errorbar(limits,position=dodge,width=0.25)+
    ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=size_x_axis, angle = 90,vjust=0.5), axis.text.y = ggplot2::element_text(size=size_y_axis))


  name_plot_v <- paste("v_occurrence_",plot_name,".pdf",sep="")

  grDevices::pdf(file.path(figure_path, name_plot_v),  width = 15,  height = 5)
  grid::pushViewport( grid::viewport(layout=grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(vdj_occurrence_plot_comp_V, vp=vplayout(1,1))
  grDevices::dev.off()

  name_plot_j <- paste("j_occurrence_",plot_name,".pdf",sep="")

  grDevices::pdf(file.path(figure_path, name_plot_j),  width = 15,  height = 5)
  grid::pushViewport( grid::viewport(layout=grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(vdj_occurrence_plot_comp_J, vp=vplayout(1,1))
  grDevices::dev.off()

  return(TRUE)

}
