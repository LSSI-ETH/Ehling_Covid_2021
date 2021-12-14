#' Barplot isotype usage
#' @param list_repertoires Patient Metadata.
#' @param figure_path Path to which figures should be saved.
#' @return TRUE (plots saved as pdfs into figure_path)
#' @examples
#' plot_isotype_usage_bar(list_repertoires=example_data, figure_path="")

barplot_isotype_usage<-function (list_repertoires, order_ids,chain = "H", sample_name, text_size = 5, width_height = c(x = 10, y = 10), figure_path = "",color_scale_category=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f')){

  sample_name<-"all"

  chain<-"H"

  bar_chart_df_plot_list<-list()
  for(i in 1:length(list_repertoires)){
    c_call_curr<-list_repertoires[[i]][!duplicated(list_repertoires[[i]]),]

    if(chain=="H"){
      isotype_occ<-sapply(1:nrow(c_call_curr),function(x) substr(as.character(c_call_curr$c_call[x]),1,4))
      isotype_occ[isotype_occ == "Miss"]<-"Missing"
      isotype_occ[isotype_occ==""] <- "Missing"
      isotype_occ <- isotype_occ[isotype_occ %in% c("IGHA","IGHD","IGHG","IGHM","IGHE","Missing")]
    }else if(chain=="L"){
      isotype_occ<-sapply(1:nrow(c_call_curr),function(x) substr(as.character(c_call_curr$c_call.light[x]),1,4))
      isotype_occ[isotype_occ == "Miss"]<-"Missing"
      isotype_occ[isotype_occ==""] <- "Missing"
      isotype_occ <- isotype_occ[isotype_occ %in% c("IGKC","IGLC","Missing")]
    }


    isotype_df_h<-reshape2::melt(100*table(isotype_occ)/sum(table(isotype_occ)))
    isotype_df_h$class<-names(list_repertoires)[i]
    names(isotype_df_h)<-c("isotype","Freq","class")

    bar_chart_df<-isotype_df_h

    if(chain=="H"){
      bar_chart_df$category<-factor(bar_chart_df$isotype,levels=c("IGHA","IGHD","IGHE","IGHG","IGHM","Missing"))
      color_scale <- color_scale_category
      names(color_scale) <- c("IGHA","IGHD","IGHG","IGHM","IGHE","Missing")
    }else if(chain=="L"){
      bar_chart_df$category<-factor(bar_chart_df$isotype,levels=c("IGKC","IGLC","Missing"))
      color_scale <- c('#8da0cb','#fc8d62','grey')
      names(color_scale) <- c("IGKC","IGLC","Missing")
    }

    bar_chart_df_plot_tmp <- bar_chart_df %>%
      arrange(desc(category)) %>%
      mutate(prop = Freq) %>%
      mutate(ypos = cumsum(prop)- 0.5*prop )

    bar_chart_df_plot_list[[i]]<-bar_chart_df_plot_tmp
  }

  c_class_freq_plot<-do.call(rbind,bar_chart_df_plot_list)
  c_class_freq_plot$prop<-NULL
  c_class_freq_plot$ypos<-NULL


  if(chain=="H"){
    c_class_freq_plot$category<-factor(c_class_freq_plot$isotype,levels=c("IGHA","IGHD","IGHE","IGHG","IGHM","Missing"))
    color_scale <- color_scale_category
    names(color_scale) <- c("IGHA","IGHD","IGHG","IGHM","IGHE","Missing")
  }else if(chain=="L"){
    c_class_freq_plot$category<-factor(c_class_freq_plot$isotype,levels=c("IGKC","IGLC","Missing"))
    color_scale <- c('#8da0cb','#fc8d62','grey')
    names(color_scale) <- c("IGKC","IGLC","Missing")
  }


  names(c_class_freq_plot)<-c("isotype","Freq","repertoire","isotype_category")

  c_class_freq_plot_merged<-merge(c_class_freq_plot,sample_overview,by="repertoire")
  c_class_freq_plot_sum<-as.data.frame(c_class_freq_plot_merged %>%
                                         group_by(category,isotype) %>%
                                         summarize(mean_freq=mean(Freq),sd=sd(Freq),se=plotrix::std.error(Freq)))

  names(c_class_freq_plot_sum)<-c("Category","Isotype","Freq","sd","se")

  c_class_freq_plot_sum$sd[is.na(c_class_freq_plot_sum$sd)]<-0
  c_class_freq_plot_sum$se[is.na(c_class_freq_plot_sum$se)]<-0
  c_class_freq_plot_sum$Isotype<-factor(c_class_freq_plot_sum$Isotype,levels=c("IGHA","IGHD","IGHG","IGHM","IGHE","Missing"))

  dodge <- position_dodge(width=0.9)
  limits <- aes(ymax = Freq + se, ymin=Freq - se) #shell for calculation of limit bars.

  isotype_usage_plot<- ggplot2::ggplot(data=c_class_freq_plot_sum, ggplot2::aes(Category, Freq,fill=Isotype))+
    ggplot2::geom_bar(stat="identity",position="dodge",colour='black') +
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=color_scale)+
    ggplot2::geom_errorbar(limits,position=dodge,width=0.25)+
    ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=16, angle = 0),axis.title.y = ggplot2::element_text(size=16),axis.text.y = ggplot2::element_text(size=16, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Isotype_usage_category_",sample_name,"_bar.pdf",sep="")),  width = 5,  height = 5)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(isotype_usage_plot, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()

  return(TRUE)


}
