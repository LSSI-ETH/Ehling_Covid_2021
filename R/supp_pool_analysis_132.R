#' Generates supplementary plots characterizing the 132 chosen sequences
#' @param chosen_132 Dataframe of 132 chosen sequences (96+36)
#' @param figure_path Path where figures are saved
#' @return TRUE (plots saved as pdfs into path_figures)
#' @examples
#' supp_pool_analysis_132(chosen_132,binder_matches,figure_path="")
#'

supp_pool_analysis_132 <- function(chosen_132,figure_path){
  pool_a<-chosen_132[chosen_132$batch.1=="Pool A",]
  pool_b<-chosen_132[chosen_132$batch.1=="Pool B",]

  list_pools<-list()
  list_pools[["Pool A"]]<-pool_a
  list_pools[["Pool B"]]<-pool_b


  v_usage_A<-reshape2::melt(100*table(pool_a$v_gene)/sum(table(pool_a$v_gene)))
  v_usage_A$category<-"Pool A"

  v_usage_B<-reshape2::melt(100*table(pool_b$v_gene)/sum(table(pool_b$v_gene)))
  v_usage_B$category<-"Pool B"

  add_to_A<-unique(as.character(v_usage_B$Var1[!(v_usage_B$Var1 %in% v_usage_A$Var1)]))
  add_to_B<-unique(as.character(v_usage_A$Var1[!(v_usage_A$Var1 %in% v_usage_B$Var1)]))




  pools_Vusage<-rbind(v_usage_A,v_usage_B)

  pools_Vusage_all<-rbind(pools_Vusage,rbind(data.frame(Var1=add_to_A,value=0,category="Pool A"),data.frame(Var1=add_to_B,value=0,category="Pool B")))

  ordered_v_genes<-levels(pools_Vusage_all$Var1)[order(levels(pools_Vusage_all$Var1))]

  pools_Vusage_all$Var1<-factor(pools_Vusage_all$Var1,levels=ordered_v_genes)

  pools_Vusage_all_combo<-pools_Vusage_all
  pools_Vusage_all_combo$category<-"Pool A+B"
  pools_Vusage_all_combo2<-as.data.frame(pools_Vusage_all_combo %>%
                                           group_by(Var1) %>%
                                           summarize(value=sum(value),category=category))


  pools_Vusage_all_combo2<-pools_Vusage_all_combo2[!duplicated(pools_Vusage_all_combo2),]
  pools_Vusage_all_combo2$value<-pools_Vusage_all_combo2$value/2

  grid_v_usage<-rbind(pools_Vusage_all,pools_Vusage_all_combo2)
  #



  v_usage_plot_combo<- ggplot2::ggplot(data=grid_v_usage, ggplot2::aes(Var1, value,fill=category))+
    ggplot2::geom_bar(stat="identity",position="dodge",colour='black') +
    ggplot2::facet_grid(.~category)+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=c("gray30","red","darkgrey"))+
    ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=5, angle = 90,vjust=0.5),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_Analysis_Vusage_pools_bar.pdf",sep="")),  width = 10,  height = 3)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(v_usage_plot_combo, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()


  #j_gene

  v_usage_A<-reshape2::melt(100*table(pool_a$j_gene)/sum(table(pool_a$j_gene)))
  v_usage_A$category<-"Pool A"

  v_usage_B<-reshape2::melt(100*table(pool_b$j_gene)/sum(table(pool_b$j_gene)))
  v_usage_B$category<-"Pool B"


  add_to_A<-unique(as.character(v_usage_B$Var1[!(v_usage_B$Var1 %in% v_usage_A$Var1)]))
  add_to_B<-unique(as.character(v_usage_A$Var1[!(v_usage_A$Var1 %in% v_usage_B$Var1)]))


  pools_Vusage<-rbind(v_usage_A,v_usage_B)
  pools_Vusage_all<-rbind(pools_Vusage,data.frame(Var1=add_to_A,value=0,category="Pool A"))

  ordered_v_genes<-levels(pools_Vusage_all$Var1)[order(levels(pools_Vusage_all$Var1))]

  pools_Vusage_all$Var1<-factor(pools_Vusage_all$Var1,levels=ordered_v_genes)
  pools_Vusage_all_combo<-pools_Vusage_all
  pools_Vusage_all_combo$category<-"Pool A+B"
  pools_Vusage_all_combo2<-as.data.frame(pools_Vusage_all_combo %>%
                                           group_by(Var1) %>%
                                           summarize(value=sum(value),category=category))

  pools_Vusage_all_combo2<-pools_Vusage_all_combo2[!duplicated(pools_Vusage_all_combo2),]
  pools_Vusage_all_combo2$value<-pools_Vusage_all_combo2$value/2

  grid_v_usage<-rbind(pools_Vusage_all,pools_Vusage_all_combo2)
  #

  v_usage_plot_combo<- ggplot2::ggplot(data=grid_v_usage, ggplot2::aes(Var1, value,fill=category))+
    ggplot2::geom_bar(stat="identity",position="dodge",colour='black') +
    ggplot2::facet_grid(.~category)+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=c("gray30","red","darkgrey"))+
    ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=8, angle = 0),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_Analysis_Jusage_pools_bar.pdf",sep="")),  width = 10,  height = 3)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(v_usage_plot_combo, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()



  #vgene light

  v_usage_A<-reshape2::melt(100*table(pool_a$v_call)/sum(table(pool_a$v_call)))
  v_usage_A$category<-"Pool A"

  v_usage_B<-reshape2::melt(100*table(pool_b$v_call)/sum(table(pool_b$v_call)))
  v_usage_B$category<-"Pool B"

  add_to_A<-unique(as.character(v_usage_B$Var1[!(v_usage_B$Var1 %in% v_usage_A$Var1)]))
  add_to_B<-unique(as.character(v_usage_A$Var1[!(v_usage_A$Var1 %in% v_usage_B$Var1)]))

  pools_Vusage<-rbind(v_usage_A,v_usage_B)

  pools_Vusage_all<-rbind(pools_Vusage,rbind(data.frame(Var1=add_to_A,value=0,category="Pool A"),data.frame(Var1=add_to_B,value=0,category="Pool B")))

  ordered_v_genes<-levels(pools_Vusage_all$Var1)[order(levels(pools_Vusage_all$Var1))]

  pools_Vusage_all$Var1<-factor(pools_Vusage_all$Var1,levels=ordered_v_genes)

  pools_Vusage_all_combo<-pools_Vusage_all
  pools_Vusage_all_combo$category<-"Pool A+B"
  pools_Vusage_all_combo2<-as.data.frame(pools_Vusage_all_combo %>%
                                           group_by(Var1) %>%
                                           summarize(value=sum(value),category=category))

  pools_Vusage_all_combo2<-pools_Vusage_all_combo2[!duplicated(pools_Vusage_all_combo2),]
  pools_Vusage_all_combo2$value<-pools_Vusage_all_combo2$value/2


  grid_v_usage<-rbind(pools_Vusage_all,pools_Vusage_all_combo2)
  #

  v_usage_plot_combo<- ggplot2::ggplot(data=grid_v_usage, ggplot2::aes(Var1, value,fill=category))+
    ggplot2::geom_bar(stat="identity",position="dodge",colour='black') +
    ggplot2::facet_grid(.~category)+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=c("gray30","red","darkgrey"))+#3rd pos))+
    ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=5, angle = 90,vjust=0.5),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_Analysis_VLusage_pools_bar.pdf",sep="")),  width = 10,  height = 3)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(v_usage_plot_combo, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()






  #j_gene (Light)

  v_usage_A<-reshape2::melt(100*table(pool_a$j_call)/sum(table(pool_a$j_call)))
  v_usage_A$category<-"Pool A"

  v_usage_B<-reshape2::melt(100*table(pool_b$j_call)/sum(table(pool_b$j_call)))
  v_usage_B$category<-"Pool B"


  add_to_A<-unique(as.character(v_usage_B$Var1[!(v_usage_B$Var1 %in% v_usage_A$Var1)]))
  add_to_B<-unique(as.character(v_usage_A$Var1[!(v_usage_A$Var1 %in% v_usage_B$Var1)]))




  pools_Vusage<-rbind(v_usage_A,v_usage_B)


  pools_Vusage_all<-rbind(pools_Vusage,data.frame(Var1=add_to_A,value=0,category="Pool A"))

  ordered_v_genes<-levels(pools_Vusage_all$Var1)[order(levels(pools_Vusage_all$Var1))]

  pools_Vusage_all$Var1<-factor(pools_Vusage_all$Var1,levels=ordered_v_genes)

  pools_Vusage_all_combo<-pools_Vusage_all
  pools_Vusage_all_combo$category<-"Pool A+B"
  pools_Vusage_all_combo2<-as.data.frame(pools_Vusage_all_combo %>%
                                           group_by(Var1) %>%
                                           summarize(value=sum(value),category=category))


  pools_Vusage_all_combo2<-pools_Vusage_all_combo2[!duplicated(pools_Vusage_all_combo2),]
  pools_Vusage_all_combo2$value<-pools_Vusage_all_combo2$value/2


  grid_v_usage<-rbind(pools_Vusage_all,pools_Vusage_all_combo2)
  #

  v_usage_plot_combo<- ggplot2::ggplot(data=grid_v_usage, ggplot2::aes(Var1, value,fill=category))+
    ggplot2::geom_bar(stat="identity",position="dodge",colour='black') +
    ggplot2::facet_grid(.~category)+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=c("gray30","red","darkgrey"))+
    ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=6, angle =0),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_Analysis_JLusage_pools_bar.pdf",sep="")),  width = 10,  height = 3)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(v_usage_plot_combo, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()







  #length_cdrh3

  cdr3_length_usage_A<-reshape2::melt(100*table(nchar(pool_a$cdr3))/sum(table(nchar(pool_a$cdr3))))
  cdr3_length_usage_A$category<-"Pool A"

  cdr3_length_usage_B<-reshape2::melt(100*table(nchar(pool_b$cdr3))/sum(table(nchar(pool_b$cdr3))))
  cdr3_length_usage_B$category<-"Pool B"

  pools_cdr3_length<-rbind(cdr3_length_usage_A,cdr3_length_usage_B)
  missing_nbs<-min(pools_cdr3_length$Var1):max(pools_cdr3_length$Var1)


  add_to_A<-unique(missing_nbs[!(missing_nbs %in% cdr3_length_usage_A$Var1)])
  add_to_B<-unique(missing_nbs[!(missing_nbs %in% cdr3_length_usage_B$Var1)])


  pools_cdr3_len_all<-rbind(pools_cdr3_length,rbind(data.frame(Var1=add_to_A,value=0,category="Pool A"),data.frame(Var1=add_to_B,value=0,category="Pool B")))

  pools_cdr3_len_all$Var1<-factor(pools_cdr3_len_all$Var1,levels=missing_nbs)

  pools_cdr3_len_all_combo<-pools_cdr3_len_all
  pools_cdr3_len_all_combo$category<-"Pool A+B"
  pools_cdr3_len_all_combo2<-as.data.frame(pools_cdr3_len_all_combo %>%
                                             group_by(Var1) %>%
                                             summarize(value=sum(value),category=category))


  pools_cdr3_len_all_combo2<-pools_cdr3_len_all_combo2[!duplicated(pools_cdr3_len_all_combo2),]
  pools_cdr3_len_all_combo2$value<-pools_cdr3_len_all_combo2$value/2



  grid_cdr3_length<-rbind(pools_cdr3_len_all,pools_cdr3_len_all_combo2)
  #



  cdrh3_len<- ggplot2::ggplot(data=grid_cdr3_length, ggplot2::aes(Var1, value,fill=category))+
    ggplot2::geom_bar(stat="identity",position="dodge",colour='black') +
    ggplot2::facet_grid(.~category)+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=c("gray30","red","darkgrey"))+
    ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=7),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_Analysis_cdrh3_len_pools_bar.pdf",sep="")),  width = 10,  height = 3)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(cdrh3_len, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()


  #length_cdrl3

  cdr3_length_usage_A<-reshape2::melt(100*table(nchar(pool_a$junction_aa.light))/sum(table(nchar(pool_a$junction_aa.light))))
  cdr3_length_usage_A$category<-"Pool A"

  cdr3_length_usage_B<-reshape2::melt(100*table(nchar(pool_b$junction_aa.light))/sum(table(nchar(pool_b$junction_aa.light))))
  cdr3_length_usage_B$category<-"Pool B"

  pools_cdr3_length<-rbind(cdr3_length_usage_A,cdr3_length_usage_B)
  missing_nbs<-min(pools_cdr3_length$Var1):max(pools_cdr3_length$Var1)


  add_to_A<-unique(missing_nbs[!(missing_nbs %in% cdr3_length_usage_A$Var1)])
  add_to_B<-unique(missing_nbs[!(missing_nbs %in% cdr3_length_usage_B$Var1)])


  pools_cdr3_len_all<-rbind(pools_cdr3_length,rbind(data.frame(Var1=add_to_A,value=0,category="Pool A"),data.frame(Var1=add_to_B,value=0,category="Pool B")))

  pools_cdr3_len_all$Var1<-factor(pools_cdr3_len_all$Var1,levels=missing_nbs)

  pools_cdr3_len_all_combo<-pools_cdr3_len_all
  pools_cdr3_len_all_combo$category<-"Pool A+B"
  pools_cdr3_len_all_combo2<-as.data.frame(pools_cdr3_len_all_combo %>%
                                             group_by(Var1) %>%
                                             summarize(value=sum(value),category=category))


  pools_cdr3_len_all_combo2<-pools_cdr3_len_all_combo2[!duplicated(pools_cdr3_len_all_combo2),]
  pools_cdr3_len_all_combo2$value<-pools_cdr3_len_all_combo2$value/2



  grid_cdr3_length<-rbind(pools_cdr3_len_all,pools_cdr3_len_all_combo2)
  #



  cdrh3_len<- ggplot2::ggplot(data=grid_cdr3_length, ggplot2::aes(Var1, value,fill=category))+
    ggplot2::geom_bar(stat="identity",position="dodge",colour='black') +
    ggplot2::facet_grid(.~category)+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=c("gray30","red","darkgrey"))+
    ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=10),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_Analysis_cdrl3_len_pools_bar.pdf",sep="")),  width = 10,  height = 3)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(cdrh3_len, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()












  hl_pairing<-chosen_132[,c("v_gene","v_call.light","batch.1")]

  table(paste(hl_pairing$v_gene,hl_pairing$v_call.light,sep="-"))

  names(hl_pairing)<-c("v_call","v_call.light","category")

  hl_pairing<-hl_pairing[hl_pairing$v_call!="",]

  hl_usage_plot_combo<- ggplot2::ggplot(data=hl_pairing, ggplot2::aes(v_call, v_call.light,colour=category))+
    ggplot2::geom_point(position=ggplot2::position_jitter(h=0.2, w=0.1),alpha = 0.8, size = 3)+
    ggplot2::theme_bw()+
    ggplot2::scale_colour_manual(values=c("gray30","darkgrey"))+
    ggplot2::labs(x = "", y = "", fill = "", colour = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=9, angle = 90,vjust=0.5),axis.title.y = ggplot2::element_text(size=8),axis.text.y = ggplot2::element_text(size=8, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_Analysis_HL_Vusage_pool_bar_selected_only.pdf",sep="")),  width = 11,  height = 6)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(hl_usage_plot_combo, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()






  hl_pairing<-chosen_132[,c("j_gene","j_call.light","batch.1")]

  table(paste(hl_pairing$j_gene,hl_pairing$j_call.light,sep="-"))

  names(hl_pairing)<-c("j_call","j_call.light","category")

  hl_pairing<-hl_pairing[hl_pairing$j_call!="",]

  hl_usage_plot_combo<- ggplot2::ggplot(data=hl_pairing, ggplot2::aes(j_call, j_call.light,colour=category))+
    ggplot2::geom_point(position=ggplot2::position_jitter(h=0.2, w=0.1),alpha = 0.8, size = 3)+
    ggplot2::theme_bw()+
    ggplot2::scale_colour_manual(values=c("gray30","darkgrey"))+
    ggplot2::labs(x = "", y = "", fill = "", colour = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, angle = 0),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_Analysis_HL_Jusage_pools_bar_selected_only.pdf",sep="")),  width = 11,  height = 6)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(hl_usage_plot_combo, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()



  hl_pairing<-chosen_132[,c("v_gene","v_call.light","reactive")]

  table(paste(hl_pairing$v_call,hl_pairing$v_call.light,sep="-"))

  names(hl_pairing)<-c("v_call","v_call.light","category")
  hl_pairing_binders_only <- hl_pairing[hl_pairing$category!="",]

  vhl_usage_plot_combo<- ggplot2::ggplot(data=hl_pairing_binders_only, ggplot2::aes(v_call, v_call.light,colour=category))+
    ggplot2::geom_point(position=ggplot2::position_jitter(h=0.2, w=0.1),alpha = 0.6, size = 3)+
    ggplot2::theme_bw()+
    ggplot2::scale_colour_manual(values=c("tomato","red3","grey","black"))+#3rd pos))+
    ggplot2::labs(x = "", y = "", fill = "", colour = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, angle = 0),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=6, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_Analysis_HL_Vusage_pools_bar_binders_only.pdf",sep="")),  width = 11,  height = 6)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(vhl_usage_plot_combo, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()







  hl_pairing<-chosen_132[,c("j_gene","j_call.light","reactive")]

  names(hl_pairing)<-c("j_call","j_call.light","category")
  hl_pairing_binders_only <- hl_pairing[hl_pairing$category!="",]

  jhl_usage_plot_combo<- ggplot2::ggplot(data=hl_pairing_binders_only, ggplot2::aes(j_call, j_call.light,colour=category))+
    ggplot2::geom_point(position=ggplot2::position_jitter(h=0.2, w=0.1),alpha = 0.6, size = 3)+
    ggplot2::theme_bw()+
    ggplot2::scale_colour_manual(values=c("tomato","red3","grey","black"))+
    ggplot2::labs(x = "", y = "", fill = "", colour = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, angle = 0),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_Analysis_Supp_Pool_AnalysisHL_Jusage_pools_bar_binders_only.pdf",sep="")),  width = 11,  height = 6)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(jhl_usage_plot_combo, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()





  #compare V_identity 37 bin

  binders_37<-chosen_132[chosen_132$reactive!="",]


  v_id_h<-binders_37[,c("v_identity","reactive")]
  v_id_h$chain<-"Heavy"
  v_id_l<-binders_37[,c("v_identity.light","reactive")]
  v_id_l$chain<-"Light"
  names(v_id_l)<-c("v_identity","reactive","chain")

  v_id_37<-rbind(v_id_h,v_id_l)



  v_id_h_132<-chosen_132[,c("v_identity","reactive")]
  v_id_h_132$chain<-"Heavy"
  v_id_l_132<-chosen_132[,c("v_identity.light","reactive")]
  v_id_l_132$chain<-"Light"
  names(v_id_l_132)<-c("v_identity","reactive","chain")

  v_id_132<-rbind(v_id_h_132,v_id_l_132)
  v_id_132[,"reactive"]<-"Selected 132"


  v_id<-rbind(v_id_37,v_id_132)

  vid_plot_combo<- ggplot2::ggplot(data=v_id, ggplot2::aes(x=reactive,y=v_identity,fill=reactive))+
    ggplot2::geom_boxplot()+
    ggplot2::facet_grid(.~chain)+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=c("tomato","red3","grey","black","white"))+
    ggplot2::labs(x = "", y = "", fill = "", colour = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, angle = 0),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_V_identity_box_binders_vs_selected.pdf",sep="")),  width = 11,  height = 6)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(vid_plot_combo, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()



  example_data_scCovid<-bind_rows(example_data[1:16])

  example_data_scCovid$c_call_short<-substr(example_data_scCovid$c_call,1,4)
  unique(example_data_scCovid$c_call_short)

  example_data_scCovid_igag<-example_data_scCovid[example_data_scCovid$c_call_short %in% c("IGHA","IGHG"),]

  v_id_h_all<-example_data_scCovid_igag[,c("v_identity","dataset_id")]
  v_id_h_all$chain<-"Heavy"
  names(v_id_h_all)<-c("v_identity","reactive","chain")


  v_id_l_all<-example_data_scCovid_igag[,c("v_identity.light","dataset_id")]
  v_id_l_all$chain<-"Light"
  names(v_id_l_all)<-c("v_identity","reactive","chain")

  v_id_all<-rbind(v_id_h_all,v_id_l_all)
  v_id_all[,"reactive"]<-"All IGHG/A"

  v_id_incl_all<-rbind(v_id,v_id_all)

  unique(v_id_incl_all$reactive)

  v_id_incl_all$reactive<-factor(v_id_incl_all$reactive,levels=c("S1", "S1:RBD", "S1;S2","S2", "Selected 132", "All IGHG/A"))

  unique(v_id_incl_all$reactive)


  vid_plot_combo_all<- ggplot2::ggplot(data=v_id_incl_all, ggplot2::aes(x=reactive,y=v_identity,fill=reactive))+
    ggplot2::geom_boxplot(show.legend = FALSE)+
    ggplot2::facet_grid(.~chain)+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=c("tomato","red3","grey","black","ghostwhite","white"))+
    ggplot2::labs(x = "", y = "V-identity (%)", fill = "", colour = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=11, angle = 0),axis.title.y = ggplot2::element_text(size=14),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_V_identity_box_binders_vs_selected_vs_all.pdf",sep="")),  width = 13,  height = 6)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(vid_plot_combo_all, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()


  vid_plot_combo_all<- ggplot2::ggplot(data=v_id_incl_all, ggplot2::aes(x=reactive,y=v_identity,fill=reactive))+
    ggplot2::geom_boxplot(show.legend = FALSE)+
    ggplot2::geom_jitter(alpha=0.4,show.legend = FALSE)+
    ggplot2::facet_grid(.~chain)+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values=c("tomato","red3","grey","black","ghostwhite","white"))+
    ggplot2::labs(x = "", y = "V-identity (%)", fill = "", colour = "") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=7, angle = 0),axis.title.y = ggplot2::element_text(size=12),axis.text.y = ggplot2::element_text(size=12, angle = 0))

  grDevices::pdf(file.path(figure_path, paste("Supp_Pool_V_identity_box_binders_vs_selected_vs_all_jitter.pdf",sep="")),  width = 8,  height = 4)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  print(vid_plot_combo_all, vp=vplayout(1,1))   # plot for row 1, column 1
  dev.off()


  return(TRUE)



}


