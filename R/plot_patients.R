#' Sideways barplot for patient timeline
#' @param patients Patient Metadata.
#' @param path_figures Path to which figures should be saved.
#' @return TRUE (plots saved as pdfs into path_figures)
#' @examples
#' plot_patients(patients=patients_metadata,path_figures="./figures")

plot_patients<-function(patients=patients_metadata,figure_path="Ehling_figs"){

  #Curate data for overview Figure
  patients<-patients[!is.na(patients$record_id),]
  patients$name<-paste(patients$fid,", age: ",patients$age,sep="")
  patients_name<-patients[order(patients$age),"name"]
  patients$name<-factor(patients$name,levels=patients_name)

  patients$symptoms_to_blood_collect<-patients$blood_collect-patients$days_ill_total
  max_y_axis<-max(patients$blood_collect)


  patients_ill<-patients[,c("name","days_ill_total","sex","Total.Sympt")]
  patients_blood<-patients[,c("name","symptoms_to_blood_collect","sex","Total.Sympt")]

  patients_ill$category<-"days_ill"
  patients_blood$category<-"days_to_blood_collect"

  names(patients_ill)<-c("name", "days", "sex", "Total.Sympt", "category")
  names(patients_blood)<-c("name", "days", "sex", "Total.Sympt", "category")
  patients_blood$sex<-"blood_collect"

  days_df<-rbind(patients_ill,patients_blood)

  days_df$sex<-factor(days_df$sex,levels=c("blood_collect","f","m"))

  #ggplot2 print of barplot
  fig1 <- ggplot2::ggplot(data=days_df, ggplot2::aes(name, days,fill=sex,label=Total.Sympt))+
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme_bw()+
    ggplot2::coord_flip()+
    ggplot2::scale_fill_manual(values=c("grey",'tomato','skyblue3'))+
    ggplot2::scale_y_continuous(breaks=seq(1:max_y_axis))+
    ggplot2::labs(x = "", y = "Days ill", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, angle = 0,vjust=0.5))+
    ggplot2::theme(legend.position="top", legend.box = "horizontal")

  grDevices::pdf(file.path(figure_path,"fig1_out.pdf"), width = 7.5,  height = 5)
  grid::pushViewport( grid::viewport(layout=grid::grid.layout(1,1, heights = grid::unit(c(0.5, 4),"null")))) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(fig1, vp=vplayout(1,1))
  grDevices::dev.off()##

  return(TRUE)
}
