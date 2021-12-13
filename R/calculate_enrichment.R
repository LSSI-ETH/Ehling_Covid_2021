
.calculate_enrichment<-function(rep_1_name,rep_2_name){

  rep_1<-enrichment_data[[rep_1_name]][c("cloneCount","cloneFraction","nSeqCDR3","aaSeqCDR3")]

  rep_1 <- as.data.frame(rep_1 %>%
                           dplyr::group_by(aaSeqCDR3) %>%
                           dplyr::summarise(cloneCount=sum(cloneCount),cloneFraction=sum(cloneFraction),nSeqCDR3=paste(nSeqCDR3,collapse="|")))


  rep_2<-enrichment_data[[rep_2_name]][c("cloneCount","cloneFraction","nSeqCDR3","aaSeqCDR3")]

  rep_2 <- as.data.frame(rep_2 %>%
                           dplyr::group_by(aaSeqCDR3) %>%
                           dplyr::summarise(cloneCount=sum(cloneCount),cloneFraction=sum(cloneFraction),nSeqCDR3=paste(nSeqCDR3,collapse="|")))



  names(rep_1)<-c("aaSeqCDR3","cloneCount_rep1","cloneFraction_rep1","nSeqCDR3_rep1")
  names(rep_2)<-c("aaSeqCDR3","cloneCount_rep2","cloneFraction_rep2","nSeqCDR3_rep2")

  rep1_vs_rep2<-merge(rep_1,rep_2,by="aaSeqCDR3")


  rep1_vs_rep2$enrichment<-log2(rep1_vs_rep2$cloneFraction_rep1/rep1_vs_rep2$cloneFraction_rep2)
  rep1_vs_rep2<-rep1_vs_rep2[order(-rep1_vs_rep2$enrichment),]

  rep1_vs_rep2<-rep1_vs_rep2[!grepl("\\_",rep1_vs_rep2$aaSeqCDR3) & !grepl("\\*",rep1_vs_rep2$aaSeqCDR3), ]

  return(rep1_vs_rep2)
}
