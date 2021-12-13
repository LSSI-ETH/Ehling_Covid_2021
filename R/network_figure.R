#' Network plot of binders and closest found neighbors
#' @param binders_all Identified binders
#' @param binder_matches Closest neighbors to identified binders from repertoire data.
#' @param figure_path Path where figures are saved
#' @return TRUE (plots saved as pdfs into path_figures)
#' @examples
#' network_figure(binders_all,binder_matches,figure_path="")
#'

network_figure <- function(binders_all,binder_matches,figure_path=""){

  ### based on a threshold percentage (here 0.2 the full datasets are searched for matches within that distance resulting in bindermatches)
  ##thresholddistance set by percent of input length
  #threshold_distance<-0.2
  #
  #list_seqs_found<-list()
  #for(i in 1:nrow(binders_all)){
  #  curr_seq<-binders_all$junction_aa[i]
  #  curr_binds<-binders_all$reactive.to[i]
  #
  #  allowed_dist<-nchar(curr_seq)*threshold_distance
  #  all_dists<-stringdist::stringdist(curr_seq,all_seqs$junction_aa,method="hamming")
  #  seqs_within_dist<-which(all_dists<=allowed_dist)
  #
  #  all_seqs_within_dist<-all_seqs[seqs_within_dist,c("junction_aa","dataset_id","v_identity","v_call")]
  #  all_seqs_within_dist[!duplicated(all_seqs_within_dist),]
  #  #all_seqs_within_dist<-all_seqs_within_dist[all_seqs_within_dist!=curr_seq]
  #  all_seqs_within_dist$reactive.to<-curr_binds
  #  list_seqs_found[[i]]<-all_seqs_within_dist#c(curr_seq,all_seqs_within_dist)
  #  #binders_all$found_closest[i]<-paste(all_seqs_within_dist,collapse="|")
  #
  #}
  #
  #binder_matches<-bind_rows(list_seqs_found)
  #binder_matches<-binder_matches[!duplicated(binder_matches),]

  binder_matches_S1<-binder_matches[binder_matches$reactive.to=="S1",]
  binder_matches_S2<-binder_matches[binder_matches$reactive.to=="S2",]
  binder_matches_S1RBD<-binder_matches[binder_matches$reactive.to=="S1:RBD",]
  binder_matches_S1S2<-binder_matches[binder_matches$reactive.to=="S1;S2",]

  s1s2<-binder_matches[,c("junction_aa","reactive.to")]


  s1s2out<-s1s2 %>%
    group_by(junction_aa) %>%
    summarize(binding=paste(unique(reactive.to),collapse="_"))

  network_input<-s1s2out

  network_input[!(network_input$junction_aa %in% binders_all$junction_aa),"binding"]<-"Neighbor"


  ldistance<-3

  ld_matrix<-stringdist::stringdistmatrix(network_input$junction_aa,network_input$junction_aa,method="hamming")
  rownames(ld_matrix)<-unlist(network_input$junction_aa)
  colnames(ld_matrix)<-unlist(network_input$junction_aa)

  adj_matrix<-ld_matrix

  similarity <- ldistance
  adj_matrix[adj_matrix<=similarity]<-1
  adj_matrix[adj_matrix>similarity]<-0

  #construct graph
  graph_similarity<-igraph::graph_from_adjacency_matrix(adj_matrix,mode=c('undirected'),diag=F)

  order_graph_nodes<-names(igraph::V(graph_similarity))

  names_nodes_list<-list()
  for(i in 1:length(order_graph_nodes)){
    curr_node<-order_graph_nodes[i]
    if(curr_node %in% binders_all$junction_aa){
      curr_mab<-binders_all[binders_all$junction_aa==curr_node,"mAb"]
      names_nodes_list[[i]]<-curr_mab
    }else{
      names_nodes_list[[i]]<-""
    }
  }

  names_nodes<-unlist(names_nodes_list)

  node_cols=c(S1="tomato",S2="black",S1_S2="grey",S1RBD='red3',Neighbor="white")

  for(i in 1:length(order_graph_nodes)){
    curr_bindr<-network_input[network_input$junction_aa==order_graph_nodes[i],"binding"]
    igraph::V(graph_similarity)$color[i] <- ifelse(curr_bindr == "S1",node_cols["S1"],ifelse(curr_bindr == "S2",node_cols["S2"],ifelse(curr_bindr=="S1:RBD",node_cols["S1RBD"],ifelse(curr_bindr=="S1;S2",node_cols["S1_S2"],node_cols["Neighbor"]))))
    igraph::V(graph_similarity)$label[i]<-names_nodes[i]#order_graph_nodes[i]
  }

  cairo_pdf(file.path(figure_path,"graph_labeled.pdf"), family="Helvetica", width = 20,  height = 20,onefile=TRUE)
  par(mfrow=c(1,1))
  plot(graph_similarity, layout=igraph::layout_nicely,
       vertex.color=igraph::V(graph_similarity)$color,
       vertex.size=4.5,
       vertex.frame.width=0.02,
       vertex.label.cex=c(1.5), #original 1.5
       vertex.label=igraph::V(graph_similarity)$label,
       vertex.label.dist=1.1,
       edge.width=1.5
  )

  legend("bottomright",legend=c("S1","S1:RBD","S2","S1_S2","Closest Neighbor"),fill=c(node_cols["S1"],node_cols["S1RBD"],node_cols["S2"],node_cols["S1_S2"],node_cols["Neighbor"]), cex = 0.7)

  dev.off()

  return(TRUE)

}

