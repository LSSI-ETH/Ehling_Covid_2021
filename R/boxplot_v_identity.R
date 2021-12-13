#' Boxplot for comparison of V_identity across samples
#' @param list_repertoires List of repertoires
#' @param rep_category metadata.
#' @return TRUE (plots saved as pdfs into path_figures)
#' @examples
#' boxplot_v_identity(list_repertoires=example_data,rep_category=sample_overview, path_figures="")

boxplot_v_identity<-function (list_repertoires, rep_category, chain = "H", color_scale = c("#d7191c",
                                                                                                    "#fdae61", "#2c7bb6", "#ffffbf", "#abd9e9"), plot_name = "",figure_path = ""){
  list_v_identity <- list()
  for (i in 1:length(list_repertoires)) {

    repertoire_name <- names(list_repertoires)[i]
    if (chain == "H") {
      ids_curr <- list_repertoires[[i]]$v_identity
      names(color_scale) <- c("IGHA","IGHD","IGHG","IGHM","IGHE","Missing")
    }
    else if (chain == "L") {
      ids_curr <- list_repertoires[[i]]$v_identity.light

    }
    v_ids <- data.frame(v_identity = ids_curr, c_call = list_repertoires[[i]]$c_call,
                        repertoire = repertoire_name)
    list_v_identity[[repertoire_name]] <- v_ids
  }
  v_id_df <- do.call(rbind, list_v_identity)
  rep_cat_min <- rep_category[rep_category$repertoire %in%
                                v_id_df$repertoire, ]
  v_id_df_cat <- merge(v_id_df, sample_overview, by = "repertoire")
  v_id_df_cat$c_call <- substr(v_id_df_cat$c_call, 1, 4)
  v_id_df_cat$c_call[v_id_df_cat$c_call=="Miss"]<-"Missing"
  v_id_df_cat$c_call[v_id_df_cat$c_call==""]<-"Missing"

  v_id_df_cat$v_identity<-100*v_id_df_cat$v_identity

  v_id_df_cat$c_call<-factor(v_id_df_cat$c_call,levels=c("IGHA","IGHD","IGHG","IGHM","IGHE","Missing"))

  v_id_df_stat<-as.data.frame(v_id_df_cat %>%
                                group_by(category,c_call) %>%
                                summarize(mean_freq=mean(v_identity),sd=sd(v_identity),se=plotrix::std.error(v_identity)))

  v_identity_plot <- ggplot2::ggplot(data = v_id_df_cat, ggplot2::aes(category, v_identity, fill = c_call)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = color_scale) +
    ggplot2::labs(x = "", y = "V-identity (%)", fill = "", colour = "Repertoire") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16, angle = 0, vjust = 0),axis.title.x = ggplot2::element_text(size=16),axis.text.y = ggplot2::element_text(size = 16, angle = 0),axis.title.y = ggplot2::element_text(size = 16, angle = 90))

  grDevices::pdf(file.path(figure_path, paste("v_identity_", plot_name, ".pdf", sep = "")), width = 5, height = 5)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 1, heights = grid::unit(c(0.5, 4), "null"))))
  vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y)
  base::print(v_identity_plot, vp = vplayout(1, 1))
  grDevices::dev.off()
}
