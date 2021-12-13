#' Sideways barplot
#' @param list_repertoires Patient Metadata.
#' @param order_ids Path to which figures should be saved.
#' @param chain Path to which figures should be saved.
#' @param sample_name Path to which figures should be saved.
#' @param text_size Path to which figures should be saved.
#' @param width_height Path to which figures should be saved.
#' @param plot_directory Path to which figures should be saved.
#' @param color_scale_category Path to which figures should be saved.
#' @return TRUE (plots saved as pdfs into path_figures)
#' @examples
#' pie_chart_isotype_usage(list_repertoires,order_ids,chain="H",sample_name)

pie_chart_isotype_usage<-function (list_repertoires, order_ids,chain = "H", sample_name, text_size = 5,
                              width_height = c(x = 15, y = 15), figure_path = "Ehling_figs",color_scale_category=c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'))
{
  pie_chart_df_plot_list <- list()
  for (i in 1:length(list_repertoires)) {
    c_call_curr <- list_repertoires[[i]][!duplicated(list_repertoires[[i]]),]
    if (chain == "H") {
      isotype_occ <- sapply(1:nrow(c_call_curr), function(x) substr(as.character(c_call_curr$c_call[x]),1, 4))
      isotype_occ[isotype_occ == "Miss"]<-"Missing"
      isotype_occ[isotype_occ == ""] <- "Missing"
      isotype_occ <- isotype_occ[isotype_occ %in% c("IGHA", "IGHD", "IGHG", "IGHM", "IGHE", "Missing")]
    }
    else if (chain == "L") {
      isotype_occ <- sapply(1:nrow(c_call_curr), function(x) substr(as.character(c_call_curr$c_call.light[x]),
                                                                    1, 4))
      isotype_occ[isotype_occ == "Miss"]<-"Missing"
      isotype_occ[isotype_occ == ""] <- "Missing"
      isotype_occ <- isotype_occ[isotype_occ %in% c("IGKC", "IGLC", "Missing")]
    }
    isotype_df_h <- reshape2::melt(100 * table(isotype_occ)/sum(table(isotype_occ)))
    isotype_df_h$class <- names(list_repertoires)[i]
    names(isotype_df_h) <- c("isotype", "Freq", "class")
    pie_chart_df <- isotype_df_h
    if (chain == "H") {
      pie_chart_df$category <- factor(pie_chart_df$isotype,
                                      levels = c("IGHA", "IGHD", "IGHE", "IGHG", "IGHM",
                                                 "Missing"))
      color_scale <- color_scale_category

      names(color_scale) <- c("IGHA", "IGHD", "IGHG", "IGHM",
                              "IGHE", "Missing")
    }
    else if (chain == "L") {
      pie_chart_df$category <- factor(pie_chart_df$isotype,
                                      levels = c("IGKC", "IGLC", "Missing"))
      color_scale <- c("#8da0cb", "#fc8d62", "grey")
      names(color_scale) <- c("IGKC", "IGLC", "Missing")
    }
    pie_chart_df_plot_tmp <- pie_chart_df %>% arrange(desc(category)) %>%
      mutate(prop = Freq) %>% mutate(ypos = cumsum(prop) -
                                       0.5 * prop)
    pie_chart_df_plot_list[[i]] <- pie_chart_df_plot_tmp
  }
  pie_chart_df_plot <- do.call(rbind, pie_chart_df_plot_list)

  pie_chart_df_plot$class<-factor(pie_chart_df_plot$class,levels=order_ids)

  pie_chart <- ggplot2::ggplot(data = pie_chart_df_plot, ggplot2::aes(x = "", y = prop, fill = category)) +
    ggplot2::geom_bar(width = 1, stat = "identity", colour = "black") +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::facet_wrap(~class) +
    ggplot2::scale_fill_manual(values = color_scale) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "", fill = "Occurrence", colour = "") +
    ggplot2::geom_text(data = pie_chart_df_plot,
                       ggplot2::aes(y = ypos, label = round(Freq, 2)), color = "black", size = text_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(2), hjust = 0), plot.margin = ggplot2::unit(c(2, 2, 2, 2), "points"), axis.text.x = ggplot2::element_text(size = 9, angle = 0), axis.text.y = ggplot2::element_text(size = 14), axis.title = ggplot2::element_text(size = 16, vjust = 0.4)) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(legend.position = "right", legend.key = ggplot2::element_blank())
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"), strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(4)), panel.grid.major.x = ggplot2::element_blank())

  grDevices::pdf(file.path(figure_path, paste("Isotype_usage_pie_", sample_name, ".pdf", sep = "")), width = width_height[["x"]], height = width_height[["y"]])
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 1), height = 1, width = 1))
  vplayout <- function(x, y) grid::viewport(layout.pos.row = x,
                                            layout.pos.col = y)
  print(pie_chart, vp = vplayout(1, 1))
  dev.off()

  return(TRUE)
}
