
#' Save diagnostic plots with the data
#'
#' @param counts count matrix!! (not data.table)
#' @param meta metadata data.table
#' @param directory directory where to save the plots to
#'
#' @return
#' @export
#'
#' @examples
save_plots<-function(counts, meta, directory){

  colors  = viridis::inferno(length(unique(meta$grn)))
  names(colors) = unique(meta$grn)

  pgrn<-scMultiSim::plot_tsne(counts, labels = meta$grn)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = 'bottom')+
    ggplot2::ylab('tSNE2')+
    ggplot2::xlab('tSNE1')+
    ggplot2::scale_color_manual(values = colors)

  pgrn

  ggplot2::ggsave(pgrn, file = file.path(directory, 'tsne_plot.pdf'), height = 15, width = 15, unit = 'cm')

  ha = ComplexHeatmap::HeatmapAnnotation(
    simulation = ComplexHeatmap::anno_simple(meta$grn, col=colors),
    cell = ComplexHeatmap::anno_simple(meta$pop),
    annotation_name_side = "left"
  )

  pdf(file.path(directory, 'expression_heatmap.pdf'))
  plot(ComplexHeatmap::Heatmap(counts, row_names_gp =grid::gpar(fontsize = 5),
                               top_annotation = ha, show_column_names = FALSE ))
  dev.off()


}
