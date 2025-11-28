# require(grn2gex)
# require(scMultiSim)
# require(data.table)
#
# library(ComplexHeatmap)
# library(circlize)
#
# # # #
# # # # require(yaml) # ty = yaml.load_file('/home/bionets-og86asub/Documents/netmap/thenetmap/dockerize-grn2gex/config.yaml')
# #
# net.dir<-'/home/bionets-og86asub/Documents/netmap/thenetmap//NetMap_LRP/data/simulation/collectri_subnetworks_2'
# collectri.file <-'/home/bionets-og86asub/Documents/netmap/thenetmap/NetMap_LRP/data/collectri.tsv'
# gex.dir<-'/home/bionets-og86asub/Documents/netmap/thenetmap/NetMap_LRP/data/simulation/collectri_simulated_data2'
#
# dir.create(gex.dir, recursive = TRUE)
#
# nr.modules <- 1
# nr.genes.per.module<-3
# disregulation.type <- 'adjust-weight-of-others'
# seed <- 11 #
#
# collectri <- loadOrDownloadCollectTRI(collectri.file)
# colnames(collectri)<- c('source', 'target', 'sad')
# graph_list<-clusterNetwork(collectri, min_nodes = 200, max_nodes = 500 )
#
# graph_list<-remove_irrelevant_networks(graph_list)
# orig_graph<-createGraph(collectri)
# for (g in graph_list){
#   save_small_subgraph(g, net.dir, orig_graph)
#  }
#
# res<-restore_original_directionality(graph_list[[4]], orig_graph)
# res2<-restore_original_directionality(graph_list[[11]], orig_graph)
# res3<-restore_original_directionality(graph_list[[5]], orig_graph)
# res4<-restore_original_directionality(graph_list[[2]], orig_graph)
# res5<-restore_original_directionality(graph_list[[3]], orig_graph)
#
# el1 <- rbind(res$egdelist, res2$egdelist, res3$egdelist, res4$egdelist)
# ll <- nrow(res$egdelist)
# ll2 <- nrow(res2$egdelist)
# ll3<- nrow(res3$egdelist)
# ll4<- nrow(res4$egdelist)
# el1$module<-c(rep(1, ll), rep(2, ll2), rep(3, ll3), rep(4, ll4))
#
#
# dataset<-create_gex_data_easy(el1, res5$egdelist, 'test', base_effect = 'standard-normal', mean = 4, sd=1)
#
# el1<-as.data.table(el)
#
# dataset<-create_gex_data(res$egdelist, res$nodelist, 'test', base_effect = 'standard-normal', mean = 4, sd=1, tf_effect ='adjust-weight-of-tf', weight_delta = 0.5, nr_grns = 3 )
#
# pgrn<-scMultiSim::plot_tsne(dataset$counts, labels = dataset$meta$grn)+ ggplot2::theme_bw()+ ggplot2::theme(legend.position = 'bottom')+
# ggplot2::ylab('tSNE2')+ ggplot2::xlab('tSNE1')
# pgrn
#
# scMultiSim::plot_gene_module_cor_heatmap(dataset$counts)
#
# plot(res$clustering, gg , vertex.size = 5, edge.arrow.size = 1, vertex.label.cex = 1)
#
#
# ha = ComplexHeatmap::HeatmapAnnotation(
#   annotation_name_side = "left"
# )
#
#
# plot(ComplexHeatmap::Heatmap(dataset$counts, row_names_gp =grid::gpar(fontsize = 5),
#                               show_column_names = FALSE ))
# gex <-dataset$counts
# # Compute the gene-gene coexpression matrix (correlation matrix)
# coexpression_matrix <- cor(t(gex), method = "pearson")
#
# # Define a color gradient for the heatmap
#
# col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
#
#
# # Create the coexpression heatmap
# heatmap <- ComplexHeatmap::Heatmap(
#   coexpression_matrix,
#   name = "Coexpression",
#   cluster_rows = TRUE,        # Cluster genes
#   cluster_columns = TRUE,     # Cluster genes
#   show_row_names = TRUE,
#   show_column_names = TRUE,
#   column_title = "Gene Coexpression",
#   row_title = "Gene Coexpression",
#   heatmap_legend_param = list(title = "Correlation")
# )
#
# # Draw the heatmap
# plot(heatmap)
#
