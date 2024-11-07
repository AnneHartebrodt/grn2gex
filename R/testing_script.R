# #require(grn2gex)
# #require(scMultiSim)
# #require(data.table)
#
# require(yaml)
#
# ty = yaml.load_file('/home/bionets-og86asub/Documents/netmap/thenetmap/dockerize-grn2gex/config.yaml')
#
# net.dir<-'/home/bionets-og86asub/Documents/netmap/thenetmap//NetMap_LRP/data/simulation/collectri_subnetworks_2'
# collectri.file <- '/home/bionets-og86asub/Documents/netmap/thenetmap/NetMap_LRP/data/public/collecTRI2/collectri_2024-11-04.tsv'
# gex.dir<-'/home/bionets-og86asub/Documents/netmap/thenetmap/NetMap_LRP/data/simulation/collectri_simulated_data2'
#
# dir.create(gex.dir, recursive = TRUE)
#
# nr.modules <- 2
# nr.genes.per.module<-2
# disregulation.type <- 'remove-tf-effect-others'
# seed <- 11
#
#
# collectri <- loadOrDownloadCollectTRI(collectri.file)
# colnames(collectri)<- c('source', 'target', 'sad')
# graph_list<-clusterNetwork(collectri )
#
# graph_list<-remove_irrelevant_networks(graph_list)
# orig_graph<-createGraph(collectri)
# for (g in graph_list){
#   save_small_subgraph(g, net.dir, orig_graph)
# }
#
# res<-restore_original_directionality(graph_list[[1]], orig_graph)
#
# dataset<-create_gex_data(res$egdelist, res$nodelist, 'test')
#
# save_generated_data(dataset$net, dataset$counts, dataset$meta, dataset$info, gex_dir = gex.dir)
