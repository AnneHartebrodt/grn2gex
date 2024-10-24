require(decoupleR)
require(OmnipathR)
require(igraph)
require(data.table)

###
# Function saves a small subgraph with a final round of clustering
###
save_small_subgraph<-function(sub, network.dir, orig.graph){


  # add number of nodes to filename
  counter <-vcount(sub)
  ## Add timestamp to avoid overwriting same sized graphs
  timestamp = as.numeric(format(Sys.time(), "%OS3")) * 1000

  local.dir<-file.path(network.dir, paste0('net_', counter, '_',timestamp))
  dir.create(local.dir, recursive = TRUE)
  # Create filename
  file.plot <- file.path(local.dir, paste0('net.pdf'))
  file.edge <- file.path(local.dir, paste0('edges.tsv'))
  file.node<-file.path(local.dir, paste0('nodes.tsv'))


  # Layout and create submodules
  lay <- layout_with_kk(sub)
  go<-cluster_fast_greedy(sub)


  sub<-subgraph(orig.graph, V(sub)$name)
  # Plot network
  pdf(file = file.plot, width = 600, height = 600 )
  plot(go, sub , vertex.size = 5, edge.arrow.size = 1, vertex.label.cex = 1)
  dev.off()

  # save edge list
  edgelist<-as.data.table(as_edgelist(sub) )
  colnames(edgelist)<-c('source', 'target')
  fwrite(edgelist, file = file.edge , sep = '\t')

  # save nodelist with community
  nodelist<-data.table(node=go$names, module.greedy=go$membership)
  fwrite(nodelist, file = file.node, sep = '\t')
}


cluster<-function(gr, network.dir, orig.graph){
  # recursive function subdividing the network
  # until a size of less than 200 nodes is reached
  # At 200 or fewer nodes the network is clustered one last
  # time to obtain modules and saved.
  groups<-cluster_fast_greedy(gr)
  group_sizes<-table(groups$membership)
  for (cl in names(group_sizes)){
    gs <- group_sizes[cl]
    sub<-subgraph(gr, which(groups$membership==cl))
    if (gs>150){
      cluster(sub, network.dir, orig.graph)
    }
    else if (gs<200 & gs>50){

      save_small_subgraph(sub, network.dir, orig.graph)
    }
  }
}


clusterCollecTRI<-function(net.dir, collectri_file){
  # Main function to cluster collectri network.
  # Check if collectri is stored on disk
  # otherwise download from web
  if (!file.exists(collectri_file)){
    # get network
    collectri <- get_collectri(organism='human', split_complexes=FALSE)
    fwrite(collectri, file = collectri_file , sep = '\t')
  }
  else{
    collectri = fread(collectri_file, sep = '\t')
  }
  ## make undirected graph (required for greedy clustering)
  collectri$mor<-1
  g<-graph_from_data_frame(collectri, directed = FALSE, vertices = NULL)
  orig.graph <- graph_from_data_frame(collectri, directed = TRUE, vertices = NULL)
  ## delete multiedges (required for greedy clustering)
  m<-which(which_multiple(g))
  g<-delete_edges(g, m)


  # Create network directory
  dir.create(net.dir, recursive = T)

  # call recursive clustering function.
  cluster(g, net.dir, orig.graph)
}

## Remove networks which contain only one or 2 hubnodes
remove_irrelevant_networks<-function(basedir){
  edgelists<-list.files(basedir, pattern = '.*edges.tsv', recursive = T)
  hedgehogs<-c()
  betweenness_list<-list()

  for( e in edgelists){
    net<-fread(file.path(basedir, e))
    g<-graph_from_data_frame(net, directed = FALSE, vertices = NULL)
    bet<- betweenness(g)
    #print(sort(bet))

    if (length(which(bet>0))<3){
      # remove networks with one gene regulating all the others
      hedgehogs<-c(hedgehogs, e)
    }
    else{
      betweenness_list[[e]]<-data.table(betweeness = bet, name = rep(e, length(bet)))
    }
  }

  edgelists<-setdiff(edgelists, hedgehogs)
  netlist<-data.table(net.filename<-c(edgelists, hedgehogs), type = c(rep('modularized', length(edgelists)), rep('hedgehog', length(hedgehogs))))
  fwrite(netlist, file = file.path(basedir, 'network_stats.tsv'), sep = '\t')

  return(edgelists)
}

color<-function(x, n = 3){
  pal<-viridis::inferno(n)
  return(pal[x])
}

## Compute spectral clustering
cluster_and_plot<-function(edgelist, basedir, n.clusters = 3){
  for (e in edgelist){
    net<-fread(file.path(basedir, e))
    g<-graph_from_data_frame(net, directed = FALSE, vertices = NULL)
    #plot(g, vertex.size = 5, edge.arrow.size = 1, vertex.label.cex = 1)

    ## Spectral clustering
    clu <- embed_laplacian_matrix(g, no=n.clusters, type = 'DAD', which = 'la')
    km <- kmeans(clu[["X"]], centers =n.clusters)

    nodes.file <- file.path(basedir, gsub('edges.tsv', 'nodes.tsv', e))
    nodes<-fread(nodes.file, sep = '\t')

    columnname<-paste0('cluster.spectral.', n.clusters)
    node.clu<-data.table(V(g)$name,km$cluster)
    colnames(node.clu)<-c('name', columnname)

    nodes<-merge(nodes, node.clu, by.x = 'node', by.y = 'name')
    fwrite(nodes, file = nodes.file, sep = '\t')

    V(g)$color<-sapply(km$cluster, function(x) color(x,n.clusters))

    file.plot <- file.path(basedir,  gsub('edges.tsv', paste0('spectral_', n.clusters, '.pdf'), e))
    pdf(file = file.plot, width = 600, height = 600 )
    plot(g, vertex.size = 5, edge.arrow.size = 1, vertex.label.cex = 1)
    dev.off()
  }
}



