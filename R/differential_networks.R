
#' Assign modules to be disregulated to each GRN
#'
#' @param node.labels A table with the name of the nodes and the module they belong to
#' @param n.grns Number of GRNs/cell populations to produce
#'
#' @return A data.table containing modules and the grn they have been assigned to
#' @export
#'
#' @examples
randomly_select_modules<-function(node.labels, n.grns){
  nr.module<-length(unique(node.labels$module.greedy))
  modules<-unique(node.labels$module.greedy)
  step<-max(1, floor(nr.module/(n.grns+1)))
  groups<-c()
  label<-c()
  for (i in 0:(n.grns-1)){
    print((i*step+1))
    groups<-c(groups, modules[(i*step+1):((i+1)*step)])
    label<-c(label, rep(i,step))
  }
  return(data.table::data.table(module=groups, grn=label))
}



#' Title
#'
#' @param g The graph to which the GRNs should be added
#' @param l The list of modules to be perturberd
#' @param node.labels A data table of nodenames and their modules.
#' @param nr.dis Number of disregulated genes per module (will be rounded down if too many)
#'
#' @return A data frame containing the module, which GRN is assigned to to it and the name of the disregulated gene
#' @export
#'
#' @examples
randomly_select_disregulated_node<-function(g, l, node.labels, nr.dis){
  collect<-list()
  regs<-which(igraph::degree(g, mode = 'out')>0)
  disregulated_regulators<-c()
  for (i in 1:nrow(l)){
    dis<-l[i, module]
    eligible<-node.labels[node %in% names(regs)& module.greedy==dis]$node
    nr.dis.c<-min(nr.dis, length(eligible))
    nr.dis.c<-max(1, nr.dis.c)
    dn<-sample(eligible, size = nr.dis.c)
    for (d in dn){
      collect[[paste0(d,i)]]<-c(l[i,], d)
    }
  }
  collect<-rbindlist(collect)
  colnames(collect)<-c('module', 'grn', 'disregulated_gene')
  return(collect)

}



#' Add a base effect to all edges to the GRN
#' At the moment this samples from a random Gaussian with mean 2.5 and standard deviation 1
#'
#' @param net The network as an edge list
#' @param disregulated_regulators The data table containing the perturbation information
#'
#' @return The network with an additional column base.effect containing the base weights
#' @export
#'
#' @examples
add_base_effect_to_grn<-function(net, disregulated_regulators){
  net$base.effect<-rnorm(nrow(net), mean = 2.5, sd = 1)
  n.interactions<-nrow(net)
  net<-net[rep(1:n.interactions, length(unique(disregulated_regulators$grn)))]
  net$grn<-rep(unique(disregulated_regulators$grn), each = n.interactions)
  return(net)
}

#' Modifify the weights to induce the differential behaviour for the networks.
#' This function modifies the weight of the edges in the mapping.
#' this means one has to be careful when working with more than two grns, because
#' the effect is only added/substracted to the selected modules
#' Good for amping up the effect of the TF
#'
#' @param net The network with source, target and base.effect columns
#' @param disregulated_regulators The mapping for the disregulated modules
#' @param weight_delta The effect to be added to the selected edges
#'
#' @return net A network with an additional column grn.effect which is a copy of the base effect but modified according
#' to the network mapping
#' @export
#'
#' @examples
modify_weight_of_edges<-function(net, disregulated_regulators, weight_delta = 0.5){
  net$grn.effect<-net$base.effect
  for (gg in unique(disregulated_regulators$grn)){
    net[(grn == gg) & (source %in% disregulated_regulators[grn==gg]$disregulated_gene),]$grn.effect<-
      net[(grn == gg) & (source %in% disregulated_regulators[grn==gg]$disregulated_gene),]$grn.effect + weight_delta
  }
  return(net)
}

#' Modifies the selected TF in the GRNS which are NOT part of the selected module.
#' Good for reducing the effect of TFs in all but the selected module
#'
#' @param net The input network with columns source, target and base.effect
#' @param disregulated_regulators The mapping with the perturbed regulators
#' @param weight_delta The weight to be added/substracted
#'
#' @return net A network with an additional column grn.effect which is a copy of the base effect but modified according
#' to the network mapping
#' @export
#'
#' @examples
modify_weight_of_other_edges<-function(net, disregulated_regulators, weight_delta = 0.5){
  net$grn.effect<-net$base.effect
  for (gg in unique(disregulated_regulators$grn)){
    net[(grn != gg) & (source %in% disregulated_regulators[grn!=gg]$disregulated_gene),]$grn.effect<-
      net[(grn != gg) & (source %in% disregulated_regulators[grn!=gg]$disregulated_gene),]$grn.effect + weight_delta
  }
  return(net)
}

#' Set the effect of the Master TFs down for all modules expect the one where it is supposed to be
#' active.
#'
#' @param net The input network with columns source, target and base.effect
#' @param disregulated_regulators The mapping with the perturbed regulators
#'
#' @return net A network with an additional column grn.effect which is a copy of the base effect but modified according
#' to the network mapping
#' @export
#'
#' @examples
remove_master_tf_effect_from_other_modules<-function(net, disregulated_regulators){
  net$grn.effect<-net$base.effect
  for (gg in unique(disregulated_regulators$grn)){
    print(gg)
    net[(grn != gg) & (source %in% disregulated_regulators[grn!=gg]$disregulated_gene),]$grn.effect<-0.01
  }
  return(net)
}

#' Set the effect of the TFs to a specific value according to their effect.
#'
#' @param net The input network with columns source, target and base.effect
#' @param disregulated_regulators The mapping with the perturbed regulators
#' @param tf_effect The weight to set the edges to
#'
#' @return net A network with an additional column grn.effect which is a copy of the base effect but modified according
#' to the network mapping
#' @export
#'
#' @examples
set_master_tf_effect<-function(net, disregulated_regulators, tf_effect=5){
  net$grn.effect<-net$base.effect
  for (gg in unique(disregulated_regulators$grn)){
    net[(grn == gg) & (source %in% disregulated_regulators[grn==gg]$disregulated_gene),]$grn.effect<-tf_effect
  }
  return(net)
}

#' Use scMultiSim to generate the data with the master network
#'
#' @param net
#' @param n.cells Number of cells per GRN
#' @param seed Random seet
#' @param tree Differentiation tree required for scMultisim (please refer to documentation)
#' @param noise Whether to add typical noise to the count data
#'
#' @return List of two elements: count matrix, and associated metadata as a data.table
#' @export
#'
#' @examples
generate_data_from_grn<-function(net, n.cells = 750, seed=11, tree = scMultiSim::Phyla1(), noise = FALSE){

  count_list<-list()
  meta_list<-list()

  for (gg in unique(net$grn)){
    # Subset grn and reorder
    sub.grn<-net[grn == gg]
    sub.grn<-sub.grn[, c('target','source', 'grn.effect')]

    results <- scMultiSim::sim_true_counts(list(
      # required options
      rand.seed = seed,
      GRN = sub.grn,
      tree = tree,
      num.cells = n.cells,
      # optional options
      num.cif = 50,
      discrete.cif = TRUE,
      cif.sigma = 0.1,
      unregulated.gene.ratio  = 0.05))

    if (noise){
      scMultiSim::add_expr_noise(
        results,
        # options go here
        alpha_mean = 1e4
      )
    }
    counts<-results$counts
    # remove additional genes, to avoid spurious clustering
    counts<-counts[sapply(rownames(counts), function(x) !startsWith(x, 'gene')),]
    meta<-results$cell_meta
    meta$grn<-gg

    count_list[[(gg+1)]]<-counts
    meta_list[[(gg+1)]]<-meta
  }

  meta<-rbindlist(meta_list)
  counts<-do.call("cbind", count_list)

  return(list(counts, meta))
}

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


#' Generate data for a list of pregenerated networks.
#'
#' @param net.dir
#' @param gex.dir
#' @param network.list
#' @param nr.modules
#' @param nr.genes.per.module
#' @param disregulation.type
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
batch_create_gex_data<-function(net.dir,
                                gex.dir,
                                network.list,
                                nr.modules = 2,
                                nr.genes.per.module = 2,
                                disregulation.type= 'tf-effect-0',
                                seed=11){
  metadata.list<-list()
  for (net.name in network.list){

    # Get names of files
    net.path<-file.path(net.dir, net.name, 'edges.tsv')
    node.path <-file.path(net.dir, net.name, 'nodes.tsv')

    # read network data
    net<-data.table::fread(net.path, sep = '\t')
    node.labels <- data.table::fread(node.path, sep = '\t')

    # create graph
    g<-igraph::graph_from_data_frame(net, directed = TRUE, vertices = NULL)

    # select modules and genes to be perturbed
    l<-randomly_select_modules(node.labels, nr.modules)
    disregulated_regulators<-randomly_select_disregulated_node(g, l,node.labels, nr.genes.per.module)

    # Add gaussian weights as a base
    net<-add_base_effect_to_grn(net, disregulated_regulators)

    # Set effects of selected regulators t 0
    if (disregulation.type ==  'remove-tf-effect-others'){
      net<-remove_master_tf_effect_from_other_modules(net, disregulated_regulators)
    }else if (disregulation.type ==  'set-tf-value'){
      net<-set_master_tf_effect(net, disregulated_regulators, tf_effect = 5.5)
    } else if (disregulation.type == 'remove-weight-from-tf'){
      net<- modify_weight_of_edges(net, disregulated_regulators , weight_delta = -0.5)
    }else if(disregulated.type == 'remove-weight-from-others'){
      net<- modify_weight_of_other_edges(net, disregulated_regulators , weight_delta = -0.5)
    }
    else
    {
      net<-remove_master_tf_effect_from_other_modules(net, disregulated_regulators)
    }

    # Generate data
    generated_data <- generate_data_from_grn(net, seed = seed)

    counts <- generated_data[[1]]
    meta<- generated_data[[2]]

    # save the data, plots and files in a directory
    current_experiment_id<- glue::glue("{net.name}_nm_{nr.modules}_gpm_{nr.genes.per.module}_{disregulation.type}_{seed}")
    gex.current.dir<-file.path(gex.dir, net.name, current_experiment_id)
    gex.plot<-file.path(gex.current.dir, 'plots')
    dir.create(gex.plot, recursive = TRUE)

    # save plots before transforming data into data.table
    save_plots(counts, meta, gex.plot)

    # save counts and metadata
    gene_names<-rownames(counts)
    counts<-data.table::as.data.table(t(counts))
    colnames(counts)<-gene_names
    data.table::fwrite(counts, file = file.path(gex.current.dir, 'gex.tsv'), sep = '\t')
    data.table::fwrite(meta, file = file.path(gex.current.dir, 'metadata.tsv'), sep = '\t')

    # save network
    data.table::fwrite(net, file = file.path(gex.current.dir, 'net.tsv'), sep = '\t')

    metadata.list[[current_experiment_id]]<-list(current_experiment_id, net.name, nr.modules, nr.genes.per.module, disregulation.type, seed, gex.plot, file.path(gex.current.dir, 'gex.tsv'), file.path(gex.current.dir, 'metadata.tsv'))
  }
  metadata.list<-rbindlist(metadata.list)
  colnames(metadata.list)<-c('experiment_id', 'base.network.name', 'nr.grns', 'nr.genes.per.grn', 'disregulation.type', 'seed', 'plot.dir', 'gex.file', 'metadata.file')
  data.table::fwrite(metadata.list, file = file.path(gex.dir, 'metadata.tsv'), sep = '\t')
  return(metadata.list)
}


