
#' Title
#'
#' @param node.labels
#' @param n.grns
#'
#' @return
#' @export
#'
#' @examples
randomly_select_modules<-function(node.labels, n.grns){
  nr.module<-length(unique(node.labels$module.greedy))
  modules<-unique(node.labels$module.greedy)
  step<-max(1, floor(nr.module/(n.grns+1)))
  print(step)
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
#' @param g
#' @param l
#' @param node.labels
#' @param nr.dis
#'
#' @return
#' @export
#'
#' @examples
randomly_select_disregulated_node<-function(g, l, node.labels, nr.dis){
  collect<-list()
  regs<-which(igraph::degree(g, mode = 'out')>0)
  disregulated_regulators<-c()
  for (i in 1:nrow(l)){
    dis<-l[i, module]
    print(dis)
    print(names(regs))
    eligible<-node.labels[node %in% names(regs)& module.greedy==dis]$node
    nr.dis.c<-min(nr.dis, length(eligible))
    nr.dis.c<-max(1, nr.dis.c)
    print(nr.dis.c)
    print(eligible)
    dn<-sample(eligible, size = nr.dis.c)
    for (d in dn){
      collect[[paste0(d,i)]]<-c(l[i,], d)
    }
  }
  collect<-rbindlist(collect)
  colnames(collect)<-c('module', 'grn', 'disregulated_gene')
  return(collect)

}



#' Title
#'
#' @param net
#' @param disregulated_regulators
#'
#' @return
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

#' Title
#'
#' @param net
#' @param disregulated_regulators
#' @param weight_delta
#'
#' @return
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

#' Title
#'
#' @param net
#' @param disregulated_regulators
#' @param weight_delta
#'
#' @return
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

#' Title
#'
#' @param net
#' @param disregulated_regulators
#'
#' @return
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

#' Title
#'
#' @param net
#' @param disregulated_regulators
#' @param tf_effect
#'
#' @return
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

#' Title
#'
#' @param net
#' @param n.cells
#' @param seed
#' @param tree
#' @param noise
#'
#' @return
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

#' Title
#'
#' @param counts
#' @param meta
#' @param directory
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


#' Title
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


