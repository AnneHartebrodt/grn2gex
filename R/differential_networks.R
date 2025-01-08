#' Add a base effect to all edges to the GRN
#' At the moment this samples from a random Gaussian with mean 2.5 and standard deviation 1
#'
#' @param net The network as an edge list
#' @param disregulated_regulators The data table containing the perturbation information
#'
#' @return The network with an additional column base_effect containing the base weights
#' @export
#'
#' @examples
add_base_effect_gaussian<-function(net, disregulated_regulators, mean = 2.5, sd = 1.0){
  net$base_effect<-rnorm(nrow(net), mean = mean, sd = sd)
  n_interactions<-nrow(net)
  net<-net[rep(1:n_interactions, length(unique(disregulated_regulators$grn)))]
  net$grn<-rep(unique(disregulated_regulators$grn), each = n_interactions)
  return(net)
}

#' Modifify the weights to induce the differential behaviour for the networks.
#' This function modifies the weight of the edges in the mapping.
#' this means one has to be careful when working with more than two grns, because
#' the effect is only added/substracted to the selected modules
#' Good for amping up the effect of the TF
#'
#' @param net The network with source, target and base_effect columns
#' @param disregulated_regulators The mapping for the disregulated modules
#' @param weight_delta The effect to be added to the selected edges
#'
#' @return net A network with an additional column grn_effect which is a copy of the base effect but modified according
#' to the network mapping
#' @export
#'
#' @examples
modify_weight_of_edges<-function(net, disregulated_regulators, weight_delta = 0.5){
  net$grn_effect<-net$base_effect
  for (gg in unique(disregulated_regulators$grn)){
    net[(grn == gg) & (source %in% disregulated_regulators[grn==gg]$disregulated_gene),]$grn_effect<-
      net[(grn == gg) & (source %in% disregulated_regulators[grn==gg]$disregulated_gene),]$grn_effect + weight_delta
  }
  return(net)
}

#' Modifies the selected TF in the GRNS which are NOT part of the selected module.
#' Good for reducing the effect of TFs in all but the selected module
#'
#' @param net The input network with columns source, target and base_effect
#' @param disregulated_regulators The mapping with the perturbed regulators
#' @param weight_delta The weight to be added/substracted
#'
#' @return net A network with an additional column grn_effect which is a copy of the base effect but modified according
#' to the network mapping
#' @export
#'
#' @examples
modify_weight_of_other_edges<-function(net, disregulated_regulators, weight_delta = 0.5){
  net$grn_effect<-net$base_effect
  for (gg in unique(disregulated_regulators$grn)){
    net[(grn != gg) & (source %in% disregulated_regulators[grn!=gg]$disregulated_gene),]$grn_effect<-
      net[(grn != gg) & (source %in% disregulated_regulators[grn!=gg]$disregulated_gene),]$grn_effect + weight_delta
  }
  return(net)
}

#' Set the effect of the Master TFs down for all modules expect the one where it is supposed to be
#' active.
#'
#' @param net The input network with columns source, target and base_effect
#' @param disregulated_regulators The mapping with the perturbed regulators
#'
#' @return net A network with an additional column grn_effect which is a copy of the base effect but modified according
#' to the network mapping
#' @export
#'
#' @examples
remove_master_tf_effect_from_other_modules<-function(net, disregulated_regulators){
  net$grn_effect<-net$base_effect
  for (gg in unique(disregulated_regulators$grn)){
    net[(grn != gg) & (source %in% disregulated_regulators[grn==gg]$disregulated_gene),]$grn_effect<-0.01
  }
  return(net)
}

#' Set the effect of the TFs to a specific value according to their effect.
#'
#' @param net The input network with columns source, target and base_effect
#' @param disregulated_regulators The mapping with the perturbed regulators
#' @param tf_effect The weight to set the edges to
#'
#' @return net A network with an additional column grn_effect which is a copy of the base effect but modified according
#' to the network mapping
#' @export
#'
#' @examples
set_master_tf_effect<-function(net, disregulated_regulators, tf_effect=5){
  if(!('base_effect' %in% colnames(net))) BBmisc::stopf('net does not have a column "base effect", please initialize first')

  net$grn_effect<-net$base_effect
  for (gg in unique(disregulated_regulators$grn)){
    net[(grn == gg) & (source %in% disregulated_regulators[grn==gg]$disregulated_gene),]$grn_effect<-tf_effect
  }
  return(net)
}

#' Use scMultiSim to generate the data with the modified network
#'
#' @param net
#' @param n_cells Number of cells per GRN
#' @param seed Random seed
#' @param tree Differentiation tree required for scMultisim (please refer to documentation)
#' @param noise Whether to add typical noise to the count data
#'
#' @return List of two elements: count matrix, and associated metadata as a data.table
#' @export
#'
#' @examples
generate_data_from_grn<-function(net, n_cells = 750, seed=11, tree = scMultiSim::Phyla1(), noise = FALSE){

  if(!('source' %in% colnames(net))) BBmisc::stopf("Column 'source' does not exist in dataframe")
  if(!('target' %in% colnames(net))) BBmisc::stopf("Column 'target' does not exist in dataframe")
  if(!('grn_effect' %in% colnames(net))) BBmisc::stopf("Column 'grn_effect' does not exist in dataframe, please initialize run 'add_base_effect_to_grn' firs, or initialize weights otherwise")
  if(!('grn' %in% colnames(net))) BBmisc::stopf('GRN colum is not contained in dataframe')
  if (nrow(net)==0) BBmisc::stopf('GRN data frame is empty')

  count_list<-list()
  meta_list<-list()

  for (gg in unique(net$grn)){
    # Subset grn and reorder
    sub.grn<-net[grn == gg]
    sub.grn<-sub.grn[, c('target','source', 'grn_effect')]

    set.seed(seed)
    results <- scMultiSim::sim_true_counts(list(
      # required options
      rand.seed = seed,
      GRN = sub.grn,
      tree = tree,
      num.cells = n_cells,
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
  meta$cell_id <- paste0(meta$cell_id, '-', meta$grn)
  counts<-do.call("cbind", count_list)

  results<-list(counts, meta)
  names(results)<-c('counts', 'meta')
  return(results)
}




#' Generate data for a list of pregenerated networks.
#'
#' @param net
#' @param node_labels
#' @param network_list
#' @param nr_modules
#' @param nr_genes_per_module
#' @param disregulation_type
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
create_gex_data<-function(net,
                          node_labels,
                          net_name,
                          nr_grns = 2,
                          nr_modules = 1,
                          nr_genes_per_module = 1,
                          base_effect = 'standard-normal',
                          mean = 2.5,
                          sd = 1.0,
                          disregulation_type= 'remove-tf-effect-others',
                          seed=11,
                          tf_effect = 5.5,
                          weight_delta = 0.5){


  if (!disregulation_type %in% c('remove-tf-effect-others',
                                 'set-tf-value',
                                 'adjust-weight-of-tf',
                                 'adjust-weight-of-others'
                                 )) {
    BBmisc::stopf('Invalid differential effect name')
  }

  # create graph
  g<-igraph::graph_from_data_frame(net, directed = TRUE, vertices = NULL)

  # select modules and genes to be perturbed
  l<-randomly_select_modules(node_labels, nr_grns, nr_modules)
  disregulated_regulators<-randomly_select_disregulated_node(g, l,node_labels, nr_genes_per_module)

  print(disregulated_regulators)
  if (base_effect == 'standard-normal'){
    # Add gaussian weights as a base
    net<-add_base_effect_gaussian(net, disregulated_regulators, mean=mean, sd = sd)
  }
  else{
  # Add gaussian weights as a base
  net<-add_base_effect_gaussian(net, disregulated_regulators, mean=mean, sd = sd)
  }
  # Set effects of selected regulators t 0
  if (disregulation_type ==  'remove-tf-effect-others'){
    net<-remove_master_tf_effect_from_other_modules(net, disregulated_regulators)
    disregulated_regulators$effect<-disregulation_type
  }else if (disregulation_type ==  'set-tf-value'){
    net<-set_master_tf_effect(net, disregulated_regulators, tf_effect = tf_effect )
    disregulated_regulators$effect<-paste0(disregulation_type, '_', tf_effect)
  } else if (disregulation_type == 'adjust-weight-of-tf'){
    net<- modify_weight_of_edges(net, disregulated_regulators , weight_delta = weight_delta)
    disregulated_regulators$effect<-paste0(disregulation_type, '_', weight_delta)
  }else if(disregulation_type == 'adjust-weight-of-others'){
    net<- modify_weight_of_other_edges(net, disregulated_regulators , weight_delta = weight_delta)
    disregulated_regulators$effect<-paste0(disregulation_type, '_', weight_delta)
  }
  else
  {
    net<-remove_master_tf_effect_from_other_modules(net, disregulated_regulators)
  }

  # Generate data
  generated_data <- generate_data_from_grn(net, seed = seed)

  counts <- generated_data[[1]]
  meta<- generated_data[[2]]
  current_experiment_id<- glue::glue("{net_name}_nm_{nr_modules}_gpm_{nr_genes_per_module}_{disregulation_type}_{seed}")

  info<-list(current_experiment_id, net_name, nr_modules, nr_genes_per_module, disregulation_type, seed)

  results<- list(net, counts, meta, info, disregulated_regulators)
  names(results)<-c('net', 'counts', 'meta', 'info', 'disregulated_regulators')
  return(results)
}


#' Assign modules to be disregulated to each GRN
#'
#' @param node_labels A table with the name of the nodes and the module they belong to
#' @param n_grns Number of GRNs/cell populations to produce
#'
#' @return A data.table containing modules and the grn they have been assigned to
#' @export
#'
#' @examples
randomly_select_modules<-function(node_labels, n_grns, nr_modules){
  # Total number of modules
  nr.module<-length(unique(node_labels$module_greedy))
  # unique module names
  modules<-unique(node_labels$module_greedy)
  max_modules_per_grn <- max(1, floor(nr.module/(n_grns+1)))
  step<-min(max_modules_per_grn, nr_modules)
  print(step)
  groups<-c()
  label<-c()
  for (i in 0:(n_grns-1)){
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
#' @param node_labels A data table of nodenames and their modules.
#' @param nr_disregulated_genes Number of disregulated genes per module (will be rounded down if too many)
#'
#' @return A data frame containing the module, which GRN is assigned to to it and the name of the disregulated gene
#' @export
#'
#' @examples
randomly_select_disregulated_node<-function(g, l, node_labels, nr_disregulated_genes){
  collect<-list()
  regs<-which(igraph::degree(g, mode = 'out')>0)
  disregulated_regulators<-c()
  for (i in 1:nrow(l)){
    dis<-l[i, module]
    duplicate_nodes<-unique(res$egdelist$source[duplicated(res$egdelist$source)])
    eligible<-node_labels[node %in% names(regs)& module_greedy==dis & node %in% duplicate_nodes]$node
    nr_disregulated_genes_c<-min(nr_disregulated_genes, length(eligible))
    nr_disregulated_genes_c<-max(1, nr_disregulated_genes_c)
    dn<-sample(eligible, size = nr_disregulated_genes_c)
    for (d in dn){
      collect[[paste0(d,i)]]<-c(l[i,], d)
    }
  }
  collect<-data.table::rbindlist(collect)
  colnames(collect)<-c('module', 'grn', 'disregulated_gene')
  return(collect)

}

