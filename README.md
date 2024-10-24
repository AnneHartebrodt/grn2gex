
# grn2gex

<!-- badges: start -->
<!-- badges: end -->

grn2gex can be used to simulate gene expression data from differential gene regulatory networks

## Installation

You can install the development version of grn2gex like so:

``` r
devtools::install_github('https://github.com/AnneHartebrodt/grn2gex')
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
require(grn2gex)
require(scMultiSim)
require(data.table)
net.dir<-'/path/to/store/your/project/data/simulation/collectri_subnetworks'
collectri <- '/path/to/store/your/project/data/public/collecTRI/collectri_2024-10-24.tsv'
gex.dir<-'/path/to/store/your/project/data/simulation/collectri_simulated_data'

dir.create(gex.dir, recursive = TRUE)
```

## Cluster network
``` r
nr.modules <- 2
nr.genes.per.module<-2
disregulation.type <- 'remove-tf-effect-others'
seed <- 11

# Download/load network and cluster 
clusterCollecTRI(net.dir, collectri)
# remove networks with only 1 hubnode 
edgelist<-remove_irrelevant_networks(net.dir)
```

## Read list of relevant networks and generate data
``` r
network.stat <-fread(file.path(net.dir, 'network_stats.tsv'))

network.stat<-network.stat[type == 'modularized']
network.list <- sapply(network.stat$V1, function(x) dirname(x))

disregulation.type = 'remove-tf-effect-others'
nl<-network.list[1]
batch_create_gex_data(net.dir,
                      gex.dir,  
                      nl,
                      disregulation.type = disregulation.type)
```
