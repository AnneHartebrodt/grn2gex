% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_modules.R
\name{clusterNetwork}
\alias{clusterNetwork}
\title{Wrapper function to cluster a network in a datatable.}
\usage{
clusterNetwork(network_dt, min_nodes = 50, max_nodes = 200)
}
\arguments{
\item{network_dt}{A data.table containing a colum 'source' and 'target' as the first and second column}
}
\value{
A list of igraph Graph objects
}
\description{
Wrapper function to cluster a network in a datatable.
}
\examples{
collectri.file <- '/tmp/collectri.tsv'
collectri <- loadOrDownloadCollectTRI(collectri.file)
colnames(collectri)<- c('source', 'target', 'sad')
netlist<-clusterNetwork(collectri )
}
