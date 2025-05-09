% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_networks.R
\name{modify_weight_of_edges}
\alias{modify_weight_of_edges}
\title{Modifify the weights to induce the differential behaviour for the networks.
This function modifies the weight of the edges in the mapping.
this means one has to be careful when working with more than two grns, because
the effect is only added/substracted to the selected modules
Good for amping up the effect of the TF}
\usage{
modify_weight_of_edges(net, disregulated_regulators, weight_delta = 0.5)
}
\arguments{
\item{net}{The network with source, target and base_effect columns}

\item{disregulated_regulators}{The mapping for the disregulated modules}

\item{weight_delta}{The effect to be added to the selected edges}
}
\value{
net A network with an additional column grn_effect which is a copy of the base effect but modified according
to the network mapping
}
\description{
Modifify the weights to induce the differential behaviour for the networks.
This function modifies the weight of the edges in the mapping.
this means one has to be careful when working with more than two grns, because
the effect is only added/substracted to the selected modules
Good for amping up the effect of the TF
}
