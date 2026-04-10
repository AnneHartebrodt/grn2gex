library(testthat)
library(igraph)
library(data.table)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_simple_network_dt <- function() {
  data.table::data.table(
    source = c("A", "B", "C", "D", "E"),
    target = c("B", "C", "D", "E", "A"),
    weight = c(1, 1, 1, 1, 1)
  )
}

# A star network: one hub (A) connected to many leaves – a "hedgehog"
make_star_network_dt <- function(n_leaves = 10) {
  data.table::data.table(
    source = "A",
    target = paste0("leaf", seq_len(n_leaves)),
    weight = 1
  )
}

# A larger network with genuine hub structure (enough betweenness diversity)
make_rich_network_dt <- function() {
  set.seed(42)
  n <- 20
  edges <- data.table::data.table(
    source = paste0("g", sample(seq_len(n), 40, replace = TRUE)),
    target = paste0("g", sample(seq_len(n), 40, replace = TRUE)),
    weight = 1
  )
  # Remove self-loops
  edges <- edges[source != target]
  unique(edges)
}

# ---------------------------------------------------------------------------
# createGraph
# ---------------------------------------------------------------------------

test_that("createGraph returns an igraph object", {
  net <- make_simple_network_dt()
  g <- createGraph(net)
  expect_true(igraph::is.igraph(g))
})

test_that("createGraph creates undirected graph by default", {
  net <- make_simple_network_dt()
  g <- createGraph(net)
  expect_false(igraph::is.directed(g))
})

test_that("createGraph creates directed graph when directed = TRUE", {
  net <- make_simple_network_dt()
  g <- createGraph(net, directed = TRUE)
  expect_true(igraph::is.directed(g))
})

test_that("createGraph removes multi-edges", {
  # Duplicate every edge so we have multi-edges
  net <- make_simple_network_dt()
  net_dup <- data.table::rbindlist(list(net, net))
  g <- createGraph(net_dup)
  expect_false(any(igraph::which_multiple(g)))
})

test_that("createGraph stops when first column is not 'source'", {
  net <- make_simple_network_dt()
  data.table::setnames(net, "source", "src")
  expect_error(createGraph(net))
})

test_that("createGraph stops when second column is not 'target'", {
  net <- make_simple_network_dt()
  data.table::setnames(net, "target", "tgt")
  expect_error(createGraph(net))
})

test_that("createGraph preserves all unique vertices", {
  net <- make_simple_network_dt()
  g <- createGraph(net)
  expected_nodes <- sort(unique(c(net$source, net$target)))
  expect_equal(sort(igraph::V(g)$name), expected_nodes)
})

# ---------------------------------------------------------------------------
# remove_irrelevant_networks
# ---------------------------------------------------------------------------

make_graph_list <- function() {
  star_dt <- make_star_network_dt(10)
  star_g   <- igraph::graph_from_data_frame(star_dt, directed = FALSE)

  rich_dt <- make_rich_network_dt()
  rich_g  <- igraph::graph_from_data_frame(rich_dt, directed = FALSE)
  rich_g  <- igraph::simplify(rich_g, remove.loops = TRUE)

  list(star_g, rich_g)
}

test_that("remove_irrelevant_networks removes hedgehog networks", {
  gl <- make_graph_list()
  relevant <- remove_irrelevant_networks(gl)
  # The star graph (1 hub) should be removed; only the rich graph remains
  expect_true(length(relevant) < length(gl))
})

test_that("remove_irrelevant_networks returns a list", {
  gl <- make_graph_list()
  relevant <- remove_irrelevant_networks(gl)
  expect_true(is.list(relevant))
})

test_that("remove_irrelevant_networks keeps networks with >2 betweenness nodes", {
  rich_dt <- make_rich_network_dt()
  rich_g  <- igraph::graph_from_data_frame(rich_dt, directed = FALSE)
  rich_g  <- igraph::simplify(rich_g, remove.loops = TRUE)
  relevant <- remove_irrelevant_networks(list(rich_g))
  expect_equal(length(relevant), 1)
})

test_that("remove_irrelevant_networks returns empty list when all are hedgehogs", {
  star_g <- igraph::graph_from_data_frame(make_star_network_dt(5), directed = FALSE)
  result <- remove_irrelevant_networks(list(star_g))
  expect_equal(length(result), 0)
})

test_that("remove_irrelevant_networks handles an empty input list", {
  result <- remove_irrelevant_networks(list())
  expect_equal(length(result), 0)
})

# ---------------------------------------------------------------------------
# subsample_and_sparsify_network
# ---------------------------------------------------------------------------

make_large_network_dt <- function() {
  set.seed(1)
  sources  <- paste0("TF", seq_len(20))
  targets  <- paste0("gene", seq_len(50))
  data.table::data.table(
    source = sample(sources, 100, replace = TRUE),
    target = sample(targets, 100, replace = TRUE),
    weight = 1
  )
}

test_that("subsample_and_sparsify_network returns a data.table", {
  net <- make_large_network_dt()
  result <- subsample_and_sparsify_network(net, num_sources = 5, max_out_degree = 3)
  expect_true(data.table::is.data.table(result))
})

test_that("subsample_and_sparsify_network respects num_sources", {
  net <- make_large_network_dt()
  result <- subsample_and_sparsify_network(net, num_sources = 5, max_out_degree = 10)
  expect_true(length(unique(result$source)) <= 5)
})

test_that("subsample_and_sparsify_network respects max_out_degree", {
  net <- make_large_network_dt()
  result <- subsample_and_sparsify_network(net, num_sources = 10, max_out_degree = 2)
  out_degrees <- table(result$source)
  expect_true(all(out_degrees <= 2))
})

test_that("subsample_and_sparsify_network stops when not enough unique sources", {
  net <- data.table::data.table(
    source = c("A", "A", "B"),
    target = c("X", "Y", "Z"),
    weight = 1
  )
  expect_error(subsample_and_sparsify_network(net, num_sources = 5, max_out_degree = 2))
})

test_that("subsample_and_sparsify_network stops when first column is not 'source'", {
  net <- make_large_network_dt()
  data.table::setnames(net, "source", "src")
  expect_error(subsample_and_sparsify_network(net, num_sources = 3, max_out_degree = 3))
})

test_that("subsample_and_sparsify_network stops when second column is not 'target'", {
  net <- make_large_network_dt()
  data.table::setnames(net, "target", "tgt")
  expect_error(subsample_and_sparsify_network(net, num_sources = 3, max_out_degree = 3))
})

test_that("subsample_and_sparsify_network removes target-as-source rows", {
  # After sparsification targets that also appear as sources are removed
  net <- make_large_network_dt()
  result <- subsample_and_sparsify_network(net, num_sources = 5, max_out_degree = 5)
  # Targets in the result should not appear as sources in the result
  expect_true(length(intersect(result$target, result$source)) == 0)
})

# ---------------------------------------------------------------------------
# cluster (recursive clustering)
# ---------------------------------------------------------------------------

make_medium_graph <- function(n = 80, seed = 7) {
  set.seed(seed)
  el <- data.table::data.table(
    source = paste0("n", sample(seq_len(n), n * 3, replace = TRUE)),
    target = paste0("n", sample(seq_len(n), n * 3, replace = TRUE))
  )
  el <- el[source != target]
  el <- unique(el)
  igraph::graph_from_data_frame(el, directed = FALSE)
}

test_that("cluster returns a list", {
  g <- make_medium_graph(80)
  result <- cluster(g, min_nodes = 10, max_nodes = 50)
  expect_true(is.list(result))
})

test_that("cluster returns graphs with at most max_nodes nodes", {
  g <- make_medium_graph(120)
  result <- cluster(g, min_nodes = 10, max_nodes = 50)
  sizes <- vapply(result, igraph::vcount, integer(1))
  expect_true(all(sizes <= 50))
})

test_that("cluster returns igraph objects", {
  g <- make_medium_graph(80)
  result <- cluster(g, min_nodes = 5, max_nodes = 40)
  for (sub in result) {
    expect_true(igraph::is.igraph(sub))
  }
})

# ---------------------------------------------------------------------------
# clusterNetwork
# ---------------------------------------------------------------------------

test_that("clusterNetwork returns a list of igraph objects", {
  net <- make_rich_network_dt()
  result <- clusterNetwork(net, min_nodes = 2, max_nodes = 8)
  expect_true(is.list(result))
  for (g in result) {
    expect_true(igraph::is.igraph(g))
  }
})

test_that("clusterNetwork stops when first column is not 'source'", {
  net <- make_rich_network_dt()
  data.table::setnames(net, "source", "from")
  expect_error(clusterNetwork(net))
})

# ---------------------------------------------------------------------------
# restore_original_directionality
# ---------------------------------------------------------------------------

test_that("restore_original_directionality returns a named list with 4 elements", {
  net <- make_rich_network_dt()
  orig_g  <- igraph::graph_from_data_frame(net, directed = TRUE)
  orig_g  <- igraph::simplify(orig_g, remove.loops = TRUE)
  undirected_g <- igraph::as.undirected(orig_g)
  result <- restore_original_directionality(undirected_g, orig_g)

  expect_true(is.list(result))
  expect_equal(length(result), 4)
  # Note: "egdelist" (not "edgelist") is the key name as defined in restore_original_directionality
  expect_true(all(c("subgraph", "clustering", "nodelist", "egdelist") %in% names(result)))
})

test_that("restore_original_directionality nodelist has expected columns", {
  net <- make_rich_network_dt()
  orig_g  <- igraph::graph_from_data_frame(net, directed = TRUE)
  orig_g  <- igraph::simplify(orig_g, remove.loops = TRUE)
  undirected_g <- igraph::as.undirected(orig_g)
  result <- restore_original_directionality(undirected_g, orig_g)

  expect_true("node" %in% colnames(result$nodelist))
  expect_true("module_greedy" %in% colnames(result$nodelist))
})

test_that("restore_original_directionality edgelist has source/target columns", {
  net <- make_rich_network_dt()
  orig_g  <- igraph::graph_from_data_frame(net, directed = TRUE)
  orig_g  <- igraph::simplify(orig_g, remove.loops = TRUE)
  undirected_g <- igraph::as.undirected(orig_g)
  result <- restore_original_directionality(undirected_g, orig_g)

  # Note: "egdelist" (not "edgelist") is the key name as defined in restore_original_directionality
  expect_true(all(c("source", "target") %in% colnames(result$egdelist)))
})
