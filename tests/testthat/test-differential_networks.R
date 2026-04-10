library(testthat)
library(data.table)
library(igraph)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# A minimal edge-list with two GRNs already assigned
make_net_with_grn <- function() {
  data.table::data.table(
    source      = c("TF1", "TF1", "TF2", "TF2", "TF3", "TF3"),
    target      = c("G1",  "G2",  "G3",  "G4",  "G5",  "G6"),
    base_effect = c(2.0,   2.5,   3.0,   2.0,   1.5,   2.0),
    grn         = c(0L,    0L,    1L,    1L,    0L,    1L)
  )
}

# Disregulation table matching the network above
make_disreg <- function() {
  data.table::data.table(
    module            = c(1L, 2L),
    grn               = c(0L, 1L),
    disregulated_gene = c("TF1", "TF2")
  )
}

# ---------------------------------------------------------------------------
# add_base_effect_gaussian
# ---------------------------------------------------------------------------

test_that("add_base_effect_gaussian adds base_effect column", {
  net <- data.table::data.table(
    source = c("TF1", "TF2"),
    target = c("G1",  "G2")
  )
  disreg <- data.table::data.table(
    grn               = c(0L, 1L),
    disregulated_gene = c("TF1", "TF2")
  )
  set.seed(1)
  result <- add_base_effect_gaussian(net, disreg, mean = 2.5, sd = 1.0)
  expect_true("base_effect" %in% colnames(result))
})

test_that("add_base_effect_gaussian replicates rows for each GRN", {
  net <- data.table::data.table(
    source = c("TF1", "TF2"),
    target = c("G1",  "G2")
  )
  disreg <- data.table::data.table(
    grn               = c(0L, 1L),
    disregulated_gene = c("TF1", "TF2")
  )
  set.seed(1)
  result <- add_base_effect_gaussian(net, disreg)
  # 2 original rows × 2 GRNs = 4 rows
  expect_equal(nrow(result), nrow(net) * length(unique(disreg$grn)))
})

test_that("add_base_effect_gaussian adds grn column", {
  net <- data.table::data.table(
    source = c("TF1"),
    target = c("G1")
  )
  disreg <- data.table::data.table(
    grn               = c(0L, 1L, 2L),
    disregulated_gene = c("TF1", "TF1", "TF1")
  )
  set.seed(1)
  result <- add_base_effect_gaussian(net, disreg)
  expect_true("grn" %in% colnames(result))
  expect_equal(sort(unique(result$grn)), c(0L, 1L, 2L))
})

test_that("add_base_effect_gaussian respects mean and sd parameters", {
  net <- data.table::data.table(
    source = paste0("TF", seq_len(500)),
    target = paste0("G",  seq_len(500))
  )
  disreg <- data.table::data.table(grn = 0L, disregulated_gene = "TF1")
  set.seed(42)
  result <- add_base_effect_gaussian(net, disreg, mean = 10.0, sd = 0.01)
  # With 500 samples and tiny sd, mean should be very close to 10
  expect_true(abs(mean(result$base_effect) - 10.0) < 0.1)
})

# ---------------------------------------------------------------------------
# modify_weight_of_edges
# ---------------------------------------------------------------------------

test_that("modify_weight_of_edges adds grn_effect column", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- modify_weight_of_edges(net, disreg, weight_delta = 0.5)
  expect_true("grn_effect" %in% colnames(result))
})

test_that("modify_weight_of_edges increases weight for selected source in correct GRN", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- modify_weight_of_edges(net, disreg, weight_delta = 0.5)

  # TF1 in grn=0 should have grn_effect = base_effect + 0.5
  affected <- result[source == "TF1" & grn == 0L]
  expect_true(all(affected$grn_effect == affected$base_effect + 0.5))
})

test_that("modify_weight_of_edges does not change unaffected rows", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- modify_weight_of_edges(net, disreg, weight_delta = 0.5)

  # TF3 is not a disregulated gene in either GRN
  unaffected <- result[source == "TF3"]
  expect_true(all(unaffected$grn_effect == unaffected$base_effect))
})

# ---------------------------------------------------------------------------
# modify_weight_of_other_edges
# ---------------------------------------------------------------------------

test_that("modify_weight_of_other_edges adds grn_effect column", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- modify_weight_of_other_edges(net, disreg, weight_delta = 0.5)
  expect_true("grn_effect" %in% colnames(result))
})

test_that("modify_weight_of_other_edges modifies TF in opposite GRNs", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- modify_weight_of_other_edges(net, disreg, weight_delta = 0.5)

  # TF1 is assigned to grn=0. In grn=1, TF1's weight should be changed.
  affected <- result[source == "TF1" & grn == 1L]
  if (nrow(affected) > 0) {
    expect_true(all(affected$grn_effect == affected$base_effect + 0.5))
  }
})

# ---------------------------------------------------------------------------
# remove_master_tf_effect_from_other_modules
# ---------------------------------------------------------------------------

test_that("remove_master_tf_effect_from_other_modules adds grn_effect column", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- remove_master_tf_effect_from_other_modules(net, disreg)
  expect_true("grn_effect" %in% colnames(result))
})

test_that("remove_master_tf_effect_from_other_modules sets cross-grn effect to 0.01", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- remove_master_tf_effect_from_other_modules(net, disreg)

  # TF1 belongs to grn=0; in grn=1 its effect should be 0.01
  cross_grn <- result[source == "TF1" & grn == 1L]
  if (nrow(cross_grn) > 0) {
    expect_true(all(cross_grn$grn_effect == 0.01))
  }
})

test_that("remove_master_tf_effect_from_other_modules preserves in-module effect", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- remove_master_tf_effect_from_other_modules(net, disreg)

  # TF1 belongs to grn=0; within grn=0 its effect should equal base_effect
  in_module <- result[source == "TF1" & grn == 0L]
  expect_true(all(in_module$grn_effect == in_module$base_effect))
})

# ---------------------------------------------------------------------------
# set_master_tf_effect
# ---------------------------------------------------------------------------

test_that("set_master_tf_effect stops when base_effect column is missing", {
  net <- data.table::data.table(
    source = "TF1",
    target = "G1",
    grn    = 0L
  )
  disreg <- make_disreg()
  expect_error(set_master_tf_effect(net, disreg))
})

test_that("set_master_tf_effect adds grn_effect column", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- set_master_tf_effect(net, disreg, tf_effect = 5.0)
  expect_true("grn_effect" %in% colnames(result))
})

test_that("set_master_tf_effect sets exact tf_effect value for selected TF in correct GRN", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- set_master_tf_effect(net, disreg, tf_effect = 9.9)

  affected <- result[source == "TF1" & grn == 0L]
  expect_true(all(affected$grn_effect == 9.9))
})

test_that("set_master_tf_effect does not alter unaffected rows", {
  net    <- make_net_with_grn()
  disreg <- make_disreg()
  result <- set_master_tf_effect(net, disreg, tf_effect = 9.9)

  unaffected <- result[source == "TF3"]
  expect_true(all(unaffected$grn_effect == unaffected$base_effect))
})

# ---------------------------------------------------------------------------
# randomly_select_modules
# ---------------------------------------------------------------------------

make_node_labels <- function(n_modules = 6, nodes_per_module = 5) {
  data.table::data.table(
    node          = paste0("gene", seq_len(n_modules * nodes_per_module)),
    module_greedy = rep(seq_len(n_modules), each = nodes_per_module)
  )
}

test_that("randomly_select_modules returns a data.table with module and grn columns", {
  set.seed(1)
  nl     <- make_node_labels()
  result <- randomly_select_modules(nl, n_grns = 2, nr_modules = 1)
  expect_true(data.table::is.data.table(result))
  expect_true(all(c("module", "grn") %in% colnames(result)))
})

test_that("randomly_select_modules returns correct number of rows (n_grns * nr_modules)", {
  set.seed(1)
  nl     <- make_node_labels(n_modules = 6)
  result <- randomly_select_modules(nl, n_grns = 3, nr_modules = 1)
  expect_equal(nrow(result), 3)
})

test_that("randomly_select_modules assigns distinct grn labels", {
  set.seed(1)
  nl     <- make_node_labels(n_modules = 6)
  result <- randomly_select_modules(nl, n_grns = 3, nr_modules = 1)
  expect_equal(length(unique(result$grn)), 3)
})

# ---------------------------------------------------------------------------
# randomly_select_disregulated_node
# ---------------------------------------------------------------------------

make_directed_graph_for_selection <- function() {
  # Create a directed graph where TFs (sources) have out-degree > 1
  el <- data.table::data.table(
    source = c("TF1", "TF1", "TF2", "TF2", "TF3", "TF3"),
    target = c("G1",  "G2",  "G3",  "G4",  "G5",  "G6")
  )
  igraph::graph_from_data_frame(el, directed = TRUE)
}

test_that("randomly_select_disregulated_node returns a data.table", {
  set.seed(1)
  nl     <- make_node_labels(n_modules = 4)
  g      <- make_directed_graph_for_selection()
  # Add TFs to the node labels
  nl2 <- data.table::data.table(
    node          = c("TF1", "TF2", "TF3", "G1", "G2", "G3", "G4", "G5", "G6"),
    module_greedy = c(1L,    2L,    3L,    1L,   1L,   2L,   2L,   3L,   3L)
  )
  l  <- data.table::data.table(module = c(1L, 2L), grn = c(0L, 1L))
  result <- randomly_select_disregulated_node(g, l, nl2, nr_disregulated_genes = 1)
  expect_true(data.table::is.data.table(result))
})

test_that("randomly_select_disregulated_node has required columns", {
  set.seed(1)
  nl2 <- data.table::data.table(
    node          = c("TF1", "TF2", "TF3", "G1", "G2", "G3", "G4", "G5", "G6"),
    module_greedy = c(1L,    2L,    3L,    1L,   1L,   2L,   2L,   3L,   3L)
  )
  g <- make_directed_graph_for_selection()
  l <- data.table::data.table(module = c(1L, 2L), grn = c(0L, 1L))
  result <- randomly_select_disregulated_node(g, l, nl2, nr_disregulated_genes = 1)
  expect_true(all(c("module", "grn", "disregulated_gene") %in% colnames(result)))
})

test_that("randomly_select_disregulated_node selected gene has out-degree > 1", {
  set.seed(42)
  nl2 <- data.table::data.table(
    node          = c("TF1", "TF2", "G1", "G2", "G3", "G4"),
    module_greedy = c(1L,    2L,    1L,   1L,   2L,   2L)
  )
  g <- make_directed_graph_for_selection()
  l <- data.table::data.table(module = 1L, grn = 0L)
  result <- randomly_select_disregulated_node(g, l, nl2, nr_disregulated_genes = 1)
  sel <- result$disregulated_gene
  out_deg <- igraph::degree(g, v = sel, mode = "out")
  expect_true(all(out_deg > 1))
})

# ---------------------------------------------------------------------------
# generate_data_from_grn – input validation only (no scMultiSim call)
# ---------------------------------------------------------------------------

test_that("generate_data_from_grn stops when 'source' column missing", {
  net <- data.table::data.table(
    target     = "G1",
    grn_effect = 1.0,
    grn        = 0L
  )
  expect_error(generate_data_from_grn(net))
})

test_that("generate_data_from_grn stops when 'target' column missing", {
  net <- data.table::data.table(
    source     = "TF1",
    grn_effect = 1.0,
    grn        = 0L
  )
  expect_error(generate_data_from_grn(net))
})

test_that("generate_data_from_grn stops when 'grn_effect' column missing", {
  net <- data.table::data.table(
    source = "TF1",
    target = "G1",
    grn    = 0L
  )
  expect_error(generate_data_from_grn(net))
})

test_that("generate_data_from_grn stops when 'grn' column missing", {
  net <- data.table::data.table(
    source     = "TF1",
    target     = "G1",
    grn_effect = 1.0
  )
  expect_error(generate_data_from_grn(net))
})

test_that("generate_data_from_grn stops when network is empty", {
  net <- data.table::data.table(
    source     = character(0),
    target     = character(0),
    grn_effect = numeric(0),
    grn        = integer(0)
  )
  expect_error(generate_data_from_grn(net))
})

# ---------------------------------------------------------------------------
# create_gex_data – input validation only
# ---------------------------------------------------------------------------

test_that("create_gex_data stops on invalid disregulation_type", {
  net <- data.table::data.table(source = "TF1", target = "G1", weight = 1)
  nl  <- data.table::data.table(node = c("TF1", "G1"), module_greedy = c(1L, 1L))
  expect_error(
    create_gex_data(net, nl, "test_net", disregulation_type = "invalid-type")
  )
})

test_that("create_gex_data stops when nr_grns exceeds unique modules", {
  net <- data.table::data.table(source = "TF1", target = "G1", weight = 1)
  nl  <- data.table::data.table(node = c("TF1", "G1"), module_greedy = c(1L, 1L))
  # Only 1 unique module but requesting 3 GRNs
  expect_error(
    create_gex_data(net, nl, "test_net", nr_grns = 3,
                    disregulation_type = "set-tf-value")
  )
})

# ---------------------------------------------------------------------------
# create_gex_data_easy – input validation only
# ---------------------------------------------------------------------------

test_that("create_gex_data_easy stops when 'module' column is missing from net", {
  net <- data.table::data.table(
    source = "TF1",
    target = "G1"
  )
  common_net <- data.table::data.table(
    source = "TF1",
    target = "G2"
  )
  expect_error(
    create_gex_data_easy(net, common_net, "test_net")
  )
})
