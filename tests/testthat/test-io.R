library(testthat)
library(data.table)
library(igraph)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Minimal count matrix: genes × cells
make_counts <- function(n_genes = 5, n_cells = 10) {
  m <- matrix(rpois(n_genes * n_cells, lambda = 10),
              nrow = n_genes, ncol = n_cells)
  rownames(m) <- paste0("TF", seq_len(n_genes))
  colnames(m) <- paste0("cell", seq_len(n_cells))
  m
}

# Matching meta data.table
make_meta <- function(n_cells = 10, n_grns = 2) {
  data.table::data.table(
    cell_id = paste0("cell", seq_len(n_cells)),
    grn     = rep(seq_len(n_grns) - 1L, length.out = n_cells),
    pop     = rep(1L, n_cells)
  )
}

# ---------------------------------------------------------------------------
# save_generated_data
# ---------------------------------------------------------------------------

test_that("save_generated_data creates the output directory", {
  tmp    <- tempfile()
  net    <- data.table::data.table(
    source     = "TF1", target = "G1",
    grn_effect = 1.0, grn = 0L
  )
  counts <- make_counts()
  meta   <- make_meta()

  # save_plots tries to use scMultiSim / ComplexHeatmap; mock it by replacing
  # save_plots temporarily
  local_mocked_bindings(
    save_plots = function(...) invisible(NULL),
    .package   = "grn2gex",
    {
      save_generated_data(net, counts, meta, tmp)
    }
  )
  expect_true(dir.exists(tmp))
})

test_that("save_generated_data writes gex.tsv", {
  tmp    <- tempfile()
  net    <- data.table::data.table(
    source     = "TF1", target = "G1",
    grn_effect = 1.0, grn = 0L
  )
  counts <- make_counts()
  meta   <- make_meta()

  local_mocked_bindings(
    save_plots = function(...) invisible(NULL),
    .package   = "grn2gex",
    {
      save_generated_data(net, counts, meta, tmp)
    }
  )
  expect_true(file.exists(file.path(tmp, "gex.tsv")))
})

test_that("save_generated_data writes metadata.tsv", {
  tmp    <- tempfile()
  net    <- data.table::data.table(
    source     = "TF1", target = "G1",
    grn_effect = 1.0, grn = 0L
  )
  counts <- make_counts()
  meta   <- make_meta()

  local_mocked_bindings(
    save_plots = function(...) invisible(NULL),
    .package   = "grn2gex",
    {
      save_generated_data(net, counts, meta, tmp)
    }
  )
  expect_true(file.exists(file.path(tmp, "metadata.tsv")))
})

test_that("save_generated_data writes net.tsv", {
  tmp    <- tempfile()
  net    <- data.table::data.table(
    source     = "TF1", target = "G1",
    grn_effect = 1.0, grn = 0L
  )
  counts <- make_counts()
  meta   <- make_meta()

  local_mocked_bindings(
    save_plots = function(...) invisible(NULL),
    .package   = "grn2gex",
    {
      save_generated_data(net, counts, meta, tmp)
    }
  )
  expect_true(file.exists(file.path(tmp, "net.tsv")))
})

test_that("save_generated_data writes disregulation_info.tsv when provided", {
  tmp    <- tempfile()
  net    <- data.table::data.table(
    source     = "TF1", target = "G1",
    grn_effect = 1.0, grn = 0L
  )
  counts    <- make_counts()
  meta      <- make_meta()
  disreg    <- data.table::data.table(
    module = 1L, grn = 0L, disregulated_gene = "TF1"
  )

  local_mocked_bindings(
    save_plots = function(...) invisible(NULL),
    .package   = "grn2gex",
    {
      save_generated_data(net, counts, meta, tmp, disregulated_info = disreg)
    }
  )
  expect_true(file.exists(file.path(tmp, "disregulation_info.tsv")))
})

test_that("save_generated_data does NOT write disregulation_info.tsv when not provided", {
  tmp    <- tempfile()
  net    <- data.table::data.table(
    source     = "TF1", target = "G1",
    grn_effect = 1.0, grn = 0L
  )
  counts <- make_counts()
  meta   <- make_meta()

  local_mocked_bindings(
    save_plots = function(...) invisible(NULL),
    .package   = "grn2gex",
    {
      save_generated_data(net, counts, meta, tmp)
    }
  )
  expect_false(file.exists(file.path(tmp, "disregulation_info.tsv")))
})

test_that("save_generated_data returns the output directory path", {
  tmp    <- tempfile()
  net    <- data.table::data.table(
    source     = "TF1", target = "G1",
    grn_effect = 1.0, grn = 0L
  )
  counts <- make_counts()
  meta   <- make_meta()

  result <- local_mocked_bindings(
    save_plots = function(...) invisible(NULL),
    .package   = "grn2gex",
    {
      save_generated_data(net, counts, meta, tmp)
    }
  )
  expect_equal(result, tmp)
})

test_that("save_generated_data gex.tsv has correct number of columns", {
  tmp    <- tempfile()
  n_genes <- 4
  counts <- make_counts(n_genes = n_genes, n_cells = 6)
  meta   <- make_meta(n_cells = 6)
  net    <- data.table::data.table(
    source     = "TF1", target = "G1",
    grn_effect = 1.0, grn = 0L
  )

  local_mocked_bindings(
    save_plots = function(...) invisible(NULL),
    .package   = "grn2gex",
    {
      save_generated_data(net, counts, meta, tmp)
    }
  )
  gex <- data.table::fread(file.path(tmp, "gex.tsv"), sep = "\t")
  # Rows = cells, columns = genes
  expect_equal(ncol(gex), n_genes)
  expect_equal(nrow(gex), 6)
})

# ---------------------------------------------------------------------------
# save_small_subgraph
# ---------------------------------------------------------------------------

make_test_orig_graph <- function() {
  el <- data.table::data.table(
    source = c("A", "B", "C", "D", "E", "F", "A", "B", "C"),
    target = c("B", "C", "D", "E", "F", "A", "D", "E", "F")
  )
  igraph::graph_from_data_frame(el, directed = TRUE)
}

make_test_subgraph <- function(orig) {
  igraph::as.undirected(orig)
}

test_that("save_small_subgraph creates a subdirectory inside network_dir", {
  tmp      <- tempfile()
  dir.create(tmp)
  orig_g   <- make_test_orig_graph()
  sub_g    <- make_test_subgraph(orig_g)

  save_small_subgraph(sub_g, tmp, orig_g)

  subdirs <- list.dirs(tmp, recursive = FALSE)
  expect_true(length(subdirs) >= 1)
})

test_that("save_small_subgraph creates edges.tsv in the subdirectory", {
  tmp      <- tempfile()
  dir.create(tmp)
  orig_g   <- make_test_orig_graph()
  sub_g    <- make_test_subgraph(orig_g)

  save_small_subgraph(sub_g, tmp, orig_g)

  edge_files <- list.files(tmp, pattern = "edges.tsv", recursive = TRUE)
  expect_true(length(edge_files) >= 1)
})

test_that("save_small_subgraph creates nodes.tsv in the subdirectory", {
  tmp      <- tempfile()
  dir.create(tmp)
  orig_g   <- make_test_orig_graph()
  sub_g    <- make_test_subgraph(orig_g)

  save_small_subgraph(sub_g, tmp, orig_g)

  node_files <- list.files(tmp, pattern = "nodes.tsv", recursive = TRUE)
  expect_true(length(node_files) >= 1)
})

test_that("save_small_subgraph creates net.pdf in the subdirectory", {
  tmp      <- tempfile()
  dir.create(tmp)
  orig_g   <- make_test_orig_graph()
  sub_g    <- make_test_subgraph(orig_g)

  save_small_subgraph(sub_g, tmp, orig_g)

  pdf_files <- list.files(tmp, pattern = "net.pdf", recursive = TRUE)
  expect_true(length(pdf_files) >= 1)
})

test_that("save_small_subgraph edges.tsv has source and target columns", {
  tmp      <- tempfile()
  dir.create(tmp)
  orig_g   <- make_test_orig_graph()
  sub_g    <- make_test_subgraph(orig_g)

  save_small_subgraph(sub_g, tmp, orig_g)

  edge_file <- list.files(tmp, pattern = "edges.tsv", recursive = TRUE, full.names = TRUE)[1]
  edges_dt  <- data.table::fread(edge_file)
  expect_true(all(c("source", "target") %in% colnames(edges_dt)))
})

test_that("save_small_subgraph nodes.tsv has node and module_greedy columns", {
  tmp      <- tempfile()
  dir.create(tmp)
  orig_g   <- make_test_orig_graph()
  sub_g    <- make_test_subgraph(orig_g)

  save_small_subgraph(sub_g, tmp, orig_g)

  node_file <- list.files(tmp, pattern = "nodes.tsv", recursive = TRUE, full.names = TRUE)[1]
  nodes_dt  <- data.table::fread(node_file)
  expect_true(all(c("node", "module_greedy") %in% colnames(nodes_dt)))
})
