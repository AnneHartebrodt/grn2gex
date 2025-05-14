suppressMessages(require(grn2gex))
suppressMessages(require(scMultiSim))
suppressMessages(require(data.table))
suppressMessages(require(optparse))
suppressMessages(library(yaml))
suppressMessages(require(SingleCellExperiment))
suppressMessages(library(sceasy))
suppressMessages(library(reticulate))

## Save as h5ad file
use_python('/opt/conda/bin/python')

use_condaenv('base')
print(py_config())
# 
# Define command line options
option_list <- list(
  make_option(
    c("-c", "--config"),
    type = "character",
    default = NULL,
    help = "Path to the YAML configuration file",
    metavar = "FILE"
  )
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the config file is provided
if (is.null(opt$config)) {
  stop("Error: Please provide a YAML config file using the --config option.")
}

# Load the YAML configuration file
print(opt$config)
config <- yaml.load_file(opt$config)


# Check if the 'data_simulation' section exists
if (!"data_simulation" %in% names(config)) {
  stop("Error: Missing 'data_simulation' section in config file.")
}

###### DATA SIMULATION #######################################################################
# Define required variables within 'data_simulation'
required_vars_simulation <- c("n_cells", "n_celltypes", "edgelist", 'nodelist', "seed", "modification_type", "n_genes_per_module", 'n_celltypes', 'n_modules')

# Check for missing variables in the 'data_simulation' section
missing_vars <- required_vars_simulation[!required_vars_simulation %in% names(config$data_simulation)]
if (length(missing_vars) > 0) {
  stop(paste("Error: Missing required variables in 'data_simulation' section:", paste(missing_vars, collapse = ", ")))
}

# # Check if 'input_grn' exists and is a valid file path
# if (!file.exists(config$data_simulation$edgelist)) {
#   stop(paste("Error: Specified file", config$data_simulation$edgelist, "does not exist."))
# }
# # Check if 'input_grn' exists and is a valid file path
# if (!file.exists(config$data_simulation$nodelist)) {
#   stop(paste("Error: Specified file", config$data_simulation$nodelist, "does not exist."))
# }

read_dataframes<-function(path, filelist){
    list_edgelist<-list()
    counter<-1
    for(e in filelist){
        edgelist<-fread(file.path(path, e), sep='\t')
        edgelist$module<-counter
        list_edgelist[[counter]]<-edgelist
        counter<-counter+1
    }
    edgelist<-rbindlist(list_edgelist)

    return(edgelist)
}


output_dir<-file.path('/usr/src/app/output')



edgelist<-read_dataframes('/usr/src/app/input',config$data_simulation$edgelist )
nodelist<-read_dataframes('/usr/src/app/input', config$data_simulation$nodelist)

commonlist<-read_dataframes('/usr/src/app/input',config$data_simulation$common_edges )

# config <- yaml.load_file("/home/bionets-og86asub/Documents/netmap/thenetmap/dockerize-grn2gex/config.yaml")
# edgelist<-fread(file.path("/home/bionets-og86asub/Documents/netmap/thenetmap/NetMap_LRP/data/simulation/collectri_subnetworks_2/net_135_44105/", config$data_simulation$edgelist), sep='\t')
# nodelist<-fread(file.path("/home/bionets-og86asub/Documents/netmap/thenetmap/NetMap_LRP/data/simulation/collectri_subnetworks_2/net_135_44105/", config$data_simulation$nodelist), sep = '\t')
# output_dir<-file.path("/home/bionets-og86asub/Documents/netmap/thenetmap/NetMap_LRP/data/simulation", config$data_simulation$dataset_id)

print('UPDATED')

dataset<-create_gex_data_easy(net = edgelist, 
                         common_net = commonlist,
                         net_name = config$data_simulation$dataset_id,
                         seed = config$data_simulation$seed,
                         base_effect = config$data_simulation$base_effect,
                         mean = config$data_simulation$mean,
                         sd = config$data_simulation$sd)

gex.dir<-file.path(output_dir,config$data_simulation$dataset_id )
gex.dir<-save_generated_data(dataset$net, dataset$counts, dataset$meta, gex_dir =gex.dir, disregulated_info = NULL)

print(gex.dir)


sce <- SingleCellExperiment(list(counts=dataset$counts), colData=dataset$meta, rowData = DataFrame(genes = rownames(dataset$counts)))

out_path <- file.path(gex.dir, 'data.h5ad')

sce <- sceasy::convertFormat(sce, from="sce", to="anndata",
                      outFile=out_path)
