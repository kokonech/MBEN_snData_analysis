#' Get Subsets
#' @param dataset_name name of dataset
#' @param subset_category category of subset which is a metadata column
#' @return subsets
#' @export
get_subsets = function (dataset_name, subset_category){
  # extract info from analysis_params
  subset_info <- analysis_params$other[[dataset_name]]$subs[[subset_category]]
  subset_info_collapsed = sapply(subset_info, function(elems){
    paste0(elems, collapse = "")
  })

  # subset the data
  subsets <- lapply(subset_info[-1], function(elems){
    subset(x = datasets[[dataset_name]], subset = !!as.symbol(subset_category) %in% unlist(elems))
  })
  names(subsets) <- paste0(dataset_name, "_", subset_info_collapsed[1], subset_info_collapsed[-1]) # all_t_clust34

  # save subset
  lapply(names(subsets), function(subset_name){
    ds_analysis_params <- read_ds_analysis_params(subset_name, datasetpath)
    path <- paste0(ds_analysis_params$paths$datapath_tmp) %>% paste0(collapse = "/") %>% paste0("/", subset_name, ".rds")
    saveRDS(subsets[subset_name], path)
  })
  subsets
}

#' get_subsets_of_dataset To add subset names via analysis_params$other$<dataset_name>$subs$<subset_category>
#' @param datasetname name of dataset
#' @param metacol for example "seurat_clusters"
#' @export
get_subsets_of_dataset <- function(datasets, datasetname, metacol){
  subsets <-  get_subsets(datasetname, metacol)
  # Add subsets to datasets
  datasets <- c(datasets, subsets)
  datasets
}
