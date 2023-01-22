#' @title read_data
#' @author Thalla
#' @param dataset_names names of datasets
#' @param datapath path to data in mounted folder, assumption1: data folder is called 'data', assumption2: 'original_name' contains the full name of the file
#' @param analysis_params params
#' @param func reader function
#' @return datasets (list)
#' @export
read_data <- function(dataset_names,datapath, analysis_params, func){
  datasets <- lapply(dataset_names, function(dataset){
    func(paste0(datapath, "/", dataset, "/data/", analysis_params$other[[dataset]]$original_name))
  })
  names(datasets) <- dataset_names
  datasets
}


#' Read Dataset Analysis Params
#' @param datasetpath path of dataset
#' @param subset_name name of subset
#' @description read in dataset specific parameters
#' @export
read_ds_analysis_params <- function(dataset_name, datasetpath){
  yaml::read_yaml(paste(datasetpath, dataset_name, "analysis_params.yaml", sep = "/"))
}


#' read_liana_methods
#' @description same as 'res'/result from liana_wrap
#' @param ds_name dataset name
#' @export
read_liana_methods <- function(dataset, clustercolname){
  config <- dataset$config
  res <- sapply(c(config$methods), function(methodname){
    read.csv(paste0(config$resultpath, "/", clustercolname, "/", methodname, ".csv"))
  })
  names(res) <- config$methods
  sapply(res, tibble::as_tibble)
}
