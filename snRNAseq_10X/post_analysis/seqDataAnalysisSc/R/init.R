#' @title save_tool_config
#' @author Thalla
#' @description Creates and saves a config file for a specific tool. Each dataset is added with the relevant paths. File is saved in the version folder.
#' @param datasets list of datasets
#' @param toolname name of tool
#' @keywords config init
#' @return config table
#' @export
save_tool_config <- function(datasets, datasetpath, toolname){
  # create_tool_config reads dataset specific analysis_params file, creates data path and specific for liana result & figure path
  # @param dataset_name name of dataset
  # @param datasetpath path of dataset
  # @return dataframe with name, datapath, resultpath, figpath
  create_tool_config <- function(dataset_name, datasetpath, toolname){
    ap <- yaml::read_yaml(paste(datasetpath, dataset_name, "analysis_params.yaml", sep = "/"))
    datapath <- paste0(ap$paths$datapath_tmp) %>% paste0(collapse = "/") %>% paste0("/", dataset_name, ".rds")
    resultpath <- ap$paths$liana_resdir  %>% paste0(collapse = "/")
    figpath <- ap$paths$liana_figdir %>% paste0(collapse = "/")
    aggpath <- ap$paths$liana_aggdir %>% paste0(collapse = "/")

    lig_rec <- if(is.null(ap[["liana"]][["lig_rec"]])) NA
    else{
      sapply(ap[["liana"]][["lig_rec"]], function(elem){
        if(length(elem) <= 1){
          elem %>% unlist
        }else{
          list(elem[[1]] %>% unlist, elem[[2]] %>% unlist)
        }
      })
    }
    sapply(c(resultpath, figpath, aggpath), function(path) dir.create(path, recursive = TRUE))
    list(dataset_name = dataset_name, datapath = datapath, resultpath = resultpath, figpath = figpath, aggpath = aggpath, base = eval(parse(text = 'ap[["liana"]][["base"]]')), methods = ap[["liana"]][["methods"]] %>% unlist, lig_rec = lig_rec)
  }

  tool_config <- do.call(rbind, lapply(names(datasets), function(dataset_name) create_tool_config(dataset_name, datasetpath, toolname)))
  write.table(tool_config, paste("../../", "config", paste0(toolname, "_config.tsv"), sep = "/"), quote=FALSE, sep='\t', row.names = FALSE)
  tool_config
}



#' create_analysis_obj
#' @description maps dataset name, data and params, reads dataset specific analysis params
#' @param sces SingleCellExperiment objects
#' @param datasetpath path to dataset
#' @return list of lists, each representing a dataset with the attributes data, params and config
#' @export
create_analysis_obj <- function(datasets, datasetpath){
  # Create Config
  config <- save_tool_config(datasets, datasetpath, "liana")
  rownames(config) <- config[, "dataset_name"] %>% unlist

  # To SCE
  sces <- sapply(datasets, function(dataset){
    if(class(dataset) == "Seurat") Seurat::as.SingleCellExperiment(dataset)
    else dataset
  })

  analysis <- lapply(names(sces), function(name){
    list("data" = sces[[name]], "params" = read_yaml(paste(datasetpath, name, "analysis_params.yaml", sep = "/")), "config" = config[name,])
  })
  names(analysis) <- names(sces)
  analysis
}


#' rename_all uses analysis params to rename metadata columns
#' @description for each dataset the metadata columns are renamed according to the information in the analysis_params
#' @param sces a list of singleCellExperiment objects
#' @param analysis a list of additional information for each dataset in sces (analysis_params and paths)
#' @export
rename_all <- function(sces, analysis){
  # rename_dataset
  # @description recursive function, each function call renames one metadata column
  # @reasoning using sapply didn't work therefore using two recursive functions
  # @param cols list of column names that shall be renamed
  # @param meta metadata table
  # @param analysis_ds the analysis object of a specific dataset
  rename_dataset <- function(cols, meta, analysis_ds){
    if(length(cols) == 0){
      meta %>% DataFrame # conversion needed for colData slot
    } else {
      col <- cols[1]
      meta <- meta %>% dplyr::rename({{ col }} := analysis_ds$analysis_params$rename[[col]])
      rename_dataset(cols[-1], meta, analysis_ds)
    }
  }
  # for each dataset do rename
  dsnames <- names(sces)
  rename <- function(dsnames, sces){
    if(length(dsnames) == 0) {sces}
    else{
      cols <- names(analysis[[dsnames[1]]]$analysis_params$rename)
      meta <- sces[[dsnames[1]]]@colData %>% data.frame
      colData(sces[[dsnames[1]]]) <- rename_dataset(cols, meta, analysis[[dsnames[1]]])
      rename(dsnames[-1], sces)
    }
  }
  rename(dsnames, sces)
}



#' @title basic
#' @author Thalla
#' @description a package is loaded or installed and loaded
#' @param package the name of the package
#' @param fromGithub if the package shall be installed from Github
#' @param fromBioconductor if the package shall be installed from Bioconductor
#' @examples loadPackage("base")
#' @examples loadPackage("ying14/yingtools2", 1)
#' @return package is loaded
#' @keywords load package install
#' @export
loadPackage <- function(package, fromGithub = F, fromBioconductor = F){
  package <- gsub(".*/(.*)","\\1", package)
  if(!require(package, character.only = T)){
    if(!fromGithub){
      if(!fromBioconductor){
        install.packages(package)
      }
      else{
        if (!requireNamespace("BiocManager", quietly = TRUE)){
          install.packages("BiocManager")
        }
        BiocManager::install(package)
      }
    }
    else{
      if(!require("devtools", character.only = T)){
        install.packages("devtools")
      }
      devtools::install_github(package)
    }
  }
  library(package, character.only = T)
}


#' @title save_tool_config
#' @author Thalla
#' @description Creates and saves a config file for a specific tool. Each dataset is added with the relevant paths. File is saved in the version folder.
#' @param datasets list of datasets
#' @param toolname name of tool
#' @keywords config init
#' @return config table
#' @export
save_tool_config <- function(datasets, datasetpath, toolname){
  # create_tool_config reads dataset specific analysis_params file, creates data path and specific for liana result & figure path
  # @param dataset_name name of dataset
  # @param datasetpath path of dataset
  # @return dataframe with name, datapath, resultpath, figpath
  create_tool_config <- function(dataset_name, datasetpath, toolname){
    ap <- yaml::read_yaml(paste(datasetpath, dataset_name, "analysis_params.yaml", sep = "/"))
    datapath <- paste0(ap$paths$datapath_tmp) %>% paste0(collapse = "/") %>% paste0("/", dataset_name, ".rds")
    resultpath <- ap$paths$liana_resdir  %>% paste0(collapse = "/")
    figpath <- ap$paths$liana_figdir %>% paste0(collapse = "/")
    aggpath <- ap$paths$liana_aggdir %>% paste0(collapse = "/")

    lig_rec <- if(is.null(ap[["liana"]][["lig_rec"]])) NA
    else{
      sapply(ap[["liana"]][["lig_rec"]], function(elem){
        if(length(elem) <= 1){
          elem %>% unlist
        }else{
          list(elem[[1]] %>% unlist, elem[[2]] %>% unlist)
        }
      })
    }
    sapply(c(resultpath, figpath, aggpath), function(path) dir.create(path, recursive = TRUE))
    list(dataset_name = dataset_name, datapath = datapath, resultpath = resultpath, figpath = figpath, aggpath = aggpath, base = eval(parse(text = 'ap[["liana"]][["base"]]')), methods = ap[["liana"]][["methods"]] %>% unlist, lig_rec = lig_rec)
  }

  tool_config <- do.call(rbind, lapply(names(datasets), function(dataset_name) create_tool_config(dataset_name, datasetpath, toolname)))
  write.table(tool_config, paste(getwd(), paste0(toolname, "_config.tsv"), sep = "/"), quote=FALSE, sep='\t', row.names = FALSE)
  tool_config
}



#' create_analysis_obj
#' @description maps dataset name, data and params, reads dataset specific analysis params
#' @param sces SingleCellExperiment objects
#' @param datasetpath path to dataset
#' @return list of lists, each representing a dataset with the attributes data, params and config
#' @export
create_analysis_obj <- function(datasets, datasetpath){
  # Create Config
  config <- save_tool_config(datasets, datasetpath, "liana")
  rownames(config) <- config[, "dataset_name"] %>% unlist

  # To SCE
  sces <- sapply(datasets, function(dataset){
    if(class(dataset) == "Seurat") Seurat::as.SingleCellExperiment(dataset)
    else dataset
  })

  analysis <- lapply(names(sces), function(name){
    list("data" = sces[[name]], "params" = read_yaml(paste(datasetpath, name, "analysis_params.yaml", sep = "/")), "config" = config[name,])
  })
  names(analysis) <- names(sces)
  analysis
}


