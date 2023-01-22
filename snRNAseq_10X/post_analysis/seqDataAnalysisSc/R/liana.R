#' run_liana_wrap
#' @param dataset element of list 'analysis'
#' @param cluster_col_name name of the column that defines the cell clusters
#' @param new (boolean), if true (re)run liana_wrap
#' @export
run_liana_wrap <- function(dataset, cluster_col_name, new){
  res <- NULL
  res_dir <- paste(dataset$config$resultpath, cluster_col_name,
                   sep = "/")
  methods <- dataset$config$methods
  is_existing <- sapply(methods, function(methodname) {
    file.exists(paste0(res_dir, "/", methodname, ".csv"))
  }) %>% sum %>% equals(length(methods))
  if (new | !is_existing) {
    dir.create(res_dir, recursive = TRUE)
    message("....Calculating....")
    res <- dataset$data %>% liana_wrap(method = methods,
                                       parallelize = TRUE, base = eval(parse(text = dataset$config$base)),
                                       idents_col = cluster_col_name)
    message("....Writing.... -")
    sapply(names(res), function(methodname) {
      write.csv(res[[methodname]], paste0(res_dir, "/",
                                          methodname, ".csv"))
    })
    message("Saved: Method Results ")
  }
  else {
    message("- Info - Results already exist. Please read them in with 'read_liana_methods'.")
  }
  res
}

#' get_liana_agg
#' @param res result from run_liana_wrap
#' @param dataset element of 'analysis' object
#' @param new <boolean> if true (re)run aggregation else read previous results in
#' @export
get_liana_agg <- function(res, dataset, new = FALSE){
  agg <- tibble::tibble()
  agg_dir <- paste0(dataset$config$resultpath, "/aggregated/")
  agg_file <- paste0(agg_dir, "aggregated_results.csv")
  dir.create(agg_dir)

  if(new | !file.exists(agg_file)){
    agg <- res %>% liana_aggregate()
    write.csv(agg, agg_file)
    message("Aggregation results are saved.")
  }else{
    message("....Reading preexisting results....")
    agg <- read.csv(agg_file) %>% tibble::as_tibble()
  }
  agg
}


#' @title run_liana_plots
#' @author Thalla
#' @description plot aggregation results and save the plots, loops over ltes, ntops, lig_rec
#' @param agg aggregated results from LIANA
#' @param dataset element of 'analysis', has data, config, ...
#' @param cluster_col_name name of clustering column that was used with LIANA so far, if a list is provided loops over list
#' @param ltes filter aggregate_rank according to lower than or equals the given value, if a list is provided loops over list
#' @param ntops number of top n results chosen for dotplot, if a list is provided loops over list
#' @usage run_liana_plots(agg, analysis$all, "seurat_clusters", lte = c(0.01, 0.05))
#' @export
run_liana_plots <- function(agg, dataset, cluster_col_name, ltes = c(0.01), ntops = c(25)){
  sapply(ltes, function(lte){
    sapply(ntops, function(ntop){
      lig_recs <- dataset$config$lig_rec
      sapply(lig_recs, function(lig_rec_elem){
        plot_lig_rec <- function(lig_rec, file_suffix = NULL){
          if(is.null(file_suffix)){
            file_suffix <- paste0("_", paste(lig_rec[[1]], collapse = ''), "_", paste(lig_rec[[2]], collapse = ''))
          }
          c(source_groups, target_groups) %<-% list(lig_rec[[1]], lig_rec[[2]])
          # Paths
          figdir <- paste0(dataset$config$figpath, "/", cluster_col_name)
          resultdir <- paste0(dataset$config$resultpath, "/", cluster_col_name)
          dir.create(figdir)
          # @function
          plot_agg <- function(agg, file_suffix){
            # depending on withFDR with correction, lte filter, dotplot, heatmap
            plot_agg_adj <- function(agg, withFDR, file_suffix){
              if(withFDR){
                agg %<>% dplyr::mutate(aggregate_rank = p.adjust(aggregate_rank, method="BH"))
                file_suffix <- paste0(file_suffix , "_bh")
              }
              agg %<>% dplyr::filter(aggregate_rank <= lte)
              file_suffix <- paste0(file_suffix, "_lte", lte) %>% stringr::str_replace("\\.", "")
              # @function dotplot
              plot_agg_dotplot <- function(agg, file_suffix){
                agg %<>% dplyr::filter(source %in% source_groups) %>% dplyr::filter(target %in% target_groups)
                if(! agg %>% nrow %>% equals(0)){
                  filename <- paste0("/dotplot", file_suffix, "_n", ntop)
                  pdf(file = paste0(figdir, filename, ".pdf"),  width = 12, height = 10)
                  agg %>% liana_dotplot(source_groups = source_groups, target_groups = target_groups, ntop = ntop) %>% plot
                  dev.off()
                  write.csv(agg, paste0(resultdir, filename, ".csv"))
                  message("-Saved-: Dotplot")
                }else{
                  message("-Info-: Considering the source and target groups there are no results to plot.")
                }
              }
              # @function heatmap
              plot_agg_heatmap <- function(agg, file_suffix){
                aggnew <- agg
                repeat{
                  agg <- aggnew
                  aggnew <- agg %>% dplyr::filter((source %in% unique(.$target)) & (target %in% unique(.$source)))
                  if(nrow(aggnew) == nrow(agg)){ break }
                }
                if(nrow(aggnew) == 0) return
                filename <- paste0("/heatmap", file_suffix)
                if(length(unique(agg$source)) >= 2 & length(unique(agg$target)) >= 2){
                  pdf(paste0(figdir, filename, ".pdf"),  width = 7, height = 6)
                  agg %>% heat_freq %>% plot
                  dev.off()
                } else message("-Info-: Not enough distinct values. No heatmap is generated.")
                write.csv(agg, paste0(resultdir, filename, ".csv"))
                message("-Saved-: Heatmap")
              }
              plot_agg_dotplot(agg, file_suffix)
              plot_agg_heatmap(agg, file_suffix)
            }
            plot_agg_adj(agg, 0, file_suffix)
            plot_agg_adj(agg, 1, file_suffix)
          }
          plot_agg(agg, file_suffix)
        }
        if(length(lig_rec_elem) == 1){ # lig and rec for the same clusters (lig_rec length is 1 (NA))
          lig_rec_levels <- levels(droplevels(factor(agg$source)))
          plot_lig_rec(list(lig_rec_levels, lig_rec_levels), "_all")
        }else{ # lig cluster and rec cluster are different
          c(e1, e2) %<-% list(lig_rec_elem[[1]] %>% droplevels, lig_rec_elem[[2]] %>% droplevels)
          plot_lig_rec(list(e1,e2))
          plot_lig_rec(list(e2, e1))
        }
      })
    })
  })
}


#' @description Heatmap showing the number of significant interactions between clusters
#' @title src_trgt_heatmap
#' @author Thalla
#' @param agg_plot aggregation data, eventually preprocessed for plotting
#' @param title title of plot
#' @param nrow for filterint table matrix: table(source, target)[1:nrow, 1:ncol]
#' @param ncol see nrow
#' @export
src_trgt_heatmap <- function(agg_plot, title, nrow, ncol, figpath){
  #coul <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "PuBu"))(25)
  getplot <- function(){
    #    heatmap(col = coul, Colv = NA, Rowv = NA, scale="column", main = title, xlab = "target", ylab = "source")
    #  legend(x="bottomright", legend=c("min", "ave", "max"), fill=grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "PuBu"))(3))
    setHook("grid.newpage", function() grid::pushViewport(grid::viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
    table(agg_plot$source, agg_plot$target)[1:nrow, 1:ncol] %>% # table(row, col)
      as.matrix %>%
      pheatmap::pheatmap(display_numbers = T, cluster_rows = FALSE, cluster_cols = FALSE, main = title)
    setHook("grid.newpage", NULL, "replace")
    grid::grid.text("target", y=-0.07, gp=grid::gpar(fontsize=16))
    grid::grid.text("source", x=-0.07, rot=90, gp=grid::gpar(fontsize=16))
  }
  pdf(file=paste0(figpath, "/", title, ".pdf"), width = 8, height = 6)
  getplot()
  dev.off()
  getplot()
}
