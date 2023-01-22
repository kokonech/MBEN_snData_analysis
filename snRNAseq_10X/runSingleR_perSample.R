# module load R/4.0.1-foss-2017a


require(SingleCellExperiment)
library(SingleR)
library(Seurat)
library(stringr)
library(patchwork)

# input data: seurat object
# example
#inFile = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/perSampleAnalysisSV3/MBEN_AK3_S3_res.20191008_1726/MBEN_AK3_object.Rdata"


args =  commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
    inFile = args[1]
}


if (length(args) > 1) {
    ctrlFile = args[2]
    ctrlName = args[3]
} else {
    # default
    ctrlFile = "/b06x-isilon/b06x-m/mbCSF/results/humanCbData/extFetal/Cerebellum_cell_types_pseudobulk.gn.log2.txt"
    ctrlName = "cerebellum_fetal" 
}

print(paste("Input data:", inFile))
print(paste("Control data :",ctrlFile))




vals <- str_split(inFile, "/")[[1]]
sId =   gsub("_object.Rdata","", vals[length(vals)])

print(paste("Load sample:", sId))

sObj <- readRDS(inFile)

# normalized input tumor data
# 10X data
counts_tumor = sObj@assays$RNA@data
input_tumor <- SingleCellExperiment(assays = list(logcounts = counts_tumor))



print(paste("Compare to ", ctrlName,  "...") )

# assumme log2 adjustement
logCbCtrl <- read.delim(ctrlFile)


cb_ref <- SummarizedExperiment( assays=list(logcounts = logCbCtrl) )
colData(cb_ref)$celltypes <-  colnames(logCbCtrl)

# analysis

pred.cb <- SingleR(test=input_tumor, ref = cb_ref , assay.type.test=1,
    labels = cb_ref$celltypes )

resName1 = gsub("_object.Rdata", paste0("_singleR_",ctrlName,"_check.pdf"),inFile) 
pdf(resName1,width=8,height=6)
plotScoreHeatmap(pred.cb)
plotDeltaDistribution(pred.cb, ncol = 3)
dev.off()

sObj@meta.data$cell_type <- pred.cb$labels
sObj@meta.data$filtered_ct <- pred.cb$pruned.labels

resName2= gsub("_object.Rdata",paste0("_UMAP.",ctrlName,"_SingleR.pdf"),inFile)
#pdf(resName2,width = 14, height = 6)
#Seurat::DimPlot(object = sObj,reduction = 'umap', group.by=c("cell_type","filtered_ct")) + plot_annotation(paste("SingleR: comparison to",ctrlName))
pdf(resName2,width = 8, height = 6)
Seurat::DimPlot(object = sObj,reduction = 'umap', group.by=c("filtered_ct")) + plot_annotation(paste("SingleR: comparison to",ctrlName))
dev.off()

#resName2b = gsub("_object.Rdata",paste0("_tSNE.",ctrlName,"_SingleR.pdf"),inFile)
#pdf(resName2b,width = 8, height = 6)
#Seurat::DimPlot(object = sObj,reduction = 'tsne', group.by=c("filtered_ct")) + plot_annotation(paste("SingleR: comparison to",ctrlName))
#dev.off()

resName3 = gsub("_object.Rdata",paste0("_SingleR_",ctrlName,"_ann.txt"),inFile)
write.table(pred.cb, resName3, sep="\t", quote=F)

print("Finished!")
