# module load R/4.0.3-foss-2020a
require(SingleCellExperiment)
library(Seurat)
library(celda)
library(stringr)

# input data: Seurat object

args =  commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
    inFile = args[1]
}


print(paste("Input data:", inFile))

vals <- str_split(inFile, "/")[[1]]
sId =   gsub("_object.Rdata","", vals[length(vals)])

print(paste("Load sample:", sId))
tumObj <- readRDS(inFile)

counts_tumor = tumObj@assays$RNA@counts
input_tumor <- SingleCellExperiment(assays = list(counts = counts_tumor))

decont_tumor <- decontX(x = input_tumor)

print("Draw figures...")

# draw initla result
umap <- reducedDim(decont_tumor, "decontX_UMAP")
fName1 <- gsub(".Rdata",".decontX_check.pdf",inFile)
pdf(fName1)
plotDimReduceCluster(x = decont_tumor$decontX_clusters,
                     dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXContamination(decont_tumor)
dev.off()



# draw via coputed initially UMAP
tumObj@meta.data$decontX_level <- decont_tumor$decontX_contamination

fName2= gsub("_object.Rdata",paste0("_UMAP.decontX_res.pdf"),inFile)
pdf(fName2,width = 8, height = 6)
FeaturePlot(tumObj, reduction = "umap",features = "decontX_level")
dev.off()


# save full result
print("Saving...")

resName <- gsub(".Rdata",".decontX.Rdata",inFile) 
saveRDS(decont_tumor,resName)

resName2 <- gsub(".Rdata",".decontX_ann.txt",inFile)
resDf <- data.frame(decontX_res = decont_tumor$decontX_contamination)
rownames(resDf) <- colnames(decont_tumor)
write.table(resDf, resName2,sep="\t",quote=F)

print("Finished!")
