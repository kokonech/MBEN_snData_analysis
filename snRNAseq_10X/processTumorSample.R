# requirements
# for Seurat: R-bundle/20180906-foss-2017a-R-3.5.1
# for UMAP:  Python/2.7.12-foss-2017a

library(Seurat)
library(dplyr)
options(stringsAsFactors=FALSE)

### FUNCTIONS

extractClusterMatrix <- function(inData, id, highlyVar = F) {
    ctypes <- inData@meta.data[,id]
    resTable <- NULL
    for (k in sort(unique(ctypes))) {
        cNames <- rownames(inData@meta.data[ ctypes == k, ] )
        print(paste("cluster:",k,"num cells:",length(cNames)))
        vals <- rowSums(inData@raw.data[,cNames])
        resTable <- cbind(resTable,vals)
    }
    colnames(resTable) <- paste0("c",sort(unique(ctypes)))
    resTable
}


geneLengths = read.table("/b06x-isilon/b06x-m/mbCSF/annotation/hg19/hg19_genes.length", header=1)

convToRpkm <- function(clusterMatrix, geneLengths) {
    ns=colSums(clusterMatrix)
    rpkmTable=NULL
    # issue with different genes
    commonGenes <- intersect(rownames(clusterMatrix) , rownames(geneLengths))
    geneLengths <- geneLengths[commonGenes, ,drop=FALSE]
    clusterMatrix <- clusterMatrix[commonGenes,]

    geneLengths <- geneLengths[row.names(clusterMatrix), ,drop=FALSE]
    for(i in 1:ncol(clusterMatrix)){
        rpkmTable = cbind(rpkmTable,
        round((clusterMatrix[,i]*10^9)/(geneLengths$Length*ns[i]), digits=3))
    }
    colnames(rpkmTable) = colnames(clusterMatrix)
    rownames(rpkmTable) = rownames(clusterMatrix)
    rpkmTable
}


##### MAIN ANALYSIS

args =  commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
    sId = args[1]
    dataPath = args[2]
}

print(paste("Sample ID:", sId))
print(paste("Input data:", dataPath))

print("Read initial data...")
cb.data <- Read10X( data.dir = dataPath )

timeStr = format(Sys.time(), "%Y%m%d_%H%M")
mainDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/perSampleAnalysis/"
resDir = paste0(mainDir,sId,"_res.", timeStr, "/")

if (!file.exists(resDir)){
    dir.create(resDir)
}

ndim = 10

print("Start analysis...")
colnames(cb.data) <- paste0(sId,"_c",1:ncol(cb.data))

cb <- CreateSeuratObject(cb.data, min.cells = 5, min.genes = 300, project = sId)

mito.genes <- grep("^MT", rownames(cb@data), value = T)
percent.mito <- Matrix::colSums(cb@raw.data[mito.genes, ])/Matrix::colSums(cb@raw.data)
cb <- AddMetaData(cb, percent.mito, "percent.mito")

print("QC report...")
pdf(paste0(resDir, sId,"_premrna_VlnPlot.pdf"), width=10, height=5)
VlnPlot(cb, c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()


# Filter
# Max limits gene, transcript for Single Cell v2: 4000,8000
cb <- FilterCells(object = cb, subset.names = c("nGene","nUMI",  "percent.mito"), low.thresholds = c(300,300,  -Inf), high.thresholds = c(4000, 8000, 0.05))


print("Normalization...")
cb <- NormalizeData(object = cb,  normalization.method = "LogNormalize",
                    scale.factor = 10000 )

cb <-FindVariableGenes(cb ,mean.function = ExpMean, dispersion.function =LogVMR,  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(cb@var.genes) # check num genes

cb <- ScaleData(  object = cb,   vars.to.regress = c("nUMI","percent.mito"))


print("Dimensioional reduction...")
cb <- RunPCA(cb, pc.genes = cb@var.genes, do.print = F, pcs.print = 1:5 )

png(paste0(resDir,sId,"_PcElbow.png"), width=600, height=400)
PCElbowPlot(cb, num.pc = 20)
dev.off()

# default: 0.75
cb <- FindClusters(cb, dims.use = 1:ndim, resolution = 0.75, print.output = F, save.SNN = T,force.recalc = T)

print("Check clusters...")
table(cb@ident)


print("tSNE...")
cb <- RunTSNE(cb, dims.use = 1:ndim, do.fast = T)
pdf(paste0(resDir,sId,"_tSNE.pdf"))
TSNEPlot(cb, do.label = T)
dev.off()

print("UMAP...")
cb <- RunUMAP(object = cb, dims.use = 1:ndim)
pdf(paste0(resDir,sId,"_UMAP.pdf"),width = 8, height = 6)
DimPlot(object = cb, reduction.use = 'umap')
dev.off()

print("Find DEGs...")
cb.markers <- FindAllMarkers(cb, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

write.table(cb.markers, paste0(resDir, sId,"_clusters_marker_genes.txt") , sep = "\t", quote = FALSE)

top10 <- cb.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

png(paste0(resDir,sId, "_top10_per_cluster_heatmap.png"), width=1200, height=1000)
DoHeatmap(object = cb, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()


print("Save result...")
saveRDS(cb, paste0(resDir,sId,"_object.Rdata"))


print("Save expression...")
clCounts <- extractClusterMatrix(cb,5)
write.table(clCounts, paste0(resDir,sId,"_cluster_expr_counts.txt"),sep="\t",quote = F)

rpkmData <- convToRpkm(clCounts,geneLengths)
write.table(rpkmData, paste0(resDir,sId,"_cluster_expr_RPKM.txt"),sep="\t",quote = F)

print("Done!")




