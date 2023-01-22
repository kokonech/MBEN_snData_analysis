# requirements
# for Seurat: module load R/3.6.0-foss-2017a
# for UMAP:  Python/2.7.12-foss-2017a

library(Seurat)
library(dplyr)
library(config)
options(stringsAsFactors=FALSE)

### FUNCTIONS

extractClusterMatrix <- function(rawTable, annData ) {
    resTable <- NULL
    clusters <- sort(unique(as.numeric(annData$seurat_clusters)) )
    print(clusters)
    for (k in clusters) {
        cNames <- rownames(annData)[as.numeric(annData$seurat_clusters) == k ]
        print(paste("cluster:",k,"num cells:",length(cNames)))
        vals <- rowSums(as.matrix(rawTable[,cNames]))
        resTable <- cbind(resTable,vals)
    }
    colnames(resTable) <- paste0("c",clusters -1) 
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


print("Read configruation:")
cfgPath = paste0("/b06x-isilon/b06x-m/mbCSF/results/smSeq2Res/configs/",sId,"_config.txt")

if (!file.exists(cfgPath)) {
    print("Using default config...")
    cfgPath = "/b06x-isilon/b06x-m/mbCSF/results/smSeq2Res/configs/default_config.txt"
}

print(cfgPath)

config <- config::get(file=cfgPath)


print("Read initial data...")
mols <- read.table(dataPath, sep = "\t") 

timeStr = format(Sys.time(), "%Y%m%d")
mainDir = "/b06x-isilon/b06x-m/mbCSF/results/smSeq2Res/perSampleAnalysisSV3/"

resDir = paste0(mainDir,sId,"_S3_res.", timeStr, "/")

if (!file.exists(resDir)){
    dir.create(resDir)
}

ndim = 10

print("Start analysis...")
#colnames(mols) <- paste0(sId,"_c",1:ncol(mol)) 

cb <- CreateSeuratObject(mols, min.cells = 5, min.features = 300, project = sId)

# save for raw counts extraction

#rawData <- cb@assays$RNA

cb[["percent.mito"]] <- PercentageFeatureSet(cb, pattern = "^MT-")


print("QC report...")
pdf(paste0(resDir, sId,"_premrna_VlnPlot.pdf"), width=10, height=5)
VlnPlot(cb, c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

print("Filter...")

cb <- subset(cb, subset = nFeature_RNA > config$min_genes  & nFeature_RNA < config$max_genes  & percent.mito < config$min_mt_perc )


print("Normalization...")
cb <- NormalizeData(cb, normalization.method = "LogNormalize", scale.factor = 10000)

# control the amount
# Q: how many genes to select?
# intial 2000, default applied 2500

cb <- FindVariableFeatures(cb, selection.method = "vst", nfeatures = 2500)

top10 <- head(VariableFeatures(cb), 10)

plot1 <- VariableFeaturePlot(cb)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(paste0(resDir, sId,"_top10_variable.pdf"), width=12, height=6 )
CombinePlots(plots = list(plot1, plot2))
dev.off()


# Q: regress variables as previously?
# possible, but novel method is suggested
all.genes <- rownames(cb)
cb <- ScaleData( cb, features = all.genes )


print("Dimensioional reduction...")
cb <- RunPCA(cb, features = VariableFeatures(object = cb))


png(paste0(resDir,sId,"_PcElbow.png"), width=600, height=400)
ElbowPlot(cb, ndims = 20)
dev.off()


print("Find clusters...")

cb <- FindNeighbors(cb, dims = 1:ndim)
cb <- FindClusters(cb, resolution = 0.5)
table(Idents(cb))


print("tSNE...")
cb <- RunTSNE(cb, dims.use = 1:ndim)
pdf(paste0(resDir,sId,"_tSNE.pdf"))
TSNEPlot(cb, pt.size = 1, label = T)
dev.off()

print("UMAP...")
cb <- RunUMAP(object = cb, dims = 1:ndim)
pdf(paste0(resDir,sId,"_UMAP.pdf"),width = 8, height = 6)
DimPlot(object = cb, pt.size=1, reduction = 'umap')
dev.off()

print("Save result...")
saveRDS(cb, paste0(resDir,sId,"_object.Rdata"))


gNames <- c("EOMES", "LMX1A", "MYC", "MYCN", "ELP4", "IMMP1L", "PAX6", "SNCAIP", "PRDM6")

gNames <- gNames[gNames %in% rownames(cb)]

pdf(paste0(resDir, sId,"_UMAP.gene_info.pdf"), width=8, height=6)
for (gName in gNames) {
  print(gName)
  print(FeaturePlot(cb, features = gName, reduction = "umap"))
}
dev.off()

print("Find DEGs...")
cb.markers <- FindAllMarkers(cb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(cb.markers, paste0(resDir, sId,"_clusters_marker_genes.txt") , sep = "\t", quote = FALSE)



if (FALSE) {


# issue with infity value
top10 <- cb.markers %>% group_by(cluster) %>% top_n(n=10 , wt=avg_logFC)

png(paste0(resDir,sId, "_top10_per_cluster_heatmap.png"), width=1200, height=1000)
DoHeatmap(object = cb, features = top10$gene) + NoLegend()
dev.off()
print("Save pseudobulk expression...")
clCounts <- extractClusterMatrix(rawData, cb@meta.data  )
write.table(clCounts, paste0(resDir,sId,"_cluster_expr_counts.txt"),sep="\t",quote = F)

rpkmData <- convToRpkm(clCounts,geneLengths)
write.table(rpkmData, paste0(resDir,sId,"_cluster_expr_RPKM.txt"),sep="\t",quote = F)

}
print("Done!")




