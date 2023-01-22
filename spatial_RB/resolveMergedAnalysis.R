library(Seurat)
library(dplyr)
setwd("~/work/spatialMB/")

dataDir = "/Users/okonechin/work/spatialMB/MB_spatial_results/"
resDir = "/Users/okonechin//work/spatialMB/MB_spatial_combined/"
# load samples list

grpName = "MBEN"
sIds = c("MB263", "MB266", "MB295", "MB299")

ob.list <- list()
orig.ids <- c()
for (sId in sIds) {
    sPath = paste0(dataDir,sId,".res.rds")
    print(sId)
    ob.list[[sId]] <- readRDS(sPath)
    orig.ids <- c(orig.ids,rep(sId,ncol(ob.list[[sId]])))
}

mb.comb <- merge(ob.list[[1]],
                y = ob.list[2: length(sIds)],
                add.cell.ids = sIds,
                project = grpName)
mb.comb@meta.data$sample <- orig.ids

# clean
mb.comb$seurat_clusters <- NULL


# method #1: combine samples with batch effect adjustment
resName = paste0(grpName,"_comb")
require(harmony)


# default ncells = 5K
mb.comb <- SCTransform(mb.comb, assay = "Spatial", verbose = FALSE,ncells = 10000)
mb.comb <- RunPCA(mb.comb, assay = "SCT", verbose = FALSE,npcs = 30, approx=FALSE)

png(paste0(resDir,resName,".PCA_elbow.png"),width=600, height = 400)
ElbowPlot(mb.comb)
dev.off()

summary(as.factor(mb.comb@meta.data$sample))
mb.comb <- harmony::RunHarmony(mb.comb, dims=1:20, 
                              group.by.vars = "sample",assay.use="SCT")

mb.comb <- FindNeighbors(mb.comb, reduction = "harmony", dims = 1:20)
mb.comb <- FindClusters(mb.comb, resolution = 0.4)
mb.comb <- RunUMAP(mb.comb, reduction = "harmony", dims = 1:20)


pdf(paste0(resDir, resName,".UMAP.pdf"), width=6, height = 4)
DimPlot(mb.comb, reduction = "umap", label = TRUE)
DimPlot(mb.comb, reduction = "umap", group.by = "sample")
dev.off()


DimPlot(mb.comb, reduction = "umap", label = TRUE, #split.by = "seurat_clusters")
        cells = mb.comb@meta.data$seurat_clusters == "22" )

gNames <- rownames(mb.comb)

pdf(paste0(resDir, resName,".UMAP.gene_expr.pdf"), width=6, height = 4)

for (gName in gNames) {
    print(gName)
    print(FeaturePlot(mb.comb, features = gName, reduction = "umap"))
}
dev.off()


pdf(paste0(resDir, resName,".gene_vln_profiles.pdf"), width=12, height = 4)
for (gName in gNames) {
    print(gName)
    print(VlnPlot(mb.comb, features = gName, combine = TRUE,pt.size=0))
}
dev.off()

markers <- FindAllMarkers(mb.comb, only.pos = TRUE)
fName = paste0(resDir, resName,".DEG_table.txt")
write.table(markers, fName, sep="\t", quote=F)

markers2 <- markers %>% group_by(cluster)

pdf(paste0(resDir, resName,".DEG_heatmap.pdf"),width=16,height = 14)
DoHeatmap(mb.comb,features = markers2$gene)
dev.off()

# inspect contamation
require(SingleCellExperiment)
require(celda)
counts_tumor = mb.comb@assays$Spatial@counts
input_tumor <- SingleCellExperiment(assays = list(counts = counts_tumor))

decont_tumor <- decontX(x = input_tumor)


saveRDS(mb.comb, paste0(resDir,resName,".object.rds"))


# markers visualization

ngList = c("NRXN3","RBFOX3","CNTN2",  # NOD
        "LAMA2", "GLI1", "TRPM3", "IGFBP7", "COL1A2", # INTER
        "AIF1", "ITGAM", "CD14", "CD163", # immune
        "AQP4", "NPAS3", "TNC", # astro
        "TULP1","MKI67", "FCGR3A") # other


fName <- paste0(resDir, resName, ".DEG_dotplot.pdf")
pdf(fName,width = 8, height = 6)
DotPlot(mb.comb, features = ngList) + RotatedAxis()
dev.off()

# rename for plotting
initClusters <- mb.comb@active.ident

# adjust
renamedClusters <- initClusters 


levels(renamedClusters) <- c("Diff_neuronal-like",
                             "Migr_CGNP-like", "Stromal",
                             "Early_CGNP-like_prolif","Early_CGNP-like",
                             "Late_GCNP-like","Vascular",
                             "Immune", "Astrocytes","Migr_CGNP-like",
                             "Diff_neuronal-like","Prolif","Immune")


mb.comb@active.ident <- renamedClusters


# manual colors control for MBEN
    
clColors <- c("#F8766D", "#00BA38","purple", # diff, Migr, stromal
              "#00BFC4", "#B79F00","coral4", # early profli, early, later
              "magenta1","powderblue", "#619CFF", "palegreen") # vascular, immune, astro, prolif

names(clColors) <- levels(mb.comb@active.ident)

pdf(paste0(resDir, resName,".UMAP_renamed.pdf"), width=10, height = 8)
#png(paste0(resDir, resName,".UMAP_renamed.v3.png"), width=1000, height = 800)
DimPlot(mb.comb, reduction = "umap", label = FALSE, cols = clColors)

dev.off()

fName <- paste0(resDir, resName, ".DEG_heatmap_renamed.pdf")
pdf(fName,width = 8, height = 6)
DotPlot(mb.comb, features = ngList) + RotatedAxis()
dev.off()

# save with assigned names
mb.comb@meta.data$assigned_cluster <- renamedClusters
saveRDS(mb.comb, paste0(resDir,resName,".object_renamed.rds"))

# source: https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html
require(ggplot2)

propDf <- NULL
cT <- summary(as.factor(mb.comb@meta.data$assigned_cluster))
for (cName in names(cT) ) {
    targAnn <- mb.comb@meta.data[ mb.comb@meta.data$assigned_cluster == cName,]
    vals <- summary(as.factor(targAnn$sample))
    #names(vals) <- gsub("main", "MB165", names(vals))
    clDf <- data.frame(type=cName ,sId = names(vals), value = vals)
    propDf <- rbind(propDf,clDf)
}

propDf$type <- factor(propDf$type,levels = names(cT))

ggplot(propDf, aes(fill=sId, y=value, x=type)) + 
    geom_bar(position="fill", stat="identity") +
    theme_minimal() +
    xlab("") + ylab("Sample cells prop")


# Method 2 : no batch effect adjustment



resName = paste0(grpName,"_merged")

# normalize?
mb.merged <- SCTransform(mb.comb, assay = "Spatial", verbose = FALSE)

mb.merged <- RunPCA(mb.merged, assay = "SCT", verbose = FALSE,npcs = 30, approx=FALSE)
ElbowPlot(mb.merged)
mb.merged <- FindNeighbors(mb.merged, reduction = "pca", dims = 1:20)
mb.merged <- FindClusters(mb.merged, resolution = 0.4)
mb.merged <- RunUMAP(mb.merged, reduction = "pca", dims = 1:20)

pdf(paste0(resDir, resName,".UMAP.pdf"), width=6, height = 4)
DimPlot(mb.merged, reduction = "umap", label = TRUE)
DimPlot(mb.merged, reduction = "umap",group.by = "sample")
dev.off()

gNames <- rownames(mb.merged)

for (gName in gNames) {
    print(gName)
    png(paste0(resDir, resName,".UMAP.",gName,".png"), width=600, height = 400)
    print(FeaturePlot(mb.merged, features = gName, reduction = "umap"))
    dev.off()
}

markers <- FindAllMarkers(mb.merged, only.pos = TRUE)
fName = paste0(resDir, resName,".DEG_table.txt")
write.table(markers, fName, sep="\t", quote=F)

markers2 <- markers %>% group_by(cluster)

pdf(paste0(resDir, resName,".DEG_heatmap.pdf"),width=16,height = 14)
DoHeatmap(mb.merged,features = markers2$gene)
dev.off()

nod = c("NRXN3","RBFOX3","CNTN2")
inter = c("LAMA2", "GLI1", "TRPM3", "IGFBP7", "COL1A2")
immune = c("AIF1", "ITGAM", "CD14", "CD163")
astro = c("AQP4", "NPAS3", "TNC")
other = c("TULP1","FCGR3A","CD37","PTPRC")

ngList <- c(nod,inter, immune, astro,other)

fName <- paste0(resDir, resName, ".target_DEG_check.pdf")
pdf(fName,width = 10, height = 6)
DotPlot(mb.comb, features = ngList) + RotatedAxis()
dev.off()

saveRDS(mb.merged, paste0(resDir,resName,".object.rds"))

## trajectory analysis

library(slingshot)

scTarg <- as.SingleCellExperiment(mb.comb, assay="SCT")
scTarg$cluster <- scTarg$seurat_clusters

# no start cluster
scTarg2 <- slingshot(scTarg, clusterLabels = 'cluster', 
                     reducedDim = 'UMAP',approx_points = 100)


require(ggplot2)
require(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

plotcol <- colors[cut(scTarg2$slingPseudotime_1, breaks=100)]
pdf(paste0(resDir, resName, "_slingshot_trajectory.clusters_UMAP.pdf"), width= 8, height=6)

plot(reducedDims(scTarg2)$UMAP, col = plotcol, pch=16, asp = 1,cex=0.5)
lines(SlingshotDataSet(scTarg2), lwd=1, col='black')
#lines(SlingshotDataSet(scTarg2), lwd=1, type = 'lineages', col = 'black')

plot(reducedDims(scTarg2)$UMAP, col = brewer.pal(9,'Set1')[scTarg2$cluster], pch=16, asp = 1,cex=0.5)
lines(SlingshotDataSet(scTarg2), lwd=1, type = 'lineages', col = 'black')

dev.off()


scTarg2 <- slingshot(scTarg, clusterLabels = 'cluster', 
                     start.clus = "0", end.clus = "6",
                     reducedDim = 'UMAP',approx_points = 100)

plotcol <- colors[cut(scTarg2$slingPseudotime_1, breaks=100)]

pdf(paste0(resDir, resName, "_slingshot_trajectory_assigned.UMAP.pdf"), width= 8, height=6)

plot(reducedDims(scTarg2)$UMAP, col = plotcol, pch=16, asp = 1,cex=0.5)
lines(SlingshotDataSet(scTarg2), lwd=1, col='black')
#lines(SlingshotDataSet(scTarg2), lwd=1, type = 'lineages', col = 'black')

plot(reducedDims(scTarg2)$UMAP, col = brewer.pal(9,'Set1')[scTarg2$cluster], pch=16, asp = 1,cex=0.5)
lines(SlingshotDataSet(scTarg2), lwd=1, type = 'lineages', col = 'black')

dev.off()

