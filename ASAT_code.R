## Analysis of single-cell RNA-seq of abdominal adipose

## Load required packages
source("scRNA_seq_functions.R")
library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(SingleCellExperiment)
library(scran)
library(SingleR)
library(scater)
library(Matrix)
library(pheatmap)
library(sfsmisc)
library(MASS)

## Load data
load("YASAT_P1.Rdata")
load("YASAT_P2.Rdata")
load("YASAT_P3.Rdata")
load("OASAT_P1.Rdata")
load("OASAT_P2.Rdata")
load("OASAT_P3.Rdata")

## Fliter data
YASAT_P1[["percent.mt"]] <- PercentageFeatureSet(YASAT_P1, pattern = "^MT-")
VlnPlot(YASAT_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(YASAT_P1$nCount_RNA)
summary(YASAT_P1$nFeature_RNA)
YASAT_P1 <- subset(YASAT_P1, subset = percent.mt < 10 & nCount_RNA < 30000 & nFeature_RNA > 500)

YASAT_P2[["percent.mt"]] <- PercentageFeatureSet(YASAT_P2, pattern = "^MT-")
VlnPlot(YASAT_P2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(YASAT_P2$nCount_RNA)
summary(YASAT_P2$nFeature_RNA)
YASAT_P2 <- subset(YASAT_P2, subset = percent.mt < 10 & nCount_RNA < 20000 & nFeature_RNA > 500)

YASAT_P3[["percent.mt"]] <- PercentageFeatureSet(YASAT_P3, pattern = "^MT-")
VlnPlot(YASAT_P3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(YASAT_P3$nCount_RNA)
summary(YASAT_P3$nFeature_RNA)
YASAT_P3 <- subset(YASAT_P3, subset = percent.mt < 10 & nCount_RNA < 30000 & nFeature_RNA > 500)

OASAT_P1[["percent.mt"]] <- PercentageFeatureSet(OASAT_P1, pattern = "^MT-")
VlnPlot(OASAT_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(OASAT_P1$nCount_RNA)
summary(OASAT_P1$nFeature_RNA)
OASAT_P1 <- subset(OASAT_P1, subset = percent.mt < 10 & nCount_RNA < 30000 & nFeature_RNA > 500)

OASAT_P2[["percent.mt"]] <- PercentageFeatureSet(OASAT_P2, pattern = "^MT-")
VlnPlot(OASAT_P2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(OASAT_P2$nCount_RNA)
summary(OASAT_P2$nFeature_RNA)
OASAT_P2 <- subset(OASAT_P2, subset = percent.mt < 10 & nCount_RNA < 20000 & nFeature_RNA > 500)

OASAT_P3[["percent.mt"]] <- PercentageFeatureSet(OASAT_P3, pattern = "^MT-")
VlnPlot(OASAT_P3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(OASAT_P3$nCount_RNA)
summary(OASAT_P3$nFeature_RNA)
OASAT_P3 <- subset(OASAT_P3, subset = percent.mt < 10 & nCount_RNA < 20000 & nFeature_RNA > 500)

## Merge datasets
ab.adipose <- merge(YASAT_P1, y = c(YASAT_P2,YASAT_P3,OASAT_P1,OASAT_P2,OASAT_P3), add.cell.ids = c( "YASAT_P1", "YASAT_P2","YASAT_P3", "OASAT_P1","OASAT_P2", "OASAT_P3"), project = "AbdomenAdipose")

## Aalysis of all cells without performing integration
adipose.Unint <- ab.adipose
adipose.Unint <- NormalizeData(adipose.Unint, normalization.method = "LogNormalize", scale.factor = 10000)
adipose.Unint <- FindVariableFeatures(adipose.Unint, selection.method = "vst", nfeatures = 2000)
adipose.Unint <- ScaleData(adipose.Unint, vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
adipose.Unint <- RunPCA(adipose.Unint, features = VariableFeatures(object = adipose.Unint), npcs = 50)
adipose.Unint <- RunTSNE(adipose.Unint, dims = 1:20)
Unint_patients_cols <- c("skyblue4", "skyblue3", "skyblue2", "sienna4", "sienna3", "sienna2")
DimPlot(adipose.Unint, group.by = "orig.ident", reduction = "pca", cols = Unint_patients_cols, combine = FALSE, pt.size = 0.1)
DimPlot(adipose.Unint, group.by = "orig.ident", reduction = "tsne", cols = Unint_patients_cols, combine = FALSE, pt.size = 0.1)

## Perform integrative analysis
adipose.list <- SplitObject(ab.adipose, split.by = "orig.ident")

for (i in names(adipose.list)) {
  adipose.list[[i]] <- SCTransform(adipose.list[[i]], vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"),verbose = TRUE)
}

adipose.features <- SelectIntegrationFeatures(object.list = adipose.list, nfeatures = 1000)
adipose.list <- PrepSCTIntegration(object.list = adipose.list, anchor.features = adipose.features, verbose = TRUE)
adipose.anchors <- FindIntegrationAnchors(object.list = adipose.list, normalization.method = "SCT", 
                                          anchor.features = adipose.features, verbose = TRUE)
adipose.integrated <- IntegrateData(anchorset = adipose.anchors, normalization.method = "SCT", verbose = TRUE)
adipose.integrated <- RunPCA(adipose.integrated, verbose = TRUE)
### adipose.integrated <- JackStraw(adipose.integrated, num.replicate = 200, dims = 50)
### adipose.integrated <- ScoreJackStraw(adipose.integrated, dims = 1:50)
### ElbowPlot(adipose.integrated,ndims =30)
### JackStrawPlot(adipose.integrated, dims = 1:50)
adipose.integrated <- RunTSNE(adipose.integrated, dims = 1:9)
adipose.integrated <- FindNeighbors(adipose.integrated, reduction = "pca", dims = 1:9)
adipose.integrated <- FindClusters(adipose.integrated, resolution = 0.5)

new.cluster.ids <- c("APC1", "APC2", "APC3", "APC4", "APC5", "APC6","IC1", "IC2", "APC7", "IC3", "SMC", "VEC")
names(new.cluster.ids) <- levels(adipose.integrated)
adipose.integrated <- RenameIdents(adipose.integrated, new.cluster.ids)
adipose.integrated$clusters <- Idents(adipose.integrated)

adipose.integrated_colorPalette <- c(RColorBrewer::brewer.pal(12,"Paired"))
DimPlot(adipose.integrated, reduction = "tsne", group.by = "clusters", cols = adipose.integrated_colorPalette, pt.size = 0.5) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3)))
DimPlot(adipose.integrated, reduction = "tsne", group.by = "clusters", split.by = "Age", cols = adipose.integrated_colorPalette, pt.size = 0.5) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3)))
DimPlot(adipose.integrated, reduction = "tsne", group.by = "Age", split.by = "Age", pt.size = 0.5) + 
  theme(legend.position = "top") 
DimPlot(adipose.integrated, reduction = "tsne", group.by = "clusters", split.by = "orig.ident", cols = adipose.integrated_colorPalette, pt.size = 0.5) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)))
DimPlot(adipose.integrated, group.by = "orig.ident", reduction = "tsne", cols = Unint_patients_cols, combine = FALSE, pt.size = 0.1)

## Analysis of APC sub-populations
DefaultAssay(adipose.integrated) <- "RNA"
adipose.integrated <- NormalizeData(adipose.integrated, verbose = TRUE)
adipose.integrated <- ScaleData(adipose.integrated)
apc.integrated <- subset(adipose.integrated, idents = paste0("APC", 1:7))
apc.integrated_colorPalette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#CAB2D6")
DimPlot(apc.integrated, reduction = "tsne", group.by = "clusters", cols = apc.integrated_colorPalette, pt.size = 0.5) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)))

## Find common and unique DEGs between old and young APC (Old vs Young)
apc.integrated$group <- paste(apc.integrated$clusters, apc.integrated$Age, sep = "_")
OA_vs_YA_APC_DEGs <- getDEGs(seurat.object = apc.integrated, ident.1 = "Old", ident.2 = "Young")
OA_vs_YA_APC_common_up_DEGs <- getCommonUpDEGs(OA_vs_YA_APC_DEGs)
OA_vs_YA_APC_common_down_DEGs <- getCommonDownDEGs(OA_vs_YA_APC_DEGs)
OA_vs_YA_APC_unique_up_DEGs <- getUniqueUpGenes(OA_vs_YA_APC_DEGs)
OA_vs_YA_APC_unique_down_DEGs <- getUniqueDownGenes(OA_vs_YA_APC_DEGs)

## plot OA_vs_YA_APC_common_and_unique_up_genes
OA_vs_YA_APC_up_DEGs <- list()
for (i in 1:length(OA_vs_YA_APC_DEGs)) {
  DEG_data <- as.data.frame(OA_vs_YA_APC_DEGs[[i]])
  OA_vs_YA_APC_up_DEGs[[i]] <- rownames(DEG_data)[which(DEG_data$avg_logFC>0)]
}

num_dec1 <- matrix(0,ncol = 7, nrow = 355)
for (i in 1:355) {
  for (j in 1:7){
    if (OA_vs_YA_APC_common_up_DEGs[i] %in% OA_vs_YA_APC_up_DEGs[[j]]) {
      num_dec1[i,j] <- 1
    }
  }
}

OA_vs_YA_APC_unique_up_DEGs_list <- unlist(OA_vs_YA_APC_unique_up_DEGs)
num_dec2 <- matrix(0, ncol = 7, nrow = 270)
for (i in 1:270) {
  for (j in 1:7){
    if (OA_vs_YA_APC_unique_up_DEGs_list[i] %in% OA_vs_YA_APC_up_DEGs[[j]]) {
      num_dec2[i,j] <- 1
    }
  }
}

num_dec <- rbind(num_dec1,num_dec2)
colnames(num_dec) <- levels(Idents(apc.integrated))
pheatmap(num_dec, cluster_rows = F, cluster_cols = F, angle_col = 0, color = colorRampPalette(c("lightgrey", "white", "#F8766C"))(100))

### plot OA_vs_YA_APC_common_and_unique_down_genes
OA_vs_YA_APC_down_DEGs <- list()
for (i in 1:length(OA_vs_YA_APC_DEGs)) {
  DEG_data <- as.data.frame(OA_vs_YA_APC_DEGs[[i]])
  OA_vs_YA_APC_down_DEGs[[i]] <- rownames(DEG_data)[which(DEG_data$avg_logFC<0)]
}

num_dec1 <- matrix(0,ncol = 7, nrow = 406)
for (i in 1:406) {
  for (j in 1:7){
    if (OA_vs_YA_APC_common_down_DEGs[i] %in% OA_vs_YA_APC_down_DEGs[[j]]) {
      num_dec1[i,j] <- 1
    }
  }
}

OA_vs_YA_APC_unique_down_DEGs_list <- unlist(OA_vs_YA_APC_unique_down_DEGs)
num_dec2 <- matrix(0,ncol = 7, nrow = 257)
for (i in 1:257) {
  for (j in 1:7){
    if (OA_vs_YA_APC_unique_down_DEGs_list[i] %in% OA_vs_YA_APC_down_DEGs[[j]]) {
      num_dec2[i,j] <- 1
    }
  }
}

num_dec <- rbind(num_dec1,num_dec2)
colnames(num_dec) <- levels(Idents(apc.integrated))
pheatmap(num_dec, cluster_rows = F, cluster_cols = F, angle_col = 0, color = colorRampPalette(c("lightgrey", "white", "#F8766C"))(100))


## Analysis of immune cells subpopulation
immune.cells.idents <- names(adipose.integrated$clusters)[c(which(adipose.integrated$clusters == "IC1"),
                                                            which(adipose.integrated$clusters == "IC2"),
                                                            which(adipose.integrated$clusters == "IC3"))]
immune.cells <- ab.adipose[,immune.cells.idents]
immune.cells.list <- SplitObject(immune.cells, split.by = "orig.ident")

for (i in names(immune.cells.list)) {
  immune.cells.list[[i]] <- SCTransform(immune.cells.list[[i]], vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"),verbose = TRUE)
}
immune.cells.features <- SelectIntegrationFeatures(object.list = immune.cells.list, nfeatures = 1000)
immune.cells.list <- PrepSCTIntegration(object.list = immune.cells.list, anchor.features = immune.cells.features, verbose = TRUE)
immune.cells.anchors <- FindIntegrationAnchors(object.list = immune.cells.list, normalization.method = "SCT", 
                                               anchor.features = immune.cells.features, k.filter = 110, verbose = TRUE)
immune.cells.integrated <- IntegrateData(anchorset = immune.cells.anchors, normalization.method = "SCT", verbose = TRUE)
immune.cells.integrated <- RunPCA(immune.cells.integrated, verbose = FALSE)
### ElbowPlot(immune.cells.integrated,ndims = 50)
### immune.cells.integrated <- JackStraw(immune.cells.integrated, num.replicate = 200, dims = 50)
### immune.cells.integrated <- ScoreJackStraw(immune.cells.integrated, dims = 1:50)
### JackStrawPlot(immune.cells.integrated, dims = 1:50)
immune.cells.integrated <- RunTSNE(immune.cells.integrated, dims = 1:10)
immune.cells.integrated <- RunUMAP(immune.cells.integrated, dims = 1:10)
immune.cells.integrated <- FindNeighbors(immune.cells.integrated, reduction = "pca", dims = 1:10)
immune.cells.integrated <- FindClusters(immune.cells.integrated, resolution = 0.4)

new.cluster.ids.ic <- c("ICS1","ICS2","ICS3","ICS4","ICS5","ICS6","ICS7","ICS8")
names(new.cluster.ids.ic) <- levels(immune.cells.integrated)
immune.cells.integrated <- RenameIdents(immune.cells.integrated, new.cluster.ids.ic)
immune.cells.integrated$clusters <- Idents(immune.cells.integrated)

immune.cells_colorPalette <- c(RColorBrewer::brewer.pal(8,"Set2"))
DimPlot(immune.cells.integrated, reduction = "tsne", group.by = "orig.ident", cols = Unint_patients_cols, pt.size = 1) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3)))
DimPlot(immune.cells.integrated, reduction = "tsne", group.by = "clusters",cols = immune.cells_colorPalette, pt.size = 1) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)))
DimPlot(immune.cells.integrated, reduction = "tsne", group.by = "clusters", split.by = "Age", cols = immune.cells_colorPalette, pt.size = 1) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)))
DimPlot(immune.cells.integrated, reduction = "tsne", group.by = "Age", split.by = "orig.ident", cols = c("skyblue3","sienna3"), pt.size = 0.5) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)))

## Find common and unique DEGs between IC clusters
```{r}
DefaultAssay(immune.cells.integrated) <- "RNA"
immune.cells.integrated <- NormalizeData(immune.cells.integrated, verbose = TRUE)
immune.cells.integrated <- ScaleData(immune.cells.integrated)
immune.cells.integrated$group <- paste(immune.cells.integrated$clusters, immune.cells.integrated$Age, sep = "_")
OA_vs_YA_IC_DEGs <- getDEGs(seurat.object = immune.cells.integrated, ident.1 = "Old", ident.2 = "Young")
OA_vs_YA_IC_common_up_DEGs <- getCommonUpDEGs(OA_vs_YA_IC_DEGs)
OA_vs_YA_IC_common_down_DEGs <- getCommonDownDEGs(OA_vs_YA_IC_DEGs)
OA_vs_YA_IC_unique_up_DEGs <- getUniqueUpGenes(OA_vs_YA_IC_DEGs)
OA_vs_YA_IC_unique_down_DEGs <- getUniqueDownGenes(OA_vs_YA_IC_DEGs)

## plot OA_vs_YA_IC_common_and_unique_up_genes
OA_vs_YA_IC_up_DEGs <- list()
for (i in 1:length(OA_vs_YA_IC_DEGs)) {
  DEG_data <- as.data.frame(OA_vs_YA_IC_DEGs[[i]])
  OA_vs_YA_IC_up_DEGs[[i]] <- rownames(DEG_data)[which(DEG_data$avg_logFC>0)]
}

num_dec1 <- matrix(0,ncol = 8, nrow = 23)
for (i in 1:23) {
  for (j in 1:8){
    if (OA_vs_YA_IC_common_up_DEGs[i] %in% OA_vs_YA_IC_up_DEGs[[j]]) {
      num_dec1[i,j] <- 1
    }
  }
}

OA_vs_YA_IC_unique_up_DEGs_list <- unlist(OA_vs_YA_IC_unique_up_DEGs)
num_dec2 <- matrix(0,ncol = 8, nrow = 112)
for (i in 1:112) {
  for (j in 1:8){
    if (OA_vs_YA_IC_unique_up_DEGs_list[i] %in% OA_vs_YA_IC_up_DEGs[[j]]) {
      num_dec2[i,j] <- 1
    }
  }
}

num_dec <- rbind(num_dec1,num_dec2)
colnames(num_dec) <- levels(Idents(immune.cells.integrated))
pheatmap(num_dec, cluster_rows = F, cluster_cols = F, angle_col = 0, color = colorRampPalette(c("lightgrey", "white", "#F8766C"))(100))

## plot OA_vs_YA_IC_common_and_unique_down_genes
OA_vs_YA_IC_down_DEGs <- list()
for (i in 1:length(OA_vs_YA_IC_DEGs)) {
  DEG_data <- as.data.frame(OA_vs_YA_IC_DEGs[[i]])
  OA_vs_YA_IC_down_DEGs[[i]] <- rownames(DEG_data)[which(DEG_data$avg_logFC<0)]
}

num_dec1 <- matrix(0,ncol = 8, nrow = 95)
for (i in 1:95) {
  for (j in 1:8){
    if (OA_vs_YA_IC_common_down_DEGs[i] %in% OA_vs_YA_IC_down_DEGs[[j]]) {
      num_dec1[i,j] <- 1
    }
  }
}

OA_vs_YA_IC_unique_down_DEGs_list <- unlist(OA_vs_YA_IC_unique_down_DEGs)
num_dec2 <- matrix(0,ncol = 8, nrow = 471)
for (i in 1:471) {
  for (j in 1:8){
    if (OA_vs_YA_IC_unique_down_DEGs_list[i] %in% OA_vs_YA_IC_down_DEGs[[j]]) {
      num_dec2[i,j] <- 1
    }
  }
}

num_dec <- rbind(num_dec1,num_dec2)
colnames(num_dec) <- levels(Idents(immune.cells.integrated))
pheatmap(num_dec, cluster_rows = F, cluster_cols = F, angle_col = 0, color = colorRampPalette(c("lightgrey", "white", "#F8766C"))(100))
