# Analysis of single-cell RNA-seq of gluteofemoral adipose

## Load required packages
source("scRNA_seq_functions.R")
library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(scater)
library(scran)
library(ggplot2)
library(reshape2)
library(dplyr)
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

## Filter data
YGSAT_P1[["percent.mt"]] <- PercentageFeatureSet(YGSAT_P1, pattern = "^MT-")
VlnPlot(YGSAT_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(YGSAT_P1$nCount_RNA)
summary(YGSAT_P1$nFeature_RNA)
YGSAT_P1 <- subset(YGSAT_P1, subset = percent.mt < 10 & nCount_RNA < 30000 & nFeature_RNA > 500)

YGSAT_P2[["percent.mt"]] <- PercentageFeatureSet(YGSAT_P2, pattern = "^MT-")
VlnPlot(YGSAT_P2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(YGSAT_P2$nCount_RNA)
summary(YGSAT_P2$nFeature_RNA)
YGSAT_P2 <- subset(YGSAT_P2, subset = percent.mt < 10 & nCount_RNA < 20000 & nFeature_RNA > 500)

YGSAT_P3[["percent.mt"]] <- PercentageFeatureSet(YGSAT_P3, pattern = "^MT-")
VlnPlot(YGSAT_P3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(YGSAT_P3$nCount_RNA)
summary(YGSAT_P3$nFeature_RNA)
YGSAT_P3 <- subset(YGSAT_P3, subset = percent.mt < 10 & nCount_RNA < 30000 & nFeature_RNA > 500)

YGSAT_P4[["percent.mt"]] <- PercentageFeatureSet(YGSAT_P4, pattern = "^MT-")
VlnPlot(YGSAT_P4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(YGSAT_P4$nCount_RNA)
summary(YGSAT_P4$nFeature_RNA)
YGSAT_P4 <- subset(YGSAT_P4, subset = percent.mt < 10 & nCount_RNA < 30000 & nFeature_RNA > 500)

OGSAT_P1[["percent.mt"]] <- PercentageFeatureSet(OGSAT_P1, pattern = "^MT-")
VlnPlot(OGSAT_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(OGSAT_P1$nCount_RNA)
summary(OGSAT_P1$nFeature_RNA)
OGSAT_P1 <- subset(OGSAT_P1, subset = percent.mt < 10 & nCount_RNA < 30000 & nFeature_RNA > 500)

OGSAT_P2[["percent.mt"]] <- PercentageFeatureSet(OGSAT_P2, pattern = "^MT-")
VlnPlot(OGSAT_P2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(OGSAT_P2$nCount_RNA)
summary(OGSAT_P2$nFeature_RNA)
OGSAT_P2 <- subset(OGSAT_P2, subset = percent.mt < 10 & nCount_RNA < 20000 & nFeature_RNA > 500)

OGSAT_P3[["percent.mt"]] <- PercentageFeatureSet(OGSAT_P3, pattern = "^MT-")
VlnPlot(OGSAT_P3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(OGSAT_P3$nCount_RNA)
summary(OGSAT_P3$nFeature_RNA)
OGSAT_P3 <- subset(OGSAT_P3, subset = percent.mt < 10 & nCount_RNA < 20000 & nFeature_RNA > 500)

OGSAT_P4[["percent.mt"]] <- PercentageFeatureSet(OGSAT_P4, pattern = "^MT-")
VlnPlot(OGSAT_P4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(OGSAT_P4$nCount_RNA)
summary(OGSAT_P4$nFeature_RNA)
OGSAT_P4 <- subset(OGSAT_P4, subset = percent.mt < 10 & nCount_RNA < 20000 & nFeature_RNA > 500)

## Merge datasets
lg.adipose<- merge(YGSAT_P1, y = c(YGSAT_P2,YGSAT_P3,YGSAT_P4,OGSAT_P1,OGSAT_P2,OGSAT_P3,OGSAT_P4), 
                   add.cell.ids = c( "YGSAT_P1","YGSAT_P2","YGSAT_P3","YGSAT_P4","OGSAT_P1","OGSAT_P2","OGSAT_P3","OGSAT_P4"), project = "lgAdipose")

## Aalysis of all cells without performing integration
adipose.Unint <- lg.adipose
adipose.Unint <- NormalizeData(adipose.Unint, normalization.method = "LogNormalize", scale.factor = 10000)
adipose.Unint <- FindVariableFeatures(adipose.Unint, selection.method = "vst", nfeatures = 2000)
adipose.Unint <- ScaleData(adipose.Unint, vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
adipose.Unint <- RunPCA(adipose.Unint, features = VariableFeatures(object = adipose.Unint), npcs = 50)
ElbowPlot(adipose.Unint, ndims = 50)
Unint_patients_cols <- c("skyblue4", "skyblue3", "skyblue2", "skyblue1", "sienna4", "sienna3", "sienna2", "sienna1")
DimPlot(adipose.Unint, group.by = "orig.ident", reduction = "pca", cols = Unint_patients_cols, combine = FALSE, pt.size = 0.1)
adipose.Unint <- RunTSNE(adipose.Unint, dims = 1:20)
DimPlot(adipose.Unint, group.by = "orig.ident", reduction = "tsne", cols = Unint_patients_cols, combine = FALSE, pt.size = 0.1)

# Perform integrative analysis
adipose.list <- SplitObject(lg.adipose, split.by = "orig.ident")
for (i in names(adipose.list)) {
  adipose.list[[i]] <- SCTransform(adipose.list[[i]], vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"),verbose = TRUE)
}

adipose.features <- SelectIntegrationFeatures(object.list = adipose.list, nfeatures = 1000)
adipose.list <- PrepSCTIntegration(object.list = adipose.list, anchor.features = adipose.features, verbose = TRUE)
adipose.anchors <- FindIntegrationAnchors(object.list = adipose.list, normalization.method = "SCT", 
                                          anchor.features = adipose.features, verbose = TRUE)
adipose.integrated <- IntegrateData(anchorset = adipose.anchors, normalization.method = "SCT", verbose = TRUE)
adipose.integrated <- RunPCA(adipose.integrated, verbose = TRUE)
adipose.integrated <- RunTSNE(adipose.integrated, dims = 1:10)
adipose.integrated <- FindNeighbors(adipose.integrated, reduction = "pca", dims = 1:10)
adipose.integrated <- FindClusters(adipose.integrated, resolution = 0.5)

new.cluster.ids <- c("APC1", "APC2", "APC3", "APC4", "APC5", "APC6","VEC","IC1", "IC2", "IC3", "IC4", "SMC")
names(new.cluster.ids) <- levels(adipose.integrated)
adipose.integrated <- RenameIdents(adipose.integrated, new.cluster.ids)
adipose.integrated$clusters <- Idents(adipose.integrated)

adipose.integrated_colorPalette <- c(RColorBrewer::brewer.pal(12,"Paired"))
DimPlot(adipose.integrated, reduction = "tsne", group.by = "clusters", split.by = "orig.ident", cols = adipose.integrated_colorPalette, pt.size = 0.5) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)))

## Violin plot of key markers
DefaultAssay(adipose.integrated) <- "RNA"
adipose.integrated <- NormalizeData(adipose.integrated, verbose = TRUE)
VlnPlot(adipose.integrated, features = c("PDGFRA", "PECAM1", "PTPRC","ACTA2"), cols = adipose.integrated_colorPalette, fill.by = "ident", stack = TRUE, flip = TRUE)







