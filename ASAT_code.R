source("scRNA_seq_functions.R")
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(scran)
library(SingleR)
library(scater)
library(Matrix)
library(pheatmap)
library(sfsmisc)
library(MASS)

###########################################################################################################################
# Pre-processing
###########################################################################################################################
# load count data
YLSAT <- Read10X(data.dir = "./data/Young_Lg_SAT/")
YASAT <- Read10X(data.dir = "./data/Young_Ab_SAT/")
OLSAT <- Read10X(data.dir = "./data/Old_Lg_SAT/")
OASAT <- Read10X(data.dir = "./data/Old_Ab_SAT/")

# load souporcell data
YLSAT_souporcell <- read.table("./souporcell/YLSAT/clusters.tsv", header = T)
YASAT_souporcell <- read.table("./souporcell/YASAT/clusters.tsv", header = T)
OLSAT_souporcell <- read.table("./souporcell/OLSAT/clusters.tsv", header = T)
OASAT_souporcell <- read.table("./souporcell/OASAT/clusters.tsv", header = T)

# split data
YLSAT_P1 <- YLSAT[,which(YLSAT_souporcell$assignment == "0")]
YLSAT_P2 <- YLSAT[,which(YLSAT_souporcell$assignment == "1")]
YLSAT_P3 <- YLSAT[,which(YLSAT_souporcell$assignment == "2")]
YLSAT_P4 <- YLSAT[,which(YLSAT_souporcell$assignment == "3")]
OLSAT_P1 <- OLSAT[,which(OLSAT_souporcell$assignment == "0")]
OLSAT_P2 <- OLSAT[,which(OLSAT_souporcell$assignment == "1")]
OLSAT_P3 <- OLSAT[,which(OLSAT_souporcell$assignment == "2")]
OLSAT_P4 <- OLSAT[,which(OLSAT_souporcell$assignment == "3")]
YASAT_P1 <- YASAT[,which(YASAT_souporcell$assignment == "0")]
YASAT_P2 <- YASAT[,which(YASAT_souporcell$assignment == "1")]
YASAT_P3 <- YASAT[,which(YASAT_souporcell$assignment == "2")]
OASAT_P1 <- OASAT[,which(OASAT_souporcell$assignment == "0")]
OASAT_P2 <- OASAT[,which(OASAT_souporcell$assignment == "1")]
OASAT_P3 <- OASAT[,which(OASAT_souporcell$assignment == "2")]

# Create Seurat Object and subset data
YASAT_P1S <- CreateSeuratObject(counts = YASAT_P1, project = "YASAT_P1", min.cells = 10, min.features = 0)
YASAT_P1S$Age <- "Young"
YASAT_P1S$Anatomic_site <- "Abdomen"
YASAT_P1S[["percent.mt"]] <- PercentageFeatureSet(YASAT_P1S, pattern = "^MT-")
VlnPlot(YASAT_P1S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(YASAT_P1S$nCount_RNA)
summary(YASAT_P1S$nFeature_RNA)
YASAT_P1S <- subset(YASAT_P1S, subset = percent.mt < 10 & nCount_RNA < 30000 & nFeature_RNA > 500)

YASAT_P2S <- CreateSeuratObject(counts = YASAT_P2, project = "YASAT_P2", min.cells = 10, min.features = 0)
YASAT_P2S$Age <- "Young"
YASAT_P2S$Anatomic_site <- "Abdomen"
YASAT_P2S[["percent.mt"]] <- PercentageFeatureSet(YASAT_P2S, pattern = "^MT-")
VlnPlot(YASAT_P2S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(YASAT_P2S$nCount_RNA)
summary(YASAT_P2S$nFeature_RNA)
YASAT_P2S <- subset(YASAT_P2S, subset = percent.mt < 10 & nCount_RNA < 20000 & nFeature_RNA > 500)

YASAT_P3S <- CreateSeuratObject(counts = YASAT_P3, project = "YASAT_P3", min.cells = 10, min.features = 0)
YASAT_P3S$Age <- "Young"
YASAT_P3S$Anatomic_site <- "Abdomen"
YASAT_P3S[["percent.mt"]] <- PercentageFeatureSet(YASAT_P3S, pattern = "^MT-")
VlnPlot(YASAT_P3S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(YASAT_P3S$nCount_RNA)
summary(YASAT_P3S$nFeature_RNA)
YASAT_P3S <- subset(YASAT_P3S, subset = percent.mt < 10 & nCount_RNA < 30000 & nFeature_RNA > 500)

OASAT_P1S <- CreateSeuratObject(counts = OASAT_P1, project = "OASAT_P1", min.cells = 10, min.features = 0)
OASAT_P1S$Age <- "Old"
OASAT_P1S$Anatomic_site <- "Abdomen"
OASAT_P1S[["percent.mt"]] <- PercentageFeatureSet(OASAT_P1S, pattern = "^MT-")
VlnPlot(OASAT_P1S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(OASAT_P1S$nCount_RNA)
summary(OASAT_P1S$nFeature_RNA)
OASAT_P1S <- subset(OASAT_P1S, subset = percent.mt < 10 & nCount_RNA < 30000 & nFeature_RNA > 500)

OASAT_P2S <- CreateSeuratObject(counts = OASAT_P2, project = "OASAT_P2", min.cells = 10, min.features = 0)
OASAT_P2S$Age <- "Old"
OASAT_P2S$Anatomic_site <- "Abdomen"
OASAT_P2S[["percent.mt"]] <- PercentageFeatureSet(OASAT_P2S, pattern = "^MT-")
VlnPlot(OASAT_P2S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(OASAT_P2S$nCount_RNA)
summary(OASAT_P2S$nFeature_RNA)
OASAT_P2S <- subset(OASAT_P2S, subset = percent.mt < 10 & nCount_RNA < 20000 & nFeature_RNA > 500)

OASAT_P3S <- CreateSeuratObject(counts = OASAT_P3, project = "OASAT_P3", min.cells = 10, min.features = 0)
OASAT_P3S$Age <- "Old"
OASAT_P3S$Anatomic_site <- "Abdomen"
OASAT_P3S[["percent.mt"]] <- PercentageFeatureSet(OASAT_P3S, pattern = "^MT-")
VlnPlot(OASAT_P3S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
summary(OASAT_P3S$nCount_RNA)
summary(OASAT_P3S$nFeature_RNA)
OASAT_P3S <- subset(OASAT_P3S, subset = percent.mt < 10 & nCount_RNA < 20000 & nFeature_RNA > 500)

ab.adipose <- merge(YASAT_P1S, y = c(YASAT_P2S,YASAT_P3S,OASAT_P1S,OASAT_P2S,OASAT_P3S), 
                    add.cell.ids = c( "YASAT_P1", "YASAT_P2","YASAT_P3", "OASAT_P1","OASAT_P2", "OASAT_P3"), project = "AbdomenAdipose")
###########################################################################################################################



###########################################################################################################################
# Data analysis of all cells without performing integration
###########################################################################################################################
adipose.Unint <- ab.adipose
adipose.Unint <- NormalizeData(adipose.Unint, normalization.method = "LogNormalize", scale.factor = 10000)
adipose.Unint <- FindVariableFeatures(adipose.Unint, selection.method = "vst", nfeatures = 2000)
adipose.Unint <- ScaleData(adipose.Unint, vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
adipose.Unint <- RunPCA(adipose.Unint, features = VariableFeatures(object = adipose.Unint), npcs = 50)
adipose.Unint <- RunTSNE(adipose.Unint, dims = 1:20)
# Define coloring
Unint_patients_cols <- c("skyblue4", "skyblue3", "skyblue2", "sienna4", "sienna3", "sienna2")
DimPlot(adipose.Unint, group.by = "orig.ident", reduction = "pca", cols = Unint_patients_cols, combine = FALSE, pt.size = 0.1)
DimPlot(adipose.Unint, group.by = "orig.ident", reduction = "tsne", cols = Unint_patients_cols, combine = FALSE, pt.size = 0.1)
###########################################################################################################################



###########################################################################################################################
# Perform integrative analysis
###########################################################################################################################
# setup the Seurat object list
adipose.list <- SplitObject(ab.adipose, split.by = "orig.ident")

# Perform SCTransform normalization separately for each dataset
for (i in names(adipose.list)) {
  adipose.list[[i]] <- SCTransform(adipose.list[[i]], vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"),verbose = TRUE)
}

# Select features for downstream integration 
adipose.features <- SelectIntegrationFeatures(object.list = adipose.list, nfeatures = 1000)
adipose.list <- PrepSCTIntegration(object.list = adipose.list, anchor.features = adipose.features, verbose = TRUE)

# Identify anchors and integrate the datasets
adipose.anchors <- FindIntegrationAnchors(object.list = adipose.list, normalization.method = "SCT", 
                                          anchor.features = adipose.features, verbose = TRUE)
adipose.integrated <- IntegrateData(anchorset = adipose.anchors, normalization.method = "SCT", verbose = TRUE)

#  Dim reduction and cluster Visualization
adipose.integrated <- RunPCA(adipose.integrated, verbose = TRUE)
DimPlot(adipose.integrated, group.by = "orig.ident", reduction = "pca",combine = FALSE, pt.size = 0.1) 
ggsave("Results_AllCellsIntegrated/PCA_integratd.pdf", height = 5, width = 6)

# Determine the number of dims for clustering
adipose.integrated <- JackStraw(adipose.integrated, num.replicate = 200, dims = 50)
adipose.integrated <- ScoreJackStraw(adipose.integrated, dims = 1:50)
ElbowPlot(adipose.integrated,ndims =30)
JackStrawPlot(adipose.integrated, dims = 1:50)

# Clustering and visualization
adipose.integrated <- RunTSNE(adipose.integrated, dims = 1:9)
adipose.integrated <- RunUMAP(adipose.integrated, dims = 1:9)
adipose.integrated <- FindNeighbors(adipose.integrated, reduction = "pca", dims = 1:9)
adipose.integrated <- FindClusters(adipose.integrated, resolution = 0.5)

# Assign cell types
new.cluster.ids <- c("APC1", "APC2", "APC3", "APC4", "APC5", "APC6","IC1", "IC2", "APC7", "IC3", "SMC", "VEC")
names(new.cluster.ids) <- levels(adipose.integrated)
adipose.integrated <- RenameIdents(adipose.integrated, new.cluster.ids)
adipose.integrated$clusters <- Idents(adipose.integrated)

# TSNE plots
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

# Plot representative cell markers
FeaturePlot(adipose.integrated, features = c("PDGFRA","PECAM1","PTPRC","ACTA2"),
            reduction = "tsne", pt.size = 0.1, order = TRUE, label.size = 2)

VlnPlot(adipose.integrated, features = c("PDGFRA", "PECAM1", "PTPRC", "ACTA2"), cols = adipose.integrated_colorPalette, 
        fill.by = "ident", stack = T, flip = TRUE) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45)) +
  NoLegend()
###########################################################################################################################



###########################################################################################################################
# Analysis of APC sub-populations
###########################################################################################################################
# subset integrated APC data
apc.integrated <- subset(adipose.integrated, idents = paste0("APC", 1:7))

# TSNE plots
apc.integrated_colorPalette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#CAB2D6")

DimPlot(apc.integrated, reduction = "tsne", group.by = "clusters", cols = apc.integrated_colorPalette, pt.size = 0.5) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)))

# Plots of known markers
DefaultAssay(apc.integrated) <- "RNA"
FeaturePlot(apc.integrated, features = c("CD55", "PI16", "SEMA3C", "DPP4"), reduction = "tsne", pt.size = 0.05, label.size = 2, ncol = 4)

FeaturePlot(apc.integrated, features = c("FABP4", "APOE", "GGT5", "DEPP1"), reduction = "tsne",  pt.size = 0.05, label.size = 2, ncol = 4)

DotPlot(apc.integrated, features = c("CYTOR","S100A16","HSPA6","HSPA1B","SFRP4","FBN1","ZFP36","JUNB","CXCL8","CXCL3","TXNIP","MGP","CCL2","IGFBP7"),  cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() + labs(x = "", y = "")

# Find DEGs between old and young APC (Old vs Young)
apc.integrated$group <- paste(apc.integrated$clusters, apc.integrated$Age, sep = "_")
OA_vs_YA_APC_DEGs <- getDEGs(seurat.object = apc.integrated, ident.1 = "Old", ident.2 = "Young")

# Calculate proportion of APCs in each cluster
nApcsInCluster <- table(Idents(apc.integrated), apc.integrated$orig.ident)
nApcsInCluster <- melt(nApcsInCluster,varnames = c("Cluster","orig.ident"))
propApcsInCluster <- as.data.frame(prop.table(table(Idents(apc.integrated), apc.integrated$orig.ident), margin = 2))
propApcsInCluster$Age <- as.factor(c(rep("Old",21),rep("Young",21)))
propApcsInCluster <- subset(propApcsInCluster,select = -Var2)
df_apc <- data_summary(propApcsInCluster, varname = "Freq", groupnames = c("Age", "Var1"))
Varx <- factor(c("APC1","APC2","APC3","APC4","APC5","APC6","APC7","APC1","APC2","APC3","APC4","APC5","APC6","APC7"), levels = c("APC1","APC2","APC3","APC4","APC5","APC6","APC7") )
df_apc$Varx <- Varx
ggplot(df_apc, aes(x = Varx, y = Freq, fill = Age))+
  geom_bar(stat = "identity", width = .5, position = position_dodge2(reverse = TRUE))+
  geom_errorbar(aes(ymin= Freq-sd, ymax = Freq+sd), width = .2, position = position_dodge2(reverse = TRUE))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  theme_classic() +
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black")) +
  theme( axis.text.x = element_text(size = 12,  color = "black", angle = 45, vjust = 0.5)) +
  labs(x="",y="Cell type distribution",title="")

# Heatmap of marker genes of APC
apc_markers <- FindAllMarkers(apc.integrated, min.pct = 0.25, only.pos = TRUE, verbose = TRUE)
apc_markers <- apc_markers[which(apc_markers$p_val_adj < 0.05),]
apc_markers_top50 <-  apc_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(apc.integrated, features = apc_markers_top50$gene) + 
  theme(text = element_text(size = 0))
###########################################################################################################################



###########################################################################################################################
# Analysis of immune cells subpopulation
###########################################################################################################################
# subset data
immune.cells.idents <- names(adipose.integrated$clusters)[c(which(adipose.integrated$clusters == "IC1"),
                                                            which(adipose.integrated$clusters == "IC2"),
                                                            which(adipose.integrated$clusters == "IC3"))]
immune.cells <- ab.adipose[,immune.cells.idents]
immune.cells.list <- SplitObject(immune.cells, split.by = "orig.ident")

# Perform SCTransform normalization separately for each dataset
for (i in names(immune.cells.list)) {
  immune.cells.list[[i]] <- SCTransform(immune.cells.list[[i]], vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt"),verbose = TRUE)
}

# Select features for downstream integration 
immune.cells.features <- SelectIntegrationFeatures(object.list = immune.cells.list, nfeatures = 1000)
immune.cells.list <- PrepSCTIntegration(object.list = immune.cells.list, anchor.features = immune.cells.features, verbose = TRUE)

# Identify anchors and integrate the datasets
immune.cells.anchors <- FindIntegrationAnchors(object.list = immune.cells.list, normalization.method = "SCT", 
                                               anchor.features = immune.cells.features, k.filter = 110, verbose = TRUE)
immune.cells.integrated <- IntegrateData(anchorset = immune.cells.anchors, normalization.method = "SCT", verbose = TRUE)

# Dim reduction and Visualization
immune.cells.integrated <- RunPCA(immune.cells.integrated, verbose = FALSE)

# Determine the number of dims for clustering
ElbowPlot(immune.cells.integrated,ndims = 50)
immune.cells.integrated <- JackStraw(immune.cells.integrated, num.replicate = 200, dims = 50)
immune.cells.integrated <- ScoreJackStraw(immune.cells.integrated, dims = 1:50)
JackStrawPlot(immune.cells.integrated, dims = 1:50)

# Clustering and visualization
DefaultAssay(immune.cells.integrated) <- "integrated"
table(immune.cells.integrated$seurat_clusters)
immune.cells_colorPalette <- c(RColorBrewer::brewer.pal(6,"Set2"))
immune.cells.integrated <- RunTSNE(immune.cells.integrated, dims = 1:10)
immune.cells.integrated <- RunUMAP(immune.cells.integrated, dims = 1:10)
immune.cells.integrated <- FindNeighbors(immune.cells.integrated, reduction = "pca", dims = 1:10)
immune.cells.integrated <- FindClusters(immune.cells.integrated, resolution = 0.4)

# Assign cell types
new.cluster.ids.ic <- c("ICS1","ICS2","ICS3","ICS4","ICS5","ICS6","ICS7","ICS8")
names(new.cluster.ids.ic) <- levels(immune.cells.integrated)
immune.cells.integrated <- RenameIdents(immune.cells.integrated, new.cluster.ids.ic)
immune.cells.integrated$clusters <- Idents(immune.cells.integrated)

# proportion of IC clusters
nIcsInCluster <- table(Idents(immune.cells.integrated), immune.cells.integrated$orig.ident)
nIcsInCluster <- melt(nIcsInCluster,varnames = c("Cluster","orig.ident"))
propIcsInCluster <- as.data.frame(prop.table(table(Idents(immune.cells.integrated), immune.cells.integrated$orig.ident), margin = 2))
propIcsInCluster$Age <- as.factor(c(rep("Old",24),rep("Young",24)))
propIcsInCluster <- subset(propIcsInCluster,select = -Var2)
write.csv(propIcsInCluster, file = "Results_ImmuneCellsIntegrated/Immune_cell_type_distr.csv")
df_ic <- data_summary(propIcsInCluster, varname = "Freq", groupnames = c("Age", "Var1"))
Varx <- factor(c("ICS1","ICS2","ICS3","ICS4","ICS5","ICS6","ICS7","ICS8","ICS1","ICS2","ICS3","ICS4","ICS5","ICS6","ICS7","ICS8"), levels = c("ICS1","ICS2","ICS3","ICS4","ICS5","ICS6","ICS7","ICS8") )
df_ic$Varx <- Varx
ggplot(df_ic, aes(x = Varx, y = Freq, fill = Age))+
  geom_bar(stat = "identity", width = .5, position = position_dodge2(reverse = TRUE))+
  geom_errorbar(aes(ymin= Freq-sd, ymax = Freq+sd), width = .2, position = position_dodge2(reverse = TRUE))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  theme_classic() +
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black")) +
  theme( axis.text.x = element_text(size = 12,  color = "black", angle = 45, vjust = 0.5)) +
  labs(x="",y="Cell type distribution",title="")

# Tsne plots
immune.cells_colorPalette <- c(RColorBrewer::brewer.pal(8,"Set2"))

DimPlot(immune.cells.integrated, reduction = "tsne", group.by = "orig.ident", cols = Unint_patients_cols, pt.size = 1) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3)))

DimPlot(immune.cells.integrated, reduction = "tsne", group.by = "clusters",cols = immune.cells_colorPalette, pt.size = 1) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)))

DimPlot(immune.cells.integrated, reduction = "tsne", group.by = "clusters", split.by = "Age", cols = immune.cells_colorPalette, pt.size = 1) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)))

DimPlot(immune.cells.integrated, reduction = "tsne", group.by = "Age", split.by = "orig.ident", cols = c("skyblue3","sienna3"), pt.size = 0.5) + 
  theme(legend.position = "top") + guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3)))

# Plot known markers
DefaultAssay(immune.cells.integrated) <- "RNA"
DotPlot(immune.cells.integrated, features = c("CD34","CD3D","CD8A","NKG7","CCL5","IL7R","CD69","CXCR4","CD83","CD68","ITGAX","ITGAM","HLA-DQA2","CD163","MRC1","FCGR3A"), group.by = c("seurat_clusters"), cols = c("lightgrey", "red"), dot.scale = 8)+ RotatedAxis() + labs(x = "", y = "")

# Heatmap of marker genes of ICS
immunecell.markers <- FindAllMarkers(immune.cells.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immunecell.markers <- immunecell.markers[which(immunecell.markers$p_val_adj < 0.05),]
ic_markers_top50 <- immunecell.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(immune.cells.integrated, features = ic_markers_top50$gene) + 
  theme(text = element_text(size = 0))
###########################################################################################################################



###########################################################################################################################
# Using SingleR to annotate immune cell data
###########################################################################################################################
hpca.se <- HumanPrimaryCellAtlasData()
bp.se <- BlueprintEncodeData(rm.NA = "rows")

pred.immune.cell <- list()
for (i in 1:8) {
  ed <- as.matrix(GetAssayData(immune.cells.integrated, slot = "counts")[, WhichCells(immune.cells.integrated, ident = as.character(0:7)[i])])
  ed <- SingleCellExperiment(assays = list(counts = ed))
  ed <- logNormCounts(ed)
  pred.immune.cell[[i]] <- SingleR(test = ed, ref = hpca.se, labels = hpca.se$label.fine)
}

pred.immune.cell.C0 <- as.data.frame(table(pred.immune.cell[[1]]$labels))
pred.immune.cell.C1 <- as.data.frame(table(pred.immune.cell[[2]]$labels))
pred.immune.cell.C2 <- as.data.frame(table(pred.immune.cell[[3]]$labels))
pred.immune.cell.C3 <- as.data.frame(table(pred.immune.cell[[4]]$labels))
pred.immune.cell.C4 <- as.data.frame(table(pred.immune.cell[[5]]$labels))
pred.immune.cell.C5 <- as.data.frame(table(pred.immune.cell[[6]]$labels))
pred.immune.cell.C6 <- as.data.frame(table(pred.immune.cell[[7]]$labels))
pred.immune.cell.C7 <- as.data.frame(table(pred.immune.cell[[8]]$labels))
###########################################################################################################################



###########################################################################################################################
#  Find common and unique DEGs between APC or IC clusters
###########################################################################################################################

# Find common DEGs between APC
OA_vs_YA_APC_common_up_DEGs <- getCommonUpDEGs(OA_vs_YA_APC_DEGs)
OA_vs_YA_APC_common_down_DEGs <- getCommonDownDEGs(OA_vs_YA_APC_DEGs)

# Find common DEGs between IC
OA_vs_YA_IC_common_up_DEGs <- getCommonUpDEGs(OA_vs_YA_IC_DEGs)
OA_vs_YA_IC_common_down_DEGs <- getCommonDownDEGs(OA_vs_YA_IC_DEGs)

# Find unique DEGs between APC
OA_vs_YA_APC_unique_up_DEGs <- getUniqueUpGenes(OA_vs_YA_APC_DEGs)
OA_vs_YA_APC_unique_down_DEGs <- getUniqueDownGenes(OA_vs_YA_APC_DEGs)

# Find unique DEGs between IC
OA_vs_YA_IC_unique_up_DEGs <- getUniqueUpGenes(OA_vs_YA_IC_DEGs)
OA_vs_YA_IC_unique_down_DEGs <- getUniqueDownGenes(OA_vs_YA_IC_DEGs)

# plot OA_vs_YA_APC_common_and_unique_up_genes
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
pdf(file = "Results_ApcCellsIntegrated/APC_DEGs_OA_vs_YA/APC_common_and_specific_up_genes.pdf")
pheatmap::pheatmap(num_dec, cluster_rows = F, cluster_cols = F, angle_col = 0, color = colorRampPalette(c("lightgrey", "white", "#F8766C"))(100))
dev.off()


# plot OA_vs_YA_APC_common_and_unique_down_genes
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
pdf(file = "Results_ApcCellsIntegrated/APC_DEGs_OA_vs_YA/APC_common_and_specific_down_genes.pdf")
pheatmap::pheatmap(num_dec, cluster_rows = F, cluster_cols = F, angle_col = 0, color = colorRampPalette(c("lightgrey", "white", "#F8766C"))(100))
dev.off()


# plot OA_vs_YA_IC_common_and_unique_up_genes
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
pdf(file = "Results_ImmuneCellsIntegrated/IC_DEGs_OA_vs_YA/IC_common_and_specific_up_genes.pdf")
pheatmap::pheatmap(num_dec, cluster_rows = F, cluster_cols = F, angle_col = 0, color = colorRampPalette(c("lightgrey", "white", "#F8766C"))(100))
dev.off()


# plot OA_vs_YA_IC_common_and_unique_down_genes
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
pdf(file = "Results_ImmuneCellsIntegrated/IC_DEGs_OA_vs_YA/IC_common_and_specific_down_genes.pdf")
pheatmap::pheatmap(num_dec, cluster_rows = F, cluster_cols = F, angle_col = 0, color = colorRampPalette(c("lightgrey", "white", "#F8766C"))(100))
dev.off()
###########################################################################################################################


