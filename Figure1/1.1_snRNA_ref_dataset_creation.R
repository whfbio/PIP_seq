# ------------------------------------------------------------------------------
# Title: Preparation of Reference Seurat Datasets for Cytospace
# Author: Haofei Wang
# Date: March 24, 2025
# Description:
# This script prepares Seurat objects for reference datasets across timepoints 
# (p0, p7, p14, p21). It reads 10X Genomics data, creates and merges Seurat objects 
# for replicates, performs quality control, and conducts PCA, clustering, and UMAP 
# visualization. These processed objects can then be used for Cytospace or downstream analysis.
#
# Dependencies:
# - Seurat
# - Harmony
# - tidyverse
# - cowplot
#
# Usage:
# Update the working directory and 10X input paths before running the script.
# ------------------------------------------------------------------------------


library(Seurat)
library(harmony)
library(tidyverse)
library(cowplot)

#set wrkdr where "ref_dataset" is to read10X files with seurat

p21_1_inputdata <- Read10X("./ref_dataset/P21__1_CCCACCACAA-AAGCGGAGGT/outs/filtered_feature_bc_matrix")
p21_2_inputdata <- Read10X("./ref_dataset/P21__2_CGGCTGGATG-GTGCTTATCA/outs/filtered_feature_bc_matrix")

p14_1_inputdata <- Read10X("./ref_dataset/P14__1_AAGGGCCGCA-GAGGAATCAG/outs/filtered_feature_bc_matrix")
p14_2_inputdata <- Read10X("./ref_dataset/P14__2_GAGAGGATAT-CCCATTTCAA/outs/filtered_feature_bc_matrix")

p7_1_inputdata  <- Read10X("./ref_dataset/batch2/P7__1_CCGGCAACTG-CGGTTTAACA/outs/filtered_feature_bc_matrix" )
p7_2_inputdata <- Read10X("./ref_dataset/batch2/P7__2_TTCACACCTT-TAGTGTACAC/outs/filtered_feature_bc_matrix" )

p0_1_inputdata  <- Read10X("./ref_dataset/P0__1_TTGCCCGTGC-AATCTCACGC/outs/filtered_feature_bc_matrix" )
p0_2_inputdata  <- Read10X("./ref_dataset/P0__2_AATGTATCCA-TAAGCTCATT/outs/filtered_feature_bc_matrix" )

#now that files are read we will create seurat object with them
p21_1 <- CreateSeuratObject(counts = p21_1_inputdata, project = "p21_1")
p21_2 <- CreateSeuratObject(counts = p21_2_inputdata, project = "p21_2")

p14_1 <- CreateSeuratObject(counts = p14_1_inputdata, project = "p14_1")
p14_2 <- CreateSeuratObject(counts = p14_2_inputdata, project = "p14_2")

p7_1  <- CreateSeuratObject(counts = p7_1_inputdata_v2, project = "p7_1_v2")
p7_2 <- CreateSeuratObject(counts  = p7_2_inputdata_v2, project = "p7_2_v2")

p0_1  <- CreateSeuratObject(counts = p0_1_inputdata, project = "p0_1")
p0_2  <- CreateSeuratObject(counts = p0_2_inputdata, project = "p0_2")
#
#
#
#
#
########merge replicates to evaluate if replicate batch effects persists
#
p21 <- merge(p21_1, y = p21_2, add.cell.ids = c("p21_1", "p21_2"), project = "p21")

p14 <- merge(p14_1, y = p14_2, add.cell.ids = c("p14_1", "p14_2"), project = "p14")

p7 <- merge(p7_1, y = p7_2, add.cell.ids = c("p7_1", "p7_2"), project = "p7")

p0  <- merge(p0_1, y = p0_2, add.cell.ids = c("p0_1", "p0_2"), project = "p0")

################# QC Section ################# 
# look at distribution of counts in cells and perform QC
p21[["percent.mt"]] <- PercentageFeatureSet(p21, pattern = "^mt-") #add mito percentages per cell

pdf("./QC_folder/p21_QC/p21_QC_pre_vs_post_filter.pdf", onefile = TRUE)
VlnPlot(p21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(p21, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p21, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

p21 <- subset(p21, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 5)


VlnPlot(p21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(p21, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p21, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

p14[["percent.mt"]] <- PercentageFeatureSet(p14, pattern = "^mt-") #add mito percentages per cell

pdf("./QC_folder/p14_QC/p14_QC_pre_vs_post_filter.pdf", onefile = TRUE)
VlnPlot(p14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(p14, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

p14 <- subset(p14, subset = nFeature_RNA > 750 & nFeature_RNA < 13000 & percent.mt < 5)

VlnPlot(p14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(p14, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

p7[["percent.mt"]] <- PercentageFeatureSet(p7, pattern = "^mt-") #add mito percentages per cell
pdf("./QC_folder/p7_QC/p7_QC_pre_vs_post_filter.pdf", onefile = TRUE)
VlnPlot(p7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(p7, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


p7 <- subset(p7, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 5)


VlnPlot(p7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(p7, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()



p0[["percent.mt"]] <- PercentageFeatureSet(p0, pattern = "^mt-") #add mito percentages per cell

pdf("./QC_folder/p0_QC/p0_QC_pre_vs_post_filter.pdf", onefile = TRUE)
VlnPlot(p0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(p0, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


p0 <- subset(p0, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 5)

VlnPlot(p0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(p0, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

#Next step is to perform clustering

pdf("./QC_folder/p21_QC/p21_Clustering_QC.pdf", onefile = TRUE)
p21 <- p21 %>% NormalizeData() %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(p21), 10)
LabelPoints(plot = VariableFeaturePlot(p21), points = top10, repel = TRUE)

all.genes <- rownames(p21)
p21 <- ScaleData(p21, features = all.genes)
p21 <- RunPCA(p21, features = VariableFeatures(object = p21))

DimPlot(p21, reduction = "pca", group.by = "orig.ident")
ElbowPlot(p21) #choose 15

p21 <- FindNeighbors(p21, dims = 1:15)
p21 <- FindClusters(p21, resolution = .8)
p21 <- RunUMAP(p21, dims = 1:15)

DimPlot(p21, reduction = "umap", label = TRUE)
DimPlot(p21, reduction = "umap", label = TRUE, group.by = "orig.ident")
dev.off()


p14 <- p14 %>% NormalizeData() %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(p14), 10)
LabelPoints(plot = VariableFeaturePlot(p14), points = top10, repel = TRUE)

all.genes <- rownames(p14)
p14 <- ScaleData(p14, features = all.genes)
p14 <- RunPCA(p14, features = VariableFeatures(object = p14))

DimPlot(p14, reduction = "pca", group.by = "orig.ident")
ElbowPlot(p14) #choose 15

p14 <- FindNeighbors(p14, dims = 1:15)
p14 <- FindClusters(p14, resolution = .8)
p14 <- RunUMAP(p14, dims = 1:15)

pdf("./QC_folder/p14_QC/p14_Clustering_QC.pdf", onefile = TRUE)
LabelPoints(plot = VariableFeaturePlot(p14), points = top10, repel = TRUE)
DimPlot(p14, reduction = "pca", group.by = "orig.ident")
ElbowPlot(p14) #choose 15
DimPlot(p14, reduction = "umap", label = TRUE)
DimPlot(p14, reduction = "umap", label = TRUE, group.by = "orig.ident")
dev.off()



p7 <- p7 %>% NormalizeData() %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(p7), 10)
pdf("top_variable_genes.pdf")
LabelPoints(plot = VariableFeaturePlot(p7), points = top10, repel = TRUE)
dev.off()


all.genes <- rownames(p7)
p7 <- ScaleData(p7, features = all.genes)
p7 <- RunPCA(p7, features = VariableFeatures(object = p7))


DimPlot(p7, reduction = "pca", group.by = "orig.ident")

ElbowPlot(p7) #choose 15

p7 <- FindNeighbors(p7, dims = 1:15)
p7 <- FindClusters(p7, resolution = 0.8)
p7 <- RunUMAP(p7, dims = 1:15)

pdf("./QC_folder/p7_QC/p7_Clustering_QC.pdf", onefile = TRUE)
LabelPoints(plot = VariableFeaturePlot(p7), points = top10, repel = TRUE)
DimPlot(p7, reduction = "pca", group.by = "orig.ident")
ElbowPlot(p7) #choose 15
DimPlot(p7, reduction = "umap", label = TRUE)
DimPlot(p7, reduction = "umap", label = TRUE, group.by = "orig.ident")
dev.off()


p0 <- p0 %>% NormalizeData() %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(p0), 10)
LabelPoints(plot = VariableFeaturePlot(p0), points = top10, repel = TRUE)

all.genes <- rownames(p0)
p0 <- ScaleData(p0, features = all.genes)
p0 <- RunPCA(p0, features = VariableFeatures(object = p0))

DimPlot(p0, reduction = "pca", group.by = "orig.ident")
ElbowPlot(p0) #choose 15

p0 <- FindNeighbors(p0, dims = 1:15)
p0 <- FindClusters(p0, resolution = 0.8)
p0 <- RunUMAP(p0, dims = 1:15)

pdf("p0_Clustering_QC.pdf", onefile = TRUE)
LabelPoints(plot = VariableFeaturePlot(p0), points = top10, repel = TRUE)
DimPlot(p0, reduction = "pca", group.by = "orig.ident")
ElbowPlot(p0) #choose 15
DimPlot(p0, reduction = "umap", label = TRUE)
DimPlot(p0, reduction = "umap", label = TRUE, group.by = "orig.ident")
dev.off()

#########################################################
#
#
#Here in this section we have to remove and document the removal of doublets
#
#
#########################################################
library(scDblFinder)

dbl.dens.p21 <- computeDoubletDensity(p21@assays$RNA@counts, dims=15)
p21$DoubletScore <- dbl.dens.p21
DimPlot(p21, split.by = "DoubletScore")

dbl.calls.p21 <- doubletThresholding(data.frame(score=dbl.dens.p21),
                                     method="griffiths", returnType="call")

summary(dbl.calls.p21)

p21$DoubletScoreOutlier <- dbl.calls.p21

pdf("./QC_folder/p21_QC/Doublets_Score_QC_p21.pdf", onefile = TRUE, width = 12)
DimPlot(p21, group.by = "seurat_clusters", label = TRUE)
FeaturePlot(p21, features = "DoubletScore", label = TRUE)
DimPlot(p21, group.by = "DoubletScoreOutlier")
DimPlot(p21, group.by ="seurat_clusters", split.by  = "DoubletScoreOutlier", label = TRUE)
VlnPlot(p21, group.by = "seurat_clusters", features = "DoubletScore" )
VlnPlot(p21, group.by = "seurat_clusters", features = "DoubletScore", split.by = "DoubletScoreOutlier" )
dev.off()


dbl.dens.p14 <- computeDoubletDensity(p14@assays$RNA@counts, dims=15)
p14$DoubletScore <- dbl.dens.p14
DimPlot(p14, split.by = "DoubletScore")

dbl.calls.p14 <- doubletThresholding(data.frame(score=dbl.dens.p14),
                                     method="griffiths", returnType="call")
summary(dbl.calls.p14)

p14$DoubletScoreOutlier <- dbl.calls.p14

pdf("./QC_folder/p14_QC/Doublets_Score_QC_p14.pdf", onefile = TRUE, width = 12)
DimPlot(p14, group.by = "seurat_clusters", label = TRUE)
FeaturePlot(p14, features = "DoubletScore", label = TRUE)
DimPlot(p14, group.by = "DoubletScoreOutlier")
DimPlot(p14, group.by ="seurat_clusters", split.by  = "DoubletScoreOutlier", label = TRUE)
VlnPlot(p14, group.by = "seurat_clusters", features = "DoubletScore" )
VlnPlot(p14, group.by = "seurat_clusters", features = "DoubletScore", split.by = "DoubletScoreOutlier" )
dev.off()


dbl.dens.p7 <- computeDoubletDensity(p7@assays$RNA@counts, dims=15)
p7$DoubletScore <- dbl.dens.p7

dbl.calls.p7 <- doubletThresholding(data.frame(score=dbl.dens.p7),
                                       method="griffiths", returnType="call")
summary(dbl.calls.p7)

p7$DoubletScoreOutlier <- dbl.calls.p7

pdf("./QC_folder/p7_QC/Doublets_Score_QC_p7.pdf", onefile = TRUE, width = 12)
DimPlot(p7, group.by = "seurat_clusters", label = TRUE)
FeaturePlot(p7, features = "DoubletScore", label = TRUE)
DimPlot(p7, group.by = "DoubletScoreOutlier")
DimPlot(p7, group.by ="seurat_clusters", split.by  = "DoubletScoreOutlier", label = TRUE)
VlnPlot(p7, group.by = "seurat_clusters", features = "DoubletScore" )
VlnPlot(p7, group.by = "seurat_clusters", features = "DoubletScore", split.by = "DoubletScoreOutlier" )
dev.off()


dbl.dens.p0 <- computeDoubletDensity(p0@assays$RNA@counts, dims=15)
p0$DoubletScore <- dbl.dens.p0
DimPlot(p0, split.by = "DoubletScore")

dbl.calls.p0 <- doubletThresholding(data.frame(score=dbl.dens.p0),
                                    method="griffiths", returnType="call")
summary(dbl.calls.p0)

p0$DoubletScoreOutlier <- dbl.calls.p0

pdf("./Doublets/Doublets_Score_QC_p0.pdf", onefile = TRUE, width = 12)
DimPlot(p0, group.by = "seurat_clusters", label = TRUE)
FeaturePlot(p0, features = "DoubletScore", label = TRUE)
DimPlot(p0, group.by = "DoubletScoreOutlier")
DimPlot(p0, group.by ="seurat_clusters", split.by  = "DoubletScoreOutlier", label = TRUE)
VlnPlot(p0, group.by = "seurat_clusters", features = "DoubletScore" )
VlnPlot(p0, group.by = "seurat_clusters", features = "DoubletScore", split.by = "DoubletScoreOutlier" )
dev.off()

#Finally save these objects
# saveRDS(p0, "p0_predoublets.rds")
# saveRDS(p7_v1, "p7_v1_predoublets.rds")
# saveRDS(p7_v2, "p7_v2_predoublets.rds")
# saveRDS(p14, "p14_predoublets.rds")
# saveRDS(p21, "p21_predoublets.rds")


#save the doublet finder output here when you can 

p21_predoublets <- readRDS("/proj/liulab/users/mdcolon/Ref_Data_Analysis_PerturbSeq/Data_Set_Creation/rObjs/Before_Doublet_Removal/p21_predoublets.rds")
p14_predoublets <- readRDS("/proj/liulab/users/mdcolon/Ref_Data_Analysis_PerturbSeq/Data_Set_Creation/rObjs/Before_Doublet_Removal/p14_predoublets.rds")

p7_predoublets <- readRDS("./rObjs/Before_Doublet_Removal/p7_new_predoublets.rds")

p0_predoublets <- readRDS("/proj/liulab/users/mdcolon/Ref_Data_Analysis_PerturbSeq/Data_Set_Creation/rObjs/Before_Doublet_Removal/p0_predoublets.rds")

#In the below section I finally perform SCTranform and remove doublets and more mitochondria...

###NOTE doublet score actually doubletThresholding() function will generate different values every time, to reproduce this analysis you can store the output that was first calculated
####fir kater and then filter based on those was determined by decreasing threshold until I see minimal clustering of data w.r.t DoubletScore
####note that doubletscore had to be made more restricteve than what the software scDblFinder package determined

##
p21    <- subset(p21_predoublets, subset = DoubletScoreOutlier =="singlet" & percent.mt < 1) #single is <= .800
p14    <- subset(p14_predoublets, subset = DoubletScore <.9 & percent.mt <.93 )
p7     <- subset(p7_predoublets,  subset = DoubletScore <.5 & percent.mt <1) 
p0     <- subset(p0_predoublets, subset = DoubletScore < 1.2 & percent.mt < 1)
##

p21 <- SCTransform(p21, method = "glmGamPoi")
p21 <- RunPCA(p21)
p21 <- RunUMAP(p21, dims = 1:15)
p21 <- FindNeighbors(p21, dims = 1:15)
p21 <- FindClusters(p21, resolution = .8)
DimPlot(p21, label = TRUE)
FeaturePlot(p21, c("percent.mt", "Ttn", "DoubletScore"))
#p21[["celltype.lv3"]] <- p21_celltype.lv3[names(p21_celltype.lv3) %in%  colnames(p21)]
#DimPlot(p21, label = TRUE, group.by = "celltype.lv3")

p21 <- subset(p21, subset = seurat_clusters != c(16)) #this cluster is in both CMs and Fib and it is most likely a doublets cluster. 
p21 <- SCTransform(p21, method = "glmGamPoi")
p21 <- RunPCA(p21)
p21 <- RunUMAP(p21, dims = 1:15)
p21 <- FindNeighbors(p21, dims = 1:30)
p21 <- FindClusters(p21, resolution = .8)
DimPlot(p21, label = TRUE)
FeaturePlot(p21, c("percent.mt", "Ttn", "DoubletScore"))            

p14 <- SCTransform(p14, method = "glmGamPoi")
p14 <- RunPCA(p14)
p14 <- RunUMAP(p14, dims = 1:15)
p14 <- FindNeighbors(p14, dims = 1:15)
p14 <- FindClusters(p14, resolution = .8)
DimPlot(p14, label = TRUE)
FeaturePlot(p14, c("percent.mt", "Ttn", "DoubletScore"))

#p14$remove <- Cells(p14) %in% names(remove_p14)
#p14$remove <- Cells(p14) %in% names(remove_p14)
p14$remove <- remove_p14
p14 <- subset(p14, subset = remove == FALSE)
p14 <- SCTransform(p14, method = "glmGamPoi")
p14 <- RunPCA(p14)
p14 <- RunUMAP(p14, dims = 1:15)
p14 <- FindNeighbors(p14, dims = 1:30)
p14 <- FindClusters(p14, resolution = 1)
DimPlot(p14, label = TRUE)
FeaturePlot(p14, c("percent.mt","Mki67" ))
FeaturePlot(p14, c("percent.mt", "Ttn", "DoubletScore"))


DimPlot(p14, group.by = "celltype.lv2")

#p7 is slightly different because there are to samples of p7 that need to be combined
#I combine them with harmony
p7 <- SCTransform(p7, method = "glmGamPoi")
p7 <- RunPCA(p7)
p7 <- RunUMAP(p7, dims = 1:15)
p7 <- FindNeighbors(p7, dims = 1:15)
p7 <- FindClusters(p7, resolution = .8)

DimPlot(p7, label = TRUE)
FeaturePlot(p7, "DoubletScore")


p0 <- SCTransform(p0, method = "glmGamPoi")
p0 <- RunPCA(p0)
p0 <- RunUMAP(p0, dims = 1:15)
p0 <- FindNeighbors(p0, dims = 1:15)
p0 <- FindClusters(p0, resolution = .8)
FeaturePlot(p0, c("percent.mt", "Ttn", "DoubletScore"))

p0  <- subset(p0, subset = seurat_clusters != c(17)) #this clusters looks like doublets between a_CM and v_CM as it lacks Tbx3(av Marker) but having markers for both 
p0 <- SCTransform(p0, method = "glmGamPoi")
p0 <- RunPCA(p0)
p0 <- RunUMAP(p0, dims = 1:15)
p0 <- FindNeighbors(p0, dims = 1:15)
p0 <- FindClusters(p0, resolution = .8)
DimPlot(p0, label = TRUE) + FeaturePlot(p0, c("DoubletScore", "Ttn"))


FeaturePlot(p21,      c("percent.mt", "DoubletScore"))
FeaturePlot(p14,      c("percent.mt", "DoubletScore"))
FeaturePlot(p7,       c("percent.mt", "DoubletScore"))
FeaturePlot(p0,       c("percent.mt", "DoubletScore"))


#we add the time_point metadata here to facilitate future plotting 
p0$time_point  <- "p0"
p7$time_point  <- "p7"
p14$time_point <- "p14"
p21$time_point <- "p21"

#I save as sct applied to signify they have been processed with SCT and filter and now only need to be manually annotated
#the finale annotated version are created in Manual_annotations script.
saveRDS(p0,  file = "./rObjs/SCT_applied/p0_sct.rds")
saveRDS(p7,  file = "./rObjs/SCT_applied/p7_sct.rds")
saveRDS(p14,  file = "./rObjs/SCT_applied/p14_sct.rds")
saveRDS(p21,  file = "./rObjs/SCT_applied/p21_sct.rds")





