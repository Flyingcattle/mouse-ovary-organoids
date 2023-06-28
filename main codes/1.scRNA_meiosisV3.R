
################################################# Seurat3.0 pipline #########################################
############################################# used for scRNA seq _ ovary-organoids #####################################


### load required packages 
library(Seurat)
library(dplyr)
library(cowplot)
library(DoubletFinder)
library(ggplot2)


### read inputdata 
CON_11.data <- Read10X(data.dir = "E:/Lab of Germ Cell Bio/10xGenomics_meiosis_V3/scRNA_meiosis/cellranger_output/GR11/filtered_feature_bc_matrix")
CON_12.data <- Read10X(data.dir = "E:/Lab of Germ Cell Bio/10xGenomics_meiosis_V3/scRNA_meiosis/cellranger_output/GR12/filtered_feature_bc_matrix")
CON_13.data <- Read10X(data.dir = "E:/Lab of Germ Cell Bio/10xGenomics_meiosis_V3/scRNA_meiosis/cellranger_output/GR13/filtered_feature_bc_matrix")


colnames(x = CON_11.data) <- paste('GR11', colnames(x = CON_11.data), sep = '_')
colnames(x = CON_12.data) <- paste('GR12', colnames(x = CON_12.data), sep = '_')
colnames(x = CON_13.data) <- paste('GR13', colnames(x = CON_13.data), sep = '_')


## Construct Seurat object
CON_11 <- CreateSeuratObject(counts = CON_11.data, project = "GR11", min.cells = 3, min.features = 200)
CON_11[["percent.mt"]]<-PercentageFeatureSet(CON_11,pattern = "^Mt")
VlnPlot(CON_11, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) ## 5x3
CON_11 <- subset(CON_11, subset =  nFeature_RNA > 2000 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >200)
CON_11 <- NormalizeData(CON_11, normalization.method = "LogNormalize", scale.factor = 10000)
CON_11 <- FindVariableFeatures(CON_11, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
CON_11$group<-"GR11"
all.genes <- rownames(CON_11)
CON_11 <- ScaleData(CON_11, features = all.genes)
CON_11@meta.data$time <- "GR11"

CON_12 <- CreateSeuratObject(counts = CON_12.data, project = "GR12", min.cells = 3, min.features = 200)
CON_12[["percent.mt"]]<-PercentageFeatureSet(CON_12,pattern = "^Mt")
VlnPlot(CON_12, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) ## 5x3
CON_12 <- subset(CON_12, subset =  nFeature_RNA > 2000 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >200)
CON_12 <- NormalizeData(CON_12, normalization.method = "LogNormalize", scale.factor = 10000)
CON_12 <- FindVariableFeatures(CON_12, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
CON_12$group<-"GR12"
all.genes <- rownames(CON_12)
CON_12 <- ScaleData(CON_12, features = all.genes)
CON_12@meta.data$time <- "GR12"

CON_13 <- CreateSeuratObject(counts = CON_13.data, project = "GR13", min.cells = 3, min.features = 200)
CON_13[["percent.mt"]]<-PercentageFeatureSet(CON_13,pattern = "^Mt")
VlnPlot(CON_13, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) ## 5x3
CON_13 <- subset(CON_13, subset =  nFeature_RNA > 2000 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >200)
CON_13 <- NormalizeData(CON_13, normalization.method = "LogNormalize", scale.factor = 10000)
CON_13 <- FindVariableFeatures(CON_13, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
CON_13$group<-"GR13"
all.genes <- rownames(CON_13)
CON_13 <- ScaleData(CON_13, features = all.genes)
CON_13@meta.data$time <- "GR13"



#################################### remove douplets cells with Doublefinder #################################
################################ https://www.jianshu.com/p/b1947c4156ad   ####################################
#########################   https://github.com/chris-mcginnis-ucsf/DoubletFinder #############################

CON_11<-RunPCA(CON_11)
CON_11<-RunUMAP(CON_11,dims = 1:10)
## pK Identification (no ground-truth)-----------------------------------------------------------------------------------------
sweep.res.list_CON_11 <- paramSweep_v3(CON_11, PCs = 1:20, sct = FALSE)
head(sweep.res.list_CON_11)
sweep.stats_CON_11 <- summarizeSweep(sweep.res.list_CON_11, GT = FALSE)
bcmvn_CON12 <- find.pK(sweep.stats_CON_11)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_CON_11 <- paramSweep_v3(CON_11, PCs = 1:20, sct = FALSE)
# gt.calls <- CON_11@meta.data[rownames(sweep.res.list_CON_11[[1]]), "GT"]
# sweep.stats_CON_11 <- summarizeSweep(sweep.res.list_CON_11, GT = TRUE, GT.calls = gt.calls)
# bcmvn_CON12 <- find.pK(sweep.stats_CON_11)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
CON_11<-FindNeighbors(CON_11,reduction ="pca",dims = 1:20 )
CON_11<-FindClusters(CON_11,resolution = 0.5)
head(CON_11@meta.data)
annotationscon_E12<-CON_11@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_E12)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.09*(length(CON_11@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_con_E11 <- doubletFinder_v3(CON_11, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(seu_con_E11@meta.data)
# seu_con_E11@meta.data$DF.classifications_0.25_0.09_1413
table(seu_con_E11$DF.classifications_0.25_0.09_501)
## Doublet Singlet 
## 501    5064

seu_con_E11@meta.data$cellfilter <- seu_con_E11@meta.data$DF.classifications_0.25_0.09_501
seu_con_E11@meta.data <-seu_con_E11@meta.data[,-10]

# Doublet Singlet 
# 523    6453
seu_con_E11@meta.data$time<- "GR11"


CON_12<-RunPCA(CON_12)
CON_12<-RunUMAP(CON_12,dims = 1:10)
## pK Identification (no ground-truth)-----------------------------------------------------------------------------------------
sweep.res.list_CON_12 <- paramSweep_v3(CON_12, PCs = 1:20, sct = FALSE)
head(sweep.res.list_CON_12)
sweep.stats_CON_12 <- summarizeSweep(sweep.res.list_CON_12, GT = FALSE)
bcmvn_CON14 <- find.pK(sweep.stats_CON_12)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_CON_12 <- paramSweep_v3(CON_12, PCs = 1:20, sct = FALSE)
# gt.calls <- CON_12@meta.data[rownames(sweep.res.list_CON_12[[1]]), "GT"]
# sweep.stats_CON_12 <- summarizeSweep(sweep.res.list_CON_12, GT = TRUE, GT.calls = gt.calls)
# bcmvn_CON12 <- find.pK(sweep.stats_CON_12)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
CON_12<-FindNeighbors(CON_12,reduction ="pca",dims = 1:20 )
CON_12<-FindClusters(CON_12,resolution = 0.5)
head(CON_12@meta.data)
annotationscon_E14<-CON_12@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_E14)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.09*(length(CON_12@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_con_E12 <- doubletFinder_v3(CON_12, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# seu_con_E12@meta.data$DF.classifications_0.25_0.09_1413
head(seu_con_E12@meta.data)

table(seu_con_E12$DF.classifications_0.25_0.09_612)
## Doublet Singlet 
##612    6189
seu_con_E12@meta.data$cellfilter <- seu_con_E12@meta.data$DF.classifications_0.25_0.09_612
seu_con_E12@meta.data <-seu_con_E12@meta.data[,-10]

seu_con_E12@meta.data$time<- "GR12"


CON_13<-RunPCA(CON_13)
CON_13<-RunUMAP(CON_13,dims = 1:10)
## pK Identification (no ground-truth)-----------------------------------------------------------------------------------------
sweep.res.list_CON_13 <- paramSweep_v3(CON_13, PCs = 1:20, sct = FALSE)
head(sweep.res.list_CON_13)
sweep.stats_CON_13 <- summarizeSweep(sweep.res.list_CON_13, GT = FALSE)
bcmvn_DEHP12 <- find.pK(sweep.stats_CON_13)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_CON_13 <- paramSweep_v3(CON_13, PCs = 1:20, sct = FALSE)
# gt.calls <- CON_13@meta.data[rownames(sweep.res.list_CON_13[[1]]), "GT"]
# sweep.stats_CON_13 <- summarizeSweep(sweep.res.list_CON_13, GT = TRUE, GT.calls = gt.calls)
# bcmvn_DEHP12 <- find.pK(sweep.stats_CON_13)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
CON_13<-FindNeighbors(CON_13,reduction ="pca",dims = 1:20 )
CON_13<-FindClusters(CON_13,resolution = 0.5)
head(CON_13@meta.data)
annotationsDEHP_E12<-CON_13@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationsDEHP_E12)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.09*(length(CON_13@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_con_E13 <- doubletFinder_v3(CON_13, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# seu_con_E13@meta.data$DF.classifications_0.25_0.09_1413
head(seu_con_E13@meta.data)
table(seu_con_E13$DF.classifications_0.25_0.09_436)
# Doublet Singlet 
#  436    4408
seu_con_E13@meta.data$cellfilter <- seu_con_E13@meta.data$DF.classifications_0.25_0.09_436
seu_con_E13@meta.data <-seu_con_E13@meta.data[,-10]
seu_con_E13@meta.data$time<- "GR13"


###################################### I need to remove the douplets  #######################################
#########################################  Extract the singlets  ############################################
head(seu_con_E12@meta.data)
DimPlot(seu_con_E11,group.by = "cellfilter") # E11
DimPlot(seu_con_E12,group.by = "cellfilter") # E12
DimPlot(seu_con_E13,group.by = "cellfilter") # E13


head(seu_con_E11@meta.data)
CON11.singlet <- subset(seu_con_E11, subset = cellfilter == 'Singlet')

head(seu_con_E12@meta.data)
CON12.singlet <- subset(seu_con_E12, subset = cellfilter == 'Singlet')

head(seu_con_E13@meta.data)
CON13.singlet <- subset(seu_con_E13, subset = cellfilter == 'Singlet')


## Here, add some UMAP plots to show the doublets dstribution is desired 

########################################### end doublefinder  ################################################
################################ https://www.jianshu.com/p/b1947c4156ad   ####################################
#########################   https://github.com/chris-mcginnis-ucsf/DoubletFinder #############################



################################################# Perform downstream integration  ########################################################
#################################################    Standard workflow    ###################################################

GRs.anchors <- FindIntegrationAnchors(object.list = list(CON11.singlet, CON12.singlet, CON13.singlet), dims = 1:20)  ## NS
GR.combined <- IntegrateData(anchorset = GRs.anchors, dims = 1:20)

saveRDS(GR.combined, file = "E:\\Lab of Germ Cell Bio\\10xGenomics_meiosis_V3\\scRNA_meiosis\\GR_V3.combined.rds")
# setwd("E:\\Lab of Germ Cell Bio\\10xGenomics_DEHP\\SeuratV3")
# save.image("GR.combined.RData")

GR.combined <- readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_meiosis_V3\\scRNA_meiosis\\GR_V3.combined.rds")

## Perform an integrated analysis
###   We use the integrated assay to jointly define cell types in stimulated/control cells,
###   and the RNA assay to define markers and cell-type specific responses.
###   You should use the integrated assay when trying to 'align' cell states that are 
###   shared across datasets (i.e. for clustering, visualization, learning pseudotime, etc.)
###   You should use the RNA assay when exploring the genes that change either across 
###   clusters, trajectories, or conditions.

DefaultAssay(GR.combined) <- "integrated"
### PCA
GR.combined <- ScaleData(GR.combined, verbose = FALSE)
GR.combined <- RunPCA(GR.combined, verbose = FALSE,npcs = 30)
ElbowPlot(GR.combined)   ## Here I used 10 pcs for downstream analysis

### visulization
GR.combined <- RunUMAP(GR.combined,reduction = "pca", dims = 1:15)
GR.combined <- FindNeighbors(GR.combined, reduction = "pca", dims = 1:15)
GR.combined <- FindClusters(GR.combined, resolution = 0.6)

### Evaluate the undesired clusters 
VlnPlot(object = GR.combined, features = c("nFeature_RNA", "nCount_RNA"))  ## doublets? 

### set the defined color 
### Dimplot to show the clusters 
my_color=c("#1D77B2","#ACC8E6", "#FD7F0D", "#FEBA73", "#2CA029", 
           "#97E083", "#D7262B", "#FD9797", "#9667BB",
           "#C4B1D3", "#8B564C", "#7859A3", "#E376BE",
           "#EEB9D1", "#929597","#4DC1B4")
             # "#4DC1B4", "#EE82EE", "#40E0D0","#DB7093")

##  color 2 : 1D77B2 ACC8E6 FD7F0D FEBA73 2CA029 97E083 D7262B FD9797 9667BB 
##            C4B1D3  8B564C  C39C98 E376BE  EEB9D1 929597 CCCECE A7A123
my_color=c("#1D77B2","#ACC8E6", "#FD7F0D", "#FEBA73", "#2CA029", 
             "#97E083", "#D7262B", "#FD9797", "#9667BB",
             "#C4B1D3", "#8B564C", "#7859A3", "#E376BE",
             "#EEB9D1", "#929597", "#CCCECE", "#A7A123","#DB7093")
my_color2=c("#6E4B9E", "#FEE500", "#C06CAB", "#3BBCA8", 
           "#89288F", "#8A9FD1", "#F37B7D", "#D51F26",
           "#D8A767", "#D24B27", "#7859A3", "#90D5E4",
           "#89C75F", "#D2691E","#00BFFF","#800000","#1E90FF","#DB7093")


DimPlot(GR.combined, reduction = "umap", cols =my_color, label = FALSE)+
                   theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20)) # 6x5
DimPlot(GR.combined, reduction = "umap",  cols =my_color, label = TRUE)+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20)) # 6x5
DimPlot(GR.combined, reduction = "umap",  cols =my_color, label = TRUE, split.by = "time")+
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20)) # 6x5


DimPlot(GR.combined, reduction = "umap", cols =colors.use, label = TRUE)+
                   theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20)) # 6x5
DimPlot(GR.combined, reduction = "umap", cols = c("#FD7F0D","#9667BB","#6495ED"), group.by = "time")
DimPlot(GR.combined, reduction = "umap")

head(GR.combined@meta.data)
table(GR.combined@meta.data$seurat_clusters)

################################################ Determine the cluster identity  ##############################################
################################################ Determine the cluster identity  ##############################################

DefaultAssay(GR.combined) <- "RNA"
# nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
# head(nk.markers)
head(GR.combined@meta.data)

### find markers for every cluster compared to all remaining cells, report only the positive ones
GR.combined@misc$markers <- FindAllMarkers(object = GR.combined, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
write.table(GR.combined@misc$markers,file='E:\\Lab of Germ Cell Bio\\10xGenomics_meiosis_V3\\scRNA_meiosis\\All_cluster_markers.txt',row.names = FALSE,quote = FALSE,sep = '\t')

## markers used 
## Germ cell: Dazl, Stra8, Sycp3             cluster 3
## Granulosa cells : Wnt4, Wnt6              cluster 0,1,2,7,9
## Mesothelial:    Lhx9 , Upk3b              cluster 4,8,10,11
## Interstitial:   Bgn  Col1a2               cluster 5,6,12
## endothelial:    Pecam1, Kdr               cluster 13
## Erythriod:      Alad  Alas2               cluster 14
## Immune :        Cd52  Car2                cluster 15

FeaturePlot(GR.combined, features = c("Ddx4"), cols = c("#C0C0C0","#FF4500"), min.cutoff = "q9") +
   theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20)) ## set axis font size



############################### ?????? 2021.5.12 Stacked violin plot ?????? ##########################################
############################### ?????? 2021.5.12 Stacked violin plot ?????? ##########################################

library(Seurat)
library(ggplot2)

modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
StackedVlnPlot(GR.combined, c('Sycp3', 'Stra8', 'Dazl', 'Pecam1', 'Kdr'), pt.size=0, cols=my_color)

## I need to order the clusters 
## markers used 
## Germ cell: Dazl, Stra8, Sycp3             cluster 3
## Granulosa cells : Wnt4, Wnt6              cluster 0,1,2,7,9
## Mesothelial:    Lhx9 , Upk3b              cluster 4,8,10,11
## Interstitial:   Bgn  Col1a2               cluster 5,6,12
## endothelial:    Pecam1, Kdr               cluster 13
## Erythriod:      Alad  Alas2               cluster 14
## Immune :        Cd52  Car2                cluster 15

# Define order of appearance of timepoints
my_levels <- c("3", "0", "1", "2", "7", "9", "4","8","10","11","5","6","12","13","14","15")
# my_levels <- c(2,4,5,12,17,0,1,13,3,9,10,6,7,11,14,15,8,16)
head(GR.combined@meta.data)

## relevel the cluster identity
GR.combined@meta.data$seurat_clusters <- factor(x = GR.combined@meta.data$integrated_snn_res.0.6, levels = my_levels)

StackedVlnPlot(GR.combined, c('Dazl', 'Ddx4', 
                              'Wnt4', 'Wnt6',
                              "Krt19","Upk3b",
                              "Bgn","Kdr",
                              "Alas2","Cd52"
), pt.size=0, cols=my_color, group.by = "seurat_clusters")


############################### ?????? 2021.5.12 Stacked violin plot ?????? ##########################################
############################### ?????? 2021.5.12 Stacked violin plot ?????? ##########################################


########################### ?????? 2020.9.3 Extract germ cells ?????? ##########################################
########################### ?????? 2020.9.3 Extract germ cells ?????? ##########################################
head(GR.combined@meta.data)
germsub = GR.combined[,GR.combined@meta.data$seurat_clusters %in% c(3)]
head(germsub@meta.data)

sce=germsub
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 

sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 0.5 )
# Look at cluster IDs of the first 5 cells
head(Idents(sce), 5)
table(sce$seurat_clusters) 
sce <- RunUMAP(sce, dims = 1:10)
DimPlot(sce, reduction = 'umap')

sce <- subset(sce, idents = c("6","7"), invert = TRUE)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 

sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 0.5 )
# Look at cluster IDs of the first 5 cells
head(Idents(sce), 5)
table(sce$seurat_clusters) 
sce <- RunUMAP(sce, dims = 1:10)
DimPlot(sce, reduction = 'umap')


genes_to_check = c('Dazl', 'Ddx4', 'Stra8', 'Sycp3',"Xist")
DotPlot(sce, group.by = 'seurat_clusters',
        features = unique(genes_to_check)) + RotatedAxis()

p1=DimPlot(sce, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
p2=DotPlot(sce, group.by = 'seurat_clusters',
           features = unique(genes_to_check)) + RotatedAxis()
library(patchwork)
p1+p2

my_color=c("#1D77B2","#ACC8E6", "#FD7F0D", "#FEBA73", "#2CA029", 
           "#97E083", "#D7262B", "#FD9797", "#9667BB",
           "#C4B1D3", "#8B564C", "#7859A3", "#E376BE",
           "#EEB9D1", "#929597", "#CCCECE", "#A7A123","#DB7093")

DimPlot(sce, reduction = "umap", cols = c("#2CA029","#9667BB","#FD7F0D"), group.by = "time")
DimPlot(sce, reduction = "umap",cols=c( "#2CA029", "#97E083", "#D7262B", "#FD9797", "#9667BB","#C4B1D3"))

FeaturePlot(sce, features = c("Stra8"), min.cutoff = "q9") +
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20)) ## set axis font size
head(sce@meta.data)
VlnPlot(sce, features = c("Cenpf"),  group.by = "time", cols = c("#2CA029","#9667BB","#FD7F0D"),
        pt.size = 0, combine = FALSE)

library(ggpubr)

VlnPlot(sce, features = c("Cenpf"),  group.by = "time", cols = c("#2CA029","#9667BB","#FD7F0D"),
        pt.size = 0, combine = FALSE)

my_comparisons <- list(c("GR11","GR12"), c("GR12", "GR13"),c("GR11", "GR13"))
my_comparisons <- list(c("GR12", "GR13"))

VlnPlot(sce, features = c("Stra8"),
        pt.size = 0, 
        group.by = "time", 
        cols = c("#2CA029","#9667BB","#FD7F0D"),
       y.max = 4 # add the y-axis maximum value - otherwise p-value hidden
) + stat_compare_means(comparisons = my_comparisons, label.y=c(3.5, 4, 4.5))

VlnPlot(sce, features = c("Stk31"),
        pt.size = 0, 
        group.by = "time", 
        cols = c("#2CA029","#9667BB","#FD7F0D"),
        y.max = 4 # add the y-axis maximum value - otherwise p-value hidden
) + stat_compare_means(comparisons = my_comparisons, label.y=c(2.5, 3, 3.5))



save.image(file="E:\\Lab of Germ Cell Bio\\10xGenomics_meiosis_V3\\scRNA_meiosis\\@@@2021.5.12_GR11_vs_GR12\\2021.5.26_Germ.RData")
########################### ?????? 2020.9.3 Extract germ cells ?????? ##########################################
########################### ?????? 2020.9.3 Extract germ cells ?????? ##########################################



