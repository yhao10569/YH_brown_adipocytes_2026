
####single sample##########
library(Seurat)
library(tidyverse)
library(scutilsR)

##d4
SIM_Matrix=Read10X(data.dir = 'data/10X_matrix/SIM_D4')
SIM=CreateSeuratObject(counts = SIM_Matrix,project = 'SIM',min.cells = 3,min.features = 200)
Idents(SIM)=factor('SIM')
SIM@meta.data$orig.ident=factor('SIM')
##QC
SIM[["percent.mt"]] = PercentageFeatureSet(SIM, pattern = "^MT-")
SIM[['percent.RPS']]=PercentageFeatureSet(SIM,pattern = '^RPS')
SIM[['percent.RPL']]=PercentageFeatureSet(SIM,pattern = '^RPL')
SIM$QC <- "high"
SIM$QC[SIM$nFeature_RNA<1000|SIM$percent.mt > 15] <- "low"
SIM <- subset(SIM, QC == "high")

SIM <- NormalizeData(SIM)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
SIM <- CellCycleScoring(SIM, s.features = s.genes, g2m.features = g2m.genes)
SIM <- scutilsR::MarkDoublets(SIM)
QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
SIM.list <- SplitObject(SIM, split.by = "orig.ident")
SIM.list <- lapply(SIM.list, QuickCluster)
clusters <- lapply(SIM.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
SIM$quick_clusters <- clusters[rownames(SIM@meta.data)]
rm(SIM.list)
gc()
SIM <- scutilsR::RemoveAmbientRNAs(SIM, split.by = "orig.ident", cluster.name = "quick_clusters")
qs::qsave(SIM, file.path('output/01_output/QC', "SIM_QC2.qs"))

SIM<- subset(E,subset =nFeature_RNA >1000 & nFeature_RNA <5000 & percent.mt>3 & percent.mt<15 & nCount_RNA <15000 )
SIM=subset(x=SIM, DF.classifications == 'Singlet')
table(SIM$orig.ident,SIM$DF.classifications)
SIM=NormalizeData(SIM) %>% 
  FindVariableFeatures() %>% 
  ScaleData(., vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA()

resolutions <- c(0.1,0.2,0.3)
SIM<- FindNeighbors(SIM, reduction = "pca", dims = 1:20, k.param = 20)
SIM <- FindClusters(SIM, resolution = resolutions)
SIM=RunUMAP(SIM,reduction = 'pca',dims = 1:30,seed.use = 32)

##annotation
Idents(SIM) <- SIM$RNA_snn_res.0.1

SIM<-RenameIdents(SIM,
                  '0'='Somite',
                  '1'='End.Pro.',
                  '2'='VE')
SIM$annotation <- Idents(SIM)
qs::qsave(SIM, file.path('output/01_output/01_singel_sample/', "SIM_downstream.qs"))

##d8
TIFL_Matrix=Read10X(data.dir = 'data/10X_matrix/TIFL_D2')
TIFL=CreateSeuratObject(counts = TIFL_Matrix,project = 'SIM',min.cells = 3,min.features = 200)
Idents(TIFL)=factor('TIFL')
TIFL@meta.data$orig.ident=factor('TIFL')
####QC
TIFL[["percent.mt"]] = PercentageFeatureSet(TIFL, pattern = "^MT-")
TIFL[['percent.RPS']]=PercentageFeatureSet(TIFL,pattern = '^RPS')
TIFL[['percent.RPL']]=PercentageFeatureSet(TIFL,pattern = '^RPL')
TIFL$QC <- "high"
TIFL$QC[TIFL$nFeature_RNA<1000 | TIFL$percent.mt > 15 | TIFL$percent.mt < 3] <- "low"
table(TIFL$orig.ident, TIFL$QC)
TIFL <- subset(TIFL, QC == "high")
TIFL <- NormalizeData(TIFL)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
TIFL <- CellCycleScoring(TIFL, s.features = s.genes, g2m.features = g2m.genes)
TIFL <- scutilsR::MarkDoublets(TIFL)
QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
TIFL.list <- SplitObject(TIFL, split.by = "orig.ident")
TIFL.list <- lapply(TIFL.list, QuickCluster)
clusters <- lapply(TIFL.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
TIFL$quick_clusters <- clusters[rownames(TIFL@meta.data)]
rm(TIFL.list)
gc()
TIFL <- scutilsR::RemoveAmbientRNAs(TIFL, split.by = "orig.ident", cluster.name = "quick_clusters")
qs::qsave(TIFL, file.path('output/01_output/QC', "TIFL_QC2.qs"))
TIFL<-subset(TIFL,subset =nFeature_RNA >1000 & nFeature_RNA <5000 & percent.mt>3 & percent.mt<15 & nCount_RNA <15000 )
TIFL <-subset(TIFL, DF.classifications  == 'Singlet')
TIFL=NormalizeData(TIFL) %>% 
  FindVariableFeatures() %>% 
  ScaleData(., vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA()

##cluster
resolutions <- c(0.1,0.2,0.3,0.5,0.8)
TIFL<- FindNeighbors(TIFL, reduction = "pca", dims = 1:20, k.param = 20)
TIFL <- FindClusters(TIFL, resolution = resolutions)
TIFL=RunUMAP(TIFL,reduction = 'pca',dims = 1:30,seed.use = 3)
#annotation
Idents(TIFL) <- TIFL$RNA_snn_res.0.5
Idents(TIFL)<-TIFL$annotation
TIFL <- RenameIdents(TIFL,
                     '0_5'='End.Pro.',
                     '1_2'='End.mural',
                     '3'='End.VE',
                     '4_6'='Somite Pro.',
                     '7'='VE')
TIFL$annotation <- Idents(TIFL)
qs::qsave(TIFL, file.path('output/01_output/01_singel_sample/', "TIFL_downstream.qs"))

##d12
Priming_Matrix=Read10X(data.dir = 'data/10X_matrix/Priming_D4')
Priming_D4=CreateSeuratObject(counts = Priming_Matrix,project = 'Priming_D4',min.cells = 3,min.features = 200)

Idents(Priming_D4)=factor('Priming_D4')
Priming_D4@meta.data$orig.ident=factor('priming_D4')
#QC
Priming_D4[["percent.mt"]] = PercentageFeatureSet(Priming_D4, pattern = "^MT-")
Priming_D4[['percent.RPS']]=PercentageFeatureSet(Priming_D4,pattern = '^RPS')
Priming_D4[['percent.RPL']]=PercentageFeatureSet(Priming_D4,pattern = '^RPL')
Priming_D4$QC <- "high"
Priming_D4$QC[Priming_D4$nFeature_RNA<1000 | Priming_D4$percent.mt > 15] <- "low"
Priming_D4 <- subset(Priming_D4, QC == "high")
Priming_D4 <- NormalizeData(Priming_D4)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
Priming_D4 <- CellCycleScoring(Priming_D4, s.features = s.genes, g2m.features = g2m.genes)
Priming_D4 <- scutilsR::MarkDoublets(Priming_D4)
QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
Priming_D4.list <- SplitObject(Priming_D4, split.by = "orig.ident")
Priming_D4.list <- lapply(Priming_D4.list, QuickCluster)
clusters <- lapply(Priming_D4.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
Priming_D4$quick_clusters <- clusters[rownames(Priming_D4@meta.data)]
rm(Priming_D4.list)
gc()
Priming_D4 <- scutilsR::RemoveAmbientRNAs(Priming_D4, split.by = "orig.ident", cluster.name = "quick_clusters")
Priming_D4 
qs::qsave(Priming_D4, file.path('output/01_output/QC', "Priming_QC2.qs"))
Priming_D4=NormalizeData(Priming_D4) %>% 
  FindVariableFeatures() %>% 
  ScaleData(., vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA()
Priming_D4=RunUMAP(Priming_D4,reduction = 'pca',dims = 1:30,seed.use = 32)
qs::qsave(Priming_D4, file.path('output/01_output/01_singel_sample/', "Priming_downstream.qs"))

##d16
BA_Matrix=Read10X(data.dir = 'data/10X_matrix/BA_D4')
BA_D4=CreateSeuratObject(counts = BA_Matrix,project = 'BA_D4',min.cells = 3,min.features = 200)

Idents(BA_D4)=factor('BA_D4')
BA_D4@meta.data$orig.ident=factor('BA_D4')
####QC
BA_D4[["percent.mt"]] = PercentageFeatureSet(BA_D4, pattern = "^MT-")
BA_D4[['percent.RPS']]=PercentageFeatureSet(BA_D4,pattern = '^RPS')
BA_D4[['percent.RPL']]=PercentageFeatureSet(BA_D4,pattern = '^RPL')
BA_D4$QC <- "high"
BA_D4$QC[BA_D4$nFeature_RNA<1000 | BA_D4$percent.mt > 15] <- "low"
BA_D4 <- subset(BA_D4, QC == "high")
BA_D4 <- NormalizeData(BA_D4)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
BA_D4 <- CellCycleScoring(BA_D4, s.features = s.genes, g2m.features = g2m.genes)
BA_D4 <- scutilsR::MarkDoublets(BA_D4)
QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
BA_D4.list <- SplitObject(BA_D4, split.by = "orig.ident")
BA_D4.list <- lapply(BA_D4.list, QuickCluster)
clusters <- lapply(BA_D4.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
BA_D4$quick_clusters <- clusters[rownames(BA_D4@meta.data)]
rm(BA_D4.list)
gc()
BA_D4 <- scutilsR::RemoveAmbientRNAs(BA_D4, split.by = "orig.ident", cluster.name = "quick_clusters")
qs::qsave(BA_D4, file.path('output/01_output/QC', "BA4_QC2.qs"))
BA_D4=NormalizeData(BA_D4) %>% 
  FindVariableFeatures() %>% 
  ScaleData(., vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA()
BA_D4=RunUMAP(BA_D4,reduction = 'pca',dims = 1:30,seed.use = 32)
qs::qsave(BA_D4, file.path('output/01_output/01_singel_sample/BA4_downstream.qs'))


##d20
BA_Matrix=Read10X(data.dir = 'data/10X_matrix/BA_D8')
BA_D8=CreateSeuratObject(counts = BA_Matrix,project = 'BA_D8',min.cells = 3,min.features = 200)

Idents(BA_D8)=factor('BA_D8')
BA_D8@meta.data$orig.ident=factor('BA_D8')
####QC
BA_D8[["percent.mt"]] = PercentageFeatureSet(BA_D8, pattern = "^MT-")
BA_D8[['percent.RPS']]=PercentageFeatureSet(BA_D8,pattern = '^RPS')
BA_D8[['percent.RPL']]=PercentageFeatureSet(BA_D8,pattern = '^RPL')
BA_D8$QC <- "high"
BA_D8$QC[BA_D8$nFeature_RNA<1000 | BA_D8$percent.mt > 15] <- "low"
BA_D8 <- subset(BA_D8, QC == "high")
BA_D8 <- NormalizeData(BA_D8)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
BA_D8<- CellCycleScoring(BA_D8, s.features = s.genes, g2m.features = g2m.genes)
BA_D8 <- scutilsR::MarkDoublets(BA_D8)

QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
BA_D8.list <- SplitObject(BA_D8, split.by = "orig.ident")
BA_D8.list <- lapply(BA_D8.list, QuickCluster)
clusters <- lapply(BA_D8.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
BA_D8$quick_clusters <- clusters[rownames(BA_D8@meta.data)]
rm(BA_D8.list)
gc()

BA_D8 <- scutilsR::RemoveAmbientRNAs(BA_D8, split.by = "orig.ident", cluster.name = "quick_clusters")
qs::qsave(BA_D8, file.path('output/01_output/QC', "BA8_QC2.qs"))
BA_D8=NormalizeData(BA_D8) %>% 
  FindVariableFeatures() %>% 
  ScaleData(., vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA()
BA_D8=RunUMAP(BA_D8,reduction = 'pca',dims = 1:30,seed.use = 32)
qs::qsave(BA_D8, file.path('output/01_output/01_singel_sample/BA8_downstream.qs'))

##d24
BA_Matrix=Read10X(data.dir = 'data/10X_matrix/BA_D12')
BA_D12=CreateSeuratObject(counts = BA_Matrix,project = 'BA_D12',min.cells = 3,min.features = 200)

Idents(BA_D12)=factor('BA_D12')
BA_D12@meta.data$orig.ident=factor('BA_D12')
####QC
BA_D12[["percent.mt"]] = PercentageFeatureSet(BA_D12, pattern = "^MT-")
BA_D12[['percent.RPS']]=PercentageFeatureSet(BA_D12,pattern = '^RPS')
BA_D12[['percent.RPL']]=PercentageFeatureSet(BA_D12,pattern = '^RPL')
BA_D12$QC <- "high"
BA_D12$QC[BA_D12$nFeature_RNA<1000 | BA_D12$percent.mt > 15] <- "low"
BA_D12 <- subset(BA_D12, QC == "high")
BA_D12 <- NormalizeData(BA_D12)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
BA_D12 <- CellCycleScoring(BA_D12, s.features = s.genes, g2m.features = g2m.genes)
BA_D12 <- scutilsR::MarkDoublets(BA_D12)
QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
BA_D12.list <- SplitObject(BA_D12, split.by = "orig.ident")
BA_D12.list <- lapply(BA_D12.list, QuickCluster)
clusters <- lapply(BA_D12.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
BA_D12$quick_clusters <- clusters[rownames(BA_D12@meta.data)]
rm(BA_D12.list)
gc()
BA_D12 <- scutilsR::RemoveAmbientRNAs(BA_D12, split.by = "orig.ident", cluster.name = "quick_clusters")
qs::qsave(BA_D12, file.path('output/01_output/QC', "BA12_QC2.qs"))
BA_D12=NormalizeData(BA_D12) %>% 
  FindVariableFeatures() %>% 
  ScaleData(., vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA()
BA_D12=RunUMAP(BA_D12,reduction = 'pca',dims = 1:30,seed.use = 32)
qs::qsave(BA_D12, file.path('output/01_output/01_singel_sample/BA12_downstream.qs'))


#######integration###########
library(Seurat)
library(harmony)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(glue)
library(bbknnR)
data1 <- qs::qread(file.path('output/01_output/01_singel_sample/SIM_downstream.qs'))
data2 <- qs::qread(file.path('output/01_output/01_singel_sample/TIFL_downstream.qs'))
data3 <- qs::qread(file.path('output/01_output/01_singel_sample/Priming_downstream.qs'))
data4 <- qs::qread(file.path('output/01_output/01_singel_sample/BA4_downstream.qs'))
data5 <- qs::qread(file.path('output/01_output/01_singel_sample/BA8_downstream.qs'))
data6 <- qs::qread(file.path('output/01_output/01_singel_sample/BA12_downstream.qs'))

data <- list(data1, data2,data3,data4,data5,data6)
E <- merge(x=data[[1]],y=data[-1])
rm(list = ls(pattern="data.*"))
E@assays$decontX = NULL
Idents(E) <- E$orig.ident
VlnPlot(E_BA, features = c("percent.mt",'percent.RPS','percent.RPL'), ncol = 3,raster=FALSE)
E <- subset(E,subset =nFeature_RNA >1000 & nFeature_RNA <5000 & percent.mt>3 & percent.mt<15 & nCount_RNA <15000 )
E<- subset(E,subset = DF.classifications == 'Singlet')

E <- NormalizeData(E)
E <- FindVariableFeatures(E, selection.method = "vst", nfeatures = 2000)
E <- ScaleData(E, vars.to.regress = c('S.Score','G2M.Score','percent.RPL','percent.RPS'))
E <- RunPCA(E, features = VariableFeatures(object = E),reduction.name = "pca")
E <- RunUMAP(E, reduction = "pca", dims = 1:30,reduction.name = "umap")
E <- RunHarmony(E,reduction.use = "pca",group.by.vars = 'orig.ident',reduction.save = "harmony",dims.use=1:50)
E <- RunUMAP(E, reduction = "harmony", dims = 1:30,reduction.name = "umap.harmony")

DimPlot(E, reduction = "umap.harmony",group.by = "orig.ident",raster=FALSE)
##using bbknnR
E=RunBBKNN(E,reduction = 'pca',run_UMAP = TRUE, UMAP_name="umap.bbknn",run_TSNE= FALSE,
           batch_key = 'orig.ident',seed = 23)
qs::qsave(BA,file = 'output/01_output/02_two_pathway_integration/Integration/01_1_E_cc_rpl.qs')

##annotation
E<-qs::qread(file = 'output/01_output/02_two_pathway_integration/Integration/01_1_E_cc_rpl.qs')
reso<-c(0.5,0.8,1,1.2,1.5,2)
E<- FindClusters(E,resolution =reso,graph.name = 'bbknn' )
Idents(E)<-E$bbknn_res.1.2

E<-FindSubCluster(E,resolution = 0.6,cluster = 5,subcluster.name = 'c5_res6',graph.name = 'bbknn')
E$c5_res6<-as.factor(E$c5_res6)
Idents(E)<-E$c5_res6
E<-RenameIdents(E,
                '0'='Somite',
                '1'='Adipogenic fibroblast',
                '2'='WNT primed precusors',
                '3'='Brown progenitors',
                '4'='Adipogenic fibroblast',
                '5_0'='SMC',
                '5_1'='Endotome',
                '5_2'='SMC',
                '5_3'='SMC',
                '5_4'='Endotome',
                '5_5'='Endotome',
                '6'='Endotome',
                '7'='SKM',
                '8'= 'Endotome',
                '9'='BAs',
                '10'='SKM',
                '11'='Neural progenitors',
                '12'='Endothelium cells')
qs::qsave(E,file = 'output/01_output/02_two_pathway_integration/Integration/01_1_E_cc_rpl.qs')

##BA subset
E<-qs::qread(file = 'output/01_output/02_two_pathway_integration/Integration/01_1_E_cc_rpl.qs')
##BA subset
E_BA<-subset(E,celltype %in% c('Somite','Endotome','WNT primed precusors','Adipogenic fibroblast','Brown progenitors','BAs'))
E_BA <- NormalizeData(E_BA)
E_BA <- FindVariableFeatures(E_BA, selection.method = "vst", nfeatures = 2000)
E_BA <- ScaleData(E_BA, vars.to.regress = c('S.Score','G2M.Score','percent.RPL','percent.RPS'))
E_BA <- RunPCA(E_BA, features = VariableFeatures(object = E_BA),reduction.name = "pca")
E_BA <- RunHarmony(E_BA,reduction.use = "pca",group.by.vars = 'orig.ident',reduction.save = "harmony",dims.use=1:50)
E_BA <- RunUMAP(E_BA, reduction = "pca", dims = 1:30,reduction.name = "umap")
E_BA <- RunUMAP(E_BA, reduction = "harmony", dims = 1:30,reduction.name = "umap.harmony")

DimPlot(E_BA, reduction = "umap.harmony",group.by = "celltype",raster=FALSE)
##using bbknnR
E_BA=RunBBKNN(E_BA,reduction = 'pca',run_UMAP = TRUE, UMAP_name="umap.bbknn",run_TSNE= FALSE,
              batch_key = 'orig.ident')
DimPlot(E_BA, reduction = "umap.bbknn",group.by = "celltype",raster=FALSE)

qs::qsave(E_BA,file = 'output/01_output/02_two_pathway_integration/Integration/01_5_EBF2_BA.qs')

#fdg
E_BA<-qs::qread(file = 'output/01_output/02_two_pathway_integration/Integration/01_5_EBF2_BA.qs')

RunFDG <- function(object, reduction, dims = NULL, reduction.name = "fr",
                   reduction.key = "FDG_", ...) {
  if (is.null(dims)) {
    dims <- 1:ncol(Embeddings(object, reduction = reduction))
  }
  object <- FindNeighbors(object, reduction = reduction, dims = dims, ...)
  g <- as(object[["RNA_snn"]], "dgCMatrix")
  g <- igraph::graph_from_adjacency_matrix(adjmatrix = g,
                                           mode = "undirected",
                                           weighted = TRUE,
                                           add.colnames = TRUE)
  fr.layout <- igraph::layout_with_fr(g)
  rownames(fr.layout) <- colnames(object)
  ## store cell embeddings to objectrat
  object[[reduction.name]] <- CreateDimReducObject(fr.layout, key = reduction.key)
  return(object)
}

set.seed(32)
cellstokeep <- sample(rownames(E_BA@meta.data),20000)
E_BA2 <- subset(E_BA, cells = cellstokeep)

E_BA2<-RunHarmony(E_BA2,,reduction.use = "pca",group.by.vars = 'orig.ident',reduction.save = "harmony",dims.use = 1:30,
                  theta = 1.5, lambda = 0.5)
set.seed(18)
E_BA2<- RunFDG(E_BA2, reduction = "harmony", dims = 1:30)
DimPlot(E_BA2, reduction = "fr", group.by = "celltype", label = T)
qs::qsave(E_BA2,file = 'output/01_output/02_two_pathway_integration/Integration/01_5_EBF2_BA_slim.qs')
E_BA3<-E_BA2[rownames(E_BA2[['RNA']]@scale.data),]
sceasy::convertFormat(E_BA3, from = "seurat", to = "anndata", main_layer = "scale.data", 
                      outFile = file.path("output/01_output/02_two_pathway_integration/Integration/01_5_EBF2_BA_scaled.h5ad"))


root.cells <- CellSelector(FeaturePlot(E_BA2, reduction = 'fr', features = 'PAX3'))
"SIM_D4-GCCGTTATT_CGCTGTAAG_GGTGGTATC"



####visualization
library(Seurat)
library(scCustomize)
library(Nebulosa)
source('Resource/color_map.R')


##annotation
E<-qs::qread(file = 'output/01_output/02_two_pathway_integration/Integration/01_1_E_cc_rpl.qs')
Idents(E)<-E$celltype
scCustomize::DimPlot_scCustom(E, figure_plot = TRUE,pt.size = 0.2,reduction = 'umap.bbknn',raster = F,colors_use = EBF2_path_1,label = F,group.by = 'celltype') 
ggsave(filename = 'output/01_output/02_two_pathway_integration/Figures/01_E_annotation.pdf',width = 9.79,height = 6.43)
ggsave(filename = 'output/01_output/02_two_pathway_integration/Figures/01_E_annotation.tiff',width = 9.79,height = 6.43)
##percentage
E$sample<-factor(E$group)
E$sample<-factor(E$sample,labels = c('Day6','Day8','Day12','Day16','Day20','Day24'))
data<-as.data.frame(table(E$sample,E$celltype))
colnames(data)<-c('Sample','CellType','Freq')
library(dplyr)
df <- data %>% 
  group_by(Sample) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = Freq/Total) %>% 
  as.data.frame()
df$CellType  <- factor(df$CellType,levels = unique(df$CellType))
library(ggplot2)
ggplot(df, aes(x = Sample, y = Percent, fill = CellType)) +
  geom_bar(position = "fill", stat="identity", color = 'white', alpha = 1, width = 0.95) +
  scale_fill_manual(values = EBF2_path_1) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14, face = "bold", angle = 45, vjust = 1,hjust = 1),  # 横轴字体
        axis.text.y = element_text(size = 14, face = "bold"),                          # 纵轴字体
        axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 12))
ggsave(filename = 'output/01_output/02_two_pathway_integration/Figures/01_E_celltype_proportion.pdf',width = 7.36,height = 7.91)
ggsave(filename = 'output/01_output/02_two_pathway_integration/Figures/01_E_celltype_proportion.tiff',width = 7.36,height = 7.91)
##
library(RColorBrewer)
scCustomize::DimPlot_scCustom(E, figure_plot = TRUE,pt.size = 0.2,reduction = 'umap.bbknn',raster = F,colors_use = RColorBrewer::brewer.pal(6, "Set2"),label = F,group.by = 'sample') 
ggsave(filename ='output/01_output/02_two_pathway_integration/Figures/01_E_sample.pdf',width = 9.79,height = 6.43)
ggsave(filename ='output/01_output/02_two_pathway_integration/Figures/01_E_sample.tiff',width = 9.79,height = 6.43)
scCustomize::DimPlot_scCustom(E,figure_plot = TRUE,pt.size = 0.2,reduction = 'umap.bbknn',raster = F,colors_use = EBF2_path_1,label = F,split.by = 'sample') 

DimPlot(E,split.by = 'sample',reduction = 'umap.bbknn',pt.size = 0.2,label = F,group.by = 'celltype',cols = EBF2_path_1)
ggsave(filename ='output/01_output/02_two_pathway_integration/Figures/01_E_sample_split.pdf',width = 18,height = 4)
ggsave(filename ='output/01_output/02_two_pathway_integration/Figures/01_E_sample_split.tiff',width = 18,height = 4)

####heatmap
markers <- c('HES7','PAX3','SIX1',#somite
             'FNDC1','FOXS1','DLX1','MYH11','GDF7','CXCL12','SYNPO',  #Endotome
             'TGFBI','PENK','ENG','SOX9','EBF2',#BA progenitor
             'THY1','HIC1','PDGFRA','OSR1','SCARA5','WNT11','SCARA3',#adipogenic fibroblast
             'CEBPD','PRRX2','FGFBP2','PPARG',#brown pre
             'CEBPA','FABP4','UCP1','PLIN1','PLIN5','ADIPOQ','KLB','TRARG1','SGK2', #brown
             'CNN1','ACTA2','TAGLN', #smooth muscle
             'MYH8','MYOD1','MYOG',  #skeleton muscle
             'SOX17','CLDN5','KDR', #endothelium
             'SOX2','HES5','WNT4') #neural 

library(dplyr)
library(Seurat) 
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(scales)


P <- DotPlot(E,features = markers)
data <- P$data
data$id<-factor(data$id)
data$id<-factor(data$id,levels = rev(levels(E$celltype)))
data$color<-EBF2_path_1[as.character(data$id)]

library(ggnewscale)
ggplot(data) +
  geom_tile(aes(x = features.plot, y = id, fill = avg.exp.scaled),
            color = NA, width = 1, height = 1) +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#FFF5EB", "#FDBB84", "#E34A33", "#B30000"))(100),
    name = "avg.exp.scaled"
  ) +
  geom_point(
    aes(x = features.plot, y = id, size = squish(pct.exp, c(0, 50))),
    shape = 21, color = "black", stroke = 0.5
  ) +scale_size_continuous(range = c(0,6), limits = c(5,50)) +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = -0.5, y = id, fill = color),
            width = 1.1, height = 1.05, color = NA) +
  scale_fill_identity(name = "EBF2 Path") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.border      = element_blank(),
    axis.title        = element_blank(),
    axis.text.x       = element_text(size = 15, angle = 90,face = 'bold',
                                     hjust = 1, vjust = .5),
    axis.text.y       = element_blank(),
    legend.text       = element_text(size = 8),
    legend.title      = element_text(size = 8),
    legend.key        = element_rect(fill = "white", colour = NA),
    legend.key.size   = unit(2.5, "mm"),
    legend.spacing    = unit(-1, "mm")
  ) +
  scale_x_discrete(expand = c(0, 0)) +  
  scale_y_discrete(expand = c(0, 0))   
ggsave(filename = 'output/01_output/02_two_pathway_integration/Figures/01_E_markers_heatmap_new.pdf',width = 10.6,height = 5)


###featureplot
marker_genes<-c('PAX3','EBF2','MYF5','PAX7','EN1','GATA6','PENK','HIC1','PPARG','PRDM16','ADIPOQ','CEBPD','UCP1')
scCustomize::FeaturePlot_scCustom(E,features = c('PAX7','MYOD1','MYOG'),order = T,colors_use = c("#FFF5EB", "#FDBB84", "#E34A33", "#B30000"),reduction = 'umap.bbknn')+theme_void() + 
  theme(legend.position = "none")


marker_genes<-c('PAX3','SOX2','PAX7','EN1',
                'PENK','GDF7','TAGLN','ERG',
                'MYH11','FOXS1','DLX1','EBF2',
                'HIC1','ZIC1','WNT11','MCAM',
                'PPARG','FABP4','NR4A1','CD34',
                'CEBPA','PLIN1','ADIPOQ','UCP1','CEBPD')
for (marker in marker_genes) {
  filename <- paste0("output/01_output/02_two_pathway_integration/Figures/01_E_", marker, ".tiff")
  tiff(file = filename, width = 258, height = 233)
  print(scCustomize::FeaturePlot_scCustom(E,features = marker,order = T,colors_use = c("#FFF5EB", "#FDBB84", "#E34A33", "#B30000"),reduction = 'umap.bbknn')+theme_void() + 
          theme(legend.position = "none"))
  dev.off()
}


####subseted BA
E_BA<-qs::qread(file = 'output/01_output/02_two_pathway_integration/Integration/01_5_EBF2_BA_slim.qs')
E_BA$celltype1<-as.factor(E_BA$celltype)
E_BA$celltype1<-factor(E_BA$celltype,labels = c('1:Somite','2:Endotome','3:WNT primed precusors','4:Adipogenic fibroblast','5:Brown progenitors',
                                                '6:BAs'))
E_BA$celltype2<-factor(E_BA$celltype2,labels = c('1','2','3','4','5','6'))
scCustomize::DimPlot_scCustom(E_BA, figure_plot = TRUE,pt.size = 0.6,reduction = 'fr',raster = F,colors_use = EBF2_path_1,label = F) 
ggsave(filename = 'output/01_output/02_two_pathway_integration/Figures/01_E_BA_FDG.pdf',width = 10.57,height = 6.36)
ggsave(filename = 'output/01_output/02_two_pathway_integration/Figures/01_E_BA_FDG.tiff',width = 10.57,height = 6.36)

scCustomize::DimPlot_scCustom(E_BA, figure_plot = TRUE,pt.size = 0.6,reduction = 'fr',raster = F,colors_use = EBF2_path_2,label = F,group.by = 'celltype1',label.size = 8) 
ggsave(filename = 'output/01_output/02_two_pathway_integration/Figures/01_E_BA_FDG_celltype1.pdf',width = 5.34,height = 4.36)
ggsave(filename = 'output/01_output/02_two_pathway_integration/Figures/01_E_BA_FDG_celltype1.tiff',width =  5.34,height = 4.36)

qs::qsave(E_BA,file = 'output/01_output/02_two_pathway_integration/Integration/01_5_EBF2_BA_slim.qs')


E_BA<-qs::qread(file = 'output/01_output/02_two_pathway_integration/Integration/01_5_EBF2_BA_slim.qs')

marker_genes<-c(
  "PAX3", "EBF2", "TGFBR2", "TGFBI", "CDH15", "PENK", "PLVAP", "MYH11",
  "DLX1", "FOXS1", "MYF5", "PAX7", "EN1", "GATA6", "HIC1", "OSR1", 
  "OSR2", "ZIC1", "PDGFRA", "PPARG", "ADIPOQ", "UCP1")
plot_density(E_BA,features = 'CD34',reduction = 'fr')+theme_void()
for (marker in marker_genes) {
  filename <- paste0("output/01_output/02_two_pathway_integration/Figures/01_E_slim", marker, ".tiff")
  tiff(file = filename, width = 258, height = 233)
  print(scCustomize::FeaturePlot_scCustom(E_BA,features = marker,order = T,colors_use = c("#FFF5EB", "#FDBB84", "#E34A33", "#B30000"),reduction = 'fr',pt.size = 0.5)+theme_void() + 
          theme(legend.position = "none"))
  dev.off()
}