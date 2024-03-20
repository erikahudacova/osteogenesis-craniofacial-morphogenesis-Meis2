#scRNA seq analysis of CNC-derived MSCs preprocessing E12.5-E13.5 by ERIKA H. 
r libraries ,message=FALSE, warning=FALSE 
library(Seurat)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyverse)
library(dplyr)
library(STutility)
library(sctransform)
library(ggraph)
library(clustree)
library(scater)
library(writexl)

#open subsetted mesenchymal cluster
seurat_mesenchymal <- readRDS("seurat_mesenchymal_clusters.rds")

#subset mitochondrial content and removal of debris cluster
seurat_mesenchymal_mt <- subset(seurat_mesenchymal, subset = percent.mt < 10)
seurat_mesenchymal_mt  <- SetIdent(seurat_mesenchymal_mt, value = "nicknames")
seurat_filtered <- subset(seurat_mesenchymal_mt , idents= "Neuro", invert = TRUE)

saveRDS(seurat_filtered,"seurat_mesenchymal_filtered.rds")

#normalization and regression of cell cycle

DefaultAssay(seurat_filtered) <- "RNA"
# Make a list of all samples
experiment.list <- SplitObject(seurat_filtered, split.by = "Sample")

# run normalization for each sample

for (i in 1:length(experiment.list)) {
  experiment.list[[i]] <- SCTransform(experiment.list[[i]], assay = "RNA", vst.flavor="v2", verbose = FALSE, return.only.var.genes = TRUE, vars.to.regress = c("G2M.Score","S.Score"))}

saveRDS(experiment.list,"seurat_mesenchymal_SCTransfromregressed.rds")

#integration across 4 conditions
experiment.list  <- readRDS("seurat_mesenchymal_SCTransfromregressed.rds")

options(future.globals.maxSize = 100 * 1024^3)

experiment.features <- SelectIntegrationFeatures(object.list = experiment.list, nfeatures = 4000)

experiment.features <- experiment.features[substr(experiment.features,1,3) != "mt-" &
                                             substr(experiment.features,1,3) != "Rpl" &
                                             substr(experiment.features,1,3) != "Rps"]

featuresforint <- experiment.features[!(experiment.features %in% cellcyclevector)] #removing cell cycle genes from exp features by subsetting vector 

experiment <- PrepSCTIntegration(object.list = experiment.list, anchor.features = experiment.features, 
                                 verbose = TRUE)

experiment.anchors <- FindIntegrationAnchors(object.list = experiment, normalization.method = "SCT", 
                                             anchor.features = experiment.features, verbose = TRUE, dims = 1:30)
experiment.integrated <- IntegrateData(anchorset = experiment.anchors, normalization.method = "SCT", 
                                       verbose = TRUE, dims = 1:30)

saveRDS(experiment.integrated,"seurat_mesenchymal_integratedregressed.rds")


experiment.integrated  <- readRDS("seurat_mesenchymal_integratedregressed.rds")

DefaultAssay(experiment.integrated) <- "integrated"
all.genes <- row.names(experiment.integrated)

# use Elbow Plot to decide some ‘optimal” number of PCs for Neighbors finding and clustering
seuratPCA <- RunPCA(experiment.integrated,
                    
                    assay = "SCT",
                    
                    npcs = 80, verbose = FALSE,
                    
                    reduction.name = "pca", features = all.genes)

ElbowPlot(seuratPCA, ndims = 100)

seurat_integrated_mesenchymal <- RunPCA(seurat_integrated_mesenchymal, npcs = 50, verbose = FALSE, features = all.genes)

experiment <- RunPCA(experiment.integrated, npcs = 50, verbose = FALSE, features = all.genes)

experiment<- RunUMAP(experiment, reduction = "pca", verbose = TRUE,features = all.genes)
experiment<- FindNeighbors(experiment, reduction = "pca", dims = 1:41, verbose = TRUE)

#clustering resolutions, bring you data after finding Neighbors and pick your resolution

experiment<- FindClusters(experiment, resolution = 0.6, verbose = FALSE,)

#normalize and scale rna
DefaultAssay(experiment) <- "RNA"
experiment <- NormalizeData(experiment, verbose = TRUE)
experiment <- ScaleData(experiment, verbose = TRUE)

#DEG marker analysis
experiment <- SetIdent(experiment,value = "integrated_snn_res.0.5")
DefaultAssay(experiment) <- "SCT"
experiment <- PrepSCTFindMarkers(experiment)
markers0.5<- FindAllMarkers(experiment,
                            min.pct = 0.25,
                            logfc.threshold = 0.2,
                            only.pos = FALSE,
                            verbose = FALSE,
                            assay = "SCT",
                            test.use = "wilcox")

write.xlsx(markers0.6,"marker0.6sct.xlsx", rowNames = TRUE)

saveRDS(experiment,file = "seurat_mesenchymal_final.rds")


