# Sergio Alías, 20231115
# Last modified 20231115

# Script for replicating mouse retina paper (integration of all samples)

library(Seurat)
library(scCustomize)
library(patchwork)
library(harmony)

dataset_loc <- "/mnt2/fscratch/users/bio_267_uma/sergioalias/NGS_projects/MOUSE/count_results"
ids <- c("WT_1",
         "WT_2",
         "rd10_1",
         "rd10_2")

seu.list <- sapply(ids, function(i){ # Loading
  d10x <- Read10X(file.path(dataset_loc,i,"cellranger_0000",i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  seu <- CreateSeuratObject(counts = d10x, project = "mouse", min.cells = 1, min.features = 1)
  seu@meta.data$samplename <- c(rep(i, nrow(seu@meta.data)))
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 20)
  seu
})

seu.list$WT_1@meta.data$genotype <- c(rep("WT", nrow(seu.list$WT_1@meta.data)))
seu.list$WT_2@meta.data$genotype <- c(rep("WT", nrow(seu.list$WT_2@meta.data)))
seu.list$rd10_1@meta.data$genotype <- c(rep("mutant", nrow(seu.list$rd10_1@meta.data)))
seu.list$rd10_2@meta.data$genotype <- c(rep("mutant", nrow(seu.list$rd10_2@meta.data)))

# Preprocessing trying to follow mouse paper

seu <- Merge_Seurat_List(list_seurat = seu.list)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunHarmony(seu, "samplename", plot_convergence = FALSE)

seu <- RunUMAP(seu, dims = 1:10, reduction = "harmony")
seu <- FindNeighbors(seu, dims = 1:10, reduction = "harmony")
seu <- FindClusters(seu, resolution = 0.6)

seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# We save the Seurat object and marker genes

setwd("outs")

saveRDS(seu, file = "seu.mouse.RDS")
saveRDS(seu.markers, file = "seu.markers.mouse.RDS")