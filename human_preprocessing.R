# Sergio Alías, 20231020
# Last modified 20231107

# Script for replicating human retina paper (integration of PCW samples)

library(Seurat)
library(DoubletFinder)
library(scCustomize)
library(patchwork)
library(harmony)

dataset_loc <- "/mnt2/fscratch/users/bio_267_uma/sergioalias/NGS_projects/MACULAR/count_results"
ids <- c("21PCW", # Figure 1 covers only 12 to 21 PCW samples
         "20PCW",
         "12PCW",
         "16PCW")

seu.list <- sapply(ids, function(i){ # Loading
  d10x <- Read10X(file.path(dataset_loc,i,"cellranger_0000",i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  seu <- CreateSeuratObject(counts = d10x, project = "macular", min.cells = 1, min.features = 1)
  seu@meta.data$samplename <- c(rep(i, nrow(seu@meta.data)))
  seu
})

seu.doublets <- lapply(seu.list, function(seu){ # Preprocessing following the paper (see the Obsidian for reference)
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 20)

    # It looks like they used standard Seurat preprocessing before DoubletFinder, so:
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    seu <- ScaleData(seu)
    seu <- RunPCA(seu)
    seu <- RunUMAP(seu, dims = 1:10)

    # Now we run DoubletFinder

    nExp <- round(ncol(seu) * 0.04)  # Expect 4% doublets
    seu <- doubletFinder_v3(seu, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
    DF.name <- colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]
    return(seu@meta.data[, DF.name])

    # We don't filter doublet cells, we store those cells to remove them later
})

seu.list <- sapply(ids, function(i){ # Loading
  seu <- seu.list[[i]]
  is.doublet <- seu.doublets[[i]]
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 20)
  seu <- seu[, is.doublet == "Singlet"]
})


# DoubletFinder recommends to perform preprocessing separately in each sample, but Harmony needs a single merged Seurat object. So:

seu <- Merge_Seurat_List(list_seurat = seu.list)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunHarmony(seu, "samplename", plot_convergence = FALSE)

# Now clustering and non linear dimensional reduction using Harmony embeddings

seu <- FindNeighbors(seu, dims = 1:10, reduction = "harmony")
seu <- FindClusters(seu, resolution = 0.6)
seu <- RunUMAP(seu, dims = 1:10, reduction = "harmony")

# And marker gene selection

seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Some metadata useful for the DE genes step:

seu$stage <- ifelse(seu$samplename %in% c("12PCW", "16PCW"), "early", "late")

# We save the Seurat object and marker genes

setwd("outs")

saveRDS(seu.markers, file = "seu.markers.human.RDS")
saveRDS(seu, file = "seu.human.RDS")