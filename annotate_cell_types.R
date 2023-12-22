# Sergio Alías, 20231020
# Last modified 20231215

# Script for adding the maual annotation to the Seurat objects (human and mouse)

library(Seurat)
library(optparse)


option_list <- list(
  make_option(c("-s", "--seu"), type = "character",
              help="Path to the Seurat object"),
  make_option(c("-c", "--celltypes"), type = "character",
              help="Path to the file with the cell types"),
)  

opt <- parse_args(OptionParser(option_list = option_list))

# Main

seu <- readRDS(opt$seu)

new.cluster.ids <- read.table(file = opt$celltypes, sep = '\t', header = FALSE)[,2]
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)

saveRDS(seu, file = opt$seu)