# Sergio Alías, 20231219
# Last modified 20231220

# Script for finding markers between WT and mutant (rd10) genotypes in mouse dataset

library(Seurat)

setwd("outs")

seu <- readRDS("seu.mouse.RDS")

seu$stage <- ifelse(seu$samplename %in% c("WT_1", "WT_2"), "WT", "mutant")

marker_list <- vector(mode = "list", length = length(levels(Idents(seu))))
names(marker_list) <- levels(Idents(seu))

seu$celltype.condition <- paste(Idents(seu), seu$stage, sep="_")
seu$celltype <- Idents(seu)
Idents(seu) <- "celltype.condition"
for (i in levels(seu$celltype)){
  ident1 <- paste0(i,"_WT")
  ident2 <- paste0(i,"_mutant")
  marker_list[[i]] <- tryCatch(
        {
            markers <- FindMarkers(seu, ident.1 = ident1, ident.2=ident2, min.pct=0.01, logfc.threshold=0.1, only.pos = FALSE)
            markers <- markers[markers$p_val_adj <= 0.05,]
            markers$gene.name <- rownames(markers)
            markers
        },
        error = function(cond) {
            message(paste("Error handling", i))
            message("Here's the original error message:")
            message(conditionMessage(cond))
            message("We'll remove it from the output")
            NULL
        }
    )
  
}
saveRDS(marker_list, file = "mouse_markers_logfc01_pval05.RDS")

cl_names <- gsub(" ", ".", names(marker_list))
names(marker_list) <- gsub("[0-9]+\\.\\.", "", cl_names)
full_table <- dplyr::bind_rows(marker_list, .id = "column_label")
genes <- sort(unique(full_table$gene.name))

gene_lines <- sapply(genes, function(g) {
  paste(unique(full_table$column_label[full_table$gene.name == g]), collapse=",")
})

gene_cell_subtypes <- data.frame(gene_lines)
write.table(gene_cell_subtypes, col.names=FALSE, file="DEgenes_mouse_WT_vs_mutant_cell_subtypes_logFC01_pval05.tsv", sep="\t", quote=FALSE)