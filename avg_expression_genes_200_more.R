# Sergio Alías, 20231114
# Last modified 20231114

# Script for getting top 200 expressed genes present in at least 50% cells per cluster

library(Seurat)
library(scCustomize)
library(optparse)


option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help="human or mouse"),
  make_option(c("-o", "--output"), type = "character",
              help="Output folder")
)  

opt <- parse_args(OptionParser(option_list = option_list))

seu <- readRDS(file.path(opt$output, paste0("seu.", opt$input, ".RDS")))

percent_stats <- Percent_Expressing(seu, features = rownames(seu), threshold = 0)
avg_expression <- AverageExpression(seu)$RNA

top200_df <- vector(mode = "list", length = length(levels(Idents(seu))))
names(top200_df) <- levels(Idents(seu))

for (i in 1:(length(levels(Idents(seu))))){ 
  clus_avg <- avg_expression[c(rownames(percent_stats[percent_stats[,i] >= 50,, drop = FALSE])), i, drop = FALSE]
  if (nrow(clus_avg) > 200){
    clus_top200 <- clus_avg[order(clus_avg[,1], decreasing = TRUE),, drop = FALSE][1:200,, drop = FALSE]
  } else {
  	clus_top200 <- clus_avg[order(clus_avg[,1], decreasing = TRUE),, drop = FALSE]
  }
  top200_df[[i]] <- clus_top200
}

names(top200_df) <- gsub("X[0-9]+\\.\\.", "", names(top200_df))

names(top200_df) <- gsub("^\\d+\\.\\s", "", names(top200_df))

top200_conv_df <- lapply(top200_df, function(vec) {
  data.frame(
    names = rownames(vec),
    values = as.numeric(vec),
    stringsAsFactors = FALSE
  )
})

full_table <- dplyr::bind_rows(top200_conv_df, .id = "column_label")
genes <- sort(unique(full_table$names))
gene_lines <- sapply(genes, function(g) {
  paste(unique(full_table$column_label[full_table$names == g]), collapse=",")
})
gene_cell_subtypes <- data.frame(gene_lines)

write.table(gene_cell_subtypes, col.names=FALSE, file=file.path(opt$output, paste0("gene_cell_subtypes_top200_", opt$input, ".tsv")), sep="\t", quote=FALSE)


