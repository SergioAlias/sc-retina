# Sergio Alías, 20230920
# Last modified 20231221

# Script for aggregating single-cell data to gene normcounts per sample and norm number of cells expressing each gene
# This script can generate ítems 3 and 5 in the Obsidian (and 4 if you remove the normalization)

library(Seurat)
library(optparse)


option_list <- list(
  make_option(c("-d", "--dataset"), type = "character",
              help="Dataset to use: 'human' or 'mouse'"),
  make_option(c("-o", "--output"), type = "character",
              help="Output to obtain: 'ncounts' or 'ncells'. Read the Obsidian for more info."),
)  

opt <- parse_args(OptionParser(option_list = option_list))


if (opt$dataset == "human"){
dataset_loc <- "/mnt2/fscratch/users/bio_267_uma/sergioalias/NGS_projects/MACULAR/count_results"
ids <- c("Un_Peri",
         "Un_Macu",
         "AMD_Peri",
         "AMD_Macu",
         "Adult_2",
         "Adult_5",
         "Adult_3",
         "Adult_1",
         "Adult_4",
         "21PCW",
         "20PCW",
         "12PCW",
         "16PCW")
} else if (opt$dataset == "mouse"){
  dataset_loc <- "/mnt2/fscratch/users/bio_267_uma/sergioalias/NGS_projects/MOUSE/count_results"
  ids <- c("WT_1",
           "WT_2",
           "rd10_1",
           "rd10_2")
}


d10x.data <- sapply(ids, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"cellranger_0000",i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

experiment.data <- do.call("cbind", d10x.data)

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "macular",
  min.cells = 1,
  min.features = 1,
  names.field = 2,
  names.delim = "\\-")


if (opt$output == "ncounts"){ # if you want normalized counts
  norm.exp <- NormalizeData(experiment.aggregate)    
} else if (opt$output == "ncells"){ # if you want normalized number of cells expressing gene
  norm.exp <- experiment.aggregate # does nothing but makes me not changing code below                  
}

col_names <- norm.exp@meta.data$orig.ident
unique_col_names <- unique(col_names)

mat <- as.matrix(norm.exp@assays$RNA@data)


result_matrix <- matrix(0, nrow = nrow(mat), ncol = length(unique_col_names))


for (i in 1:length(unique_col_names)) {
  col_name <- unique_col_names[i]
  col_indices <- which(col_names == col_name)
  if (opt$output == "ncounts"){ # if you want normalized counts
    result_matrix[, i] <- rowSums(mat[, col_indices, drop = FALSE])    
  } else if (opt$output == "ncells"){ # if you want normalized number of cells expressing gene
    result_matrix[, i] <- apply(mat, 1, function(row, col_indices){sum(row[col_indices] != 0) / length(col_indices)}, col_indices = col_indices)                   
  }
}

rownames(result_matrix) <- rownames(mat)
colnames(result_matrix) <- unique_col_names


if (opt$output == "ncounts"){
  if (opt$dataset == "human"){
    filename <- "norm_counts_per_sample.csv"
  } else if (opt$dataset == "mouse"){
    filename <- "norm_counts_per_sample_MOUSE.csv"
  }    
} else if (opt$output == "ncells"){ 
  if (opt$dataset == "human"){
    filename <- "normalized_num_cells_expressing_per_sample.csv"
  } else if (opt$dataset == "mouse"){
    filename <- "normalized_num_cells_expressing_per_sample_MOUSE.csv"
  }                      
}

write.csv(result_matrix, file = filename) 
print("Finished! :D")