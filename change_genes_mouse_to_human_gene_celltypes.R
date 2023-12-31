# Este script también sirve para el fichero de DEgenes

setwd("outs")

ee <- read.csv("gene_cell_subtypes_top200_mouse.tsv", sep = "\t", header=FALSE)

rownames(ee) <- ee$V1

qq <- orthogene::convert_orthologs(as.matrix(ee), input_species="mouse", output_species="human", agg_fun="sum", gene_output="columns")

qq <- qq[,-1, drop = FALSE]

write.table(qq, file="human_ort_gene_cell_subtypes_top200_mouse.tsv", quote=FALSE, sep="\t", col.names=FALSE)
