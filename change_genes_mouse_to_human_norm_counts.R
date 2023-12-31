setwd("outs")

ee <- read.csv("normalized_num_cells_expressing_per_sample_MOUSE.csv", row.names=1)

qq <- orthogene::convert_orthologs(as.matrix(ee), input_species="mouse", output_species="human", agg_fun="sum", gene_output="columns")

write.table(qq, file="human_ort_norm_num_cells.csv", quote=FALSE, sep="\t")
