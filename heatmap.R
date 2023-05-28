#Import the data 
DEGs_comparison_expressio_counts <- read.csv("/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/log2_2/analisi_definitiu_m/fulla/comparison_far_central/DEGs_comparison_south_expressio_counts.csv")
#Select the columns of interest
X <- DEGs_comparison_expressio_counts[,2:7]

# Re-elaborate dataframe scaling
mat3=scale(X ,scale=T, center = T)
m <- data.matrix(mat3, rownames.force = NA)

#Genrate heatmap
library(pheatmap)
heatmap <- pheatmap(m, scale="row", Kmeans_K = NA, 
                    clustering_method = "ward.D2", breaks = NA, cluster_row = TRUE,
                    cluster_cols = T, main = "DEGs exclusive south", 
                    treeheight_col = 25, treeheight_row = 25,
                    cellwidth = 20, show_rownames = F,
                    labels_col = c("ST TREAT", "ST TREAT", "ST TREAT", "ST CONT","ST CONT", "ST CONT"),
                    color = colorRampPalette(c("#430c55", "white", "#234F1E"))(50))
