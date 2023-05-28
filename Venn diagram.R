#Import the files
all_north <- read.csv("/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/log2_2/analisi_definitiu_m/fulla/all_north/SigGenes_all_north_.csv")
north_inland <- read.csv("/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/log2_2/analisi_definitiu_m/fulla/north_inland/SigGenes_north_inland.csv")
north_costal <- read.csv("/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/log2_2/analisi_definitiu_m/fulla/north_costal/SigGenes_north_costal_noout.csv")
head(north_costal)

library(tidyverse)
library(VennDiagram)
#create the array of each list
AA <- all_north %>% 
  unlist()
BB <- north_inland %>% 
  unlist()
CC <- north_costal %>% 
  unlist()

#perform function to create the venn diagram
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
total_list <- list(A=AA, B=BB, C = CC) 
display_venn(total_list)

# Further customization
display_venn(
  total_list,
  category.names = c("All north", "North inlad", "North coastal"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#00a087", "#4c5488", "#ff745f"),
  # Numbers
  cex = 3,
  fontface = "italic",
  # Set names
  cat.cex =  2.5,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

#To extract the common list
ol <- calculate.overlap(total_list)
ol$a3


