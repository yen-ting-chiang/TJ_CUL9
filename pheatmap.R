setwd("C:/Users/dannyj/TDocuments/Rproject/TJ_CUL9")
getwd()
library(dplyr)
library(pheatmap)
library(colorRampPalette)
library(RColorBrewer)

oncogenic_signature_for_heatmap <- 
  read.csv(file = "oncogenic_signature_for_heatmap.csv", 
           header = T)
oncogenic_signature_for_heatmap_tmp = 
  oncogenic_signature_for_heatmap[,c(3, 4, 5)]
row.names(oncogenic_signature_for_heatmap_tmp) = 
  oncogenic_signature_for_heatmap[,2]

pheatmap(oncogenic_signature_for_heatmap_tmp,
         color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                 name = "YlGnBu")))(100),
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color = NA,
         na_col = "white")

# install.packages("colorRampPalette")
# install.packages("brewer.pal")
# install.packages("pheatmap")