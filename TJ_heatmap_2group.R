setwd("C:/R_project/TJ_CUL9")
getwd()
library(dplyr)
library(readr)
# install.packages("read_delim")

filenames <- list.files(pattern="gsea_report", full.names=TRUE)
ldf <- lapply(filenames, read.csv, sep = "\t")
res <- lapply(ldf, summary)
names(ldf) <- substr(filenames, 19,100)
names(res) <- substr(filenames, 19,100)

D15_GO_BP <- 
  bind_rows(ldf[["RBD_D15_GO_BP.tsv"]], 
            ldf[["RBM_D15_GO_BP.tsv"]])
D24_GO_BP <- 
  bind_rows(ldf[["RBD_D24_GO_BP.tsv"]], 
            ldf[["RBM_D24_GO_BP.tsv"]])

D15_GO_BP_filtered <- D15_GO_BP %>% 
  filter(NOM.p.val<=0.05) %>% 
  filter(FDR.q.val<=0.25)
D24_GO_BP_filtered <- D24_GO_BP %>% 
  filter(NOM.p.val<=0.05) %>% 
  filter(FDR.q.val<=0.25)

GO_BP_filtered_join_1 = full_join(D15_GO_BP_filtered,
                                 D24_GO_BP_filtered,
                                 by = c("NAME" = "NAME"),
                                 suffix = c(".D15", ".D24"),
                                 keep = FALSE)

GO_BP_for_heatmap <- GO_BP_filtered_join_1 %>% 
  arrange(desc(NES.D24)) %>% 
  select(NAME, NES.D15, NES.D24) %>% 
  filter(is.na(NES.D24) == FALSE)

write.csv(GO_BP_for_heatmap, 
          file = "GO_BP_for_heatmap.csv")
