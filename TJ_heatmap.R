setwd("C:/Users/dannyj/Documents/Rproject/TJ_CUL9")
getwd()
library(dplyr)
library(readr)
library(fread)
install.packages("read_delim")

filenames <- list.files(pattern="*.tsv", full.names=TRUE)
ldf <- lapply(filenames, read.csv, sep = "\t")
res <- lapply(ldf, summary)
names(ldf) <- substr(filenames, 19,100)
names(res) <- substr(filenames, 19,100)



D0_kegg <- 
  bind_rows(ldf[["RBD_D0_kegg.tsv"]], 
            ldf[["RBM_D0_kegg.tsv"]])
D15_kegg <- 
  bind_rows(ldf[["RBD_D15_kegg.tsv"]], 
            ldf[["RBM_D15_kegg.tsv"]])
D24_kegg <- 
  bind_rows(ldf[["RBD_D24_kegg.tsv"]], 
            ldf[["RBM_D24_kegg.tsv"]])
D0_transcription_factor_targets <- 
  bind_rows(ldf[["RBD_D0_transcription_factor_targets.tsv"]], 
            ldf[["RBM_D0_transcription_factor_targets.tsv"]])
D15_transcription_factor_targets <- 
  bind_rows(ldf[["RBD_D15_transcription_factor_targets.tsv"]], 
            ldf[["RBM_D15_transcription_factor_targets.tsv"]])
D24_transcription_factor_targets <- 
  bind_rows(ldf[["RBD_D24_transcription_factor_targets.tsv"]], 
            ldf[["RBM_D24_transcription_factor_targets.tsv"]])

D0_oncogenic_signature <- 
  bind_rows(ldf[["RBD_D0_oncogenic_signature.tsv"]], 
            ldf[["RBM_D0_oncogenic_signature.tsv"]])
D15_oncogenic_signature <- 
  bind_rows(ldf[["RBD_D15_oncogenic_signature.tsv"]], 
            ldf[["RBM_D15_oncogenic_signature.tsv"]])
D24_oncogenic_signature <- 
  bind_rows(ldf[["RBD_D24_oncogenic_signature.tsv"]], 
            ldf[["RBM_D24_oncogenic_signature.tsv"]])

D0_kegg_filtered <- D0_kegg %>% 
  filter(NOM.p.val<=0.05) %>% 
  filter(FDR.q.val<=0.25)
D15_kegg_filtered <- D15_kegg %>% 
  filter(NOM.p.val<=0.05) %>% 
  filter(FDR.q.val<=0.25)
D24_kegg_filtered <- D24_kegg %>% 
  filter(NOM.p.val<=0.05) %>% 
  filter(FDR.q.val<=0.25)

D0_transcription_factor_targets_filtered <- 
  D0_transcription_factor_targets %>% 
  filter(NOM.p.val<=0.005) %>% 
  filter(FDR.q.val<=0.025)
D15_transcription_factor_targets_filtered <- 
  D15_transcription_factor_targets %>% 
  filter(NOM.p.val<=0.005) %>% 
  filter(FDR.q.val<=0.025)
D24_transcription_factor_targets_filtered <- 
  D24_transcription_factor_targets %>% 
  filter(NOM.p.val<=0.005) %>% 
  filter(FDR.q.val<=0.025)

D0_oncogenic_signature_filtered <- D0_oncogenic_signature %>% 
  filter(NOM.p.val<=0.05) %>% 
  filter(FDR.q.val<=0.25)
D15_oncogenic_signature_filtered <- D15_oncogenic_signature %>% 
  filter(NOM.p.val<=0.05) %>% 
  filter(FDR.q.val<=0.25)
D24_oncogenic_signature_filtered <- D24_oncogenic_signature %>% 
  filter(NOM.p.val<=0.05) %>% 
  filter(FDR.q.val<=0.25)


kegg_filtered_join_1 = full_join(D0_kegg_filtered,
                   D15_kegg_filtered,
                   by = c("NAME" = "NAME"),
                   suffix = c(".D0", ".D15"),
                   keep = FALSE)
kegg_filtered_join_2 = full_join(kegg_filtered_join_1,
                        D24_kegg_filtered,
                        by = c("NAME" = "NAME"),
                        keep = FALSE)

transcription_factor_targets_filtered_join_1 = 
  full_join(D0_transcription_factor_targets_filtered,
            D15_transcription_factor_targets_filtered,
            by = c("NAME" = "NAME"),
            suffix = c(".D0", ".D15"),
            keep = FALSE)

transcription_factor_targets_filtered_join_2 = 
  full_join(transcription_factor_targets_filtered_join_1,
            D24_transcription_factor_targets_filtered,
            by = c("NAME" = "NAME"),
            keep = FALSE)

oncogenic_signature_filtered_join_1 = 
  full_join(D0_oncogenic_signature_filtered,
            D15_oncogenic_signature_filtered,
            by = c("NAME" = "NAME"),
            suffix = c(".D0", ".D15"),
            keep = FALSE)

oncogenic_signature_filtered_join_2 = 
  full_join(oncogenic_signature_filtered_join_1,
            D24_oncogenic_signature_filtered,
            by = c("NAME" = "NAME"),
            keep = FALSE)


kegg_for_heatmap <- kegg_filtered_join_2 %>% 
  arrange(desc(NES)) %>% 
  select(NAME, NES.D0, NES.D15, NES) %>% 
  filter(is.na(NES) == FALSE)

transcription_factor_targets_for_heatmap <- 
  transcription_factor_targets_filtered_join_2 %>% 
  arrange(desc(NES)) %>% 
  select(NAME, NES.D0, NES.D15, NES) %>% 
  filter(is.na(NES) == FALSE)

oncogenic_signature_for_heatmap <- 
  oncogenic_signature_filtered_join_2 %>% 
  arrange(desc(NES)) %>% 
  select(NAME, NES.D0, NES.D15, NES) %>% 
  filter(is.na(NES) == FALSE)

write.csv(kegg_for_heatmap, 
          file = "kegg_for_heatmap.csv")
write.csv(transcription_factor_targets_for_heatmap, 
          file = "transcription_factor_targets_for_heatmap.csv")
write.csv(oncogenic_signature_for_heatmap, 
          file = "oncogenic_signature_for_heatmap.csv")
