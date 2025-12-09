#Transcript expression levels summarized per gene in 1206 cell lines of which 1132 are classified as cancer cell lines (see column Cancer cell lines in the Cell lines table below). The tab-separated file includes Ensembl gene identifier, analysed cell line, transcripts per million ("TPM"), protein-coding transcripts per million ("pTPM") and normalized expression ("nTPM"). The data is based on The Human Protein Atlas version 25.0 and Ensembl version 109.
#https://www.proteinatlas.org/humanproteome/cell+line/data#cell_lines

library(readr)
setwd("S:/Ye/Tool Dev/cell line NGS proteomic source/")
con <- unz("S:/Ye/Tool Dev/cell line NGS proteomic source/rna_celline.tsv.zip", "rna_celline.tsv")
dat <- read_tsv(con)

head(dat)

getwd()
dat_H4 <- dat %>%
  filter(`Cell line` %in% "H4")

write.table(dat_H4, "H4_RNA_proteinaltlas.txt", quote = F, sep = "\t", row.names = F)

#neuron target
GOI <- c("GRN", "SCN1A", "SHANK3", "GBA1", "FGF21", "SGLT2", "MSH3", "LDLR", "STXBP1", "KL", "KLB")

dat_H4_GOI <- dat_H4 %>%
  filter(`Gene name` %in% GOI) %>%
  arrange(desc(TPM))

write.table(dat_H4_GOI, "H4_RNA_proteinaltlas_GOI.txt", quote = F, sep = "\t", row.names = F)

close(con)
