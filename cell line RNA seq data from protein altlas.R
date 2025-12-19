#Transcript expression levels summarized per gene in 1206 cell lines of which 1132 are classified as cancer cell lines (see column Cancer cell lines in the Cell lines table below). The tab-separated file includes Ensembl gene identifier, analysed cell line, transcripts per million ("TPM"), protein-coding transcripts per million ("pTPM") and normalized expression ("nTPM"). The data is based on The Human Protein Atlas version 25.0 and Ensembl version 109.
#https://www.proteinatlas.org/humanproteome/cell+line/data#cell_lines
#be aware that many cell lines are not included in this database, such as HK-2

#mouse cell line seq data from ENCODE
#https://www.encodeproject.org/rnaget-report/?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.biosample_ontology.classification=cell+line

library(readr)
library(tidyverse)
setwd("S:/Ye/Tool Dev/cell line NGS proteomic source/")

#human cell line dataset
con <- unz("S:/Ye/Tool Dev/cell line NGS proteomic source/rna_celline.tsv.zip", "rna_celline.tsv")
dat <- read_tsv(con)

head(dat)

getwd()
dat_H4 <- dat %>%
  filter(`Cell line` %in% "H4")

write.table(dat_H4, "H4_RNA_proteinaltlas.txt", quote = F, sep = "\t", row.names = F)

#neuron target
GOI <- c("GRN", "SCN1A", "SHANK3", "GBA1", "FGF21", "SLC5A2", "MSH3", "LDLR", "STXBP1", "KL", "KLB")

dat_H4_GOI <- dat_H4 %>%
  filter(`Gene name` %in% GOI) %>%
  arrange(desc(TPM)) %>%
  rowwise() %>%
  mutate(percen2top = (100*TPM/max(.$TPM)))


head(dat_H4_GOI)


write.table(dat_H4_GOI, "H4_RNA_proteinaltlas_GOI.txt", quote = F, sep = "\t", row.names = F)

#GOI <- "KL"

dat_GOI <-  dat %>%
  select(2:4) %>%
  filter(`Gene name` %in% GOI) %>%
  filter(TPM > 0) %>%
  group_by(`Gene name`) %>%
  slice_max(order_by = TPM, n = 20, with_ties = FALSE) %>%
  mutate(cell_TPM = paste0(`Cell line`, "(", TPM, ")")) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  select(cell_TPM, `Gene name`, rank) %>%
  pivot_wider(
    names_from = `Gene name`,
    values_from = cell_TPM
  )

head(dat_GOI)

write.table(dat_GOI, "top20 expression in 1206 cell lines GOI.txt", quote = F, sep = "\t", row.names = F)


#close(con)

#mouse cell line database

dat_ms <- read_tsv("rna_expression_report_2025_12_12_9h_29m.tsv", skip = 1)
head(dat_ms)

dat_ms_GOI <- dat_ms %>%
  filter(`Gene symbol` == "Kl")
