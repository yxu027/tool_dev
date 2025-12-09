# install.packages(c("biomaRt", "dplyr"))
library(biomaRt)
library(dplyr)

check_exon_frame_biomart <- function(transcript_id,
                                     dataset = "hsapiens_gene_ensembl") {
  library(biomaRt)
  library(dplyr)
  
  ensembl <- useEnsembl(
    biomart = "ensembl",
    dataset = dataset,
    host = "https://www.ensembl.org"
  )
  
  # ---- Fetch gene + protein info ----
  gene_info <- getBM(
    attributes = c("ensembl_transcript_id",
                   "external_gene_name",
                   "ensembl_peptide_id",
                   "uniprotswissprot"),
    filters = "ensembl_transcript_id",
    values = transcript_id,
    mart = ensembl
  )
  
  gene_name <- unique(gene_info$external_gene_name)[1]
  pep_id    <- unique(gene_info$ensembl_peptide_id)[1]
  uniprot   <- unique(gene_info$uniprotswissprot)[1]
  
  # ---- Fetch exon info ----
  attrs <- c(
    "ensembl_transcript_id", "ensembl_exon_id", "rank",
    "phase", "end_phase", "cds_start", "cds_end",
    "chromosome_name", "strand", "exon_chrom_start", "exon_chrom_end",
    "transcript_biotype"
  )
  
  res <- getBM(
    attributes = attrs,
    filters = "ensembl_transcript_id",
    values = transcript_id,
    mart = ensembl
  )
  if (nrow(res) == 0) stop("No exon results. Check transcript ID or species dataset.")
  
  df <- res %>%
    arrange(rank, exon_chrom_start) %>%
    mutate(
      phase      = suppressWarnings(as.integer(phase)),
      end_phase  = suppressWarnings(as.integer(end_phase)),
      cds_bp     = ifelse(is.na(cds_start) | is.na(cds_end), 0L,
                          as.integer(cds_end) - as.integer(cds_start) + 1L),
      coding_exon = ((phase >= 0 | end_phase >= 0) & cds_bp > 0),
      starts_on_codon = (phase == 0 & coding_exon),
      ends_on_codon   = (end_phase == 0 & coding_exon)
    ) %>%
    arrange(rank) %>%
    mutate(
      prev_end_phase = lag(end_phase),
      next_phase     = lead(phase),
      prev_coding    = lag(coding_exon),
      next_coding    = lead(coding_exon)
    )
  
  crit_len_frameshift   <- df$coding_exon & (df$cds_bp %% 3 != 0)
  can_eval_neighbor     <- df$coding_exon & df$prev_coding & df$next_coding &
    !is.na(df$prev_end_phase) & !is.na(df$next_phase) &
    df$prev_end_phase >= 0 & df$next_phase >= 0
  crit_neighbor_incompat <- can_eval_neighbor & (df$prev_end_phase != df$next_phase)
  
  df <- df %>%
    mutate(
      affects_frame_if_skipped = crit_len_frameshift | crit_neighbor_incompat,
      affects_frame_reason = case_when(
        crit_len_frameshift & crit_neighbor_incompat ~ "len%3!=0 + neighbor_phase_mismatch",
        crit_len_frameshift ~ "len%3!=0",
        crit_neighbor_incompat ~ "neighbor_phase_mismatch",
        TRUE ~ NA_character_
      ),
      gene_name = gene_name,
      peptide_id = pep_id,
      uniprot_id = uniprot
    )
  
  # ---- Save CSV with gene name in filename ----
  out_csv <- paste0(gene_name, "_", transcript_id, "_exon_frames.csv")
  utils::write.csv(df, out_csv, row.names = FALSE)
  
  message("Saved: ", out_csv)
  
  list(
    full_table = df,
    prf_affecting_exons = df %>% filter(affects_frame_if_skipped)
  )
}

## Example run:
# out <- check_exon_frame_biomart("ENST00000266088")
# head(out$full_table)


## Example:
out <- check_exon_frame_biomart("ENST00000265081")
View(out$full_table)
View(out$prf_affecting_exons)

cat(length(out$prf_affecting_exons$rank), " of ", length(out$full_table$rank),
    " exons would shift the ORF if spliced out (exon ranks: ",
    paste(out$prf_affecting_exons$rank, collapse = ", "), ").\n", sep = "")


# SGLT2, ENST00000330498
# SGLT1, ENST00000266088
# MSH3, ENST00000265081


ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = "https://www.ensembl.org"
)
listAttributes(ensembl)[1:20, ]
listAttributes(ensembl)[grep("exon", listAttributes(ensembl)$name), ]
