library(tidyverse)

coms <- dplyr::bind_rows(
  read.table("raw_data/blast/blastn_tabular_coiA.txt"),
  read.table("raw_data/blast/blastn_tabular_comCDE.txt")
) %>% 
  dplyr::rename_with(~ c("qseqid", "temporary_sseqid", "pident", "length", "mismatch",
                     "gapopen", "qstart", "qend", "sstart", "send",
                     "evalue", "bitscore")) %>% 
  dplyr::mutate(sseqid = temporary_sseqid) %>% 
  tidyr::separate(temporary_sseqid, into = c("temporary_ID", "temporary_contig"), sep = "contig") %>%
  dplyr::mutate(ID = paste0(substr(temporary_ID, 1, 8), "_ukhsa_assembly_trimmed_500bp_contigs"),
                contig = paste0(ID, "_", temporary_contig)) %>% 
  dplyr::select(-temporary_ID, -temporary_contig)


# Generate gene list for BLAST against NCBI servers and databases:
# https://www.ncbi.nlm.nih.gov/books/NBK569856/

