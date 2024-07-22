library(tidyverse)

coms <- dplyr::bind_rows(
  read.table("raw_data/blast/blastn_tabular_coiA.txt"),
  read.table("raw_data/blast/blastn_tabular_comCDE.txt")
) %>% 
  dplyr::rename_with(~ c("qseqid", "sseqid", "pident", "length", "mismatch",
                     "gapopen", "qstart", "qend", "sstart", "send",
                     "evalue", "bitscore"))
