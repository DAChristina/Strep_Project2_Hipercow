library(tidyverse)
library(ggtree)
library(epitools)

# load data
data <- read.csv("raw_data/gubbins/n703/phandango_microreact_check/microreact_tre_names.csv")

gene_data_protl <- read.csv("raw_data/panaroo/gene_data_protl.csv", head = F)
colnames(gene_data_protl) <- c("sequence", "protein", "prot_length")

name <- "01435762_ukhsa_assembly_trimmed_500bp_contigs"

coiA_samples <- read.csv("raw_data/panaroo/protein_samples/coiA_samples.csv", head = F)
comA_samples <- read.csv("raw_data/panaroo/protein_samples/comA_samples.csv", head = F)
comC_samples <- read.csv("raw_data/panaroo/protein_samples/comC_samples.csv", head = F)
comD_samples <- read.csv("raw_data/panaroo/protein_samples/comD_samples.csv", head = F)
comE_samples <- read.csv("raw_data/panaroo/protein_samples/comE_samples.csv", head = F)


GNAT_samples <- read.csv("raw_data/panaroo/protein_samples/GNAT_acetyltransferase_samples.csv", head = F)
blpA2_samples <- read.csv("raw_data/panaroo/protein_samples/blpA2_samples.csv", head = F)
bsaA_samples <- read.csv("raw_data/panaroo/protein_samples/bsaA_samples.csv", head = F)
cinA_samples <- read.csv("raw_data/panaroo/protein_samples/cinA_samples.csv", head = F)
folP_samples <- read.csv("raw_data/panaroo/protein_samples/folP_samples.csv", head = F)
hysA_samples <- read.csv("raw_data/panaroo/protein_samples/hysA_samples.csv", head = F)
liaF_samples <- read.csv("raw_data/panaroo/protein_samples/liaF_samples.csv", head = F)

# Focused on competence proteins ###############################################
coiA_samples <- as.data.frame(t(coiA_samples)) %>% 
  magrittr::set_colnames(c("all_protein_coiA")) %>% 
  dplyr::filter(!is.na(all_protein_coiA),
                all_protein_coiA != "group_843" & all_protein_coiA != "hypothetical protein") %>% 
  dplyr::mutate(seq_name_coiA = substr(all_protein_coiA, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_coiA), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_coiA" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_coiA) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_coiA = prot_length) %>% 
  dplyr::ungroup()

comA_samples <- as.data.frame(t(comA_samples)) %>% 
  magrittr::set_colnames(c("all_protein_comA")) %>% 
  dplyr::filter(!is.na(all_protein_comA),
                all_protein_comA != "msbA~~~lagD_3~~~lagD_1" & all_protein_comA != "msbA;lagD_3;lagD_1" & 
                  all_protein_comA != "Lipid A export ATP-binding/permease protein MsbA;Lactococcin-G-processing and transport ATP-binding protein LagD") %>% 
  dplyr::mutate(seq_name_comA = substr(all_protein_comA, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_comA), sep = ";") %>%
  dplyr::left_join(gene_data_protl, by = c("all_protein_comA" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_comA) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_comA = prot_length) %>% 
  dplyr::ungroup()

# For comC it's better to separate the df into comC and comC21
comC_samples_1 <- as.data.frame(t(comC_samples)) %>% 
  magrittr::set_colnames(c("all_protein_comC", "all_protein_comC21")) %>% 
  dplyr::select(all_protein_comC) %>% 
  dplyr::filter(!is.na(all_protein_comC),
                all_protein_comC != "comC" & all_protein_comC != "Type 4 prepilin-like proteins leader peptide-processing enzyme") %>% 
  dplyr::mutate(seq_name_comC = substr(all_protein_comC, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_comC), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_comC" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_comC) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_comC = prot_length) %>% 
  dplyr::ungroup()

comC21_samples <- as.data.frame(t(comC_samples)) %>% 
  magrittr::set_colnames(c("all_protein_comC", "all_protein_comC21")) %>% 
  dplyr::select(all_protein_comC21) %>% 
  dplyr::filter(!is.na(all_protein_comC21),
                all_protein_comC21 != "comC" & all_protein_comC21 != "Type 4 prepilin-like proteins leader peptide-processing enzyme") %>% 
  dplyr::mutate(seq_name_comC21 = substr(all_protein_comC21, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_comC21), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_comC21" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_comC21) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_comC21 = prot_length) %>% 
  dplyr::ungroup()

comD_samples <- as.data.frame(t(comD_samples)) %>% 
  magrittr::set_colnames(c("all_protein_comD")) %>% 
  dplyr::filter(!is.na(all_protein_comD),
                all_protein_comD != "group_89" & all_protein_comD != "hypothetical protein") %>% 
  dplyr::mutate(seq_name_comD = substr(all_protein_comD, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_comD), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_comD" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_comD) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_comD = prot_length) %>% 
  dplyr::ungroup()

comE_samples <- as.data.frame(t(comE_samples)) %>% 
  magrittr::set_colnames(c("all_protein_comE")) %>% 
  dplyr::filter(!is.na(all_protein_comE),
                all_protein_comE != "ComEC" & all_protein_comE != "comEC" & all_protein_comE != "ComE operon protein 3") %>% 
  dplyr::mutate(seq_name_comE = substr(all_protein_comE, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_comE), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_comE" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_comE) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_comE = prot_length) %>% 
  dplyr::ungroup()

# Focused on other "additional" proteins! ######################################
# GNAT family N-acetyltransferase (Acetylation linked to AMR)
GNAT_samples <- as.data.frame(t(GNAT_samples)) %>% 
  magrittr::set_colnames(c("all_protein_GNAT")) %>% 
  dplyr::filter(!is.na(all_protein_GNAT),
                all_protein_GNAT != "group_363" & all_protein_GNAT != "hypothetical protein") %>% 
  dplyr::mutate(seq_name_GNAT = substr(all_protein_GNAT, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_GNAT), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_GNAT" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_GNAT) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_GNAT = prot_length) %>% 
  dplyr::ungroup()

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6016807/ (previously misidentified as comA)
blpA2_samples <- as.data.frame(t(blpA2_samples)) %>% 
  magrittr::set_colnames(c("all_protein_blpA2")) %>% 
  dplyr::filter(!is.na(all_protein_blpA2),
                all_protein_blpA2 != "lagD_1" & all_protein_blpA2 != "lagD_2" & all_protein_blpA2 != "lagD_1~~~lagD_2" & 
                  all_protein_blpA2 != "hypothetical protein" & all_protein_blpA2 != "Lactococcin-G-processing and transport ATP-binding protein LagD") %>% 
  dplyr::mutate(seq_name_blpA2 = substr(all_protein_blpA2, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_blpA2), sep = ";") %>%
  dplyr::left_join(gene_data_protl, by = c("all_protein_blpA2" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_blpA2) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_blpA2 = prot_length) %>% 
  dplyr::ungroup()

# glutathione peroxidase (Related to oxidative stress resistance systems (and metal ion toxicity))
bsaA_samples <- as.data.frame(t(bsaA_samples)) %>% 
  magrittr::set_colnames(c("all_protein_bsaA")) %>% 
  dplyr::filter(!is.na(all_protein_bsaA),
                all_protein_bsaA != "bsaA" & all_protein_bsaA != "Glutathione peroxidase BsaA") %>% 
  dplyr::mutate(seq_name_bsaA = substr(all_protein_bsaA, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_bsaA), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_bsaA" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_bsaA) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_bsaA = prot_length) %>% 
  dplyr::ungroup()

# competence/damage-inducible protein A (Natural competence)
cinA_samples <- as.data.frame(t(cinA_samples)) %>% 
  magrittr::set_colnames(c("all_protein_cinA")) %>% 
  dplyr::filter(!is.na(all_protein_cinA),
                all_protein_cinA != "cinA" & all_protein_cinA != "Putative competence-damage inducible protein") %>% 
  dplyr::mutate(seq_name_cinA = substr(all_protein_cinA, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_cinA), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_cinA" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_cinA) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_cinA = prot_length) %>% 
  dplyr::ungroup()

# dihydropteroate synthase (targeted by Sulfonamides, linked to AMR)
folP_samples <- as.data.frame(t(folP_samples)) %>% 
  magrittr::set_colnames(c("all_protein_folP")) %>% 
  dplyr::filter(!is.na(all_protein_folP),
                all_protein_folP != "sulA" & all_protein_folP != "Dihydropteroate synthase") %>% 
  dplyr::mutate(seq_name_folP = substr(all_protein_folP, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_folP), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_folP" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_folP) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_folP = prot_length) %>% 
  dplyr::ungroup()

# cell wall-active antibiotics response protein LiaF (Relate to sense cell envelope stress; studied in NESp)
liaF_samples <- as.data.frame(t(liaF_samples)) %>% 
  magrittr::set_colnames(c("all_protein_liaF")) %>% 
  dplyr::filter(!is.na(all_protein_liaF),
                all_protein_liaF != "sulA" & all_protein_liaF != "Dihydropteroate synthase") %>% 
  dplyr::mutate(seq_name_liaF = substr(all_protein_liaF, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_liaF), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_liaF" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_liaF) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_liaF = prot_length) %>% 
  dplyr::ungroup()

# LPXTG-anchored hyaluronate lyase (Surface enzyme linked to spreading throughout host tissue)
hysA_samples <- as.data.frame(t(hysA_samples)) %>% 
  magrittr::set_colnames(c("all_protein_hysA")) %>% 
  dplyr::filter(!is.na(all_protein_hysA),
                all_protein_hysA != "group_924" & all_protein_hysA != "Hyaluronate lyase") %>% 
  dplyr::mutate(seq_name_hysA = substr(all_protein_hysA, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein_hysA), sep = ";") %>% 
  dplyr::left_join(gene_data_protl, by = c("all_protein_hysA" = "protein")) %>% 
  dplyr::select(-sequence) %>% 
  dplyr::group_by(seq_name_hysA) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::rename(protl_hysA = prot_length) %>% 
  dplyr::ungroup()


# Create combined dataframe! ###################################################
# ggtree failed to load a tree from one huge dataframe.
combined_data <- data %>% 
  dplyr::left_join(coiA_samples, by = c("tre.tip.label" = "seq_name_coiA")) %>% 
  dplyr::left_join(comA_samples, by = c("tre.tip.label" = "seq_name_comA")) %>% 
  dplyr::left_join(comC_samples_1, by = c("tre.tip.label" = "seq_name_comC")) %>% 
  dplyr::left_join(comC21_samples, by = c("tre.tip.label" = "seq_name_comC21")) %>% 
  dplyr::left_join(comD_samples, by = c("tre.tip.label" = "seq_name_comD")) %>% 
  dplyr::left_join(comE_samples, by = c("tre.tip.label" = "seq_name_comE")) %>% 
  # Other proteins
  dplyr::left_join(GNAT_samples, by = c("tre.tip.label" = "seq_name_GNAT")) %>% 
  dplyr::left_join(blpA2_samples, by = c("tre.tip.label" = "seq_name_blpA2")) %>% 
  dplyr::left_join(bsaA_samples, by = c("tre.tip.label" = "seq_name_bsaA")) %>% 
  dplyr::left_join(cinA_samples, by = c("tre.tip.label" = "seq_name_cinA")) %>% 
  dplyr::left_join(folP_samples, by = c("tre.tip.label" = "seq_name_folP")) %>% 
  dplyr::left_join(hysA_samples, by = c("tre.tip.label" = "seq_name_hysA")) %>% 
  dplyr::left_join(liaF_samples, by = c("tre.tip.label" = "seq_name_liaF")) %>% 
  glimpse()


# BLAST Result for coiA and comCDE #############################################
blast_coiA <- read.table("raw_data/blast/blastn_tabular_coiA.txt") %>% 
  dplyr::rename_with(~ c("qseqid", "temporary_sseqid", "pident", "length", "mismatch",
                         "gapopen", "qstart", "qend", "sstart", "send",
                         "evalue", "bitscore")) %>% 
  dplyr::mutate(sseqid = temporary_sseqid) %>% 
  tidyr::separate(temporary_sseqid, into = c("temporary_ID", "temporary_contig"), sep = "contig") %>%
  dplyr::mutate(ID = paste0(substr(temporary_ID, 1, 8), "_ukhsa_assembly_trimmed_500bp_contigs"),
                contig = paste0(ID, "_", temporary_contig)) %>% 
  dplyr::group_by(ID) %>%
  dplyr::slice(which.max(pident)) %>% # Precaution: combined data below is based on max(pident) (maximum percentage of identification)
  dplyr::ungroup() %>% 
  dplyr::select(-temporary_ID, -temporary_contig, -sseqid) %>% 
  dplyr::rename_with(~ paste0(., "_blast_coiA"), -ID)

blast_comCDE <- read.table("raw_data/blast/blastn_tabular_comCDE.txt") %>% 
  dplyr::rename_with(~ c("qseqid", "temporary_sseqid", "pident", "length", "mismatch",
                         "gapopen", "qstart", "qend", "sstart", "send",
                         "evalue", "bitscore")) %>% 
  dplyr::mutate(sseqid = temporary_sseqid) %>% 
  tidyr::separate(temporary_sseqid, into = c("temporary_ID", "temporary_contig"), sep = "contig") %>%
  dplyr::mutate(ID = paste0(substr(temporary_ID, 1, 8), "_ukhsa_assembly_trimmed_500bp_contigs"),
                contig = paste0(ID, "_", temporary_contig)) %>% 
  dplyr::group_by(ID) %>%
  dplyr::slice(which.max(pident)) %>% # Precaution: combined data below is based on max(pident) (maximum percentage of identification)
  dplyr::ungroup() %>% 
  dplyr::select(-temporary_ID, -temporary_contig, -sseqid) %>% 
  dplyr::rename_with(~ paste0(., "_blast_comCDE"), -ID)

combined_data <- combined_data %>% 
  dplyr::left_join(blast_coiA, by = c("tre.tip.label" = "ID")) %>% 
  dplyr::left_join(blast_comCDE, by = c("tre.tip.label" = "ID"))

# Precaution: combined data below is based on max(pident) from BLAST result:
representative_clades <- combined_data %>% 
  dplyr::filter(tre.tip.label %in% c("01474021_ukhsa_assembly_trimmed_500bp_contigs", # clade1
                                     "01474041_ukhsa_assembly_trimmed_500bp_contigs", # clade2
                                     "01474105_ukhsa_assembly_trimmed_500bp_contigs") # clade3
  )
################################################################################
# tre_gubbins <- BactDating::loadGubbins("raw_data/gubbins/n703/n703_") # gubbins output
tre_BD <- read_rds("outputs/genomics/choosen_n703/method_strictgamma_1e6/mcmc_bacdating.rds") # BactDating output

# calculate branch point
branchp <- as.data.frame(tre_BD$CI) %>% 
  dplyr::mutate(branchp = (V1+V2)/2) %>% 
  view()

tre_names <- as.data.frame(tre_BD$tree$tip.label)
tre_names <- dplyr::left_join(tre_names, combined_data, by = c("tre_BD$tree$tip.label" = "tre.tip.label")) %>% 
  # dplyr::filter(!is.na(clade)) %>% 
  dplyr::select(-ID) %>% 
  dplyr::rename(ID = 'tre_BD$tree$tip.label') %>% 
  dplyr::mutate(current.region.name = ifelse(is.na(current.region.name), "Unknown", current.region.name),
                ageGroup7 = ifelse(is.na(ageGroup7), "Unknown", ageGroup7),
                ageGroup2 = ifelse(is.na(ageGroup2), "Unknown", ageGroup2))

write.csv(tre_names, "raw_data/tree_inputs_final.csv", row.names = FALSE)

tre_names <- read.csv("raw_data/tree_inputs_final.csv")
# Visual inspection for every prot_l in tre_names
tre_names_long <- tre_names %>%
  tidyr::pivot_longer(cols = starts_with("protl_"), names_to = "proteins", values_to = "AA_length")

# Create histogram with facet_wrap
# png("pictures/genomics/protein_length_histograms.png", width = 24, height = 20, unit = "cm", res = 1200)
ggplot(tre_names_long, aes(x = AA_length)) +
  geom_histogram(binwidth = 5, fill = "deepskyblue3", alpha = 0.7) +
  facet_wrap(~ proteins, scales = "free_x") +
  scale_y_continuous(trans = "log1p") +
  labs(title = "Histograms Amino Acid Length",
       x = "Sequence",
       y = "Count") +
  theme_bw()
# dev.off()


# Specified to mutations in coiA and comCDE (from BLAST output)
# grouped by proteins
# https://stackoverflow.com/questions/74880935/grouped-barplot-with-sd-bars-from-two-different-groups-with-ggplot2
# Shapiro-Wilk test for normality
# Data are clearly not following criteria for normal distribution

# <tobecontinued>
stat_KW <- kruskal.test(gene_length ~ clade, data = tre_names)
stat_DB <- FSA::dunnTest(gene_length ~ clade, data = tre_names, method = "bonferroni")
stat_KW <- kruskal.test(length_blast_coiA ~ clade, data = tre_names)



tobetested <- tre_names %>% 
  dplyr::select(contains(c("mismatch_blast", "length_blast"))) %>% 
  names()

stat_result <- list()

for (col in tobetested) {
  tre_names <- tre_names %>% 
    dplyr::filter(!is.na(clade))
  
  kruskal_test <- kruskal.test(tre_names[[col]] ~ tre_names$clade)
  dunn_test <- FSA::dunnTest(tre_names[[col]], tre_names$clade, method = "bonferroni")
  
  stat_result[[col]] <- list(
    kruskal = kruskal_test,
    dunn = dunn_test$res)
}

# Display results
stat_result


tre_blast_summary <- tre_names %>% 
  dplyr::filter(!is.na(clade)) %>% 
  dplyr::group_by(clade) %>% 
  dplyr::summarise(mean_l_coiA = mean(length_blast_coiA),
                   mean_m_coiA = mean(mismatch_blast_coiA),
                   mean_l_comCDE = mean(length_blast_comCDE),
                   mean_m_comCDE = mean(mismatch_blast_comCDE),
                   sd_l_coiA = sd(length_blast_coiA),
                   sd_m_coiA = sd(mismatch_blast_coiA),
                   sd_l_comCDE = sd(length_blast_comCDE),
                   sd_m_comCDE = sd(mismatch_blast_comCDE)) %>% 
  glimpse() %>% 
  tidyr::pivot_longer(cols = -clade, names_to = c("stats", "blast", "genes"), names_sep = "_") %>%
  pivot_wider(names_from = stats, values_from = value) %>% 
  dplyr::mutate(blast = case_when(blast == "l" ~ "Gene sequence length",
                                  blast == "m" ~ "Gene sequence mismatch")) %>% 
  glimpse()

# Conduct pairwise t-tests
pairwise_tests <- tre_blast_summary %>%
  dplyr::group_by(genes, blast) %>%
  dplyr::do(broom::tidy(aov(mean ~ clade, data = .))) %>%
  glimpse()

# Format the results for annotation
annotation_data <- pairwise_tests %>%
  dplyr::filter(term == "clade") %>%
  dplyr::mutate(label = case_when(p.value < 0.001 ~ "***",
                                  p.value < 0.01 ~ "**",
                                  p.value < 0.05 ~ "*",
                                  T ~ "")) %>% 
  glimpse()

png("pictures/genomics/sequence_blast_mean_within_clades.png", width = 24, height = 12, unit = "cm", res = 1200)
pic <- ggplot(tre_blast_summary, aes(x = clade, y = mean, fill = genes)) +
  geom_bar(stat = "identity", position = position_dodge(), color = NA) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), 
                position = position_dodge(0.9), width = 0.25) +
  facet_wrap(~ blast, scales = "free_y") +
  labs(title = "Length and Mismatch of coiA and comCDE Across GPSC31 Clades",
       x = "Clades of GPSC31",
       y = "Mean Value",
       fill = "Genes") +
  theme_bw()
pic +
  geom_text(data = annotation_data, aes(x = 1.5, y = max(tre_blast_summary$mean + tre_blast_summary$sd) + 0.1, label = label), 
            position = position_dodge(width = 0.9), color = "red")

dev.off()



# Add annotations to the plot

  
  
# Calculate mean difference between clades (don't think this is necessary):
mean_diff <- tre_blast_summary %>%
dplyr::filter(blast == "Gene sequence length") %>% 
dplyr::select(clade, genes, mean) %>%
tidyr::pivot_wider(names_from = clade, values_from = mean) %>%
dplyr::mutate(diff_clade1_clade2 = clade1 - clade2,
              diff_clade1_clade3 = clade1 - clade3,
              diff_clade2_clade3 = clade2 - clade3)
mean_diff

mean_diff_long <- mean_diff %>%
  tidyr::pivot_longer(cols = starts_with("diff"), names_to = "comparison", values_to = "mean_diff") %>% 
  glimpse()

ggplot(mean_diff_long, aes(x = comparison, y = mean_diff, fill = genes)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  facet_wrap(~ genes, scales = "free_y") +
  labs(title = "Mean Differences Between Clades",
       x = "Comparison",
       y = "Mean Difference",
       fill = "Type") +
  theme_minimal()






# Focused on competence genes tree #############################################
# Nice guideline: https://yulab-smu.top/treedata-book/chapter7.html
ggtree_clades <- ggtree(tre_BD$tree,
                        mrsd = 2014-07-11) %<+%
  tre_names +
  geom_tippoint(aes(color=clade))
ggtree_clades

ggtree_ageGroup2 <- ggtree(tre_BD$tree,
                           mrsd = 2014-07-11) %<+%
  tre_names +
  geom_tippoint(aes(color=ageGroup2))
ggtree_ageGroup2

ggtree_ageGroup5 <- ggtree(tre_BD$tree,
                           mrsd = 2014-07-11) %<+%
  tre_names +
  geom_tippoint(aes(color=ageGroup))
ggtree_ageGroup5

ggtree_ageGroup7 <- ggtree(tre_BD$tree,
                           mrsd = 2014-07-11) %<+%
  tre_names +
  geom_tippoint(aes(color=ageGroup7))
ggtree_ageGroup7

ggtree_vacc <- ggtree(tre_BD$tree,
                      mrsd = 2014-07-11) %<+%
  tre_names +
  geom_tippoint(aes(color=vacc))
ggtree_vacc

ggtree_mlst <- ggtree(tre_BD$tree,
                      mrsd = 2014-07-11) %<+%
  tre_names +
  geom_tippoint(aes(color=MLST_ST))
ggtree_mlst

ggtree_others <- ggtree(tre_BD$tree,
                        mrsd = 2014-07-11) %<+%
  tre_names +
  geom_tippoint(aes(color=resistance_smx)) +
  geom_nodepoint(color="#b5e521", alpha=1/4, size=5)
ggtree_others


# FINAL TREE XD ################################################################
# Colour annotations!
# Based on col_map from 1_data_preparation.R XD
col_map <- c(
  # Vaccination
  "Pre-PCV7" = "lightblue1",
  "PCV7" = "gray70",
  "PCV13" = "gray20",
  # 5 age bands 
  "<5" = "indianred4",
  "5-18" = "orange",
  "19-30" = "seagreen4",
  "31-64" = "steelblue",
  "65+" = "purple3",
  "Unknown" = "white",
  # 7 age bands
  "<2" = "indianred4", 
  "2-4" = "indianred2", 
  "5-14" = "orange",
  "15-30" = "seagreen1", # Edit the Age-band into 15-30 & 31-44
  "31-44" = "seagreen4", # Edit the Age-band into 15-30 & 31-44
  "45-64" = "steelblue",
  "65+" = "purple3",
  # 2 age bands
  "children" = "indianred2",
  "adults" = "darkblue",
  # Regions (From north to south)
  "North West" = "indianred4",
  "North East" = "steelblue",
  "Yorkshire and The Humber" = "seagreen4",
  "East Midlands" = "purple3",
  "West Midlands" = "orange",
  "East of England" = "indianred2",
  "London" = "seagreen1",
  "South East" = "deepskyblue",
  "South West" = "gold1",
  # Cases vs sequenced
  "Serotype 1 Case" = "gray75",
  "Sequenced" = "deepskyblue3",
  "Meningitis" = "green",
  "30 Day Death" = "maroon",
  # MLST
  "227" = "steelblue",
  "228" = "green",
  "306" = "gold1",
  "3446" = "indianred4",
  "3918" = "purple3",
  "-" = "white"
)

# For the clade group
temp_clade <- tre_names %>% select(c("ID", "clade"))
temp_clade <- aggregate(.~clade, temp_clade, FUN=paste, collapse=",")
clades <- lapply(temp_clade$ID, function(x){unlist(strsplit(x,split=","))})
names(clades) <- temp_clade$clade

tr <- groupOTU(tre_BD$tree, clades, "Clade")
Clade <- NULL
ggtree_clade <- ggtree(tr=tr, layout="fan", mrsd = 2014-07-11,
                       open.angle=15, size=0.75, aes(colour=Clade)) +
  scale_colour_manual(
    name="GPSC31 Clades",
    values=c("gray75","steelblue","darkgreen","red"),
    labels=c("","Clade 1", "Clade 2", "Clade 3"),
    guide=guide_legend(keywidth=0.8,
                       keyheight=0.8,
                       order=1,
                       override.aes=list(linetype=c("0"=NA,
                                                    "Clade1"=1,
                                                    "Clade2"=1,
                                                    "Clade3"=1
                       )
                       )
    )
  ) + 
  ggnewscale::new_scale_colour()
ggtree_clade

# For Vaccination era
ggtree_vacc <- ggtree_clade %<+%
  tre_names +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=vacc),
    width=10,
    offset=0.01
  ) +
  scale_fill_manual(
    name="Vaccination Era",
    values=c(col_map),
    labels=c("Pre-PCV7", "PCV7", "PCV13"),
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=3, order=3)
  ) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  )
ggtree_vacc

# For region
ggtree_region <- ggtree_vacc %<+%
  tre_names +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=tre_names$current.region.name),
    width=10,
    offset=0.05
  ) +
  scale_fill_manual(
    name="Region (Northen to Southern)",
    values=c(col_map),
    breaks = c("North West", "North East", "Yorkshire and The Humber",
               "East Midlands", "West Midlands", "East of England",
               "London", "South East", "South West", "Unknown"),
    labels =  c("North West", "North East", "Yorkshire and The Humber",
                "East Midlands", "West Midlands", "East of England",
                "London", "South East", "South West", "Unknown"),
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=3)
  ) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  )
ggtree_region

# For 7 age groups
ggtree_ageGroup7 <- ggtree_region %<+%
  tre_names +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=tre_names$ageGroup7),
    width=10,
    offset=0.05
  ) +
  scale_fill_manual(
    name="Demographic Groups (7)",
    values=c(col_map),
    breaks=c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown"),
    labels=c("<2", "2-4", "5-14", "15-30", "31-44", "45-64", "65+", "Unknown"),
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=3)
  ) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  )
ggtree_ageGroup7

# For 2 age groups
png("pictures/genomics/tree_epiData.png", width = 24, height = 12, unit = "cm", res = 1200)
ggtree_ageGroup2 <- ggtree_ageGroup7 %<+%
  tre_names +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=tre_names$ageGroup2),
    width=10,
    offset=0.05
  ) +
  scale_fill_manual(
    name="Demographic Groups (2)",
    values=c(col_map),
    breaks = c("children", "adults", "Unknown"),
    labels = c("Children (< 15)", "Adults", "Unknown"),
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=3)
  ) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  )
ggtree_ageGroup2
dev.off()

# For MLST
png("pictures/genomics/tree_epiData_plus_MLST.png", width = 24, height = 12, unit = "cm", res = 1200)
ggtree_MLST <- ggtree_ageGroup2 %<+%
  tre_names +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=as.character(tre_names$MLST_ST)),
    width=15,
    offset=0.05
  ) +
  scale_fill_manual(
    name="Sequence Type",
    values=c(col_map),
    breaks = c("227", "228", "306", "3446", "3918", "-"),
    labels = c("ST227", "ST228", "ST306", "ST3446", "ST3918", "Unknown"),
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=3)
  ) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  )
ggtree_MLST
dev.off()



################################################################################
# ggtree failed to load a tree from one huge df, better to separate info into some df
# between tree, general info, and other data (e.g. those loaded in geom_facet())

# Focused only from gathered proteins with different length based on histogram inspection
gathered_protl <- tre_names %>% 
  dplyr::select(ID, protl_coiA, protl_comA, protl_comD, protl_blpA2, protl_folP, protl_hysA, resistance_smx) %>% 
  dplyr::rename(pl_coiA = protl_coiA,
                pl_comA = protl_comA,
                pl_comD = protl_comD,
                pl_blpA2 = protl_blpA2,
                pl_folP = protl_folP,
                pl_hysA = protl_hysA,
                res_smx = resistance_smx)


prot_map <- c(
  # coiA
  "602" = "yellow",
  "605" = "gold",
  "954" = "deepskyblue4",
  # comCDE
  "819" = "yellow",
  "2210" = "deepskyblue4"
)

png("pictures/genomics/protein_length_competence_round.png", width = 24, height = 12, unit = "cm", res = 1200)
gfacet_coms_l1 <- ggtree_clade %<+%
  tre_names +
  # ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=as.character(tre_names$qend_blast_coiA)),
    width=20,
    offset=0.05
  ) +
  scale_fill_manual(
    values = c(prot_map),
    name = "End of coimCDE and coiA Alignment",
    breaks = c("602", "605", "954"),
    labels = c("Disrupted", "Disrupted", "Normal"),
    # labels = c("602 b (disrupted)", "605 b (disrupted)", "954 b (normal)"),
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=1, order=3)
    ) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  )
gfacet_coms_l1

gfacet_coms_l2 <- gfacet_coms_l1 %<+%
  tre_names +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=as.character(tre_names$qend_blast_comCDE)),
    width=20,
    offset=0.15
  ) +
  scale_fill_manual(
    values = c(prot_map),
    name = "End of comCDE Alignment",
    breaks = c(1200, 2000),
    labels = c("1.2 kb (disrupted)", "2 kb (normal)"),
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=3)
    ) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  )
gfacet_coms_l2
dev.off()

png("pictures/genomics/protein_length_others.png", width = 24, height = 12, unit = "cm", res = 1200)
gfacet_others_l <- ggtree_vacc +
  geom_facet(panel = "Length of blpA2", data = gathered_protl, geom = geom_col, 
             aes(x = protl_blpA2,
                 colour = res_smx, fill = res_smx
             ), orientation = 'y', width = .6) +
  geom_facet(panel = "Length of folP", data = gathered_protl, geom = geom_col, 
             aes(x = protl_folP,
                 colour = res_smx, fill = res_smx
             ), orientation = 'y', width = .6) +
  geom_facet(panel = "Length of hysA", data = gathered_protl, geom = geom_col, 
             aes(x = protl_hysA,
                 colour = res_smx, fill = res_smx
             ), orientation = 'y', width = .6) +
  theme_tree2(legend.position=c(.05, .85))
gfacet_others_l
dev.off()


# Focused on BLAST results for coiA and comCDE
gathered_blast_seq <- tre_names %>% 
  dplyr::select(ID, length_blast_coiA, mismatch_blast_coiA, length_blast_comCDE, mismatch_blast_comCDE) %>% 
  dplyr::rename(l_coiA = length_blast_coiA,
                m_coiA  = mismatch_blast_coiA,
                l_comCDE = length_blast_comCDE,
                m_comCDE = mismatch_blast_comCDE)

png("pictures/genomics/protein_length_competence_blast_length.png", width = 24, height = 12, unit = "cm", res = 1200)
gfacet_blast_l <- ggtree_vacc +
  geom_facet(panel = "Length of coiA genomic sequence", data = gathered_blast_seq, geom = geom_col, 
             aes(x = l_coiA), orientation = 'y', width = .6) +
  geom_facet(panel = "Length of comCDE genomic sequence", data = gathered_blast_seq, geom = geom_col, 
             aes(x = l_comCDE), orientation = 'y', width = .6) +
  theme_tree2(legend.position=c(.05, .85))
gfacet_blast_l
dev.off()

png("pictures/genomics/protein_length_competence_blast_mismatch.png", width = 24, height = 12, unit = "cm", res = 1200)
gfacet_blast_m <- ggtree_vacc +
  geom_facet(panel = "Length of coiA genomic sequence mismatch", data = gathered_blast_seq, geom = geom_col, 
             aes(x = m_coiA), orientation = 'y', width = .6) +
  geom_facet(panel = "Length of comCDE genomic sequence mismatch", data = gathered_blast_seq, geom = geom_col, 
             aes(x = m_comCDE), orientation = 'y', width = .6) +
  theme_tree2(legend.position=c(.05, .85))
gfacet_blast_m
dev.off()













# Trial viz loop!
protl_cols <- grep("^protl_", colnames(combined_data), value = T)
for (c in protl_cols) {
  # Generate plots!
  png(file = paste("pictures/genomics/", c, ".png", sep = ""), width = 24, height = 12, unit = "cm", res = 1200)
  tree_proteins <-
    ggtree_vacc +
    geom_facet(panel = "Product", data = combined_data, geom = geom_col, 
               aes(x = c, color = c, 
                   fill = c), orientation = 'y', width = .6) +
    theme_tree2(legend.position=c(.05, .85)) #+
  #ggtitle(paste("Plot ", c))
  tree_proteins
  dev.off()
}

filtered_coiA <- combined_data %>% 
  dplyr::filter(!is.na(protl_coiA)) %>% 
  dplyr::select(tre.tip.label, protl_coiA)



ggtree_vacc +
  geom_facet(panel = "Trait", data = combined_data, geom = geom_col, 
             aes(x = serotype, colour = vacc, 
                 fill = vacc), orientation = 'y', width = .6) +
  theme_tree2(legend.position=c(.05, .85))





# Trial gheatmap
gheatmap(trial_ggtree, tre_comD_clade,
         offset=5, width=0.5, font.size=3, colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("clade1", "clade2", "clade3"), 
                    values=c("steelblue", "firebrick", "darkgreen"), name="Protein Length")
# Breaks based on hist(tre_comD_length$prot_length)

# Trial geom_facet
trial_ggtree+
  geom_facet(panel = "Trait", data = tre_comD_length, geom = geom_col, 
             aes(x = prot_length, color = prot_length, 
                 fill = prot_length), orientation = 'y', width = .6) +
  theme_tree2(legend.position=c(.05, .85))


# Focused on folP tree #########################################################
tre_folP <- dplyr::left_join(tre_names, folP_combined, by = c("tre$tip.label" = "seq_name"))
tre_folP_length <- tre_folP %>% 
  dplyr::select(`tre$tip.label`, prot_length)
tre_folP_clade <- tre_folP %>% 
  dplyr::select(`tre$tip.label`, clade)

# Nice guideline: https://yulab-smu.top/treedata-book/chapter7.html
trial_ggtree <- ggtree(tre,
                       mrsd = 2014-07-11)
trial_ggtree

# Trial gheatmap
gheatmap(trial_ggtree, tre_folP_clade,
         offset=5, width=0.5, font.size=3, colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("clade1", "clade2", "clade3"), 
                    values=c("steelblue", "firebrick", "darkgreen"), name="Protein Length")
# Breaks based on hist(tre_folP_length$prot_length)

# Trial geom_facet
trial_ggtree+
  geom_facet(panel = "Trait", data = tre_folP_length, geom = geom_col, 
             aes(x = prot_length, color = prot_length, 
                 fill = prot_length), orientation = 'y', width = .6) +
  theme_tree2(legend.position=c(.05, .85))
