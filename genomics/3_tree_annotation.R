library(tidyverse)
library(ggtree)

# load data
gene_data_protl <- read.csv("raw_data/panaroo/gene_data_protl.csv", head = F)
colnames(gene_data_protl) <- c("sequence", "protein", "prot_length")

name <- "01435762_ukhsa_assembly_trimmed_500bp_contigs"

# Focused on comD ##############################################################
comD_samples <- read.csv("raw_data/panaroo/comD_samples.csv", head = F)
comD_samples <- as.data.frame(t(comD_samples)) %>% 
  magrittr::set_colnames(c("all_protein")) %>% 
  dplyr::filter(!is.na(all_protein),
                all_protein != "group_89" & all_protein != "hypothetical protein") %>% 
  dplyr::mutate(seq_name = substr(all_protein, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein), sep = ";")

# filter out all protein list in gene_data_protl
comD_combined <- dplyr::left_join(comD_samples, gene_data_protl, by = c("all_protein" = "protein")) %>% 
  dplyr::group_by(seq_name) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(all_protein, seq_name, prot_length)

# Focused on folP ##############################################################
folP_samples <- read.csv("raw_data/panaroo/folP_samples.csv", head = F)
folP_samples <- as.data.frame(t(folP_samples)) %>% 
  magrittr::set_colnames(c("all_protein")) %>% 
  dplyr::filter(!is.na(all_protein),
                all_protein != "sulA" & all_protein != "Dihydropteroate synthase") %>% 
  dplyr::mutate(seq_name = substr(all_protein, 1, nchar(name))) %>% 
  tidyr::separate_rows(c(all_protein), sep = ";")

# filter out all protein list in gene_data_protl
folP_combined <- dplyr::left_join(folP_samples, gene_data_protl, by = c("all_protein" = "protein")) %>% 
  dplyr::group_by(seq_name) %>%
  dplyr::slice(which.max(prot_length)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(all_protein, seq_name, prot_length)

################################################################################
tre_gubbins <- BactDating::loadGubbins("raw_data/gubbins/n703/n703_") # gubbins output
tre_BD <- read_rds("outputs/genomics/choosen_n703/method_strictgamma_1e6/mcmc_bacdating.rds") # BactDating output

data <- read.csv("raw_data/gubbins/n703/phandango_microreact_check/microreact_tre_names.csv")


tre_names <- as.data.frame(tre$tip.label) #%>% 
# rename(ID = 'tre$tip.label')
tre_names$ID <- substr(tre_BD$tree$tip.label, 1, 8)
tre_names <- dplyr::left_join(tre_names, data, by = c("tre$tip.label" = "tre.tip.label")) %>% 
  dplyr::filter(!is.na("clade"))

# Focused on comD tree #########################################################
tre_comD <- dplyr::left_join(tre_names, comD_combined, by = c("tre$tip.label" = "seq_name"))
tre_comD_length <- tre_comD %>% 
  dplyr::select(`tre$tip.label`, prot_length)
tre_comD_clade <- tre_comD %>% 
  dplyr::select(`tre$tip.label`, clade)

# Nice guideline: https://yulab-smu.top/treedata-book/chapter7.html
trial_ggtree <- ggtree(tre_BD$tree,
                       mrsd = 2014-07-11)
trial_ggtree

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
