# devtools::install_github("xavierdidelot/BactDating")

library(ape)
library(tidyverse)
library(BactDating)
library(readxl)
library(coda)
library(socialmixr)

source("global/all_function.R") # Collected functions stored here!

# Treedater analysis
# More info: https://github.com/emvolz/treedater
# library(treedater)
# tre_GPSC31 <- ape::read.tree("raw_data/gubbins/GPSC_31/n712_.node_labelled.final_tree.tre")
# tre_GPSC2 <- ape::read.tree("raw_data/gubbins/GPSC_2/n17_.node_labelled.final_tree.tre")

# sts as sample times (in year-decimal) in vector format
# sts_df <- read.csv("raw_data/temporary_microreact_check.csv") %>% 
#   dplyr::mutate(Earliest.specimen.date = as.Date(Earliest.specimen.date)) %>% 
#   dplyr::filter(!is.na(microreact_ID)) %>% 
#   dplyr::select(microreact_ID, Earliest.specimen.date, year)

# sts <- setNames(sts_df$year, sts_df$microreact_ID)

# Run treedater!
# WARNING! TAKE A VERY LONG TIME TO FINISH!
# treedater_GPSC2 <- treedater::dater(tre = tre_GPSC2, sts = sts, omega0 = NA)
# 
# saveRDS(treedater_GPSC2, file = "outputs/genomics/GPSC2/treedater_GPSC2.rds")
# 
# # Trial load:
# treedater_GPSC2 <- readRDS("outputs/genomics/GPSC2/treedater_GPSC2.rds")

# see help(package='BactDating') for more info
# Time-scaled tree with BactDating
# https://xavierdidelot.github.io/BactDating/articles/Staph.html

# Result is Strict Clock. That's it.

run_bactdating <- function(nbIts){
  # 1. Data wrangling ##########################################################
  # library(data.table)
  data <- read.csv("raw_data/serotype1_UKHSA_imperial_date_age_region_MOLIS_sequenced_postThesis_cleaned.csv") %>% 
    dplyr::mutate(gubbins_ID = stringr::str_remove(assembly_name, ".fasta"),
                  Earliest.specimen.date = as.Date(Earliest.specimen.date))
  
  tre <- BactDating::loadGubbins("raw_data/gubbins/GPSC_2/n17_")
  
  ##############################################################################
  tre_names <- as.data.frame(tre$tip.label) #%>% 
  # tre_names$ID <- substr(tre$tip.label, 1, 8)
  tre_names <- dplyr::left_join(tre_names, data, by = c("tre$tip.label" = "gubbins_ID"))
  
  # 2. BactDating ##############################################################
  # https://xavierdidelot.github.io/BactDating/articles/yourData.html
  d <- cbind(tre_names$year, tre_names$year+1) # 2-d matrix according to the articles above
  
  set.seed(0)
  
  dir.create("outputs/genomics/GPSC2", FALSE, TRUE)
  dir.create("pictures/genomics/GPSC2", FALSE, TRUE)
  
  res_pr <- BactDating::bactdate(tre,d,nbIts=nbIts, # Put 1e6 or 1e10 on hipercow
                                 model = "strictgamma", # "arc",
                                 showProgress = T)
  
  saveRDS(res_pr, "outputs/genomics/GPSC2/mcmc_bacdating_GPSC2.rds")
  
  # Figures!
  png("pictures/genomics/GPSC2/tree_treeCI.png", width = 24, height = 12, unit = "cm", res = 1200)
  par(mfrow = c(1,1), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
  plot(res_pr,'treeCI',show.tip.label = F)
  dev.off()
  
  png("pictures/genomics/GPSC2/tree_trace1.png", width = 24, height = 12, unit = "cm", res = 1200)
  par(mfrow = c(1,1), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
  plot(res_pr,'trace')
  dev.off()
  
  # MCMC analysis
  mcmc_result <- BactDating::as.mcmc.resBactDating(res_pr)
  
  # Calculating ESS & Acceptance Rate
  calc_ess <- ess_calculation(mcmc_result)
  write.csv(calc_ess, "outputs/genomics/GPSC2/calc_ess.csv", row.names = TRUE)
  
  # Figures! (still failed, margin error)
  png("pictures/genomics/GPSC2/tree_trace2.png", width = 24, height = 12, unit = "cm", res = 1200)
  par(mfrow = c(1,1), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
  pmcmc_trace(mcmc_result)
  dev.off()
  
  Sys.sleep(10) # wait 10 secs
  
  # Rooted tree!
  rooted_tree <- BactDating::initRoot(tre,d[,1]) # Incompatible dimensions because of d as matrix of (74,2)
  saveRDS(rooted_tree, "outputs/genomics/GPSC2/rooted_tree.rds")
  
  png("pictures/genomics/GPSC2/tree_rootedtree.png", width = 24, height = 12, unit = "cm", res = 1200)
  par(mfrow = c(1,1), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
  plot(rooted_tree, show.tip.label = F)
  dev.off()
  
  Sys.sleep(10) # wait 10 secs
  
  png("pictures/genomics/GPSC2/tree_roottotip.png", width = 24, height = 12, unit = "cm", res = 1200)
  par(mfrow = c(1,1), mar = c(3, 3, 2, 2), mgp = c(1.7, 0.7, 0), bty = "n")
  res_roottotip <- BactDating::roottotip(rooted_tree,d[,1])
  saveRDS(res_roottotip, "outputs/genomics/GPSC2/res_roottotip.rds")
  dev.off()
  
}

