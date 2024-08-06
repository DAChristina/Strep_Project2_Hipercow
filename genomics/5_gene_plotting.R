library(tidyverse)
library(ggtree)
library(ape)c
library(ggstance)
library(magrittr)
library(vcd)
tre_names <- read.csv("raw_data/tree_inputs_final.csv")
gpsc31.tre <- ape::read.tree("raw_data/gubbins/n703/n703_.node_labelled.final_tree.tre")
gpsc31.tre$tip.label<-gsub("-[12]$","",gpsc31.tre$tip.label)
gg_31 <- ggtree(gpsc31.tre)

library(gggenes)
library(cowplot)

comcde_genes.df <-
  data.frame(
    "start" = c(1,750,2096),
    "end" = c(753,2075,2221),
    "strand" = c(-1,-1,-1),
    "name" = c("comE","comD","comC"),
    "genome" = "gpsc31"
  )

comCDE_plot <-
  ggplot(comcde_genes.df,
         aes(xmin = start,
             xmax = end,
             y = genome,
             forward = strand,
             label = name)) +
  geom_gene_arrow(fill = "steelblue") +
  geom_gene_label(align = "centre",
                  colour = "white") +
  ylab("") +
  #theme_genes() +
  #theme(axis.text.y=element_blank())
  xlim(c(1,2221)) +
  theme_void()

comCDE_lower <-
  cowplot::plot_grid(plotlist = list(NULL,comCDE_plot))

comcde.df <- tre_names %>%
  dplyr::select(ID, MLST_ST, clade, contains("com")) %>% 
  dplyr::rename(Start = qstart_blast_comCDE,
                End = qend_blast_comCDE) %>% 
  dplyr::mutate(Length = End - Start + 1)

gpsc31_comCDE <-
  gg_31 + geom_facet(panel = "comCDE", data = comcde.df, geom = geom_errorbarh, 
                     mapping=aes(x = Start, xmin = Start, xmax = End, colour = Length)) +
  xlim_expand(c(1,2221), "comCDE") +
  theme_void() +
  theme(legend.position = "top")

gpsc31_comCDE_analysis <-
  cowplot::plot_grid(plotlist = list(gpsc31_comCDE,comCDE_lower),
                     ncol = 1,
                     rel_heights = c(0.95,0.05))

ggsave(gpsc31_comCDE_analysis,
       file = "~/Documents/12F_analysis/com_gene_analysis/gpsc31_comCDE_analysis.pdf",
       height = 8,
       width = 8)

gg_31 + geom_facet(panel = "pipR", data = pipr.df, geom = geom_errorbarh, 
                   mapping=aes(x = Start, xmin = Start, xmax = End, colour = Length))

gpsc32_comCDE <-
  gg_32 + geom_facet(panel = "comCDE", data = comcde.df, geom = geom_errorbarh, 
                     mapping=aes(x = Start, xmin = Start, xmax = End, colour = Length)) +
  xlim_expand(c(1,2221), "comCDE") +
  theme_void() +
  theme(legend.position = "top")

gpsc32_comCDE_analysis <-
  cowplot::plot_grid(plotlist = list(gpsc32_comCDE,comCDE_lower),
                     ncol = 1,
                     rel_heights = c(0.95,0.05))

ggsave(gpsc32_comCDE_analysis,
       file = "~/Documents/12F_analysis/com_gene_analysis/GPSC32_comCDE_analysis.pdf",
       height = 8,
       width = 8)

gg_32 + geom_facet(panel = "pipR", data = pipr.df, geom = geom_errorbarh, 
                   mapping=aes(x = Start, xmin = Start, xmax = End, colour = Length))

gpsc55_comCDE <-
  gg_55 + geom_facet(panel = "comCDE", data = comcde.df, geom = geom_errorbarh, 
                     mapping=aes(x = Start, xmin = Start, xmax = End, colour = Length)) +
  xlim_expand(c(1,2221), "comCDE") +
  theme_void() +
  theme(legend.position = "top")

gpsc55_comCDE_analysis <-
  cowplot::plot_grid(plotlist = list(gpsc55_comCDE,comCDE_lower),
                     ncol = 1,
                     rel_heights = c(0.95,0.05))

ggsave(gpsc55_comCDE_analysis,
       file = "~/Documents/12F_analysis/com_gene_analysis/GPSC55_comCDE_analysis.pdf",
       height = 8,
       width = 8)

gg_55 + geom_facet(panel = "pipR", data = pipr.df, geom = geom_errorbarh, 
                   mapping=aes(x = Start, xmin = Start, xmax = End, colour = Length))

##########

# ICE locus affected by recombination in the GPSC55 reference
# locus 271137..313269

rec.df <-
  read.table(file = "~/Documents/12F_analysis/project_gubbins/gubbins_55.recombination_predictions.gff",
             #skip = 1,
             sep = "\t",
             header = FALSE)

rec_count <- rep(0,times = max(rec.df[,5]))

for (i in 1:nrow(rec.df)) {
  for (j in rec.df[i,4]:rec.df[i,5]) {
    rec_count[j] = rec_count[j] + 1
  }
}

rec_levels.df <-
  data.frame(
    "Position" = 1:length(rec_count),
    "Recombination" =  rec_count
  )

ggplot(rec_levels.df,
       aes(x = Position, y = Recombination)) +
  geom_line() +
  xlim(c(271137,313269)) +
  geom_vline(xintercept = 26992+271137) +
  ylab("Number of recombination events affecting base")

library(genoPlotR)

min_seq_id<-75
max_seq_id<-100

range_bounds <- round(seq(from=min_seq_id,to=max_seq_id,length.out=10))
reds <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", 
          "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")

process_comparison_file <- function(comparison,range_bounds,reds) {
  
  ref.query.comparison.df<-read.table(comparison, comment.char="")
  colnames(ref.query.comparison.df) <- c("name1", "name2", "per_id", "aln_len", "mism", 
                                         "gaps", "start1", "end1", "start2", "end2", "e_value", 
                                         "bit_score")
  ref.query.comparison.df <- ref.query.comparison.df[, c(colnames(ref.query.comparison.df)[c(7:10, 1:6, 11:12)])]
  ref.query.comparison.df<-ref.query.comparison.df[order(ref.query.comparison.df$aln_len),]
  ref.query.comparison.actOrder<-as.comparison(ref.query.comparison.df)
  
  # Reverse comparison file from ACT order
  ref.query.comparison<-ref.query.comparison.actOrder
  ref.query.comparison$start1<-ref.query.comparison.actOrder$start2
  ref.query.comparison$end1<-ref.query.comparison.actOrder$end2
  ref.query.comparison$start2<-ref.query.comparison.actOrder$start1
  ref.query.comparison$end2<-ref.query.comparison.actOrder$end1
  
  # Add colours
  col.df <-
    data.frame(
      "rounded_id" = range_bounds[2:10],
      "col" = reds
    )
  
  ref.query.comparison %<>%
    dplyr::filter(per_id >= min(range_bounds)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(rounded_id = col.df$rounded_id[which.min(abs(col.df$rounded_id - per_id))]) %>%
    ungroup() %>%
    dplyr::left_join(col.df, by = c("rounded_id" = "rounded_id"))
  
  return(as.comparison(ref.query.comparison))
}

gpsc55_ice_locus <- read_dna_seg_from_embl(
  "~/Documents/12F_analysis/ICE_analysis/GPSC55_ICE_locus.embl",
  tagsToParse = c("CDS"))
gpsc31_ice_locus <- read_dna_seg_from_embl(
  "~/Documents/12F_analysis/ICE_analysis/gpsc31_ICE_locus.embl",
  tagsToParse = c("CDS"))

gpsc55_ice_locus$gene_type[gpsc55_ice_locus$feature=="CDS"]<-"arrows"
gpsc31_ice_locus$gene_type[gpsc31_ice_locus$feature=="CDS"]<-"arrows"

gpsc55_ice_annotation<-auto_annotate(gpsc55_ice_locus,
                                     locus_tag_pattern=NULL,
                                     names=gpsc55_ice_locus$gene,
                                     rot=45)
gpsc55_ice_locus$fill<-gpsc55_ice_locus$col
gpsc31_ice_annotation<-auto_annotate(gpsc31_ice_locus,
                                     locus_tag_pattern=NULL,
                                     names=gpsc31_ice_locus$gene,
                                     rot=45)
gpsc31_ice_locus$fill<-gpsc31_ice_locus$col

seq_names <- c("GPS_IL_27095","GPS_ZA_1974")

ice_locus_rec_levels.df <-
  rec_levels.df %>%
  dplyr::filter(Position >= 271137 & Position <= 313269) %>%
  dplyr::mutate(Position = Position - 271137 + 1)

gpsc55_seg_plot <-
  seg_plot(func = linesGrob,
           args=list(x=ice_locus_rec_levels.df$Position,
                     y=ice_locus_rec_levels.df$Recombination))

ice_locus_comparison <- process_comparison_file("~/Documents/12F_analysis/ICE_analysis/GPSC55_ICE_locus.dna.gpsc31_ICE_locus.dna.crunch",
                                                range_bounds,
                                                reds)

# Start saving
output_fn<-"~/Documents/12F_analysis/ICE_analysis/12F_ICE_comparison.pdf"
pdf(file=output_fn,width=6,height=5)

plot_gene_map(list(gpsc55_ice_locus,
                   gpsc31_ice_locus),
              comparisons = list(ice_locus_comparison),
              seg_plots = list(gpsc55_seg_plot,NULL),
              annotations = list(gpsc55_ice_annotation,
                                 gpsc31_ice_annotation),
              offsets = c(0,(max(gpsc55_ice_locus$end)-max(gpsc31_ice_locus$end))/2),
              annotation_cex = 1,
              annotation_height = 0.0,
              dna_seg_label_cex = 1.0,
              dna_seg_labels = seq_names,
              top_annotation_banner = 3.0)

# Generate legend
range_bounds <- round(seq(from=min_seq_id,to=max_seq_id,length.out=10))
reds <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", 
          "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")

legend_text<-NULL
for (i in 1:9) {
  legend_text<-c(legend_text,
                 paste0(range_bounds[i],"-",range_bounds[i+1]))
}

suppressWarnings(
  invisible(capture.output(
    pushViewport(viewport("left",
                          height = .9,
                          width = .9,
    ))
  ))
)

grid_legend(
  x = "left",
  labels = legend_text,
  title = "Identity (%)",
  col = reds,
  size = 2,
  pch = 15,
  draw = TRUE,
  frame = FALSE,
  vgap = 0,
  just = "left",
  inset = c(0,0),
  gp=gpar(cex=.85)
)

# Save
dev.off()

##########

seq_names <- c("GPS_IL_27095","467530_H18014004907")

process_sequence_file <- function(fn, tags = c("CDS")) {
  seq_locus <- read_dna_seg_from_embl(fn, tagsToParse = tags)
  seq_locus$gene_type[seq_locus$feature=="CDS"]<-"arrows"
  seq_annotation<-auto_annotate(seq_locus,
                                locus_tag_pattern=NULL,
                                names=seq_locus$gene,
                                rot=45)
  seq_locus$fill<-seq_locus$col
  return(list(seq_locus,seq_annotation))
}

com_ref_locus <- process_sequence_file("~/Documents/12F_analysis/com_gene_analysis/comparisons/gpsc31_comCDE.embl")
com_467530_locus <- process_sequence_file("~/Documents/12F_analysis/com_gene_analysis/comparisons/467530_H18014004907_comCDE.embl",
                                          tags = c("CDS","repeat_region","misc_feature"))

com_467530_locus[[1]] %<>%
  dplyr::mutate(col = dplyr::case_when(
    gene == "IRL" ~ "green",
    gene == "IRR" ~ "green4",
    TRUE ~ col
  )
  ) %>%
  dplyr::mutate(fill = col) %>%
  #dplyr::filter(feature != "CDS_intron_pseudo") %>%
  dplyr::mutate(gene_type = dplyr::case_when(
    feature == "CDS_pseudo" & end == 63 ~ "arrows",
    feature == "CDS_pseudo" & end == 1919 ~ "arrows",
    TRUE ~ gene_type
  )
  )

com_ref_locus[[2]] <-
  as.annotation(
    data.frame(
      # "x1" = c(1,750,2096),
      # "x2" = c(753,2075,2221),
      "x1" = c(377,1412.5,2158.5),
      "x2" = c(377,1412.5,2158.5),
      "text" = c("comE","comD","comC"),
      "color" = "black",
      "rot" = 45
    )
  )

com_locus_comparison <- process_comparison_file("~/Documents/12F_analysis/com_gene_analysis/comparisons/gpsc31_comCDE.dna.467530_H18014004907_comCDE.fa.crunch",
                                                range_bounds,
                                                reds)

# Start saving
output_fn<-"~/Documents/12F_analysis/com_gene_analysis/com_467530_comparison.pdf"
pdf(file=output_fn,width=6,height=5)

plot_gene_map(list(com_ref_locus[[1]],
                   com_467530_locus[[1]]),
              comparisons = list(com_locus_comparison),
              annotations = list(com_ref_locus[[2]],
                                 NULL),
              #offsets = c(0,(max(gpsc55_ice_locus$end)-max(gpsc31_ice_locus$end))/2),
              annotation_cex = 1,
              annotation_height = 0.0,
              dna_seg_label_cex = 1.0,
              dna_seg_labels = seq_names,
              top_annotation_banner = 3.0)

# Generate legend
range_bounds <- round(seq(from=min_seq_id,to=max_seq_id,length.out=10))
reds <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", 
          "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")

legend_text<-NULL
for (i in 1:9) {
  legend_text<-c(legend_text,
                 paste0(range_bounds[i],"-",range_bounds[i+1]))
}

suppressWarnings(
  invisible(capture.output(
    pushViewport(viewport("left",
                          height = .9,
                          width = .9,
    ))
  ))
)

grid_legend(
  x = "left",
  labels = legend_text,
  title = "Identity (%)",
  col = reds,
  size = 2,
  pch = 15,
  draw = TRUE,
  frame = FALSE,
  vgap = 0,
  just = "left",
  inset = c(0,0),
  gp=gpar(cex=.85)
)

# Save
dev.off()

#############

com_517608_locus <- process_sequence_file("~/Documents/12F_analysis/com_gene_analysis/comparisons/517608_H18142019706-2_comCDE.embl",
                                          tags = c("CDS","repeat_region","misc_feature"))

second_com_locus_comparison <- process_comparison_file("~/Documents/12F_analysis/com_gene_analysis/comparisons/gpsc31_comCDE.dna.517608_H18142019706-2_comCDE.dna.crunch",
                                                       range_bounds,
                                                       reds)

second_output_fn<-"~/Documents/12F_analysis/com_gene_analysis/com_517608_comparison.pdf"
pdf(file=second_output_fn,width=6,height=5)

second_seq_names<-c(seq_names[1],"517608_H18142019706")

plot_gene_map(list(com_ref_locus[[1]],
                   com_517608_locus[[1]]),
              comparisons = list(second_com_locus_comparison),
              annotations = list(com_ref_locus[[2]],
                                 NULL),
              offsets = c(0,(max(com_ref_locus[[1]]$end)-max(com_517608_locus[[1]]$end))/2),
              annotation_cex = 1,
              annotation_height = 0.0,
              dna_seg_label_cex = 1.0,
              dna_seg_labels = second_seq_names,
              top_annotation_banner = 3.0)

# Generate legend
range_bounds <- round(seq(from=min_seq_id,to=max_seq_id,length.out=10))
reds <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", 
          "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")

legend_text<-NULL
for (i in 1:9) {
  legend_text<-c(legend_text,
                 paste0(range_bounds[i],"-",range_bounds[i+1]))
}

suppressWarnings(
  invisible(capture.output(
    pushViewport(viewport("left",
                          height = .9,
                          width = .9,
    ))
  ))
)

grid_legend(
  x = "left",
  labels = legend_text,
  title = "Identity (%)",
  col = reds,
  size = 2,
  pch = 15,
  draw = TRUE,
  frame = FALSE,
  vgap = 0,
  just = "left",
  inset = c(0,0),
  gp=gpar(cex=.85)
)

# Save
dev.off()

#############

com_689707_locus <- process_sequence_file("~/Documents/12F_analysis/com_gene_analysis/comparisons/689707_H19074023606-1_comCDE.embl",
                                          tags = c("CDS","repeat_region","misc_feature"))

com_689707_locus[[1]] %<>%
  dplyr::mutate(col = dplyr::case_when(
    gene == "IRL" ~ "green",
    gene == "IRR" ~ "green4",
    TRUE ~ col
  )
  ) %>%
  dplyr::mutate(fill = col) %>%
  #dplyr::filter(feature != "CDS_intron_pseudo") %>%
  dplyr::mutate(gene_type = dplyr::case_when(
    feature == "CDS" & end == 935 ~ "exons",
    TRUE ~ gene_type
  )
  )

third_com_locus_comparison <- process_comparison_file("~/Documents/12F_analysis/com_gene_analysis/comparisons/gpsc31_comCDE.dna.689707_H19074023606-1_comCDE.fa.crunch",
                                                      range_bounds,
                                                      reds)

third_output_fn<-"~/Documents/12F_analysis/com_gene_analysis/com_689707_comparison.pdf"
pdf(file=third_output_fn,width=6,height=5)

third_seq_names<-c(seq_names[1],"689707_H19074023606")

plot_gene_map(list(com_ref_locus[[1]],
                   com_689707_locus[[1]]),
              comparisons = list(third_com_locus_comparison),
              annotations = list(com_ref_locus[[2]],
                                 NULL),
              offsets = c(0,(max(com_ref_locus[[1]]$end)-max(com_689707_locus[[1]]$end))/2),
              annotation_cex = 1,
              annotation_height = 0.0,
              dna_seg_label_cex = 1.0,
              dna_seg_labels = third_seq_names,
              top_annotation_banner = 3.0)

# Generate legend
range_bounds <- round(seq(from=min_seq_id,to=max_seq_id,length.out=10))
reds <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", 
          "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")

legend_text<-NULL
for (i in 1:9) {
  legend_text<-c(legend_text,
                 paste0(range_bounds[i],"-",range_bounds[i+1]))
}

suppressWarnings(
  invisible(capture.output(
    pushViewport(viewport("left",
                          height = .9,
                          width = .9,
    ))
  ))
)

grid_legend(
  x = "left",
  labels = legend_text,
  title = "Identity (%)",
  col = reds,
  size = 2,
  pch = 15,
  draw = TRUE,
  frame = FALSE,
  vgap = 0,
  just = "left",
  inset = c(0,0),
  gp=gpar(cex=.85)
)

# Save
dev.off()

#############

com_777372_locus <- process_sequence_file("~/Documents/12F_analysis/com_gene_analysis/comparisons/777372_H19288054206-2_comCDE.embl",
                                          tags = c("CDS","repeat_region","misc_feature"))

com_777372_locus[[1]] %<>%
  dplyr::mutate(col = dplyr::case_when(
    gene == "IRL" ~ "green",
    gene == "IRR" ~ "green4",
    TRUE ~ col
  )
  ) %>%
  dplyr::mutate(fill = col) %>%
  #dplyr::filter(feature != "CDS_intron_pseudo") %>%
  dplyr::mutate(gene_type = dplyr::case_when(
    feature == "CDS" & end == 2246 ~ "exons",
    TRUE ~ gene_type
  )
  )

fourth_com_locus_comparison <- process_comparison_file("~/Documents/12F_analysis/com_gene_analysis/comparisons/gpsc31_comCDE.dna.777372_H19288054206-2_comCDE.fa.crunch",
                                                       range_bounds,
                                                       reds)

fourth_output_fn<-"~/Documents/12F_analysis/com_gene_analysis/com_777372_comparison.pdf"
pdf(file=fourth_output_fn,width=6,height=5)

fourth_seq_names<-c(seq_names[1],"777372_H19288054206")

plot_gene_map(list(com_ref_locus[[1]],
                   com_777372_locus[[1]]),
              comparisons = list(fourth_com_locus_comparison),
              annotations = list(com_ref_locus[[2]],
                                 NULL),
              offsets = c(0,(max(com_ref_locus[[1]]$end)-max(com_777372_locus[[1]]$end))/2),
              annotation_cex = 1,
              annotation_height = 0.0,
              dna_seg_label_cex = 1.0,
              dna_seg_labels = fourth_seq_names,
              top_annotation_banner = 3.0)

# Generate legend
range_bounds <- round(seq(from=min_seq_id,to=max_seq_id,length.out=10))
reds <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", 
          "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")

legend_text<-NULL
for (i in 1:9) {
  legend_text<-c(legend_text,
                 paste0(range_bounds[i],"-",range_bounds[i+1]))
}

suppressWarnings(
  invisible(capture.output(
    pushViewport(viewport("left",
                          height = .9,
                          width = .9,
    ))
  ))
)

grid_legend(
  x = "left",
  labels = legend_text,
  title = "Identity (%)",
  col = reds,
  size = 2,
  pch = 15,
  draw = TRUE,
  frame = FALSE,
  vgap = 0,
  just = "left",
  inset = c(0,0),
  gp=gpar(cex=.85)
)

# Save
dev.off()


#############

clade.list<-c("506090_H18106005606-1.fa","680082_H19048054409-1.fa","457406_H17502016106-1.fa","516354_H18132048706-1.fa","435328_H17420042303-1.fa","447693_H17470016406-1.fa","564307_H18262032806-2.fa","579904_H18304025806-1.fa","448288_H17472039206-2.fa","456766_H17500014006-1.fa","678111_H19042016806-1.fa","441564_H17444022006-1.fa","773375_H19274042506-2.fa","520183_H18146028206-1.fa","642720_H18466033806-2.fa","456249_H17496004206-2.fa","459989_H17508005506-1.fa","838214_H19454015407-1.fa","220448_H16046042510-1.fa","220449_H16052071010-1.fa","872204_H20024010606-2.fa","848345_H19478008206-1.fa","441572_H17444023006-1.fa","669357_H19022053806-2.fa","853095_H19490006106-2.fa","467530_H18014004907-1.fa","470371_H18022024106-1.fa","773223_H19272023006-1.fa","589176_H18326024005-2.fa","636516_H18448023906-2.fa","848408_H19482050306-1.fa","888044_H20068035906-1.fa","619126_H18396022006-1.fa","417343_H17364003309-1.fa","838225_H19454055407-1.fa","646235_H18474032306-1.fa","747884_H19208019310-1.fa","526243_H18158051106-1.fa","667355_H19020016406-1.fa","661897_H18520017710-1.fa","665069_H18528050308-1.fa","779162_H19292026306-2.fa","471050_H18026023606-2.fa","886798_H20066042706-2.fa","482144_H18038039506-1.fa","861153_H19512012706-1.fa","526234_H18158025606-1.fa","408621_H17336011005-1.fa","444132_H17456047606-1.fa","462693_H17516015306-1.fa","444455_H17456047507-2.fa","560704_H18256031306-2.fa","511687_H18122037206-1.fa","660129_H18514030806-2.fa","869355_H20018122906-2.fa","429124_H17402034402-2.fa","1016800.fa","444499_H17458032806-2.fa","869249_H20018012506-1.fa","523826_H18152023806-2.fa","692657_H19082052306-1.fa","471092_H18026037306-2.fa","694520_H19088003706-2.fa","524256_H18148025011-1.fa","692611_H19082030306-1.fa","724688_H19144025806-2.fa","1142515_H21112083807-1.fa","886747_H20062080006-2.fa","822728_H19414020606-2.fa","448291_H17472041006-2.fa","619132_H18396022706-1.fa","619093_H18392038006-1.fa","506447_H18110032406-2.fa","663173_H19010051307-1.fa","537116_H18192018506-1.fa","461404_H17512022206-2.fa","950383_H20218054006-2.fa","984473_H20364064106-1.fa","495576_H18086029006-2.fa","1013481.fa","838201_H19452057107-1.fa","448287_H17472039106-2.fa","789534_H19316042508-2.fa","560779_H18252042306-1.fa","839968_H19456031308-2.fa","857501_H19502066106-1.fa","452431_H17484027006-1.fa","454005_H17492029706-2.fa","447694_H17470016506-1.fa","758107_H19236025306-1.fa","1017804.fa","1019854.fa","1019854_H160520723-2.fa")

blast.df <-
  read.table("/Users/ncrouche/Documents/12F_analysis/com_gene_analysis/com_gene_blast.blast.out", header = F)

colnames(blast.df) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

dist_between_matches.df <-
  blast.df %>%
  dplyr::filter(sseqid %in% clade.list) %>%
  dplyr::group_by(sseqid) %>%
  dplyr::summarise(Max = max(sstart),
                   Min = min(sstart)
  ) %>%
  dplyr::mutate(Range = Max - Min)

