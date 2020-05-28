library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")

source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11488_cheese_rind/42238_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 2516 taxa and 361 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11488_cheese_rind/11488_20180802-123901.txt', header=TRUE, sep='\t')

############ GLOBAL VARIABLE ##############
STUDY_ID = 11488
SUMMARY = "cheese_rinds"
WRITE = FALSE
##########################################

show_plot(otu, df_sample, "sample_name", "continent", "TSS", "geo_loc", TRUE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "continent", "GMPR", "food", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "continent", "GMPR2", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "continent", "CSS", "geo_loc", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "pasteurized", "TSS", "treatment", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "pasteurized", "GMPR2", "treatment", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "pasteurized", "CSS", "treatment", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "animal", "TSS", "", TRUE, write = WRITE, log = TRUE)
show_plot(otu, df_sample, "sample_name", "animal", "GMPR2", "", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "animal", "CSS", "", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "rindtype", "TSS", "", TRUE, write = WRITE, log = TRUE)
show_plot(otu, df_sample, "sample_name", "rindtype", "GMPR2", "", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "rindtype", "CSS", "", TRUE, write = WRITE)

