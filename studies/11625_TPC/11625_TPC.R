library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11625_TPC/44930_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 3690 taxa and 160 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)


df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11625_TPC/11625_20180418-110120.txt', header=TRUE, sep='\t')

############ GLOBAL VARIABLE ##############
STUDY_ID = 11625
SUMMARY = "TPC"
WRITE = F
##########################################

show_plot(otu, df_sample[df_sample$day==35,], "sample_name", "treatment", "TSS", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$day==35,], "sample_name", "treatment", "GMPR", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$day==35,], "sample_name", "treatment", "GMPR2", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$day==35,], "sample_name", "treatment", "CSS", "chemicals", TRUE, write = WRITE)
# 
# show_plot(otu, df_sample, "sample_name", "treatment", "TSS", "chemicals", TRUE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "treatment", "GMPR", "chemicals", TRUE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "treatment", "GMPR2", "chemicals", TRUE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "treatment", "CSS", "chemicals", TRUE, write = WRITE)

