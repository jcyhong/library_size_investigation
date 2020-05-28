library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11653_macaque_vacc/47024_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 3911 taxa and 70 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11653_macaque_vacc/11653_20180301-112907.txt', header=TRUE, sep='\t')
df_sample$vacc = sapply(strsplit(as.character(df_sample$description), "[.]"),`[[`, 1)
############ GLOBAL VARIABLE ##############
STUDY_ID = 11653
SUMMARY = "macaque_vacc"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "host_life_stage", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_life_stage", "GMPR2", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_life_stage", "CSS", "demographics", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "vacc", "TSS", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "vacc", "GMPR2", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "vacc", "CSS", "chemicals", TRUE, write = WRITE)

