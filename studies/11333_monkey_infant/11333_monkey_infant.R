library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")

source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11333_monkey_infant/47171_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 3539 taxa and 58 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11333_monkey_infant/11333_20180418-110057.txt', header=TRUE, sep='\t')

############ GLOBAL VARIABLE ##############
STUDY_ID = 11333
SUMMARY = "monkey_infant"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "assigned_dietary_group", "TSS", "food", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "assigned_dietary_group", "GMPR2", "food", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "assigned_dietary_group", "CSS", "food", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "diet", "TSS", "food", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diet", "GMPR2", "food", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diet", "CSS", "food", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "age", "TSS", "time", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "age", "GMPR2", "time", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "age", "CSS", "time", TRUE, write = WRITE)

