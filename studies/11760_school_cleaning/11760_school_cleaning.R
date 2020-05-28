library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11760_school_cleaning/53224_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 1242 taxa and 26 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11760_school_cleaning/11760_20180514-123250.txt', header=TRUE, sep='\t')
df_sample$time = sapply(strsplit(as.character(df_sample$description), ","),`[[`, 4)
############ GLOBAL VARIABLE ##############
STUDY_ID = 11760
SUMMARY = "school_cleaning"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "host_subject_id", "TSS", "geo_loc", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_subject_id", "GMPR2", "geo_loc", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_subject_id", "CSS", "geo_loc", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "time", "TSS", "time", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "time", "GMPR2", "time", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "time", "CSS", "time", TRUE, write = WRITE)


