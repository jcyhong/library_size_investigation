library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11336_multiple_sclerosis/47173_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 4068 taxa and 54 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11336_multiple_sclerosis/11336_20170825-194934.txt', header=TRUE, sep='\t')
df_sample = df_sample[df_sample$host_common_name=="human",]
############ GLOBAL VARIABLE ##############
STUDY_ID = 11336
SUMMARY = "multiple_sclerosis"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "geo_loc_name", "TSS", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "geo_loc_name", "GMPR2", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "geo_loc_name", "CSS", "geo_loc", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "disease_state", "TSS", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "disease_state", "GMPR2", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "disease_state", "CSS", "disease", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "sex", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "GMPR2", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "CSS", "demographics", TRUE, write = WRITE)

