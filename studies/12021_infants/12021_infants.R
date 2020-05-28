library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(GMPR)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")


path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/12021_infants/64777_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 2052 rows and 294 samples (cols)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/12021_infants/12021_20181203-132340.txt', header=TRUE, sep='\t')


############ GLOBAL VARIABLE ##############
STUDY_ID = 12021
SUMMARY = "infants"
WRITE = T
##########################################
show_plot(otu, df_sample, "sample_name", "diet", "TSS", "food", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diet", "GMPR2", "food", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diet", "CSS", "food", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "gender", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "gender", "GMPR2", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "gender", "CSS", "demographics", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "host_age", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_age", "GMPR2", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_age", "CSS", "demographics", FALSE, write = WRITE)
