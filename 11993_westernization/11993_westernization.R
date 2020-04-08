library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(GMPR)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11993_westernization/62384_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 4834 taxa and 441 samples
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]


df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11993_westernization/11993_20181011-082238.txt', header=TRUE, sep='\t')

############ GLOBAL VARIABLE ##############
STUDY_ID = 11993
SUMMARY = "westernization"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "age_range", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "age_range", "GMPR", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "age_range", "CSS", "demographics", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "bmi_class", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "bmi_class", "GMPR", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "bmi_class", "CSS", "demographics", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "city", "TSS", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "city", "GMPR", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "city", "CSS", "geo_loc", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "medicament", "TSS", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "medicament", "GMPR", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "medicament", "CSS", "disease", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "sex", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "GMPR", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "CSS", "demographics", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "stool_consistency", "TSS", "body_site", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "stool_consistency", "GMPR", "body_site", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "stool_consistency", "CSS", "body_site", FALSE, write = WRITE)

